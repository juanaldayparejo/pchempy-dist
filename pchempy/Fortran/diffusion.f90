MODULE diffusion
 
  CONTAINS
  
!==========================================================================================================

  subroutine solve_ode_rosenbrock(nlay,ngas,h,T,N,mmol,A,s,B,&
    physrates,typelbc,valuelbc,typeubc,valueubc,fix_species,dt,Nnew,err)

    !  Function to calculate the density at the next time step using a second-order Rosenbrock solver
    !
    !  Using the method from Verwer et al. (1999) and Hobbs et al. (2019)

    !use omp_lib
    implicit none
    

    !Inputs
    real, intent(in) :: h(nlay)                           !Altitude (m)
    real, intent(in) :: T(nlay)                           !Temperature (K)
    double precision, intent(in) :: N(nlay,ngas)          !Number density of each species (m-3)
    real, intent(in) :: mmol(ngas)                        !Molecular weight (uma)
    real, intent(in) :: A(ngas),s(ngas),B(ngas)           !Coefficients for molecular diffusion
    double precision, intent(in) :: physrates(nlay,ngas)  !Production/loss rates from other processes (such as chemistry) (m-3 s-1)
    double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)     !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
    integer, intent(in) :: typelbc(ngas),typeubc(ngas)    !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
    integer, intent(in) :: nlay,ngas
    integer, intent(in) :: fix_species(nlay,ngas)         !Flag to indicate if density must be fixed (if 1)
    double precision, intent(in) :: dt                                !Timestep (s)

    !Local 
    integer :: igas,ilay,jlay
    double precision :: Jmat(nlay,nlay,ngas),N0(nlay)
    double precision :: b1(nlay,ngas),g1(nlay,ngas),gi(nlay)
    double precision :: g2(nlay,ngas),b2(nlay,ngas),a2(nlay,ngas)
    double precision :: Mmat(nlay,nlay),gamma
    real :: K(nlay)                             !Eddy diffusion coefficient (m2 s-1)
    real :: D(nlay,ngas)                        !Molecular diffuson coefficient (m2 s-1)
    real :: mmean(nlay)                         !Mean molecular weight (uma)
    real :: scaleH0(nlay),scaleH(nlay,ngas)     !Scale height (m)
    double precision :: ksi(nlay,ngas),klsi(nlay,ngas),ksim1(nlay,ngas),klsim1(nlay,ngas)
    real :: delz
    integer :: code, indx(nlay)

    integer :: num_threads

    external dgesv

    !Output
    double precision, intent(out) :: Nnew(nlay,ngas)
    double precision, intent(out) :: err(nlay,ngas)

    !Performing calculations to characterise the atmosphere
    !############################################################

    !Calculating the total density (m-3)
    do ilay=1,nlay
      N0(ilay) = 0.d0
      do igas=1,ngas
        N0(ilay) = N0(ilay) + N(ilay,igas)
      enddo
    enddo

    delz = h(2) - h(1)  !W

    !Calculating the Eddy diffusion coefficient in each layer
    call calc_Keddy(nlay,h,N0,K)

    !Calculating the molecular diffusion coefficient in each layer
    call calc_Dmoldiff(nlay,ngas,N0,T,A,s,D)

    !Calculating the mean molecular weight
    call calc_mmean(nlay,ngas,N,mmol,mmean)

    !Calculating the scale height
    call calc_scaleH_mean(nlay,h,T,mmean,scaleH0)

    !Calculating the species-dependent scale height
    call calc_scaleH(nlay,ngas,h,T,mmol,scaleH)


    !Calculating the diffusion parameters
    !################################################################

    !Calculating the diffusion coefficients to fill the jacobian 
    call calc_diffusion_coeff(nlay,ngas,h,T,scaleH0,scaleH,K,D,B,ksi,klsi,ksim1,klsim1)

    !Calculating the Jacobian matrix
    call calc_jacobian_diffusion(nlay,ngas,ksi,klsi,ksim1,klsim1,typelbc,valuelbc,typeubc,valueubc,fix_species,Jmat)

    !do ilay=1,nlay
    !  do jlay=1,nlay
    !    print*,Jmat(ilay,jlay,1)
    !  enddo
    !enddo

    !Evaluate the jacobian function at the current density
    call eval_fn_diffusion(nlay,ngas,delz,N,Jmat,physrates,typelbc,valuelbc,typeubc,valueubc,fix_species,b1)

    !do ilay=1,nlay
    !  print*,b1(ilay,1)
    !enddo

    !Rosenbrock solver (see Hobbs et al., 2019)
    !################################################################

    !Build matrix (I-gamma*dt*J) and solve for g1
    gamma = 1.d0 + 1.d0/sqrt(2.d0)
    

    !!$OMP PARALLEL PRIVATE(igas,ilay,jlay,Mmat,gi,indx,code) SHARED(gamma,dt,Jmat,g1,nlay,b1,ngas) NUM_THREADS(8)
    !!$OMP DO 

    do igas=1,ngas

      !print*,igas

      !Building matrix
      do ilay=1,nlay
        do jlay=1,nlay
          if(ilay.eq.jlay)then
            Mmat(ilay,jlay) = 1.0d0 - gamma*dt*Jmat(ilay,jlay,igas)
          else
            Mmat(ilay,jlay) = -gamma*dt*Jmat(ilay,jlay,igas)
          endif
        enddo
      enddo

      !Solve for g1   (I-gamma*dt*J) * g1 = b1
      gi(:) = b1(:,igas)
      call dgesv(nlay,1,Mmat,nlay,indx,gi,nlay,code)

      g1(:,igas) = gi(:)

    enddo
    !!$OMP END DO 
    !!$OMP END PARALLEL

    !Calculate b2 = F(N+dt*g1)
    !!$OMP PARALLEL PRIVATE(igas,ilay) SHARED(nlay,ngas,N,dt,g1,a2)
    !!$OMP DO 
    do igas=1,ngas
      do ilay=1,nlay
        a2(ilay,igas) = N(ilay,igas)+dt*g1(ilay,igas)
      enddo
    enddo
    !!$OMP END DO 
    !!$OMP END PARALLEL

    call eval_fn_diffusion(nlay,ngas,delz,a2,Jmat,physrates,typelbc,valuelbc,typeubc,valueubc,fix_species,b2)

    !Building matrix (I-gamma*dt*J) and solve for g2
    !!$OMP PARALLEL PRIVATE(igas,ilay,jlay,Mmat,gi,indx,code) SHARED(gamma,dt,Jmat,g1,nlay,b2,ngas,g2) NUM_THREADS(8)
    !!$OMP DO 
    do igas=1,ngas

      !Building matrix
      do ilay=1,nlay
        do jlay=1,nlay
          if(ilay.eq.jlay)then
            Mmat(ilay,jlay) = 1.0d0 - gamma*dt*Jmat(ilay,jlay,igas)
          else
            Mmat(ilay,jlay) = -gamma*dt*Jmat(ilay,jlay,igas)
          endif
        enddo
      enddo

      !Solve for g2   (I-gamma*dt*J) * g2 = b1
      gi(:) = b2(:,igas) - 2.d0*g1(:,igas)
      call dgesv(nlay,1,Mmat,nlay,indx,gi,nlay,code)

      g2(:,igas) = gi(:)

    enddo
    !!$OMP END DO 
    !!$OMP END PARALLEL

    !Finding the new density
    do igas=1,ngas
      do ilay=1,nlay
        Nnew(ilay,igas) = N(ilay,igas) + 1.5*dt*g1(ilay,igas) + 0.5*dt*g2(ilay,igas)

        if(Nnew(ilay,igas).le.1.0d-30)then
          Nnew(ilay,igas) = 0.d0
        endif
      enddo
    enddo


    !Computing the local error
    do igas=1,ngas
      do ilay=1,nlay
        err(ilay,igas) = abs(Nnew(ilay,igas) - (N(ilay,igas) + dt*g1(ilay,igas)))
      enddo
    enddo

  end subroutine

!==========================================================================================================

  subroutine solve_ode_crancknicholson(nlay,ngas,h,T,N,mmol,A,s,B,&
    physrates,typelbc,valuelbc,typeubc,valueubc,fix_species,dt,Nnew)

    !  Function to calculate the density at the next time step using a second-order Rosenbrock solver
    !
    !  Using the method from Verwer et al. (1999) and Hobbs et al. (2019)

    !use omp_lib
    implicit none
    

    !Inputs
    real, intent(in) :: h(nlay)                           !Altitude (m)
    real, intent(in) :: T(nlay)                           !Temperature (K)
    double precision, intent(in) :: N(nlay,ngas)          !Number density of each species (m-3)
    real, intent(in) :: mmol(ngas)                        !Molecular weight (uma)
    real, intent(in) :: A(ngas),s(ngas),B(ngas)           !Coefficients for molecular diffusion
    double precision, intent(in) :: physrates(nlay,ngas)  !Production/loss rates from other processes (such as chemistry) (m-3 s-1)
    double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)     !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
    integer, intent(in) :: typelbc(ngas),typeubc(ngas)    !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
    integer, intent(in) :: nlay,ngas
    integer, intent(in) :: fix_species(nlay,ngas)         !Flag to indicate if density must be fixed (if 1)
    double precision, intent(in) :: dt                                !Timestep (s)

    !Local 
    integer :: igas,ilay,jlay
    double precision :: Jmat(nlay,nlay,ngas),N0(nlay)
    double precision :: b1(nlay,ngas),g1(nlay,ngas),gi(nlay)
    double precision :: g2(nlay,ngas),b2(nlay,ngas),a2(nlay,ngas)
    double precision :: Mmat(nlay,nlay),gamma
    real :: K(nlay)                             !Eddy diffusion coefficient (m2 s-1)
    real :: D(nlay,ngas)                        !Molecular diffuson coefficient (m2 s-1)
    real :: mmean(nlay)                         !Mean molecular weight (uma)
    real :: scaleH0(nlay),scaleH(nlay,ngas)     !Scale height (m)
    double precision :: ksi(nlay,ngas),klsi(nlay,ngas),ksim1(nlay,ngas),klsim1(nlay,ngas)
    real :: delz
    integer :: code, indx(nlay)

    integer :: num_threads

    external dgesv

    !Output
    double precision, intent(out) :: Nnew(nlay,ngas)


    !Performing calculations to characterise the atmosphere
    !############################################################

    !Calculating the total density (m-3)
    do ilay=1,nlay
      N0(ilay) = 0.d0
      do igas=1,ngas
        N0(ilay) = N0(ilay) + N(ilay,igas)
      enddo
    enddo

    delz = h(2) - h(1)  !W

    !Calculating the Eddy diffusion coefficient in each layer
    call calc_Keddy(nlay,h,N0,K)

    !Calculating the molecular diffusion coefficient in each layer
    call calc_Dmoldiff(nlay,ngas,N0,T,A,s,D)

    !Calculating the mean molecular weight
    call calc_mmean(nlay,ngas,N,mmol,mmean)

    !Calculating the scale height
    call calc_scaleH_mean(nlay,h,T,mmean,scaleH0)

    !Calculating the species-dependent scale height
    call calc_scaleH(nlay,ngas,h,T,mmol,scaleH)


    !Calculating the diffusion parameters
    !################################################################

    !Calculating the diffusion coefficients to fill the jacobian 
    call calc_diffusion_coeff(nlay,ngas,h,T,scaleH0,scaleH,K,D,B,ksi,klsi,ksim1,klsim1)

    !Calculating the Jacobian matrix
    call calc_jacobian_diffusion(nlay,ngas,ksi,klsi,ksim1,klsim1,typelbc,valuelbc,typeubc,valueubc,fix_species,Jmat)

    !Evaluating the function to get b
    call calc_matrix_crancknicholson(nlay,ngas,N,physrates,ksi,klsi,ksim1,klsim1,typelbc,valuelbc,typeubc,valueubc,&
                                     fix_species,dt,delz,Jmat,b1)


    !Solving the ODE to get the new density
    !!$OMP PARALLEL PRIVATE(igas,ilay,jlay,Mmat,gi,indx,code) SHARED(gamma,dt,Jmat,g1,nlay,b1,ngas) NUM_THREADS(8)
    !!$OMP DO 
    do igas=1,ngas

      !Building matrix
      do ilay=1,nlay
        do jlay=1,nlay
          if(ilay.eq.jlay)then
            Mmat(ilay,jlay) = Jmat(ilay,jlay,igas)
          else
            Mmat(ilay,jlay) = Jmat(ilay,jlay,igas)
          endif
        enddo
      enddo

      gi(:) = b1(:,igas)
      call dgesv(nlay,1,Mmat,nlay,indx,gi,nlay,code)

      Nnew(:,igas) = gi(:)

    enddo
    !!$OMP END DO 
    !!$OMP END PARALLEL

  end subroutine

!==========================================================================================================
  
    subroutine calc_matrix_crancknicholson(nlay,ngas,N,physrates,ksi,klsi,ksim1,klsim1,&
                                          typelbc,valuelbc,typeubc,valueubc,&
                                          fix_species,deltat,delz,Jmat,b)

      !Function to calculate the jacobian matrix for diffusion
  
      !use omp_lib
  
  
      !Inputs
      double precision, intent(in) :: ksi(nlay,ngas)
      double precision, intent(in) :: klsi(nlay,ngas)
      double precision, intent(in) :: ksim1(nlay,ngas)
      double precision, intent(in) :: klsim1(nlay,ngas)
      double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)   !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
      double precision, intent(in) :: N(nlay,ngas)                    !Number density of each species (m-3)
      double precision, intent(in) :: physrates(nlay,ngas)            !Production/loss rates from other processes (such as chemistry)
      integer, intent(in) :: typelbc(ngas),typeubc(ngas)              !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
      integer, intent(in) :: nlay,ngas
      integer, intent(in) :: fix_species(nlay,ngas)                   !Flag to indicate if density must be fixed (if 1)
      double precision, intent(in) :: deltat                          !Timestep (s)
      real, intent(in) :: delz                                        !Layer size (m)
  
      !Local 
      integer :: igas,ilay
  
      !Output
      double precision, intent(out) :: Jmat(nlay,nlay,ngas)
      double precision, intent(out) :: b(nlay,ngas)
  
      
  
      !Filling the Jacobian matrix
      !!$OMP PARALLEL PRIVATE(igas,ilay,jlay) &
      !!$OMP SHARED(Jmat,b,N,physrates,typelbc,ksi,klsi,ksim1,klsim1,typeubc,deltat,delz) &
      !!$OMP NUM_THREADS(8)
      !!$OMP DO 
      do igas=1,ngas
  
        !Initialising the matrix
        do ilay=1,nlay
          do jlay=1,nlay
            Jmat(ilay,jlay,igas) = 0.d0
          enddo
          b(ilay,igas) = 0.d0
        enddo
  
        !Inbetween layers
        !##############################################
  
        do ilay=2,nlay-1
  
          Jmat(ilay,ilay,igas) = 1.d0 + deltat/2.d0 * (klsi(ilay,igas)+ksim1(ilay,igas))
          Jmat(ilay,ilay-1,igas) = -deltat/2.d0*klsim1(ilay,igas)
          Jmat(ilay,ilay+1,igas) = -deltat/2.d0*ksi(ilay,igas)

          b(ilay,igas) = N(ilay,igas) * (1.d0 - deltat/2.d0 * (klsi(ilay,igas)+ksim1(ilay,igas))) + &
                         N(ilay-1,igas) * deltat/2.d0*klsim1(ilay,igas) + &
                         N(ilay+1,igas) * deltat/2.d0*ksi(ilay,igas) + physrates(ilay,igas) * deltat
  
        enddo

        !Lower boundary
        !##############################################
  
        ilay = 1
  
        if(typelbc(igas).eq.1)then !Fixed density

          Jmat(ilay,ilay,igas) = 1.d0
          Jmat(ilay+1,ilay,igas) = 0.d0
          b(ilay,igas) = valuelbc(igas)

          !print*,igas,typelbc(igas),valuelbc(igas)
  
        elseif(typelbc(igas).eq.2)then   !Fixed flux

          Jmat(ilay,ilay,igas) = 1.0 + deltat/2.d0 * klsi(ilay,igas)
          Jmat(ilay,ilay+1,igas) = -deltat/2.d0 * ksi(ilay,igas)

          b(ilay,igas) = N(ilay,igas) * (1.0-deltat/2.d0*klsi(ilay,igas)) + &
                         N(ilay+1,igas) * (deltat/2.d0*ksi(ilay,igas)) + &
                         deltat/delz*valuelbc(igas) + physrates(ilay,igas) * deltat
  
        endif
  
  
        !Upper boundary
        !##############################################
  
        ilay = nlay
  
        if(typeubc(igas).eq.1)then !Fixed density
  
          Jmat(ilay,ilay,igas) = 1.d0

          b(ilay,igas) = valueubc(igas)
  
        elseif(typeubc(igas).eq.2)then   !Fixed flux
  
          Jmat(ilay,ilay,igas) = 1.d0 + deltat/2.d0 * ksim1(ilay,igas)
          Jmat(ilay,ilay-1,igas) = -deltat/2.d0*klsim1(ilay,igas)

          b(ilay,igas) = N(ilay,igas) * ( 1.d0 - deltat/2.d0*ksim1(ilay,igas) ) + &
                         N(ilay-1,igas) * (deltat/2.*klsim1(ilay,igas)) - &
                         deltat/delz*valueubc(igas) + physrates(ilay,igas) * deltat
  
        elseif(typeubc(igas).eq.3)then   !Fixed velocity
  
          Jmat(ilay,ilay,igas) = (1.d0 + deltat/2.d0/delz *valueubc(igas) + deltat/2.d0 * ksim1(ilay,igas))
          Jmat(ilay,ilay-1,igas) = (-deltat/2.d0*klsim1(ilay,igas))

          b(ilay,igas) = N(ilay,igas) * ( 1.d0 - deltat/2.d0/delz * valueubc(igas) - deltat/2.d0*ksim1(ilay,igas) ) + &
                         N(ilay-1,igas) * (deltat/2.*klsim1(ilay,igas)) + &
                         physrates(ilay,igas) * deltat
          
        endif

  
      enddo
      !!$OMP END DO 
      !!$OMP END PARALLEL
  
      !Re-computing the Jacobian matrix if some species is fixed
      do igas=1,ngas
        do ilay=1,nlay
          if(fix_species(ilay,igas).eq.1)then
            !do jlay=1,nlay
            !  Jmat(ilay,jlay,igas) = 1.d0
            !  b(ilay,igas) = N(ilay,igas)
            !enddo
            Jmat(ilay,:,igas) = 0.d0
            Jmat(ilay,ilay,igas) = 1.d0
            b(ilay,igas) = N(ilay,igas)
          endif
        enddo
      enddo
  
  end subroutine

!==========================================================================================================
  
  subroutine calc_jacobian_diffusion(nlay,ngas,ksi,klsi,ksim1,klsim1,typelbc,valuelbc,typeubc,valueubc,fix_species,Jmat)

    !Function to calculate the jacobian matrix for diffusion

    use omp_lib


    !Inputs
    double precision, intent(in) :: ksi(nlay,ngas)
    double precision, intent(in) :: klsi(nlay,ngas)
    double precision, intent(in) :: ksim1(nlay,ngas)
    double precision, intent(in) :: klsim1(nlay,ngas)
    double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)   !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
    integer, intent(in) :: typelbc(ngas),typeubc(ngas)  !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
    integer, intent(in) :: nlay,ngas
    integer, intent(in) :: fix_species(nlay,ngas)       !Flag to indicate if density must be fixed (if 1)

    !Local 
    integer :: igas,ilay

    !Output
    double precision, intent(out) :: Jmat(nlay,nlay,ngas)

    

    !Filling the Jacobian matrix
    !!$OMP PARALLEL PRIVATE(igas,ilay,jlay) SHARED(Jmat,typelbc,ksi,klsi,ksim1,klsim1,typeubc) NUM_THREADS(8)
    !!$OMP DO 
    do igas=1,ngas

      !Initialising the matrix
      do ilay=1,nlay
        do jlay=1,nlay
          Jmat(ilay,jlay,igas) = 0.d0
        enddo
      enddo


      !Lower boundary
      !##############################################

      ilay = 1

      if(typelbc(igas).eq.1)then !Fixed density

        Jmat(ilay,ilay,igas) = 0.d0

      elseif(typelbc(igas).eq.2)then   !Fixed flux

        Jmat(ilay,ilay,igas) = -klsi(ilay,igas)
        Jmat(ilay,ilay+1,igas) = ksi(ilay,igas)

      endif


      !Upper boundary
      !##############################################

      ilay = nlay

      if(typeubc(igas).eq.1)then !Fixed density

        Jmat(ilay,ilay,igas) = 0.d0

      elseif(typeubc(igas).eq.2)then   !Fixed flux

        Jmat(ilay,ilay,igas) = - ksim1(ilay,igas)
        Jmat(ilay,ilay-1,igas) = klsim1(ilay,igas)

      elseif(typeubc(igas).eq.3)then   !Fixed velocity

        Jmat(ilay,ilay,igas) = (- ksim1(ilay,igas))
        Jmat(ilay,ilay-1,igas) = (klsim1(ilay,igas))
        
      endif

      !Inbetween layers
      !##############################################

      do ilay=2,nlay-1

        Jmat(ilay,ilay,igas) = - (klsi(ilay,igas)+ksim1(ilay,igas))
        Jmat(ilay,ilay-1,igas) = klsim1(ilay,igas)
        Jmat(ilay,ilay+1,igas) = ksi(ilay,igas)

      enddo

    enddo
    !!$OMP END DO 
    !!$OMP END PARALLEL

    !Re-computing the Jacobian matrix if some species is fixed
    do igas=1,ngas
      do ilay=1,nlay
        if(fix_species(ilay,igas).eq.1)then
          Jmat(ilay,:,igas) = 0.d0
        endif
      enddo
    enddo

  end subroutine
  
!==========================================================================================================
  
  subroutine eval_fn_diffusion(nlay,ngas,delz,N,Jmat,physrates,typelbc,valuelbc,typeubc,valueubc,fix_species,b)
  
    !Function to evaluate the diffusion function (-dphi/dz + P - L) for a given density state
    !
    !This is used for solving the continuity equation, in which we want to solve:
    !(I+deltat*gamma*J) * N_t+1 = fn(N_t)
  
  
    !Inputs
    real, intent(in) :: delz                              !Layer width (m) assumed to be constant
    double precision, intent(in) :: N(nlay,ngas)          !Number density of each species (m-3)
    double precision, intent(in) :: Jmat(nlay,nlay,ngas)  !Jacobian matrix
    double precision, intent(in) :: physrates(nlay,ngas)              !Production/loss rates from other processes (such as chemistry)
    double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)     !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
    integer, intent(in) :: typelbc(ngas),typeubc(ngas)    !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
    integer, intent(in) :: nlay,ngas
    integer, intent(in) :: fix_species(nlay,ngas)         !Flag to indicate if density must be fixed (if 1)

    !Local 
    integer :: igas,ilay
    double precision :: FNi(nlay,1),Ji(nlay,nlay),Ni(nlay,1),FN(nlay,ngas)

    !Output
    double precision, intent(out) :: b(nlay,ngas)


    !First we multiply the Jacobian matrix times N
    do igas=1,ngas
      do ilay=1,nlay
        do jlay=1,nlay
          Ji(ilay,jlay) = Jmat(ilay,jlay,igas)
        enddo
        Ni(ilay,1) = N(ilay,igas)
      enddo

      FNi = MATMUL(Ji,Ni)
      !call mmul(nlay,nlay,1,Ji,Ni,Fni)

      do ilay=1,nlay
        FN(ilay,igas) = FNi(ilay,1)
      enddo
    enddo

    !Now we include the boundary conditions and the rates from other processes
    !!$OMP DO 
    do igas=1,ngas  

      !Initialising array
      do ilay=1,nlay
        b(ilay,igas) = 0.d0
      enddo

      !Lower boundary
      !##############################################

      ilay = 1

      if(typelbc(igas).eq.1)then !Fixed density

        b(ilay,igas) = FN(ilay,igas) + valuelbc(igas)

      elseif(typelbc(igas).eq.2)then   !Fixed flux

        b(ilay,igas) =  FN(ilay,igas) + valuelbc(igas)/delz + physrates(ilay,igas)

      endif


      !Upper boundary
      !##############################################

      ilay = nlay


      if(typeubc(igas).eq.1)then !Fixed density

        b(ilay,igas) = FN(ilay,igas) + valueubc(igas)

      elseif(typeubc(igas).eq.2)then   !Fixed flux

        b(ilay,igas) = FN(ilay,igas) - valueubc(igas)/delz + physrates(ilay,igas)

      elseif(typeubc(igas).eq.3)then   !Fixed velocity

        b(ilay,igas) = FN(ilay,igas) - N(ilay,igas) * valueubc(igas)/delz + physrates(ilay,igas)
        
      endif


      !Inbetween layers
      !##############################################

      do ilay=2,nlay-1

        b(ilay,igas) = FN(ilay,igas) + physrates(ilay,igas)

      enddo

    enddo
    !!$OMP END DO

    !Re-computing the function if some species is fixed
    do igas=1,ngas
      do ilay=1,nlay
        if(fix_species(ilay,igas).eq.1)then
          b(ilay,igas) = FN(ilay,igas)
        endif
      enddo
    enddo
    
  end subroutine
  
!==========================================================================================================
  
  subroutine calc_diffusion_coeff(nlay,ngas,h,T,scaleH0,scaleH,K,D,B,ksi,klsi,ksim1,klsim1)
  
    !Function to calculate the diffusion coefficients that fill the Jacobian matrix
  
  
    !Inputs
    real, intent(in) :: h(nlay)            !Altitude (m)
    real, intent(in) :: T(nlay)            !Temperature (K)
    real, intent(in) :: scaleH0(nlay)      !Mean scale height (m)
    real, intent(in) :: scaleH(nlay,ngas)  !Scale height of each species (m)
    real, intent(in) :: K(nlay)            !Eddy diffusion coefficient (m2 s-1)
    real, intent(in) :: D(nlay,ngas)       !Molecular diffusion coefficient (m2 s-1)
    real, intent(in) :: B(ngas)            !Molecular thermal diffusion coefficient of each gas
    integer, intent(in) :: nlay,ngas


    !Local
    integer :: ilay
    real :: K_i,T_i,H0_i,H_i
    real :: delz
  
    !Output
    ! Coefficients k_i*, kl_i*, k_(i-1)*, kl_(i-1)* in each layer (s-1)
    double precision, intent(out) :: ksi(nlay,ngas)
    double precision, intent(out) :: klsi(nlay,ngas)
    double precision, intent(out) :: ksim1(nlay,ngas)
    double precision, intent(out) :: klsim1(nlay,ngas)
  
  
    delz = h(2) - h(1)   !Layer width (m) assumed to be constant

    do ilay=1,nlay
      
      !Coefficients between layers i and i-1
      if(ilay.ne.1)then

        K_i = (K(ilay)+K(ilay-1))/2.
        T_i = (T(ilay)+T(ilay-1))/2.
        H0_i = (scaleH0(ilay)+scaleH0(ilay-1))/2.

        do igas=1,ngas

          D_i = (D(ilay,igas)+D(ilay-1,igas))/2.
          H_i = (scaleH(ilay,igas)+scaleH(ilay-1,igas))/2.

          ksim1(ilay,igas) = (D_i+K_i)/delz/delz 
          klsim1(ilay,igas) = (K_i/delz/delz) * (1.0 - delz/H0_i - (T(ilay)-T(ilay-1))/T_i ) + &
                              (D_i/delz/delz) * (1.0 - delz/H_i - (1.0+B(igas))*(T(ilay)-T(ilay-1))/T_i )

        enddo

      endif

      !Coefficients between layers i and i+1
      if(ilay.ne.nlay)then

        K_i = (K(ilay)+K(ilay+1))/2.
        T_i = (T(ilay)+T(ilay+1))/2.
        H0_i = (scaleH0(ilay)+scaleH0(ilay+1))/2.

        do igas=1,ngas

          D_i = (D(ilay,igas)+D(ilay+1,igas))/2.
          H_i = (scaleH(ilay,igas)+scaleH(ilay+1,igas))/2.

          ksi(ilay,igas) = (D_i+K_i)/delz/delz
          klsi(ilay,igas) = (K_i/delz/delz) * (1.0 - delz/H0_i - (T(ilay+1)-T(ilay))/T_i ) + &
                            (D_i/delz/delz) * (1.0 - delz/H_i - (1.0+B(igas))*(T(ilay+1)-T(ilay))/T_i )

        enddo

      endif

    end do
      
  end subroutine
  
!==========================================================================================================
  
  subroutine calc_Keddy(nlay,h,N0,K)
  
    !Routine to calculate the Eddy diffusion coefficient at each altitude level in (m2 s-1)
      
    !We assume :
      
    !      K = 1.0e6                          (z<=60 km)     [cm2 s-1]
    !      K = 2.013 * (1./n(z))**0.5         (z>60km)       [cm2 s-1]
      
    !Inputs
    real, intent(in) :: h(nlay)               !Altitude (m)
    double precision, intent(in) :: N0(nlay)  !Number density (m-3)
    integer, intent(in) :: nlay
  
    !Local
    integer :: ilay
  
    !Output
    real, intent(out) :: K(nlay)  !Eddy diffusion coefficient (m2 s-1)
  
  
    do ilay=1,nlay
        
      if (h(ilay)/1.0e3.LT.60.) then
        K(ilay) = 1.0e6 * 1.0e-4  !m2 s-1
      else
        K(ilay) = 2.0e13 * sqrt(1./(N0(ilay)*1.0e-6)) * 1.0e-4  !m2 s-1
      endif
      
    end do
      
  end subroutine
      
!============================================================================================================
  
  subroutine calc_Dmoldiff(nlay,ngas,N0,T,A,s,D)
  
    ! Routine to calculate the molecular diffusion coefficients for each species at each level
  
    ! The main citation for this is Hunten (1973). It is calculated as:
  
    ! D_i = A_i * temp**s_i / numdens
    
    ! The diffusion coefficients of H and H2 in CO2 are provided in this reference.
  
    ! For the rest of the species, we set A and s to 1.0 and 0.75 (Cangi et al. 2020)
        
    ! Inputs
    double precision, intent(in) :: N0(nlay)    !Number density (m-3)
    real, intent(in) :: T(nlay)                 !Temperature (K)
    integer, intent(in) :: nlay,ngas
    real, intent(in) :: A(ngas),s(ngas)         !Coefficients
        
    !Local
    integer :: ilay,igas
        
    !Output
    real, intent(out) :: D(nlay,ngas)  !Eddy diffusion coefficient (m2 s-1)
        
        
    do ilay=1,nlay
      do igas=1,ngas
        D(ilay,igas) = A(igas) * T(ilay)**s(igas) / (N0(ilay)*1.0e-6) * 1.0e-4  !m2 s-1
      enddo
    end do
  
  end subroutine
  
!============================================================================================================
    
  subroutine calc_mmean(nlay,ngas,N,mmol,mmean)
    
    !Routine to calculate the mean molecular weight of the atmosphere
        
    ! Inputs
    double precision, intent(in) :: N(nlay,ngas)   !Number density of each species (m-3)
    real, intent(in) :: mmol(ngas)     !Molecular weight (uma)
    integer, intent(in) :: nlay,ngas
        
    !Local
    integer :: ilay,igas
    real :: sum1,sum2
        
    !Output
    real, intent(out) :: mmean(nlay)   !Mean molecular weight (uma)
        
    do ilay=1,nlay
      sum1 = 0.0
      sum2 = 0.0
      do igas=1,ngas
        sum1 = sum1 + N(ilay,igas)
        sum2 = sum2 + N(ilay,igas) * mmol(igas)
      enddo
      mmean(ilay) = sum2 / sum1
    enddo
    
  end subroutine
    
!============================================================================================================
    
  subroutine calc_scaleH_mean(nlay,h,T,mmean,scaleH0)
  
    ! Routine to calculate the scale height of the atmosphere
      
    include 'const.h'
  
    ! Inputs
    real, intent(in) :: h(nlay)         !Altitude (m)
    real, intent(in) :: T(nlay)         !Temperature (K)
    real, intent(in) :: mmean(nlay)     !Molecular weight (uma)
    integer, intent(in) :: nlay
        
    !Local
    integer :: ilay
    real :: grav(nlay)
        
    !Output
    real, intent(out) :: scaleH0(nlay)   !Scale height (m)
        
    do ilay=1,nlay
      grav(ilay) = Mplanet * G / (h(ilay)+Rplanet)**2.0
      scaleH0(ilay) = k_B * T(ilay) / ( (mmean(ilay)/N_A/1.0e3) * grav(ilay) )
    enddo
  
  end subroutine
    
!============================================================================================================
    
  subroutine calc_scaleH(nlay,ngas,h,T,mmol,scaleH)
  
    ! Routine to calculate the scale height of the atmosphere for each species
      
    include 'const.h'
  
    ! Inputs
    integer, intent(in) :: nlay,ngas
    real, intent(in) :: h(nlay)         !Altitude (m)
    real, intent(in) :: T(nlay)         !Temperature (K)
    real, intent(in) :: mmol(ngas)      !Molecular weight (uma)
  
        
    !Local
    integer :: ilay,igas
    real :: grav(nlay)
        
    !Output
    real, intent(out) :: scaleH(nlay,ngas)   !Scale height (m)
        
    do ilay=1,nlay
      do igas=1,ngas
        grav(ilay) = Mplanet * G / (h(ilay)+Rplanet)**2.0
        scaleH(ilay,igas) = k_B * T(ilay) / ( (mmol(igas)/N_A/1.0e3) * grav(ilay) )
      enddo
    enddo
  
  end subroutine

!============================================================================================================
    
  subroutine calc_tscale_diffusion(nlay,ngas,h,T,N,mmol,A,s,B,&
    typelbc,valuelbc,typeubc,valueubc,fix_species,t_diff)
  
    ! Routine to calculate the characteristic timescale for diffusion using
    !
    ! t_diff = N / (dphi/dz)
    !
    ! from Vuitton et al. (2019)

      
    include 'const.h'
  
    ! Inputs
    real, intent(in) :: h(nlay)                           !Altitude (m)
    real, intent(in) :: T(nlay)                           !Temperature (K)
    double precision, intent(in) :: N(nlay,ngas)          !Number density of each species (m-3)
    real, intent(in) :: mmol(ngas)                        !Molecular weight (uma)
    real, intent(in) :: A(ngas),s(ngas),B(ngas)           !Coefficients for molecular diffusion
    double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)     !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
    integer, intent(in) :: typelbc(ngas),typeubc(ngas)    !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
    integer, intent(in) :: nlay,ngas
    integer, intent(in) :: fix_species(nlay,ngas)         !Flag to indicate if density must be fixed (if 1)
    
    
    !Local
    double precision :: N0(nlay)
    real :: K(nlay)                             !Eddy diffusion coefficient (m2 s-1)
    real :: D(nlay,ngas)                        !Molecular diffuson coefficient (m2 s-1)
    real :: mmean(nlay)                         !Mean molecular weight (uma)
    real :: scaleH0(nlay),scaleH(nlay,ngas)     !Scale height (m)
    double precision :: dphi(nlay,ngas)
    double precision :: ksi(nlay,ngas),klsi(nlay,ngas),ksim1(nlay,ngas),klsim1(nlay,ngas)
  

    !Output
    double precision, intent(out) :: t_diff(nlay,ngas)    !Diffusion timescale (s)
        
    !Calculating the total density (m-3)
    do ilay=1,nlay
      N0(ilay) = 0.d0
      do igas=1,ngas
        N0(ilay) = N0(ilay) + N(ilay,igas)
      enddo
    enddo

    delz = h(2) - h(1)  !W

    !Calculating the Eddy diffusion coefficient in each layer
    call calc_Keddy(nlay,h,N0,K)

    !Calculating the molecular diffusion coefficient in each layer
    call calc_Dmoldiff(nlay,ngas,N0,T,A,s,D)

    !Calculating the mean molecular weight
    call calc_mmean(nlay,ngas,N,mmol,mmean)

    !Calculating the scale height
    call calc_scaleH_mean(nlay,h,T,mmean,scaleH0)

    !Calculating the species-dependent scale height
    call calc_scaleH(nlay,ngas,h,T,mmol,scaleH)

    !Calculating the diffusion coefficients to fill the jacobian 
    call calc_diffusion_coeff(nlay,ngas,h,T,scaleH0,scaleH,K,D,B,ksi,klsi,ksim1,klsim1)

    !Calculating the divergence of the flux dphi/dz
    do igas=1,ngas
  
      !Initialising the matrix
      do ilay=1,nlay
        dphi(ilay,igas) = 0.d0
      enddo


      !Lower boundary
      !##############################################

      ilay = 1

      if(typelbc(igas).eq.1)then !Fixed density

        dphi(ilay,igas) = 0.d0

      elseif(typelbc(igas).eq.2)then   !Fixed flux

        dphi(ilay,igas) = N(ilay,igas) * klsi(ilay,igas) - &
                       N(ilay+1,igas) * ksi(ilay,igas) - &
                       valuelbc(igas)/delz

      endif


      !Upper boundary
      !##############################################

      ilay = nlay

      if(typeubc(igas).eq.1)then !Fixed density

        dphi(ilay,igas) = 0.d0

      elseif(typeubc(igas).eq.2)then   !Fixed flux

        dphi(ilay,igas) = N(ilay,igas) * ksim1(ilay,igas) - &
                       N(ilay-1,igas) * klsim1(ilay,igas) + &
                       valueubc(igas)/delz

      elseif(typeubc(igas).eq.3)then   !Fixed velocity

        dphi(ilay,igas) = N(ilay,igas) * ksim1(ilay,igas) - &
                          N(ilay-1,igas) * klsim1(ilay,igas) + &
                          valueubc(igas)/delz*N(ilay,igas)
        
      endif

      !Inbetween layers
      !##############################################

      do ilay=2,nlay-1

        dphi(ilay,igas) = N(ilay,igas) * (klsi(ilay,igas)+ksim1(ilay,igas)) - &
                       N(ilay-1,igas) * klsim1(ilay,igas) - &
                       N(ilay+1,igas) * ksi(ilay,igas) 

      enddo

    enddo

    !Setting the timescale to 0 if a species is fixed
    do igas=1,ngas
      do ilay=1,nlay
        if(fix_species(ilay,igas).eq.1)then
          do jlay=1,nlay
            dphi(ilay,igas) = 0.d0
          enddo
        endif
      enddo
    enddo

    !Calculating the characteristic diffusion timescale
    do igas=1,ngas
      do ilay=1,nlay
        t_diff(ilay,igas) = abs( N(ilay,igas) / dphi(ilay,igas) )
      enddo
    enddo

  end subroutine


!============================================================================================================
  
  subroutine mmul(IDM1,IDM2,IDM3,AM1,AM2,ANS)

    !Routine to perform a matrix multiplication

    IMPLICIT DOUBLE PRECISION (A-H,O-Z)

    integer,intent(in) :: IDM1,IDM2,IDM3
    double precision,intent(in) :: AM1(IDM1,IDM2),AM2(IDM2,IDM3)

    integer :: i,j,k
    double precision :: AIJ

    double precision,intent(out) :: ANS(IDM1,IDM3)

    do j=1,IDM3
      do i=1,IDM1
        AIJ = 0.0D0
        do K=1,IDM2
          AIJ = AIJ + AM1(I,K)*AM2(K,J)
        enddo
        ANS(i,j)=AIJ
      enddo
    enddo

  end subroutine
    

END MODULE diffusion