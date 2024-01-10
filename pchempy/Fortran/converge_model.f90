!==============================================================================================================================================

subroutine converge_model_mars_rosenbrock(nlay,ngas,h,P,T,N,gasID,isoID,mmol,A,s,B,&
    typelbc,valuelbc,typeubc,valueubc,fix_species,sza,tauvis,dist_sun,dtmin,dtmax,&
    chemistry_inc,nitrogen_chemistry,c13_chemistry,o18_chemistry,Nnew,time,dtcur)

    !  Function to calculate the density at the next time step using a second-order Rosenbrock solver
    !
    !  Using the method from Verwer et al. (1999) and Hobbs et al. (2019)

    use diffusion
    use chemistry
    use photolysis
    use mars_chemistry

    implicit none

    !Inputs
    real, intent(in) :: h(nlay)                           !Altitude (m)
    double precision, intent(in) :: P(nlay)               !Pressure (Pa)
    real, intent(in) :: T(nlay)                           !Temperature (K)
    double precision, intent(in) :: N(nlay,ngas)          !Number density of each species (m-3)
    real, intent(in) :: mmol(ngas)                        !Molecular weight (uma)
    real, intent(in) :: A(ngas),s(ngas),B(ngas)           !Coefficients for molecular diffusion
    double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)     !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
    integer, intent(in) :: typelbc(ngas),typeubc(ngas)    !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
    integer, intent(in) :: gasID(ngas),isoID(ngas)        !Gas and Isotope IDs
    integer, intent(in) :: nlay,ngas
    integer, intent(in) :: fix_species(nlay,ngas)         !Flag to indicate if density must be fixed (if 1)
    real, intent(in) :: sza,tauvis,dist_sun               !Solar zenith angle (degrees), dust opacity, Sun-planet distance (AU)
    double precision, intent(in) :: dtmin                 !Timestep to start the simulation (s)
    double precision, intent(in) :: dtmax                 !Timestep at which to convert the simulation (s)
    logical, intent(in) :: chemistry_inc                  !Flag to include chemistry
    logical, intent(in) :: nitrogen_chemistry             !Flag to include the nitrogen chemistry
    logical, intent(in) :: c13_chemistry                  !Flag to include the (13C) chemistry
    logical, intent(in) :: o18_chemistry                  !Flag to include the (18O) chemistry

    !Local 
    double precision :: c(nlay,ngas),err(nlay,ngas),Jmat(ngas,ngas),N0(nlay)
    double precision :: M_chem(nlay,ngas,ngas),M_diff(nlay,nlay,ngas)
    double precision :: A_tri(nlay,ngas,ngas),B_tri(nlay,ngas,ngas),C_tri(nlay,ngas,ngas)
    double precision :: A1_tri(nlay,ngas,ngas),B1_tri(nlay,ngas,ngas),C1_tri(nlay,ngas,ngas)
    double precision :: R(nlay,ngas),X(nlay,ngas),g1(nlay,ngas),g2(nlay,ngas),b1(nlay,ngas)
    integer :: ilay,igas,jgas,nreactions,n_photolysis
    integer :: ngas_phot, gasID_phot(ngas), isoID_phot(ngas)
    double precision :: Nold(nlay,ngas),Ncur(nlay,ngas)
    double precision :: ccur(nlay,ngas),cnew(nlay,ngas)
    double precision :: cdens(nlay,ngas)
    real :: K(nlay)                             !Eddy diffusion coefficient (m2 s-1)
    real :: D(nlay,ngas)                        !Molecular diffuson coefficient (m2 s-1)
    real :: mmean(nlay)                         !Mean molecular weight (uma)
    real :: scaleH0(nlay),scaleH(nlay,ngas)     !Scale height (m)
    double precision :: ksi(nlay,ngas),klsi(nlay,ngas),ksim1(nlay,ngas),klsim1(nlay,ngas)
    double precision :: gamma
    integer :: iter,iterc
    real :: delz

    !all reactions
    double precision, allocatable :: rrates(:,:)
    integer, allocatable :: rtype(:)
    integer, allocatable :: ns(:),sID(:,:),sISO(:,:)
    integer, allocatable :: npr(:),pID(:,:),pISO(:,:)
    real, allocatable :: sf(:,:),pf(:,:)
    integer, allocatable :: sID_pos(:,:),pID_pos(:,:)

    !photolysis reactions
    double precision, allocatable :: rratesp(:,:)
    integer, allocatable :: rtypep(:)
    integer, allocatable :: nsp(:),sIDp(:,:),sISOp(:,:)
    integer, allocatable :: nprp(:),pIDp(:,:),pISOp(:,:)
    real, allocatable :: sfp(:,:),pfp(:,:)

    !convergence
    logical :: converged
    real :: rtol = 1.0e-4 !relative tolerance
    real :: atol = 0.05   !absolute tolerance
    integer :: max_iter = 5000
    double precision :: mindt = 1.0d-4

    double precision :: dtnew
    real :: lerr !local error
    real :: coef,e , egas(ngas)

    !Flags
    logical :: photolysis_inc


    !Output
    double precision, intent(out) :: Nnew(nlay,ngas)
    double precision, intent(out) :: time,dtcur



    !INITIALISATIONS
    !################################################################################

    !Defining if photolysis needs to be computed
    if(sza.le.135.0)then
        photolysis_inc = .true.
    endif

    !Calculating the number of reactions
    call number_reactions(nitrogen_chemistry,c13_chemistry,o18_chemistry,photolysis_inc,&
    nreactions,n_photolysis)

    !Allocating arrays
    allocate(rrates(nlay,nreactions))
    allocate(rtype(nreactions))
    allocate(ns(nreactions),sID(2,nreactions),sISO(2,nreactions),sf(2,nreactions))
    allocate(npr(nreactions),pID(2,nreactions),pISO(2,nreactions),pf(2,nreactions))
    allocate(sID_pos(2,nreactions),pID_pos(2,nreactions))

    !Loading the chemical reactions and calculating the reaction rates
    call load_reactions(nlay,ngas,gasID,isoID,h,P,T,N,nitrogen_chemistry,c13_chemistry,o18_chemistry,photolysis_inc,&
    nreactions,n_photolysis,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf)

    !Computing the photolysis rates for a first iteration to define reactions
    if(photolysis_inc)then

        !Allocating arrays
        allocate(rratesp(nlay,n_photolysis))
        allocate(rtypep(n_photolysis))
        allocate(nsp(n_photolysis),sIDp(2,n_photolysis),sISOp(2,n_photolysis),sfp(2,n_photolysis))
        allocate(nprp(n_photolysis),pIDp(2,n_photolysis),pISOp(2,n_photolysis),pfp(2,n_photolysis))


        !loading the gases to be photolysed
        call load_reactions_photolysis(ngas,nitrogen_chemistry,c13_chemistry,o18_chemistry,&
                                        ngas_phot,gasID_phot,isoID_phot)

        !calculating the column density of each species in each layer
        !we assume the density is constant in each layer
        delz = h(2) - h(1)   !Assuming all layers have the same height
        cdens(:,:) = 0.d0
        do ilay=1,nlay
            do igas=1,ngas
                cdens(ilay,igas) = N(ilay,igas) * delz
            enddo
        enddo

        !calculating the photolysis rates 
        call photolysis_online(nlay, ngas, ngas_phot, gasID, isoID, gasID_phot, isoID_phot, &
                                h, T, cdens, sza, tauvis, dist_sun, n_photolysis, &
                                rtypep, nsp, sIDp, sISOp, sfp, nprp, pIDp, pISOp, pfp, rratesp)

        !mapping the photolysis reactions into the overall reaction array
        rtype(1:n_photolysis) = rtypep(1:n_photolysis)
        ns(1:n_photolysis) = nsp(1:n_photolysis)
        sID(:,1:n_photolysis) = sIDp(:,1:n_photolysis)
        sISO(:,1:n_photolysis) = sISOp(:,1:n_photolysis)
        sf(:,1:n_photolysis) = sfp(:,1:n_photolysis)
        npr(1:n_photolysis) = nprp(1:n_photolysis)
        pID(:,1:n_photolysis) = pIDp(:,1:n_photolysis)
        pISO(:,1:n_photolysis) = pISOp(:,1:n_photolysis)
        pf(:,1:n_photolysis) = pfp(:,1:n_photolysis)
        rrates(:,1:n_photolysis) = rratesp(:,1:n_photolysis)

    endif

    !Locating the sources/products in atmospheric gas array
    call locate_gas_reactions(ngas,gasID,isoID,nreactions, ns, sID, sISO, npr, pID, pISO,&
    sID_pos,pID_pos)



    !CONVERGENCE LOOP
    !################################################################################

    !Initialising parameters
    Nold(:,:) = N(:,:)
    Ncur(:,:) = N(:,:)
    dtcur = dtmin
    time = 0.d0
    e = 0.0
    converged = .false.

    iter = 1
    iterc = 1
    do while(converged.eqv..false.)

        !Making sure the timestep is not higher than dt
        if(dtcur.ge.dtmax)then
            converged = .true.
            dtcur = dtmax
        endif

        !Printing values just once every 100 iterations
        if(iterc.eq.1)then
            print*,'Iteration',iter
            print*,'dt',dtcur
            print*,'e',e
            do igas=1,ngas
                print*,'gas',gasID(igas),isoID(igas),'err',egas(igas)
            enddo
        endif

        !CHEMISTRY
        !############################################################################

        if(chemistry_inc)then


            !Loading the chemical reactions and calculating the reaction rates
            call load_reactions(nlay,ngas,gasID,isoID,h,P,T,Ncur,nitrogen_chemistry,c13_chemistry,o18_chemistry,photolysis_inc,&
            nreactions,n_photolysis,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf)

            !Computing the photolysis rates
            if(photolysis_inc)then

                !loading the gases to be photolysed
                call load_reactions_photolysis(ngas,nitrogen_chemistry,c13_chemistry,o18_chemistry,&
                                                ngas_phot,gasID_phot,isoID_phot)

                !calculating the column density of each species in each layer
                !we assume the density is constant in each layer
                delz = h(2) - h(1)   !Assuming all layers have the same height
                cdens(:,:) = 0.d0
                do ilay=1,nlay
                    do igas=1,ngas
                        cdens(ilay,igas) = Ncur(ilay,igas) * delz
                    enddo
                enddo

                !calculating the photolysis rates 
                call photolysis_online(nlay, ngas, ngas_phot, gasID, isoID, gasID_phot, isoID_phot, &
                                        h, T, cdens, sza, tauvis, dist_sun, n_photolysis, &
                                        rtypep, nsp, sIDp, sISOp, sfp, nprp, pIDp, pISOp, pfp, rratesp)


                !mapping the photolysis reactions into the overall reaction array
                rtype(1:n_photolysis) = rtypep(1:n_photolysis)
                ns(1:n_photolysis) = nsp(1:n_photolysis)
                sID(:,1:n_photolysis) = sIDp(:,1:n_photolysis)
                sISO(:,1:n_photolysis) = sISOp(:,1:n_photolysis)
                sf(:,1:n_photolysis) = sfp(:,1:n_photolysis)
                npr(1:n_photolysis) = nprp(1:n_photolysis)
                pID(:,1:n_photolysis) = pIDp(:,1:n_photolysis)
                pISO(:,1:n_photolysis) = pISOp(:,1:n_photolysis)
                pf(:,1:n_photolysis) = pfp(:,1:n_photolysis)
                !rrates(:,1:n_photolysis) = rratesp(:,1:n_photolysis) / 2.!we divide by factor of two to account for night
                rrates(:,1:n_photolysis) = rratesp(:,1:n_photolysis)  !Mean dayside conditions

            endif

            !Calculating the Jacobian for the chemistry
            M_chem(:,:,:) = 0.d0
            ccur(:,:) = Ncur(:,:) * 1.0d-6

            do ilay=1,nlay
                call calc_jacobian_chemistry(nlay, ngas, ilay, ccur,  &
                nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, Jmat)
                M_chem(ilay,:,:) = Jmat(:,:)
            enddo

            !Fixing the Jacobian for specified species
            do ilay=1,nlay
                do igas=1,ngas
                    if(fix_species(ilay,igas).eq.1)then
                        M_chem(ilay,igas,:) = 0.d0
                    endif
                enddo
            enddo

            !Fixing the Jacobian if the boundary conditions specify a fixed density
            do igas=1,ngas
                if(typelbc(igas).eq.1)then
                    M_chem(1,igas,:) = 0.d0
                endif

                if(typeubc(igas).eq.1)then
                    M_chem(nlay,igas,:) = 0.d0
                endif                
            enddo

        else

            M_chem(:,:,:) = 0.d0

        endif

        !DIFFUSION
        !############################################################################

        !Calculating the Jacobian for the diffusion

        !Performing calculations to characterise the atmosphere
        !############################################################

        !Calculating the total density (m-3)
        do ilay=1,nlay
            N0(ilay) = 0.d0
            do igas=1,ngas
                N0(ilay) = N0(ilay) + Ncur(ilay,igas)
            enddo
        enddo
        
        delz = h(2) - h(1)  !W
        
        !Calculating the Eddy diffusion coefficient in each layer
        call calc_Keddy(nlay,h,N0,K)
        
        !Calculating the molecular diffusion coefficient in each layer
        call calc_Dmoldiff(nlay,ngas,N0,T,A,s,D)
        
        !Calculating the mean molecular weight
        call calc_mmean(nlay,ngas,Ncur,mmol,mmean)
    
        !Calculating the scale height
        call calc_scaleH_mean(nlay,h,T,mmean,scaleH0)
    
        !Calculating the species-dependent scale height
        call calc_scaleH(nlay,ngas,h,T,mmol,scaleH)
    
    
        !Calculating the diffusion parameters
        !################################################################
    
        !Calculating the diffusion coefficients to fill the jacobian 
        call calc_diffusion_coeff(nlay,ngas,h,T,scaleH0,scaleH,K,D,B,ksi,klsi,ksim1,klsim1)
    
        !Calculating the Jacobian matrix
        call calc_jacobian_diffusion(nlay,ngas,ksi,klsi,ksim1,klsim1,typelbc,valuelbc,typeubc,valueubc,fix_species,M_diff)


        !SOLVING CHEMISTRY+DIFFUSION SYSTEM
        !############################################################################


        !Calculating the block tridiagonal matrix (filling the A,B and C blocks)
        A_tri(:,:,:) = 0.d0
        B_tri(:,:,:) = 0.d0
        C_tri(:,:,:) = 0.d0

        do ilay=1,nlay

            !Diagonal
            B_tri(ilay,:,:) = M_chem(ilay,:,:)
            do igas=1,ngas
                B_tri(ilay,igas,igas) = B_tri(ilay,igas,igas) + M_diff(ilay,ilay,igas)
            enddo

            !Subdiagonal
            if(ilay.gt.1)then
                do igas=1,ngas
                    A_tri(ilay,igas,igas) = A_tri(ilay,igas,igas) + M_diff(ilay,ilay-1,igas)
                enddo
            endif

            !Superdiagonal
            if(ilay.lt.nlay)then
                do igas=1,ngas
                    C_tri(ilay,igas,igas) = C_tri(ilay,igas,igas) + M_diff(ilay,ilay+1,igas)
                enddo
            endif

        enddo

        

        !Solving the system using a second-order Rosenbrock solver
        !==============================================================

        ccur(:,:) = Ncur(:,:) * 1.0d-6
        gamma = 1.d0 + 1.d0/sqrt(2.d0)

        !Calculating the matrix (I - gamma*dt*J) and right hand side at f(N)

        do ilay=1,nlay
            do igas=1,ngas
                do jgas=1,ngas
                    if(igas.eq.jgas)then
                        B1_tri(ilay,igas,jgas) = 1.d0 - dtcur * gamma * B_tri(ilay,igas,jgas)
                    else
                        B1_tri(ilay,igas,jgas) = - dtcur * gamma * B_tri(ilay,igas,jgas)
                    endif
                enddo
            enddo

            A1_tri(ilay,:,:) = -dtcur * gamma * A_tri(ilay,:,:)
            C1_tri(ilay,:,:) = -dtcur * gamma * C_tri(ilay,:,:)

        enddo

        call eval_fn_system(ngas,nlay,A_tri,B_tri,C_tri,ccur,typelbc,valuelbc,typeubc,valueubc,fix_species,delz,R)


        !Solving for g1
        call blktri(ngas,nlay,A1_tri,B1_tri,C1_tri,R,g1)

        !Calculating b1
        call eval_fn_system(ngas,nlay,A_tri,B_tri,C_tri,(ccur+dtcur*g1),typelbc,valuelbc,typeubc,valueubc,fix_species,delz,R)
        R(:,:) = R(:,:) - 2.d0*g1(:,:)

        !Solving for g2
        call blktri(ngas,nlay,A1_tri,B1_tri,C1_tri,R,g2)

        !Calculating new density
        cnew(:,:) = ccur(:,:) + 1.5d0 * dtcur * g1(:,:) + 0.5d0 * dtcur * g2(:,:)
        Nnew(:,:) = cnew(:,:) * 1.0d6
        !cnew(:,:) = ccur(:,:) + 1.d0 * dtcur * g1(:,:)
        !Nnew(:,:) = cnew(:,:) * 1.0d6


        !CONVERGENCE ASSESMENT
        !############################################################################

        !Defining the local error
        !============================

        err(:,:) = abs(cnew(:,:) - (ccur(:,:) + dtcur * g1(:,:))) * 1.0d6   !m-3

        e = 0.0
        egas(:) = 0.0
        do ilay=1,nlay
            do igas=1,ngas
                lerr = dabs( (err(ilay,igas)) &
                            / (Ncur(ilay,igas)*rtol+atol) )
                if(lerr.gt.e)then
                    e = lerr
                endif
                if(lerr.gt.egas(igas))then
                    egas(igas) = lerr
                endif
            enddo
        enddo

        !do igas=1,ngas
        !    print*,egas(igas)
        !enddo

        if(e.le.0.0)e=0.1


        !Computing next timestep
        !========================

        !Timestep correction
        coef = max(0.1, min(2.0,1.0/sqrt(e)))
        !dtnew = max(dtcur,dtcur*coef)
        dtnew = max(mindt,dtcur*coef)

        !print*,dtcur,e,dtnew
        !pause


        !Updating densities and times if error is lower than limit
        !============================================================

        !print*,dtnew

        if(e.le.1.1)then

        !    Nold(:,:) = Ncur(:,:)
            Ncur(:,:) = Nnew(:,:)

            time = time + dtcur
        !    dtold = dtcur
            dtcur = dtnew

        else

            if(dtcur.eq.mindt)then
                Ncur(:,:) = Nnew(:,:)
                time = time + dtcur
            endif

            dtcur = dtnew

        endif


        if(iter.eq.max_iter)then
            converged = .true.
        endif

        iter = iter + 1
        iterc = iterc + 1

        if(iterc.eq.101)iterc=1


    enddo

end subroutine

!==============================================================================================================================================

subroutine blktri(M,N,A,B,C,R,X)
    !     .......THIS PROGRAM SOLVES A TRI-BLOCK-DIAGONAL MATRIX.
    !      .......NDIM = DIMENSION OF THE A,B AND C IN THE MAIN PROGRAM
    !      .......MDIM = DIMENSION OF BLOCKS IN MAIN PROGRAM
    !      .......M = SIZE OF BLOCKS
    !      .......N = NUMBER OF BLOCKS
    !      .......A,B,C = M*M MATRICES; A GOES FROM 2 TO N, B GOES
    !      .......        FROM 1 TO N AND C GOES FROM 1 TO N-1.
    !      .......X = RESULT VECTOR
    !      .......R = RIGHT HAND SIDE VECTOR
    !     C

    !Input
    integer, intent(in) :: M, N
    double precision, intent(inout) :: A(N,M,M),B(N,M,M),C(N,M,M)
    double precision, intent(inout) :: R(N,M)

    !Local
    integer :: MM1,JP1,I,J,MP
    double precision :: T,S

    !Output
    double precision, intent(out) :: X(N,M)

    !.....FORWARD ELIMINATION-BLOCKS
    do I=1,N

        !.....GAUSS JORDAN FOR B(I)
        MM1 = M - 1

        do J=1,MM1

            !.....NORMALIZE PIVOT ROW
            JP1 = J+1
            do K=JP1,M
                T=1.d0/B(I,J,J)
                B(I,J,K) = T*B(I,J,K)
            enddo

            do K=1,M
                if(I.NE.M)then
                    C(I,J,K) = T*C(I,J,K)
                endif
            enddo

            R(I,J) = T*R(I,J)

            !.....LOWER ELIMINATION IN B
            do K=JP1,M
                do KK=JP1,M
                    B(I,K,KK) = B(I,K,KK)-B(I,K,J)*B(I,J,KK)
                enddo

                do KK=1,M
                    if(I.NE.M)then
                        C(I,K,KK) = C(I,K,KK)-B(I,K,J)*C(I,J,KK)
                    endif
                enddo

                R(I,K) = R(I,K)-B(I,K,J)*R(I,J)

            enddo

        enddo


        !.....UPPER (GAUSS JORDAN) ELIMINATION IN B

        T= 1.d0/B(I,M,M)

        do KK=1,M
            C(I,M,KK) = T*C(I,M,KK)
        enddo
        R(I,M) = T*R(I,M)

        MM1=M-1

        do J=1,MM1

            MP = M-J+1
            MPM1=MP-1

            do K=1,MPM1

                MR = MP-K

                do KK=1,M
                    if(I.NE.N)then
                        C(I,MR,KK) = C(I,MR,KK)-B(I,MR,MP)*C(I,MP,KK)
                    endif
                enddo

                R(I,MR) = R(I,MR)-B(I,MR,MP)*R(I,MP)

            enddo

        enddo

        !.....B(I) IS NOW THE UNIT MATRIX
        !.....ELIMINATE A(I+1)
        if(I.NE.N)then

            do J=1,M
                do K=1,M
                    do KK=1,M
                        B(I+1,K,KK) = B(I+1,K,KK)-A(I+1,K,J)*C(I,J,KK)
                    enddo
                    R(I+1,K) = R(I+1,K)-A(I+1,K,J)*R(I,J)
                enddo
            enddo

        endif

    enddo

    !.....BACK SUBSTITUTION
    do K=1,M
        X(N,K) = R(N,K)
    enddo

    NM1 = N-1
    do J=1,NM1
        JB=N-J
        do K=1,M
            KR = M - K + 1
            S=0.d0
            do KK=1,M
                S = C(JB,KR,KK)*X(JB+1,KK)+S
            enddo
            X(JB,KR) = R(JB,KR)-S
        enddo
    enddo

end subroutine


!==============================================================================================================================================

subroutine eval_fn_system(ngas,nlay,A,B,C,N,typelbc,valuelbc,typeubc,valueubc,fix_species,delz,X)
    ! This program evaluates the right hand side of the continuity equation using the
    ! current density array

    ! A, B and C are the sub-diagonal, diagonal and superdiagonal arrays of the transport/chemistry jacobian matrix
    ! which is of the form of a block tridiagonal matrix

    ! N are the densities of each species at each layer

    ! The rest of the parameters define the boundary conditions and the species to be fixed

    !Input
    integer, intent(in) :: nlay, ngas
    double precision, intent(in) :: A(nlay,ngas,ngas),B(nlay,ngas,ngas),C(nlay,ngas,ngas)
    double precision, intent(in) :: N(nlay,ngas)          !Density (cm-3)
    double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)     !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
    real, intent(in) :: delz                  !Layer altitude (m)
    integer, intent(in) :: typelbc(ngas),typeubc(ngas)    !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
    integer, intent(in) :: fix_species(nlay,ngas)         !Flag to indicate if density must be fixed (if 1)

    !Local
    integer :: igas,jgas,ilay,jlay
    double precision :: FN(nlay,ngas)

    !Output
    double precision, intent(out) :: X(nlay,ngas)


    ! Bottom layer
    !###################################################################################

    ilay = 1

    FN(ilay,:) = MATMUL(B(ilay,:,:),N(ilay,:)) + MATMUL(C(ilay,:,:),N(ilay+1,:))

    do igas=1,ngas

        if(typelbc(igas).eq.1)then  !Fixed density

            X(ilay,igas) = FN(ilay,igas) + valuelbc(igas) * 1.0d-6

        elseif(typelbc(igas).eq.2)then  !Fixed flux

            X(ilay,igas) = FN(ilay,igas) + (valuelbc(igas)*1.0d-4)/(delz*1.0d2)  

        endif
    enddo

    ! Top layer
    !###################################################################################

    ilay = nlay

    FN(ilay,:) = MATMUL(A(ilay,:,:),N(ilay-1,:)) + MATMUL(B(ilay,:,:),N(ilay,:))

    do igas=1,ngas

        if(typeubc(igas).eq.1)then  !Fixed density

            X(ilay,igas) = FN(ilay,igas) + valueubc(igas) * 1.0d-6

        elseif(typeubc(igas).eq.2)then  !Fixed flux

            X(ilay,igas) = FN(ilay,igas) - (valueubc(igas)*1.0d-4)/(delz*1.0d2)  

        elseif(typeubc(igas).eq.3)then  !Fixed velocity

            X(ilay,igas) = FN(ilay,igas) - N(ilay,igas) * (valueubc(igas)*1.0d2)/(delz*1.0d2)

        endif
    enddo


    ! In-between layers
    !###################################################################################

    do ilay=2,nlay-1

        X(ilay,:) = MATMUL(A(ilay,:,:),N(ilay-1,:)) + MATMUL(B(ilay,:,:),N(ilay,:)) + MATMUL(C(ilay,:,:),N(ilay+1,:))

    enddo


end subroutine