MODULE chemistry

    CONTAINS

!===========================================================================================================================

    subroutine solve_chemistry_asis(nlay,ngas,nreactions,rtype,ns,sID_pos,sf,&
        npr,pID_pos,pf,rrates,N,dt,chemrates)

        !Routine to calculate the values of the chemical Jacobian matrix 

        implicit none

        ! input
        integer, intent(in) :: nlay                      ! number of atmospheric layers
        integer, intent(in) :: ngas                      ! number of gases in the atmosphere
        integer, intent(in) :: nreactions                ! number of reactions
        integer, intent(in) :: rtype(nreactions)         ! reaction type (0,1,2)
        integer, intent(in) :: ns(nreactions)            ! number of source species in each reaction
        integer, intent(in) :: sID_pos(2,nreactions)     ! position of the source species in the density array
        real, intent(in) :: sf(2,nreactions)             ! number of molecules of each source species
        integer, intent(in) :: npr(nreactions)           ! number of product species in each reaction
        integer, intent(in) :: pID_pos(2,nreactions)     ! position of the product species in the density array
        real, intent(in) :: pf(2,nreactions)             ! number of molecules of each product species
        double precision, intent(in) :: rrates(nlay,nreactions) !Reaction rates (s-1)
        double precision, intent(in) :: N(nlay,ngas)     ! number density of each species (m-3)
        double precision, intent(in) :: dt               ! timestep (s)

        ! local
        double precision :: mat(ngas,ngas),mat1(ngas,ngas) !jacobian matrix
        double precision :: prod(ngas),loss(ngas)          !production/loss terms
        integer :: ilay,igas
        double precision :: cnew(ngas),N0(nlay),c(nlay,ngas),cold(ngas),cini(nlay,ngas)
        integer :: code, indx(nlay)
        double precision :: time,dt_guess,dt_corrected

        ! output
        double precision, intent(out) :: chemrates(nlay,ngas) ! chemistry tendencies (cm-3 s-1)

        external dgesv

        ! initialisations
        chemrates(:,:) = 0.d0

        !Converting the density to cm-3
        c(:,:) = N(:,:) * 1.0d-6 
        cini(:,:) = N(:,:) * 1.0d-6

        ! Calculating the total atmospheric density (cm-3)
        do ilay=1,nlay
            N0(ilay) = 0.d0
            do igas=1,ngas
              N0(ilay) = N0(ilay) + c(ilay,igas)
            enddo
        enddo

        ! Solving the system in each layer
        !!$OMP DO
        do ilay=1,nlay

            !Initialisations
            time = 0.d0
            dt_guess = dt
            cold(:) = c(ilay,:)

            !Internal loop for the chemistry
            do while (time < dt)


                !Filling the jacobian matrix
                call fill_chem_matrix(nlay, ngas, ilay, c,  &
                nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, mat1)

                !Calculating the production and loss terms
                call calc_prodloss_chemistry(nlay, ngas, ilay, c,  &
                nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, prod, loss)

                !Adaptative evaluation of the sub time step
                call define_dt_chemistry(ngas, dt_guess, cold, c(ilay,:), N0(ilay), prod, loss, dt_corrected)

                if (time + dt_corrected > dt) then
                    dt_corrected = dt - time
                end if

                !Calculating the matrix (I+dt*J) to solve the system (I+dt*J) * N(t+1) = N(t)
                mat(:,:) = mat1(:,:)*dt_corrected
                do igas = 1,ngas
                    mat(igas,igas) = 1.d0 + mat(igas,igas)
                end do

                !Solving the system
                cnew(:) = c(ilay,:)
                call dgesv(ngas,1,mat,ngas,indx,cnew,ngas,code)

                !eliminate small values
                where (cnew(:)/N0(ilay) < 1.e-30)
                    cnew(:) = 0.d0
                end where

                !update concentrations
                cold(:) = c(ilay,:)
                c(ilay,:) = cnew(:)

                time = time + dt_corrected
                dt_guess = dt_corrected

            enddo

        !calculating the tendencies
        chemrates(ilay,:) = (c(ilay,:) - cini(ilay,:))/dt

        enddo
        !!$OMP END DO


    end subroutine

!===========================================================================================================================

    subroutine solve_ode_chemistry(nlay,ngas,nreactions,rtype,ns,sID_pos,sf,&
        npr,pID_pos,pf,rrates,N,dt,chemrates)

        !Routine to calculate the values of the chemical Jacobian matrix 

        implicit none

        ! input
        integer, intent(in) :: nlay                      ! number of atmospheric layers
        integer, intent(in) :: ngas                      ! number of gases in the atmosphere
        integer, intent(in) :: nreactions                ! number of reactions
        integer, intent(in) :: rtype(nreactions)         ! reaction type (0,1,2)
        integer, intent(in) :: ns(nreactions)            ! number of source species in each reaction
        integer, intent(in) :: sID_pos(2,nreactions)     ! position of the source species in the density array
        real, intent(in) :: sf(2,nreactions)             ! number of molecules of each source species
        integer, intent(in) :: npr(nreactions)           ! number of product species in each reaction
        integer, intent(in) :: pID_pos(2,nreactions)     ! position of the product species in the density array
        real, intent(in) :: pf(2,nreactions)             ! number of molecules of each product species
        double precision, intent(in) :: rrates(nlay,nreactions) !Reaction rates (s-1)
        double precision, intent(in) :: N(nlay,ngas)     ! number density of each species (m-3)
        double precision, intent(in) :: dt               ! timestep (s)

        ! local
        double precision :: mat(ngas,ngas),mat1(ngas,ngas) !jacobian matrix
        integer :: ilay,igas
        double precision :: cnew(ngas),N0(nlay),c(nlay,ngas)
        integer :: code, indx(nlay)

        ! output
        double precision, intent(out) :: chemrates(nlay,ngas) ! chemistry tendencies (cm-3 s-1)

        external dgesv

        ! initialisations
        chemrates(:,:) = 0.d0

        !Converting the density to cm-3
        c(:,:) = N(:,:) * 1.0d-6 


        ! Calculating the total atmospheric density (cm-3)
        do ilay=1,nlay
            N0(ilay) = 0.d0
            do igas=1,ngas
              N0(ilay) = N0(ilay) + c(ilay,igas)
            enddo
        enddo

        ! Solving the system in each layer
        !!$OMP DO
        do ilay=1,nlay

            !Filling the jacobian matrix
            call fill_chem_matrix(nlay, ngas, ilay, c,  &
            nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, mat1)

            !Calculating the matrix (I+dt*J) to solve the system (I+dt*J) * N(t+1) = N(t)
            mat(:,:) = mat1(:,:)*dt
            do igas = 1,ngas
               mat(igas,igas) = 1.d0 + mat(igas,igas)
            end do

            !Solving the system
            cnew(:) = c(ilay,:)
            call dgesv(ngas,1,mat,ngas,indx,cnew,ngas,code)

            !eliminate small values
            !where (cnew(:)/N0(ilay) < 1.e-30)
            !    cnew(:) = 0.d0
            !end where

            !calculating the tendencies
            chemrates(ilay,:) = (cnew(:) - c(ilay,:))/dt

            !print*,chemrates(ilay,:)
            !pause

        enddo
        !!$OMP END DO


    end subroutine

!======================================================================

    subroutine solve_ode_chemistry_rosenbrock(nlay,ngas,nreactions,rtype,ns,sID_pos,sf,&
        npr,pID_pos,pf,rrates,N,dt,chemrates,err)

        !Routine to calculate the values of the chemical Jacobian matrix 

        implicit none

        ! input
        integer, intent(in) :: nlay                      ! number of atmospheric layers
        integer, intent(in) :: ngas                      ! number of gases in the atmosphere
        integer, intent(in) :: nreactions                ! number of reactions
        integer, intent(in) :: rtype(nreactions)         ! reaction type (0,1,2)
        integer, intent(in) :: ns(nreactions)            ! number of source species in each reaction
        integer, intent(in) :: sID_pos(2,nreactions)     ! position of the source species in the density array
        real, intent(in) :: sf(2,nreactions)             ! number of molecules of each source species
        integer, intent(in) :: npr(nreactions)           ! number of product species in each reaction
        integer, intent(in) :: pID_pos(2,nreactions)     ! position of the product species in the density array
        real, intent(in) :: pf(2,nreactions)             ! number of molecules of each product species
        double precision, intent(in) :: rrates(nlay,nreactions) !Reaction rates (s-1)
        double precision, intent(in) :: N(nlay,ngas)     ! number density of each species (m-3)
        double precision, intent(in) :: dt               ! timestep (s)

        ! local
        double precision :: Jmat(ngas,ngas),Mmat(ngas,ngas) !jacobian matrix
        integer :: ilay,igas,jgas
        double precision :: gamma
        double precision :: cnew(ngas),N0(nlay),c(nlay,ngas)
        double precision :: b1(ngas),g1(ngas),b2(ngas),g2(ngas),a2(ngas)
        integer :: code, indx(ngas)

        ! output
        double precision, intent(out) :: chemrates(nlay,ngas) ! chemistry tendencies (cm-3 s-1)
        double precision, intent(out) :: err(nlay,ngas)       ! local error (cm-3)

        external dgesv

        ! initialisations
        chemrates(:,:) = 0.d0

        !Converting the density to cm-3
        c(:,:) = N(:,:) * 1.0d-6 


        ! Calculating the total atmospheric density (cm-3)
        do ilay=1,nlay
            N0(ilay) = 0.d0
            do igas=1,ngas
              N0(ilay) = N0(ilay) + c(ilay,igas)
            enddo
        enddo

        ! Solving the system in each layer
        !!$OMP DO
        !!$OMP PARALLEL PRIVATE(ilay,igas,jgas,Jmat,b1,Mmat,g1,indx,code,a2,b2,g2,cnew) &
        !!$OMP SHARED(ngas,nlay,c,nreactions,rtype,ns,sID_pos, sf, npr, pID_pos, pf, rrates) &
        !!$OMP SHARED(gamma,dt,N0,err,chemrates) &
        !!$OMP NUM_THREADS(50)
        !!$OMP DO 

        do ilay=1,nlay

            !Filling the jacobian matrix
            call calc_jacobian_chemistry(nlay, ngas, ilay, c,  &
            nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, Jmat)

            !Calculating b1
            b1 = MATMUL(Jmat,c(ilay,:))

            !print*,'b1'
            !print*,b1(:)

            !Build matrix (I-gamma*dt*J) and solve for g1
            gamma = 1.d0 + 1.d0/sqrt(2.d0)

            !Building matrix
            do igas=1,ngas
                do jgas=1,ngas
                    if(igas.eq.jgas)then
                        Mmat(igas,jgas) = 1.0d0 - gamma*dt*Jmat(igas,jgas)
                    else
                        Mmat(igas,jgas) = -gamma*dt*Jmat(igas,jgas)
                    endif
                enddo
            enddo

            !Solve for g1   (I-gamma*dt*J) * g1 = b1
            g1(:) = b1(:)
            call dgesv(ngas,1,Mmat,ngas,indx,g1,ngas,code)


            !Calculate b2 = F(N+dt*g1)
            do igas=1,ngas
                a2(igas) = c(ilay,igas)+dt*g1(igas)
            enddo
            b2 = MATMUL(Jmat,a2(:))

            !print*,'b2'
            !print*,b2(:)

            !Solve for g2   (I-gamma*dt*J) * g2 = b2 - 2*g1
            g2(:) = b2(:) - 2.d0*g1(:)
            call dgesv(ngas,1,Mmat,ngas,indx,g2,ngas,code)

            !print*,'g2'
            !print*,g2(:)
      
            !Finding the new density
            do igas=1,ngas
                cnew(igas) = c(ilay,igas) + 1.5*dt*g1(igas) + 0.5*dt*g2(igas)
            enddo

            !print*,cnew(:)

            !eliminate small values
            where (cnew(:)/N0(ilay) < 1.e-30)
                cnew(:) = 0.d0
            end where


            !Computing the local error
            do igas=1,ngas
                  err(ilay,igas) = abs(cnew(igas) - (c(ilay,igas) + dt*g1(igas)))
            enddo

            !calculating the tendencies
            chemrates(ilay,:) = (cnew(:) - c(ilay,:))/dt

            !do igas=1,ngas
            !    chemrates(ilay,igas) = 1.5*g1(igas) + 0.5*g2(igas)
            !enddo

            !print*,'chemrates',ilay
            !print*,chemrates(ilay,:)
            !pause

        enddo
        !!$OMP END DO
        !!$OMP END DO 
        !!$OMP END PARALLEL

    end subroutine solve_ode_chemistry_rosenbrock

!==========================================================================================================================

    subroutine fill_chem_matrix(nlay, ngas, ilay, c,  &
        nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, mat)

        !Routine to calculate the values of the chemical Jacobian matrix 

        implicit none

        ! input
        integer, intent(in) :: ilay    ! level index
        integer, intent(in) :: ngas    ! number of species in the chemistry
        integer, intent(in) :: nlay    ! number of atmospheric layers

        integer, intent(in) :: nreactions  !number of reactions
        integer, intent(in) :: rtype(nreactions)   !reaction tpye
        double precision, intent(in) :: rrates(nlay,nreactions)  !reaction rates
        integer, intent(in) :: ns(nreactions),sID_pos(2,nreactions)  !source species
        integer, intent(in) :: npr(nreactions),pID_pos(2,nreactions) !product species
        real, intent(in) :: sf(2,nreactions),pf(2,nreactions) !number of molecules for each source/product
        double precision, intent(in) :: c(nlay,ngas) !Number density of each species (cm-3)

        ! local
        integer :: ir
        integer :: ind_phot_2,ind_phot_4,ind_phot_6
        integer :: ind_3_2,ind_3_4,ind_3_6
        integer :: ind_4_2,ind_4_4,ind_4_6,ind_4_8
        double precision :: eps, eps_4


        ! output
        double precision, intent(out) :: mat(ngas,ngas) !matrix


        ! initialisations
        mat(:,:) = 0.d0


        do ir=1,nreactions

            if(rtype(ir).eq.1)then

                ! photodissociations
                ! or reactions a + c -> b + c
                ! or reactions a + ice -> b + c
                !################################################################################

                ind_phot_2 = sID_pos(1,ir)
                ind_phot_4 = pID_pos(1,ir)
                ind_phot_6 = pID_pos(2,ir)

                mat(ind_phot_2,ind_phot_2) = mat(ind_phot_2,ind_phot_2) + sf(1,ir)*rrates(ilay,ir)

                if(npr(ir).eq.1)then
                    mat(ind_phot_4,ind_phot_2) = mat(ind_phot_4,ind_phot_2) - pf(1,ir)*rrates(ilay,ir)
                elseif(npr(ir).eq.2)then
                    mat(ind_phot_4,ind_phot_2) = mat(ind_phot_4,ind_phot_2) - pf(1,ir)*rrates(ilay,ir)
                    mat(ind_phot_6,ind_phot_2) = mat(ind_phot_6,ind_phot_2) - pf(2,ir)*rrates(ilay,ir)
                endif

            elseif(rtype(ir).eq.2)then

                ! reactions a + a -> b + c
                !################################################################################

                ind_3_2 = sID_pos(1,ir)
                ind_3_4 = pID_pos(1,ir)
                ind_3_6 = pID_pos(2,ir)

                mat(ind_3_2,ind_3_2) = mat(ind_3_2,ind_3_2) + sf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)

                if(npr(ir).eq.1)then
                    mat(ind_3_4,ind_3_2) = mat(ind_3_4,ind_3_2) - pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)
                elseif(npr(ir).eq.2)then
                    mat(ind_3_4,ind_3_2) = mat(ind_3_4,ind_3_2) - pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)
                    mat(ind_3_6,ind_3_2) = mat(ind_3_6,ind_3_2) - pf(2,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)
                endif

            elseif(rtype(ir).eq.3)then
                ! reactions a + b -> c + d
                !################################################################################        

                eps = 1.d-10

                ind_4_2 = sID_pos(1,ir)
                ind_4_4 = sID_pos(2,ir)
                ind_4_6 = pID_pos(1,ir)
                ind_4_8 = pID_pos(2,ir)

                eps_4 = abs(c(ilay,ind_4_2))/(abs(c(ilay,ind_4_2)) + abs(c(ilay,ind_4_4)) + eps)
                eps_4 = min(eps_4,1.d0)

                mat(ind_4_2,ind_4_2) = mat(ind_4_2,ind_4_2) + sf(1,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                mat(ind_4_2,ind_4_4) = mat(ind_4_2,ind_4_4) + sf(1,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                mat(ind_4_4,ind_4_2) = mat(ind_4_4,ind_4_2) + sf(2,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                mat(ind_4_4,ind_4_4) = mat(ind_4_4,ind_4_4) + sf(2,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)

                if(npr(ir).eq.1)then
                    mat(ind_4_6,ind_4_2) = mat(ind_4_6,ind_4_2) - pf(1,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                    mat(ind_4_6,ind_4_4) = mat(ind_4_6,ind_4_4) - pf(1,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                elseif(npr(ir).eq.2)then
                    mat(ind_4_6,ind_4_2) = mat(ind_4_6,ind_4_2) - pf(1,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                    mat(ind_4_6,ind_4_4) = mat(ind_4_6,ind_4_4) - pf(1,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                    mat(ind_4_8,ind_4_2) = mat(ind_4_8,ind_4_2) - pf(2,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                    mat(ind_4_8,ind_4_4) = mat(ind_4_8,ind_4_4) - pf(2,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                endif

            else

                print*,'error, reaction type must be 1,2 or 3'
                print*,'reaction',ir,'. type',rtype(ir)
                stop

            endif

        enddo

    end subroutine fill_chem_matrix


!===============================================================================================================================================

    subroutine calc_jacobian_chemistry(nlay, ngas, ilay, c,  &
        nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, Jmat)

        !Routine to calculate the values of the chemical Jacobian matrix 

        implicit none

        ! input
        integer, intent(in) :: ilay    ! level index
        integer, intent(in) :: ngas    ! number of species in the chemistry
        integer, intent(in) :: nlay    ! number of atmospheric layers

        integer, intent(in) :: nreactions  !number of reactions
        integer, intent(in) :: rtype(nreactions)   !reaction tpye
        double precision, intent(in) :: rrates(nlay,nreactions)  !reaction rates
        integer, intent(in) :: ns(nreactions),sID_pos(2,nreactions)  !source species
        integer, intent(in) :: npr(nreactions),pID_pos(2,nreactions) !product species
        real, intent(in) :: sf(2,nreactions),pf(2,nreactions) !number of molecules for each source/product
        double precision, intent(in) :: c(nlay,ngas) !Number density of each species (cm-3)

        ! local
        integer :: ir
        integer :: ind_phot_2,ind_phot_4,ind_phot_6
        integer :: ind_3_2,ind_3_4,ind_3_6
        integer :: ind_4_2,ind_4_4,ind_4_6,ind_4_8
        double precision :: eps, eps_4


        ! output
        double precision, intent(out) :: Jmat(ngas,ngas) !matrix


        ! initialisations
        Jmat(:,:) = 0.d0


        do ir=1,nreactions

            if(rtype(ir).eq.1)then

                ! photodissociations
                ! or reactions a + c -> b + c
                ! or reactions a + ice -> b + c
                !################################################################################

                ind_phot_2 = sID_pos(1,ir)
                ind_phot_4 = pID_pos(1,ir)
                ind_phot_6 = pID_pos(2,ir)

                Jmat(ind_phot_2,ind_phot_2) = Jmat(ind_phot_2,ind_phot_2) - sf(1,ir)*rrates(ilay,ir)

                if(npr(ir).eq.1)then
                    Jmat(ind_phot_4,ind_phot_2) = Jmat(ind_phot_4,ind_phot_2) + pf(1,ir)*rrates(ilay,ir)
                elseif(npr(ir).eq.2)then
                    Jmat(ind_phot_4,ind_phot_2) = Jmat(ind_phot_4,ind_phot_2) + pf(1,ir)*rrates(ilay,ir)
                    Jmat(ind_phot_6,ind_phot_2) = Jmat(ind_phot_6,ind_phot_2) + pf(2,ir)*rrates(ilay,ir)
                endif

            elseif(rtype(ir).eq.2)then

                ! reactions a + a -> b + c
                !################################################################################

                ind_3_2 = sID_pos(1,ir)
                ind_3_4 = pID_pos(1,ir)
                ind_3_6 = pID_pos(2,ir)

                Jmat(ind_3_2,ind_3_2) = Jmat(ind_3_2,ind_3_2) - sf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)

                if(npr(ir).eq.1)then
                    Jmat(ind_3_4,ind_3_2) = Jmat(ind_3_4,ind_3_2) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)
                elseif(npr(ir).eq.2)then
                    Jmat(ind_3_4,ind_3_2) = Jmat(ind_3_4,ind_3_2) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)
                    Jmat(ind_3_6,ind_3_2) = Jmat(ind_3_6,ind_3_2) + pf(2,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)
                endif

            elseif(rtype(ir).eq.3)then
                ! reactions a + b -> c + d
                !################################################################################        

                eps = 1.d-10

                ind_4_2 = sID_pos(1,ir)
                ind_4_4 = sID_pos(2,ir)
                ind_4_6 = pID_pos(1,ir)
                ind_4_8 = pID_pos(2,ir)

                eps_4 = abs(c(ilay,ind_4_2))/(abs(c(ilay,ind_4_2)) + abs(c(ilay,ind_4_4)) + eps)
                eps_4 = min(eps_4,1.d0)

                Jmat(ind_4_2,ind_4_2) = Jmat(ind_4_2,ind_4_2) - sf(1,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                Jmat(ind_4_2,ind_4_4) = Jmat(ind_4_2,ind_4_4) - sf(1,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                Jmat(ind_4_4,ind_4_2) = Jmat(ind_4_4,ind_4_2) - sf(2,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                Jmat(ind_4_4,ind_4_4) = Jmat(ind_4_4,ind_4_4) - sf(2,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)

                if(npr(ir).eq.1)then
                    Jmat(ind_4_6,ind_4_2) = Jmat(ind_4_6,ind_4_2) + pf(1,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                    Jmat(ind_4_6,ind_4_4) = Jmat(ind_4_6,ind_4_4) + pf(1,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                elseif(npr(ir).eq.2)then
                    Jmat(ind_4_6,ind_4_2) = Jmat(ind_4_6,ind_4_2) + pf(1,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                    Jmat(ind_4_6,ind_4_4) = Jmat(ind_4_6,ind_4_4) + pf(1,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                    Jmat(ind_4_8,ind_4_2) = Jmat(ind_4_8,ind_4_2) + pf(2,ir)*rrates(ilay,ir)*(1. - eps_4)*c(ilay,ind_4_4)
                    Jmat(ind_4_8,ind_4_4) = Jmat(ind_4_8,ind_4_4) + pf(2,ir)*rrates(ilay,ir)*eps_4*c(ilay,ind_4_2)
                endif

            else

                print*,'error, reaction type must be 1,2 or 3'
                print*,'reaction',ir,'. type',rtype(ir)
                stop

            endif

        enddo

    end subroutine calc_jacobian_chemistry

!================================================================================================================================

    subroutine calc_prodloss_chemistry(nlay, ngas, ilay, c,  &
        nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates, prod, loss)

        !Routine to calculate the values of the chemical Jacobian matrix 

        implicit none

        ! input
        integer, intent(in) :: ilay    ! level index
        integer, intent(in) :: ngas    ! number of species in the chemistry
        integer, intent(in) :: nlay    ! number of atmospheric layers

        integer, intent(in) :: nreactions  !number of reactions
        integer, intent(in) :: rtype(nreactions)   !reaction tpye
        double precision, intent(in) :: rrates(nlay,nreactions)  !reaction rates
        integer, intent(in) :: ns(nreactions),sID_pos(2,nreactions)  !source species
        integer, intent(in) :: npr(nreactions),pID_pos(2,nreactions) !product species
        real, intent(in) :: sf(2,nreactions),pf(2,nreactions) !number of molecules for each source/product
        double precision, intent(in) :: c(nlay,ngas) !Number density of each species (cm-3)

        ! local
        integer :: ir
        integer :: ind_phot_2,ind_phot_4,ind_phot_6
        integer :: ind_3_2,ind_3_4,ind_3_6
        integer :: ind_4_2,ind_4_4,ind_4_6,ind_4_8
        double precision :: eps, eps_4


        ! output
        double precision, intent(out) :: prod(ngas),loss(ngas) !production and loss


        ! initialisations
        prod(:) = 0.d0
        loss(:) = 0.d0


        do ir=1,nreactions

            if(rtype(ir).eq.1)then

                ! photodissociations
                ! or reactions a + c -> b + c
                ! or reactions a + ice -> b + c
                !################################################################################

                ind_phot_2 = sID_pos(1,ir)
                ind_phot_4 = pID_pos(1,ir)
                ind_phot_6 = pID_pos(2,ir)

                loss(ind_phot_2) = loss(ind_phot_2) + sf(1,ir)*rrates(ilay,ir)

                if(npr(ir).eq.1)then
                    prod(ind_phot_4) = prod(ind_phot_4) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_phot_2)
                elseif(npr(ir).eq.2)then
                    prod(ind_phot_4) = prod(ind_phot_4) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_phot_2)
                    prod(ind_phot_6) = prod(ind_phot_6) + pf(2,ir)*rrates(ilay,ir)*c(ilay,ind_phot_2)
                endif

            elseif(rtype(ir).eq.2)then

                ! reactions a + a -> b + c
                !################################################################################

                ind_3_2 = sID_pos(1,ir)
                ind_3_4 = pID_pos(1,ir)
                ind_3_6 = pID_pos(2,ir)

                loss(ind_3_2) = loss(ind_3_2) + sf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)

                if(npr(ir).eq.1)then
                    prod(ind_3_4) = prod(ind_3_4) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)*c(ilay,ind_3_2)
                elseif(npr(ir).eq.2)then
                    prod(ind_3_4) = prod(ind_3_4) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)*c(ilay,ind_3_2)
                    prod(ind_3_6) = prod(ind_3_6) + pf(2,ir)*rrates(ilay,ir)*c(ilay,ind_3_2)*c(ilay,ind_3_2)
                endif

            elseif(rtype(ir).eq.3)then

                ! reactions a + b -> c + d
                !################################################################################        

                eps = 1.d-10

                ind_4_2 = sID_pos(1,ir)
                ind_4_4 = sID_pos(2,ir)
                ind_4_6 = pID_pos(1,ir)
                ind_4_8 = pID_pos(2,ir)

                eps_4 = abs(c(ilay,ind_4_2))/(abs(c(ilay,ind_4_2)) + abs(c(ilay,ind_4_4)) + eps)
                eps_4 = min(eps_4,1.d0)

                loss(ind_4_2) = loss(ind_4_2) + sf(1,ir)*rrates(ilay,ir)*c(ilay,ind_4_4)
                loss(ind_4_4) = loss(ind_4_4) + sf(2,ir)*rrates(ilay,ir)*c(ilay,ind_4_2)

                if(npr(ir).eq.1)then
                    prod(ind_4_6) = prod(ind_4_6) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_4_2)*c(ilay,ind_4_4)
                elseif(npr(ir).eq.2)then
                    prod(ind_4_6) = prod(ind_4_6) + pf(1,ir)*rrates(ilay,ir)*c(ilay,ind_4_2)*c(ilay,ind_4_4)
                    prod(ind_4_8) = prod(ind_4_8) + pf(2,ir)*rrates(ilay,ir)*c(ilay,ind_4_2)*c(ilay,ind_4_4)
                endif

            else

                print*,'error, reaction type must be 1,2 or 3'
                print*,'reaction',ir,'. type',rtype(ir)
                stop

            endif

        enddo

    end subroutine calc_prodloss_chemistry

!================================================================================================================================

    subroutine define_dt_chemistry(ngas, dt_guess, cold, ccur, N0, prod, loss, dt_corrected)

        !================================================================
        ! iterative evaluation of the appropriate time step dtnew
        ! according to curvature criterion based on
        ! e = 2 Rtol [r Cn+1 -(1-r) Cn + Cn-1 ]/[(1+r) Cn]
        ! with r = (tn - tn-1)/(tn+1 - tn)
        !================================================================

        implicit none

        ! input
        integer, intent(in) :: ngas    ! number of species in the chemistry

        double precision, intent(in) :: dt_guess !First guess of time step
        double precision, intent(in) :: cold(ngas),ccur(ngas) !Number density of each species (cm-3)
        double precision, intent(in) :: N0 !Total number density (cm-3)
        double precision, intent(in) :: prod(ngas),loss(ngas)  !Production and loss terms

        ! local
        double precision :: dttest,atol,ratio,e,es,coef
        double precision :: cnew(ngas)
        integer :: iter,igas

        ! parameters

        real (kind = 8), parameter :: dtmin   = 1.e-3    ! minimum time step (s)
        real (kind = 8), parameter :: vmrtol  = 1.e-11   ! absolute tolerance on vmr
        real (kind = 8), parameter :: rtol    = 0.05     ! rtol recommended value : 0.1-0.02
        integer,         parameter :: niter   = 3        ! number of iterations
        real (kind = 8), parameter :: coefmax = 2.
        real (kind = 8), parameter :: coefmin = 0.1
        

        !output
        double precision, intent(out) :: dt_corrected


        dttest = dt_guess
        atol = vmrtol*N0

        do iter=1,niter

            !fast semi-imlicit method
            do igas=1,ngas
                cnew(igas) = (ccur(igas) + prod(igas)*dttest)/(1.d0 + loss(igas)*dttest)
            enddo

            ratio = dt_guess/dttest

            ! e : local error indicator
            e = 0.

            do igas = 1,ngas
                es = 2.*abs((ratio*cnew(igas) - (1. + ratio)*ccur(igas) + cold(igas))   &
                        /(1. + ratio)/max(ccur(igas)*rtol,atol))

                if (es > e) then
                    e = es
                end if
            end do

            ! timestep correction
            coef = max(coefmin, min(coefmax,0.8/sqrt(e)))

            dttest = max(dtmin,dttest*coef)
            dttest = min(dt_guess,dttest)

        enddo

        dt_corrected = dttest


    end subroutine define_dt_chemistry

    !======================================================================

    subroutine locate_gas_reactions(ngas,gasID,isoID,nreactions, ns, sID, sISO, npr, pID, pISO,&
        sID_pos,pID_pos)

        !Routine to find the location of the source/products in each reaction in the Gas ID array
        !defining the gases in the atmosphere 

        implicit none

        ! input
        integer, intent(in) :: ngas                      ! number of species in the chemistry
        integer, intent(in) :: gasID(ngas),isoID(ngas)   ! ID of the gases in the atmosphere

        integer, intent(in) :: nreactions  !number of reactions
        integer, intent(in) :: ns(nreactions),sID(2,nreactions),sISO(2,nreactions)  !source species
        integer, intent(in) :: npr(nreactions),pID(2,nreactions),pISO(2,nreactions) !product species


        !local
        integer :: ir,j,igas,igasx

        !output
        integer, intent(out) :: sID_pos(2,nreactions),pID_pos(2,nreactions)


        !Looping through each source
        do ir=1,nreactions

            do j=1,ns(ir)
                igasx = 0
                do igas=1,ngas
                    if((sID(j,ir).eq.gasID(igas)).and.(sISO(j,ir).eq.isoID(igas)))then
                        sID_pos(j,ir) = igas
                        igasx = 1
                    endif
                enddo
                if(igasx.eq.0)then
                    print*,'error, reaction involves a gas not present in the atmosphere (source)'
                    print*,'reaction',ir,nreactions
                    print*,'gasID,isoID',sID(j,ir),sISO(j,ir)
                    stop
                endif
            enddo

            do j=1,npr(ir)
                igasx = 0
                do igas=1,ngas
                    if((pID(j,ir).eq.gasID(igas)).and.(pISO(j,ir).eq.isoID(igas)))then
                        pID_pos(j,ir) = igas
                        igasx = 1
                    endif
                enddo
                if(igasx.eq.0)then
                    print*,'error, reaction involves a gas not present in the atmosphere (product)'
                    print*,'reaction',ir,nreactions
                    print*,'gasID,isoID',pID(j,ir),pISO(j,ir)
                    stop
                endif
            enddo

        enddo

end subroutine



END MODULE chemistry