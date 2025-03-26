
MODULE mars_c13_diffusion
    
    real, parameter :: G = 6.67199976e-11         !m3 kg-1 s-2
    real, parameter :: k_B = 1.38065e-23     !m2 kg s-2 K-1
    real, parameter :: N_A = 6.02214e+23     !mol-1
    
    real,parameter :: Rplanet = 3397200.0    !Radius (m)
    real,parameter :: Mplanet = 6.4169e+23    !Mass (kg)
    
CONTAINS

!==========================================================================================================
  
    subroutine calc_jacobian_diffusion(nlay,ngas,ksi,klsi,ksim1,klsim1,typelbc,valuelbc,typeubc,valueubc,fix_species,Jmat)

        !Function to calculate the jacobian matrix for diffusion

        !use omp_lib


        !Inputs
        double precision, intent(in) :: ksi(nlay,ngas)                  !Parameters from the diffusion equation (see Alday et al., 2023)
        double precision, intent(in) :: klsi(nlay,ngas)
        double precision, intent(in) :: ksim1(nlay,ngas)
        double precision, intent(in) :: klsim1(nlay,ngas)
        double precision, intent(in) :: valuelbc(ngas),valueubc(ngas)   !Boundary conditions (in m-3 ; m-2 s-1 ; m s-1)
        integer, intent(in) :: typelbc(ngas),typeubc(ngas)              !Type of boundary conditions (1 - Fixed density ; 2 - Fixed flux ; 3 - Fixed velocity)
        integer, intent(in) :: nlay,ngas
        integer, intent(in) :: fix_species(nlay,ngas)                   !Flag to indicate if density must be fixed (if 1)

        !Local 
        integer :: igas,ilay

        !Output
        double precision, intent(out) :: Jmat(nlay,nlay,ngas)

        
        !Filling the Jacobian matrix
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

        !Re-computing the Jacobian matrix if some species is fixed
        do igas=1,ngas
        do ilay=1,nlay
            if(fix_species(ilay,igas).eq.1)then
            Jmat(ilay,:,igas) = 0.d0
            endif
        enddo
        enddo

    end subroutine calc_jacobian_diffusion

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
        real, intent(out) :: D(nlay,ngas)           !Eddy diffusion coefficient (m2 s-1)
            
            
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
        real, intent(in) :: mmol(ngas)                 !Molecular weight (uma)
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
    
END MODULE mars_c13_diffusion

    