
MODULE mars_c13_chemical_network

    integer, save :: ngas_phot = 8         !number of gases to be photolysed

CONTAINS

    
!============================================================================================================================

    subroutine number_reactions(photolysis_inc,nreactions,n_photolysis)
    
        !Function to calculate the number of reactions in the chemical network

        logical, intent(in) :: photolysis_inc

        integer, intent(out) :: nreactions,n_photolysis
        
        
        !Determining the number of reactions depending on the flags
        !#####################################################################################

        !Photolysis
        n_photolysis = 0
        if(photolysis_inc) n_photolysis = 12

        !Chemistry
        nreactions = 31 + n_photolysis
    
    end subroutine number_reactions
    
    
!=======================================================================================================

    subroutine load_reactions(nlay,ngas,gasID,isoID,h,P,T,N,nreactions,n_photolysis,&
        rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf)

        !Loading the chemical network

        !Each reaction is defined by:

            !rtype :: Reaction type
            !   1 =     a + hv ---> b + c   or   a + c ---> b + c
            !   2 =     a + a ---> b + c
            !   3 =     a + b ---> c + d

            !ns :: Number of source species in the reaction (either 1 or 2)
            !sf(ns) :: Number of molecules of each source species 
            !sID(ns) :: Gas ID of each source species
            !sISO(ns) :: Isotope ID of each source species
            !npr :: Number of product species in the reaction
            !pf(npr) :: Number of molecules of each product
            !pID(ns) :: Gas ID of each product species
            !pISO(ns) :: Isotope ID of each product species

            !rrates(nlay) :: Reaction rates

        use reactions
        
        !Inputs
        integer, intent(in) :: nlay,ngas,nreactions,n_photolysis
        double precision, intent(in) :: P(nlay)           !Pressure (Pa)
        double precision, intent(in) :: N(nlay,ngas)      !Density of each species (m-3)
        integer, intent(in) :: gasID(ngas),isoID(ngas)    !Gas and isotopes IDs
        real, intent(in) :: h(nlay),T(nlay)               !Altitude (m) and Temperature (K)

        !Local 
        double precision :: dens(nlay),co2(nlay),o2(nlay),o(nlay),n2(nlay)
        integer :: ilay 
        character (len = 100) :: ref
        
        !Output
        double precision, intent(out) :: rrates(nlay,nreactions)
        integer, intent(out) :: rtype(nreactions)
        integer, intent(out) :: ns(nreactions),sID(2,nreactions),sISO(2,nreactions)
        integer, intent(out) :: npr(nreactions),pID(2,nreactions),pISO(2,nreactions)
        real, intent(out) :: sf(2,nreactions),pf(2,nreactions)
        
        
        !Initial calculations
        !##########################################################################################
        
        !Calculating the total atmospheric density in cm-3
        do ilay=1,nlay
            dens(ilay) = 0.d0
            do igas=1,ngas
                dens(ilay) = dens(ilay) + N(ilay,igas) * 1.0d-6
            enddo
        enddo
        
        !Calculating the number density of certain species (the ones as ambient gas in three-body reactions)
        do igas=1,ngas

            if((gasID(igas).eq.2).and.(isoID(igas).eq.0))then
                co2(:) = N(:,igas) * 1.0d-6
                
            elseif((gasID(igas).eq.7).and.(isoID(igas).eq.0))then
                o2(:) = N(:,igas) * 1.0d-6

            elseif((gasID(igas).eq.22).and.(isoID(igas).eq.0))then
                n2(:) = N(:,igas) * 1.0d-6

            elseif((gasID(igas).eq.45).and.(isoID(igas).eq.0))then
                o(:) = N(:,igas) * 1.0d-6

            endif

        enddo
        
        !Initialising the arrays
        rrates(:,:) = 0.d0
        ns(:) = 0
        npr(:) = 0
        sID(:,:) = 0
        pID(:,:) = 0
        sISO(:,:) = 0
        pISO(:,:) = 0
        sf(:,:) = 0.0
        pf(:,:) = 0.0
        rtype(:) = 0 
        
        !Starting to include reactions after the photolysis ones
        i0 = n_photolysis
    
        !O + O2 + CO2 -> O3 + CO2
        i0 = i0 + 1
        call reaction0001(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O + O + CO2 -> O2 + CO2
        i0 = i0 + 1
        call reaction0002(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O + O3 -> O2 + O2
        i0 = i0 + 1
        call reaction0003(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O(1D) + CO2 -> O + CO2
        i0 = i0 + 1
        call reaction0004(nlay,p(:),t(:),co2(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O(1D) + H2O -> OH + OH
        i0 = i0 + 1
        call reaction0005(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O(1D) + H2 -> OH + H
        i0 = i0 + 1
        call reaction0006(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O(1D) + O2 -> O + O2
        i0 = i0 + 1
        call reaction0007(nlay,p(:),t(:),o2(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O(1D) + O3  -> O2 + O2
        i0 = i0 + 1
        call reaction0008(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O(1D) + O3  -> O2 + O + O
        i0 = i0 + 1
        call reaction0009(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O + HO2 -> OH + O2
        i0 = i0 + 1
        call reaction0010(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O + OH -> O2 + H
        i0 = i0 + 1
        call reaction0011(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !H + O3 -> OH + O2
        i0 = i0 + 1
        call reaction0012(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !H + HO2 -> OH + OH
        i0 = i0 + 1
        call reaction0013(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !H + HO2 -> H2 + O2
        i0 = i0 + 1
        call reaction0014(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !H + HO2 -> H2O + O
        i0 = i0 + 1
        call reaction0015(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !OH + HO2 -> H2O + O2
        i0 = i0 + 1
        call reaction0016(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !HO2 + HO2 -> H2O2 + O2
        i0 = i0 + 1
        call reaction0017(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !OH + H2O2 -> H2O + HO2
        i0 = i0 + 1
        call reaction0018(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !OH + H2 -> H2O + H
        i0 = i0 + 1
        call reaction0019(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !H + O2 + CO2 -> HO2 + CO2
        i0 = i0 + 1
        call reaction0020(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O + H2O2 -> OH + HO2
        i0 = i0 + 1
        call reaction0021(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !OH + OH -> H2O + O
        i0 = i0 + 1
        call reaction0022(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !OH + O3 -> HO2 + O2
        i0 = i0 + 1
        call reaction0023(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !HO2 + O3 -> OH + O2 + O2
        i0 = i0 + 1
        call reaction0024(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !HO2 + HO2 + CO2 -> H2O2 + O2 + CO2
        i0 = i0 + 1
        call reaction0025(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !OH + OH + CO2 -> H2O2 + CO2
        i0 = i0 + 1
        call reaction0026(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !H + H + CO2 -> H2 + CO2
        i0 = i0 + 1
        call reaction0027(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !OH + CO -> CO2 + H
        i0 = i0 + 1
        call reaction0040(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !O + CO + M -> CO2 + M
        i0 = i0 + 1
        call reaction0041(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        
        !13-C reactions
        !##############################################
        
        !OH + (13C)O -> (13C)O2 + H
        i0 = i0 + 1
        call reaction0040_13c(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

            
        !O + (13C)O + M -> (13C)O2 + M
        i0 = i0 + 1
        call reaction0041_13c(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)
            
        end subroutine load_reactions
    
!============================================================================================================================

    subroutine load_reactions_photolysis(ngas_phot,gasID_phot,isoID_phot)

        !Loading the photolysis reactions 

        !Inputs
        integer, intent(in) :: ngas_phot  !Number of gases to be photolysed

        !Outputs
        integer, intent(out) :: gasID_phot(ngas_phot),isoID_phot(ngas_phot)  !Gases to be photolysed
            
    
                gasID_phot(1) = 7   !O2
                isoID_phot(1) = 0   
        
        
                gasID_phot(2) = 2   !CO2
                isoID_phot(2) = 0   
        
        
                gasID_phot(3) = 3   !O3
                isoID_phot(3) = 0   
        
        
                gasID_phot(4) = 1   !H2O
                isoID_phot(4) = 0   
        
        
                gasID_phot(5) = 25   !H2O2
                isoID_phot(5) = 0   
        
        
                gasID_phot(6) = 44   !HO2
                isoID_phot(6) = 0   
        
        
                gasID_phot(7) = 39   !H2
                isoID_phot(7) = 0   
        
        
                gasID_phot(8) = 2   !(13C)O2
                isoID_phot(8) = 2   
        
        
    end subroutine load_reactions_photolysis
                
    
    
END MODULE mars_c13_chemical_network
    