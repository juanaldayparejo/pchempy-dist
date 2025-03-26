MODULE mars_chemistry

    CONTAINS

!============================================================================================================================

    subroutine number_reactions(nitrogen_chemistry,c13_chemistry,o18_chemistry,photolysis,&
                                nreactions,n_photolysis)

        !Function to calculate the number of reactions in the chemistry

        logical, intent(in) :: nitrogen_chemistry, c13_chemistry, o18_chemistry, photolysis

        integer, intent(out) ::nreactions,n_photolysis


        !Determining the number of reactions depending on the flags
        !#####################################################################################

        !Photolysis
        n_photolysis = 0
        if(photolysis)then

            n_photolysis = 10

            if(nitrogen_chemistry)then
                n_photolysis = n_photolysis + 3
            endif

            if(c13_chemistry)then
                n_photolysis = n_photolysis + 2
            endif

            if(o18_chemistry)then
                n_photolysis = n_photolysis + 4
            endif

        endif

        !Chemistry
        nreactions = 29 + n_photolysis

        if(nitrogen_chemistry)then
            nreactions = nreactions + 12
        endif

        if(c13_chemistry)then
            nreactions = nreactions + 2

            if(nitrogen_chemistry)then
                nreactions = nreactions + 1
            endif
        endif

        if(o18_chemistry)then
            nreactions = nreactions + 3

            if(nitrogen_chemistry)then
                nreactions = nreactions + 1
            endif
        endif

    end subroutine

!============================================================================================================================

    subroutine load_reactions(nlay,ngas,gasID,isoID,h,P,T,N,nitrogen_chemistry,c13_chemistry,o18_chemistry,photolysis,&
        nreactions,n_photolysis,rrates,rtype,ns,sID,sISO,sf,npr,pID,pISO,pf)

        !Loading the reactions relevant for the atmosphere of Mars

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
        double precision, intent(in) :: h(nlay),T(nlay)   !Altitude (m) and Temperature (K)
        logical, intent(in) :: nitrogen_chemistry, c13_chemistry, o18_chemistry, photolysis


        !Local 
        double precision :: dens(nlay),co2(nlay),o2(nlay),o(nlay),n2(nlay)
        integer :: ilay
        double precision :: ak0, ak1, xpo, rate, rate1, rate2
        double precision :: k0,kinf,kint,kf,kca
        double precision :: o18ratio_oh   !(18O)/(16O) ratio in OH to tweak the model 
        character (len = 100) :: ref


        !Output
        double precision, intent(out) :: rrates(nlay,nreactions)
        integer, intent(out) :: rtype(nreactions)
        integer, intent(out) :: ns(nreactions),sID(2,nreactions),sISO(2,nreactions)
        integer, intent(out) :: npr(nreactions),pID(2,nreactions),pISO(2,nreactions)
        real, intent(out) :: sf(2,nreactions),pf(2,nreactions)


    
        !Input to tweak the O18 chemistry
        o18ratio_oh = 1.048

        !Adding the reactions
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

        !We leave space for the first reactions to be filled with the photolysis

        !O + O2 + CO2 -> O3 + CO2
        i0 = n_photolysis + 1
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
        call reaction0004(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        !O(1D) + H2O -> OH + OH
        i0 = i0 + 1
        call reaction0005(nlay,p(:),t(:),dens(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)



        !==========================================================================================
        !---b002: o(1d) + h2o  -> oh + oh

        !   jpl 2006

        i0 = i0 + 1

        b002(:) = 1.63e-10*exp(60./t(:))

        rrates(:,i0) = b002(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 133
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 1
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 1
        pID(1,i0) = 13
        pISO(1,i0) = 0
        pf(1,i0) = 2.0


        !==========================================================================================
        !---b003: o(1d) + h2  -> oh + h

        !   jpl 2011

        i0 = i0 + 1

        b003(:) = 1.2e-10

        rrates(:,i0) = b003(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 133
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 39
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 13
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 48
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !==========================================================================================
        !---b004: o(1d) + o2  -> o + o2

        !   jpl 2006

        i0 = i0 + 1

        b004(:) = 3.3e-11*exp(55./t(:))

        rrates(:,i0) = b004(:)*o2(:)
        rtype(i0) = 1

        ns(i0) = 1
        sID(1,i0) = 133
        sISO(1,i0) = 0
        sf(1,i0) = 1.0

        npr(i0) = 1
        pID(1,i0) = 45
        pISO(1,i0) = 0
        pf(1,i0) = 1.0


        !==========================================================================================
        !---b005: o(1d) + o3  -> o2 + o2

        !   jpl 2003

        i0 = i0 + 1

        b005(:) = 1.2e-10

        rrates(:,i0) = b005(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 133
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 3
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 1
        pID(1,i0) = 7
        pISO(1,i0) = 0
        pf(1,i0) = 2.0


        !==========================================================================================
        !---b006: o(1d) + o3  -> o2 + o + o

        !   jpl 2003


        i0 = i0 + 1

        b006(:) = 1.2e-10

        rrates(:,i0) = b006(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 133
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 3
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 7
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 45
        pISO(2,i0) = 0
        pf(2,i0) = 2.0


        !==========================================================================================
        !HYDROGEN REACTIONS
        !==========================================================================================

        !=========================================================================================
        !---c001: o + ho2 -> oh + o2

        !   jpl 2003

        i0 = i0 + 1

        c001(:) = 3.0e-11*exp(200./t(:))

        rrates(:,i0) = c001(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 45
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 44
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 7
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 13
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !--- c002: o + oh -> o2 + h

        !   jpl 2011

        i0 = i0 + 1

        c002(:) = 1.8e-11*exp(180./t(:))

        rrates(:,i0) = c002(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 45
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 13
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 7
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 48
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c003: h + o3 -> oh + o2

        !   jpl 2003

        i0 = i0 + 1

        c003(:) = 1.4e-10*exp(-470./t(:))

        rrates(:,i0) = c003(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 48
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 3
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 13
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 7
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c004: h + ho2 -> oh + oh

        !   jpl 2006


        i0 = i0 + 1

        c004(:) = 7.2e-11

        rrates(:,i0) = c004(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 48
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 44
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 1
        pID(1,i0) = 13
        pISO(1,i0) = 0
        pf(1,i0) = 2.0


        !=========================================================================================
        !---c005: h + ho2 -> h2 + o2

        !   jpl 2006

        i0 = i0 + 1

        c005(:) = 6.9e-12

        rrates(:,i0) = c005(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 48
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 44
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 39
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 7
        pISO(2,i0) = 0
        pf(2,i0) = 1.0

        !=========================================================================================
        !---c006: h + ho2 -> h2o + o

        !   jpl 2006


        i0 = i0 + 1

        c006(:) = 1.6e-12

        rrates(:,i0) = c006(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 48
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 44
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 1
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 45
        pISO(2,i0) = 0
        pf(2,i0) = 1.0

        !=========================================================================================
        !---c007: oh + ho2 -> h2o + o2

        !   jpl 2003

        !   canty et al., grl, 2006 suggest to increase this rate
        !   by 20%. not done here.
 
        i0 = i0 + 1

        c007(:) = 4.8e-11*exp(250./t(:))

        rrates(:,i0) = c007(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 13
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 44
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 1
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 7
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c008: ho2 + ho2 -> h2o2 + o2

        !   jpl 2015

        !       c008(:) = 3.0e-13*exp(460./t(:))

        !   christensen et al., grl, 13, 2002
        
        !       c008(:) = 1.5e-12*exp(19./t(:))


        i0 = i0 + 1

        c008(:) = 3.0e-13*exp(460./t(:))

        rrates(:,i0) = c008(:)
        rtype(i0) = 2

        ns(i0) = 1
        sID(1,i0) = 44
        sISO(1,i0) = 0
        sf(1,i0) = 2.0

        npr(i0) = 2
        pID(1,i0) = 25
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 7
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c009: oh + h2o2 -> h2o + ho2

        !   jpl 2006


        i0 = i0 + 1

        c009(:) = 1.8e-12

        rrates(:,i0) = c009(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 13
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 25
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 1
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 44
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c010: oh + h2 -> h2o + h

        !   jpl 2006

        i0 = i0 + 1

        c010(:) = 2.8e-12*exp(-1800./t(:))

        rrates(:,i0) = c010(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 13
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 39
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 1
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 48
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c011: h + o2 + co2 -> ho2 + co2

        !   jpl 2011
        !   co2/n2 efficiency as a third body = 2.4
        !   from ashman and haynes, 27th symposium on combustion, 1998.

        i0 = i0 + 1

        do ilay = 1,nlay
            !ak0 = 3.1*2.4*4.4e-32*(t(ilev)/300.)**(-1.3) ! FL li et al 2017
            ak0 = 2.4*4.4e-32*(t(ilay)/300.)**(-1.3)
            ak1 = 7.5e-11*(t(ilay)/300.)**(0.2)
            
            rate = (ak0*dens(ilay))/(1. + ak0*dens(ilay)/ak1)
            xpo = 1./(1. + dlog10((ak0*dens(ilay))/ak1)**2)
            c011(ilay) = rate*0.6**xpo
        end do


        rrates(:,i0) = c011(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 48
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 7
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 1
        pID(1,i0) = 44
        pISO(1,i0) = 0
        pf(1,i0) = 1.0


        !=========================================================================================
        !---c012: o + h2o2 -> oh + ho2

        !   jpl 2003


        i0 = i0 + 1

        c012(:) = 1.4e-12*exp(-2000./t(:))

        rrates(:,i0) = c012(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 45
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 25
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 13
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 44
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c013: oh + oh -> h2o + o

        !   jpl 2006


        i0 = i0 + 1

        c013(:) = 1.8e-12

        rrates(:,i0) = c013(:)
        rtype(i0) = 2

        ns(i0) = 1
        sID(1,i0) = 13
        sISO(1,i0) = 0
        sf(1,i0) = 2.0

        npr(i0) = 2
        pID(1,i0) = 1
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 45
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c014: oh + o3 -> ho2 + o2

        !   jpl 2003

        i0 = i0 + 1

        c014(:) = 1.7e-12*exp(-940./t(:))

        rrates(:,i0) = c014(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 13
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 3
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 44
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 7
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c015: ho2 + o3 -> oh + o2 + o2

        !   jpl 2003

        i0 = i0 + 1

        c015(:) = 1.0e-14*exp(-490./t(:))

        rrates(:,i0) = c015(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 44
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 3
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 13
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 7
        pISO(2,i0) = 0
        pf(2,i0) = 2.0


        !=========================================================================================
        !---c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2

        !   jpl 2011


        i0 = i0 + 1

        c016(:) = 2.5*2.1e-33*exp(920./t(:))*dens(:)

        rrates(:,i0) = c016(:)
        rtype(i0) = 2

        ns(i0) = 1
        sID(1,i0) = 44
        sISO(1,i0) = 0
        sf(1,i0) = 2.0

        npr(i0) = 2
        pID(1,i0) = 25
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 7
        pISO(2,i0) = 0
        pf(2,i0) = 1.0


        !=========================================================================================
        !---c017: oh + oh + co2 -> h2o2 + co2

        !   jpl 2003

        i0 = i0 + 1

        do ilay = 1,nlay
            ak0 = 2.5*6.9e-31*(t(ilay)/300.)**(-1.0)
            ak1 = 2.6e-11*(t(ilay)/300.)**(0.0)
   
            rate = (ak0*dens(ilay))/(1. + ak0*dens(ilay)/ak1)
            xpo = 1./(1. + dlog10((ak0*dens(ilay))/ak1)**2)
            c017(ilay) = rate*0.6**xpo
         end do

        rrates(:,i0) = c017(:)
        rtype(i0) = 2

        ns(i0) = 1
        sID(1,i0) = 13
        sISO(1,i0) = 0
        sf(1,i0) = 2.0

        npr(i0) = 1
        pID(1,i0) = 25
        pISO(1,i0) = 0
        pf(1,i0) = 1.0


        !=========================================================================================
        !---c018: h + h + co2 -> h2 + co2

        !   baulch et al., 2005

        i0 = i0 + 1

        c018(:) = 2.5*1.8e-30*(t(:)**(-1.0))*dens(:)

        rrates(:,i0) = c018(:)
        rtype(i0) = 2

        ns(i0) = 1
        sID(1,i0) = 48
        sISO(1,i0) = 0
        sf(1,i0) = 2.0

        npr(i0) = 1
        pID(1,i0) = 39
        pISO(1,i0) = 0
        pf(1,i0) = 1.0


        !==========================================================================================
        !NITROGEN REACTIONS
        !==========================================================================================

        if(nitrogen_chemistry)then

            !=========================================================================================
            !---  d001: no2 + o -> no + o2

            !     jpl 2006


            i0 = i0 + 1

            d001(:) = 5.1e-12*exp(210./t(:))

            rrates(:,i0) = d001(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 10
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 45
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 8
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(2,i0) = 7
            pISO(2,i0) = 0
            pf(2,i0) = 1.0

            !=========================================================================================
            !---  d002: no + o3 -> no2 + o2

            !     jpl 2006


            i0 = i0 + 1

            d002(:) = 3.0e-12*exp(-1500./t(:))

            rrates(:,i0) = d002(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 8
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 3
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 10
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(2,i0) = 7
            pISO(2,i0) = 0
            pf(2,i0) = 1.0


            !=========================================================================================
            !---  d003: no + ho2 -> no2 + oh

            !     jpl 2011


            i0 = i0 + 1

            d003(:) = 3.3e-12*exp(270./t(:))

            rrates(:,i0) = d003(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 8
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 44
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 10
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(2,i0) = 13
            pISO(2,i0) = 0
            pf(2,i0) = 1.0


            !=========================================================================================
            !---  d004: n + no -> n2 + o

            !     jpl 2011


            i0 = i0 + 1

            d004(:) = 2.1e-11*exp(100./t(:))

            rrates(:,i0) = d004(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 134
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 8
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 22
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(2,i0) = 45
            pISO(2,i0) = 0
            pf(2,i0) = 1.0


            !=========================================================================================
            !---  d005: n + o2 -> no + o

            !     jpl 2011

            i0 = i0 + 1

            d005(:) = 1.5e-11*exp(-3600./t(:))

            rrates(:,i0) = d005(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 134
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 7
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 8
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(2,i0) = 45
            pISO(2,i0) = 0
            pf(2,i0) = 1.0


            !=========================================================================================
            !---  d006: no2 + h -> no + oh

            !     jpl 2011


            i0 = i0 + 1

            d006(:) = 4.0e-10*exp(-340./t(:))

            rrates(:,i0) = d006(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 10
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 48
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 8
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(2,i0) = 13
            pISO(2,i0) = 0
            pf(2,i0) = 1.0



            !=========================================================================================
            !---  d007: n + o -> no

            i0 = i0 + 1

            d007(:) = 2.8e-17*(300./t(:))**0.5

            rrates(:,i0) = d007(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 134
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 45
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 1
            pID(1,i0) = 8
            pISO(1,i0) = 0
            pf(1,i0) = 1.0


            !=========================================================================================
            !---  d008: n + ho2 -> no + oh

            !     brune et al., j. chem. phys., 87, 1983    


            i0 = i0 + 1

            d008(:) = 2.19e-11

            rrates(:,i0) = d008(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 134
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 44
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 8
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(1,i0) = 13
            pISO(1,i0) = 0
            pf(1,i0) = 1.0

            !=========================================================================================
            !---  d009: n + oh -> no + h

            !     atkinson et al., j. phys. chem. ref. data, 18, 881, 1989


            i0 = i0 + 1

            d009(:) = 3.8e-11*exp(85./t(:))

            rrates(:,i0) = d009(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 134
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 13
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 8
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(1,i0) = 48
            pISO(1,i0) = 0
            pf(1,i0) = 1.0


            !=========================================================================================
            !---  d010: n(2d) + o  -> n + o

            !     herron, j. phys. chem. ref. data, 1999


            i0 = i0 + 1

            d010(:) = 3.3e-12*exp(-260./t(:))

            rrates(:,i0) = d010(:)*o(:)
            rtype(i0) = 1
    
            ns(i0) = 1
            sID(1,i0) = 135
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
    
            npr(i0) = 1
            pID(1,i0) = 134
            pISO(1,i0) = 0
            pf(1,i0) = 1.0

            !=========================================================================================
            !---  d011: n(2d) + n2  -> n + n2

            !     herron, j. phys. chem. ref. data, 1999

            i0 = i0 + 1

            d011(:) = 1.7e-14

            rrates(:,i0) = d011(:)*n2(:)
            rtype(i0) = 1
    
            ns(i0) = 1
            sID(1,i0) = 135
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
    
            npr(i0) = 1
            pID(1,i0) = 134
            pISO(1,i0) = 0
            pf(1,i0) = 1.0


            !=========================================================================================
            !---  d012: n(2d) + co2  -> no + co

            !     herron, j. phys. chem. ref. data, 1999


            i0 = i0 + 1

            d012(:) = 3.6e-13

            rrates(:,i0) = d012(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 135
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 2
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 8
            pISO(1,i0) = 0
            pf(1,i0) = 1.0
            pID(1,i0) = 5
            pISO(1,i0) = 0
            pf(1,i0) = 1.0

        endif
        

        !==========================================================================================
        !CARBON REACTIONS
        !==========================================================================================

        !=========================================================================================
        !---  e001: oh + co -> co2 + h

        !     jpl 2003

        !     e001(:) = 1.5e-13*(1 + 0.6*press(:)/1013.)

        !     mccabe et al., grl, 28, 3135, 2001

        !     e001(:) = 1.57e-13 + 3.54e-33*dens(:)

        !     jpl 2015

        !do ilay = 1,nlay

            !branch 1 : oh + co -> h + co2
            
        !    rate1 = 1.5e-13*(t(ilay)/300.)**(0.0)
            
            !branch 2 : oh + co + m -> hoco + m
            
        !    ak0 = 5.9e-33*(t(ilay)/300.)**(-1.0)
        !    ak1 = 1.1e-12*(t(ilay)/300.)**(1.3)
        !    rate2 = (ak0*dens(ilay))/(1. + ak0*dens(ilay)/ak1)
        !    xpo = 1./(1. + dlog10((ak0*dens(ilay))/ak1)**2)
            
        !    e001(ilay) = rate1 + rate2*0.6**xpo
        !end do

        i0 = i0 + 1 

        !     jpl 2019

        !     For the sake of simplicity, it is assumed that the association yield (kf)
        !     gives the same product as the chemical activation yield (kca).
        !     Thus the only products are h + co2. There is no production of hoco.

        do ilay=1,nlay

            !association
            
            k0 = 2.5*6.9e-33*(298./t(ilay))**(2.1)
            kinf = 1.1e-12*(298./t(ilay))**(-1.3)
            
            kf = (kinf*k0*dens(ilay)/(kinf + k0*dens(ilay)))  &
                *0.6**(1. + (log10(k0*dens(ilay)/kinf))**2.)**(-1.0)
            
            !chemical activation
            
            kint = 1.85e-13*exp(-65./t(ilay))
            
            kca = kint*(1. - kf/kinf)
            
            !total : association + chemical activation
            
            e001(ilay) = kf + kca
            
        end do

        if(o18_chemistry)then
            !Tweaking the model to account for constant (18O)/(16O) ratio in OH
            !This reaction is for (16O)H
            rrates(:,i0) = e001(:)/(1.d0+o18ratio_oh*2005.20d-6)

        else
            rrates(:,i0) = e001(:)
        endif
        
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 13
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 5
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 2
        pID(1,i0) = 2
        pISO(1,i0) = 0
        pf(1,i0) = 1.0
        pID(2,i0) = 48
        pISO(2,i0) = 0
        pf(2,i0) = 1.0

        !=========================================================================================
        !---  e002: o + co + m -> co2 + m

        !     tsang and hampson, 1986.

        i0 = i0 + 1

        e002(:) = 2.5*6.5e-33*exp(-2184./t(:))*dens(:)

        rrates(:,i0) = e002(:)
        rtype(i0) = 3

        ns(i0) = 2
        sID(1,i0) = 45
        sISO(1,i0) = 0
        sf(1,i0) = 1.0
        sID(2,i0) = 5
        sISO(2,i0) = 0
        sf(2,i0) = 1.0

        npr(i0) = 1
        pID(1,i0) = 2
        pISO(1,i0) = 0
        pf(1,i0) = 1.0

        !=========================================================================================

        !==========================================================================================
        !CARBON-13 REACTIONS
        !==========================================================================================

        if(c13_chemistry)then

            !=========================================================================================
            !---g001: oh + (13c)o -> (13c)o2 + h


            !Assumed to be the same as e001, but with a 
            !fractionation factor from Stevens et al. (1980)
    
            !A polynomial function is fit to capture the pressure-
            !dependence of the fractionation
    
            !k13 / k12 = 1.00638 - 1.693e-5*press(hPa) + 4.6968e-9 * press(hPa)**2.
            !k13 / k12 = 1.00638 - 1.693e-7*press(Pa) + 4.6968e-13 * press(Pa)**2.

            i0 = i0 + 1

            g001(:) = e001(:)*(1.00638 - 1.693e-7*P(:) + 4.6968e-13*P(:)**2.)

            rrates(:,i0) = g001(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 13
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 5
            sISO(2,i0) = 2
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 2
            pISO(1,i0) = 2
            pf(1,i0) = 1.0
            pID(2,i0) = 48
            pISO(2,i0) = 0
            pf(2,i0) = 1.0


            !=========================================================================================
            !---  g002: o + (13c)o + m -> (13c)o2 + m

            !     we assume it is the same rate as e002

            i0 = i0 + 1
    
            rrates(:,i0) = e002(:) * 1.0
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 45
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 5
            sISO(2,i0) = 2
            sf(2,i0) = 1.0
    
            npr(i0) = 1
            pID(1,i0) = 2
            pISO(1,i0) = 2
            pf(1,i0) = 1.0


            if(nitrogen_chemistry)then
                !=========================================================================================
                !---  g003: n(2d) + (13c)o2  -> no + (13c)o

                !     we assume it is the same rate as d012


                i0 = i0 + 1
    
                rrates(:,i0) = d012(:)
                rtype(i0) = 3
        
                ns(i0) = 2
                sID(1,i0) = 135
                sISO(1,i0) = 0
                sf(1,i0) = 1.0
                sID(2,i0) = 2
                sISO(2,i0) = 2
                sf(2,i0) = 1.0
        
                npr(i0) = 2
                pID(1,i0) = 8
                pISO(1,i0) = 0
                pf(1,i0) = 1.0
                pID(1,i0) = 5
                pISO(1,i0) = 2
                pf(1,i0) = 1.0

            endif


        endif


        !=========================================================================================

        !==========================================================================================
        !OXYGEN-18 REACTIONS
        !==========================================================================================

        if(o18_chemistry)then

            !=========================================================================================
            !---h001: oh + c(18o) -> (18o)(12c)(16o) + h


            !Assumed to be the same as e001, but with a 
            !fractionation factor from Stevens et al. (1980)
    
            !A polynomial function is fit to capture the pressure-
            !dependence of the fractionation
    
            !k18 / k16 = 1.01195 - 4.443404e-8*press(Pa) + 1.938355e-13 * press(Pa)**2.

            i0 = i0 + 1

            h001(:) = e001(:)*(1.01195 - 4.443404e-8*P(:) + 1.938355e-13*P(:)**2.)

            !rrates(:,i0) = h001(:)
            !We assume that the (18O)/(16O) ratio is constant and we adapt the reaction rate for it
            rrates(:,i0) = h001(:)/(1.d0+o18ratio_oh*2005.20d-6)

            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 13
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 5
            sISO(2,i0) = 3
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 2
            pISO(1,i0) = 3
            pf(1,i0) = 1.0
            pID(2,i0) = 48
            pISO(2,i0) = 0
            pf(2,i0) = 1.0


            !=========================================================================================
            !---  h002: o + (12c)(18o) + m -> (18o)(12c)(16o) + m

            !     we assume it is the same rate as e002

            i0 = i0 + 1
    
            rrates(:,i0) = e002(:) * 1.0
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 45
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 5
            sISO(2,i0) = 3
            sf(2,i0) = 1.0
    
            npr(i0) = 1
            pID(1,i0) = 2
            pISO(1,i0) = 3
            pf(1,i0) = 1.0


            !=========================================================================================
            !---h001_2: (18o)h + c(16o) -> (18o)(12c)(16o) + h

            !Assumed to be the same as e001, but tweaking the model to compute the rate from (16O)H.
            !The tweak is to assume the rate is lower by a factor of (18O)/(16O), which represents the
            !O isotopic ratio in OH (here we assume it follows the value from Curiosity in CO2 = 1.046 VSMOW)

            i0 = i0 + 1

            h001_2(:) = e001(:)*(o18ratio_oh*2005.20d-6)/(1.d0+o18ratio_oh*2005.20d-6)

            rrates(:,i0) = h001_2(:)
            rtype(i0) = 3
    
            ns(i0) = 2
            sID(1,i0) = 13
            sISO(1,i0) = 0
            sf(1,i0) = 1.0
            sID(2,i0) = 5
            sISO(2,i0) = 0
            sf(2,i0) = 1.0
    
            npr(i0) = 2
            pID(1,i0) = 2
            pISO(1,i0) = 3
            pf(1,i0) = 1.0
            pID(2,i0) = 48
            pISO(2,i0) = 0
            pf(2,i0) = 1.0


            if(nitrogen_chemistry)then
                !=========================================================================================
                !---  h003: n(2d) + (18o)(12c)(16o)  -> no + (12c)(18o)

                !     we assume it is the same rate as d012


                i0 = i0 + 1
    
                rrates(:,i0) = d012(:)
                rtype(i0) = 3
        
                ns(i0) = 2
                sID(1,i0) = 135
                sISO(1,i0) = 0
                sf(1,i0) = 1.0
                sID(2,i0) = 2
                sISO(2,i0) = 3
                sf(2,i0) = 1.0
        
                npr(i0) = 2
                pID(1,i0) = 8
                pISO(1,i0) = 0
                pf(1,i0) = 1.0
                pID(1,i0) = 5
                pISO(1,i0) = 3
                pf(1,i0) = 1.0

            endif


        endif


        !=========================================================================================


    end subroutine

!============================================================================================================================

    subroutine load_reactions_photolysis(ngas,nitrogen_chemistry,c13_chemistry,o18_chemistry,&
        ngas_phot,gasID_phot,isoID_phot)

        !Loading the photolysis reactions relevant for the atmosphere of Mars

        !Inputs
        integer, intent(in) :: ngas
        logical, intent(in) :: nitrogen_chemistry, c13_chemistry, o18_chemistry

        !local
        integer :: ilast

        !Outputs
        integer, intent(out) :: ngas_phot    !Number of gases to be photolysed
        integer, intent(out) :: gasID_phot(ngas),isoID_phot(ngas)  !Gases to be photolysed


        ngas_phot = 7 

        gasID_phot(1) = 7   !O2
        isoID_phot(1) = 0   

        gasID_phot(2) = 2   !CO2
        isoID_phot(2) = 0  

        gasID_phot(3) = 3   !O3
        isoID_phot(3) = 0  

        gasID_phot(4) = 1   !H2O
        isoID_phot(4) = 0  

        gasID_phot(5) = 25  !H2O2
        isoID_phot(5) = 0  

        gasID_phot(6) = 44  !HO2
        isoID_phot(6) = 0  

        gasID_phot(7) = 39  !H2
        isoID_phot(7) = 0  

        ilast = 7

        if(nitrogen_chemistry)then

            ngas_phot = ngas_phot + 3

            ilast = ilast + 1
            gasID_phot(ilast) = 8  !NO
            isoID_phot(ilast) = 0

            ilast = ilast + 1
            gasID_phot(ilast) = 10  !NO2
            isoID_phot(ilast) = 0

            ilast = ilast + 1
            gasID_phot(ilast) = 22  !N2
            isoID_phot(ilast) = 0

        endif


        if(c13_chemistry)then

            ngas_phot = ngas_phot + 1

            ilast = ilast + 1
            gasID_phot(ilast) = 2  !(13C)O2
            isoID_phot(ilast) = 2
        endif

        if(o18_chemistry)then

            ngas_phot = ngas_phot + 1

            ilast = ilast + 1
            gasID_phot(ilast) = 2  !(18O)(12C)(16O)
            isoID_phot(ilast) = 3
        endif

    end subroutine


END MODULE mars_chemistry
