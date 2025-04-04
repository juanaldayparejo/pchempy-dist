MODULE photolysis


    !PHOTOLYSIS REACTIONS INCLUDED IN THE SCHEME
    !##################################################################################

    !Defining the gases whose photolysis can be determined in the current setup
    integer, save :: n_gas_phot_inc = 19

    integer, parameter, dimension(19) :: gasID_phot_inc = (/&
    7, 2, 3, 1, 25, 44, 39, 8, 10, 22, & !O2,CO2,O3,H2O,H2O2,HO2,H2,NO,NO2,N2
    2, &                                 !13C - CO2
    7, 2, 3, 1, 25, 44, 8, 10/)          !18O - O2,CO2,O3,H2O,H2O2,HO2,NO,NO2

    integer, parameter, dimension(19) :: isoID_phot_inc = (/&
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &     !O2,CO2,O3,H2O,H2O2,HO2,H2,NO,NO2,N2
    2, &                                !13C - CO2
    2, 3, 2, 2, 2, 2, 3, 3/)            !18O - O2,CO2,O3,H2O,H2O2,HO2,NO,NO2

    !Number of reactions associated with each photolysis
    integer, parameter, dimension(19) :: n_phot_inc = (/&
    2, 2, 2, 1, 1, 1, 1, 1, 1, 1, &    !O2,CO2,O3,H2O,H2O2,HO2,H2,NO,NO2,N2,
    2, &                               !13C - CO2
    2, 2, 2, 1, 1, 1, 1, 1/)           !18O - O2,CO2,O3,H2O,H2O2,HO2,NO,NO2

    CONTAINS

!==============================================================================================================================================

    subroutine photolysis_online(nlay, ngas, ngas_phot, gasID, isoID, gasID_phot, isoID_phot, &
        h, T, cdens, sza, tau, dist_sun, n_phot, &
        rtype, ns, sID, sISO, sf, npr, pID, pISO, pf, rrates)

        ! Function to calculate the photolysis rates of a given set of reactions

        use photolysis_mod

        implicit none

        ! input
        integer, intent(in) :: nlay                         ! number of layers in atmosphere
        integer, intent(in) :: ngas                         ! number of species in the chemistry
        integer, intent(in) :: gasID(ngas),isoID(ngas)      ! ID of the gases in the atmosphere
        integer, intent(in) :: ngas_phot                    ! number of active species
        integer, intent(in) :: gasID_phot(ngas_phot),isoID_phot(ngas_phot)   ! ID of the active gases 
        real, intent(in) :: h(nlay)                         ! Altitude (m)
        real, intent(in) :: T(nlay)                         ! Temperature (K)
        double precision, intent(in) :: cdens(nlay,ngas)    ! column density of each species (m-2)
        real, intent(in) :: sza,tau,dist_sun                ! Solar zenith angle (degrees), dust opacity, Sun-planet distance (AU)
        integer, intent(in) :: n_phot                       ! number of photolysis reactions (calculated before based on particular chemistry)


        ! local
        integer :: ID_pos(ngas_phot)
        integer :: i,igas,igasx,n_phot1,i0,ilay,iw
        integer :: j_o2_o, j_o2_o1d
        integer :: j_18o2_o, j_18o2_o1d_1, j_18o2_o1d_2
        integer :: j_co2_o, j_co2_o1d
        integer :: j_13co2_o, j_13co2_o1d
        integer :: j_18co2_o_1, j_18co2_o_2, j_18co2_o1d_1, j_18co2_o1d_2
        integer :: j_o3_o, j_o3_o1d
        integer :: j_18o3_o_1, j_18o3_o_2, j_18o3_o1d_1, j_18o3_o1d_2
        integer :: j_h2o,j_h2o2,j_ho2,j_h2
        integer :: j_18h2o,j_18h2o2,j_18ho2_1,j_18ho2_2
        integer :: j_no,j_no2,j_n2
        integer :: j_18no,j_18no2_1,j_18no2_2
        integer :: mopt
        logical :: xsloaded
        double precision :: colincg(nlay,ngas),vmr(nlay,ngas),colinc(nlay)
        real, dimension(nlay,nw,ngas_phot) :: dtgas
        real, dimension(nlay,nw,n_phot) :: sj
        real, dimension(nlay,nw) :: dagas
        real, dimension(nlay) :: alt
        real, dimension(nlay,nw) :: dtaer       
        real, dimension(nlay,nw) :: omaer                        
        real, dimension(nlay,nw) :: gaer
        real, dimension(nlay,nw) :: dtcld                   
        real, dimension(nlay,nw) :: omcld                
        real, dimension(nlay,nw) :: gcld  
        real, dimension(nlay,nw) :: dtrl
        integer, dimension(0:nlay) :: nid
        real, dimension(0:nlay,nlay) :: dsdh
        real, dimension(nw) :: fmars 
        real :: factor
        real, dimension(nlay) :: edir, edn, eup        
        real, dimension(nlay) :: fdir, fdn, fup     
        real, dimension(nlay) :: saflux
        real :: deltaj  

        ! output
        integer, intent(out) :: rtype(n_phot)                ! reaction type - always 1 for photolysis
        integer, intent(out) :: ns(n_phot),npr(n_phot)       ! number of source and products
        integer, intent(out) :: sID(2,n_phot),sISO(2,n_phot) ! ID of the source molecules
        integer, intent(out) :: pID(2,n_phot),pISO(2,n_phot) ! ID of the product molecules
        real, intent(out) :: sf(2,n_phot),pf(2,n_phot)       ! number of molecules per source/product
        double precision, intent(out) :: rrates(nlay,n_phot) ! reaction rates (s-1)



        !reading all the cross sections for the different species
        !####################################################################

        mopt = 1

        xsloaded = allocated(wc)   !Flag to see if arrays have been already allocated (i.e. read)
        if(xsloaded .eqv. .false.) call init_photolysis(ngas_phot,gasID_phot,isoID_phot,mopt)

        !Converting the column densities to cm-2 and calculating the mixing ratios
        !####################################################################

        !Converting the column density to cm-2
        colincg(:,:) = cdens(:,:) * 1.0d-4 

        ! Calculating the total column density (cm-2) and the mixing ratios
        do ilay=1,nlay
            colinc(ilay) = 0.d0
            do igas=1,ngas
              colinc(ilay) = colinc(ilay) + colincg(ilay,igas)
            enddo

            do igas=1,ngas
                vmr(ilay,igas) = colincg(ilay,igas) / colinc(ilay)
            enddo

        enddo

        !Converting the altitudes to km
        alt(:) = h(:) * 1.0e-3


        !Locating the position of the active gases 
        !################################################

        do i=1,ngas_phot
            igasx = 0
            do igas=1,ngas
                if((gasID(igas).eq.gasID_phot(i)).and.(isoID(igas).eq.isoID_phot(i)))then
                    ID_pos(i) = igas
                    igasx = 1
                endif
            enddo
            if(igasx.eq.0)then
                print*,'error, photolysis involves a gas not present in the atmosphere'
                print*,'gas',gasID_phot(i),isoID_phot(i)
                stop
            endif
        enddo


        !Cheking that all the active gases have an associated photolysis implemented
        !#####################################################################################

        do i=1,ngas_phot
            igasx = 0
            do igas=1,n_gas_phot_inc
                if((gasID_phot_inc(igas).eq.gasID_phot(i)).and.(isoID_phot_inc(igas).eq.isoID_phot(i)))then
                    igasx = 1
                endif
            enddo
            if(igasx.eq.0)then
                print*,'error, the photolysis of a certain gas has not been included yet'
                print*,'gas',gasID_phot(i),isoID_phot(i)
                stop
            endif
        enddo


        !Defining the number of reactions and the reactions types for each photolysis
        !and calculating the optical depth for each gas
        !######################################################################################

        !initialising parameters
        rtype(1:n_phot) = 0
        ns(:) = 0
        sID(:,:) = 0
        sISO(:,:) = 0
        sf(:,:) = 0.0
        npr(:) = 0
        pID(:,:) = 0
        pISO(:,:) = 0
        pf(:,:) = 0.0

        n_phot1 = 0
        i0 = 0
        do i=1,ngas_phot

            !O2 PHOTOLYSIS
            if( (gasID_phot(i).eq.7).and.(isoID_phot(i).eq.0) )then     !O2  

                n_phot1 = n_phot1 + 2

                !O2 + hv --> O + O
                i0 = i0 + 1
                j_o2_o = i0

                !O2 + hv --> O + O(1D)
                i0 = i0 + 1
                j_o2_o1d = i0

                call seto2(n_phot, nlay, nw, wc, mopt, T, xso2_150, xso2_200,& 
                           xso2_250, xso2_300, yieldo2, j_o2_o, j_o2_o1d, &
                           colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                           rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !(18O)(16O) PHOTOLYSIS
            elseif( (gasID_phot(i).eq.7).and.(isoID_phot(i).eq.2) )then     !(18O)(16O)  

                n_phot1 = n_phot1 + 3

                !(18O)(16O) + hv --> (18O) + O
                i0 = i0 + 1
                j_18o2_o = i0

                !(18O)(16O) + hv --> (18O) + O(1D)
                i0 = i0 + 1
                j_18o2_o1d_1 = i0

                !(18O)(16O) + hv --> O + 18O(1D)
                i0 = i0 + 1
                j_18o2_o1d_2 = i0

                call set18o2(n_phot, nlay, nw, wc, mopt, T, &
                           xs18o2_150, xs18o2_200, xs18o2_250, xs18o2_300, yieldo2, &
                           j_18o2_o, j_18o2_o1d_1, j_18o2_o1d_2, &
                           colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                           rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !CO2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.2).and.(isoID_phot(i).eq.0) )then  !CO2

                n_phot1 = n_phot1 + 2

                !CO2 + hv --> CO + O
                i0 = i0 + 1
                j_co2_o = i0

                !CO2 + hv --> CO + O(1D)
                i0 = i0 + 1
                j_co2_o1d = i0

                call setco2(n_phot, nlay, nw, wc, T, xsco2_195, xsco2_295, &
                            xsco2_370, yieldco2, j_co2_o, j_co2_o1d, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj,&
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !(13C)O2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.2).and.(isoID_phot(i).eq.2) )then  !(13C)O2

                n_phot1 = n_phot1 + 2

                !(13C)O2 + hv --> (13C)O + O
                i0 = i0 + 1
                j_13co2_o = i0

                !(13C)O2 + hv --> (13C)O + O(1D)
                i0 = i0 + 1
                j_13co2_o1d = i0

                call set13co2(n_phot, nlay, nw, wc, T, xs13co2_195, xs13co2_295, &
                            xs13co2_370, yieldco2, j_13co2_o, j_13co2_o1d, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj,&
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !(18O)(12C)(16O) PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.2).and.(isoID_phot(i).eq.3) )then  !(18O)(12C)(16O)

                n_phot1 = n_phot1 + 4

                !(18O)(12C)(16O) + hv --> (12C)(18O) + O
                i0 = i0 + 1
                j_18co2_o_1 = i0

                !(12C)(16O)2 + hv --> (12C)(16O) + 18O
                i0 = i0 + 1
                j_18co2_o_2 = i0

                !(18O)(12C)(16O) + hv --> (12C)(18O) + O(1D)
                i0 = i0 + 1
                j_18co2_o1d_1 = i0

                !(18O)(12C)(16O) + hv --> (12C)(16O) + 18O(1D)
                i0 = i0 + 1
                j_18co2_o1d_2 = i0

                call set18co2(n_phot, nlay, nw, wc, T, xs18co2_195, xs18co2_295, &
                            xs18co2_370, yieldco2, j_18co2_o_1, j_18co2_o_2, j_18co2_o1d_1, j_18co2_o1d_2, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj,&
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !O3 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.3).and.(isoID_phot(i).eq.0) )then  !O3

                n_phot1 = n_phot1 + 2

                !O3 + hv --> O2 + O
                i0 = i0 + 1
                j_o3_o = i0

                !O3 + hv --> O2 + O(1D)
                i0 = i0 + 1
                j_o3_o1d = i0

                call seto3(n_phot, nlay, nw, wc, T, xso3_218, xso3_298, &
                            j_o3_o, j_o3_o1d, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj,&
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)


            !(18O)O2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.3).and.(isoID_phot(i).eq.2) )then  !(18O)O2

                n_phot1 = n_phot1 + 4

                !(18O)O2 + hv --> (18O)O + O
                i0 = i0 + 1
                j_18o3_o_1 = i0

                !(18O)O2 + hv --> O2 + (18O)
                i0 = i0 + 1
                j_18o3_o_2 = i0

                !(18O)O2 + hv --> (18O)O + O(1D)
                i0 = i0 + 1
                j_18o3_o1d_1 = i0

                !(18O)O2 + hv --> O2 + 18O(1D)
                i0 = i0 + 1
                j_18o3_o1d_2 = i0

                call set18o3(n_phot, nlay, nw, wc, T, xso3_218, xso3_298, &
                            j_18o3_o_1, j_18o3_o_2, j_18o3_o1d_1, j_18o3_o1d_2, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj,&
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)


            !H2O PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.1).and.(isoID_phot(i).eq.0) )then  !H2O

                n_phot1 = n_phot1 + 1

                !H2O + hv --> OH + H
                i0 = i0 + 1
                j_h2o = i0

                call seth2o(n_phot, nlay, nw, wc, T, xsh2o, j_h2o, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)


            !H2(18O) PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.1).and.(isoID_phot(i).eq.2) )then  !H2(18O)

                n_phot1 = n_phot1 + 1

                !H2(18O) + hv --> (18O)H + H
                i0 = i0 + 1
                j_18h2o = i0

                call set18h2o(n_phot, nlay, nw, wc, T, xs18h2o, j_18h2o, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)


            !H2O2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.25).and.(isoID_phot(i).eq.0) )then  !H2O2

                n_phot1 = n_phot1 + 1

                !H2O2 + hv --> OH + OH
                i0 = i0 + 1
                j_h2o2 = i0

                call seth2o2(n_phot, nlay, nw, wc, T, xsh2o2, j_h2o2,&
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !H2(18O)O PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.25).and.(isoID_phot(i).eq.2) )then  !H2(18O)O

                n_phot1 = n_phot1 + 1

                !H2(18O)O + hv --> (18O)H + OH
                i0 = i0 + 1
                j_18h2o2 = i0

                call set18h2o2(n_phot, nlay, nw, wc, T, xs18h2o2, j_18h2o2,&
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)


            !HO2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.44).and.(isoID_phot(i).eq.0) )then  !HO2

                n_phot1 = n_phot1 + 1

                !HO2 + hv --> OH + O
                i0 = i0 + 1
                j_ho2 = i0

                call setho2(n_phot, nlay, nw, wc, T, xsho2, j_ho2, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)


            !H(18O)O PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.44).and.(isoID_phot(i).eq.2) )then  !H(18O)O

                n_phot1 = n_phot1 + 2

                !H(18O)O + hv --> (18O)H + O
                i0 = i0 + 1
                j_18ho2_1 = i0

                !H(18O)O + hv --> OH + (18O)
                i0 = i0 + 1
                j_18ho2_2 = i0

                call set18ho2(n_phot, nlay, nw, wc, T, xs18ho2, j_18ho2_1, j_18ho2_2, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)


            !H2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.39).and.(isoID_phot(i).eq.0) )then  !H2

                n_phot1 = n_phot1 + 1

                !H2 + hv --> H + H
                i0 = i0 + 1
                j_h2 = i0

                call seth2(n_phot, nlay, nw, wc, T, xsh2, yieldh2, j_h2, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !NO PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.8).and.(isoID_phot(i).eq.0) )then  !NO

                n_phot1 = n_phot1 + 1

                !NO + hv --> N + O
                i0 = i0 + 1
                j_no = i0

                call setno(n_phot, nlay, nw, wc, T, xsno, yieldno, j_no, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !N(18O) PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.8).and.(isoID_phot(i).eq.3) )then  !N(18O)

                n_phot1 = n_phot1 + 1

                !N(18O) + hv --> N + (18O)
                i0 = i0 + 1
                j_18no = i0

                call set18no(n_phot, nlay, nw, wc, T, xs18no, yieldno, j_18no, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !NO2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.10).and.(isoID_phot(i).eq.0) )then  !NO2

                n_phot1 = n_phot1 + 1

                !NO2 + hv --> NO + O
                i0 = i0 + 1
                j_no2 = i0

                call setno2(n_phot, nlay, nw, wc, T, xsno2, xsno2_220,&
                            xsno2_294, yldno2_248, yldno2_298, j_no2,&
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj,&
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !N(18O)(16O) PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.10).and.(isoID_phot(i).eq.3) )then  !NO(18O)

                n_phot1 = n_phot1 + 2

                !NO(18O) + hv --> N(18O) + O
                i0 = i0 + 1
                j_18no2_1 = i0

                !NO(18O) + hv --> NO + (18O)
                i0 = i0 + 1
                j_18no2_2 = i0

                call set18no2(n_phot, nlay, nw, wc, T, xs18no2, xs18no2_220,&
                            xs18no2_294, yldno2_248, yldno2_298, &
                            j_18no2_1, j_18no2_2,&
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj,&
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            !N2 PHOTOLYSIS
            elseif ( (gasID_phot(i).eq.22).and.(isoID_phot(i).eq.0) )then  !N2

                n_phot1 = n_phot1 + 1

                !N2 + hv --> N + N
                i0 = i0 + 1
                j_n2 = i0

                call setn2(n_phot, nlay, nw, wc, T, xsn2, yieldn2, j_n2, &
                            colinc, vmr(:,ID_pos(i)), dtgas(:,:,i), sj, &
                            rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

            endif

        enddo

        if(n_phot1.ne.n_phot)then
            print*,'error in photolysis_online. n_phot1 is not equal to n_phot'
            print*,'n_phot1 = ',n_phot1
            print*,'n_phot = ',n_phot
            stop
        endif

        !Calculating the total gas optical depth
        !######################################################################################

        dagas(:,:) = 0.

        do ilay = 1,nlay
           do iw = 1,nw-1
              do i = 1,ngas_phot
                 dagas(ilay,iw) = dagas(ilay,iw) + dtgas(ilay,iw,i)
              end do
           end do
        end do

        !Set rayleigh optical depth (CO2 dominated atmosphere)
        !######################################################################################

        call setairlay(nlay, nw, wl, wc, colinc, dtrl)

        !Set aerosol properties and optical depth
        !######################################################################################

        call setaer(nlay,alt,tau,nw,dtaer,omaer,gaer)

        !Set cloud properties and optical depth
        !######################################################################################

        call setcld(nlay,nw,dtcld,omcld,gcld)

        !Slant path lengths in spherical geometry
        !######################################################################################

        call sphers(nlay,alt,sza,dsdh,nid)

        !solar flux at planet
        !######################################################################################

        factor = (1./dist_sun)**2.
        do iw = 1,nw-1
           fmars(iw) = f(iw)*factor
        end do

        !Calculate the photolysis rates for each reaction
        !######################################################################################

        rrates(:,1:n_phot) = 0.d0

        do iw = 1,nw-1

            !monochromatic radiative transfer. outputs are:
            !normalized irradiances     edir(nlayer), edn(nlayer), eup(nlayer) 
            !normalized actinic fluxes  fdir(nlayer), fdn(nlayer), fup(nlayer)
            !where 
            !dir = direct beam, dn = down-welling diffuse, up = up-welling diffuse

            call rtlink(nlay, nw, iw, albedo(iw), sza, dsdh, nid, dtrl,&
                       dagas, dtcld, omcld, gcld, dtaer, omaer, gaer,&
                       edir, edn, eup, fdir, fdn, fup)


            !spherical actinic flux

            do ilay = 1,nlay
                saflux(ilay) = fmars(iw)*(fdir(ilay) + fdn(ilay) + fup(ilay))
            end do

            !photolysis rate integration

            do i = 1,n_phot
                do ilay = 1,nlay
                    deltaj = saflux(ilay)*sj(ilay,iw,i)
                    rrates(ilay,i) = rrates(ilay,i) + deltaj*(wu(iw)-wl(iw))
                end do
            end do

        end do ! iw

        !eliminate small values
        where (rrates(:,1:n_phot) < 1.e-30)
            rrates(:,1:n_phot) = 0.d0
        end where

    end subroutine


!==============================================================================================================================================

    subroutine seto2(n_phot, nlayer, nw, wc, mopt, tlay, xso2_150, &
                      xso2_200, xso2_250, xso2_300, yieldo2, j_o2_o, &
                      j_o2_o1d, colinc, rm, dt, sj, &
                      rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the O2 temperature-dependent cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: mopt                                          
        integer :: j_o2_o, j_o2_o1d                             
        real, dimension(nw)     :: wc                            
        real, dimension(nw)     :: xso2_150,xso2_200,xso2_250,xso2_300
        real, dimension(nw)     :: yieldo2                       
        real, dimension(nlayer) :: tlay             !temperature             
        double precision, dimension(nlayer) :: rm               !mixing ratio of O2      
        double precision, dimension(nlayer) :: colinc           !column density (cm-2)              

        !     output:

        real, dimension(nlayer,nw)    :: dt         !optical depth             
        real, dimension(nlayer,nw,n_phot) :: sj     !cross sections * yield                

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product


        !     local:

        integer :: ilev, iw
        real    :: temp
        real    :: xso2, factor

        !     correction by factor if low-resolution in schumann-runge bands

        if (mopt == 1) then
            factor = 1.
        else if (mopt == 2) then
            factor = 0.8
        end if

        !     calculate temperature dependance

        do ilev = 1,nlayer
            temp = max(tlay(ilev),150.)
            temp = min(temp, 300.)

            do iw = 1, nw-1
            if (tlay(ilev) > 250.) then
                xso2 = xso2_250(iw) + (xso2_300(iw) - xso2_250(iw)) &
                                    /(300. - 250.)*(temp - 250.) 
            else if (tlay(ilev) > 200.) then
                xso2 = xso2_200(iw) + (xso2_250(iw) - xso2_200(iw)) & 
                                    /(250. - 200.)*(temp - 200.) 
            else
                xso2 = xso2_150(iw) + (xso2_200(iw) - xso2_150(iw)) &
                                    /(200. - 150.)*(temp - 150.) 
            end if

            if (wc(iw) > 180. .and. wc(iw) < 200.) then
                xso2 = xso2*factor
            end if

        !     optical depth

            dt(ilev,iw) = colinc(ilev)*rm(ilev)*xso2

        !     production of o(1d) for wavelengths shorter than 175 nm

            if (wc(iw) >= 175.) then
                sj(ilev,iw,j_o2_o)   = xso2*yieldo2(iw)
                sj(ilev,iw,j_o2_o1d) = 0.
            else
                sj(ilev,iw,j_o2_o)   = 0.
                sj(ilev,iw,j_o2_o1d) = xso2*yieldo2(iw)
            end if

            end do

        end do

        !     defining the reaction types and molecules involved

        rtype(j_o2_o) = 1
                
        ns(j_o2_o) = 1
        sID(1,j_o2_o) = 7
        sISO(1,j_o2_o) = 0
        sf(1,j_o2_o) = 1.0

        npr(j_o2_o) = 1
        pID(1,j_o2_o) = 45
        pISO(1,j_o2_o) = 0
        pf(1,j_o2_o) = 2.0


        rtype(j_o2_o1d) = 1
                
        ns(j_o2_o1d) = 1
        sID(1,j_o2_o1d) = 7
        sISO(1,j_o2_o1d) = 0
        sf(1,j_o2_o1d) = 1.0

        npr(j_o2_o1d) = 2
        pID(1,j_o2_o1d) = 45
        pISO(1,j_o2_o1d) = 0
        pf(1,j_o2_o1d) = 1.0 
        pID(2,j_o2_o1d) = 133
        pISO(2,j_o2_o1d) = 0
        pf(2,j_o2_o1d) = 1.0


    end subroutine seto2

!==============================================================================================================================================

    subroutine set18o2(n_phot, nlayer, nw, wc, mopt, tlay, xs18o2_150, &
                      xs18o2_200, xs18o2_250, xs18o2_300, yieldo2, j_18o2_o, &
                      j_18o2_o1d_1, j_18o2_o1d_2, colinc, rm, dt, sj, &
                      rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the (18O)(16O) temperature-dependent cross-sections and tau       =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: mopt                                          
        integer :: j_18o2_o, j_18o2_o1d_1, j_18o2_o1d_2                             
        real, dimension(nw)     :: wc                            
        real, dimension(nw)     :: xs18o2_150,xs18o2_200,xs18o2_250,xs18o2_300
        real, dimension(nw)     :: yieldo2                       
        real, dimension(nlayer) :: tlay             !temperature             
        double precision, dimension(nlayer) :: rm               !mixing ratio of O2      
        double precision, dimension(nlayer) :: colinc           !column density (cm-2)              

        !     output:

        real, dimension(nlayer,nw)    :: dt         !optical depth             
        real, dimension(nlayer,nw,n_phot) :: sj     !cross sections * yield                

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product


        !     local:

        integer :: ilev, iw
        real    :: temp
        real    :: xs18o2, factor

        !     correction by factor if low-resolution in schumann-runge bands

        if (mopt == 1) then
            factor = 1.
        else if (mopt == 2) then
            factor = 0.8
        end if

        !     calculate temperature dependance

        do ilev = 1,nlayer
            temp = max(tlay(ilev),150.)
            temp = min(temp, 300.)

            do iw = 1, nw-1
            if (tlay(ilev) > 250.) then
                xs18o2 = xs18o2_250(iw) + (xs18o2_300(iw) - xs18o2_250(iw)) &
                                    /(300. - 250.)*(temp - 250.) 
            else if (tlay(ilev) > 200.) then
                xs18o2 = xs18o2_200(iw) + (xs18o2_250(iw) - xs18o2_200(iw)) & 
                                    /(250. - 200.)*(temp - 200.) 
            else
                xs18o2 = xs18o2_150(iw) + (xs18o2_200(iw) - xs18o2_150(iw)) &
                                    /(200. - 150.)*(temp - 150.) 
            end if

            if (wc(iw) > 180. .and. wc(iw) < 200.) then
                xs18o2 = xs18o2*factor
            end if

        !     optical depth

            dt(ilev,iw) = colinc(ilev)*rm(ilev)*xs18o2

        !     production of o(1d) for wavelengths shorter than 175 nm

            if (wc(iw) >= 175.) then
                sj(ilev,iw,j_18o2_o)   = xs18o2*yieldo2(iw)
                sj(ilev,iw,j_18o2_o1d_1) = 0.
                sj(ilev,iw,j_18o2_o1d_2) = 0.
            else
                sj(ilev,iw,j_18o2_o)   = 0.
                sj(ilev,iw,j_18o2_o1d_1) = xs18o2*yieldo2(iw)/2.0
                sj(ilev,iw,j_18o2_o1d_2) = xs18o2*yieldo2(iw)/2.0
            end if

            end do

        end do

        !     defining the reaction types and molecules involved

        rtype(j_18o2_o) = 1
                
        ns(j_18o2_o) = 1
        sID(1,j_18o2_o) = 7
        sISO(1,j_18o2_o) = 2
        sf(1,j_18o2_o) = 1.0

        npr(j_18o2_o) = 2
        pID(1,j_18o2_o) = 45
        pISO(1,j_18o2_o) = 0
        pf(1,j_18o2_o) = 1.0
        pID(2,j_18o2_o) = 45
        pISO(2,j_18o2_o) = 2
        pf(2,j_18o2_o) = 1.0


        rtype(j_18o2_o1d_1) = 1
                
        ns(j_18o2_o1d_1) = 1
        sID(1,j_18o2_o1d_1) = 7
        sISO(1,j_18o2_o1d_1) = 2
        sf(1,j_18o2_o1d_1) = 1.0

        npr(j_18o2_o1d_1) = 2
        pID(1,j_18o2_o1d_1) = 45
        pISO(1,j_18o2_o1d_1) = 2
        pf(1,j_18o2_o1d_1) = 1.0 
        pID(2,j_18o2_o1d_1) = 133
        pISO(2,j_18o2_o1d_1) = 0
        pf(2,j_18o2_o1d_1) = 1.0


        rtype(j_18o2_o1d_2) = 1
                
        ns(j_18o2_o1d_2) = 1
        sID(1,j_18o2_o1d_2) = 7
        sISO(1,j_18o2_o1d_2) = 2
        sf(1,j_18o2_o1d_2) = 1.0

        npr(j_18o2_o1d_2) = 2
        pID(1,j_18o2_o1d_2) = 45
        pISO(1,j_18o2_o1d_2) = 0
        pf(1,j_18o2_o1d_2) = 1.0 
        pID(2,j_18o2_o1d_2) = 133
        pISO(2,j_18o2_o1d_2) = 2
        pf(2,j_18o2_o1d_2) = 1.0


    end subroutine set18o2

!==============================================================================================================================================

    subroutine setco2(n_phot, nlayer, nw, wc, tlay, xsco2_195, xsco2_295, &
                    xsco2_370, yieldco2, j_co2_o, j_co2_o1d, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the CO2 temperature-dependent cross-sections and optical depth    =*
        !-----------------------------------------------------------------------------*
   
        implicit none
   
        !    input:
   
        integer :: n_phot                                              
        integer :: nlayer                                         
        integer :: nw                                              
        integer :: j_co2_o, j_co2_o1d                              
        real, dimension(nw)     :: wc                             
        real, dimension(nw)     :: xsco2_195, xsco2_295, xsco2_370 
        real, dimension(nw)     :: yieldco2                        
        real, dimension(nlayer) :: tlay                            
        double precision, dimension(nlayer) :: rm                              
        double precision, dimension(nlayer) :: colinc                           
   
        !    output:
   
        real, dimension(nlayer,nw)    :: dt                        
        real, dimension(nlayer,nw,n_phot) :: sj                        

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !   local:

        integer :: extrapol
        integer :: i, l
        real    :: temp, sco2
   
        !     extrapol  = 0         no extrapolation below 195 k
        !     extrapol  = 1         extrapolation below 195 k
        extrapol  = 0

        do i = 1, nlayer
        if (extrapol == 1) then
            temp = tlay(i)             
        else
            temp = max(tlay(i), 195.)
        end if
        temp = min(temp, 370.)
        do l = 1, nw-1
            if (temp <= 295.) then
                if (xsco2_195(l) /= 0. .and. xsco2_295(l) /= 0.) then
                    sco2 =  alog(xsco2_195(l)) &
                     + (alog(xsco2_295(l)) - alog(xsco2_195(l))) &
                       /(295. - 195.)*(temp - 195.)
                    sco2 = exp(sco2)
                else
                    sco2 = 0.
                end if
            else
                if (xsco2_295(l) /= 0. .and. xsco2_370(l) /= 0.) then
                    sco2 =  alog(xsco2_295(l)) &
                     + (alog(xsco2_370(l)) - alog(xsco2_295(l))) &
                       /(370. - 295.)*(temp - 295.)
                    sco2 = exp(sco2)
                else
                    sco2 = 0.
                end if
            end if

            !           optical depth

            dt(i,l) = colinc(i)*rm(i)*sco2
    
            !           production of o(1d) for wavelengths shorter than 167 nm

            if (wc(l) >= 167.) then
                sj(i,l,j_co2_o) = sco2*yieldco2(l)
                sj(i,l,j_co2_o1d) = 0.
            else
                sj(i,l,j_co2_o) = 0.
                sj(i,l,j_co2_o1d) = sco2*yieldco2(l)
            end if
        end do
        end do

        !defining the reaction types and molecules involved

        rtype(j_co2_o) = 1
            
        ns(j_co2_o) = 1
        sID(1,j_co2_o) = 2
        sISO(1,j_co2_o) = 0
        sf(1,j_co2_o) = 1.0

        npr(j_co2_o) = 2
        pID(1,j_co2_o) = 5
        pISO(1,j_co2_o) = 0
        pf(1,j_co2_o) = 1.0
        pID(2,j_co2_o) = 45
        pISO(2,j_co2_o) = 0
        pf(2,j_co2_o) = 1.0


        rtype(j_co2_o1d) = 1
                
        ns(j_co2_o1d) = 1
        sID(1,j_co2_o1d) = 2
        sISO(1,j_co2_o1d) = 0
        sf(1,j_co2_o1d) = 1.0

        npr(j_co2_o1d) = 2
        pID(1,j_co2_o1d) = 5
        pISO(1,j_co2_o1d) = 0
        pf(1,j_co2_o1d) = 1.0 
        pID(2,j_co2_o1d) = 133
        pISO(2,j_co2_o1d) = 0
        pf(2,j_co2_o1d) = 1.0

    end subroutine setco2

!==============================================================================================================================================

    subroutine seto3(n_phot, nlayer, nw, wc, tlay, xso3_218, xso3_298, &
                    j_o3_o, j_o3_o1d, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the O3 temperature dependent cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_o3_o, j_o3_o1d                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xso3_218, xso3_298            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        !
        integer :: ilev, iw
        real :: temp
        real, dimension(nw) :: xso3(nw)
        real, dimension(nw) :: qy1d              
        real :: q1, q2, a1, a2, a3

        do ilev = 1, nlayer
        temp = max(tlay(ilev), 218.)
        temp = min(temp,298.)
        do iw = 1, nw-1
            xso3(iw) = xso3_218(iw) + (xso3_298(iw) - xso3_218(iw)) &
                                    /(298. - 218.) *(temp - 218.) 

        !     optical depth

            dt(ilev,iw) = colinc(ilev)*rm(ilev)*xso3(iw)

        end do

        !     calculate quantum yield for o(1d) production (jpl 2006)

        temp = max(tlay(ilev),200.)
        temp = min(temp,320.)
        do iw = 1, nw-1
            if (wc(iw) <= 306.) then
                qy1d(iw) = 0.90
            else if (wc(iw) > 306. .and. wc(iw) < 328.) then
                q1 = 1.
                q2 = exp(-825.518/(0.695*temp))
                a1 = (304.225 - wc(iw))/5.576
                a2 = (314.957 - wc(iw))/6.601
                a3 = (310.737 - wc(iw))/2.187
                qy1d(iw) = (q1/(q1 + q2))*0.8036*exp(-(a1*a1*a1*a1)) &
                        + (q2/(q1 + q2))*8.9061*(temp/300.)**2. &
                        *exp(-(a2*a2)) &
                        + 0.1192*(temp/300.)**1.5*exp(-(a3*a3)) &
                        + 0.0765
            else if (wc(iw) >= 328. .and. wc(iw) <= 340.) then
                qy1d(iw) = 0.08
            else
                qy1d(iw) = 0.
            endif
        end do
        do iw = 1, nw-1
            sj(ilev,iw,j_o3_o)   = xso3(iw)*(1. - qy1d(iw)) 
            sj(ilev,iw,j_o3_o1d) = xso3(iw)*qy1d(iw)
        end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_o3_o) = 1
            
        ns(j_o3_o) = 1
        sID(1,j_o3_o) = 3
        sISO(1,j_o3_o) = 0
        sf(1,j_o3_o) = 1.0

        npr(j_o3_o) = 2
        pID(1,j_o3_o) = 7
        pISO(1,j_o3_o) = 0
        pf(1,j_o3_o) = 1.0
        pID(2,j_o3_o) = 45
        pISO(2,j_o3_o) = 0
        pf(2,j_o3_o) = 1.0


        rtype(j_o3_o1d) = 1
                
        ns(j_o3_o1d) = 1
        sID(1,j_o3_o1d) = 3
        sISO(1,j_o3_o1d) = 0
        sf(1,j_o3_o1d) = 1.0

        npr(j_o3_o1d) = 2
        pID(1,j_o3_o1d) = 7
        pISO(1,j_o3_o1d) = 0
        pf(1,j_o3_o1d) = 1.0 
        pID(2,j_o3_o1d) = 133
        pISO(2,j_o3_o1d) = 0
        pf(2,j_o3_o1d) = 1.0

    end subroutine seto3

!==============================================================================================================================================

    subroutine set18o3(n_phot, nlayer, nw, wc, tlay, xs18o3_218, xs18o3_298, &
                    j_18o3_o_1, j_18o3_o_2, j_18o3_o1d_1, j_18o3_o1d_2, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the O3 temperature dependent cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_18o3_o_1, j_18o3_o_2, j_18o3_o1d_1, j_18o3_o1d_2                             
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xs18o3_218, xs18o3_298            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        !
        integer :: ilev, iw
        real :: temp
        real, dimension(nw) :: xs18o3(nw)
        real, dimension(nw) :: qy1d              
        real :: q1, q2, a1, a2, a3

        do ilev = 1, nlayer
        temp = max(tlay(ilev), 218.)
        temp = min(temp,298.)
        do iw = 1, nw-1
            xs18o3(iw) = xs18o3_218(iw) + (xs18o3_298(iw) - xs18o3_218(iw)) &
                                    /(298. - 218.) *(temp - 218.) 

        !     optical depth

            dt(ilev,iw) = colinc(ilev)*rm(ilev)*xs18o3(iw)

        end do

        !     calculate quantum yield for o(1d) production (jpl 2006)

        temp = max(tlay(ilev),200.)
        temp = min(temp,320.)
        do iw = 1, nw-1
            if (wc(iw) <= 306.) then
                qy1d(iw) = 0.90
            else if (wc(iw) > 306. .and. wc(iw) < 328.) then
                q1 = 1.
                q2 = exp(-825.518/(0.695*temp))
                a1 = (304.225 - wc(iw))/5.576
                a2 = (314.957 - wc(iw))/6.601
                a3 = (310.737 - wc(iw))/2.187
                qy1d(iw) = (q1/(q1 + q2))*0.8036*exp(-(a1*a1*a1*a1)) &
                        + (q2/(q1 + q2))*8.9061*(temp/300.)**2. &
                        *exp(-(a2*a2)) &
                        + 0.1192*(temp/300.)**1.5*exp(-(a3*a3)) &
                        + 0.0765
            else if (wc(iw) >= 328. .and. wc(iw) <= 340.) then
                qy1d(iw) = 0.08
            else
                qy1d(iw) = 0.
            endif
        end do
        do iw = 1, nw-1
            sj(ilev,iw,j_18o3_o_1)   = xs18o3(iw)*(1. - qy1d(iw))/2.0 
            sj(ilev,iw,j_18o3_o_2)   = xs18o3(iw)*(1. - qy1d(iw))/2.0 
            sj(ilev,iw,j_18o3_o1d_1) = xs18o3(iw)*qy1d(iw)/2.0
            sj(ilev,iw,j_18o3_o1d_2) = xs18o3(iw)*qy1d(iw)/2.0
        end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_18o3_o_1) = 1
            
        ns(j_18o3_o_1) = 1
        sID(1,j_18o3_o_1) = 3
        sISO(1,j_18o3_o_1) = 2
        sf(1,j_18o3_o_1) = 1.0

        npr(j_18o3_o_1) = 2
        pID(1,j_18o3_o_1) = 7
        pISO(1,j_18o3_o_1) = 2
        pf(1,j_18o3_o_1) = 1.0
        pID(2,j_18o3_o_1) = 45
        pISO(2,j_18o3_o_1) = 0
        pf(2,j_18o3_o_1) = 1.0


        rtype(j_18o3_o_2) = 1
            
        ns(j_18o3_o_2) = 1
        sID(1,j_18o3_o_2) = 3
        sISO(1,j_18o3_o_2) = 2
        sf(1,j_18o3_o_2) = 1.0

        npr(j_18o3_o_2) = 2
        pID(1,j_18o3_o_2) = 7
        pISO(1,j_18o3_o_2) = 0
        pf(1,j_18o3_o_2) = 1.0
        pID(2,j_18o3_o_2) = 45
        pISO(2,j_18o3_o_2) = 2
        pf(2,j_18o3_o_2) = 1.0


        rtype(j_18o3_o1d_1) = 1
                
        ns(j_18o3_o1d_1) = 1
        sID(1,j_18o3_o1d_1) = 3
        sISO(1,j_18o3_o1d_1) = 2
        sf(1,j_18o3_o1d_1) = 1.0

        npr(j_18o3_o1d_1) = 2
        pID(1,j_18o3_o1d_1) = 7
        pISO(1,j_18o3_o1d_1) = 2
        pf(1,j_18o3_o1d_1) = 1.0 
        pID(2,j_18o3_o1d_1) = 133
        pISO(2,j_18o3_o1d_1) = 0
        pf(2,j_18o3_o1d_1) = 1.0


        rtype(j_18o3_o1d_2) = 1
                
        ns(j_18o3_o1d_2) = 1
        sID(1,j_18o3_o1d_2) = 3
        sISO(1,j_18o3_o1d_2) = 2
        sf(1,j_18o3_o1d_2) = 1.0

        npr(j_18o3_o1d_2) = 2
        pID(1,j_18o3_o1d_2) = 7
        pISO(1,j_18o3_o1d_2) = 0
        pf(1,j_18o3_o1d_2) = 1.0 
        pID(2,j_18o3_o1d_2) = 133
        pISO(2,j_18o3_o1d_2) = 2
        pf(2,j_18o3_o1d_2) = 1.0


    end subroutine set18o3

!==============================================================================================================================================

    subroutine seth2o(n_phot, nlayer, nw, wc, tlay, xsh2o, j_h2o, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the H2O cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_h2o                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xsh2o            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_h2o)   = xsh2o(iw) 
                sj(ilev,iw,j_h2o) = xsh2o(iw)

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xsh2o(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_h2o) = 1
            
        ns(j_h2o) = 1
        sID(1,j_h2o) = 1   !H2O
        sISO(1,j_h2o) = 0
        sf(1,j_h2o) = 1.0

        npr(j_h2o) = 2
        pID(1,j_h2o) = 13  !OH
        pISO(1,j_h2o) = 0
        pf(1,j_h2o) = 1.0
        pID(2,j_h2o) = 48  !H
        pISO(2,j_h2o) = 0
        pf(2,j_h2o) = 1.0

    end subroutine seth2o

!==============================================================================================================================================

    subroutine set18h2o(n_phot, nlayer, nw, wc, tlay, xs18h2o, j_18h2o, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the H2(18O) cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_18h2o                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xs18h2o            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_18h2o)   = xs18h2o(iw) 

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xs18h2o(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_18h2o) = 1
            
        ns(j_18h2o) = 1
        sID(1,j_18h2o) = 1   !H2O
        sISO(1,j_18h2o) = 2
        sf(1,j_18h2o) = 1.0

        npr(j_18h2o) = 2
        pID(1,j_18h2o) = 13  !OH
        pISO(1,j_18h2o) = 2
        pf(1,j_18h2o) = 1.0
        pID(2,j_18h2o) = 48  !H
        pISO(2,j_18h2o) = 0
        pf(2,j_18h2o) = 1.0

    end subroutine set18h2o

!==============================================================================================================================================

    subroutine seth2o2(n_phot, nlayer, nw, wc, tlay, xsh2o2, j_h2o2,&
                        colinc, rm, dt, sj, &
                        rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the h2o2 temperature dependent cross-sections and optical depth   =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                           
        integer :: j_h2o2                                        
        real, dimension(nw)     :: wc                            
        real, dimension(nw)     :: xsh2o2                        
        real, dimension(nlayer) :: tlay                          
        double precision, dimension(nlayer) :: rm                           
        double precision, dimension(nlayer) :: colinc                     

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj       
        
        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:

        integer :: ilev, iw
        real    :: a0, a1, a2, a3, a4, a5, a6, a7
        real    :: b0, b1, b2, b3, b4
        real    :: lambda, suma, sumb, chi, temp, xs

        A0 = 6.4761E+04            
        A1 = -9.2170972E+02        
        A2 = 4.535649              
        A3 = -4.4589016E-03        
        A4 = -4.035101E-05         
        A5 = 1.6878206E-07
        A6 = -2.652014E-10
        A7 = 1.5534675E-13

        B0 = 6.8123E+03
        B1 = -5.1351E+01
        B2 = 1.1522E-01
        B3 = -3.0493E-05
        B4 = -1.0924E-07

        !     temperature dependance: jpl 2006

        do ilev = 1,nlayer
            temp = min(max(tlay(ilev),200.),400.)            
            chi = 1./(1. + exp(-1265./temp))
            do iw = 1, nw-1
                if ((wc(iw) >= 260.) .and. (wc(iw) < 350.)) then
                lambda = wc(iw)
                sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda + &
                              A4)*lambda +A3)*lambda + A2)*lambda + &
                              A1)*lambda + A0
                sumB = (((B4*lambda + B3)*lambda + B2)*lambda + &
                           B1)*lambda + B0
                xs = (chi*sumA + (1. - chi)*sumB)*1.e-21
                sj(ilev,iw,j_h2o2) = xs
                else
                sj(ilev,iw,j_h2o2) = xsh2o2(iw)
                end if

                !     optical depth

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*sj(ilev,iw,j_h2o2)
            end do
        end do

        !defining the reaction types and molecules involved

        rtype(j_h2o2) = 1
            
        ns(j_h2o2) = 1
        sID(1,j_h2o2) = 25   !H2O2
        sISO(1,j_h2o2) = 0
        sf(1,j_h2o2) = 1.0

        npr(j_h2o2) = 1
        pID(1,j_h2o2) = 13  !OH
        pISO(1,j_h2o2) = 0
        pf(1,j_h2o2) = 2.0

    end subroutine seth2o2

!==============================================================================================================================================

    subroutine set18h2o2(n_phot, nlayer, nw, wc, tlay, xs18h2o2, j_18h2o2,&
                        colinc, rm, dt, sj, &
                        rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the h2o2 temperature dependent cross-sections and optical depth   =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                           
        integer :: j_18h2o2                                        
        real, dimension(nw)     :: wc                            
        real, dimension(nw)     :: xs18h2o2                        
        real, dimension(nlayer) :: tlay                          
        double precision, dimension(nlayer) :: rm                           
        double precision, dimension(nlayer) :: colinc                     

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj       
        
        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:

        integer :: ilev, iw
        real    :: a0, a1, a2, a3, a4, a5, a6, a7
        real    :: b0, b1, b2, b3, b4
        real    :: lambda, suma, sumb, chi, temp, xs

        A0 = 6.4761E+04            
        A1 = -9.2170972E+02        
        A2 = 4.535649              
        A3 = -4.4589016E-03        
        A4 = -4.035101E-05         
        A5 = 1.6878206E-07
        A6 = -2.652014E-10
        A7 = 1.5534675E-13

        B0 = 6.8123E+03
        B1 = -5.1351E+01
        B2 = 1.1522E-01
        B3 = -3.0493E-05
        B4 = -1.0924E-07

        !     temperature dependance: jpl 2006

        do ilev = 1,nlayer
            temp = min(max(tlay(ilev),200.),400.)            
            chi = 1./(1. + exp(-1265./temp))
            do iw = 1, nw-1
                if ((wc(iw) >= 260.) .and. (wc(iw) < 350.)) then
                lambda = wc(iw)
                sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda + &
                              A4)*lambda +A3)*lambda + A2)*lambda + &
                              A1)*lambda + A0
                sumB = (((B4*lambda + B3)*lambda + B2)*lambda + &
                           B1)*lambda + B0
                xs = (chi*sumA + (1. - chi)*sumB)*1.e-21
                sj(ilev,iw,j_18h2o2) = xs
                else
                sj(ilev,iw,j_18h2o2) = xs18h2o2(iw)
                end if

                !     optical depth

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*sj(ilev,iw,j_18h2o2)
            end do
        end do

        !defining the reaction types and molecules involved

        rtype(j_18h2o2) = 1
            
        ns(j_18h2o2) = 1
        sID(1,j_18h2o2) = 25   !H2(18O)O
        sISO(1,j_18h2o2) = 2
        sf(1,j_18h2o2) = 1.0

        npr(j_18h2o2) = 2
        pID(1,j_18h2o2) = 13  !(18O)H
        pISO(1,j_18h2o2) = 2
        pf(1,j_18h2o2) = 1.0
        pID(2,j_18h2o2) = 13  !OH
        pISO(2,j_18h2o2) = 0
        pf(2,j_18h2o2) = 1.0

    end subroutine set18h2o2

!==============================================================================================================================================

    subroutine setho2(n_phot, nlayer, nw, wc, tlay, xsho2,j_ho2, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the HO2 cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_ho2                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xsho2            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_ho2)   = xsho2(iw) 
                sj(ilev,iw,j_ho2) = xsho2(iw)

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xsho2(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_ho2) = 1
            
        ns(j_ho2) = 1
        sID(1,j_ho2) = 44   !HO2
        sISO(1,j_ho2) = 0
        sf(1,j_ho2) = 1.0

        npr(j_ho2) = 2
        pID(1,j_ho2) = 13  !OH
        pISO(1,j_ho2) = 0
        pf(1,j_ho2) = 1.0
        pID(2,j_ho2) = 45  !O
        pISO(2,j_ho2) = 0
        pf(2,j_ho2) = 1.0

    end subroutine setho2

!==============================================================================================================================================

    subroutine set18ho2(n_phot, nlayer, nw, wc, tlay, xs18ho2,j_18ho2_1, j_18ho2_2, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the H(18O)O cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_18ho2_1, j_18ho2_2                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xs18ho2            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_18ho2_1) = xs18ho2(iw)/2.0 
                sj(ilev,iw,j_18ho2_2) = xs18ho2(iw)/2.0

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xs18ho2(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_18ho2_1) = 1
            
        ns(j_18ho2_1) = 1
        sID(1,j_18ho2_1) = 44   !H(18O)O
        sISO(1,j_18ho2_1) = 2
        sf(1,j_18ho2_1) = 1.0

        npr(j_18ho2_1) = 2
        pID(1,j_18ho2_1) = 13  !(18O)H
        pISO(1,j_18ho2_1) = 2
        pf(1,j_18ho2_1) = 1.0
        pID(2,j_18ho2_1) = 45  !O
        pISO(2,j_18ho2_1) = 0
        pf(2,j_18ho2_1) = 1.0


        rtype(j_18ho2_2) = 1
            
        ns(j_18ho2_2) = 1
        sID(1,j_18ho2_2) = 44   !H(18O)O
        sISO(1,j_18ho2_2) = 2
        sf(1,j_18ho2_2) = 1.0

        npr(j_18ho2_2) = 2
        pID(1,j_18ho2_2) = 13  !OH
        pISO(1,j_18ho2_2) = 0
        pf(1,j_18ho2_2) = 1.0
        pID(2,j_18ho2_2) = 45  !18O
        pISO(2,j_18ho2_2) = 2
        pf(2,j_18ho2_2) = 1.0

    end subroutine set18ho2

!==============================================================================================================================================

    subroutine seth2(n_phot, nlayer, nw, wc, tlay, xsh2, yieldh2, j_h2, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the HO2 cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_h2                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xsh2,yieldh2            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_h2)   = xsh2(iw) * yieldh2(iw) 

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xsh2(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_h2) = 1
            
        ns(j_h2) = 1
        sID(1,j_h2) = 39   !H2
        sISO(1,j_h2) = 0
        sf(1,j_h2) = 1.0

        npr(j_h2) = 1
        pID(1,j_h2) = 48  !H
        pISO(1,j_h2) = 0
        pf(1,j_h2) = 2.0

    end subroutine seth2

!==============================================================================================================================================

    subroutine setno2(n_phot, nlayer, nw, wc, tlay, xsno2, xsno2_220,&
                       xsno2_294, yldno2_248, yldno2_298, j_no2,&
                       colinc, rm, dt, sj,&
                       rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the no2 temperature-dependent cross-sections and optical depth    =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_no2                                         
        real, dimension(nw)     :: wc                            
        real, dimension(nw) :: xsno2, xsno2_220, xsno2_294       
        real, dimension(nw) :: yldno2_248, yldno2_298            
        real, dimension(nlayer) :: tlay                          
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                    

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj       
        
        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:

        integer :: ilev, iw
        real    :: temp, qy

        !     temperature dependance: jpl 2006

        do ilev = 1,nlayer
            temp = max(220.,min(tlay(ilev),294.))
            do iw = 1, nw - 1
                if (wc(iw) < 238.) then
                sj(ilev,iw,j_no2) = xsno2(iw)
                else
                sj(ilev,iw,j_no2) = xsno2_220(iw) &
                                    + (xsno2_294(iw) - xsno2_220(iw)) &
                                    /(294. - 220.)*(temp - 220.)
                end if

        !     optical depth

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*sj(ilev,iw,j_no2)
            end do
        end do

        !     quantum yield: jpl 2006

        do ilev = 1,nlayer
            temp = max(248.,min(tlay(ilev),298.))
            do iw = 1, nw - 1
                qy = yldno2_248(iw) + (yldno2_298(iw) - yldno2_248(iw)) &
                                    /(298. - 248.)*(temp - 248.)
                sj(ilev,iw,j_no2) = sj(ilev,iw,j_no2)*qy 
            end do
        end do

        !defining the reaction types and molecules involved

        rtype(j_no2) = 1
            
        ns(j_no2) = 1
        sID(1,j_no2) = 10   !NO2
        sISO(1,j_no2) = 0
        sf(1,j_no2) = 1.0

        npr(j_no2) = 2
        pID(1,j_no2) = 8    !NO
        pISO(1,j_no2) = 0
        pf(1,j_no2) = 1.0
        pID(2,j_no2) = 45   !O
        pISO(2,j_no2) = 0
        pf(2,j_no2) = 1.0


    end subroutine setno2

!==============================================================================================================================================

    subroutine set18no2(n_phot, nlayer, nw, wc, tlay, xs18no2, xs18no2_220,&
                       xs18no2_294, yldno2_248, yldno2_298, j_18no2_1, j_18no2_2,&
                       colinc, rm, dt, sj,&
                       rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the NO(18O) temp-dependent cross-sections and optical depth       =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_18no2_1, j_18no2_2                                 
        real, dimension(nw)     :: wc                            
        real, dimension(nw) :: xs18no2, xs18no2_220, xs18no2_294       
        real, dimension(nw) :: yldno2_248, yldno2_298            
        real, dimension(nlayer) :: tlay                          
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                    

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj       
        
        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:

        integer :: ilev, iw
        real    :: temp, qy

        !     temperature dependance: jpl 2006

        do ilev = 1,nlayer
            temp = max(220.,min(tlay(ilev),294.))
            do iw = 1, nw - 1
                if (wc(iw) < 238.) then
                 sj(ilev,iw,j_18no2_1) = xs18no2(iw)/2.0
                 sj(ilev,iw,j_18no2_2) = xs18no2(iw)/2.0
                else
                 sj(ilev,iw,j_18no2_1) = (xs18no2_220(iw) &
                                    + (xs18no2_294(iw) - xs18no2_220(iw)) &
                                    /(294. - 220.)*(temp - 220.))/2.0
                 sj(ilev,iw,j_18no2_2) = (xs18no2_220(iw) &
                                    + (xs18no2_294(iw) - xs18no2_220(iw)) &
                                    /(294. - 220.)*(temp - 220.))/2.0
                end if

        !     optical depth

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*(sj(ilev,iw,j_18no2_1)+sj(ilev,iw,j_18no2_2))
            end do
        end do

        !     quantum yield: jpl 2006

        do ilev = 1,nlayer
            temp = max(248.,min(tlay(ilev),298.))
            do iw = 1, nw - 1
                qy = yldno2_248(iw) + (yldno2_298(iw) - yldno2_248(iw)) &
                                    /(298. - 248.)*(temp - 248.)
                sj(ilev,iw,j_18no2_1) = sj(ilev,iw,j_18no2_1)*qy 
                sj(ilev,iw,j_18no2_2) = sj(ilev,iw,j_18no2_2)*qy
            end do
        end do

        !defining the reaction types and molecules involved

        rtype(j_18no2_1) = 1
            
        ns(j_18no2_1) = 1
        sID(1,j_18no2_1) = 10   !NO(18O)
        sISO(1,j_18no2_1) = 3
        sf(1,j_18no2_1) = 1.0

        npr(j_18no2_1) = 2
        pID(1,j_18no2_1) = 8    !N(18O)
        pISO(1,j_18no2_1) = 3
        pf(1,j_18no2_1) = 1.0
        pID(2,j_18no2_1) = 45   !O
        pISO(2,j_18no2_1) = 0
        pf(2,j_18no2_1) = 1.0


        rtype(j_18no2_2) = 1
            
        ns(j_18no2_2) = 1
        sID(1,j_18no2_2) = 10   !NO(18O)
        sISO(1,j_18no2_2) = 3
        sf(1,j_18no2_2) = 1.0

        npr(j_18no2_2) = 2
        pID(1,j_18no2_2) = 8    !NO
        pISO(1,j_18no2_2) = 0
        pf(1,j_18no2_2) = 1.0
        pID(2,j_18no2_2) = 45   !(18O)
        pISO(2,j_18no2_2) = 2
        pf(2,j_18no2_2) = 1.0


    end subroutine set18no2

!==============================================================================================================================================

    subroutine setno(n_phot, nlayer, nw, wc, tlay, xsno, yieldno, j_no, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the NO cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_no                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xsno,yieldno            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_no)   = xsno(iw) * yieldno(iw) 

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xsno(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_no) = 1
            
        ns(j_no) = 1
        sID(1,j_no) = 8   !NO
        sISO(1,j_no) = 0
        sf(1,j_no) = 1.0

        npr(j_no) = 2
        pID(1,j_no) = 47   !N
        pISO(1,j_no) = 0
        pf(1,j_no) = 1.0
        pID(2,j_no) = 45   !O
        pISO(2,j_no) = 0
        pf(2,j_no) = 1.0

    end subroutine setno

!==============================================================================================================================================

    subroutine set18no(n_phot, nlayer, nw, wc, tlay, xs18no, yieldno, j_18no, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the N(18O) cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_18no                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xs18no,yieldno            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_18no)   = xs18no(iw) * yieldno(iw) 

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xs18no(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_18no) = 1
            
        ns(j_18no) = 1
        sID(1,j_18no) = 8   !NO
        sISO(1,j_18no) = 3
        sf(1,j_18no) = 1.0

        npr(j_18no) = 2
        pID(1,j_18no) = 47   !N
        pISO(1,j_18no) = 0
        pf(1,j_18no) = 1.0
        pID(2,j_18no) = 45   !O
        pISO(2,j_18no) = 2
        pf(2,j_18no) = 1.0

    end subroutine set18no

!==============================================================================================================================================

    subroutine setn2(n_phot, nlayer, nw, wc, tlay, xsn2, yieldn2, j_n2, &
                    colinc, rm, dt, sj,&
                    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the N2 cross-sections and optical depth     =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                            
        integer :: nlayer                                        
        integer :: nw                                            
        integer :: j_n2                              
        real, dimension(nw)     :: wc                          
        real, dimension(nw)     :: xsn2,yieldn2            
        real, dimension(nlayer) :: tlay                         
        double precision, dimension(nlayer) :: rm                            
        double precision, dimension(nlayer) :: colinc                       

        !     output:

        real, dimension(nlayer,nw)    :: dt                      
        real, dimension(nlayer,nw,n_phot) :: sj                     

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:
        
        integer :: ilev, iw            

        do ilev = 1, nlayer

            !calculate the cross sections and optical depth

            do iw = 1, nw-1
                sj(ilev,iw,j_n2)   = xsn2(iw) * yieldn2(iw) 

                dt(ilev,iw) = colinc(ilev)*rm(ilev)*xsn2(iw)
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_n2) = 1
            
        ns(j_n2) = 1
        sID(1,j_n2) = 22   !N2
        sISO(1,j_n2) = 0
        sf(1,j_n2) = 1.0

        npr(j_n2) = 1
        pID(1,j_n2) = 47  !N
        pISO(1,j_n2) = 0
        pf(1,j_n2) = 2.0

    end subroutine setn2

!==============================================================================================================================================

    subroutine set13co2(n_phot,nlayer,nw,wc,tlay,xs13co2_195,xs13co2_295,& 
                          xs13co2_370, yieldco2, j_13co2_o, j_13co2_o1d, &
                          colinc, rm, dt, sj, &
                          rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the (13C)O2 temperature-dependent cross-sections and optical depth    =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                              
        integer :: nlayer                                         
        integer :: nw                                              
        integer :: j_13co2_o, j_13co2_o1d                              
        real, dimension(nw)     :: wc                             
        real, dimension(nw)     :: xs13co2_195, xs13co2_295, xs13co2_370 
        real, dimension(nw)     :: yieldco2                        
        real, dimension(nlayer) :: tlay                            
        double precision, dimension(nlayer) :: rm                              
        double precision, dimension(nlayer) :: colinc                           

        !     output:

        real, dimension(nlayer,nw)    :: dt                        
        real, dimension(nlayer,nw,n_phot) :: sj                        

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:

        integer :: extrapol
        integer :: i, l
        real    :: temp, sco2

        !     extrapol  = 0         no extrapolation below 195 k
        !     extrapol  = 1         extrapolation below 195 k

        extrapol  = 0

        do i = 1, nlayer
            if (extrapol == 1) then
                temp = tlay(i)             
            else
                temp = max(tlay(i), 195.)
            end if
            temp = min(temp, 370.)
            do l = 1, nw-1
                if (temp <= 295.) then
                    if (xs13co2_195(l) /= 0. .and. xs13co2_295(l) /= 0.) then
                        sco2 =  alog(xs13co2_195(l)) &
                        + (alog(xs13co2_295(l)) - alog(xs13co2_195(l))) &
                        /(295. - 195.)*(temp - 195.)
                        sco2 = exp(sco2)
                    else
                        sco2 = 0.
                    end if
                else
                    if (xs13co2_295(l) /= 0. .and. xs13co2_370(l) /= 0.) then
                        sco2 =  alog(xs13co2_295(l)) &
                        + (alog(xs13co2_370(l)) - alog(xs13co2_295(l))) &
                        /(370. - 295.)*(temp - 295.)
                        sco2 = exp(sco2)
                    else
                        sco2 = 0.
                    end if
                end if

                !           optical depth

                dt(i,l) = colinc(i)*rm(i)*sco2
        
                !           production of o(1d) for wavelengths shorter than 167 nm

                if (wc(l) >= 167.) then
                    sj(i,l,j_13co2_o) = sco2*yieldco2(l)
                    sj(i,l,j_13co2_o1d) = 0.
                else
                    sj(i,l,j_13co2_o) = 0.
                    sj(i,l,j_13co2_o1d) = sco2*yieldco2(l)
                end if
            end do
        end do
   
        !defining the reaction types and molecules involved

        rtype(j_13co2_o) = 1
            
        ns(j_13co2_o) = 1
        sID(1,j_13co2_o) = 2
        sISO(1,j_13co2_o) = 2
        sf(1,j_13co2_o) = 1.0

        npr(j_13co2_o) = 2
        pID(1,j_13co2_o) = 5
        pISO(1,j_13co2_o) = 2
        pf(1,j_13co2_o) = 1.0
        pID(2,j_13co2_o) = 45
        pISO(2,j_13co2_o) = 0
        pf(2,j_13co2_o) = 1.0


        rtype(j_13co2_o1d) = 1
                
        ns(j_13co2_o1d) = 1
        sID(1,j_13co2_o1d) = 2
        sISO(1,j_13co2_o1d) = 2
        sf(1,j_13co2_o1d) = 1.0

        npr(j_13co2_o1d) = 2
        pID(1,j_13co2_o1d) = 5
        pISO(1,j_13co2_o1d) = 2
        pf(1,j_13co2_o1d) = 1.0 
        pID(2,j_13co2_o1d) = 133
        pISO(2,j_13co2_o1d) = 0
        pf(2,j_13co2_o1d) = 1.0


    end subroutine set13co2

!==============================================================================================================================================

    subroutine set18co2(n_phot,nlayer,nw,wc,tlay,xs18co2_195,xs18co2_295,& 
                          xs18co2_370, yieldco2, j_18co2_o_1, j_18co2_o_2, j_18co2_o1d_1, j_18co2_o1d_2, &
                          colinc, rm, dt, sj, &
                          rtype, ns, sID, sISO, sf, npr, pID, pISO, pf)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set up the (18O)(12C)(16O) temperature-dependent cross-sections and optical depth    =*
        !-----------------------------------------------------------------------------*

        implicit none

        !     input:

        integer :: n_phot                                              
        integer :: nlayer                                         
        integer :: nw                                              
        integer :: j_18co2_o_1, j_18co2_o_2, j_18co2_o1d_1, j_18co2_o1d_2                              
        real, dimension(nw)     :: wc                             
        real, dimension(nw)     :: xs18co2_195, xs18co2_295, xs18co2_370 
        real, dimension(nw)     :: yieldco2                        
        real, dimension(nlayer) :: tlay                            
        double precision, dimension(nlayer) :: rm                              
        double precision, dimension(nlayer) :: colinc                           

        !     output:

        real, dimension(nlayer,nw)    :: dt                        
        real, dimension(nlayer,nw,n_phot) :: sj                        

        integer :: rtype(n_phot)                    ! reaction type - always 1 for photolysis
        integer :: ns(n_phot),npr(n_phot)           ! number of source and products
        integer :: sID(2,n_phot),sISO(2,n_phot)     ! ID of the source molecules
        integer :: pID(2,n_phot),pISO(2,n_phot)     ! ID of the product molecules
        real :: sf(2,n_phot),pf(2,n_phot)           ! number of molecules per source/product

        !     local:

        integer :: extrapol
        integer :: i, l
        real    :: temp, sco2

        !     extrapol  = 0         no extrapolation below 195 k
        !     extrapol  = 1         extrapolation below 195 k

        extrapol  = 0

        do i = 1, nlayer
            if (extrapol == 1) then
                temp = tlay(i)             
            else
                temp = max(tlay(i), 195.)
            end if
            temp = min(temp, 370.)
            do l = 1, nw-1
                if (temp <= 295.) then
                    if (xs18co2_195(l) /= 0. .and. xs18co2_295(l) /= 0.) then
                        sco2 =  alog(xs18co2_195(l)) &
                        + (alog(xs18co2_295(l)) - alog(xs18co2_195(l))) &
                        /(295. - 195.)*(temp - 195.)
                        sco2 = exp(sco2)
                    else
                        sco2 = 0.
                    end if
                else
                    if (xs18co2_295(l) /= 0. .and. xs18co2_370(l) /= 0.) then
                        sco2 =  alog(xs18co2_295(l)) &
                        + (alog(xs18co2_370(l)) - alog(xs18co2_295(l))) &
                        /(370. - 295.)*(temp - 295.)
                        sco2 = exp(sco2)
                    else
                        sco2 = 0.
                    end if
                end if

                !           optical depth

                dt(i,l) = colinc(i)*rm(i)*sco2
        
                !           production of o(1d) for wavelengths shorter than 167 nm

                if (wc(l) >= 167.) then
                    sj(i,l,j_18co2_o_1) = sco2*yieldco2(l)/2.
                    sj(i,l,j_18co2_o_2) = sco2*yieldco2(l)/2.
                    sj(i,l,j_18co2_o1d_1) = 0.
                    sj(i,l,j_18co2_o1d_2) = 0.
                else
                    sj(i,l,j_18co2_o_1) = 0.
                    sj(i,l,j_18co2_o_2) = 0.
                    sj(i,l,j_18co2_o1d_1) = sco2*yieldco2(l)/2.
                    sj(i,l,j_18co2_o1d_2) = sco2*yieldco2(l)/2.
                end if
            end do
        end do
   
        !defining the reaction types and molecules involved

        !(18O)(12C)(16O) + hv  ---->  (12C)(18O) + (16O)   
        rtype(j_18co2_o_1) = 1
            
        ns(j_18co2_o_1) = 1
        sID(1,j_18co2_o_1) = 2
        sISO(1,j_18co2_o_1) = 3
        sf(1,j_18co2_o_1) = 1.0

        npr(j_18co2_o_1) = 2
        pID(1,j_18co2_o_1) = 5
        pISO(1,j_18co2_o_1) = 3
        pf(1,j_18co2_o_1) = 1.0
        pID(2,j_18co2_o_1) = 45
        pISO(2,j_18co2_o_1) = 0
        pf(2,j_18co2_o_1) = 1.0


        !(18O)(12C)(16O) + hv   ---->    (12C)(16O) + 18O  !For now we do not produce 18O
        rtype(j_18co2_o_2) = 1
            
        ns(j_18co2_o_2) = 1
        sID(1,j_18co2_o_2) = 2
        sISO(1,j_18co2_o_2) = 3
        sf(1,j_18co2_o_2) = 1.0

        npr(j_18co2_o_2) = 2
        pID(1,j_18co2_o_2) = 5
        pISO(1,j_18co2_o_2) = 0
        pf(1,j_18co2_o_2) = 1.0
        pID(2,j_18co2_o_2) = 45
        pISO(2,j_18co2_o_2) = 0
        pf(2,j_18co2_o_2) = 1.0


        !(18O)(12C)(16O) + hv   ---->    (12C)(18O) + O(1D)
        rtype(j_18co2_o1d_1) = 1
                
        ns(j_18co2_o1d_1) = 1
        sID(1,j_18co2_o1d_1) = 2
        sISO(1,j_18co2_o1d_1) = 3
        sf(1,j_18co2_o1d_1) = 1.0

        npr(j_18co2_o1d_1) = 2
        pID(1,j_18co2_o1d_1) = 5
        pISO(1,j_18co2_o1d_1) = 3
        pf(1,j_18co2_o1d_1) = 1.0 
        pID(2,j_18co2_o1d_1) = 133
        pISO(2,j_18co2_o1d_1) = 0
        pf(2,j_18co2_o1d_1) = 1.0


        !(18O)(12C)(16O) + hv   ---->    (12C)(16O) + 18O(1D)  
        rtype(j_18co2_o1d_2) = 1
                
        ns(j_18co2_o1d_2) = 1
        sID(1,j_18co2_o1d_2) = 2
        sISO(1,j_18co2_o1d_2) = 3
        sf(1,j_18co2_o1d_2) = 1.0

        npr(j_18co2_o1d_2) = 2
        pID(1,j_18co2_o1d_2) = 5
        pISO(1,j_18co2_o1d_2) = 0
        pf(1,j_18co2_o1d_2) = 1.0 
        pID(2,j_18co2_o1d_2) = 133
        pISO(2,j_18co2_o1d_2) = 0
        pf(2,j_18co2_o1d_2) = 1.0



    end subroutine set18co2

!==============================================================================================================================================

    subroutine setairlay(nlev, nw, wl, wc, colinc, dtrl)
    
        !*-----------------------------------------------------------------------------*
        !*=  PURPOSE:                                                                 =*
        !*=  computes rayleigh optical depth                                          =*
        !*-----------------------------------------------------------------------------*
        
        implicit none
    
        !     input:
    
        integer, intent(in) :: nlev, nw
    
        real, intent(in)   :: wl(nw), wc(nw)             
        double precision, intent(in) :: colinc(nlev)
    
        !     output:
            
        real, intent(out) :: dtrl(nlev,nw)  
    
        !     local:
    
        real :: nu
        real, dimension(nw) :: srayl
        integer :: ilev, iw
    
        !     compute column increments 
    
        do iw = 1, nw - 1
    
            !        co2 rayleigh cross-section
            !        ityaksov et al., chem. phys. lett., 462, 31-34, 2008
    
            nu = 1./(wc(iw)*1.e-7)
            srayl(iw) = 1.78e-26*nu**(4. + 0.625)
            srayl(iw) = srayl(iw)*1.e-20 ! cm2
    
            do ilev = 1, nlev
                dtrl(ilev,iw) = colinc(ilev)*srayl(iw)   ! cm2
            end do
        end do
            
    end subroutine setairlay

!==============================================================================================================================================

    subroutine setaer(nlayer,alt,tau,nw,dtaer,omaer,gaer)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set aerosol properties for each specified altitude layer.  Properties    =*
        !=  may be wavelength dependent.                                             =*
        !-----------------------------------------------------------------------------*
        
        implicit none

        !     input

        integer :: nlayer                              
        integer :: nw                                 
        real, dimension(nlayer) :: alt                 
        real :: tau                                    

        !     output 

        real, dimension(nlayer,nw) :: dtaer             
        real, dimension(nlayer,nw) :: omaer            
        real, dimension(nlayer,nw) :: gaer             

        !     local

        integer :: ilay, iw
        real, dimension(nlayer) :: aer                                
        real :: omega, g, scaleh, gamma
        real :: dz, tautot, q0

        omega  = 0.622 ! single scattering albedo : wolff et al.(2010) at 258 nm
        g      = 0.88  ! asymmetry factor : mateshvili et al. (2007) at 210 nm
        scaleh = 10.   ! scale height (km)
        gamma  = 0.03  ! conrath parameter

        dtaer(:,:) = 0.
        omaer(:,:) = 0.
        gaer(:,:)  = 0.

        !     optical depth profile:

        tautot = 0.
        do ilay = 1, nlayer-1
            dz = alt(ilay+1) - alt(ilay)
            tautot = tautot + exp(gamma*(1. - exp(alt(ilay)/scaleh)))*dz
        end do
        
        q0 = tau/tautot
        do ilay = 1, nlayer-1
            dz = alt(ilay+1) - alt(ilay)
            dtaer(ilay,:) = q0*exp(gamma*(1. - exp(alt(ilay)/scaleh)))*dz
            omaer(ilay,:) = omega
            gaer(ilay,:)  = g
        end do
        
    end subroutine setaer

!==============================================================================================================================================

    subroutine setcld(nlayer,nw,dtcld,omcld,gcld)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Set cloud properties for each specified altitude layer.  Properties      =*
        !=  may be wavelength dependent.                                             =*
        !-----------------------------------------------------------------------------*
        
        implicit none

        !input

        integer :: nlayer                              
        integer :: nw                                 

        !output 

        real, dimension(nlayer,nw) :: dtcld            
        real, dimension(nlayer,nw) :: omcld            
        real, dimension(nlayer,nw) :: gcld             

        !local

        integer :: ilay, iw

        !dtcld : optical depth
        !omcld : single scattering albedo
        !gcld  : asymmetry factor

        do ilay = 1, nlayer - 1
            do iw = 1, nw - 1
            dtcld(ilay,iw) = 0.       ! no clouds for the moment
            omcld(ilay,iw) = 0.99
            gcld(ilay,iw)  = 0.85
            end do
        end do
        
    end subroutine setcld

!==============================================================================================================================================

    subroutine sphers(nlev, z, zen, dsdh, nid)

        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Calculate slant path over vertical depth ds/dh in spherical geometry.    =*
        !=  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model =*
        !=  for computing the radiation field available for photolysis and heating   =*
        !=  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)  =*
        !-----------------------------------------------------------------------------*
        !=  PARAMETERS:                                                              =*
        !=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
        !=            grid                                                           =*
        !=  Z       - REAL, specified altitude working grid (km)                  (I)=*
        !=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
        !=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
        !=            when travelling from the top of the atmosphere to layer i;     =*
        !=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
        !=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
        !=            travelling from the top of the atmosphere to layer i;          =*
        !=            NID(i), i = 0..NZ-1                                            =*
        !-----------------------------------------------------------------------------*
        
        implicit none

        ! input

        integer, intent(in) :: nlev
        real, dimension(nlev), intent(in) :: z
        real, intent(in) :: zen

        ! output

        INTEGER, intent(out) :: nid(0:nlev)
        REAL, intent(out) :: dsdh(0:nlev,nlev)

        ! more program constants

        REAL re, ze(nlev)
        REAL  dr
        real radius
        parameter (radius = 3393.)

            ! local 

        real :: pi, zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
        integer :: i, j, k, id, nlay

        REAL zd(0:nlev-1)

        !-----------------------------------------------------------------------------

        pi = acos(-1.0)
        dr = pi/180.
        zenrad = zen*dr

        ! number of layers:

        nlay = nlev - 1

        ! include the elevation above sea level to the radius of Mars:

        re = radius + z(1)

        ! correspondingly z changed to the elevation above Mars surface:

        DO k = 1, nlev
            ze(k) = z(k) - z(1)
        END DO

        ! inverse coordinate of z

        zd(0) = ze(nlev)
        DO k = 1, nlay
        zd(k) = ze(nlev - k)
        END DO

        ! initialise dsdh(i,j), nid(i)

        nid(:) = 0.
        dsdh(:,:) = 0.

        ! calculate ds/dh of every layer

        do i = 0,nlay
            rpsinz = (re + zd(i))*sin(zenrad)
    
        IF ( (zen .GT. 90.0) .AND. (rpsinz .LT. re) ) THEN
            nid(i) = -1
        ELSE

        ! Find index of layer in which the screening height lies

            id = i 
            if (zen > 90.) then
                do j = 1,nlay
                    IF( (rpsinz .LT. ( zd(j-1) + re ) ) .AND. &
                       (rpsinz .GE. ( zd(j) + re )) ) id = j
                end do
            end if
    
            do j = 1,id
                sm = 1.0
                IF (j .EQ. id .AND. id .EQ. i .AND. zen .GT. 90.0) &
                    sm = -1.0
    
                rj = re + zd(j-1)
                rjp1 = re + zd(j)
    
                dhj = zd(j-1) - zd(j)
    
                ga = rj*rj - rpsinz*rpsinz
                gb = rjp1*rjp1 - rpsinz*rpsinz

                ga = max(ga, 0.)
                gb = max(gb, 0.)
    
                IF (id.GT.i .AND. j.EQ.id) THEN
                    dsj = sqrt(ga)
                ELSE
                    dsj = sqrt(ga) - sm*sqrt(gb)
                END IF
                dsdh(i,j) = dsj/dhj
            end do
            nid(i) = id
            end if
        end do ! i = 0,nlay
        
    end subroutine sphers

!==============================================================================================================================================

    subroutine rtlink(nlev, nw, iw, ag, zen, dsdh, nid, dtrl, &
                          dagas, dtcld, omcld, gcld, dtaer, omaer, gaer, &
                          edir, edn, eup, fdir, fdn, fup)
   
        implicit none

        ! input

        integer, intent(in) :: nlev, nw, iw              
        REAL, intent(in) :: ag
        REAL, intent(in) :: zen
        REAL, intent(in) :: dsdh(0:nlev,nlev)
        INTEGER, intent(in) :: nid(0:nlev)

        REAL, intent(in) :: dtrl(nlev,nw)
        REAL, intent(in) :: dagas(nlev,nw)
        REAL, intent(in) :: dtcld(nlev,nw), omcld(nlev,nw), gcld(nlev,nw)
        REAL, intent(in) :: dtaer(nlev,nw), omaer(nlev,nw), gaer(nlev,nw)

        ! output

        REAL, intent(out) :: edir(nlev), edn(nlev), eup(nlev)
        REAL, intent(out) :: fdir(nlev), fdn(nlev), fup(nlev)

        ! local:

        REAL dt(nlev), om(nlev), g(nlev)
        REAL dtabs,dtsct,dscld,dsaer,dacld,daaer
        INTEGER i, ii
        real, parameter :: largest = 1.e+36

        ! specific two ps2str

        REAL ediri(nlev), edni(nlev), eupi(nlev)
        REAL fdiri(nlev), fdni(nlev), fupi(nlev)

        logical, save :: delta = .true.

        !_______________________________________________________________________

        ! initialize:

        do i = 1, nlev
            fdir(i) = 0.
            fup(i)  = 0.
            fdn(i)  = 0.
            edir(i) = 0.
            eup(i)  = 0.
            edn(i)  = 0.
        end do

        do i = 1, nlev - 1
            dscld = dtcld(i,iw)*omcld(i,iw)
            dacld = dtcld(i,iw)*(1.-omcld(i,iw))

            dsaer = dtaer(i,iw)*omaer(i,iw)
            daaer = dtaer(i,iw)*(1.-omaer(i,iw))

            dtsct = dtrl(i,iw) + dscld + dsaer
            dtabs = dagas(i,iw) + dacld + daaer 

            dtabs = amax1(dtabs,1./largest)
            dtsct = amax1(dtsct,1./largest)

        !     invert z-coordinate:

            ii = nlev - i
            dt(ii) = dtsct + dtabs
            om(ii) = dtsct/(dtsct + dtabs)
            IF(dtsct .EQ. 1./largest) om(ii) = 1./largest
            g(ii) = (gcld(i,iw)*dscld + &
                        gaer(i,iw)*dsaer)/dtsct
        end do

        print*,dt
        print*,om
        print*,g

        !     call rt routine:

        call ps2str(nlev, zen, ag, dt, om, g, &
                dsdh, nid, delta, &
                fdiri, fupi, fdni, ediri, eupi, edni)

        !     output (invert z-coordinate)

        do i = 1, nlev
            ii = nlev - i + 1
            fdir(i) = fdiri(ii)
            fup(i) = fupi(ii)
            fdn(i) = fdni(ii)
            edir(i) = ediri(ii)
            eup(i) = eupi(ii)
            edn(i) = edni(ii)
        end do
   
    end subroutine rtlink

!==============================================================================================================================================

    subroutine ps2str(nlev,zen,rsfc,tauu,omu,gu,&
                          dsdh, nid, delta,&
                          fdr, fup, fdn, edr, eup, edn)
   
        !-----------------------------------------------------------------------------*
        !=  PURPOSE:                                                                 =*
        !=  Solve two-stream equations for multiple layers.  The subroutine is based =*
        !=  on equations from:  Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.=*
        !=  It contains 9 two-stream methods to choose from.  A pseudo-spherical     =*
        !=  correction has also been added.                                          =*
        !-----------------------------------------------------------------------------*
        !=  PARAMETERS:                                                              =*
        !=  NLEVEL  - INTEGER, number of specified altitude levels in the working (I)=*
        !=            grid                                                           =*
        !=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
        !=  RSFC    - REAL, surface albedo at current wavelength                  (I)=*
        !=  TAUU    - REAL, unscaled optical depth of each layer                  (I)=*
        !=  OMU     - REAL, unscaled single scattering albedo of each layer       (I)=*
        !=  GU      - REAL, unscaled asymmetry parameter of each layer            (I)=*
        !=  DSDH    - REAL, slant path of direct beam through each layer crossed  (I)=*
        !=            when travelling from the top of the atmosphere to layer i;     =*
        !=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
        !=  NID     - INTEGER, number of layers crossed by the direct beam when   (I)=*
        !=            travelling from the top of the atmosphere to layer i;          =*
        !=            NID(i), i = 0..NZ-1                                            =*
        !=  DELTA   - LOGICAL, switch to use delta-scaling                        (I)=*
        !=            .TRUE. -> apply delta-scaling                                  =*
        !=            .FALSE.-> do not apply delta-scaling                           =*
        !=  FDR     - REAL, contribution of the direct component to the total     (O)=*
        !=            actinic flux at each altitude level                            =*
        !=  FUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
        !=            the total actinic flux at each altitude level                  =*
        !=  FDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
        !=            the total actinic flux at each altitude level                  =*
        !=  EDR     - REAL, contribution of the direct component to the total     (O)=*
        !=            spectral irradiance at each altitude level                     =*
        !=  EUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
        !=            the total spectral irradiance at each altitude level           =*
        !=  EDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
        !=            the total spectral irradiance at each altitude level           =*
        !-----------------------------------------------------------------------------*

        implicit none

        ! input:

        INTEGER, intent(in) :: nlev
        REAL, intent(in) :: zen, rsfc
        REAL, intent(in) :: tauu(nlev), omu(nlev), gu(nlev)
        REAL, intent(in) :: dsdh(0:nlev,nlev)
        INTEGER, intent(in) :: nid(0:nlev)
        LOGICAL, intent(in) :: delta

        ! output:

        REAL, intent(out) :: fup(nlev),fdn(nlev),fdr(nlev)
        REAL, intent(out) :: eup(nlev),edn(nlev),edr(nlev)

        ! local:

        REAL tausla(0:nlev), tauc(0:nlev)
        REAL mu2(0:nlev), mu, sum

        ! internal coefficients and matrix

        REAL lam(nlev),taun(nlev),bgam(nlev)
        REAL e1(nlev),e2(nlev),e3(nlev),e4(nlev)
        REAL cup(nlev),cdn(nlev),cuptn(nlev),cdntn(nlev)
        REAL mu1(nlev)
        INTEGER row
        REAL a(2*nlev),b(2*nlev),d(2*nlev),e(2*nlev),y(2*nlev)

        ! other:

        REAL pifs, fdn0
        REAL gi(nlev), omi(nlev), tempg
        REAL f, g, om
        REAL gam1, gam2, gam3, gam4
        real, parameter :: largest = 1.e+36
        real, parameter :: precis = 1.e-7

        ! For calculations of Associated Legendre Polynomials for GAMA1,2,3,4
        ! in delta-function, modified quadrature, hemispheric constant,
        ! Hybrid modified Eddington-delta function metods, p633,Table1.
        ! W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
        ! W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440, 
        ! uncomment the following two lines and the appropriate statements further
        ! down.
        !     REAL YLM0, YLM2, YLM4, YLM6, YLM8, YLM10, YLM12, YLMS, BETA0,
        !    >     BETA1, BETAn, amu1, subd

        REAL expon, expon0, expon1, divisr, temp, up, dn
        REAL ssfc
        INTEGER nlayer, mrows, lev

        INTEGER i, j

        ! Some additional program constants:

        real pi, dr
        REAL eps
        PARAMETER (eps = 1.E-3)
        !_______________________________________________________________________

        ! MU = cosine of solar zenith angle
        ! RSFC = surface albedo
        ! TAUU =  unscaled optical depth of each layer
        ! OMU  =  unscaled single scattering albedo
        ! GU   =  unscaled asymmetry factor
        ! KLEV = max dimension of number of layers in atmosphere
        ! NLAYER = number of layers in the atmosphere
        ! NLEVEL = nlayer + 1 = number of levels

        ! initial conditions:  pi*solar flux = 1;  diffuse incidence = 0

        pifs = 1.      
        fdn0 = 0.

        nlayer = nlev - 1

        pi = acos(-1.)
        dr = pi/180.
        mu = COS(zen*dr)

        !************* compute coefficients for each layer:
        ! GAM1 - GAM4 = 2-stream coefficients, different for different approximations
        ! EXPON0 = calculation of e when TAU is zero
        ! EXPON1 = calculation of e when TAU is TAUN
        ! CUP and CDN = calculation when TAU is zero
        ! CUPTN and CDNTN = calc. when TAU is TAUN
        ! DIVISR = prevents division by zero

        do j = 0, nlev
            tauc(j) = 0.
            tausla(j) = 0.
            mu2(j) = 1./SQRT(largest)
        end do

        IF (.NOT. delta) THEN
            DO i = 1, nlayer
                gi(i) = gu(i)
                omi(i) = omu(i)
                taun(i) = tauu(i)
            END DO
        ELSE 

            ! delta-scaling. Have to be done for delta-Eddington approximation, 
            ! delta discrete ordinate, Practical Improved Flux Method, delta function,
            ! and Hybrid modified Eddington-delta function methods approximations

            DO i = 1, nlayer
                f = gu(i)*gu(i)
                gi(i) = (gu(i) - f)/(1 - f)
                omi(i) = (1 - f)*omu(i)/(1 - omu(i)*f)       
                taun(i) = (1 - omu(i)*f)*tauu(i)
            END DO
        END IF

        ! calculate slant optical depth at the top of the atmosphere when zen>90.
        ! in this case, higher altitude of the top layer is recommended.

        IF (zen .GT. 90.0) THEN
            IF (nid(0) .LT. 0) THEN
                tausla(0) = largest
            ELSE
                sum = 0.0
                DO j = 1, nid(0)
                    sum = sum + 2.*taun(j)*dsdh(0,j)
                END DO
                tausla(0) = sum 
            END IF
        END IF

        DO 11, i = 1, nlayer
            g = gi(i)
            om = omi(i)
            tauc(i) = tauc(i-1) + taun(i)

            ! stay away from 1 by precision.  For g, also stay away from -1

            tempg = AMIN1(abs(g),1. - precis)
            g = SIGN(tempg,g)
            om = AMIN1(om,1.-precis)

            ! calculate slant optical depth
                
            IF (nid(i) .LT. 0) THEN
                tausla(i) = largest
            ELSE
                sum = 0.0
                DO j = 1, MIN(nid(i),i)
                    sum = sum + taun(j)*dsdh(i,j)
                END DO
                DO j = MIN(nid(i),i)+1,nid(i)
                    sum = sum + 2.*taun(j)*dsdh(i,j)
                END DO
                tausla(i) = sum 
                IF (tausla(i) .EQ. tausla(i-1)) THEN
                    mu2(i) = SQRT(largest)
                ELSE
                    mu2(i) = (tauc(i)-tauc(i-1))/(tausla(i)-tausla(i-1)) 
                    mu2(i) = SIGN( AMAX1(ABS(mu2(i)),1./SQRT(largest)),&
                              mu2(i) )
                END IF
            END IF

            !** the following gamma equations are from pg 16,289, Table 1
            !** save mu1 for each approx. for use in converting irradiance to actinic flux

            ! Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

            !c       gam1 =  (7. - om*(4. + 3.*g))/4.
            !c       gam2 = -(1. - om*(4. - 3.*g))/4.
            !c       gam3 = (2. - 3.*g*mu)/4.
            !c       gam4 = 1. - gam3
            !c       mu1(i) = 0.5

            !* quadrature (Liou, 1973, JAS, 30, 1303-1326; 1974, JAS, 31, 1473-1475):

            !c          gam1 = 1.7320508*(2. - om*(1. + g))/2.
            !c          gam2 = 1.7320508*om*(1. - g)/2.
            !c          gam3 = (1. - 1.7320508*g*mu)/2.
            !c          gam4 = 1. - gam3
            !c          mu1(i) = 1./sqrt(3.)
        
            !* hemispheric mean (Toon et al., 1089, JGR, 94, 16287):

            gam1 = 2. - om*(1. + g)
            gam2 = om*(1. - g)
            gam3 = (2. - g*mu)/4.
            gam4 = 1. - gam3
            mu1(i) = 0.5

            !* PIFM  (Zdunkovski et al.,1980, Conrib.Atmos.Phys., 53, 147-166):
            !c         GAM1 = 0.25*(8. - OM*(5. + 3.*G))
            !c         GAM2 = 0.75*OM*(1.-G)
            !c         GAM3 = 0.25*(2.-3.*G*MU)
            !c         GAM4 = 1. - GAM3
            !c         mu1(i) = 0.5

            !* delta discrete ordinates  (Schaller, 1979, Contrib.Atmos.Phys, 52, 17-26):
            !c         GAM1 = 0.5*1.7320508*(2. - OM*(1. + G))
            !c         GAM2 = 0.5*1.7320508*OM*(1.-G)
            !c         GAM3 = 0.5*(1.-1.7320508*G*MU)
            !c         GAM4 = 1. - GAM3
            !c         mu1(i) = 1./sqrt(3.)

            !* Calculations of Associated Legendre Polynomials for GAMA1,2,3,4
            !* in delta-function, modified quadrature, hemispheric constant,
            !* Hybrid modified Eddington-delta function metods, p633,Table1.
            !* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
            !* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440
            !c      YLM0 = 2.
            !c      YLM2 = -3.*G*MU
            !c      YLM4 = 0.875*G**3*MU*(5.*MU**2-3.)
            !c      YLM6=-0.171875*G**5*MU*(15.-70.*MU**2+63.*MU**4)
            !c     YLM8=+0.073242*G**7*MU*(-35.+315.*MU**2-693.*MU**4
            !c    *+429.*MU**6)
            !c     YLM10=-0.008118*G**9*MU*(315.-4620.*MU**2+18018.*MU**4
            !c    *-25740.*MU**6+12155.*MU**8)
            !c     YLM12=0.003685*G**11*MU*(-693.+15015.*MU**2-90090.*MU**4
            !c    *+218790.*MU**6-230945.*MU**8+88179.*MU**10)
            !c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
            !c      YLMS=0.25*YLMS
            !c      BETA0 = YLMS
            !c
            !c         amu1=1./1.7320508
            !c      YLM0 = 2.
            !c      YLM2 = -3.*G*amu1
            !c      YLM4 = 0.875*G**3*amu1*(5.*amu1**2-3.)
            !c      YLM6=-0.171875*G**5*amu1*(15.-70.*amu1**2+63.*amu1**4)
            !c     YLM8=+0.073242*G**7*amu1*(-35.+315.*amu1**2-693.*amu1**4
            !c    *+429.*amu1**6)
            !c     YLM10=-0.008118*G**9*amu1*(315.-4620.*amu1**2+18018.*amu1**4
            !c    *-25740.*amu1**6+12155.*amu1**8)
            !c     YLM12=0.003685*G**11*amu1*(-693.+15015.*amu1**2-90090.*amu1**4
            !c    *+218790.*amu1**6-230945.*amu1**8+88179.*amu1**10)
            !c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
            !c      YLMS=0.25*YLMS
            !c      BETA1 = YLMS
            !c
            !c         BETAn = 0.25*(2. - 1.5*G-0.21875*G**3-0.085938*G**5
            !c    *-0.045776*G**7)


            !* Hybrid modified Eddington-delta function(Meador and Weaver,1980,JAS,37,630):
            !c         subd=4.*(1.-G*G*(1.-MU))
            !c         GAM1 = (7.-3.*G*G-OM*(4.+3.*G)+OM*G*G*(4.*BETA0+3.*G))/subd
            !c         GAM2 =-(1.-G*G-OM*(4.-3.*G)-OM*G*G*(4.*BETA0+3.*G-4.))/subd
            !c         GAM3 = BETA0
            !c         GAM4 = 1. - GAM3
            !c         mu1(i) = (1. - g*g*(1.- mu) )/(2. - g*g)

            !*****
            !* delta function  (Meador, and Weaver, 1980, JAS, 37, 630):
            !c         GAM1 = (1. - OM*(1. - beta0))/MU
            !c         GAM2 = OM*BETA0/MU
            !c         GAM3 = BETA0
            !c         GAM4 = 1. - GAM3
            !c         mu1(i) = mu
            !*****
            !* modified quadrature (Meador, and Weaver, 1980, JAS, 37, 630):
            !c         GAM1 = 1.7320508*(1. - OM*(1. - beta1))
            !c         GAM2 = 1.7320508*OM*beta1
            !c         GAM3 = BETA0
            !c         GAM4 = 1. - GAM3
            !c         mu1(i) = 1./sqrt(3.)

            !* hemispheric constant (Toon et al., 1989, JGR, 94, 16287):
            !c         GAM1 = 2.*(1. - OM*(1. - betan))
            !c         GAM2 = 2.*OM*BETAn
            !c         GAM3 = BETA0
            !c         GAM4 = 1. - GAM3
            !c         mu1(i) = 0.5

            !*****

            !* lambda = pg 16,290 equation 21
            !* big gamma = pg 16,290 equation 22
            !* if gam2 = 0., then bgam = 0. 

            lam(i) = sqrt(gam1*gam1 - gam2*gam2)

            IF (gam2 .NE. 0.) THEN
                bgam(i) = (gam1 - lam(i))/gam2
            ELSE
                bgam(i) = 0.
            END IF

            expon = EXP(-lam(i)*taun(i))

            !* e1 - e4 = pg 16,292 equation 44
            
            e1(i) = 1. + bgam(i)*expon
            e2(i) = 1. - bgam(i)*expon
            e3(i) = bgam(i) + expon
            e4(i) = bgam(i) - expon

            !* the following sets up for the C equations 23, and 24
            !* found on page 16,290
            !* prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
            !* which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

            expon0 = EXP(-tausla(i-1))
            expon1 = EXP(-tausla(i))
                
            divisr = lam(i)*lam(i) - 1./(mu2(i)*mu2(i))
            temp = AMAX1(eps,abs(divisr))
            divisr = SIGN(temp,divisr)

            up = om*pifs*((gam1 - 1./mu2(i))*gam3 + gam4*gam2)/divisr
            dn = om*pifs*((gam1 + 1./mu2(i))*gam4 + gam2*gam3)/divisr
        
            ! cup and cdn are when tau is equal to zero
            ! cuptn and cdntn are when tau is equal to taun

            cup(i) = up*expon0
            cdn(i) = dn*expon0
            cuptn(i) = up*expon1
            cdntn(i) = dn*expon1

        11 CONTINUE

        !***************** set up matrix ******
        !ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

        ssfc = rsfc*mu*EXP(-tausla(nlayer))*pifs

        !MROWS = the number of rows in the matrix

        mrows = 2*nlayer     
                
        !the following are from pg 16,292  equations 39 - 43.
        !set up first row of matrix:

        i = 1
        a(1) = 0.
        b(1) = e1(i)
        d(1) = -e2(i)
        e(1) = fdn0 - cdn(i)

        row=1

        !set up odd rows 3 thru (MROWS - 1):

        i = 0
        DO 20, row = 3, mrows - 1, 2
            i = i + 1
            a(row) = e2(i)*e3(i) - e4(i)*e1(i)
            b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
            d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
            e(row) = e3(i)*(cup(i + 1) - cuptn(i)) + &
                e1(i)*(cdntn(i) - cdn(i + 1))
        20  CONTINUE

        !* set up even rows 2 thru (MROWS - 2): 

        i = 0
        DO 30, row = 2, mrows - 2, 2
            i = i + 1
            a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
            b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
            d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
            e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) - &
                (cdn(i + 1) - cdntn(i))*e4(i + 1)
        30  CONTINUE

        !* set up last row of matrix at MROWS:

        row = mrows
        i = nlayer
        
        a(row) = e1(i) - rsfc*e3(i)
        b(row) = e2(i) - rsfc*e4(i)
        d(row) = 0.
        e(row) = ssfc - cuptn(i) + rsfc*cdntn(i)

        !* solve tri-diagonal matrix:

        call tridiag(a, b, d, e, y, mrows)

        !**** unfold solution of matrix, compute output fluxes:

        row = 1 
        lev = 1
        j = 1
        
        !* the following equations are from pg 16,291  equations 31 & 32

        fdr(lev) = EXP( -tausla(0) )
        edr(lev) = mu * fdr(lev)
        edn(lev) = fdn0
        eup(lev) =  y(row)*e3(j) - y(row + 1)*e4(j) + cup(j)
        fdn(lev) = edn(lev)/mu1(lev)
        fup(lev) = eup(lev)/mu1(lev)

        DO 60, lev = 2, nlayer + 1
            fdr(lev) = EXP(-tausla(lev-1))
            edr(lev) =  mu *fdr(lev)
            edn(lev) =  y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
            eup(lev) =  y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
            fdn(lev) = edn(lev)/mu1(j)
            fup(lev) = eup(lev)/mu1(j)

            row = row + 2
            j = j + 1
        60  CONTINUE

    end subroutine ps2str

!==============================================================================================================================================

    subroutine tridiag(a,b,c,r,u,n)

        !_______________________________________________________________________
        ! solves tridiagonal system.  From Numerical Recipies, p. 40
        !_______________________________________________________________________

        IMPLICIT NONE

        ! input:

        INTEGER n
        REAL a, b, c, r
        DIMENSION a(n),b(n),c(n),r(n)

        ! output:

        REAL u
        DIMENSION u(n)

        ! local:

        INTEGER j

        REAL bet, gam
        DIMENSION gam(n)
        !_______________________________________________________________________


        IF (b(1) .EQ. 0.) STOP 1001
        bet   = b(1)
        u(1) = r(1)/bet
        DO 11, j = 2, n   
            gam(j) = c(j - 1)/bet
            bet = b(j) - a(j)*gam(j)
            IF (bet .EQ. 0.) STOP 2002 
            u(j) = (r(j) - a(j)*u(j - 1))/bet
        11 CONTINUE
        DO 12, j = n - 1, 1, -1  
            u(j) = u(j) - gam(j + 1)*u(j + 1)
        12 CONTINUE
        !   
        
    end subroutine tridiag

END MODULE photolysis
