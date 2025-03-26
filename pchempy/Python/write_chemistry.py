import numpy as np
import sys,os
from pchempy import *

#Set of routines to write the photochemical model online

##################################################################################################
##################################################################################################
#                                     CHEMICAL NETWORK
##################################################################################################
##################################################################################################


def write_chemistry_fortran(runname,gasID_phot,isoID_phot,reactions,reactions_c13=None):
    
    """
        FUNCTION NAME : write_chemistry_fortran()
        
        DESCRIPTION : Routine to write the fortran program with the reaction network.
        
        INPUTS :
        
            gasID_phot(ngas_phot) :: ID of the gases to be photolysed
            isoID_phot(ngas_phot) :: ID of the isotopes to be photolysed
            reactions(nr) :: Array containing the reaction numbers to include in the network
            
        OPTIONAL INPUTS:
        
        OUTPUTS :
        
        CALLING SEQUENCE:
        
            write_chemistry_fortran(ngas_phot,gasID_phot,isoID_phot,reactions)
        
        MODIFICATION HISTORY : Juan Alday (15/12/2023)
        
    """
    
    from pchempy.dict import reaction_info,gas_info
    
    #Reading the ID of the gases that can be photolysed in the model
    gasID_phot_inc = pchemf.photolysis.gasid_phot_inc
    isoID_phot_inc = pchemf.photolysis.isoid_phot_inc
    n_phot_inc = pchemf.photolysis.n_phot_inc
    
    #Checking that all gases whose photolysis is inteded to be included
    #can indeed be photolised in the photochemical model
    #Counting the number of photolysis reactions
    ngas_phot = len(gasID_phot)
    n_phot = 0
    for i in range(ngas_phot):
        
        iinc = False
        for j in range(len(gasID_phot_inc)):

            if( (gasID_phot[i]==gasID_phot_inc[j]) & (isoID_phot[i]==isoID_phot_inc[j]) ):
                iinc = True
                n_phot = n_phot + n_phot_inc[j]
                
        if iinc==False:
            print('gasID = ',gasID_phot[i],'isoID = ',isoID_phot[i])
            sys.exit('error :: There is a gas whose photolysis is not included in the model')
            
    
    #Counting the number of other reactions
    nreactionsx = len(reactions)
    nreactions = nreactionsx + n_phot
    
    nreactionsx_c13 = 0
    if reactions_c13 is not None:
        nreactionsx_c13 = len(reactions_c13)
        nreactions = nreactions + nreactionsx_c13
    
    
    #Writing the first part of the module
    fortran_code = f"""
MODULE {runname}_chemical_network

    integer, save :: ngas_phot = {ngas_phot}         !number of gases to be photolysed

CONTAINS

    """
    
    #Writing the number of reactions 
    loadreactions_code = f"""
!============================================================================================================================

    subroutine number_reactions(photolysis_inc,nreactions,n_photolysis)
    
        !Function to calculate the number of reactions in the chemical network

        logical, intent(in) :: photolysis_inc

        integer, intent(out) :: nreactions,n_photolysis
        
        
        !Determining the number of reactions depending on the flags
        !#####################################################################################

        !Photolysis
        n_photolysis = 0
        if(photolysis_inc) n_photolysis = {n_phot}

        !Chemistry
        nreactions = {nreactionsx+nreactionsx_c13} + n_photolysis
    
    end subroutine number_reactions
    
    """
    
    #Writing the reaction network
    loadreactions_code += f"""
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
    """
    
    for ir in range(nreactionsx):
        
        #Getting the reaction name from reaction dictionary
        reaction_name = reaction_info[str(reactions[ir])]["name"]
        nr = str(reactions[ir]).zfill(4)
        
        thirdbody = reaction_info[str(reactions[ir])]["thirdbody"]
        if thirdbody=='':
            thirdbody = 'dens'
        
        
        loadreactions_code += f"""
        !{reaction_name}
        i0 = i0 + 1
        call reaction{nr}(nlay,p(:),t(:),{thirdbody}(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

        """
        
        
    if reactions_c13 is not None:
        
        loadreactions_code += f"""
        !13-C reactions
        !##############################################
        """
        
        for ir in range(nreactionsx_c13):
            
            #Getting the reaction name from reaction dictionary
            reaction_name = reaction_c13_info[str(reactions_c13[ir])]["name"]
            nr = str(reactions_c13[ir]).zfill(4)
            
            thirdbody = reaction_c13_info[str(reactions_c13[ir])]["thirdbody"]
            if thirdbody=='':
                thirdbody = 'dens'
            
            
            loadreactions_code += f"""
        !{reaction_name}
        i0 = i0 + 1
        call reaction{nr}_13c(nlay,p(:),t(:),{thirdbody}(:),&
            rrates(:,i0),rtype(i0),ns(i0),sID(:,i0),sISO(:,i0),sf(:,i0),npr(i0),pID(:,i0),pISO(:,i0),pf(:,i0),ref)

            """
        
        
    loadreactions_code += f"""
        end subroutine load_reactions
    """
        
    
    #Writing the photolysis reactions
    loadphotolysis_code = f"""
!============================================================================================================================

    subroutine load_reactions_photolysis(ngas_phot,gasID_phot,isoID_phot)

        !Loading the photolysis reactions 

        !Inputs
        integer, intent(in) :: ngas_phot  !Number of gases to be photolysed

        !Outputs
        integer, intent(out) :: gasID_phot(ngas_phot),isoID_phot(ngas_phot)  !Gases to be photolysed
            
    """
    
    for igas in range(ngas_phot):
        
        if isoID_phot[igas]==0:
            gasname = gas_info[str(gasID_phot[igas])]["name"]
        else:
            gasname = gas_info[str(gasID_phot[igas])]["isotope"][str(isoID_phot[igas])]["name"]
        
        loadphotolysis_code += f"""
                gasID_phot({igas+1}) = {gasID_phot[igas]}   !{gasname}
                isoID_phot({igas+1}) = {isoID_phot[igas]}   
        
        """
    
    loadphotolysis_code += f"""
    end subroutine load_reactions_photolysis
                
    """
    
    
    
    #Writing the last part of the module
    final_code = f"""
    
END MODULE {runname}_chemical_network
    """
    
    fortran_code += loadreactions_code
    fortran_code += loadphotolysis_code
    fortran_code += final_code
    
    # Open the file in write mode and write the Fortran code
    with open(runname+'_chemical_network.f90', 'w') as f:
        f.write(fortran_code)
            
            
