import h5py
import numpy as np
import sys,os
from pchempy import *


#Set of routines to work with the photolysis scheme of the photochemical code

def calculate_photolysis_rates(h,P,T,N,gasID,isoID,sza=0.0,tauvis=0.0,dist_sun=1.5):
    
    """
        FUNCTION NAME : calculate_photolysis_rates()
        
        DESCRIPTION : Calculate the photolysis rates in each layer of the atmosphere
        
        INPUTS :
        
            h(nlay) :: Altitude of each layer (m)
            p(nlay) :: Pressure of each layer (Pa)
            T(nlay) :: Temperature of each layer (K)
            N(nlay,ngas) :: Number density of each gas in each layer (m-3)
            gasID(ngas) :: Gas ID
            isoID(ngas) :: Isotope ID
            
        OPTIONAL INPUTS:
         
            sza :: Solar zenith angle
            tauvis :: Dust optical depth
            dist_sun :: Distance to parent star (AU)
        
        OUTPUTS :

            rrates
        
        CALLING SEQUENCE:
        
            calibrate_diffor_acs(DataDir,Observation,DiffOrder)
        
        MODIFICATION HISTORY : Juan Alday (15/12/2023)
        
    """
    
    #Reading the ID of the gases that can be photolysed in the model
    gasID_phot_inc = pchemf.photolysis.gasid_phot_inc
    isoID_phot_inc = pchemf.photolysis.isoid_phot_inc
    n_phot_inc = pchemf.photolysis.n_phot_inc
    
    
    #Selecting which gases of the ones in the input can be photolysed
    gasIDact = []
    isoIDact = []
    n_photolysis = 0
    for i in range(len(gasID)):
        for j in range(len(gasID_phot_inc)):
            if( (gasID[i]==gasID_phot_inc[j]) & (isoID[i]==isoID_phot_inc[j]) ):
                gasIDact.append(gasID[i])
                isoIDact.append(isoID[i])
                n_photolysis = n_photolysis + n_phot_inc[j]
                
    gasIDact = np.array(gasIDact)
    isoIDact = np.array(isoIDact)
    
    #Calculating the column density of each layer
    delz = h[1] - h[0]   #Assuming all layers have the same height
    cdens = N * delz     #m-2
    
    #Calculating the photolysis rates
    rtype, ns, sID, sISO, sf, npr, pID, pISO, pf, rrates = \
        pchemf.photolysis.photolysis_online(nlay=len(h), ngas=len(gasID), ngas_phot=len(gasIDact), \
        gasid=gasID, isoid=isoID, gasid_phot=gasIDact, isoid_phot=isoIDact, \
        h=h, t=T, cdens=cdens, sza=sza, tau=tauvis, dist_sun=dist_sun, n_phot=n_photolysis)
        
    return rtype, ns, sID, sISO, sf, npr, pID, pISO, pf, rrates
        
    
    


    