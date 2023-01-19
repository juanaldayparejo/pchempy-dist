######################################################################
#                  1D photochemical model for Mars
######################################################################

import sys,os,time
from pchempy import *

#Defining inputs
######################################################################

planet = 'Mars'

#Input/Output file
#--------------------------

runname = input('runname :: ')

#Flags
#--------------

diffusion = True
chemistry = True
nitrogen_chemistry = False
c13_chemistry = True
o18_chemistry = False


#Photolysis
#--------------------------

sza = 60.         #Solar zenith angle (degrees)
dist_sun = 1.52    #AU
tauvis = 0.0               # Dust optical depth

#Fixing species
#--------------------------

#Fixing the H2O profile from 0 to 50 km
gasID_fix = np.array([1])   #H2O
isoID_fix = np.array([0])   #H2O

hmin_fix = np.array([0.])   #Altitude boundaries at which to fix the profile
hmax_fix = np.array([50.0]) 


# get the start time
st = time.time()


#Initialising the atmospheric properties from the HDF5 file (last step)
########################################################################

#Reading the HDF5 file
hlay,Tlay,gasID,isoID,Nlay1,ztime = read_hdf5(runname)

nlay = Nlay1.shape[0]
ngas = Nlay1.shape[1]
nt = len(ztime)

#Using the last timestep as the initial condition for the model
Nlay = np.zeros((nlay,ngas),dtype='float64')
Nlay[:,:] = Nlay1[:,:,nt-1]

#Calculating VMRs and total density
N0lay = np.sum(Nlay,axis=1)
VMRlay = np.transpose(np.transpose(Nlay)/N0lay)

#Calculating the pressure
Play = N0lay * phys_const['k_B'] * Tlay


#Performing initial calculations
#########################################################################

#Calculate the molecular weight of each species (g mol-1 or uma)
mmol = calc_mmol(gasID,isoID)

mmean = calc_mmean(Nlay,mmol)
print('Total mass :: ',np.sum(N0lay*mmean))

#Reading the coefficients for molecular diffusion
A,s,B = read_moldiff_params(gasID,isoID,planet)

#Reading the boundary conditions from the dictionary
typeubc,valueubc = calc_upper_bc(gasID,isoID,planet)
typelbc,valuelbc = calc_lower_bc(gasID,isoID,planet)

#Locating the gases to be fixed
##########################################################################

fix_species = np.zeros((nlay,ngas),dtype='int32')
for igasx in range(len(gasID_fix)):

    il = np.where( (hlay>=hmin_fix[igasx]*1000.) & (hlay<=hmax_fix[igasx]*1000.) )
    il = il[0]

    ifix = 0
    for igas in range(ngas):
        if ( (gasID_fix[igasx]==gasID[igas]) & (isoID_fix[igasx]==isoID[igas]) ):

            fix_species[il,igas] = 1
            ifix = 1
        
    if ifix==0:
        sys.exit('error :: The gas species to be fixed cannot be found in gas list')


#Iterating over different timescales to find convergence
#########################################################################

#Initialising some other parameters
Nold = np.zeros(Nlay.shape)
Nold[:,:] = Nlay[:,:]
N0layold = np.zeros(N0lay.shape)
N0layold[:] = N0lay[:]

  
#We make the profiles to converge at several timescales 
tscale = np.arange(-3.,8.,1.)
ntscale = len(tscale)

xtime = 0.0
for itscale in range(ntscale):

    dt = 10.**(1.0*tscale[itscale])
    dtmin = dt
    dtmax = dt * 10.
    #dtmax = dt

    converged=False
    iter = 0
    while converged==False:

        print('Converging in timescale :: 10**'+str(tscale[itscale])+' s')
        
        #Calling the fortran model (it will stop after X iterations even if it not converged, this way we keep better track of the convergence)
        Nnew,ztime,dtcur = pchemf.converge_model_mars_rosenbrock(hlay,Play,Tlay,Nlay,gasID,isoID,mmol,A,s,B,typelbc,valuelbc,typeubc,valueubc,fix_species,sza,tauvis,dist_sun,dtmin,dtmax,chemistry,nitrogen_chemistry,c13_chemistry,o18_chemistry)

        #Updating the densities of each species
        N0lay = np.sum(Nnew,axis=1)
        VMRlay = np.transpose(np.transpose(Nnew)/N0lay)
        Nlay[:,:] = Nnew[:,:]

        #Assessing whether convergence has been found for this timescale
        tolt = 0.01
        if( (dtcur>=(1.0-tolt)*dtmax) ):
            converged=True
        else:
            dtmin = dtcur

        #Maximum numbver of iterations in a given timescale
        if iter==20:
            converged=True

        iter = iter + 1


        #WRITING OUTPUTS
        #############################################################################

        #Writing data to HDF5 file
        hf = h5py.File(runname+'.h5','r+')
        hf['N'].resize((hf['N'].shape[2] + 1), axis=2)
        hf['N'][:,:,hf['N'].shape[2]-1] = Nlay[:,:]
        hf['time'].resize((hf['time'].shape[1] + 1), axis=1)
        hf['time'][:,hf['time'].shape[1] - 1] = np.array([ztime])[:,np.newaxis]
        hf.close()

# get the end time
et = time.time()

# get the execution time
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
