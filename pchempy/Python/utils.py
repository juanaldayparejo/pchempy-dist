import h5py
import numpy as np
from pchempy import *

'''
Set of routines needed to run the 1D photochemical model

read_hdf5 :: Read the HDF5 file with the correct format
write_hdf5 :: Write the HDF5 file with the correct format

calc_mmol :: Calculate the molecular weight of the gases in the atmosphere
calc_mmean :: Calculate the mean molecular weight of the atmosphere

calc_upper_bc :: Calculate the upper boundary conditions from the dictionary
calc_lower_bc :: Calculate the lower boundary conditions from the dictionary


'''


######################################################################################

def read_hdf5(runname):
    '''
    Function to read the output file of a simulation

    Outputs
    -------

    h(nlay) :: Altitude (m)
    T(nlay) :: Temperature (K)
    gasID(ngas) :: ID of the gases
    isoID(ngas) :: ID of the isotopes
    time(nt) :: Sol number 
    N(nlay,ngas,nt) :: Number density of each species (m-3)
    '''

    f = h5py.File(runname+'.h5','r')

    h = np.array(f['h'])
    gasID = np.array(f['gasID'],dtype='int32')
    isoID = np.array(f['isoID'],dtype='int32')
    T = np.array(f['T'])
    N = np.array(f['N'])
    time = np.array(f['time'])

    time = np.resize(time,(time.shape[1]))


    return h,T,gasID,isoID,N,time

######################################################################################

def calc_mmol(gasID,isoID):

    '''
    Routine to calculate the molecular weight of each species (g mol-1)

    Inputs
    ------
    
    gasID(ngas) :: Gas ID
    isoID(ngas) :: Isotopologue ID

    '''

    ngas = len(gasID)
    mmol = np.zeros(ngas)
    for i in range(ngas):

        if isoID[i]!=0:
            mmol1 = gas_info[str(gasID[i])]['isotope'][str(isoID[i])]['mass']
        else:
            mmol1 = gas_info[str(gasID[i])]['mmw']

        mmol[i] = mmol1

    return mmol

######################################################################################

def calc_mmean(vmr,mmol):

    '''
    Routine to calculate the mean molecular weight of the atmosphere at each level (g mol-1)

    Inputs
    ------
    
    vmr(nh,ngas) :: volume mixing ratio of each species at each altitude
    mmol(ngas) :: Molecular weight of each of the species included

    '''

    mmean = np.sum(vmr * mmol,axis=1)/np.sum(vmr,axis=1)

    return mmean

######################################################################################

def read_moldiff_params(gasID,isoID,planet):
    '''
    Function to read the parameters defining the molecular diffusion coefficient
    from the python dictionary (A,s).

    The molecular diffusion coefficient is then calculated using:

        D_i = A_i * temp**s_i / numdens

    The thermal diffusion coefficient B is also read in this function.

    '''

    ngas = len(gasID)

    A = np.zeros(ngas) 
    s = np.zeros(ngas)
    B = np.zeros(ngas)

    for i in range(ngas):

        if diffusion_coeff[planet].get(str(gasID[i])) is not None:
            A[i] = diffusion_coeff[planet][str(gasID[i])]['A']
            s[i] = diffusion_coeff[planet][str(gasID[i])]['s']
            B[i] = diffusion_coeff[planet][str(gasID[i])]['Btherm']
        else:
            A[i] = 1.0e17
            s[i] = 0.75
            B[i] = 0.0

    return A,s,B


######################################################################################

def calc_upper_bc(gasID,isoID,planet):

    '''
    Routine to read the upper boundary conditions from the dictionary.

    If a given species is not present in the dictionary then it is assumed that
    the upper boundary condition is given by a fixed flux of 0.0

    Inputs
    ------

    gasID :: Gas ID
    isoID :: Isotope ID
    planet :: Planet name (e.g., 'Mars','Venus')

    '''

    ngas = len(gasID)

    type = np.zeros(ngas,dtype='int32')
    value = np.zeros(ngas)

    for i in range(ngas):

        if upper_bc[planet].get(str(gasID[i])) is not None:
            type[i] = upper_bc[planet][str(gasID[i])]['type']
            value[i] = upper_bc[planet][str(gasID[i])]['value']
        else:
            type[i] = 2
            value[i] = 0.0

    return type,value

######################################################################################

def calc_lower_bc(gasID,isoID,planet):

    '''
    Routine to read the lower boundary conditions from the dictionary.

    If a given species is not present in the dictionary then it is assumed that
    the lower boundary condition is given by a fixed flux of 0.0

    Inputs
    ------

    gasID :: Gas ID
    isoID :: Isotope ID
    planet :: Planet name (e.g., 'Mars','Venus')

    '''

    ngas = len(gasID)

    type = np.zeros(ngas,dtype='int32')
    value = np.zeros(ngas)

    for i in range(ngas):

        if lower_bc[planet].get(str(gasID[i])) is not None:
            type[i] = lower_bc[planet][str(gasID[i])]['type']
            value[i] = lower_bc[planet][str(gasID[i])]['value']
        else:
            type[i] = 2
            value[i] = 0.0

    return type,value

######################################################################################