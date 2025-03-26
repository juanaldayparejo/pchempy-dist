import h5py
import numpy as np
import sys,os
from pchempy import *

'''
Set of routines needed to run the 1D photochemical model

read_hdf5 :: Read the HDF5 file with the correct format
write_hdf5 :: Write the HDF5 file with the correct format

calc_mmol :: Calculate the molecular weight of the gases in the atmosphere
calc_mmean :: Calculate the mean molecular weight of the atmosphere

calc_upper_bc :: Calculate the upper boundary conditions from the dictionary
calc_lower_bc :: Calculate the lower boundary conditions from the dictionary

read_profiles_mars :: Read input profiles for the atmosphere of Mars

ini_layer :: Given some vertical profiles, initialise the properties of each layer

adjust_vmr :: Adjust the VMRs so that they add up to 1

file_lines :: Read the number of lines in a file

'''

######################################################################################

def pchempy_path():
    import os
    pchempy_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')
    return pchempy_path

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
    time(nt) :: Time of simulation in each period (seconds)
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

def write_ini_hdf5(runname,h,T,gasID,isoID,N):
    '''
    Function to write the initial HDF5 file for a simulation

    Inputs
    -------

    runname :: Simulation run name
    h(nh) :: Altitude (km)
    T(nh) :: Temperature (K)
    gasID(ngas) :: Gas IDs
    isoID(ngas) :: Isotope IDs
    N(nlay,ngas) :: Number density of each species (m-3)
    '''

    nlay = len(h)
    ngas = len(gasID)

    hf = h5py.File(runname+'.h5','w')
    hf.create_dataset('h', data=h)
    hf.create_dataset('T', data=T)
    hf.create_dataset('gasID', data=gasID)
    hf.create_dataset('isoID', data=isoID)
    hf.create_dataset('N', data=N[:,:, np.newaxis], maxshape=(nlay,ngas,None))
    ztime = np.array([0.])
    hf.create_dataset('time', data=ztime[:,np.newaxis], maxshape=(1,None))
    hf.close()

#############################################################################

def read_gasname(gasID,isoID):
    '''
    Routine to read the name of the gases from the python dictionary

    Inputs
    ------
    gasID(ngas) :: Gas ID
    isoID(ngas) :: Isotope ID
    '''

    ngas = len(gasID)
    gasname = ['']*ngas
    for i in range(ngas):

        if isoID[i]!=0:
            gasname[i] = gas_info[str(gasID[i])]['isotope'][str(isoID[i])]['name']
        else:
            gasname[i] = gas_info[str(gasID[i])]['name']

    return gasname

#############################################################################

def read_gaslabel(gasID,isoID):
    '''
    Routine to read the name of the gases from the python dictionary
    with the LaTeX notation to include in plots

    Inputs
    ------
    gasID(ngas) :: Gas ID
    isoID(ngas) :: Isotope ID
    '''

    ngas = len(gasID)
    gasname = ['']*ngas
    for i in range(ngas):

        if isoID[i]!=0:
            gasname[i] = gas_info[str(gasID[i])]['isotope'][str(isoID[i])]['label']
        else:
            gasname[i] = gas_info[str(gasID[i])]['label']

    return gasname

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

def read_profiles_mars(datasrc=srcpath+'Data/Profiles/Mars/',nitrogen_chemistry=False,MakePlot=False):
    '''

    Function to initialise the atmospheric profiles for a Mars simulation

    Outputs
    --------

    zh(nh) :: Altitude (m)
    temp(nh) :: Temperature (K)
    press(nh) :: Pressure (Pa)
    numdens(nh) :: Number density (m-3)
    vmr(nh,ngas) :: Volume mixing ratio of each species
    mmol(ngas) :: Molecular weight of each species (amu)

    '''
    
    filename1 = datasrc+'atmosfera_LMD_may.dat'
    filename2 = datasrc+'atmosfera_LMD_min.dat'
    filename3 = datasrc+'atmosfera_LMD_nitr.dat'

    nlines = file_lines(filename1)


    #Reading the first file
    nh = nlines - 2     #Up to 200 km
    f = open(filename1,'r')
    header = f.readline().split()

    if nitrogen_chemistry==True:
        ngas = 18
        gasID = np.zeros(ngas,dtype='int32')   #ID of the gases to include
        isoID = np.zeros(ngas,dtype='int32')   #ID of the isotopes to include (0 is no differentiation between isotopes)
        gasnames = ['$^{12}$CO$_2$','Ar','N$_2$','O$_2$','$^{12}$CO','O','H$_2$','H','OH','HO$_2$','H$_2$O','H$_2$O$_2$','O(1D)','O$_3$','N','NO','NO$_2$','N(2D)']
        gasID[:] = [2,76,22,7,5,45,39,48,13,44,1,25,133,3,134,8,10,135]
        isoID[:] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    else:
        ngas = 14
        gasID = np.zeros(ngas,dtype='int32')   #ID of the gases to include
        isoID = np.zeros(ngas,dtype='int32')   #ID of the isotopes to include (0 is no differentiation between isotopes)
        gasnames = ['$^{12}$CO$_2$','Ar','N$_2$','O$_2$','$^{12}$CO','O','H$_2$','H','OH','HO$_2$','H$_2$O','H$_2$O$_2$','O(1D)','O$_3$']
        gasID[:] = [2,76,22,7,5,45,39,48,13,44,1,25,133,3]
        isoID[:] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    vmr = np.zeros((nh,ngas))
    zh = np.zeros(nh)
    temp = np.zeros(nh)
    press = np.zeros(nh)
    numdens = np.zeros(nh)


    #Reading the profiles (Mayor species)
    for i in range(nh):

        s = f.readline().split()
        zh[i] = float(s[0])*1.0e3       #Altitude (m)
        temp[i] = float(s[1])           #Temperature (K)
        press[i] = float(s[2])*100.     #Pressure (Pa)
        numdens[i] = float(s[3])*1.0e6  #Number density (molecules m-3)
        vmr[i,0] = float(s[4])          #CO2
        vmr[i,1] = float(s[5])          #Ar
        vmr[i,2] = float(s[6])          #N2
        vmr[i,3] = float(s[7])          #O2
        vmr[i,4] = float(s[8])          #CO
        vmr[i,5] = float(s[9])          #O
        vmr[i,6] = float(s[10])         #H2


    #Reading the second file (Minor species)
    f = open(filename2,'r')
    header = f.readline().split()
    for i in range(nh):

        #Reading the profiles
        s = f.readline().split()
        vmr[i,7] = float(s[1])  #H
        vmr[i,8] = float(s[2])  #OH
        vmr[i,9] = float(s[3])  #HO2
        vmr[i,10] = float(s[4])  #H2O
        vmr[i,11] = float(s[5])  #H2O2
        vmr[i,12] = float(s[6])  #O1D
        vmr[i,13] = float(s[7])  #O3

    ilast = 13

    #Reading the third file (Nitrogen species)
    if nitrogen_chemistry==True:
        f = open(filename3,'r')
        header = f.readline().split()
        for i in range(nh):

            #Reading the profiles
            s = f.readline().split()
            vmr[i,ilast+1] = float(s[1])  #N
            vmr[i,ilast+2] = float(s[2])  #NO
            vmr[i,ilast+3] = float(s[3])  #NO2
            vmr[i,ilast+4] = float(s[4])  #N2D

        ilast = ilast + 4

    #Making sure that the included VMRs add up to 1.0
    vmr = adjust_vmr(vmr)

    #Making plot
    if MakePlot==True:

        fig,(ax1,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(15,6))

        #Mayor species
        ix = 0
        for i in range(7):
            ax1.plot(vmr[:,ix],zh/1.0e3,label=gasnames[ix])
            ix = ix + 1

        #Minor species
        for i in range(7):
            ax2.plot(vmr[:,ix],zh/1.0e3,label=gasnames[ix])
            ix = ix + 1

        #Nitrogen species
        if nitrogen_chemistry==True:
            for i in range(4):
                ax3.plot(vmr[:,ix],zh/1.0e3,label=gasnames[ix])
                ix = ix + 1

        #(13C) species
        for i in range(2):
            ax4.plot(vmr[:,ix],zh/1.0e3,label=gasnames[ix])
            ix = ix + 1

        ax1.set_ylabel('Altitude (km)')
        ax2.set_ylabel('Altitude (km)')
        ax3.set_ylabel('Altitude (km)')
        ax4.set_ylabel('Altitude (km)')
        ax1.set_xlabel('Volume mixing ratio')
        ax2.set_xlabel('Volume mixing ratio')
        ax3.set_xlabel('Volume mixing ratio')
        ax4.set_xlabel('Volume mixing ratio')
        #ax1.set_xscale('log')
        ax2.set_xscale('log')
        ax3.set_xscale('log')
        ax4.set_xscale('log')
        ax1.legend()
        ax2.legend()
        ax3.legend()
        ax4.legend()
        ax1.grid()
        ax2.grid()
        ax3.grid()
        ax4.grid()
        plt.tight_layout()


        fig,(ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,6))

        ax1.plot(temp,zh/1.0e3)
        ax2.plot(press,zh/1.0e3)
        ax3.plot(numdens,zh/1.0e3)

        ax1.set_ylabel('Altitude (km)')
        ax2.set_ylabel('Altitude (km)')
        ax3.set_ylabel('Altitude (km)')
        ax1.set_xlabel('Temperature (K)')
        ax2.set_xlabel('Pressure (Pa)')
        ax3.set_xlabel('Number density (m$^{-3}$)')
        ax2.set_xscale('log')
        ax3.set_xscale('log')
        ax1.grid()
        ax2.grid()
        ax3.grid()
        plt.tight_layout()

        plt.show()

    return zh,temp,press,numdens,vmr,gasID,isoID

###############################################################################################

def ini_layer(h,temp,numdens,vmr,hlay):
    '''

    Function to initialise the atmospheric profiles

    Outputs
    --------

    zh(nh) :: Altitude of the boundaries of each layer (m)
    temp(nh) :: Temperature at the boundaries of each layer (K)
    press(nh) :: Pressure at the boundaries of each layer (Pa)
    numdens(nh) :: Number density at the boundaries of each layer (m-3)
    vmr(nh,ngas) :: Volume mixing ratio of each species at the boundaries of each layer

    '''

    from scipy.interpolate import interp1d

    nh = vmr.shape[0]
    ngas = vmr.shape[1]

    #Calculating the properties of each layer, assuming our profiles are the boundaries of the layer
    nlay = nh - 1
    Tlay = np.zeros(nlay)
    N0lay = np.zeros(nlay)
    VMRlay = np.zeros((nlay,ngas))

    #Interpolating temperatures
    f = interp1d(h,temp)
    Tlay = f(hlay)

    #Interpolating densities
    f = interp1d(h,np.log(numdens))
    N0lay = np.exp(f(hlay))

    #Interpolating VMRs
    f = interp1d(h,vmr,axis=0)
    VMRlay = f(hlay)

    Play = N0lay * phys_const['k_B'] * Tlay

    return hlay,Tlay,Play,N0lay,VMRlay


###############################################################################################

def adjust_vmr(vmr):

    '''
    Routine to adjust the VMR values so that they add up to 1.0 at all altitudes

    Inputs
    ------

    vmr(nh,ngas) :: volume mixing ratio of each species at each altitude

    '''

    nh = vmr.shape[0]
    ngas = vmr.shape[1]

    #Making sure that the sum of the volume mixing ratios add up to 1.0
    #dominent species = 1.0 - sum(all other species)
    for i in range(nh):

        vmrs = vmr[i,:]
        jmax = np.argmax(vmrs)

        vmr_rest = 0.0
        vmr_dom = 0.0
        for j in range(ngas):
            if j!=jmax:
                vmr_rest = vmr_rest + vmr[i,j]
            else:
                vmr_dom = vmr[i,j]

        vmr[i,jmax] = 1.0 - vmr_rest

    return vmr

###############################################################################################

def file_lines(fname):

    """
    FUNCTION NAME : file_lines()

    DESCRIPTION : Returns the number of lines in a given file

    INPUTS : 
 
        fname :: Name of the file

    OPTIONAL INPUTS: none
            
    OUTPUTS : 
 
        nlines :: Number of lines in file

    CALLING SEQUENCE:

        nlines = file_lines(fname)

    MODIFICATION HISTORY : Juan Alday (29/04/2019)

    """

    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


    