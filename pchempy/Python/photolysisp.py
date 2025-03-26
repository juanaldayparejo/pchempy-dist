import h5py
import numpy as np
import sys,os
from numba import jit
from pchempy import *

###############################################################################################################################

#Set of routines to work with the photolysis scheme of the photochemical code
def write_xs_hdf5(filename,wave,temp,xs,nreactions,sID,sISO,npr,pID,pISO,pf,branching_ratios):
    """
        FUNCTION NAME : write_xs_hdf5()
        
        DESCRIPTION : Write the cross sections and branching ratios into HDF5 file
        
        INPUTS :
        
            filename :: Name of the output .h5 file
            wave(nwave) :: Wavelength (nm)
            temp(ntemp) :: Temperature at which the cross sections are tabulated (K)
            xs(nwave,ntemp) :: Cross sections (cm2)
            nreactions :: Number of reactions associated with this photolysis
            sID :: ID of the gas that is photolysed
            sISO :: ID of the isotope that is photolysed
            npr(nreactions) :: Number of products in each reaction
            pID(npr,nreactions) :: ID of the photolysis products in each reaction
            pISO(npr,nreactions) :: Isotope ID of the photolysis products in each reaction
            pf(npr,nreactions) :: Number of molecules of a given product produced in each reaction
            branching_ratios(nwave,nreactions) :: Branching ratios for each of the reactions
            
        OPTIONAL INPUTS: none
        
        OUTPUTS :

            Output HDF5 file
        
        CALLING SEQUENCE:
        
            write_xs_hdf5(filename,wave,temp,xs,nreactions,sID,sISO,npr,pID,pISO,pf,branching_ratios)
        
        MODIFICATION HISTORY : Juan Alday (15/12/2023)
        
    """
    
    import h5py

    f = h5py.File(filename+'.h5','a')
    
    dset = f.create_dataset('WAVELENGTH',data=wave)
    dset.attrs['title'] = "Wavelength (NWAVE)"
    dset.attrs['units'] = 'nm'

    dset = f.create_dataset('TEMPERATURE',data=temp)
    dset.attrs['title'] = "Temperature at which the cross sections are tabulated (NTEMP)"
    dset.attrs['units'] = 'K'

    dset = f.create_dataset('CROSS_SECTIONS',data=xs)
    dset.attrs['title'] = "Cross sections at the different temperatures (NWAVE,NTEMP)"
    dset.attrs['units'] = 'cm2'

    dset = f.create_dataset('NREACTIONS',data=nreactions)
    dset.attrs['title'] = "Number of reactions associated with this photolysis"

    dset = f.create_dataset('sID',data=sID)
    dset.attrs['title'] = "ID of the gas that is photolysed"
    
    dset = f.create_dataset('sISO',data=sISO)
    dset.attrs['title'] = "Isotope ID of the gas that is photolysed"
    
    dset = f.create_dataset('NPRODUCTS',data=npr)
    dset.attrs['title'] = "Number of products in each reaction (NREACTIONS)"
    
    dset = f.create_dataset('pID',data=pID)
    dset.attrs['title'] = "ID of the photolysis products in each reaction (NPRODUCTS,NREACTIONS)"
    
    dset = f.create_dataset('pISO',data=pISO)
    dset.attrs['title'] = "Isotope ID of the photolysis products in each reaction (NPRODUCTS,NREACTIONS)"
    
    dset = f.create_dataset('pf',data=pf)
    dset.attrs['title'] = "Number of molecules of a given product produced in each reaction (NPRODUCTS,NREACTIONS)"
    
    dset = f.create_dataset('BRANCHING_RATIOS',data=branching_ratios)
    dset.attrs['title'] = "Branching ratios for each of the reactions (NWAVE,NREACTIONS)"
    
    f.close()

###############################################################################################################################

def read_xs_hdf5(filename):
    """
        FUNCTION NAME : read_xs_hdf5()
        
        DESCRIPTION : Read the cross sections and branching ratios into HDF5 file
        
        INPUTS :
        
            filename :: Name of the input .h5 file
            
        OPTIONAL INPUTS: none
        
        OUTPUTS :

            wave(nwave) :: Wavelength (nm)
            temp(ntemp) :: Temperature at which the cross sections are tabulated (K)
            xs(nwave,ntemp) :: Cross sections (cm2)
            nreactions :: Number of reactions associated with this photolysis
            sID :: ID of the gas that is photolysed
            sISO :: ID of the isotope that is photolysed
            npr(nreactions) :: Number of products in each reaction
            pID(npr,nreactions) :: ID of the photolysis products in each reaction
            pISO(npr,nreactions) :: Isotope ID of the photolysis products in each reaction
            pf(npr,nreactions) :: Number of molecules of a given product produced in each reaction
            branching_ratios(nwave,nreactions) :: Branching ratios for each of the reactions
        
        CALLING SEQUENCE:
        
            wave,temp,xs,nreactions,sID,sISO,npr,pID,pISO,pf,branching_ratios = read_xs_hdf5(filename)
        
        MODIFICATION HISTORY : Juan Alday (15/12/2023)
        
    """
    
    f = h5py.File(filename,'r')

    wave = np.array(f.get('WAVELENGTH'))
    temp = np.array(f.get('TEMPERATURE'))
    xs = np.array(f.get('CROSS_SECTIONS'))
    nreactions = np.int32(f.get('NREACTIONS'))
    sID = np.array(f.get('sID'),dtype='int32')
    sISO = np.array(f.get('sISO'),dtype='int32')
    npr = np.array(f.get('NPRODUCTS'),dtype='int32')
    pID = np.array(f.get('pID'),dtype='int32')
    pISO = np.array(f.get('pISO'),dtype='int32')
    pf = np.array(f.get('pf'))
    branching_ratios= np.array(f.get('BRANCHING_RATIOS'))
    
    f.close()
    
    return wave,temp,xs,nreactions,sID,sISO,npr,pID,pISO,pf,branching_ratios

###############################################################################################################################

def read_solflux_hdf5(filename):
    """
        FUNCTION NAME : read_solflux_hdf5()
        
        DESCRIPTION : Function to read the solar flux (at 1AU) from an HDF5 file
        
        INPUTS :

            filename :: Name of the HDF5 file
            
        OUTPUTS :
        
            wave(nwave) :: Wavelength (nm)
            solflux(nwave) :: Solar flux at 1 AU (W cm-2 nm-1)
        
        CALLING SEQUENCE:
        
            wave,solflux = read_solflux_hdf5(filename)
        
        MODIFICATION HISTORY : Juan Alday (13/10/2024)
        
    """
    
    f = h5py.File(filename,'r')
    wave = np.array(f.get('WAVELENGTH'))
    solflux = np.array(f.get('SOLFLUX'))
    f.close()
    
    return wave,solflux

###############################################################################################################################

def write_solflux_hdf5(wave,solflux,filename):
    """
        FUNCTION NAME : write_solflux_hdf5()
        
        DESCRIPTION : Function to write the solar flux (at 1AU) to an HDF5 file
        
        INPUTS :

            wave(nwave) :: Wavelength (nm)
            solflux(nwave) :: Solar flux at 1 AU (W cm-2 nm-1)
            
        OUTPUTS : HDF5 file
        
        CALLING SEQUENCE:
        
            write_solflux_hdf5(wave,solflux,filename)
        
        MODIFICATION HISTORY : Juan Alday (13/10/2024)
        
    """

    import h5py

    f = h5py.File(filename+'.h5','a')
    
    dset = f.create_dataset('WAVELENGTH',data=wave)
    dset.attrs['title'] = "Wavelength (NWAVE)"
    dset.attrs['units'] = 'nm'

    dset = f.create_dataset('SOLFLUX',data=solflux)
    dset.attrs['title'] = "Solar flux at 1 AU"
    dset.attrs['units'] = 'W cm-2 nm-1'

    f.close()

###############################################################################################################################

@jit(nopython=True)
def interp_xs_temp(temp,xs,tlay):
    """
        FUNCTION NAME : interp_xs_temp()
        
        DESCRIPTION : Interpolate the cross sections to a set of temperatures
        
        INPUTS :
        
            temp(ntemp) :: Temperatures at which the cross sections are stored (K)
            xs(nwave,ntemp) :: Cross sections at the different temperatures (cm2)
            templay(ntemplay) :: Temperatures of the atmospheric layers (K)
            
        OUTPUTS :

            xslay(nwave,ntemp) :: Interpolated cross sections at the temperature of the atmospheric layers (K)
        
        CALLING SEQUENCE:
        
            xslay = interp_xs_temp(temp,xs,tlay)
        
        MODIFICATION HISTORY : Juan Alday (13/10/2024)
        
    """
    
    ntemp = len(temp)
    nlay = len(tlay)
    nwave = xs.shape[0]
    
    xslay = np.zeros((nwave,nlay))
    
    #If there is only one temperature we cannot interpolate
    if ntemp == 1:
        for ilay in range(nlay):
            xslay[:,ilay] = xs[:,0]
    else:
        
        for ilay in range(nlay):
            
            t1 = tlay[ilay]
            
            #Finding the indices of the nearest points in the temperature array
            if t1 < np.min(temp):
                t1 = np.min(temp)
            if t1 > np.max(temp):
                t1 = np.max(temp)

            it = np.searchsorted(temp, t1) - 1
            if it < 0:
                it = 0
            if it >= ntemp - 1:
                it = len(temp) - 2

            #Calculating the interpolation factor
            u = (t1 - temp[it]) / (temp[it + 1] - temp[it])   
            klo = xs[:,it]
            khi = xs[:,it+1]
            
            #Interpolating cross sections
            igood = np.where((klo>0.0) & (khi>0.0))[0]
            xslay[igood,ilay] = np.exp((1.-u) * np.log(klo[igood]) + u * np.log(khi[igood]))
    
            ibad = np.where((klo <= 0.0) & (khi <= 0.0) )[0]
            xslay[ibad,ilay] = (1.-u) * klo[ibad] + u * khi[ibad]

    return xslay

@jit(nopython=True)
def wlgrid_lmd():
    """
        FUNCTION NAME : wlgrid_lmd()
        
        DESCRIPTION : Function to define the wavelength grid used by the LMD-GCM (good for testing purposes)
        
                       high-resolution mode (3789 intervals)

                       0-108 nm :  1.0  nm
                       108-124 nm :  0.1  nm
                       124-175 nm :  0.5  nm
                       175-205 nm :  0.01 nm
                       205-365 nm :  0.5  nm
                       365-850 nm :  5.0  nm
        
        INPUTS : None
            
        OUTPUTS :

            nbin :: Number of bins
            wl(nbin) :: Lower wavelength of each bin (nm)
            wu(nbin) :: Upper wavelength of each bin (nm)
            wc(nbin) :: Central wavelength of each bin (nm)
        
        CALLING SEQUENCE:
        
            wl,wu,wc = wlgrid_lmd()
        
        MODIFICATION HISTORY : Juan Alday (13/10/2024)
        
    """
    
    nbin = 3788
    
    wl = np.zeros(nbin)
    wu = np.zeros(nbin)
    wc = np.zeros(nbin)
    
    #1nm intervals from 0 to 108 nm
    kw = 0
    wincr = 1.
    for iw in range(108):
        wl[kw] = iw
        wu[kw] = wl[kw] + wincr
        wc[kw] = (wl[kw] + wu[kw])/2.
        kw += 1
        
    #0.1nm intervals from 108 to 124
    wincr = 0.1
    for iw in range(1080,1240):
        wl[kw] = iw/10.
        wu[kw] = wl[kw] + wincr
        wc[kw] = (wl[kw] + wu[kw])/2.
        kw += 1
        
    #0.5nm intervals from 124 to 175 nm
    wincr = 0.5
    for iw in range(1240, 1750, 5):
        wl[kw] = iw/10.
        wu[kw] = wl[kw] + wincr
        wc[kw] = (wl[kw] + wu[kw])/2.
        kw += 1

    #0.01nm intervals from 175 to 205 nm
    wincr = 0.01
    for iw in range(17500,20500):
        wl[kw] = iw/100.
        wu[kw] = wl[kw] + wincr
        wc[kw] = (wl[kw] + wu[kw])/2.
        kw += 1
        
    #0.5nm intervals from 205 to 365 nm
    wincr = 0.5
    for iw in range(2050,3650,5):
        wl[kw] = iw/10.
        wu[kw] = wl[kw] + wincr
        wc[kw] = (wl[kw] + wu[kw])/2.
        kw += 1

    #5nm intervals from 365 to 855 nm
    wincr = 5.0
    for iw in range(365,855,5):
        wl[kw] = iw
        wu[kw] = wl[kw] + wincr
        wc[kw] = (wl[kw] + wu[kw])/2.
        kw += 1
        
    return wl,wu,wc

###############################################################################################################################

@jit(nopython=True)
def bin_data(ng, xg, n, x, y):
    """
        FUNCTION NAME : bin_data()
        
        DESCRIPTION : Map input data given on discrete points onto a set of target bins via linear interpolation.
                      The average value in each target bin is found by averaging the trapezoidal area under
                      the input data curve, constructed by linearly connecting the discrete input values.
        
        INPUTS :
        
            ng :: Number of bins in the target grid
            xg(ng) :: Target grid (e.g., wavelength grid); bin i is defined as [xg[i], xg[i+1]]
            n :: Number of points in the input grid
            x(n) :: Grid on which input data is defined
            y(n) :: Input y-data corresponding to the x-data
            
        OUTPUTS :

            yg(ng) :: Interpolated y-values for the target bins
        
        RAISES:
        
            ValueError :: If the input x or xg is not sorted or if the xg values are outside the range of x.
        
        CALLING SEQUENCE:
        
            yg = bin_data(ng, xg, n, x, y)
        
        MODIFICATION HISTORY : Juan Alday (13/10/2024)
        
    """

    # Initialize the output yg array
    yg = np.zeros(ng)

    # Check if the input x data is sorted
    for i in range(1, n):
        if x[i] <= x[i-1]:
            raise ValueError("Error: input data x is not sorted")

    # Check if the target grid xg is sorted
    for i in range(1, ng):
        if xg[i] <= xg[i-1]:
            raise ValueError("Error: target grid xg is not sorted")

    # Check if the xg values are within the range of x
    if x[0] > xg[0] or x[-1] < xg[-1]:
        raise ValueError("Error: Data do not span the grid. Use ADDPNT to expand data.")

    # Initialize variables
    jstart = 0
    ngintv = ng - 1

    # Loop over target intervals
    for i in range(ngintv):
        area = 0.0
        xgl = xg[i]
        xgu = xg[i+1]

        k = jstart

        # Discard points before the first grid interval and after the last grid interval
        while k <= n - 2:
            
            if x[k+1] <= xgl:
                jstart = k + 1
                k += 1
                continue

            if x[k] >= xgu:
                break

            # Compute the x-coordinates for interpolation
            a1 = max(x[k], xgl)
            a2 = min(x[k+1], xgu)

            if x[k+1] == x[k]:
                darea = 0.0
            else:
                slope = (y[k+1] - y[k]) / (x[k+1] - x[k])
                b1 = y[k] + slope * (a1 - x[k])
                b2 = y[k] + slope * (a2 - x[k])
                darea = (a2 - a1) * (b2 + b1) / 2.0

            area += darea
            k += 1

        # Calculate the average y value for the interval
        yg[i] = area / (xgu - xgl)

    return yg

###############################################################################################################################

@jit(nopython=True)
def sphers(nlev, z, zen, radius=3393.0):
    """
    FUNCTION NAME : sphers()

    DESCRIPTION : Calculate the slant path over vertical depth (ds/dh) in spherical geometry.
                  Based on the method by A. Dahlback and K. Stamnes (1991) to compute the radiation 
                  field for photolysis and heating at twilight.

    INPUTS :

        nlev :: Number of specified altitude levels in the grid
        z(nlev) :: Specified altitude working grid (km)
        zen :: Solar zenith angle (degrees)

    OPTIONAL INPUTS :

        radius :: Radius of the planet (km)

    OUTPUTS :

        dsdh(nlev + 1, nlev) :: Slant path of the direct beam through each layer
        nid(nlev + 1) :: Number of layers crossed by the direct beam

    CALLING SEQUENCE:

        dsdh, nid = sphers(nlev, z, zen)

    """
    # Convert SZA to radians
    zenrad = zen / 180.0 * np.pi

    # Number of layers
    nlay = nlev - 1

    # Altitude adjusted to the elevation above Mars' surface
    re = radius + z[0]

    # Correspondingly, z is adjusted to the elevation above Mars' surface
    ze = z - z[0]

    # Inverse coordinates of z
    zd = np.zeros(nlev)
    zd[0] = ze[nlev - 1]
    for k in range(1, nlay + 1):
        zd[k] = ze[nlev - 1 - k]

    # Initialize outputs: slant path and number of crossed layers
    dsdh = np.zeros((nlev + 1, nlev), dtype=np.float64)
    nid = np.zeros(nlev + 1, dtype=np.int32)

    # Calculate ds/dh for each layer
    for i in range(nlay + 1):
        rpsinz = (re + zd[i]) * np.sin(zenrad)

        if zen > 90.0 and rpsinz < re:
            nid[i] = -1
        else:
            id = i
            if zen > 90.0:
                for j in range(1, nlay + 1):
                    if (rpsinz < (zd[j - 1] + re)) and (rpsinz >= (zd[j] + re)):
                        id = j

            for j in range(1, id + 1):
                sm = 1.0
                if j == id and id == i and zen > 90.0:
                    sm = -1.0

                rj = re + zd[j - 1]
                rjp1 = re + zd[j]
                dhj = zd[j - 1] - zd[j]

                ga = rj**2 - rpsinz**2
                gb = rjp1**2 - rpsinz**2

                ga = max(ga, 0.0)
                gb = max(gb, 0.0)

                if id > i and j == id:
                    dsj = np.sqrt(ga)
                else:
                    dsj = np.sqrt(ga) - sm * np.sqrt(gb)

                dsdh[i, j - 1] = dsj / dhj

            nid[i] = id

    return dsdh, nid

###############################################################################################################################

@jit(nopython=True)
def setdust_mars(nlay,alt,tau,nwave):
    """
        FUNCTION NAME : setdust_mars()
        
        DESCRIPTION : Set dust properties for each specified altitude layer
        
        INPUTS :
        
            nlay :: Number of layers
            alt(nlay) :: Altitude of each layer (km)
            tau :: Dust column optical depth
            nwave :: Number of wavelengths
            
        OPTINAL INPUTS : None
            
        OUTPUTS :

            dtau(nlay,nwave) :: Dust optical depth in each altitude layer
            omega(nlay,nwave) :: Single scattering albedo of each layer
            g(nlay,nwave) :: Asymmetry factor in each layer
            
        CALLING SEQUENCE:
        
            dtau,omega,g = setdust_mars(nlay,alt,tau,nwave)
        
        MODIFICATION HISTORY : Juan Alday (15/10/2024)
        
    """
    
    omegax  = 0.622 #single scattering albedo : wolff et al.(2010) at 258 nm
    gx     = 0.88   #asymmetry factor : mateshvili et al. (2007) at 210 nm
    scaleh = 10.    #scale height (km)
    gamma  = 0.03   #conrath parameter

    dtau = np.zeros((nlay,nwave))
    omega = np.zeros((nlay,nwave))
    g = np.zeros((nlay,nwave))
    
    #Optical depth profile
    tautot = 0.
    for ilay in range(nlay-1):
        dz = alt[ilay+1] - alt[ilay]
        tautot += np.exp(gamma*(1.0-np.exp(alt[ilay]/scaleh)))*dz

    q0 = tau/tautot
    for ilay in range(nlay-1):
        dz = alt[ilay+1] - alt[ilay]
        dtau[ilay,:] = q0 * np.exp(gamma*(1. - np.exp(alt[ilay]/scaleh)))*dz
        omega[ilay,:] = omegax
        g[ilay,:] = gx 

    return dtau,omega,g

###############################################################################################################################

@jit(nopython=True)
def setray_mars(coldens,wave):
    """
        FUNCTION NAME : setray_mars()
        
        DESCRIPTION : Calculate the Rayleigh optical depth in each layer for a CO2 dominated atmosphere
                      CO2 rayleigh cross-section from ityaksov et al., chem. phys. lett., 462, 31-34, 2008
        
        INPUTS :
        
            coldens(nlay) :: Atmospheric column density in each layer (m-2)
            wave(nwave) :: Wavelength (nm)
            
        OPTINAL INPUTS : None
            
        OUTPUTS :

            dtaur(nlay,nwave) :: Rayleigh optical depth in each layer

        CALLING SEQUENCE:
        
            dtaur = setray_mars(coldens,wave)
        
        MODIFICATION HISTORY : Juan Alday (15/10/2024)
        
    """
    
    nlay = len(coldens)
    nwave = len(wave)
    
    dtaur = np.zeros((nlay,nwave))

    #CO2 Rayleigh cross section
    nu = 1./(wave*1.e-7)  #cm-1
    srayl = 1.78e-26*nu**(4. + 0.625)
    srayl *= 1.0e-20  #cm2

    for i in range(nlay):
        dtaur[i,:] = coldens[i] * 1.0e-4 * srayl[:]

    return dtaur


###############################################################################################################################

@jit(nopython=True)
def ps2str(nlev, zen, rsfc, tauu, omu, gu, dsdh, nid, delta):
    """
    Translated and optimized version of the Fortran subroutine ps2str.
    """
    # Constants
    largest = 1.e+36
    precis = 1.e-7
    eps = 1.e-3  # From PARAMETER (eps = 1.E-3)
    pifs = 1.0
    fdn0 = 0.0
    nlayer = nlev - 1
    pi = np.arccos(-1.0)
    dr = pi / 180.0
    mu = np.cos(zen * dr)
    
    # Initialize outputs
    fup = np.zeros(nlev)
    fdn = np.zeros(nlev)
    fdr = np.zeros(nlev)
    eup = np.zeros(nlev)
    edn = np.zeros(nlev)
    edr = np.zeros(nlev)
    
    # Local variables
    tauc = np.zeros(nlev + 1)
    tausla = np.zeros(nlev + 1)
    mu2 = np.full(nlev + 1, 1.0 / np.sqrt(largest))
    
    # Internal coefficients and matrix
    lam = np.zeros(nlev)
    taun = np.zeros(nlev)
    bgam = np.zeros(nlev)
    e1 = np.zeros(nlev)
    e2 = np.zeros(nlev)
    e3 = np.zeros(nlev)
    e4 = np.zeros(nlev)
    cup = np.zeros(nlev)
    cdn = np.zeros(nlev)
    cuptn = np.zeros(nlev)
    cdntn = np.zeros(nlev)
    mu1 = np.zeros(nlev)
        
    # Prepare input arrays
    gi = np.zeros(nlev)
    omi = np.zeros(nlev)
    
    # Apply delta scaling if needed
    if not delta:
        for i in range(nlayer):
            gi[i] = gu[i]
            omi[i] = omu[i]
            taun[i] = tauu[i]
    else:
        for i in range(nlayer):
            f = gu[i] * gu[i]
            gi[i] = (gu[i] - f) / (1. - f)
            omi[i] = (1.0 - f) * omu[i] / (1.0 - omu[i] * f)
            taun[i] = (1.0 - omu[i] * f) * tauu[i]
    
    # Calculate slant optical depth at the top of the atmosphere when zen > 90
    if zen > 90.0:
        if nid[0] < 0:
            tausla[0] = largest
        else:
            sum_val = 0.0
            for j in range(nid[0]):
                sum_val += 2.0 * taun[j] * dsdh[0, j]
            tausla[0] = sum_val
    
    # Main loop over layers
    for i in range(nlayer):
        
        g = gi[i]
        om = omi[i]
        tauc[i+1] = tauc[i] + taun[i]
        
        # Stay away from 1 by precision
        tempg = min(abs(g), 1.0 - precis)
        g = np.copysign(tempg, g)
        om = min(om, 1.0 - precis)
        
        # Calculate slant optical depth
        if nid[i+1] < 0:
            tausla[i+1] = largest
        else:
            sum_val = 0.0
            min_nid_i = min(nid[i+1], i+1)
            for j in range(min_nid_i):
                sum_val += taun[j] * dsdh[i+1, j]
            for j in range(min_nid_i, nid[i+1]):
                sum_val += 2.0 * taun[j] * dsdh[i+1, j]
            tausla[i+1] = sum_val
            if tausla[i+1] == tausla[i]:
                mu2[i+1] = np.sqrt(largest)
            else:
                denom = tausla[i+1] - tausla[i]
                if denom == 0.0:
                    denom = 1e-12  # Avoid division by zero
                mu2[i+1] = (tauc[i+1] - tauc[i]) / denom
                mu2[i+1] = np.copysign(max(abs(mu2[i+1]), 1.0 / np.sqrt(largest)), mu2[i+1])
        
        # Hemispheric mean approximation
        gam1 = 2.0 - om * (1.0 + g)
        gam2 = om * (1.0 - g)
        gam3 = (2.0 - g * mu) / 4.0
        gam4 = 1.0 - gam3
        mu1[i] = 0.5
        
        # lambda and bgam
        lam[i] = np.sqrt(gam1 * gam1 - gam2 * gam2)
        if gam2 != 0.0:
            bgam[i] = (gam1 - lam[i]) / gam2
        else:
            bgam[i] = 0.0
        
        expon = np.exp(-lam[i] * taun[i])
        
        e1[i] = 1.0 + bgam[i] * expon
        e2[i] = 1.0 - bgam[i] * expon
        e3[i] = bgam[i] + expon
        e4[i] = bgam[i] - expon
        
        # Prevent division by zero
        mu2i_sq = mu2[i+1] * mu2[i+1]
        if mu2i_sq == 0.0:
            mu2i_sq = 1e-12  # Avoid division by zero
        divisr = lam[i] * lam[i] - 1.0 / mu2i_sq
        temp = max(eps, abs(divisr))
        divisr = np.copysign(temp, divisr)
        
        up = om * pifs * ((gam1 - 1.0 / mu2[i+1]) * gam3 + gam4 * gam2) / divisr
        dn = om * pifs * ((gam1 + 1.0 / mu2[i+1]) * gam4 + gam2 * gam3) / divisr
        
        expon0 = np.exp(-tausla[i])
        expon1 = np.exp(-tausla[i+1])
        
        cup[i] = up * expon0
        cdn[i] = dn * expon0
        cuptn[i] = up * expon1
        cdntn[i] = dn * expon1
        
    # Set up matrix
    ssfc = rsfc * mu * np.exp(-tausla[nlayer]) * pifs
    mrows = 2 * nlayer
    
    # Initialize matrix arrays (indices from 1 to mrows)
    a = np.zeros(mrows)
    b = np.zeros(mrows)
    d = np.zeros(mrows)
    e = np.zeros(mrows)
    y = np.zeros(mrows)
    
    # Set up first row of matrix
    i = 0
    a[0] = 0.0
    b[0] = e1[i]
    d[0] = -e2[i]
    e[0] = fdn0 - cdn[i]
    
    # Set up odd rows 3 thru (MROWS - 1)
    i = 0
    for row in range(2, mrows-1, 2):
        i += 1
        a[row] = e2[i-1] * e3[i-1] - e4[i-1] * e1[i-1]
        b[row] = e1[i-1] * e1[i] - e3[i-1] * e3[i]
        d[row] = e3[i-1] * e4[i] - e1[i-1] * e2[i]
        e[row] = e3[i-1] * (cup[i] - cuptn[i-1]) + e1[i-1] * (cdntn[i-1] - cdn[i])
    
    # Set up even rows 2 thru (MROWS - 2)
    i = 0
    for row in range(1, mrows-2, 2):
        i += 1
        a[row] = e2[i] * e1[i-1] - e3[i-1] * e4[i]
        b[row] = e2[i-1] * e2[i] - e4[i-1] * e4[i]
        d[row] = e1[i] * e4[i] - e2[i] * e3[i]
        e[row] = (cup[i] - cuptn[i-1]) * e2[i] - (cdn[i] - cdntn[i-1]) * e4[i]
    
    # Set up last row of matrix at MROWS
    row = mrows-1
    i = nlayer-1
    a[row] = e1[i] - rsfc * e3[i]
    b[row] = e2[i] - rsfc * e4[i]
    d[row] = 0.0
    e[row] = ssfc - cuptn[i] + rsfc * cdntn[i]
    
    # Solve tri-diagonal matrix
    y = tridiag(a, b, d, e, mrows)
    
    # Unfold solution of matrix, compute output fluxes
    row = 0
    lev = 0  # In Fortran lev starts from 1, in Python we start from 0
    j = 0
    
    fdr[lev] = np.exp(-tausla[0])
    edr[lev] = mu * fdr[lev]
    edn[lev] = fdn0
    eup[lev] = y[row] * e3[j] - y[row + 1] * e4[j] + cup[j]
    fdn[lev] = edn[lev] / mu1[j]
    fup[lev] = eup[lev] / mu1[j]
    
    for lev in range(1, nlayer + 1):
        fdr[lev] = np.exp(-tausla[lev])
        edr[lev] = mu * fdr[lev]
        edn[lev] = y[row] * e3[j] + y[row + 1] * e4[j] + cdntn[j]
        eup[lev] = y[row] * e1[j] + y[row + 1] * e2[j] + cuptn[j]
        fdn[lev] = edn[lev] / mu1[j]
        fup[lev] = eup[lev] / mu1[j]
        
        row += 2
        j += 1
    
    # Return outputs
    return fdr, fup, fdn, edr, eup, edn

###############################################################################################################################

@jit(nopython=True)
def tridiag(a, b, c, r, n):
    """
    Solves a tridiagonal system of equations using the Thomas algorithm.
    Optimized with Numba.
    """
    
    u = np.zeros(n)
    
    if b[0] == 0.0: 
        raise ValueError('error in tridiag (bet == 0)')
    bet = b[0]
    u[0] = r[0]/bet
    
    gam = np.zeros(n)
    for j in range(1,n):
        gam[j] = c[j-1]/bet
        bet = b[j] - a[j]*gam[j]
        if bet == 0.0:
            raise ValueError('error in tridiag (bet == 0) - second case')
        u[j] = (r[j] - a[j]*u[j-1])/bet

    # Back substitution
    for i in range(n-1):
        j = n-2-i
        u[j] = u[j] - gam[j+1]*u[j+1]
    
    return u

###############################################################################################################################

@jit(nopython=True)
def rtlink(nlev, nw, iw, ag, zen, dsdh, nid, dtrl, dagas, dtcld, omcld, gcld,
           dtaer, omaer, gaer):
    """
    Translated version of the Fortran subroutine rtlink.

    Parameters:
    nlev : int
        Number of levels.
    nw : int
        Number of wavelengths.
    iw : int
        Index for wavelength.
    ag : float
        Surface albedo.
    zen : float
        Solar zenith angle in degrees.
    dsdh : 2D array of shape (nlev + 1, nlev)
        Slant path of direct beam through each layer (indices from 0 to nlev).
    nid : array of length nlev + 1
        Number of layers crossed by the direct beam (indices from 0 to nlev).
    dtrl : 2D array of shape (nlev, nw)
        Rayleigh scattering optical depth per layer.
    dagas : 2D array of shape (nlev, nw)
        Absorption optical depth per layer (gas).
    dtcld, omcld, gcld : 2D arrays of shape (nlev, nw)
        Optical depth, single scattering albedo, asymmetry factor for clouds.
    dtaer, omaer, gaer : 2D arrays of shape (nlev, nw)
        Optical depth, single scattering albedo, asymmetry factor for aerosols.

    Returns:
    edir, edn, eup : arrays of length nlev
        Direct, downwelling, and upwelling spectral irradiance.
    fdir, fdn, fup : arrays of length nlev
        Direct, downwelling, and upwelling actinic flux.
    """
    # Initialize outputs
    edir = np.zeros(nlev)
    edn = np.zeros(nlev)
    eup = np.zeros(nlev)
    fdir = np.zeros(nlev)
    fdn = np.zeros(nlev)
    fup = np.zeros(nlev)

    largest = 1.e+36
    delta = True  # Equivalent to Fortran's "logical, save :: delta = .true."

    # Initialize local variables
    dt = np.zeros(nlev)
    om = np.zeros(nlev)
    g = np.zeros(nlev)

    # Loop over i = 1 to nlev - 1 (Fortran indices)
    for i in range(nlev-1):  # Python indices from 1 to nlev - 1

        # Calculate scattering and absorption components
        dscld = dtcld[i, iw] * omcld[i, iw]
        dacld = dtcld[i, iw] * (1.0 - omcld[i, iw])

        dsaer = dtaer[i, iw] * omaer[i, iw]
        daaer = dtaer[i, iw] * (1.0 - omaer[i, iw])

        dtsct = dtrl[i, iw] + dscld + dsaer
        dtabs = dagas[i, iw] + dacld + daaer

        dtabs = max(dtabs, 1.0 / largest)
        dtsct = max(dtsct, 1.0 / largest)

        # Invert z-coordinate
        ii = nlev - 2 - i 

        dt[ii] = dtsct + dtabs
        if (dtsct + dtabs) != 0.0:
            om[ii] = dtsct / (dtsct + dtabs)
        else:
            om[ii] = 0.0  # Avoid division by zero

        if dtsct == 1.0 / largest:
            om[ii] = 1.0 / largest

        if dtsct != 0.0:
            gcld_term = gcld[i, iw] * dscld
            gaer_term = gaer[i, iw] * dsaer
            g[ii] = (gcld_term + gaer_term) / dtsct
        else:
            g[ii] = 0.0  # Avoid division by zero

    # Call ps2str
    fdr, fup_temp, fdn_temp, edr, eup_temp, edn_temp = ps2str(
        nlev, zen, ag, dt, om, g, dsdh, nid, delta
    )

    # Output (invert z-coordinate)
    for i in range(nlev):
        ii = nlev -1 - i  

        fdir[i] = fdr[ii]
        fup[i] = fup_temp[ii]
        fdn[i] = fdn_temp[ii]
        edir[i] = edr[ii]
        eup[i] = eup_temp[ii]
        edn[i] = edn_temp[ii]

    return edir, edn, eup, fdir, fdn, fup

###############################################################################################################################

@jit(nopython=True)
def photolysis_rates(hlay,gasID,isoID,Nlay,wl,wu,wc,gasactID,isoactID,xs,solflux,planet='Mars',zen=45.,tau_aero=1.,radius=3393.,galb=0.3,dist_sun=1.5,plots=False):
    
    """
        FUNCTION NAME : photolysis_rates()
        
        DESCRIPTION : Calculate the photolysis rates in each layer of the atmosphere and for each of the involved reactions
        
        INPUTS :
        
            h(nlay) :: Altitude of each layer (m)
            gasID(ngas) :: Gas ID of each gas in the atmosphere
            isoID(ngas) :: Isotope ID of each gas in the atmosphere
            N(nlay,ngas) :: Number density of each gas in each layer (m-3)
            wl,wu,wc(nwave) :: Lower, upper and central wavelengths of each wavelength bin (nm)
            gasactID(ngasact) :: Gas ID of the gas that is photolised in each reaction
            isoactID(ngasact) :: Isotope ID of the gas that is photolised in each reaction
            xs(nwave,nlay,nreactions) :: Photolysis cross sections in each layer for each reaction (cm2)
            solflux(nwave) :: Solar flux at 1AU (W cm-2 nm-1)

        OPTIONAL INPUTS:
         
            zen :: Solar zenith angle
            tau_aero :: Dust optical depth
            radius :: Radius of the planet (km)
            galb :: Ground albedo
            dist_sun :: Distance to parent star (AU)
        
        OUTPUTS :

            rrates(nlay,nreactions) :: Photolysis rates for each reaction at each layer (s-1)
        
        CALLING SEQUENCE:
        
            rrates = photolysis_rates(hlay,gasID,isoID,Nlay,wl,wu,wc,gasactID,isoactID,xs,solflux)
        
        MODIFICATION HISTORY : Juan Alday (15/12/2023)
        
    """
    
    nwave = len(wl)
    nreactions = xs.shape[2]
    nlay = len(hlay)
    
    N0lay = np.sum(Nlay,axis=1)   #m-3
    delh = hlay[1] - hlay[0]      #m (assumed to be constant)
    coldens = N0lay * delh        #m-2
    coldens_i = Nlay * delh       #m-2
    
    #Calculate the gas optical depths
    dtau_gas_i = np.zeros((nlay,nwave,nreactions))
    for ir in range(nreactions):
        
        #Locating the gas in the atmosphere
        igas = np.where( (gasID==gasactID[ir]) & (isoID==isoactID[ir]) )[0]
        
        if len(igas)<=0:
            print('gasID',gasID)
            print('isoID',isoID)
            print('gasactID',gasactID)
            print('isoactID',isoactID)
            raise ValueError('ID of active gas not found in the atmosphere')
        
        #Calculating the optical depth for this particular reaction
        igas = igas[0]
        for ilay in range(nlay):
            dtau_gas_i[ilay,:,ir] = coldens_i[ilay,igas] * 1.0e-4 * xs[ilay,:,ir]
    
    #Calculating the total gas optical depth
    dtau_gas = np.sum(dtau_gas_i,axis=2)  #(nlay,nwave)
    
    
    if planet == 'Mars':
        
        #Calculate the aerosol optical depth
        dtau_aero,omega_aero,g_aero = setdust_mars(nlay,hlay/1.0e3,tau_aero,nwave)  #(nlay,nwave)
        
        #Calculating the Rayleigh optical depth
        dtau_ray = setray_mars(coldens,wc)   #(nlay,nwave)
        
        #Calculating the cloud optical depth
        dtau_cld = np.zeros(dtau_aero.shape)
        omega_cld = np.zeros(dtau_aero.shape)
        g_cld = np.zeros(dtau_aero.shape)
        
    else:
        
        dtau_aero = np.zeros(dtau_gas.shape)
        omega_aero = np.zeros(dtau_gas.shape)
        g_aero = np.zeros(dtau_gas.shape)
        
        dtau_ray = np.zeros(dtau_gas.shape)
        
        dtau_cld = np.zeros(dtau_gas.shape)
        omega_cld = np.zeros(dtau_gas.shape)
        g_cld = np.zeros(dtau_gas.shape)
    
    #Calculating the paths
    dsdh, nid = sphers(nlay, hlay/1.0e3, zen, radius=radius)

    #Convert solar flux to to photon.s-1.nm-1.cm-2
    sol = solflux * wc *5.039e11
    
    #Calculating the solar flux at planet
    factor = (1./dist_sun)**2.
    fplanet = sol * factor
    
    #Calculating the rates in each wavelength
    rrates = np.zeros((nlay,nreactions))
    for iw in range(nwave):
        
        #Calculating the radiative transfer
        edir, edn, eup, fdir, fdn, fup = rtlink(nlay, nwave, iw, galb, zen, dsdh, nid, dtau_ray, dtau_gas, dtau_cld, omega_cld, g_cld, dtau_aero, omega_aero, g_aero)
    
        #Calculating the spherical actinic flux
        saflux = np.zeros(nlay)
        for ilay in range(nlay):
            saflux[ilay] = fplanet[iw] * (fdir[ilay] + fdn[ilay] + fup[ilay])
        
        #Photolysis rate integration    
        for ir in range(nreactions):
            for ilay in range(nlay):
                deltaj = saflux[ilay]*xs[ilay,iw,ir]
                rrates[ilay,ir] += deltaj*(wu[iw]-wl[iw])
                
    return rrates

###############################################################################################################################

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
        
    
    


    