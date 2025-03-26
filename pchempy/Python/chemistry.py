import h5py
import numpy as np
from numba import njit

#############################################################################################################################

@jit(nopython=True)
def calc_jacobian_chemistry(nlay, ngas, ilay, c, nreactions, rtype, ns, sID_pos, sf, npr, pID_pos, pf, rrates):
    """
    Optimized routine to calculate the values of the chemical Jacobian matrix.

    Parameters:
    -----------
    nlay :: Number of atmospheric layers.
    ngas :: Number of gas species.
    ilay :: Level index at which to calculate the Jacobian matrix.
    c(nlay,ngas) :: Number density of each species.
    nreactions :: Number of reactions.
    rtype(nreactions) :: Reaction types.
    ns(nreactions) :: Number of source species.
    sID_pos(2,nreactions) :: Position indices of source species in the gasID array.
    sf(2,nreactions) :: Number of molecules for each source.
    npr(nreactions) :: Number of product species.
    pID_pos(2,nreactions) :: Position indices of product species in the gasID array.
    pf(2,nreactions) :: Number of molecules for each product.
    rrates(nlay,nreactions) :: Reaction rates (nlay, nreactions).

    Returns:
    --------
    Jmat(ngas,ngas) :: Jacobian matrix of chemical species.
    """

    # Initialize the Jacobian matrix with zeros
    Jmat = np.zeros((ngas, ngas), dtype=np.float64)

    eps = 1e-10

    for ir in range(nreactions):
        
        if rtype[ir] == 1:
            # photodissociations
            # or reactions a + c -> b + c
            # or reactions a + ice -> b + c
            !################################################################################
            
            ind_phot_2 = sID_pos[0, ir]
            ind_phot_4 = pID_pos[0, ir] 
            ind_phot_6 = pID_pos[1, ir]

            Jmat[ind_phot_2, ind_phot_2] -= sf[0, ir] * rrates[ilay, ir]

            if npr[ir] == 1:
                Jmat[ind_phot_4, ind_phot_2] += pf[0, ir] * rrates[ilay, ir]
            elif npr[ir] == 2:
                Jmat[ind_phot_4, ind_phot_2] += pf[0, ir] * rrates[ilay, ir]
                Jmat[ind_phot_6, ind_phot_2] += pf[1, ir] * rrates[ilay, ir]

        elif rtype[ir] == 2:
            # Reactions a + a -> b + c
            ################################################################################
            
            ind_3_2 = sID_pos[0, ir]
            ind_3_4 = pID_pos[0, ir]
            ind_3_6 = pID_pos[1, ir]

            Jmat[ind_3_2, ind_3_2] -= sf[0, ir] * rrates[ilay, ir] * c[ilay, ind_3_2]

            if npr[ir] == 1:
                Jmat[ind_3_4, ind_3_2] += pf[0, ir] * rrates[ilay, ir] * c[ilay, ind_3_2]
            elif npr[ir] == 2:
                Jmat[ind_3_4, ind_3_2] += pf[0, ir] * rrates[ilay, ir] * c[ilay, ind_3_2]
                Jmat[ind_3_6, ind_3_2] += pf[1, ir] * rrates[ilay, ir] * c[ilay, ind_3_2]

        elif rtype[ir] == 3:
            # Reactions a + b -> c + d
            ################################################################################
            
            ind_4_2 = sID_pos[0, ir]
            ind_4_4 = sID_pos[1, ir]
            ind_4_6 = pID_pos[0, ir]
            ind_4_8 = pID_pos[1, ir]

            eps_4 = abs(c[ilay, ind_4_2]) / (abs(c[ilay, ind_4_2]) + abs(c[ilay, ind_4_4]) + eps)
            eps_4 = min(eps_4, 1.0)

            Jmat[ind_4_2, ind_4_2] -= sf[0, ir] * rrates[ilay, ir] * (1.0 - eps_4) * c[ilay, ind_4_4]
            Jmat[ind_4_2, ind_4_4] -= sf[0, ir] * rrates[ilay, ir] * eps_4 * c[ilay, ind_4_2]
            Jmat[ind_4_4, ind_4_2] -= sf[1, ir] * rrates[ilay, ir] * (1.0 - eps_4) * c[ilay, ind_4_4]
            Jmat[ind_4_4, ind_4_4] -= sf[1, ir] * rrates[ilay, ir] * eps_4 * c[ilay, ind_4_2]

            if npr[ir] == 1:
                Jmat[ind_4_6, ind_4_2] += pf[0, ir] * rrates[ilay, ir] * (1.0 - eps_4) * c[ilay, ind_4_4]
                Jmat[ind_4_6, ind_4_4] += pf[0, ir] * rrates[ilay, ir] * eps_4 * c[ilay, ind_4_2]
            elif npr[ir] == 2:
                Jmat[ind_4_6, ind_4_2] += pf[0, ir] * rrates[ilay, ir] * (1.0 - eps_4) * c[ilay, ind_4_4]
                Jmat[ind_4_6, ind_4_4] += pf[0, ir] * rrates[ilay, ir] * eps_4 * c[ilay, ind_4_2]
                Jmat[ind_4_8, ind_4_2] += pf[1, ir] * rrates[ilay, ir] * (1.0 - eps_4) * c[ilay, ind_4_4]
                Jmat[ind_4_8, ind_4_4] += pf[1, ir] * rrates[ilay, ir] * eps_4 * c[ilay, ind_4_2]

        else:
            raise ValueError(f"Error: Reaction type must be 1, 2, or 3. Reaction {ir}, type {rtype[ir]}")

    return Jmat

#############################################################################################################################

@jit(nopython=True)
def locate_gas_reactions(ngas, gasID, isoID, nreactions, ns, sID, sISO, npr, pID, pISO):
    """
    Routine to find the location of the sources/products in each reaction
    in the Gas ID array defining the gases in the atmosphere.

    Inputs:
    -------
    ngas :: Number of gas species in the atmosphere.
    gasID(ngas) :: Array of gas IDs present in the atmosphere.
    isoID(ngas) :: Array of isotope IDs corresponding to the gases.
    nreactions :: Number of reactions.
    ns(nreactions) :: Number of sources in each reaction.
    sID(2,nreactions) :: Array of source gas IDs in each reaction.
    sISO(2,nreactions) :: Array of source isotope IDs in each reaction.
    npr(nreactions) :: Number of products in each reaction.
    pID(2,nreactions) :: Array of product gas IDs in each reaction.
    pISO(2,nreactions) :: Array of product isotope IDs in each reaction.

    Outputs:
    --------
    sID_pos(2,nreactions) :: Array indicating the positions of source gases in the gasID array.
    pID_pos(2,nreactions) :: Array indicating the positions of product gases in the gasID array.
    """

    # Initialize output arrays
    sID_pos = np.zeros((2, nreactions), dtype=np.int32)
    pID_pos = np.zeros((2, nreactions), dtype=np.int32)

    # Loop through each reaction
    for ir in range(nreactions):
        # Process source gases
        for j in range(ns[ir]):
            igasx = 0
            for igas in range(ngas):
                if sID[j, ir] == gasID[igas] and sISO[j, ir] == isoID[igas]:
                    sID_pos[j, ir] = igas
                    igasx = 1
                    break
            if igasx == 0:
                raise ValueError(f"Error: Reaction {ir+1}/{nreactions} involves a gas not present in the atmosphere (source). "
                                 f"GasID: {sID[j, ir]}, IsoID: {sISO[j, ir]}")

        # Process product gases
        for j in range(npr[ir]):
            igasx = 0
            for igas in range(ngas):
                if pID[j, ir] == gasID[igas] and pISO[j, ir] == isoID[igas]:
                    pID_pos[j, ir] = igas
                    igasx = 1
                    break
            if igasx == 0:
                raise ValueError(f"Error: Reaction {ir+1}/{nreactions} involves a gas not present in the atmosphere (product). "
                                 f"GasID: {pID[j, ir]}, IsoID: {pISO[j, ir]}")

    return sID_pos, pID_pos

