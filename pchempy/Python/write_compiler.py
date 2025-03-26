import numpy as np
import sys,os
from pchempy import *

#Set of routines to write the photochemical model online

##################################################################################################
##################################################################################################
#                                         COMPILER
##################################################################################################
##################################################################################################

def write_compiler_fortran(runname,model_name='pmodel'):
    
    """
        FUNCTION NAME : write_compiler_fortran()
        
        DESCRIPTION : Routine to write the compiler of the fortran code
        
        INPUTS :
        
            runname :: Name of the simulation run
            
        OPTIONAL INPUTS:
        
        OUTPUTS :
        
        CALLING SEQUENCE:
        
            write_compiler_fortran(runname)
        
        MODIFICATION HISTORY : Juan Alday (15/12/2023)
        
    """
    
    from pchempy.Python.utils import pchempy_path

    comp = "f2py -c -m pmodel "
    
    #Shared routines with pchempy
    pchem_path = pchempy_path()
    comp += pchem_path+'Fortran/photolysis.f90 '
    comp += pchem_path+'Fortran/photolysis_mod.f90 '
    comp += pchem_path+'Fortran/chemistry.f90 '
    comp += pchem_path+'Fortran/reactions.f90 '
    
    #Specific routines to this model
    comp += runname+"_diffusion.f90 "
    comp += runname+"_chemical_network.f90 "
    comp += runname+"_converge.f90 "
    
    #Flags
    comp += "-llapack --fcompiler='gfortran' --opt='-fcheck=bounds'"
 
    f = open(model_name+'.comp','w')
    f.write(comp)
    f.close()    

