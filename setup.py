#from setuptools import setup
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
import os

#Installation path
install_path = os.path.dirname(os.path.abspath(__file__))

#Writing installation path into path.h file
f = open(install_path+'/pchempy/Fortran/datapath.h','w')
def1str = 'character(len=100), parameter :: lp1='
def2str = 'character(len=100), parameter :: lp2='
def3str = 'character(len=100), parameter :: lp3='
defstr = 'character(len=300), parameter :: datapath=trim(lp1)//trim(lp2)//trim(lp3)'
defpath = install_path+'/pchempy/Data/'
iin = int(len(defpath)/3)
f.write(def1str+'"'+defpath[0:iin]+'"\n')
f.write(def2str+'"'+defpath[iin:2*iin]+'"\n')
f.write(def3str+'"'+defpath[2*iin:len(defpath)]+'"\n')
f.write(defstr)
f.close()
sys.exit()

#Compiling fortran code
ext1 = Extension(name='pchempy.pchemf',
                 sources=['pchempy/Fortran/converge_model.f90','pchempy/Fortran/diffusion.f90','pchempy/Fortran/chemistry.f90','pchempy/Fortran/photolysis_mod.f90','pchempy/Fortran/photolysis.f90',\
                            'pchempy/Fortran/mars_chemistry.f90'],  #Mars specific case
                 libraries=['lapack'],
                 f2py_options=[''],
                 )

#Setting up the requirements to other python libraries
setup(name='pchempy',
      version='1.0.0',
      description='Python 1D photochemical model for planetary atmospheres.',
      packages=['pchempy'],
      install_requires=['numpy','matplotlib','h5py'],
      ext_modules=[ext1],
      )
