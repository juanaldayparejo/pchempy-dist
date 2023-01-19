#from setuptools import setup
from numpy.distutils.core import setup
from numpy.distutils.core import Extension

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
