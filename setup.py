from __future__ import division, absolute_import, print_function
from setuptools import find_packages


try:
    from Cython.Build import cythonize
except:
    raise ImportWarning("Cython not installed")

try:
    from numpy import get_include
    from numpy.distutils.core import setup, Extension
except:
    raise ImportWarning("Numpy not installed")

def get_extensions():
    numpy_include = get_include()
    extensions = [
        # C files
        Extension('gravpy.cython.sie.siec', ['gravpy/cython/sie/siec.pyx'], include_dirs=[numpy_include]),
        Extension('gravpy.cython.alpha.alphac', ['gravpy/cython/alpha/alphac.pyx'], include_dirs=[numpy_include]),
        Extension('gravpy.cython.nfw.nfwc', ['gravpy/cython/nfw/nfwc.pyx'], include_dirs=[numpy_include]),
        # Fortran files
        Extension('gravpy.fortran.sie.sief', ['gravpy/fortran/sie/sief.f90']),
        Extension('gravpy.fortran.alpha.alphaf', ['gravpy/fortran/alpha/alphaf.f90', 'gravpy/fortran/alpha/asa103.f90', 'gravpy/fortran/alpha/hyp.f']),
        Extension('gravpy.fortran.nfw.nfwf', ['gravpy/fortran/nfw/nfwf.f90', 'gravpy/fortran/nfw/difsub.f90']),
    ]

    return extensions

def setup_gravpy():

    # from numpy.distutils.core import Extension as npext
    # f_modules = [npext(name= 'gravpy.fortran.sie.sief', sources = ['gravpy/fortran/sie/sief.f90'])]
    modules = cythonize(get_extensions()) # + f_modules
    
    setup(name='gravpy',
          version='0.1',
          description='A general gravitational lens solver written in python',
          author='Sourabh Cheedella',
          author_email='cheedella.sourabh@gmail.com',
          packages=find_packages(),
          install_requires=['cython','numpy','scipy','matplotlib'],
          ext_modules= modules,

    )

if __name__ == '__main__':
    setup_gravpy()




