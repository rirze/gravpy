from __future__ import division, absolute_import, print_function
from setuptools import setup, Extension
from Cython.Build import cythonize

def get_extensions():
    from numpy import get_include
    numpy_include = get_include()
    extensions = [
        Extension('gravpy.cython.alpha.alphac', ['gravpy/cython/alpha/alphac.pyx'],
                  include_dirs=[numpy_include]),
        Extension('gravpy.cython.sie.siec', ['gravpy/cython/sie/siec.pyx'],
                  include_dirs=[numpy_include]),
    ]

    return extensions

def setup_gravpy():

    
    from numpy.distutils.misc_util import Configuration
    config_dict = Configuration('gravpy',parent_package=None,top_path=None)
    
    setup(name='gravpy',
          version='0.1',
          description='A general gravitational lens solver written in python',
          author='Sourabh Cheedella',
          author_email='cheedella.sourabh@gmail.com',
          packages=['gravpy']
          install_requires=['cython','numpy','scipy','matplotlib',],
          ext_modules= cythonize(get_extensions()),
    )

if __name__ == '__main__':
    setup_gravpy()




