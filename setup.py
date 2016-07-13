from __future__ import division, absolute_import, print_function
from setuptools import setup, Extension
from Cython.Build import cythonize
    
def setup_gravpy():

    import numpy
    from numpy.distutils.misc_util import Configuration
    config_dict = Configuration('gravpy',parent_package=None,top_path=None)

    numpy_include = numpy.get_include()
    extensions = [
        Extension('gravpy.cython.alpha.alphac', ['gravpy/cython/alpha/alphac.pyx'],
                  include_dirs=[numpy_include]),
        Extension('gravpy.cython.sie.siec', ['gravpy/cython/sie/csie.pyx'],
                  include_dirs=[numpy_include]),
    ]
    
    setup(name='gravpy',
          version='0.1',
          description='A general gravitational lens solver written in python',
          author='Sourabh Cheedella',
          author_email='cheedella.sourabh@gmail.com',
          install_requires=['numpy','scipy','matplotlib','cython'],
          ext_modules= cythonize(extensions),
    )

if __name__ == '__main__':
    setup_gravpy()




