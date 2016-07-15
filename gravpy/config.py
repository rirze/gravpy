from __future__ import absolute_import

from gravpy.python.sie import siepy, sienumba
from gravpy.python.alpha import alphapy

from gravpy.cython.alpha import alphac
from gravpy.cython.sie import siec


try:
    # f2py modules imported below, clunky pathing due to fortran modules, full
    # function path is fmodel.fmodel.routine
    from gravpy.fortran.sie.sief import sief
    from gravpy.fortran.alpha.alphaf import alphaf
    from gravpy.fortran.nfw.nfwf import nfwf
except:
    pass


sie_spherical = sienumba.spherical
sie_ellipitical = sienumba.elliptical

alpha_plummer = alphapy.plummer
alpha_general = alphac.general

nfw_general = None #nfwf.nfw

