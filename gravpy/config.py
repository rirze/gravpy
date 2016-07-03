from __future__ import absolute_import

from gravpy.python.sie import siepy
from gravpy.python.alpha import alphapy

use_fortran_modules = True
if use_fortran_modules:
    # f2py modules imported below, clunky pathing due to fortran modules, full
    # function path is fmodel.fmodel.routine
    from gravpy.fortran.sie.sief import sief
    from gravpy.fortran.alpha.alphaf import alphaf
    from gravpy.fortran.nfw.nfwf import nfwf

sie_spherical = sief.spherical
sie_ellipitical = sief.elliptical

alpha_plummer = alphapy.plummer
alpha_general = alphaf.general

nfw_general = nfwf.nfw

