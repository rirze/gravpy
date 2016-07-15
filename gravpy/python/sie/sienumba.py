import numba
import gravpy.python.sie.siepy as siepy

spherical = numba.jit("f8[:,:](f8[:],f8[:],f8[:])")(siepy.spherical)

elliptical = numba.jit("f8[:,:](f8[:],f8[:],f8[:])")(siepy.elliptical)
