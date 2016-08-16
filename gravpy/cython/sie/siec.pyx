cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt,log,atan,atanh

DTYPE   = np.float64
ctypedef np.float64_t DTYPE_t 

def elliptical(np.ndarray[DTYPE_t] x, np.ndarray[DTYPE_t] y, list modelargs_list):
    
    cdef np.ndarray[DTYPE_t] modelargs = np.asarray(modelargs_list, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] out1 = np.empty((x.size, 6), dtype=DTYPE)
    cdef int i
    for i in range(x.size):
        elliptical_single_eval(x[i],y[i],modelargs,out1[i])

    return np.transpose(out1)


def spherical(np.ndarray[DTYPE_t] x, np.ndarray[DTYPE_t] y, list modelargs_list):
    
    cdef np.ndarray[DTYPE_t, ndim=1] modelargs = np.asarray(modelargs_list, dtype=DTYPE)
    cdef int i
    cdef np.ndarray out1 = np.empty((x.size, 6), dtype=DTYPE)
    for i in range(x.size):
        spherical_single_eval(x[i],y[i],modelargs,out1[i])

    return np.transpose(out1)


@cython.cdivision(True)    
cdef void elliptical_single_eval(DTYPE_t x, DTYPE_t y, np.ndarray[DTYPE_t, ndim=1] modelargs, np.ndarray[DTYPE_t, ndim =1] output):
    cdef DTYPE_t b  = modelargs[0]
    cdef DTYPE_t x0 = modelargs[1]
    cdef DTYPE_t y0 = modelargs[2]
    cdef DTYPE_t e  = modelargs[3]
    cdef DTYPE_t te = modelargs[4]
    cdef DTYPE_t s  = modelargs[5]

    cdef DTYPE_t x2  = x**2
    cdef DTYPE_t y2  = y**2
    cdef DTYPE_t s2  = s**2
    cdef DTYPE_t q   = 1.0-e
    cdef DTYPE_t q2  = q**2
    cdef DTYPE_t om  = 1.0-q2
    cdef DTYPE_t rt  = sqrt(om)
    cdef DTYPE_t psi = sqrt(q2*(s2+x2)+y2)
    cdef DTYPE_t psis= psi + s

    cdef DTYPE_t phix = b*q/rt *atan(rt*x/psis)
    cdef DTYPE_t phiy = b*q/rt *atanh(rt*y/(psi+s*q2))

    cdef DTYPE_t invDenom = 1/(psi*(om*x2+psis**2))
    cdef DTYPE_t phixx = b*q*(psi*psis-q2*x2)*invDenom
    cdef DTYPE_t phiyy = b*q*(x2+s*psis)*invDenom
    cdef DTYPE_t phixy = -b*q*x*y*invDenom

    cdef DTYPE_t pot = b*q*s*(-0.5*log(psis**2+om*x2) + log(s*(1.0+q)) ) + x*phix+y*phiy
    
    output[0] = pot
    output[1] = phix
    output[2] = phiy
    output[3] = phixx
    output[4] = phiyy
    output[5] = phixy


@cython.cdivision(True)    
cdef void spherical_single_eval(DTYPE_t x, DTYPE_t y, np.ndarray[DTYPE_t, ndim=1] modelargs, np.ndarray[DTYPE_t, ndim =1] output):
    cdef DTYPE_t b  = modelargs[0]
    cdef DTYPE_t x0 = modelargs[1]
    cdef DTYPE_t y0 = modelargs[2]
    cdef DTYPE_t e  = modelargs[3]
    cdef DTYPE_t te = modelargs[4]
    cdef DTYPE_t s  = modelargs[5]
    
    cdef DTYPE_t r = sqrt(x**2+y**2)
    cdef DTYPE_t rad = sqrt(r**2+s**2)
    cdef DTYPE_t sprad = s + rad
    cdef DTYPE_t invDenom = 1/(rad*sprad**2)
    
    cdef DTYPE_t pot = b * (rad-s*(1+log(sprad/(2*s))))
    cdef DTYPE_t phix = b * x / sprad
    cdef DTYPE_t phiy = b * y / sprad
    cdef DTYPE_t phixx = b * (s*sprad + y**2) * invDenom
    cdef DTYPE_t phiyy = b * (s*sprad + x**2) * invDenom
    cdef DTYPE_t phixy = -b*x*y *invDenom
    
    output[0] = pot
    output[1] = phix
    output[2] = phiy
    output[3] = phixx
    output[4] = phiyy
    output[5] = phixy


