cimport cython
import numpy as np
cimport numpy as np
import scipy.integrate as scipyint
import scipy.special as sp
from libc.math cimport sqrt,log

DTYPE   = np.float64
ctypedef np.float64_t DTYPE_t 

@cython.boundscheck(False)
@cython.wraparound(False)
def general(np.ndarray[DTYPE_t] x, np.ndarray[DTYPE_t] y, list modelargs_list):
    
    cdef np.ndarray[DTYPE_t] modelargs = np.asarray(modelargs_list, dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] out1 = np.empty((x.size, 6), dtype=DTYPE)
    cdef int i
    for i in range(x.size):
        single_eval(x[i],y[i],modelargs,out1[i])

    return np.transpose(out1)

cdef void single_eval(DTYPE_t x, DTYPE_t y, np.ndarray[DTYPE_t, ndim=1] modelargs, np.ndarray[DTYPE_t] out):
    cdef DTYPE_t b  = modelargs[0]
    cdef DTYPE_t e  = modelargs[3]
    cdef DTYPE_t s  = modelargs[5]
    cdef DTYPE_t a  = modelargs[6]
    
    cdef DTYPE_t r = sqrt(x*x+y*y)
    cdef DTYPE_t cost = x/r if r>0. else 1.0
    cdef DTYPE_t sint = y/r if r>0. else 0.0
    cdef DTYPE_t front = (1-e)/2.0*(b)**(2.0-a)

    cdef DTYPE_t pot, phix, phiy, phixx, phiyy, phixy
    cdef DTYPE_t phir, phir_r, phirr
    cdef DTYPE_t p0,p1,p2,p3,p4,p5

    if abs(a) > 0.01:
        front *= a

    if abs(e) < 1.0e-6:
        phir   = front*alpha0phir(r,s,a)
        phir_r = front*pow(s,a-2.0) if r==0.0 else phir/r
        phirr  = -phir_r + 2.0*front*(s*s+r*r)**(0.5*a-1.0)
        phix   = phir_r*x
        phiy   = phir_r*y
        phixx  = phir_r*sint*sint + phirr*cost*cost
        phiyy  = phir_r*cost*cost + phirr*sint*sint
        phixy  = (phirr-phir_r)*sint*cost

        if abs(s) <= 2.0*r:
            pot  = sp.hyp2f1(-0.5*a,-0.5*a,1.0-0.5*a,-(s*s)/(r*r))
            pot *= b*b*(r/b)**a/(a*a)
            pot -= b*b*(s/b)**a/a*(log(r/s)+0.5*(0.577216-sp.digamma(-0.5*a)))
            if abs(a) > 0.01:
                pot *= a

        else:
            pot = 0.0 if r==0.0 else front*scipyint.quad(alpha0phir,0.0,r,args=(s,a))[0]

    else:

        p0,p1,p2,p3,p4,p5 = alpha_integral(0.0,1.0,r,e,sint,cost,s,a)
        pot  = front*(p0/2.)
        phix = front*(p1*x)
        phiy = front*(p2*y)
        phixx= front*(p1+2.*x*x*p3)
        phiyy= front*(p2+2.*y*y*p4)
        phixy= front*   (2.*x*y*p5)

    out[0] = pot
    out[1] = phix
    out[2] = phiy
    out[3] = phixx
    out[4] = phiyy
    out[5] = phixy
    return

cdef DTYPE_t alpha0phir(r,s,a):
    if r/s < 1.0e-4:
        return r*s**(a-2.0)
    elif a==0.0:
        return log(1.0+r*r/(s*s))/r
    else:
        return 2.0/(a*r)*((s*s+r*r)**(a/2.0) - s**a)
    
cpdef alpha_integral(lower,upper,r,e,sint,cost,s,a):
    # return [scipyint.quad(alpha_integrand,lower,upper,args=(r,e,sint,cost,s,a,num))[0] for num in range(6)]
    return [scipyint.quad(function,lower,upper,args=(r,e,sint,cost,s,a,num))[0] for num in range(6)]

cdef DTYPE_t function(DTYPE_t u, DTYPE_t r, DTYPE_t e, DTYPE_t sint, DTYPE_t cost, DTYPE_t s, DTYPE_t a, int ind):
    cdef DTYPE_t q  = 1.0-e
    cdef DTYPE_t t0 = 1.0-(1.0-q*q)*u
    cdef DTYPE_t t1 = 1.0/sqrt(t0)
    cdef DTYPE_t t3 = t1/t0
    cdef DTYPE_t t5 = t3/t0
    cdef DTYPE_t t6 = (cost*cost+sint*sint/(1.0 - (1.0-q*q)*u))*r*r
    cdef DTYPE_t t7 = t6*u

    if ind==0:
        mphiu = sqrt(t7)*alpha0phir(np.sqrt(t7),s,a)/u if u!=0.0 else t6*s**(a-2.0)
        return t1*mphiu

    cdef DTYPE_t k  = (s*s+t7)**(a/2.0-1.0)
    if ind == 1:
        return t1*k
    if ind == 2:
        return t3*k

    cdef DTYPE_t kp = k*(a/2.0-1.0)/(s*s+t7)
    if ind == 3:
        return t1*kp*u
    if ind == 4:
        return t5*kp*u
    if ind == 5:
        return t3*kp*u

    else:
        return ind

    






    
    
