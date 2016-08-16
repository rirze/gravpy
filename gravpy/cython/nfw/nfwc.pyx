cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt,log,atan,atanh
#from libc.math cimport abs as abs
import scipy.integrate as scipyint

DTYPE   = np.float64
ctypedef np.float64_t DTYPE_t 

@cython.boundscheck(False)
@cython.wraparound(False)
def nfw(np.ndarray[DTYPE_t] x, np.ndarray[DTYPE_t] y, list modelargs_list):
    
    cdef np.ndarray[DTYPE_t, ndim=1] modelargs = np.asarray(modelargs_list, dtype=DTYPE)
    cdef int i
    cdef np.ndarray[DTYPE_t, ndim=2] out1 = np.empty((x.size, 6), dtype=DTYPE)
    
    for i in range(x.size):
        nfw_single_eval(x[i],y[i],modelargs,out1[i])

    # cdef np.ndarray[DTYPE_t, ndim=2] absout1 = np.abs(out1)
    # out1[np.logical_or(absout1 > 10.0, absout1 < 1.0e-16)] = 0.0
    return np.transpose(out1)

#@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef void nfw_single_eval(DTYPE_t x, DTYPE_t y, np.ndarray[DTYPE_t] modelargs, np.ndarray[DTYPE_t] output):
    cdef DTYPE_t b  = modelargs[0]
    cdef DTYPE_t e  = modelargs[3]
    cdef DTYPE_t s  = modelargs[5]

    cdef DTYPE_t x2  = x*x
    cdef DTYPE_t y2  = y*y
    cdef DTYPE_t r2  = x2+y2
    cdef DTYPE_t r   = sqrt(r2)

    cdef DTYPE_t cost, sint
    if r > 0.0:
        sint = y/r
        cost = x/r
    else:
        
        output = np.array([0.0]*6)
        return
    
    cdef DTYPE_t front = (1.0 - e)*b
    cdef DTYPE_t xx, dx, phir_r, phirr, phix, phiy, phixx ,phiyy, phixy, t1, pot, t2

    cdef np.ndarray[DTYPE_t] temparr = np.empty(6, dtype=DTYPE)
    cdef int i
    if abs(e) < 1.0e-6:
        xx = r/s
        dx = xx - 1.0
        phir_r = front*nfw0phir(xx)/xx

        if abs(dx) < 1.0e-4:
            phirr = -phir_r + 4.0*front*(1.0/3.0-0.4*dx+13.0/35.0*dx*dx-20.0/63.0*dx*dx*dx)
        else:
            phirr = -phir_r + 4.0*front*(1.0-nfwFfunc(xx))/(xx*xx-1.0)

        phix  = phir_r*x
        phiy  = phir_r*y
        phixx = phir_r*sint*sint + phirr*cost*cost
        phiyy = phir_r*cost*cost + phirr*sint*sint
        phixy = (phirr-phir_r)*sint*cost
        
        t1 = log(xx/2.0)

        if xx < 1.0e-2:
            pot = -0.5*xx*xx*(t1+(1.0+3.0*t1)*xx*xx/8.0+(3.0/32.0+5.0/24.0*t1)*xx*xx*xx*xx)
        elif xx <= 1.0:
            t2  = atanh(sqrt(1.0-xx*xx))
            pot = t1*t1-t2*t2
        else:
            t2  = atanh(sqrt(xx*xx-1.0))
            pot = t1*t1+t2*t2
            
        pot = pot * 2*front*s*s

    else:
        for i in range(6):
            # temparr[i] = scipyint.quad(nfw_integrand, 1.0e-8, 1.0, args=(r,sint,cost,e,s,i))[0]
            temparr[i] = scipyint.quad(nfw_integrand, 0.0, 1.0, args=(r,sint,cost,e,s,i), full_output=1)[0]
        pot  = front*(temparr[0]/2.)
        phix = front*(temparr[1]*x)
        phiy = front*(temparr[2]*y)
        phixx= front*(temparr[1]+2.*x*x*temparr[3])
        phiyy= front*(temparr[2]+2.*y*y*temparr[4])
        phixy= front*(2.*x*y*temparr[5])


    output[0] = pot
    output[1] = phix
    output[2] = phiy
    output[3] = phixx
    output[4] = phiyy
    output[5] = phixy

    return

cdef DTYPE_t extreme_value_check(DTYPE_t num):
    if abs(num) > 10.0:
        print("trimmed inf value")
        return np.nan
    elif abs(num) < 1.0e-16:
        return 0.0
    else:
        return num

#@cython.cdivision(True)    
cdef DTYPE_t nfw0phir(DTYPE_t x):

    if x == 0.0:
        return 0.0
    elif x < 1.0e-4:
        return -2.0*(x*(1.0+0.75*x*x)*log(0.5*x)+0.5*x*(1.0+0.875*x*x))
    else:
        return 4.0/x*(log(0.5*x)+nfwFfunc(x))



#@cython.cdivision(True)    
cdef DTYPE_t nfwFfunc(DTYPE_t x):

    cdef DTYPE_t l2x, dx, tmp
    dx = x*x - 1.0
    if (abs(x)<1.0e-2): 
        l2x  = log(2.0/x)
        res  = l2x
        res  = res + x*x*(0.5*l2x-0.25)
        res  = res + x*x*x*x*(0.375*l2x-0.21875)
        
    elif (abs(dx)<1.0e-2):
        res = 1.0
        res = res - dx/3.0
        res = res + dx*dx/5.0
        res = res - dx*dx*dx/7.0
        res = res + dx*dx*dx*dx/9.0
        res = res - dx*dx*dx*dx*dx/11.0
    elif (x < 1.0):
        tmp = sqrt(1.0 -x*x)
        res = atanh(tmp)/tmp
    else:
        tmp = sqrt(x*x-1.0)
        res = atan(tmp)/tmp

    return res

#@cython.cdivision(True)    
cdef DTYPE_t nfw_integrand(DTYPE_t x, DTYPE_t r, DTYPE_t sint, DTYPE_t cost, DTYPE_t e, DTYPE_t s, int n):
#def nfw_integrand(x, r, sint, cost, e, s, int n):
    cdef DTYPE_t q,t0,t1,t3,t5,kn#,mphiu,k,kp

    q  = 1.0 - e
    t0 = 1.0 - (1.0-q*q)*x
    t1 = 1.0/sqrt(t0)
    t3 = t1/t0
    t5 = t3/t0

    cdef int n_k
    if n == 0:
        n_k = 0
    elif n == 1 or n == 2:
        n_k = 1
    else:
        n_k = 2
        
    kn = ks(x, r, sint, cost, e, s, n_k)
    
    if n==0:
        return t1*kn
    elif n==1:
        return t1*kn
    elif n==2:
        return t3*kn
    elif n==3:
        return t1*kn*x
    elif n==4:
        return t5*kn*x
    elif n==5:
        return t3*kn*x        
    
#@cython.cdivision(True)    
cdef DTYPE_t ks(DTYPE_t u, DTYPE_t r, DTYPE_t sint, DTYPE_t cost, DTYPE_t e, DTYPE_t s, int n):
    cdef DTYPE_t q,t0,t1,t2,x2,dx,phir,F

    q  = 1.0-e
    t0 = cost*cost+sint*sint/(1.0-(1.0-q*q)*u)
    t1 = t0*r*r
    t2 = t1*u
    x2 = t2/(s*s)
    dx = x2-1.0

    if n == 0:
        x2 = t2/(s*s)
        phir = s*nfw0phir(sqrt(x2))

        return sqrt(t2)*phir/u

    if abs(dx) < 1.0e-4:
        if n == 1:
            k  = 2.0/3.0-0.4*dx+2.0/7.0*dx*dx-2.0/9.0*dx*dx*dx
            return k
        if n == 2:
            kp = -0.4+4.0/7.0*dx-2.0/3.0*dx*dx+8.0/11.0*dx*dx*dx;
            return kp

    else:
        F  = nfwFfunc(sqrt(x2))
        if n == 1:
            k  = 2.0*(1.0-F)/dx
            return k
        if n == 2:
            kp = (3.0*x2*F-2.0*x2-1.0)/(x2*s*s*dx*dx)
            return kp

        

