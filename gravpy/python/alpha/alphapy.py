import numpy as np
import scipy.special as sp
import scipy.integrate as scipyint

import itertools
import warnings
import ctypes
import pkg_resources
path  = '/home/chronos/Research/lens-solver/gravpy/cython/alpha/libintegrand.so'
lib = ctypes.CDLL(path)
lib.alpha_integrand.restype = ctypes.c_double
lib.alpha_integrand.argtypes= (ctypes.c_int, ctypes.c_double)

def plummer(x,y,modelargs):
    '''Calculation for the alpha=-1 case'''
    b,x0,y0,e,te,s  = modelargs
    
    x2 = x*x
    y2 = y*y
    s2 = s*s
    q  = 1.0-e
    q2 = q*q
    om = 1.0-q2
    rt = np.sqrt(om)
    front = b**3*q/s
    psi = np.sqrt(q2*(s2+x2)+y2)
    psi2 = psi*psi
    psi3 = psi2*psi
    
    invDenom = 1/(psi*(om+x2+(psi+s)**2))
    phix = front*x*(psi+s*q2)*invDenom
    phiy = front*y*(psi+s)   *invDenom
    
    tmp = (psi+s)**2+om*x2

    # phixx
    phixx = -2.0*front*x2*((psi+q2*s)/tmp)**2/psi2 + front*(psi2 *(psi+q2*s)-q2**2*x2*s)/(psi3*tmp)
    phiyy = -2.0*front*y2*((psi+   s)/tmp)**2/psi2 + front*(psi2 *(psi   +s)-      y2*s)/(psi3*tmp)
    phixy = -front*x*y*(2.0*(psi+s)*(psi+q2*s)/(psi2*tmp*tmp) + q2*s/(psi3*tmp))
    pot   = 0.5*front*np.log((psi*s)**2+om*x2) - front*np.log(s*(1.0+np.fabs(q)))
    
    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

def general(x,y,modelargs):
    b,x0,y0,e,te,s,a  = modelargs
    
    def single_eval(x,y):
        r = np.sqrt(x*x+y*y)
        cost = x/r if r>0. else 1.0
        sint = y/r if r>0. else 0.0
        front = (1-e)/2.0*np.power(b,2.0-a)
        
        if abs(a) > 0.01:
            front *= a
            
        if abs(e) < 1.0e-6:
            phir   = front*alpha0phir(r,s,a)
            phir_r = front*np.power(s,a-2.0) if r==0.0 else phir/r
            phirr  = -phir_r + 2.0*front*np.power(s*s+r*r,0.5*a-1.0)
            phix   = phir_r*x
            phiy   = phir_r*y
            phixx  = phir_r*sint*sint + phirr*cost*cost
            phiyy  = phir_r*cost*cost + phirr*sint*sint
            phixy  = (phirr-phir_r)*sint*cost
                
            if abs(s) <= 2.0*r:
                warnings.warn("Calling hypergeo and gamma functions")
                pot  = sp.hyp2f1(-0.5*a,-0.5*a,1.0-0.5*a,-(s*s)/(r*r))
                pot *= b*b*np.power(r/b,a)/(a*a)
                pot -= b*b*np.power(s/b,a)/a*(np.log(r/s)+0.5*(0.577216-sp.digamma(-0.5*a)))
                if abs(a) > 0.01:
                    pot *= a
                        
            else:
                pot = 0.0 if r==0.0 else front*scipyint.quad(alpha0phir,0.0,r,args=(s,a))[0]
                        
        else:
            warnings.warn("Evaluating integrals")
            p0,p1,p2,p3,p4,p5 = c_alpha_integral(0.0,1.0,r,e,sint,cost,s,a)
            pot  = front*(p0/2.)
            phix = front*(p1*x)
            phiy = front*(p2*y)
            phixx= front*(p1+2.*x*x*p3)
            phiyy= front*(p2+2.*y*y*p4)
            phixy= front*   (2.*x*y*p5)
                                        
        return (pot,phix,phiy,phixx,phiyy,phixy)

    return np.transpose([single_eval(x,y) for x,y in itertools.izip(x,y)])


def alpha0phir(r,s,a):
    if r/s < 1.0e-4:
        return r*np.power(s,a-2.0)
    elif a==0.0:
        return np.log(1.0+r*r/(s*s))/r
    else:
        return 2.0/(a*r)*(np.power(s*s+r*r,a/2.0) - np.power(s,a))


def alpha0phir_for(r,s,a):
    return [alpha0phir(r_i, s, a) for r_i in r]

def alpha_integral(lower,upper,r,e,sint,cost,s,a):
    def function(u):
        q  = 1.0-e
        t0 = 1.0-(1.0-q*q)*u
        t1 = 1.0/np.sqrt(t0)
        t3 = t1/t0
        t5 = t3/t0
        t6 = (cost*cost+sint*sint/(1.0 - (1.0-q*q)*u))*r*r
        t7 = t6*u

        
        # mphiu = t6*np.power(s,a-2.0) if u==0.0 else np.sqrt(t7)*alpha0phir(np.sqrt(t7),s,a)/u
        mphiu = np.where( u != 0.0, np.sqrt(t7)*alpha0phir_for(np.sqrt(t7),s,a)/u, t6*np.power(s,a-2.0))
        k  = np.power(s*s+t7,a/2.0-1.0)
        kp = k*(a/2.0-1.0)/(s*s+t7)

        
        out = np.squeeze([t1*mphiu, t1*k, t3*k, t1*kp*u, t5*kp*u, t3*kp*u]).T
        #print out.shape
        return out

    #return [scipyint.quad(function,lower,upper,args=(num))[0] for num in range(6)]
    val, err = cubature(function, 1, 6, [lower], [upper], adaptive='p', vectorized=True)
    return val


def c_alpha_integral(lower,upper,r,e,sint,cost,s,a):
    return [scipyint.quad(lib.alpha_integrand,lower,upper,args=(r,e,sint,cost,s,a,num))[0]
            for num in range(6)]
