import numpy as np


#modelargs: (major) radius, x-center position, y-center position, ellipticity, ellipticity angle, core radius

        
def elliptical(x,y,modelargs):
    b,x0,y0,e,te,s  = modelargs[:6]
    
    x2  = x*x
    y2  = y*y
    s2  = s*s
    q   = 1.0-e
    q2  = q*q
    om  = 1.0-q2
    rt  = np.sqrt(om)
    psi = np.sqrt(q2*(s2+x2)+y2) 
    psis= psi + s


    phix = b*q/rt *np.arctan(rt*x/psis)
    phiy = b*q/rt *np.arctanh(rt*y/(psi+s*q2))

    invDenom = b*q/(psi*(om*x2+psis*psis))
    phixx = (psi*psis-q2*x2)*invDenom
    phiyy = (x2+s*psis)*invDenom
    phixy = -x*y*invDenom

    pot = b*q*s*(-0.5*np.log(psis*psis+om*x2) + np.log(s*(1.0+q)) ) + x*phix+y*phiy
    
    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

def spherical(x,y,modelargs):
    b,x0,y0,e,te,s  = modelargs[:6]

    rad = np.sqrt(x*x+y*y+s*s)
    sprad = s + rad
    invDenom = b/(rad*sprad*sprad)

    pot = b * (rad-s*(1+np.log(sprad/(2*s))))
    phix = b * x / sprad
    phiy = b * y / sprad
    phixx = (s*sprad + y*y) * invDenom
    phiyy = (s*sprad + x*x) * invDenom
    phixy = -x*y *invDenom

    return np.array((pot,phix,phiy,phixx,phiyy,phixy))

