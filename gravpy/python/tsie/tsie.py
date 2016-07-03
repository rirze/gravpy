from math import sqrt
import numpy as np
from numpy import pi

def tsis(x,y,modelargs):
    b,x0,y0,e,te,s,a = modelargs

    r2 = x*x + y*y
    r = sqrt(r2)

    if r == 0.0:
        return [float('nan') for i in range(6)]

    cost = x/r
    sint = y/r

    if abs(a) < 0.01:
        if r < 0.01*a:
            ptmp = r*(1.0-r/(2.0*pi*a)*(1.0+r*r/(24.0*a*a)+r*r*r*r/(120.0*a*a*a*a*a)))
            atmp = (1.0-r/(pi*a)*(1.0+r*r/(12.0*a*a)+r*r*r*r/(40.0*a*a*a*a)))
            ktmp = 1.0/(2.0*r)*(1.0-2.0*r/(pi*a)*(1.0+r*r/(6.0*a*a)+3.0*r*r*r*r/(40.0*a*a*a*a)))
            gtmp = 1.0/(2.0*r)*(1.0+r*r*r/(6.0*pi*a*a*a)*(1.0+0.6*r*r/(a*a)))
        elif r < a: 
            dum  = sqrt(a*a/rsq-1.0)
            ptmp = 2.0/pi*(2.0*a-2.0*r*dum+r*atan(dum)+a*log((a+r*dum)/(2.0*a)))
            atmp = 2.0/(pi*r)*(a-r*dum+r*atan(dum))
            ktmp = 1.0/(pi*r)*atan(dum)
            gtmp = 1.0/(pi*rsq)*(2.0*a-2.0*r*dum+r*atan(dum))
        else: 
            ptmp = 2.0*a/pi*(2.0+log(0.5*r/a))
            atmp = 2.0*a/(pi*r)
            ktmp = 0.0
            gtmp = 2.0*a/(pi*rsq)
        
        ptmp -= 2.0*a/pi*(2.0-log(2.0*a))

    elif abs(a-1) < 0.01:
        rs   = 2/pi*a
        xi   = r/rs
        xi2  = xi*xi
        xi3  = xi2*xi
        xi4  = xi2*xi2
        xi5  = xi3*xi2
        xi6  = xi4*xi2
        xi8  = xi6*xi2
        xi10 = xi8*xi2
        xi12 = xi10*xi2
        op   = 1 + xi2
        rt   = sqrt(op)
        if xi > 1.0e-3:
            ptmp = rs*(1 + xi - rt + log((1 + rt)/2))
            atmp = (1 + xi - rt)/xi
            ktmp = (1/xi - 1/rt)/(2*rs)
            gtmp = (2 + xi - 2/rt - xi2/rt)/(2*rs*xi2)
        else: 
            ptmp = rs*xi - (rs*xi2)/4 + (rs*xi4)/32
            atmp = 1 - xi/2 + xi3/8
            ktmp = -1/(2*rs) + 1/(2*rs*xi) + xi2/(4*rs)
            gtmp = 1/(2*rs*xi) - xi2/(8*rs)
            
    elif abs(a-2) < 0.01:
        rs   = 4/pi*a
        xi   = r/rs
        xi2  = xi*xi
        xi3  = xi2*xi
        xi4  = xi2*xi2
        xi5  = xi3*xi2
        xi6  = xi4*xi2
        xi8  = xi6*xi2
        xi10 = xi8*xi2
        xi12 = xi10*xi2
        op   = 1 + xi2
        rt   = sqrt(op)
        
        if xi>1.0e-3:
            ptmp = (rs*(2 + 2*xi - 2*rt + log((1 + rt)/2)))/2
            atmp = (2 - (2*xi)/rt + (1 - 1/rt)/xi)/2
            ktmp = (-3*xi - 2*xi3 + 2*rt + 2*xi2*rt)/
            (4*rs*xi*op*rt)
            gtmp = (-2*xi4 + 2*xi*rt + 2*xi3*rt + 
                    2*(-1 + rt) + xi2*(-3 + 2*rt))/
            (4*rs*xi2*op*rt)
        else:
            ptmp = rs*xi - (3*rs*xi2)/8 + (5*rs*xi4)/64
            atmp = 1 - (3*xi)/4 + (5*xi3)/16
            ktmp = -3/(4*rs) + 1/(2*rs*xi) + (5*xi2)/(8*rs)
            gtmp = 1/(2*rs*xi) - (5*xi2)/(16*rs)

    elif abs(a-3.0) < 0.01:
        rs   = 16/(3*pi)*a
        xi   = r/rs
        xi2  = xi*xi
        xi3  = xi2*xi
        xi4  = xi2*xi2
        xi5  = xi3*xi2
        xi6  = xi4*xi2
        xi8  = xi6*xi2
        xi10 = xi8*xi2
        xi12 = xi10*xi2
        op   = 1 + xi2
        rt   = sqrt(op)
        if (xi>1.0e-3): 
            ptmp = (rs*(7 + 8*xi + 1/rt - 8*rt - log(8) +
                        3*log(1 + rt)))/8
            atmp = (3 + 8*xi + 1/(op*rt) + 4/rt - 8*rt)/
            (8*xi)
            ktmp = (-15*xi - 20*xi3 - 8*xi5 + 8*rt + 
                    16*xi2*rt + 8*xi4*rt)/(16*rs*xi*op*op*rt)
            gtmp = -(6 + 15*xi2 + 20*xi4 + 8*xi6 - 6*op*op*rt - 
                     8*xi*op*op*rt)/(16*rs*xi2*op*op*rt)
        else:
            ptmp = rs*xi - (15*rs*xi2)/32 + (35*rs*xi4)/256
            atmp = 1 - (15*xi)/16 + (35*xi3)/64
            ktmp = -15/(16*rs) + 1/(2*rs*xi) + (35*xi2)/(32*rs)
            gtmp = 1/(2*rs*xi) - (35*xi2)/(64*rs)
            
  
 
    elif (fabs(al-4.0)<0.01):
        rs   = 32/(5*pi)*a
        xi   = r/rs
        xi2  = xi*xi
        xi3  = xi2*xi
        xi4  = xi2*xi2
        xi5  = xi3*xi2
        xi6  = xi4*xi2
        xi8  = xi6*xi2
        xi10 = xi8*xi2
        xi12 = xi10*xi2
        op   = 1 + xi2
        rt   = sqrt(op)
        if (xi>1.0e-3): 
            ptmp = -(rs*(38 + 48*xi4 - 38*rt - 48*xi*rt - 
                         48*xi3*rt + rt*log(32768) + 
                         xi2*(87 - 38*rt + rt*log(32768)) - 
                         15*op*rt*log(1 + rt)))/(48*op*rt)
            atmp = (5 + 16*xi + 1/(op*op*rt) + 2/(op*rt) + 
                    8/rt - 16*rt)/(16*xi)
            ktmp = (16 - (xi*(35 + 70*xi2 + 56*xi4 + 16*xi6))/(op*op*op*rt))/
            (32*rs*xi)
            gtmp = -(10 + 35*xi2 + 70*xi4 + 56*xi6 + 16*xi8 - 10*op*op*op*rt - 
                     16*xi*op*op*op*rt)/(32*rs*xi2*op*op*op*rt)
        else:
            ptmp = rs*xi - (35*rs*xi2)/64 + (105*rs*xi4)/512
            atmp = 1 - (35*xi)/32 + (105*xi3)/128
            ktmp = -35/(32*rs) + 1/(2*rs*xi) + (105*xi2)/(64*rs)
            gtmp = 1/(2*rs*xi) - (105*xi2)/(128*rs)
            
  
 
    elif (fabs(al-5.0)<0.01):
        rs   = 256/(35*pi)*a
        xi   = r/rs
        xi2  = xi*xi
        xi3  = xi2*xi
        xi4  = xi2*xi2
        xi5  = xi3*xi2
        xi6  = xi4*xi2
        xi8  = xi6*xi2
        xi10 = xi8*xi2
        xi12 = xi10*xi2
        op   = 1 + xi2
        rt   = sqrt(op)
        if (xi>1.0e-3): 
            ptmp = (rs*(-281 - 965*xi2 - 1065*xi4 - 384*xi6 + 281*rt + 
                        384*xi*rt + 562*xi2*rt + 
                        768*xi3*rt + 281*xi4*rt + 
                        384*xi5*rt - 105*op*op*rt*log(2*rs) + 
                        105*op*op*rt*log(rs*(1 + rt))))/(384*op*op*rt)
            atmp = (128*rs*xi - (rs*(35 + 280*xi2 + 560*xi4 + 448*xi6 + 128*xi8 - 
                                     35*op*op*op*rt))/(op*op*op*rt))/(128*rs*xi)
            ktmp = (128 - (xi*(315 + 840*xi2 + 1008*xi4 + 576*xi6 + 128*xi8))/
                    (op*op*op*op*rt))/(256*rs*xi)
            gtmp = -(70 + 315*xi2 + 840*xi4 + 1008*xi6 + 576*xi8 + 128*xi10 - 
                     70*op*op*op*op*rt - 128*xi*op*op*op*op*rt)/
            (256*rs*xi2*op*op*op*op*rt)
        else:
            ptmp = rs*xi - (315*rs*xi2)/512 + (1155*rs*xi4)/4096
            atmp = 1 - (315*xi)/256 + (1155*xi3)/1024
            ktmp = -315/(256*rs) + 1/(2*rs*xi) + (1155*xi2)/(512*rs)
            gtmp = 1/(2*rs*xi) - (1155*xi2)/(1024*rs)
            
  
 
    elif (fabs(al-6.0)<0.01):
        rs   = 512/(63*pi)*a
        xi   = r/rs
        xi2  = xi*xi
        xi3  = xi2*xi
        xi4  = xi2*xi2
        xi5  = xi3*xi2
        xi6  = xi4*xi2
        xi8  = xi6*xi2
        xi10 = xi8*xi2
        xi12 = xi10*xi2
        op   = 1 + xi2
        rt   = sqrt(op)
        if (xi>1.0e-3): 
            ptmp = (rs*(-315*log(2*rs) + (-878 - 4018*xi2 - 6650*xi4 - 4795*xi6 - 
                                          1280*xi8 + 878*op*op*op*rt + 1280*xi*op*op*op*rt + 
                                          315*op*op*op*rt*log(rs*(1 + rt)))/(op*op*op*rt)))/
            1280
            atmp = (256*rs*xi - (rs*(63 + 630*xi2 + 1680*xi4 + 2016*xi6 + 
                                     1152*xi8 + 256*xi10 - 63*op*op*op*op*rt))/(op*op*op*op*rt))/
            (256*rs*xi)
            ktmp = (256 - (xi*(693 + 2310*xi2 + 3696*xi4 + 3168*xi6 + 1408*xi8 + 
                               256*xi10))/(op*op*op*op*op*rt))/(512*rs*xi)
            gtmp = -(126 + 693*xi2 + 2310*xi4 + 3696*xi6 + 3168*xi8 + 1408*xi10 + 
                     256*xi12 - 126*op*op*op*op*op*rt - 256*xi*op*op*op*op*op*rt)/
            (512*rs*xi2*op*op*op*op*op*rt)
        else:
            ptmp = rs*xi - (693*rs*xi2)/1024 + (3003*rs*xi4)/8192
            atmp = 1 - (693*xi)/512 + (3003*xi3)/2048
            ktmp = -693/(512*rs) + 1/(2*rs*xi) + (3003*xi2)/(1024*rs)
            gtmp = 1/(2*rs*xi) - (3003*xi2)/(2048*rs)
            
    else:
        raise ValueError('unknown value of truncation power in tsis model')

    ptmp *= b
    atmp *= b
    ktmp *= b
    gtmp *= b

    phiarr = np.zeros(6)
    
    phiarr[ 1] = ptmp
    phiarr[ 2] = atmp*cost
    phiarr[ 3] = atmp*sint
    phiarr[ 4] = ktmp - gtmp*(cost*cost-sint*sint)
    phiarr[ 5] = ktmp + gtmp*(cost*cost-sint*sint)
    phiarr[ 6] = -gtmp*2.0*cost*sint

    return phiarr

