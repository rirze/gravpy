
module nfwf
  implicit none
  
contains

  subroutine nfw(res,x,y,modelargs,len)
    integer, intent(in) :: len
    real(8), intent(in) :: x(len)
    real(8), intent(in) :: y(len)
    real(8), intent(in) :: modelargs(6)
    real(8), intent(out):: res(6,len)
    !!f2py integer intent(hide), depend(x) :: len = shape(x,0)

    integer :: i
    
    do i=1,len
       call single_eval(res(:,i),x(i),y(i),modelargs)
    enddo

  end subroutine nfw

  subroutine single_eval(res,x,y,modelargs)
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    real(8), intent(in) :: modelargs(6)
    real(8), intent(out):: res(6)

    real(8) :: b,e,s
    real(8) :: r,sint,cost,front,xx,dx,phir_r,phirr,phix,phiy,phixx,phiyy,phixy,pot,t1,t2
    real(8) :: temparr(6)
    
    b  = modelargs(1)
    e  = modelargs(4)
    s  = modelargs(6)

    r = sqrt(x*x+y*y)

    if (r > 0.0d0) then
       cost = x/r
       sint = y/r
    else
       res = (/ 0.,0.,0.,0.,0.,0. /) / r ! a list of nan-s
       return
    endif

    front = (1.-e)*b

    if (abs(e) < 1.0e-6) then
       xx = r/s
       dx = xx-1.0
       phir_r = front*nfw0phir(xx)/xx

       if (abs(dx) < 1.0e-4) then
          phirr = -phir_r + 4.0*front*(1.0/3.0-0.4*dx+13.0/35.0*dx*dx-20.0/63.0*dx*dx*dx)
       else
          phirr = -phir_r + 4.0*front*(1.0-nfwFfunc(xx))/(xx*xx-1.0)
       endif

       phix  = phir_r*x
       phiy  = phir_r*y
       phixx = phir_r*sint*sint + phirr*cost*cost
       phiyy = phir_r*cost*cost + phirr*sint*sint
       phixy = (phirr-phir_r)*sint*cost

       t1 = log(xx/2)
       if (xx < 1.0e-2) then
          pot = -0.5*xx*xx*(t1+(1.0+3.0*t1)*xx*xx/8.0+(3.0/32.0+5.0/24.0*t1)*xx*xx*xx*xx)
       else if (xx <= 1.0) then
          t2   = atanh(sqrt(1.0-xx*xx))
          pot = t1*t1-t2*t2
       else
          t2   = atanh(sqrt(xx*xx-1.0))
          pot = t1*t1+t2*t2
       endif

       pot = pot * 2*front*s*s

    else
       temparr = nfw_integral(1.0d-8,1.0d0,r,sint,cost,e,s)

       pot  = front*(temparr(1)/2.)
       phix = front*(temparr(2)*x)
       phiy = front*(temparr(3)*y)
       phixx= front*(temparr(2)+2.*x*x*temparr(4))
       phiyy= front*(temparr(3)+2.*y*y*temparr(5))
       phixy= front*(2.*x*y*temparr(6))

    endif

    res = (/ pot,phix,phiy,phixx,phiyy,phixy /)

  end subroutine single_eval

  function nfw0phir(x) result(res)
    real(8), intent(in) :: x
    real(8) :: res
        
    if (x==0.0d0) then
       res = 0.0d0
       return
    endif

    if (x < 1.0e-4) then
       res = -2.0*(x*(1.0+0.75*x*x)*log(0.5*x)+0.5*x*(1.0+0.875*x*x))
    else
       res = 4.0/x*(log(0.5*x)+nfwFfunc(x))
    endif

  end function nfw0phir

  function nfwFfunc(x) result(res)
    real(8), intent(in) :: x
    real(8) :: res
    
    real(8) :: l2x,dx,tmp

    dx = x*x - 1.0
    if (abs(x)<1.0e-2) then
       l2x  = log(2.0/x)
       res  = l2x
       res  = res + x*x*(0.5*l2x-0.25)
       res  = res + x*x*x*x*(0.375*l2x-0.21875)

    else if (abs(dx)<1.0e-2) then
       res = 1.0
       res = res - dx/3.0
       res = res + dx*dx/5.0
       res = res - dx*dx*dx/7.0
       res = res + dx*dx*dx*dx/9.0
       res = res - dx*dx*dx*dx*dx/11.0
    else if (x < 1.0) then
       tmp = sqrt(1.0 -x*x)
       res = atanh(tmp)/tmp
    else
       tmp = sqrt(x*x-1.0)
       res = atan(tmp)/tmp
    endif

  end function nfwFfunc

  function nfw_integral(a,b,r,sint,cost,e,s) result(integral)
    !use cubpack                    
    real(8), intent(in) :: a,b,r,sint,cost,e,s
    real(8)  :: integral(6)
    integer  :: ndim,ncomp,maxeval,key,nregions,neval,fail
    integer  :: rgtype(1)
    real(8)  :: epsrel,epsabs
    real(8)  :: limits(1,2,1),abserr(6)
    
    ndim = 1
    ncomp = 6
    nregions = 1
    epsrel = 1d-3
    epsabs = 1d-9
    maxeval = 100000
    key = 0
    limits(1,1,1) = a
    limits(1,2,1) = b
    rgtype(1) = 1
    fail = -1

    call integrateN(nfwIntegrand, a, b, integral)
    !call cubatr(ndim,ncomp,nfwIntegrand,nregions,limits,rgtype,integral,abserr)!,&
    !      key=key, maxpts=maxeval,neval=neval, ifail=fail)
    !
    !call integrateN(nfwIntegrand,a,b,integral)
    !print *, 'Fail: ', fail
    !print *, 'NumEval: ', neval
    !print *, 'Key: ', key
    
  contains
    !function nfwIntegrand(ncomp,x) result(res)
    subroutine nfwIntegrand(x1,y,res) 
      !integer, intent(in) :: ncomp
      !real(8), intent(in) :: x(:)
      real(8) :: x1,y(6)
      real(8) :: res(6)
      !real(8) :: res(ncomp)
      
      real(8) :: q,t0,t1,t3,t5,mphiu,k,kp
      
      ! x1 = x(1)
      q  = 1.0 - e
      t0 = 1.0 - (1.0-q*q)*x1
      t1 = 1.0/sqrt(t0)
      t3 = t1/t0
      t5 = t3/t0
      
      call nfwkap(x1,mphiu,k,kp)
            
      res = (/ t1*mphiu,t1*k,t3*k,t1*kp*x1,t5*kp*x1,t3*kp*x1 /)
      
    end subroutine nfwIntegrand
    !end function nfwIntegrand

    subroutine nfwkap(u,mphiu,k,kp)
      real(8), intent(in) :: u
      real(8), intent(out):: mphiu,k,kp

      real(8) :: q,t0,t1,t2,x2,dx,phir,F

      q  = 1.0-e
      t0 = cost*cost+sint*sint/(1.0-(1.0-q*q)*u)
      t1 = t0*r*r
      t2 = t1*u
      x2 = t2/(s*s)
      dx = x2-1.0

      phir = s*nfw0phir(sqrt(x2))

      mphiu = sqrt(t2)*phir/u

      if (abs(dx) < 1.0e-4) then
         k  = 2.0/3.0-0.4*dx+2.0/7.0*dx*dx-2.0/9.0*dx*dx*dx
         kp = -0.4+4.0/7.0*dx-2.0/3.0*dx*dx+8.0/11.0*dx*dx*dx;

      else
         F = nfwFfunc(sqrt(x2))
         k  = 2.0*(1.0-F)/dx
         kp = (3.0*x2*F-2.0*x2-1.0)/(x2*s*s*dx*dx)
      endif

    end subroutine nfwkap
    
  end function nfw_integral

end module nfwf
      
    
    
    
