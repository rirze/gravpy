subroutine integrateN(func,a,b,ans)
  real*8, intent(in) :: a,b
  real*8, intent(out) :: ans(6)
  
  real*8 :: eps
  integer:: nmax,nsteps,nbad
  nmax = 5000
  nsteps = 0
  nbad = 0
  eps = 1.0d-9
  ans = 0.0d0
  
  ! print *,"got to wrapper"
  call difsub(func,a,b,ans,eps,nsteps,nbad,nmax)
  ! print*, nsteps

end subroutine integrateN

subroutine difsub(df,a,b,y,eps,neval,nbad,nmax)
  implicit none

  integer :: n 
  integer:: nsteps,neval,nbad,nmax
  real*8 :: a,e,h,hhh,one,t4,absh,haim,hsign,t1,tsave
  real*8 :: b,eps,hgus,t2,two,dum,errmax,hh,t3,zr
  real*8 :: y(6),ysave(6),dysave(6),ysi6g(6),ydub(6),y2(6)
  real*8 :: y3(6),y1(6),dydub(6),dy1(6),ysing(6)
  logical:: done,flag11

  integer :: j
  
  ! print *,"got to difsub"
  n = 6
  nsteps = 0
  done = .FALSE.
  
  hgus = (b-a)*0.01
  
  hsign = merge(1.0d0, -1.0d0, b-a > 0.0)
  
  haim = dabs(hgus)
  
  ysave = y

  tsave = a
  
  do
     !print *,"got to do loop"
     call df(tsave,ysave,dysave)
     ! print *, "make first func call"
     neval = neval + 1
     flag11 = .TRUE.

     do while (flag11)
        h = merge(haim,-haim,hsign.ge.0.0)
        if (haim .gt. dabs(b-tsave)) then
           h = b - tsave
           done = .TRUE.
           !print *, "set done to true"
        endif
        !     HERE TAKE A STEP GIVING YSING,YDUB,DYDUB. BEGIN.
        hh=0.5*h
        hhh=0.25*h
        t1=tsave+hhh
        t2=t1+hhh
        t3=t2+hhh
        t4=tsave+h

        y2=ysave+hh*dysave
        call df(t2,y2,dy1)
        neval = neval + 1
        
        y3=ysave+hh*dy1
        y2=y2+2.0*y3

        call df(t2,y3,dy1)
        neval = neval + 1
        y3 = ysave+h*dy1
        y2 = y2+y3

        call df(t4,y3,dy1)
        neval = neval + 1
        ysing = (y2-ysave+hh*dy1)/3.0
        y2    = ysave+hhh*dysave

        call df(t1,y2,dy1)
        neval = neval + 1
        y3 = ysave+hhh*dy1
        y2 = y2+2.0*y3

        call df(t1,y3,dy1)
        neval = neval + 1
        y3 = ysave+hh*dy1
        y2 = y2+y3

        call df(t2,y3,dy1)
        neval = neval + 1
        ydub = (y2-ysave+hhh*dy1)/3.0

        call df(t2,ydub,dydub)
        neval = neval + 1
        y2 = ydub+hhh*dydub

        call df(t3,y2,dy1)
        neval = neval + 1
        y3 = ydub+hhh*dy1
        y2 = y2+2.0*y3

        call df(t3,y3,dy1)
        neval = neval + 1
        y3 = ydub+hh*dy1
        y2 = y2+y3

        call df(t4,y3,dydub)
        neval = neval + 1
        errmax = 0.0d0
        absh = dabs(h)

        ydub = (y2-ydub+hhh*dydub)/3.0

        do j=1,N
           dum = merge(1.0d0,dabs(ydub(j)),dum .lt. 1.0d0)
           
           e = dabs(ydub(j)-ysing(j))/(15.0*dum)
           if (e .gt. errmax) then
              errmax = e
           endif
        enddo
        
        errmax = errmax/eps

        !print *, "err: ",errmax
        
        if (errmax .le. 0.0d0) then
           haim = absh*4.0
           flag11 = .FALSE.
           !print *,"error less than 0"
        else if (errmax .le. 1.0d0) then
           haim = absh*(errmax**-0.2)*0.90
           flag11 = .FALSE.
           !print *,"error less than eps"
        else
           haim = absh*(errmax**-0.25)*0.90
           nbad = nbad+1
           done = .FALSE.
           if (nsteps+nbad .ge. nmax) then
              print*, "nmax exceeded in while flag loop"
              return
           endif

           flag11 = .TRUE.

        endif

     enddo

     ! print *, "exited flag loop"
     nsteps = nsteps+1
     tsave = tsave + h
     if (done) then
        !print *,"got to done loop"
        !a = tsave
        !hgus = haim
        !print *,"returning y"
        
        y = (16.0*ydub-ysing)*6.666666666666d-2
                
        return
     else
        if (nsteps+nbad .ge. nmax) then
           print*, "nmax exceeded in infinite do loop"
           return
        endif
        ysave = (16.0*ydub-ysing)*6.6666666666666d-2
     endif

  enddo
end subroutine difsub
