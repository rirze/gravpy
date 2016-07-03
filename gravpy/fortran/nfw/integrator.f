      SUBROUTINE integrateN(FUNC,A,B,ANS)
      
      real*8, intent(in) :: A,B
      real*8, intent(out) :: ANS(6)
      
      
      real*8 :: eps,hgus
      integer:: nmax,nsteps,nbad
      nmax = 5000
      nsteps = 0
      nbad = 0
      hgus = 0.0
      eps = 1.0d-6
      
      call DIFSUB(FUNC,A,B,ANS,eps,hgus,nsteps,nbad,nmax)
      print*, nsteps
      
      end SUBROUTINE integrateN


      SUBROUTINE DIFSUB(DF,A,B,Y,EPS,HGUS,NSTEPS,NBAD,NMAX)
c      implicit none
      implicit real*8 (a-h,o-z)
C
C     Note that after a call, the routine is ready to continue integrating
C     from the point it left off.
C     DF = SUBROUTINE DF(X,Y,DY). Computes rhs of DY = F(X,Y). Must be
C          declared external in calling program. Y and DY must be dimensioned.
C     NN = No. of eqns. If NN > 50, change dimension statement.
C     A  = Starting value of X. Updated to B on return.
C     B  = Final value of X. B < A is acceptable.
C     Y  = Solution vector. Must be initialized.
C     EPS= Accuracy parameter. It is a relative error tolerance per step for
C          Y's greater than 1, an absolute error tolerance for those less than
C          1. (see coding near statement 84).
C     HGUS=Initial guess for stepsize. If set to 0, uses default of .01*(B-A).
C          Returns last stepsize used.
C     NSTEPS= No. of steps taken. (Initialize to zero.)
C     NBAD= No. of steps that had to be retaken. (initialize to zero.)
C     NMAX= Max. no. of allowed steps (NSTEPS + NBAD).
C
      
c      integer NN, NSTEPS, NBAD, NMAX
      real*8 A,E,H,HHH,ONE,T4,ABSH,HAIM,HSIGN,T1,TSAVE
      real*8 B,EPS,HGUS,T2,TWO,DUM,ERRMAX,HH,T3,ZR
c      real*8 DY1,DYDUB,DYSAVE,Y,Y2,Y3,YDUB,YSAVE,YSING
      real*8 :: Y(6),YSAVE(6),DYSAVE(6),YSING(6),YDUB(6),Y2(6)
      real*8 :: Y3(6),Y1(6),DYDUB(6),DY1(6)
      LOGICAL END
      DATA ZR/0.0/,ONE/1.0/,TWO/2.0/
      END=.FALSE.
      N=6
      IF(HGUS.EQ.ZR)HGUS=(B-A)*0.01
c CRK HERE 01/21/98 for portability
      IF (B-A.ge.0.0) THEN
        HSIGN = ONE
      ELSE
        HSIGN = -ONE
      ENDIF
c      HSIGN=SIGN(ONE,B-A)
c end CRK HERE 01/21/98
      HAIM=DABS(HGUS)
      DO 10 J=1,N
   10 YSAVE(J)=Y(J)
      TSAVE=A
   12 CALL DF(TSAVE,YSAVE,DYSAVE)
c CRK HERE 01/21/98 for portability
   11 CONTINUE
      IF (HSIGN.ge.0.0) THEN
        H = HAIM
      ELSE
        H = -HAIM
      ENDIF
c   11 H=SIGN(HAIM,HSIGN)
c end CRK HERE 01/21/98
      IF(HAIM.GT.DABS(B-TSAVE))THEN
         H=B-TSAVE
         END=.TRUE.
      ENDIF
C HERE TAKE A STEP GIVING YSING,YDUB,DYDUB. BEGIN.
      HH=0.5*H
      HHH=0.25*H
      T1=TSAVE+HHH
      T2=T1+HHH
      T3=T2+HHH
      T4=TSAVE+H
      DO 61 J=1,N
   61 Y2(J)=YSAVE(J)+HH*DYSAVE(J)
      CALL DF(T2,Y2,DY1)
      DO 62 J=1,N
      Y3(J)=YSAVE(J)+HH*DY1(J)
   62 Y2(J)=Y2(J)+Y3(J)*2.
      CALL DF(T2,Y3,DY1)
      DO 63 J=1,N
      Y3(J)=YSAVE(J)+H*DY1(J)
   63 Y2(J)=Y2(J)+Y3(J)
      CALL DF(T4,Y3,DY1)
      DO 64 J=1,N
      YSING(J)=(Y2(J)-YSAVE(J)+HH*DY1(J))/3.0
   64 Y2(J)=YSAVE(J)+HHH*DYSAVE(J)
      CALL DF(T1,Y2,DY1)
      DO 72 J=1,N
      Y3(J)=YSAVE(J)+HHH*DY1(J)
   72 Y2(J)=Y2(J)+Y3(J)*2.
      CALL DF(T1,Y3,DY1)
      DO 73 J=1,N
      Y3(J)=YSAVE(J)+HH*DY1(J)
   73 Y2(J)=Y2(J)+Y3(J)
      CALL DF(T2,Y3,DY1)
      DO 74 J=1,N
   74 YDUB(J)=(Y2(J)-YSAVE(J)+HHH*DY1(J))/3.0
      CALL DF(T2,YDUB,DYDUB)
      DO 81 J=1,N
   81 Y2(J)=YDUB(J)+HHH*DYDUB(J)
      CALL DF(T3,Y2,DY1)
      DO 82 J=1,N
      Y3(J)=YDUB(J)+HHH*DY1(J)
   82 Y2(J)=Y2(J)+Y3(J)*2.0
      CALL DF(T3,Y3,DY1)
      DO 83 J=1,N
      Y3(J)=YDUB(J)+HH*DY1(J)
   83 Y2(J)=Y2(J)+Y3(J)
      CALL DF(T4,Y3,DYDUB)
      ERRMAX=ZR
      ABSH=DABS(H)
      DO 84 J=1,N
      YDUB(J)=(Y2(J)-YDUB(J)+HHH*DYDUB(J))/3.0
C END OF TAKE A STEP BLOCK.
      DUM=DMAX1(ONE,DABS(YDUB(J)))
      E=DABS(YDUB(J)-YSING(J))/(15.*DUM)
      IF(E.GT.ERRMAX)ERRMAX=E
   84 CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.LE.ZR)THEN
         HAIM=ABSH*4.0
      ELSEIF(ERRMAX.LE.ONE)THEN
         HAIM=ABSH*ERRMAX**(-0.2)*.90
      ELSE
         HAIM=ABSH*ERRMAX**(-0.25)*.90
         NBAD=NBAD+1
         END=.FALSE.
         IF(NSTEPS+NBAD.GE.NMAX)GO TO 99
         GO TO 11
      ENDIF
C  STEP SUCCEEDED,.
      NSTEPS=NSTEPS+1
      TSAVE=TSAVE+H
      IF(END)THEN
         A=TSAVE
         HGUS=HAIM
         DO 91 J=1,N
   91    Y(J)=(16.*YDUB(J)-YSING(J))*6.666666666666666e-2
         RETURN
      ELSE
         IF(NSTEPS+NBAD.GE.NMAX)GO TO 99
         DO 31 J=1,N
   31    YSAVE(J)=(16.*YDUB(J)-YSING(J))*6.666666666666666e-2
         GO TO 12
      ENDIF
c CRK HERE for Fortran -> C interface
c   99 call crkerrC('NMAX EXCEEDED')
   99 PAUSE'NMAX EXCEEDED'
      END
