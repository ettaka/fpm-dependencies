C EXAMPLE PROBLEM
C
C The program below uses ZVODE to solve the following system of 2 ODEs:
C dw/dt = -i*w*w*z, dz/dt = i*z; w(0) = 1/2.1, z(0) = 1; t = 0 to 2*pi.
C Solution: w = 1/(z + 1.1), z = exp(it).  As z traces the unit circle,
C w traces a circle of radius 10/2.1 with center at 11/2.1.
C For convenience, Main passes RPAR = (imaginary unit i) to FEX and JEX.
C
      PROGRAM TESTZVODE
      EXTERNAL FEX, JEX
      DOUBLE COMPLEX Y(2), ZWORK(24), RPAR, WTRU, ERR
      DOUBLE PRECISION ABERR, AEMAX, ATOL, RTOL, RWORK(22), T, TOUT
      DIMENSION IWORK(32)
      NEQ = 2
      Y(1) = 1.0D0/2.1D0
      Y(2) = 1.0D0
      T = 0.0D0
      DTOUT = 0.1570796326794896D0
      TOUT = DTOUT
      ITOL = 1
      RTOL = 1.D-9
      ATOL = 1.D-8
      ITASK = 1
      ISTATE = 1
      IOPT = 0
      LZW = 24
      LRW = 22
      LIW = 32
      MF = 21
      RPAR = DCMPLX(0.0D0,1.0D0)
      AEMAX = 0.0D0
      WRITE(6,10)
  10  FORMAT('   t',11X,'w',26X,'z')
      DO 40 IOUT = 1,40
        CALL ZVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
     1             ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
        WTRU = 1.0D0/DCMPLX(COS(T) + 1.1D0, SIN(T))
        ERR = Y(1) - WTRU
        ABERR = ABS(DREAL(ERR)) + ABS(DIMAG(ERR))
        AEMAX = MAX(AEMAX,ABERR)
        WRITE(6,20) T, DREAL(Y(1)),DIMAG(Y(1)), DREAL(Y(2)),DIMAG(Y(2))
  20    FORMAT(F9.5,2X,2F12.7,3X,2F12.7)
        IF (ISTATE .LT. 0) THEN
          WRITE(6,30) ISTATE
  30      FORMAT(//'***** Error halt.  ISTATE =',I3)
          STOP
          ENDIF
  40    TOUT = TOUT + DTOUT
      WRITE(6,50) IWORK(11), IWORK(12), IWORK(13), IWORK(20),
     1            IWORK(21), IWORK(22), IWORK(23), AEMAX
  50  FORMAT(/' No. steps =',I4,'   No. f-s =',I5,
     1        '   No. J-s =',I4,'   No. LU-s =',I4/
     2        ' No. nonlinear iterations =',I4/
     3        ' No. nonlinear convergence failures =',I4/
     4        ' No. error test failures =',I4/
     5        ' Max. abs. error in w =',D10.2)
      STOP
      END
 
      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
      DOUBLE COMPLEX Y(NEQ), YDOT(NEQ), RPAR
      DOUBLE PRECISION T
      YDOT(1) = -RPAR*Y(1)*Y(1)*Y(2)
      YDOT(2) = RPAR*Y(2)
      RETURN
      END
 
      SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      DOUBLE COMPLEX Y(NEQ), PD(NRPD,NEQ), RPAR
      DOUBLE PRECISION T
      PD(1,1) = -2.0D0*RPAR*Y(1)*Y(2)
      PD(1,2) = -RPAR*Y(1)*Y(1)
      PD(2,2) = RPAR
      RETURN
      END
