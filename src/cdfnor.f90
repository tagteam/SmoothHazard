      SUBROUTINE cdfnor(which,p,q,x,mean,sd,status,bound)
!**********************************************************************
!
!      SUBROUTINE CDFNOR( WHICH, P, Q, X, MEAN, SD, STATUS, BOUND )
!               Cumulative Distribution Function
!               NORmal distribution
!
!
!                              Function
!
!
!     Calculates any one parameter of the normal
!     distribution given values for the others.
!
!
!                              Arguments
!
!
!     WHICH  --> Integer indicating  which of the  next  parameter
!     values is to be calculated using values  of the others.
!     Legal range: 1..4
!               iwhich = 1 : Calculate P and Q from X,MEAN and SD
!               iwhich = 2 : Calculate X from P,Q,MEAN and SD
!               iwhich = 3 : Calculate MEAN from P,Q,X and SD
!               iwhich = 4 : Calculate SD from P,Q,X and MEAN
!                    INTEGER WHICH
!
!     P <--> The integral from -infinity to X of the normal density.
!            Input range: (0,1].
!                    DOUBLE PRECISION P
!
!     Q <--> 1-P.
!            Input range: (0, 1].
!            P + Q = 1.0.
!                    DOUBLE PRECISION Q
!
!     X < --> Upper limit of integration of the normal-density.
!             Input range: ( -infinity, +infinity)
!                    DOUBLE PRECISION X
!
!     MEAN <--> The mean of the normal density.
!               Input range: (-infinity, +infinity)
!                    DOUBLE PRECISION MEAN
!
!     SD <--> Standard Deviation of the normal density.
!             Input range: (0, +infinity).
!                    DOUBLE PRECISION SD
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q .ne. 1
!                    INTEGER STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
!
!                              Method
!
!
!
!
!     A slightly modified version of ANORM from
!
!     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
!     Package of Special Function Routines and Test Drivers"
!     acm Transactions on Mathematical Software. 19, 22-32.
!
!     is used to calulate the  cumulative standard normal distribution.
!
!     The rational functions from pages  90-95  of Kennedy and Gentle,
!     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
!     starting values to Newton's Iterations which compute the inverse
!     standard normal.  Therefore no  searches  are necessary for  any
!     parameter.
!
!     For X < -15, the asymptotic expansion for the normal is used  as
!     the starting value in finding the inverse standard normal.
!     This is formula 26.2.12 of Abramowitz and Stegun.
!
!
!                              Note
!
!
!      The normal density is proportional to
!      exp( - 0.5 * (( X - MEAN)/SD)**2)
!
!
!**********************************************************************
!     .. Parameters ..
!     ..
      implicit none	
!     .. Scalar Arguments ..
      DOUBLE PRECISION bound,mean,p,q,sd,x
      INTEGER status,which
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION z,pq
!     ..
!     .. External Functions ..

      DOUBLE PRECISION dinvnr,spmpar
      EXTERNAL dinvnr,spmpar
!     ..
!     .. External Subroutines ..
      EXTERNAL cumnor
!     ..
!     .. Executable Statements ..
!
!     Check arguments
!
      status = 0
      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 4.0D0
   20 status = -1
      RETURN

   30 IF (which.EQ.1) GO TO 70
!
!     P
!
      IF (.NOT. ((p.LE.0.0D0).OR. (p.GT.1.0D0))) GO TO 60
      IF (.NOT. (p.LE.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = 1.0D0
   50 status = -2
      RETURN

   60 CONTINUE
   70 IF (which.EQ.1) GO TO 110
!
!     Q
!
      IF (.NOT. ((q.LE.0.0D0).OR. (q.GT.1.0D0))) GO TO 100
      IF (.NOT. (q.LE.0.0D0)) GO TO 80
      bound = 0.0D0
      GO TO 90

   80 bound = 1.0D0
   90 status = -3
      RETURN

  100 CONTINUE
  110 IF (which.EQ.1) GO TO 150
!
!     P + Q
!
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.(3.0D0*spmpar(1))))	GO TO 140
      IF (.NOT. (pq.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      GO TO 130

  120 bound = 1.0D0
  130 status = 3
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 170
!
!     SD
!
      IF (.NOT. (sd.LE.0.0D0)) GO TO 160
      bound = 0.0D0
      status = -6
      RETURN

  160 CONTINUE
!
!     Calculate ANSWERS
!
  170 IF ((1).EQ. (which)) THEN
!
!     Computing P
!
          z = (x-mean)/sd
          CALL cumnor(z,p,q)

      ELSE IF ((2).EQ. (which)) THEN
!
!     Computing X
!
          z = dinvnr(p,q)
          x = sd*z + mean

      ELSE IF ((3).EQ. (which)) THEN
!
!     Computing the MEAN
!
          z = dinvnr(p,q)
          mean = x - sd*z

      ELSE IF ((4).EQ. (which)) THEN
!
!     Computing SD
!
          z = dinvnr(p,q)
          sd = (x-mean)/z
      END IF

      RETURN

      END SUBROUTINE



      DOUBLE PRECISION FUNCTION dinvnr(p,q)
!**********************************************************************
!
!     DOUBLE PRECISION FUNCTION DINVNR(P,Q)
!     Double precision NoRmal distribution INVerse
!
!
!                              Function
!
!
!     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
!     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
!
!
!                              Arguments
!
!
!     P --> The probability whose normal deviate is sought.
!                    P is DOUBLE PRECISION
!
!     Q --> 1-P
!                    P is DOUBLE PRECISION
!
!
!                              Method
!
!
!     The  rational   function   on  page 95    of Kennedy  and  Gentle,
!     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
!     value for the Newton method of finding roots.
!
!
!                              Note
!
!
!     If P or Q .lt. machine EPS returns +/- DINVNR(EPS)
!
!**********************************************************************
      implicit none
!     .. Parameters ..
      INTEGER maxit
      PARAMETER (maxit=100)
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0D-13)
      DOUBLE PRECISION r2pi
      PARAMETER (r2pi=0.3989422804014326D0)
      DOUBLE PRECISION nhalf
      PARAMETER (nhalf=-0.5D0)
!     ..
!     .. Scalar Arguments ..
      DOUBLE PRECISION p,q
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION strtx,xcur,cum,ccum,pp,dx
      INTEGER i
      LOGICAL qporq
!     ..
!     .. External Functions ..
      DOUBLE PRECISION stvaln
      EXTERNAL stvaln
!     ..
!     .. External Subroutines ..
      EXTERNAL cumnor
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION dennor,x

      dennor(x) = r2pi*exp(nhalf*x*x)
!     ..
!     .. Executable Statements ..
!
!     FIND MINIMUM OF P AND Q
!
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 10
      pp = p
      GO TO 20

   10 pp = q
!
!     INITIALIZATION STEP
!
   20 strtx = stvaln(pp)
      xcur = strtx
!
!     NEWTON INTERATIONS
!
      DO 30,i = 1,maxit
          CALL cumnor(xcur,cum,ccum)
          dx = (cum-pp)/dennor(xcur)
          xcur = xcur - dx
          IF (abs(dx/xcur).LT.eps) GO TO 40
   30 CONTINUE
      dinvnr = strtx
!
!     IF WE GET HERE, NEWTON HAS FAILED
!
      IF (.NOT.qporq) dinvnr = -dinvnr
      RETURN
!
!     IF WE GET HERE, NEWTON HAS SUCCEDED
!
   40 dinvnr = xcur
      IF (.NOT.qporq) dinvnr = -dinvnr
      RETURN

      END



      INTEGER FUNCTION ipmpar(i)
!-----------------------------------------------------------------------
!
!     IPMPAR PROVIDES THE INTEGER MA!HINE CONSTANTS FOR THE COMPUTER
!     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
!
!  INTEGERS.
!
!     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
!
!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
!
!     IPMPAR(1) = A, THE BASE.
!
!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
!
!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS.
!
!     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
!
!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
!
!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
!
!     IPMPAR(4) = B, THE BASE.
!
!  SINGLE-PRECISION
!
!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
!
!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
!
!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION
!
!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
!
!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
!
!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
!
!-----------------------------------------------------------------------
!
!     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE
!     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM
!     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN
!     COLUMN 1.)
!
!-----------------------------------------------------------------------
!
!     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
!     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
!     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
!     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
!
!-----------------------------------------------------------------------
      implicit none
!     .. Scalar Arguments ..
      INTEGER i
!     ..
!     .. Local Arrays ..
      INTEGER imach(10)
!     ..
!     .. Data statements ..
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!     DATA IMACH( 1) /   2 /
!     DATA IMACH( 2) /  31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /  16 /
!     DATA IMACH( 5) /   6 /
!     DATA IMACH( 6) / -64 /
!     DATA IMACH( 7) /  63 /
!     DATA IMACH( 8) /  14 /
!     DATA IMACH( 9) / -64 /
!     DATA IMACH(10) /  63 /
!
!     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
!     PC 7300, AND AT&T 6300.
!
!     DATA IMACH( 1) /     2 /
!     DATA IMACH( 2) /    31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /     2 /
!     DATA IMACH( 5) /    24 /
!     DATA IMACH( 6) /  -125 /
!     DATA IMACH( 7) /   128 /
!     DATA IMACH( 8) /    53 /
!     DATA IMACH( 9) / -1021 /
!     DATA IMACH(10) /  1024 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   33 /
!     DATA IMACH( 3) / 8589934591 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   24 /
!     DATA IMACH( 6) / -256 /
!     DATA IMACH( 7) /  255 /
!     DATA IMACH( 8) /   60 /
!     DATA IMACH( 9) / -256 /
!     DATA IMACH(10) /  255 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   39 /
!     DATA IMACH( 3) / 549755813887 /
!     DATA IMACH( 4) /    8 /
!     DATA IMACH( 5) /   13 /
!     DATA IMACH( 6) /  -50 /
!     DATA IMACH( 7) /   76 /
!     DATA IMACH( 8) /   26 /
!     DATA IMACH( 9) /  -50 /
!     DATA IMACH(10) /   76 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
!
!     DATA IMACH( 1) /      2 /
!     DATA IMACH( 2) /     39 /
!     DATA IMACH( 3) / 549755813887 /
!     DATA IMACH( 4) /      8 /
!     DATA IMACH( 5) /     13 /
!     DATA IMACH( 6) /    -50 /
!     DATA IMACH( 7) /     76 /
!     DATA IMACH( 8) /     26 /
!     DATA IMACH( 9) / -32754 /
!     DATA IMACH(10) /  32780 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS OPERATING SYSTEM).
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   48 /
!     DATA IMACH( 3) / 281474976710655 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   48 /
!     DATA IMACH( 6) / -974 /
!     DATA IMACH( 7) / 1070 /
!     DATA IMACH( 8) /   95 /
!     DATA IMACH( 9) / -926 /
!     DATA IMACH(10) / 1070 /
!
!     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS/VE OPERATING SYSTEM).
!
!     DATA IMACH( 1) /     2 /
!     DATA IMACH( 2) /    63 /
!     DATA IMACH( 3) / 9223372036854775807 /
!     DATA IMACH( 4) /     2 /
!     DATA IMACH( 5) /    48 /
!     DATA IMACH( 6) / -4096 /
!     DATA IMACH( 7) /  4095 /
!     DATA IMACH( 8) /    96 /
!     DATA IMACH( 9) / -4096 /
!     DATA IMACH(10) /  4095 /
!
!     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!
!     DATA IMACH( 1) /     2 /
!     DATA IMACH( 2) /    63 /
!     DATA IMACH( 3) / 9223372036854775807 /
!     DATA IMACH( 4) /     2 /
!     DATA IMACH( 5) /    47 /
!     DATA IMACH( 6) / -8189 /
!     DATA IMACH( 7) /  8190 /
!     DATA IMACH( 8) /    94 /
!     DATA IMACH( 9) / -8099 /
!     DATA IMACH(10) /  8190 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   15 /
!     DATA IMACH( 3) / 32767 /
!     DATA IMACH( 4) /   16 /
!     DATA IMACH( 5) /    6 /
!     DATA IMACH( 6) /  -64 /
!     DATA IMACH( 7) /   63 /
!     DATA IMACH( 8) /   14 /
!     DATA IMACH( 9) /  -64 /
!     DATA IMACH(10) /   63 /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   23 /
!     DATA IMACH( 3) / 8388607 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   23 /
!     DATA IMACH( 6) / -127 /
!     DATA IMACH( 7) /  127 /
!     DATA IMA!H( 8) /   38 /
!     DATA IMACH( 9) / -127 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
!     AND DPS 8/70 SERIES.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   35 /
!     DATA IMACH( 3) / 34359738367 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   27 /
!     DATA IMACH( 6) / -127 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   63 /
!     DATA IMACH( 9) / -127 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   15 /
!     DATA IMACH( 3) / 32767 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   23 /
!     DATA IMACH( 6) / -128 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   39 /
!     DATA IMACH( 9) / -128 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   15 /
!     DATA IMACH( 3) / 32767 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   23 /
!     DATA IMACH( 6) / -128 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   55 /
!     DATA IMACH( 9) / -128 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE HP 9000.
!
!     DATA IMACH( 1) /     2 /
!     DATA IMACH( 2) /    31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /     2 /
!     DATA IMACH( 5) /    24 /
!     DATA IMACH( 6) /  -126 /
!     DATA IMACH( 7) /   128 /
!     DATA IMACH( 8) /    53 /
!     DATA IMACH( 9) / -1021 /
!     DATA IMACH(10) /  1024 /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
!     5/7/9 AND THE SEL SYSTEMS 85/86.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /   16 /
!     DATA IMACH( 5) /    6 /
!     DATA IMACH( 6) /  -64 /
!     DATA IMACH( 7) /   63 /
!     DATA IMACH( 8) /   14 /
!     DATA IMACH( 9) /  -64 /
!     DATA IMACH(10) /   63 /
!
!     MACHINE CONSTANTS FOR THE IBM PC.
!
!      DATA imach(1)/2/
!      DATA imach(2)/31/
!      DATA imach(3)/2147483647/
!      DATA imach(4)/2/
!      DATA imach(5)/24/
!      DATA imach(6)/-125/
!      DATA imach(7)/128/
!      DATA imach(8)/53/
!      DATA imach(9)/-1021/
!      DATA imach(10)/1024/
!
!     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
!     MA!FORTRAN II.
!
!     DATA IMACH( 1) /     2 /
!     DATA IMACH( 2) /    31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /     2 /
!     DATA IMACH( 5) /    24 /
!     DATA IMACH( 6) /  -125 /
!     DATA IMACH( 7) /   128 /
!     DATA IMACH( 8) /    53 /
!     DATA IMACH( 9) / -1021 /
!     DATA IMACH(10) /  1024 /
!
!     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   24 /
!     DATA IMACH( 6) / -127 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   56 /
!     DATA IMACH( 9) / -127 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   35 /
!     DATA IMACH( 3) / 34359738367 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   27 /
!     DATA IMACH( 6) / -128 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   54 /
!     DATA IMACH( 9) / -101 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   35 /
!     DATA IMACH( 3) / 34359738367 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   27 /
!     DATA IMACH( 6) / -128 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   62 /
!     DATA IMACH( 9) / -128 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGER ARITHMETI!.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   24 /
!     DATA IMACH( 6) / -127 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   56 /
!     DATA IMACH( 9) / -127 /
!     DATA IMACH(10) /  127 /
!
!     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
!
!     DATA IMACH( 1) /     2 /
!     DATA IMACH( 2) /    31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /     2 /
!     DATA IMACH( 5) /    24 /
!     DATA IMACH( 6) /  -125 /
!     DATA IMACH( 7) /   128 /
!     DATA IMACH( 8) /    53 /
!     DATA IMACH( 9) / -1021 /
!     DATA IMACH(10) /  1024 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
!     SERIES (MIPS R3000 PROCESSOR).
!
!     DATA IMACH( 1) /     2 /
!     DATA IMACH( 2) /    31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /     2 /
!     DATA IMACH( 5) /    24 /
!     DATA IMACH( 6) /  -125 /
!     DATA IMACH( 7) /   128 /
!     DATA IMACH( 8) /    53 /
!     DATA IMACH( 9) / -1021 /
!     DATA IMACH(10) /  1024 /
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!
      DATA IMACH( 1) /     2 /
      DATA IMACH( 2) /    31 /
      DATA IMACH( 3) / 2147483647 /
      DATA IMACH( 4) /     2 /
      DATA IMACH( 5) /    24 /
      DATA IMACH( 6) /  -125 /
      DATA IMACH( 7) /   128 /
      DATA IMACH( 8) /    53 /
      DATA IMACH( 9) / -1021 /
      DATA IMACH(10) /  1024 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   35 /
!     DATA IMACH( 3) / 34359738367 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   27 /
!     DATA IMACH( 6) / -128 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   60 /
!     DATA IMACH( 9) /-1024 /
!     DATA IMACH(10) / 1023 /
!
!     MACHINE CONSTANTS FOR THE VAX 11/780.
!
!     DATA IMACH( 1) /    2 /
!     DATA IMACH( 2) /   31 /
!     DATA IMACH( 3) / 2147483647 /
!     DATA IMACH( 4) /    2 /
!     DATA IMACH( 5) /   24 /
!     DATA IMACH( 6) / -127 /
!     DATA IMACH( 7) /  127 /
!     DATA IMACH( 8) /   56 /
!     DATA IMACH( 9) / -127 /
!     DATA IMACH(10) /  127 /
!
      ipmpar = imach(i)
      RETURN

      END


      DOUBLE PRECISION FUNCTION spmpar(i)
!-----------------------------------------------------------------------
!
!     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
!     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
!     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
!     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
!     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
!
!        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
!
!        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
!
!        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
!
!-----------------------------------------------------------------------
!     WRITTEN BY
!        ALFRED H. MORRIS, JR.
!        NAVAL SURFACE WARFARE CENTER
!        DAHLGREN VIRGINIA
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
!     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
!     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
!-----------------------------------------------------------------------
      implicit none
!     .. Scalar Arguments ..
      INTEGER i
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION b,binv,bm1,one,w,z
      INTEGER emax,emin,ibeta,m
!     ..
!     .. External Functions ..
      INTEGER ipmpar
      EXTERNAL ipmpar
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC dble
!     ..
!     .. Executable Statements ..
!
      IF (i.GT.1) GO TO 10
      b = ipmpar(4)
      m = ipmpar(8)
      spmpar = b** (1-m)
      RETURN
!
   10 IF (i.GT.2) GO TO 20
      b = ipmpar(4)
      emin = ipmpar(9)
      one = dble(1)
      binv = one/b
      w = b** (emin+2)
      spmpar = ((w*binv)*binv)*binv
      RETURN
!
   20 ibeta = ipmpar(4)
      m = ipmpar(8)
      emax = ipmpar(10)
!
      b = ibeta
      bm1 = ibeta - 1
      one = dble(1)
      z = b** (m-1)
      w = ((z-one)*b+bm1)/ (b*z)
!
      z = b** (emax-2)
      spmpar = ((w*z)*b)*b
      RETURN

      END




      SUBROUTINE cumnor(arg,result,ccum)
!**********************************************************************
!
!     SUBROUINE CUMNOR(X,RESULT,CCUM)
!
!
!                              Function
!
!
!     Computes the cumulative  of    the  normal   distribution,   i.e.,
!     the integral from -infinity to x of
!          (1/sqrt(2*pi)) exp(-u*u/2) du
!
!     X --> Upper limit of integration.
!                                        X is DOUBLE PRECISION
!
!     RESULT <-- Cumulative normal distribution.
!                                        RESULT is DOUBLE PRECISION
!
!     CCUM <-- Compliment of Cumulative normal distribution.
!                                        CCUM is DOUBLE PRECISION
!
!
!     Renaming of function ANORM from:
!
!     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
!     Package of Special Function Routines and Test Drivers"
!     acm Transactions on Mathematical Software. 19, 22-32.
!
!     with slight modifications to return ccum and to deal with
!     machine constants.
!
!**********************************************************************
!
!
! Original Comments:
!------------------------------------------------------------------
!
! This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!   The main computation evaluates near-minimax approximations
!   derived from those in "Rational Chebyshev approximations for
!   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
!   This transportable program uses rational functions that
!   theoretically approximate the normal distribution function to
!   at least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants.
!
!   MIN   = smallest machine representable number.
!
!   EPS   = argument below which anorm(x) may be represented by
!           0.5  and above which  x*x  will not underflow.
!           A conservative value is the largest machine number X
!           such that   1.0 + X = 1.0   to machine precision.
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns  ANORM = 0     for  ARG .LE. XLOW.
!
!
! Intrinsic functions required are:
!
!     ABS, AINT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics and Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: March 15, 1992
!
!------------------------------------------------------------------
      implicit none
      INTEGER i
      DOUBLE PRECISION a,arg,b,c,d,del,eps,half,p,one,q,result,sixten,&
                       temp,sqrpi,thrsh,root32,x,xden,xnum,y,xsq,zero,&
                       min,ccum
      DIMENSION a(5),b(4),c(9),d(8),p(6),q(5)
!------------------------------------------------------------------
!  External Function
!------------------------------------------------------------------
      DOUBLE PRECISION spmpar
      EXTERNAL spmpar
!------------------------------------------------------------------
!  Mathematical constants
!
!  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
!  THRSH is the argument for which anorm = 0.75.
!------------------------------------------------------------------
      DATA one,half,zero,sixten/1.0D0,0.5D0,0.0D0,1.60D0/,&
           sqrpi/3.9894228040143267794D-1/,thrsh/0.66291D0/,&
           root32/5.656854248D0/
!------------------------------------------------------------------
!  Coefficients for approximation in first interval
!------------------------------------------------------------------
      DATA a/2.2352520354606839287D00,1.6102823106855587881D02,&
           1.0676894854603709582D03,1.8154981253343561249D04,&
           6.5682337918207449113D-2/
      DATA b/4.7202581904688241870D01,9.7609855173777669322D02,&
           1.0260932208618978205D04,4.5507789335026729956D04/
!------------------------------------------------------------------
!  Coefficients for approximation in second interval
!------------------------------------------------------------------
      DATA c/3.9894151208813466764D-1,8.8831497943883759412D00,&
           9.3506656132177855979D01,5.9727027639480026226D02,&
           2.4945375852903726711D03,6.8481904505362823326D03,&
           1.1602651437647350124D04,9.8427148383839780218D03,&
           1.0765576773720192317D-8/
      DATA d/2.2266688044328115691D01,2.3538790178262499861D02,&
           1.5193775994075548050D03,6.4855582982667607550D03,&
           1.8615571640885098091D04,3.4900952721145977266D04,&
           3.8912003286093271411D04,1.9685429676859990727D04/
!------------------------------------------------------------------
!  Coefficients for approximation in third interval
!------------------------------------------------------------------
      DATA p/2.1589853405795699D-1,1.274011611602473639D-1,&
           2.2235277870649807D-2,1.421619193227893466D-3,&
           2.9112874951168792D-5,2.307344176494017303D-2/
      DATA q/1.28426009614491121D00,4.68238212480865118D-1,&
           6.59881378689285515D-2,3.78239633202758244D-3,&
           7.29751555083966205D-5/
!------------------------------------------------------------------
!  Machine dependent constants
!------------------------------------------------------------------
      eps = spmpar(1)*0.5D0
      min = spmpar(2)
!------------------------------------------------------------------
      x = arg
      y = abs(x)
      IF (y.LE.thrsh) THEN
!------------------------------------------------------------------
!  Evaluate  anorm  for  |X| <= 0.66291
!------------------------------------------------------------------
          xsq = zero
          IF (y.GT.eps) xsq = x*x
          xnum = a(5)*xsq
          xden = xsq
          DO 10 i = 1,3
              xnum = (xnum+a(i))*xsq
              xden = (xden+b(i))*xsq
   10     CONTINUE
          result = x* (xnum+a(4))/ (xden+b(4))
          temp = result
          result = half + temp
          ccum = half - temp
!------------------------------------------------------------------
!  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
!------------------------------------------------------------------
      ELSE IF (y.LE.root32) THEN
          xnum = c(9)*y
          xden = y
          DO 20 i = 1,7
              xnum = (xnum+c(i))*y
              xden = (xden+d(i))*y
   20     CONTINUE
          result = (xnum+c(8))/ (xden+d(8))
          xsq = aint(y*sixten)/sixten
          del = (y-xsq)* (y+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF
!------------------------------------------------------------------
!  Evaluate  anorm  for |X| > sqrt(32)
!------------------------------------------------------------------
      ELSE
          result = zero
          xsq = one/ (x*x)
          xnum = p(6)*xsq
          xden = xsq
          DO 30 i = 1,4
              xnum = (xnum+p(i))*xsq
              xden = (xden+q(i))*xsq
   30     CONTINUE
          result = xsq* (xnum+p(5))/ (xden+q(5))
          result = (sqrpi-result)/y
          xsq = aint(x*sixten)/sixten
          del = (x-xsq)* (x+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF

      END IF

      IF (result.LT.min) result = 0.0D0
      IF (ccum.LT.min) ccum = 0.0D0
!------------------------------------------------------------------
!  Fix up for negative argument, erf, etc.
!------------------------------------------------------------------
!----------Last card of ANORM ----------
      END




      DOUBLE PRECISION FUNCTION stvaln(p)
!
!**********************************************************************
!
!     DOUBLE PRECISION FUNCTION STVALN(P)
!                    STarting VALue for Neton-Raphon
!                calculation of Normal distribution Inverse
!
!
!                              Function
!
!
!     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
!     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
!
!
!                              Arguments
!
!
!     P --> The probability whose normal deviate is sought.
!                    P is DOUBLE PRECISION
!
!
!                              Method
!
!
!     The  rational   function   on  page 95    of Kennedy  and  Gentle,
!     Statistical Computing, Marcel Dekker, NY , 1980.
!
!**********************************************************************
!
      implicit none
!     .. Scalar Arguments ..
      DOUBLE PRECISION p
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION sign,y,z
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION xden(5),xnum(5)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION devlpl
      EXTERNAL devlpl
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC dble,log,sqrt
!     ..
!     .. Data statements ..
      DATA xnum/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,&
           -0.204231210245D-1,-0.453642210148D-4/
      DATA xden/0.993484626060D-1,0.588581570495D0,0.531103462366D0,&
           0.103537752850D0,0.38560700634D-2/
!     ..
!     .. Executable Statements ..
      IF (.NOT. (p.LE.0.5D0)) GO TO 10
      sign = -1.0D0
      z = p
      GO TO 20

   10 sign = 1.0D0
      z = 1.0D0 - p
   20 y = sqrt(-2.0D0*log(z))
      stvaln = y + devlpl(xnum,5,y)/devlpl(xden,5,y)
      stvaln = sign*stvaln
      RETURN

      END



      DOUBLE PRECISION FUNCTION devlpl(a,n,x)
!**********************************************************************
!
!     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)
!              Double precision EVALuate a PoLynomial at X
!
!
!                              Function
!
!
!     returns
!          A(1) + A(2)*X + ... + A(N)*X**(N-1)
!
!
!                              Arguments
!
!
!     A --> Array of coefficients of the polynomial.
!                                        A is DOUBLE PRECISION(N)
!
!     N --> Length of A, also degree of polynomial - 1.
!                                        N is INTEGER
!
!     X --> Point at which the polynomial is to be evaluated.
!                                        X is DOUBLE PRECISION
!
!**********************************************************************
      implicit none
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER n
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION a(n)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION term
      INTEGER i
!     ..
!     .. Executable Statements ..
      term = a(n)
      DO 10,i = n - 1,1,-1
          term = a(i) + term*x
   10 CONTINUE
      devlpl = term
      RETURN

      END
