      REAL FUNCTION PLKINT( WNUMLO, WNUMHI, T )

c        Computes Planck function integrated between two wavenumbers
c
c  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval
c
c           WNUMHI : Upper wavenumber
c
c           T      : Temperature (K)
c
c  OUTPUT : PLKINT : Integrated Planck function ( Watts/sq m )
c                      = Integral (WNUMLO to WNUMHI) of
c                        2h c**2  nu**3 / ( EXP(hc nu/kT) - 1)
c                        (where h=Plancks constant, c=speed of
c                         light, nu=wavenumber, T=temperature,
c                         and k = Boltzmann constant)
c
c  Reference : Specifications of the Physical World: New Value
c                 of the Fundamental Constants, Dimensions/N.B.S.,
c                 Jan. 1974
c
c  Method :  For WNUMLO close to WNUMHI, a Simpson-rule quadrature
c            is done to avoid ill-conditioning; otherwise
c
c            (1)  For WNUMLO or WNUMHI small,
c                 integral(0 to WNUMLO/HI) is calculated by expanding
c                 the integrand in a power series and integrating
c                 term by term;
c
c            (2)  Otherwise, integral(WNUMLO/HI to INFINITY) is
c                 calculated by expanding the denominator of the
c                 integrand in powers of the exponential and
c                 integrating term by term.
c
c  Accuracy :  At least 6 significant digits, assuming the
c              physical constants are infinitely accurate
c
c  ERRORS WHICH ARE NOT TRAPPED:
c
c      * power or exponential series may underflow, giving no
c        significant digits.  This may or may not be of concern,
c        depending on the application.
c
c      * Simpson-rule special case is skipped when denominator of
c        integrand will cause overflow.  In that case the normal
c        procedure is used, which may be inaccurate if the
c        wavenumber limits (WNUMLO, WNUMHI) are close together.
c
c  LOCAL VARIABLES
c
c        A1,2,... :  Power series coefficients
c        C2       :  h * c / k, in units cm*K (h = Plancks constant,
c                      c = speed of light, k = Boltzmann constant)
c        D(I)     :  Exponential series expansion of integral of
c                       Planck function from WNUMLO (i=1) or WNUMHI
c                       (i=2) to infinity
c        EPSIL    :  Smallest number such that 1+EPSIL .GT. 1 on
c                       computer
c        EX       :  EXP( - V(I) )
c        EXM      :  EX**M
c        MMAX     :  No. of terms to take in exponential series
c        MV       :  Multiples of V(I)
c        P(I)     :  Power series expansion of integral of
c                       Planck function from zero to WNUMLO (I=1) or
c                       WNUMHI (I=2)
c        PI       :  3.14159...
c        SIGMA    :  Stefan-Boltzmann constant (W/m**2/K**4)
c        SIGDPI   :  SIGMA / PI
c        SMALLV   :  Number of times the power series is used (0,1,2)
c        V(I)     :  C2 * (WNUMLO(I=1) or WNUMHI(I=2)) / temperature
c        VCUT     :  Power-series cutoff point
c        VCP      :  Exponential series cutoff points
c        VMAX     :  Largest allowable argument of EXP function
c
c   Called by- DISORT
c   Calls- R1MACH, ERRMSG
c ----------------------------------------------------------------------

c     .. Parameters ..

      REAL      A1, A2, A3, A4, A5, A6
      PARAMETER ( A1 = 1. / 3., A2 = -1. / 8., A3 = 1. / 60.,
     &          A4 = -1. / 5040., A5 = 1. / 272160.,
     &          A6 = -1. / 13305600. )
c     ..
c     .. Scalar Arguments ..

      REAL      WNUMHI, WNUMLO, T
c     ..
c     .. Local Scalars ..

      INTEGER   I, K, M, MMAX, N, SMALLV
      REAL      C2, CONC, DEL, EPSIL, EX, EXM, HH, MV, OLDVAL, PI,
     &          SIGDPI, SIGMA, VAL, VAL0, VCUT, VMAX, VSQ, X
c     ..
c     .. Local Arrays ..

      REAL      D( 2 ), P( 2 ), V( 2 ), VCP( 7 )
c     ..
c     .. External Functions ..

      REAL      R1MACH
      EXTERNAL  R1MACH
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG2
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, ASIN, LOG, MOD
c     ..
c     .. Statement Functions ..

      REAL      PLKF
c     ..
      SAVE      PI, CONC, VMAX, EPSIL, SIGDPI

      DATA      C2 / 1.438786 / , SIGMA / 5.67032E-8 / , VCUT / 1.5 / ,
     &          VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      DATA      PI / 0.0 /

c     .. Statement Function definitions ..

      PLKF( X ) = X**3 / ( EXP( X ) - 1 )
c     ..

      IF( PI .EQ. 0.0 ) THEN

         PI     = 2.*ASIN( 1.0 )
         VMAX   = LOG( R1MACH( 2 ) )
         EPSIL  = R1MACH( 4 )
         SIGDPI = SIGMA / PI
         CONC   = 15. / PI**4

      END IF


      IF( T.LT.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LT.0. )
     &    CALL ERRMSG2('PLKINT--temperature or wavenums. wrong',.TRUE.)


      IF( T .LT. 1.E-4 ) THEN

         PLKINT = 0.0
         RETURN

      END IF


      V( 1 ) = C2*WNUMLO / T
      V( 2 ) = C2*WNUMHI / T

      IF( V( 1 ).GT.EPSIL .AND. V( 2 ).LT.VMAX .AND.
     &    ( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.E-2 ) THEN

c                          ** Wavenumbers are very close.  Get integral
c                          ** by iterating Simpson rule to convergence.

         HH     = V( 2 ) - V( 1 )
         OLDVAL = 0.0
         VAL0   = PLKF( V( 1 ) ) + PLKF( V( 2 ) )

         DO 20 N = 1, 10

            DEL  = HH / ( 2*N )
            VAL  = VAL0

            DO 10 K = 1, 2*N - 1
               VAL  = VAL + 2*( 1 + MOD( K,2 ) )*
     &                      PLKF( V( 1 ) + K*DEL )
   10       CONTINUE

            VAL  = DEL / 3.*VAL
            IF( ABS( ( VAL - OLDVAL ) / VAL ).LE.1.E-6 ) GO TO  30
            OLDVAL = VAL

   20    CONTINUE

         CALL ERRMSG2( 'PLKINT--Simpson rule didnt converge',.FALSE.)

   30    CONTINUE

         PLKINT = SIGDPI * T**4 * CONC * VAL

         RETURN

      END IF

c                          *** General case ***
      SMALLV = 0

      DO 60 I = 1, 2

         IF( V( I ).LT.VCUT ) THEN
c                                   ** Use power series
            SMALLV = SMALLV + 1
            VSQ    = V( I )**2
            P( I ) = CONC*VSQ*V( I )*( A1 +
     &               V( I )*( A2 + V( I )*( A3 + VSQ*( A4 + VSQ*( A5 +
     &               VSQ*A6 ) ) ) ) )

         ELSE
c                      ** Use exponential series
            MMAX  = 0
c                                ** Find upper limit of series
   40       CONTINUE
            MMAX  = MMAX + 1

            IF( V(I) .LT. VCP( MMAX ) ) GO TO  40

            EX     = EXP( - V(I) )
            EXM    = 1.0
            D( I ) = 0.0

            DO 50 M = 1, MMAX
               MV     = M*V( I )
               EXM    = EX*EXM
               D( I ) = D( I ) + EXM*( 6.+ MV*( 6.+ MV*( 3.+ MV ) ) )
     &                  / M**4
   50       CONTINUE

            D( I ) = CONC*D( I )

         END IF

   60 CONTINUE

c                              ** Handle ill-conditioning
      IF( SMALLV.EQ.2 ) THEN
c                                    ** WNUMLO and WNUMHI both small
         PLKINT = P( 2 ) - P( 1 )

      ELSE IF( SMALLV.EQ.1 ) THEN
c                                    ** WNUMLO small, WNUMHI large
         PLKINT = 1.- P( 1 ) - D( 2 )

      ELSE
c                                    ** WNUMLO and WNUMHI both large
         PLKINT = D( 1 ) - D( 2 )

      END IF

      PLKINT = SIGDPI * T**4 * PLKINT

      RETURN
      END


C  From http://www.netlib.org/port/ May 29, 2007
      REAL FUNCTION R1MACH(I)
      INTEGER I
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  R1MACH(5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C     needs to be (2) for AUTODOUBLE, HARRIS SLASH 6, ...
      INTEGER SC
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      REAL RMACH(5)
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
      INTEGER J, K, L, T3E(3)
      DATA T3E(1) / 9777664 /
      DATA T3E(2) / 5323660 /
      DATA T3E(3) / 46980 /
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES,
C  INCLUDING AUTO-DOUBLE COMPILERS.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /, SC/987/
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /, SC/987/
C
      IF (SC .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH(1) = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      SMALL(1) .EQ. 1117925532
     *          .AND. SMALL(2) .EQ. -448790528) THEN
*              *** IEEE BIG ENDIAN ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2146435071
               LARGE(2) = -1
               RIGHT(1) = 1017118720
               RIGHT(2) = 0
               DIVER(1) = 1018167296
               DIVER(2) = 0
               LOG10(1) = 1070810131
               LOG10(2) = 1352628735
            ELSE IF ( SMALL(2) .EQ. 1117925532
     *          .AND. SMALL(1) .EQ. -448790528) THEN
*              *** IEEE LITTLE ENDIAN ***
               SMALL(2) = 1048576
               SMALL(1) = 0
               LARGE(2) = 2146435071
               LARGE(1) = -1
               RIGHT(2) = 1017118720
               RIGHT(1) = 0
               DIVER(2) = 1018167296
               DIVER(1) = 0
               LOG10(2) = 1070810131
               LOG10(1) = 1352628735
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*              *** VAX WITH D_FLOATING ***
               SMALL(1) = 128
               SMALL(2) = 0
               LARGE(1) = -32769
               LARGE(2) = -1
               RIGHT(1) = 9344
               RIGHT(2) = 0
               DIVER(1) = 9472
               DIVER(2) = 0
               LOG10(1) = 546979738
               LOG10(2) = -805796613
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2147483647
               LARGE(2) = -1
               RIGHT(1) = 856686592
               RIGHT(2) = 0
               DIVER(1) = 873463808
               DIVER(2) = 0
               LOG10(1) = 1091781651
               LOG10(2) = 1352628735
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
         ELSE
            RMACH(1) = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*              *** IEEE ***
               SMALL(1) = 8388608
               LARGE(1) = 2139095039
               RIGHT(1) = 864026624
               DIVER(1) = 872415232
               LOG10(1) = 1050288283
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*              *** VAX ***
               SMALL(1) = 128
               LARGE(1) = -32769
               RIGHT(1) = 13440
               DIVER(1) = 13568
               LOG10(1) = 547045274
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               LARGE(1) = 2147483647
               RIGHT(1) = 990904320
               DIVER(1) = 1007681536
               LOG10(1) = 1091781651
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
*              *** CONVEX C-1 ***
               SMALL(1) = 8388608
               LARGE(1) = 2147483647
               RIGHT(1) = 880803840
               DIVER(1) = 889192448
               LOG10(1) = 1067065499
            ELSE
               DO 10 L = 1, 3
                  J = SMALL(1) / 10000000
                  K = SMALL(1) - 10000000*J
                  IF (K .NE. T3E(L)) GO TO 20
                  SMALL(1) = J
 10               CONTINUE
*              *** CRAY T3E ***
               CALL I1MCRA(SMALL, K, 16, 0, 0)
               CALL I1MCRA(LARGE, K, 32751, 16777215, 16777215)
               CALL I1MCRA(RIGHT, K, 15520, 0, 0)
               CALL I1MCRA(DIVER, K, 15536, 0, 0)
               CALL I1MCRA(LOG10, K, 16339, 4461392, 10451455)
               GO TO 30
 20            CALL I1MCRA(J, K, 16405, 9876536, 0)
               IF (SMALL(1) .NE. J) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
*              *** CRAY 1, XMP, 2, AND 3 ***
               CALL I1MCRA(SMALL(1), K, 8195, 8388608, 1)
               CALL I1MCRA(LARGE(1), K, 24574, 16777215, 16777214)
               CALL I1MCRA(RIGHT(1), K, 16338, 8388608, 0)
               CALL I1MCRA(DIVER(1), K, 16339, 8388608, 0)
               CALL I1MCRA(LOG10(1), K, 16383, 10100890, 8715216)
               END IF
            END IF
 30      SC = 987
         END IF
*     SANITY CHECK
      IF (RMACH(4) .GE. 1.0) STOP 776
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'R1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      R1MACH = RMACH(I)
      RETURN
 9010 FORMAT(/' Adjust autodoubled R1MACH by getting data'/
     *' appropriate for your machine from D1MACH.')
 9020 FORMAT(/' Adjust R1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
* /* C source for R1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*float r1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return FLT_MIN;
*	  case 2: return FLT_MAX;
*	  case 3: return FLT_EPSILON/FLT_RADIX;
*	  case 4: return FLT_EPSILON;
*	  case 5: return log10((double)FLT_RADIX);
*	  }
*	fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
*	exit(1); return 0; /* else complaint of missing return value */
*}
      END
      SUBROUTINE I1MCRA(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END


      SUBROUTINE  ErrMsg2( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error

      LOGICAL       FATAL, MsgLim
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( 0, '(/,2A,/)' )  'Error,  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( 0, '(/,2A,/)' )  'Warning,  ', MESSAG
      ELSE
         WRITE ( 0,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     &   'They will no longer be printed  <<<<<<<', // )
      END
