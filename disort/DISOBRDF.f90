module m_tenstr_disort_disobrdf
  use m_tenstr_disort_errpack, only: ErrMsg
  use m_tenstr_disort_bdref, only: bdref
  use m_tenstr_disort_rdi1mach, only: D1MACH

  contains
      SUBROUTINE DISOBRDF(NSTR, USRANG, NUMU, UMU, &
     &           FBEAM, UMU0, LAMBER, ALBEDO, ONLYFL,&
     &           RHOQ, RHOU, EMUST, BEMST, DEBUG,&
     &           NPHI, PHI, PHI0, BDR_BEAM_ANALYTIC,   &
     &           BRDF_TYPE, BRDF_ARG, NMUG )

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! This program calculates the Fourier components of the BRDF
! specified in function BDREF in disobrdf.f and prepares it
! to use in DISORT.
!
! ** Version 3 upgrades:
!      
!     1) Prepare gaussian quadrature, weight and azimuth cosine series
!        before the azimuth loop to compute BRDF fourier components 
!      
!     2) Compute the analytic surface bidirectional reflectance 
!        at user angle (see DISOTEST3.f variable NUMU, UMU), 
!        which is used for new intensity correction in DISORT3, 
!        for more details see DISORT3 paper Section 3.5
!
!     3) Removed all fix-dimension symbolic variables  
!      
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! +-------------------------------------------------------------------+
!
! +-------------------------------------------------------------------+
!
!  NSTR        Stream No.
!
!  USRANG      Flag that determine whether there's user angle output
!          
!  NUMU        Number of userangle
!
!  UMU(IQ)     User defined polar angles
!
!  FBEAM       Incident Solar constant
!
!  UMU0        Solar zenith angle
!
!  LAMBER      Flag that determines whether there's Lambertian lower
!               boundary
!
!  ALBEDO                      Lambertian albedo 
!
!  ONLYFL                      Flag that decides whether to only output 
!
!  NUMU                        Number of user polar angles
!
!  RHOQ(MI,0:MI,0:NAZZ)        Quadrature fourier expanded BRDF: rho^m *
!                              1st index: output polar angle 
!                              2nd index: input polar angle, 0 is direct
!                              3rd index: azimuth index
!
!  RHOU(NUMU,0:MI,0:NAZZ)     User defined fourier expanded BRDF: rho^m 
!                              same index as RHOQ
!
!  EMUST(NUMU)                Directional emissivity at user angles
!
!  BEMST(MI)                   Directional emissivity at quadrature angl
!
!
! +-------------------------------------------------------------------+


!      USE PARAMETERS 
!      INTEGER  MAXCLY, MAXMOM, MAXPHI, MAXULV, MAXUMU, MAXCMU
!     ..
!     .. Scalar Arguments ..
      INTEGER   MI, NAZZ 
      LOGICAL   LAMBER, ONLYFL, USRANG
      LOGICAL   DEBUG
      INTEGER   NSTR, NUMU
      REAL      ALBEDO, FBEAM, UMU0
      INTEGER   BRDF_TYPE, NPHI 
      REAL      PHI0, PHI(NPHI)
      INTEGER   NMUG

!     ..
!     .. Array Arguments ..
      REAL      UMU(NUMU)
      REAL      BRDF_ARG(4)

!     ..
!     .. Local Scalars ..
      INTEGER   IQ, IU, J, K, MAZIM, NAZ, NN
      REAL      DELM0, PI

!     ..
!     .. Local Arrays ..
      REAL      BDR(NSTR/2, 0:NSTR/2), BEM(NSTR/2),&
     &          CMU(NSTR), CWT(NSTR), EMU(NUMU),&
     &          RMU(NUMU, 0:NSTR/2) 
      REAL      RHOQ(NSTR/2, 0:NSTR/2, 0:NSTR-1)
      REAL      RHOU(NSTR, 0:NSTR/2, 0:NSTR-1),&
     &          EMUST(NUMU), BEMST(NSTR/2)
      REAL      BDR_BEAM_ANALYTIC(NUMU,NPHI), DPHI
      REAL      GMU(NMUG), GWT(NMUG), COSMP(0:NSTR-1,NMUG/2) 

!     ..
!     .. External Functions ..
      !REAL      R1MACH
      !EXTERNAL  R1MACH

!     ..
!     .. External Subroutines ..
      !EXTERNAL   SURFAC2, ZEROIT2

!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS, ASIN, COS, FLOAT, LEN, MAX, SQRT

!     ..

      IF(DEBUG .EQV. .TRUE.) THEN
        PRINT *, '' 
        PRINT *, "Performing accurate BRDF calculation.."
      ENDIF

      PI     = 2.*ASIN( 1.0 )

!     ** Perform various setup operations

!     ** Calculate computational polar angle cosines and associated quad
!     ** weights for Gaussian quadrature on the interval (0,1) (upward).
      NN = NSTR / 2

!      IF(.NOT.LAMBER) THEN
      
      CALL QGAUSN2(NN, CMU, CWT)
!     ** Downward (neg) angles and weights
      DO 80 IQ = 1, NN
        CMU(IQ + NN) = -CMU(IQ)
        CWT(IQ + NN) =  CWT(IQ)
   80 CONTINUE

!     ** Version 3
!     ** Preparation before Fourier loop: COSMP GMU, GWT
      CALL QGAUSN2(NMUG/2, GMU, GWT)

      DO 10 K = 1, NMUG/2
        GMU(K + NMUG/2) = -GMU(K)
        GWT(K + NMUG/2) =  GWT(K)
        COSMP(0,K)      = 1.0;
        IF(.NOT.LAMBER .AND. .NOT.ONLYFL) THEN
          DO J = 1, NSTR-1
            COSMP(J,K) = COS(J*PI*GMU(K))
          ENDDO
        ENDIF
!       Initialize zeroth order separately to avoid random error when ON
!       11-27-2017         
!       IF(.NOT.LAMBER .AND. .NOT.ONLYFL) THEN
!         DO J = 0, NSTR-1
!           COSMP(J,K) = COS(J*PI*GMU(K))
!         ENDDO
!       ENDIF
10    CONTINUE

!     ========  BEGIN LOOP OVER AZIMUTH  ========

      NAZ = NSTR - 1

!     ** Azimuth-independent case
      IF( FBEAM.EQ.0.0 .OR. ABS(1.-UMU0).LT.1.E-5 .OR. ONLYFL .OR.&
     &   ( NUMU.EQ.1 .AND. ABS(1.-UMU(1)).LT.1.E-5 ) .OR.&
     &   ( NUMU.EQ.1 .AND. ABS(1.+UMU(1)).LT.1.E-5 ) .OR.&
     &   ( NUMU.EQ.2 .AND. ABS(1.+UMU(1)).LT.1.E-5 .AND.&
     &     ABS(1.-UMU(NUMU)).LT.1.E-5 ) ) THEN     
!     &     ABS(1.-UMU(2)).LT.1.E-5 ) ) THEN
        NAZ = 0
      ENDIF

      DO 180 MAZIM = 0, NAZ
       
        WRITE(6,990,ADVANCE='NO') MAZIM
        IF(MAZIM .EQ. NAZ) THEN
          WRITE(6,*) ""
        ENDIF

        IF(MAZIM .EQ. 0) THEN
          DELM0  = 1.0
        ELSE
          DELM0  = 0.0
        ENDIF

        MI = NSTR / 2
        NAZZ = NSTR - 1
        CALL SURFAC2(ALBEDO, DELM0, CMU, FBEAM, LAMBER, MI, MAZIM,&
     &               NUMU, NN, NUMU, ONLYFL, PI, UMU, UMU0,&
     &               USRANG, BDR, EMU, BEM, RMU, &
     &               RHOQ, RHOU, EMUST, BEMST, NAZZ, DEBUG,&
     &               NSTR, NMUG, GMU, GWT, COSMP, &
     &               BRDF_TYPE, BRDF_ARG)

180   CONTINUE

!     ** Version 3
!     ** Compute analytic surface bidirectional reflectance at user angl
      IF(.NOT. LAMBER) THEN
        DO 190 IU = 1, NUMU
          DO 200 J = 1, NPHI
            IF(UMU(IU) .LE. 0.) THEN
              BDR_BEAM_ANALYTIC(IU, J) = -1.0
            ELSE
              DPHI = (PHI(J) - PHI0)*PI/180.
              BDR_BEAM_ANALYTIC(IU, J) = BDREF(UMU(IU), UMU0, DPHI,&
     &                                          BRDF_TYPE, BRDF_ARG) 
            ENDIF
200       CONTINUE
190     CONTINUE
      ENDIF
      
!      ENDIF

!      PRINT *, "NSTR=", NSTR, "UMU0=",UMU0
!      OPEN( UNIT = 58, FILE = 'BRDF.EM' )
!      OPEN( UNIT = 59, FILE = 'BRDF.RHOQ' )
!      OPEN( UNIT = 60, FILE = 'BRDF.RHOU' )
!      OPEN( UNIT = 61, FILE = 'CMU.dat' )
!      DO IQ = 1, NN
!         WRITE(58,578) IQ, BEMST(IQ), EMUST(IQ)
!      ENDDO
!     RHOQ(OUT,IN,NAZZ)
!     For each angle OUT, IQ, you have several angles IN, JQ
!       DO IQ=1,NN
!         WRITE(61,581) CMU(IQ)
!       ENDDO
!       DO IQ=1,NN
!         DO JQ=0,NN
!           IF(JQ.EQ.0) THEN
!             PRINT*, UMU0,IQ,JQ,(RHOQ(IQ,JQ,IQ2),IQ2=0,NSTR-1)
!           ELSE                  
!             PRINT *, CMU(JQ),IQ,JQ,(RHOQ(IQ,JQ,IQ2),IQ2=0,NSTR-1)
!           ENDIF
!             WRITE(59,580) (RHOQ(IQ,JQ,IQ2),IQ2=0,NSTR-1)
!         ENDDO
!       ENDDO 
!       DO IU=1,NUMU
!         DO JQ=1,NN
!           WRITE(60,582) UMU(IU),IU,JQ,(RHOU(IU,JQ,IQ2),IQ2=1,NSTR-1)
!           WRITE(60,580) (RHOU(IU,JQ,IQ2),IQ2=0,NSTR-1)
!         ENDDO
!       ENDDO
!578     FORMAT(I4,27E16.8)
!579     FORMAT(I2,I2,27E12.4)
!580     FORMAT(500E12.4)
!581     FORMAT(1E16.8)
!582     FORMAT(E12.4,I2,1X,I2,500E12.4)

990   FORMAT(I3)

      RETURN
      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --------------------------------------------------------------------- 
      SUBROUTINE SURFAC2(ALBEDO, DELM0, CMU, FBEAM, LAMBER, MI, MAZIM,&
     &                   MXUMU, NN, NUMU, ONLYFL, PI, UMU, UMU0,&
     &                   USRANG, BDR, EMU, BEM, RMU, &
     &                   RHOQ, RHOU, EMUST, BEMST, NAZZ, DEBUG,&
     &                   NSTR, NMUG, GMU, GWT, COSMP, &
     &                   BRDF_TYPE, BRDF_ARG)

!       Computes user's surface bidirectional properties, STWL(41)
!
!   I N P U T     V A R I A B L E S:
!
!       CMU    :  Computational polar angle cosines (Gaussian)
!
!       DELM0  :  Kronecker delta, delta-sub-m0
!
!       MAZIM  :  Order of azimuthal component
!
!       NN     :  Order of Double-Gauss quadrature (NSTR/2)
!
!       (Remainder are 'DISORT' input variables)
!
!    O U T P U T     V A R I A B L E S:
!
!       BDR :  Fourier expansion coefficient of surface bidirectional
!                 reflectivity (computational angles)
!
!       RMU :  Surface bidirectional reflectivity (user angles)
!
!       BEM :  Surface directional emissivity (computational angles)
!
!       EMU :  Surface directional emissivity (user angles)
!
!    I N T E R N A L     V A R I A B L E S:
!
!       DREF   :  Directional reflectivity
!
!       NMUG   :  Number of angle cosine quadrature points on (-1,1)
!                 for integrating bidirectional reflectivity to get
!                 directional emissivity (it is necessary to use a
!                 quadrature set distinct from the computational angles,
!                 because the computational angles may not be dense
!                 enough -- i.e. 'NSTR' may be too small-- to give an
!                 accurate approximation for the integration).
!
!       GMU    :  The 'NMUG' angle cosine quadrature points on (0,1)
!
!       GWT    :  The 'NMUG' angle cosine quadrature weights on (0,1)
!
!   Called by- DISORT
!   Calls- QGAUSN, BDREF, ZEROIT
!+---------------------------------------------------------------------+

!     .. Parameters ..
      INTEGER   NMUG, NSTR

!     ..
!     .. Scalar Arguments ..
      LOGICAL   LAMBER, ONLYFL, USRANG
      LOGICAL   DEBUG
      INTEGER   MAZIM, MI, MXUMU, NN, NUMU
      REAL      ALBEDO, DELM0, FBEAM, PI, UMU0
      INTEGER   BRDF_TYPE

!     ..
!     .. Array Arguments ..
      REAL      BDR( NN, 0:NN ), BEM( NN ), CMU( * ), EMU( NUMU ),&
     &          RMU( NUMU, 0:NN ), UMU( * ), &
     &          RHOQ(MI,0:MI,0:NAZZ), RHOU(MXUMU,0:MI,0:NAZZ),&
     &          EMUST(MXUMU), BEMST(MI)
      REAL      BRDF_ARG(4)

!     ..
!     .. Local Scalars ..
      INTEGER   IQ, IU, JG, JQ, K
      REAL      DREF, SUM

!     ..
!     .. Local Arrays ..
      REAL      GMU( NMUG ), GWT( NMUG ), COSMP(0:NSTR-1, NMUG/2)

!     ..
!     .. External Functions ..
      !REAL      BDREF
      !EXTERNAL  BDREF

!     ..
!     .. External Subroutines ..
      !EXTERNAL  QGAUSN2, ZEROIT2

!     ..
!     .. Intrinsic Functions ..
      INTRINSIC COS

!     ..

      IF(DEBUG .EQV. .TRUE.) THEN
        OPEN( UNIT = 57, FILE = 'BRDF.OUTPUT' )
      ENDIF

      IF(MAZIM .GT. 2*NN) THEN
        CALL ERRMSG( 'MXSTR TOO&
                     &LOW - INCREASE MXSTR OR DECREASE NSTR', .TRUE. )
      ENDIF

      CALL ZEROIT2(BDR, NN*(NN+1))
      CALL ZEROIT2(BEM, NN)

!     ** Compute Fourier expansion coefficient of surface bidirectional
!     ** reflectance at computational angles Eq. STWL (41).

      IF(LAMBER .AND. MAZIM .EQ. 0) THEN

        DO 30 IQ = 1, NN
          BEM(IQ) = 1.0 - ALBEDO
          DO 20 JQ = 0, NN
            BDR(IQ, JQ) = ALBEDO
20        CONTINUE
30      CONTINUE

!     ** Version 3: Compute BRDF accurately and efficiently
      ELSEIF(.NOT.LAMBER) THEN

        DO 70 IQ = 1, NN
          DO 50 JQ = 1, NN
            SUM  = 0.0
            DO 40 K = 1, NMUG/2
              SUM = SUM + GWT(K)&
     &             * BDREF(CMU(IQ), CMU(JQ), PI*GMU(K),&
     &                      BRDF_TYPE, BRDF_ARG)&
     &             * COSMP(MAZIM, K)
40          CONTINUE

!         ** Version 3: removed 0.5 and added PI
          BDR(IQ, JQ)         = (2. - DELM0) * SUM
          RHOQ(IQ, JQ, MAZIM) = BDR(IQ, JQ) * PI
!          BDR( IQ, JQ ) = 0.5 * ( 2. - DELM0 ) * SUM
!          RHOQ( IQ, JQ, MAZIM ) = BDR( IQ, JQ )
50        CONTINUE

          IF(FBEAM .GT. 0.0) THEN
            SUM  = 0.0
            DO 60 K = 1, NMUG/2
              SUM = SUM + GWT(K)&
     &             * BDREF(CMU(IQ), UMU0,PI*GMU(K),&
     &                      BRDF_TYPE, BRDF_ARG) &
     &             * COSMP(MAZIM, K)
60          CONTINUE
!           ** Version 3: removed 0.5 and added PI.
            BDR(IQ, 0)         = (2. - DELM0) * SUM
            RHOQ(IQ, 0, MAZIM) = BDR(IQ, 0) * PI
          ENDIF
70      CONTINUE


!       ** For 0th azimuth component, integrate bidirectional reflectivi
!       ** reflection polar angle cosines -CMU- and incident angle cosin
!       ** to get directional emissivity at computational angle cosines 
        IF(MAZIM .EQ. 0) THEN

          DO 100 IQ = 1, NN
            DREF  = 0.0
            DO 90 JG = 1, NMUG
              SUM  = 0.0
              DO 80 K = 1, NMUG / 2
                SUM = SUM + GWT(K) * GMU(K)&
     &               * BDREF(CMU(IQ), GMU(K),&
     &                        PI*GMU(JG), BRDF_TYPE, BRDF_ARG)
80            CONTINUE
            DREF = DREF + GWT(JG)*SUM
90          CONTINUE
          BEM(IQ)   = 1.0 - DREF
          BEMST(IQ) = BEM(IQ)
100       CONTINUE
        ENDIF

      ENDIF

!     ** Compute Fourier expansion coefficient of surface bidirectional
!     ** reflectance at user angles Eq. STWL (41).
      IF( .NOT.ONLYFL .AND. USRANG ) THEN

        CALL ZEROIT2(EMU, NUMU)
        CALL ZEROIT2(RMU, NUMU*(NSTR/2+1))

        DO 170 IU = 1, NUMU

          IF(UMU(IU) .GT. 0.0) THEN

            IF(LAMBER .AND. MAZIM .EQ. 0) THEN
              DO 110 IQ = 0, NN
                RMU(IU, IQ) = ALBEDO
110           CONTINUE
              EMU(IU) = 1.0 - ALBEDO

            ELSEIF(.NOT.LAMBER) THEN
              DO 130 IQ = 1, NN
                SUM  = 0.0
                DO 120 K = 1, NMUG/2
                  SUM = SUM + GWT(K)&
     &                 * BDREF(UMU(IU), CMU(IQ),&
     &                          PI*GMU(K), BRDF_TYPE, BRDF_ARG) &
     &                 * COSMP(MAZIM, K)
120             CONTINUE
!               ** Version 3: removed 0.5 and added PI.
                RMU( IU, IQ ) =  ( 2. - DELM0 ) * SUM
                RHOU( IU, IQ, MAZIM ) = RMU( IU, IQ )*PI
130           CONTINUE

              IF(FBEAM .GT. 0.0) THEN
                SUM = 0.0
                DO 140 K = 1, NMUG/2
                  SUM = SUM + GWT(K)&
     &                 *BDREF(UMU(IU), UMU0,&
     &                         PI*GMU(K), BRDF_TYPE, BRDF_ARG) &
     &                 *COSMP( MAZIM, K )
140             CONTINUE
!               ** Version 3: removed 0.5 and added PI.
                RMU(IU, 0)         = (2. - DELM0)*SUM
                RHOU(IU, 0, MAZIM) = RMU(IU, 0)*PI
              ENDIF

!             ** For 0th azimuth component, integrate bidirectional refl
!             ** at reflection angle cosines -UMU- and incident angle co
!             ** -GMU- to get directional emissivity at user angle cosin
!             ** -UMU-.
              IF(MAZIM .EQ. 0) THEN
                DREF  = 0.0
                DO 160 JG = 1, NMUG
                  SUM  = 0.0
                  DO 150 K = 1, NMUG / 2
                    SUM = SUM + GWT(K)*GMU(K)&
     &                   *BDREF(UMU(IU), GMU(K), &
     &                           PI*GMU(JG), BRDF_TYPE, BRDF_ARG)
150               CONTINUE
                  DREF  = DREF + GWT( JG ) * SUM
160             CONTINUE

                EMU(IU)   = 1.0 - DREF
                EMUST(IU) = EMU ( IU)
              ENDIF
            ENDIF
          ENDIF
170     CONTINUE

      ENDIF

      IF(DEBUG .EQV. .TRUE.) THEN
        IF(MAZIM.EQ.0) THEN
          WRITE(57,577) ( BEM(IQ), IQ=1,NN )
          WRITE(57,577) ( EMU(IQ), IQ=1,NN )
        ENDIF
        DO 571 IQ = 1, NN 
          WRITE(57,577) (RHOQ( IQ, JQ, MAZIM), JQ=0, NN) 
571     CONTINUE
        DO 572 IQ = 1, NUMU
          WRITE(57,577) (RHOU( IQ, JQ, MAZIM), JQ=0, NN)
572     CONTINUE
      ENDIF
                 
577   FORMAT(500E16.8)
     
      RETURN
      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --------------------------------------------------------------------- 
      SUBROUTINE QGAUSN2( M, GMU, GWT )

!     Compute weights and abscissae for ordinary Gaussian quadrature
!     on the interval (0,1);  that is, such that
!
!       sum(i=1 to M) ( GWT(i) f(GMU(i)) )
!
!     is a good approximation to
!
!       integral(0 to 1) ( f(x) dx )
!
!   INPUT :    M       order of quadrature rule
!
!   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
!             GWT(I)   array of weights (I = 1 TO M)
!
!   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
!                   Integration, Academic Press, New York, pp. 87, 1975
!
!   METHOD:  Compute the abscissae as roots of the Legendre
!            polynomial P-sub-M using a cubically convergent
!            refinement of Newton's method.  Compute the
!            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!            that Newton's method can very easily diverge; only a
!            very good initial guess can guarantee convergence.
!            The initial guess used here has never led to divergence
!            even for M up to 1000.
!
!   ACCURACY:  relative error no better than TOL or computer
!              precision (machine epsilon), whichever is larger
!
!   INTERNAL VARIABLES:
!
!    ITER      : number of Newton Method iterations
!    MAXIT     : maximum allowed iterations of Newton Method
!    PM2,PM1,P : 3 successive Legendre polynomials
!    PPR       : derivative of Legendre polynomial
!    P2PRI     : 2nd derivative of Legendre polynomial
!    TOL       : convergence criterion for Legendre poly root iteration
!    X,XI      : successive iterates in cubically-convergent version
!                of Newtons Method (seeking roots of Legendre poly.)
!
!   Called by- DREF, SETDIS, SURFAC
!   Calls- D1MACH, ERRMSG
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..
      INTEGER M

!     ..
!     .. Array Arguments ..
      REAL GMU(M), GWT(M)

!     ..
!     .. Local Scalars ..
      INTEGER ITER, K, LIM, MAXIT, NN, NP1
      DOUBLE PRECISION CONA, PI, T
      DOUBLE PRECISION EN, NNP1, ONE, P, P2PRI, PM1, PM2, PPR, PROD,&
     &                 TMP, TOL, TWO, X, XI

!     ..
!     .. External Functions ..
      !DOUBLE PRECISION D1MACH
      !EXTERNAL  D1MACH

!     ..
!     .. External Subroutines ..
      !EXTERNAL  ERRMSG

!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS, ASIN, COS, FLOAT, MOD, TAN

!     ..

      SAVE      PI, TOL
      DATA      PI / 0.D0 / , MAXIT / 1000 / , ONE / 1.D0 / ,&
     &          TWO / 2.D0 /


      IF(PI .EQ. 0.D0) THEN
        PI   = 2.D0*DASIN( 1.D0 )
        TOL  = 10.*D1MACH( 4 )
      ENDIF

      IF(M.LT.1) THEN
        CALL ERRMSG('QGAUSN2 -- Bad value for M',.True.)
      ENDIF

      IF(M.EQ.1) THEN
        GMU(1) = 0.5
        GWT(1) = 1.0
        RETURN
      ENDIF

      EN   = DBLE(M)
      NP1  = M + 1
      NNP1 = DBLE(M*NP1)
      CONA = DBLE( M - 1 ) / ( 8*M**3 ) 

      LIM  = M / 2

      DO 30 K = 1, LIM

!       ** Initial guess for k-th root of Legendre polynomial, from
!       ** Davis/Rabinowitz (2.7.3.3a).
        T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
        X  = DCOS( T + CONA / DTAN( T ) )
        ITER = 0

!       ** Upward recurrence for Legendre polynomials
10      CONTINUE
        ITER = ITER + 1
        PM2  = ONE
        PM1  = X

        P = 0D0
        DO 20 NN = 2, M
          P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / NN
          PM2  = PM1
          PM1  = P
20      CONTINUE

!       ** Newton Method
        TMP    = ONE / ( ONE - X**2 )
        PPR    = EN*( PM2 - X*P )*TMP
        P2PRI  = ( TWO*X*PPR - NNP1*P )*TMP
        XI     = X - ( P / PPR )*( ONE +&
     &           ( P / PPR )*P2PRI / ( TWO*PPR ) )

!       ** Check for convergence
        IF( DABS( XI - X ).GT.TOL ) THEN

          IF( ITER.GT.MAXIT ) THEN
            CALL ERRMSG('QGAUSN2 - max iteration count',.True.)
          ENDIF

          X = XI
          GOTO 10
        ENDIF

!       ** Iteration finished--calculate weights, abscissae for (-1,1)
        GMU( K ) = - REAL( X )
        GWT( K ) = REAL( TWO / ( TMP*( EN*PM2 )**2 ) )
        GMU( NP1 - K ) = -GMU( K )
        GWT( NP1 - K ) = GWT( K )
30    CONTINUE

!     ** Set middle abscissa and weight for rules of odd order
      IF(MOD( M,2 ) .NE. 0) THEN
        GMU(LIM + 1) = 0.0
        PROD = ONE
        DO 40 K = 3, M, 2
          PROD   = PROD * K / ( K - 1 )
40      CONTINUE
        GWT(LIM + 1) = REAL( TWO / PROD**2 )
      ENDIF

!     ** Convert from (-1,1) to (0,1)
      DO 50 K = 1, M
        GMU(K) = 0.5*GMU(K) + 0.5
        GWT(K) = 0.5*GWT(K)
50    CONTINUE

      RETURN
      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! --------------------------------------------------------------------- 
      SUBROUTINE ZEROIT2( A, LENGTH )

!     Zeros a real array A having LENGTH elements
!
!   Called by- DISORT, ALBTRN, SOLVE1, SURFAC, SETMTX, SOLVE0, FLUXES
! --------------------------------------------------------------------

!     .. Scalar Arguments ..
      INTEGER   LENGTH

!     ..
!     .. Array Arguments ..
      REAL      A( LENGTH )

!     ..
!     .. Local Scalars ..
      INTEGER   L

!     ..

      DO 10 L = 1, LENGTH
        A(L) = 0.0
10    CONTINUE

      RETURN
      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end module
