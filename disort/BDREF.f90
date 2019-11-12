! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! $Rev: 55 $ $Date: 2014-12-3112:16:59 -0500 (Wed, 31 Dec 2014) $
! FORTRAN 77
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
module m_tenstr_disort_bdref
  use m_tenstr_disort_errpack, only: ErrMsg
  contains
      REAL FUNCTION BDREF(MU, MUP, DPHI,&
     &                     BRDF_TYPE, BRDF_ARG)

!     Supplies surface bi-directional reflectivity.
!
!     NOTE 1: Bidirectional reflectivity in DISORT is defined
!             by Eq. 39 in STWL.
!     NOTE 2: Both MU and MU0 (cosines of reflection and incidence
!             angles) are positive.
!
!  INPUT:
!
!    MU     : Cosine of angle of reflection (positive)
!
!    MUP    : Cosine of angle of incidence (positive)
!
!    DPHI   : Difference of azimuth angles of incidence and reflection
!                (radians)
!
!  LOCAL VARIABLES:
!
!    IREF   : bidirectional reflectance options
!             1 - Hapke's BDR model
!             2 - Cox-Munk BDR model
!             3 - RPV BDR model
!             4 - Ross-Li BDR model
!
!    B0     : empirical factor to account for the finite size of
!             particles in Hapke's BDR model
!
!    B      : term that accounts for the opposition effect
!             (retroreflectance, hot spot) in Hapke's BDR model
!
!    CTHETA : cosine of phase angle in Hapke's BDR model
!
!    GAMMA  : albedo factor in Hapke's BDR model
!
!    H0     : H( mu0 ) in Hapke's BDR model
!
!    H      : H( mu ) in Hapke's BDR model
!
!    HH     : angular width parameter of opposition effect in Hapke's
!             BDR model
!
!    P      : scattering phase function in Hapke's BDR model
!
!    THETA  : phase angle (radians); the angle between incidence and
!             reflection directions in Hapke's BDR model
!
!    W      : single scattering albedo in Hapke's BDR model
!
!
!   Called by- DREF, SURFAC
! +-------------------------------------------------------------------+
!     .. Scalar Arguments ..
      REAL      DPHI, MU, MUP, BRDF_ARG(4)
      INTEGER   BRDF_TYPE
!     ..
!     .. Local Scalars ..
      INTEGER   IREF
      REAL      B0, H0, HH, W
      REAL      PWS, REFRAC_INDEX, BDREF_F
      REAL      PI
      REAL      RHO0, KAPPA, G  
      REAL      K_ISO, K_VOL, K_GEO, ALPHA0 
      LOGICAL   DO_SHADOW
!     ..
!     .. External Subroutines ..
      !EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC COS, SQRT
!     ..

      PI   = 2.*ASIN(1.)

      IREF = BRDF_TYPE

!     ** 1. Hapke BRDF
      IF ( IREF.EQ.1 ) THEN

!       ** Hapke's BRDF model (times Pi/Mu0) (Hapke, B., Theory of refle
!       ** and emittance spectroscopy, Cambridge University Press, 1993,
!       ** 8.89 on page 233. Parameters are from Fig. 8.15 on page 231, 
!       ** for w.)

        B0 = BRDF_ARG(1) !1.0
        HH = BRDF_ARG(2) !0.06
        W  = BRDF_ARG(3) !0.6

        CALL BRDF_HAPKE(MUP, MU, DPHI,&
     &                  B0, HH, W, PI,&
     &                  BDREF)

!     ** 2. Cox-Munk BRDF
      ELSEIF(IREF.EQ.2) THEN

!        PRINT *, "Calling oceabrdf"

        PWS          =  BRDF_ARG(1)
        REFRAC_INDEX =  BRDF_ARG(2)

        IF(BRDF_ARG(3) .EQ. 1) THEN
          DO_SHADOW = .TRUE.
        ELSEIF(BRDF_ARG(3) .EQ. 0) THEN
          DO_SHADOW = .FALSE.
        ELSE
          PRINT *, "ERROR SHADOW ARGUMENTS"
        ENDIF

        CALL OCEABRDF2(DO_SHADOW,&
     &                 REFRAC_INDEX, PWS, &
     &                 MUP, MU, DPHI,&
     &                 BDREF_F)

        BDREF = BDREF_F

!     ** 3. RPV BRDF
      ELSEIF(IREF .EQ. 3) THEN

        RHO0  =  BRDF_ARG(1) !0.027
        KAPPA =  BRDF_ARG(2) !0.647
        G     =  BRDF_ARG(3) !-0.169   !asymmetry factor for HG
        H0    =  BRDF_ARG(4) !0.100

        CALL BRDF_RPV(MUP, MU, DPHI,&
     &                RHO0, KAPPA, G, H0,&
     &                BDREF_F)

        BDREF = BDREF_F

!     ** 4. Ross-Li BRDF
      ELSEIF(IREF .EQ. 4) THEN
        
        K_ISO  = BRDF_ARG(1)   !0.200
        K_VOL  = BRDF_ARG(2)   !0.020
        K_GEO  = BRDF_ARG(3)   !0.300
        ALPHA0 = 1.5*pi/180.

        CALL BRDF_ROSSLI(MUP, MU, DPHI,&
     &                   K_ISO, K_VOL, K_GEO,&
     &                   ALPHA0,&
     &                   BDREF_F)

        BDREF = BDREF_F

        IF(BDREF .LT. 0.00) THEN
          BDREF = 0.00
        ENDIF

      ELSE

        CALL ERRMSG( 'BDREF--Need to supply surface BDRF model',&
     &                 .TRUE.)

      ENDIF

      RETURN
      END FUNCTION
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! +--------------------------------------------------------------------
      SUBROUTINE BRDF_HAPKE ( MUP, MU, DPHI,&
     &                        B0, HH, W, PI,&
     &                        BRDF )

! +--------------------------------------------------------------------
! Hapke "Theory of Reflectance and Emittance Spectroscopy" Chapter 10, P
! Eq. (10.2).
! Version 3 fix: definition of phase angle / scattering angle see DISORT
! paper Eqs. (25-26).
! +--------------------------------------------------------------------
      IMPLICIT NONE
      REAL MUP, MU, DPHI
      REAL B0, HH, W, PI
      REAL BRDF
      REAL CALPHA, ALPHA, P, B, H0, GAMMA, H

      CALPHA = MU * MUP - (1.-MU**2)**.5 * (1.-MUP**2)**.5&
     &         * COS( DPHI )

      ALPHA = ACOS( CALPHA )

      P     = 1. + 0.5 * CALPHA

      B     = B0 * HH / ( HH + TAN( ALPHA/2.) )

      GAMMA = SQRT( 1. - W )
      H0   = ( 1. + 2.*MUP ) / ( 1. + 2.*MUP * GAMMA )
      H    = ( 1. + 2.*MU ) / ( 1. + 2.*MU * GAMMA )

!     ** Version 3: add factor PI
      BRDF = W / (4.*PI) / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )
!     BRDF = W / 4. / (MU+MUP) * ( (1.+B)* P + H0 * H - 1.0 )

      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! +--------------------------------------------------------------------
      SUBROUTINE BRDF_RPV(MU_I, MU_R, DPHI,&
     &                    RHO0, KAPPA, G_HG, H0,&
     &                    BRDF)

! +--------------------------------------------------------------------
! DISORT Version 3: RPV BRDF
!   Input:
!
!   MU_I:  absolute cosine of incident polar angle (positive)
!   MU_R:  absolute cosine of reflected polar angle (positive)
!   DPHI:  relative azimuth to incident vector; (pi - dphi), sun-view re
!          azimuth sun located at phi = 180, while incident solar beam l
!          at phi = 0
!   RHO0:  RPV BRDF parameter, control reflectance
!   KAPPA: PRV BRDF parameter, control anisotropy
!   G:     RPV BRDF parameter, H-G asymmetry factor
!   H0:    RPV BRDF parameter, control hot spot (back scattering directi
!
!   Output:
!
!   BRDF:  RPV BRDF
! +--------------------------------------------------------------------
      IMPLICIT NONE
      REAL MU_I, MU_R, DPHI
      REAL RHO0, KAPPA, G_HG, H0
      REAL BRDF
      REAL PI
      REAL COS_ALPHA
      REAL SIN_I, SIN_R, TAN_I, TAN_R
      REAL G_SQ, G, F

      PI    = 2.*ASIN(1.)

      SIN_I = SQRT(1. - MU_I*MU_I)
      SIN_R = SQRT(1. - MU_R*MU_R)
      TAN_I = SIN_I/MU_I
      TAN_R = SIN_R/MU_R

      COS_ALPHA = MU_I*MU_R - SIN_I*SIN_R&
     & *COS(DPHI)

      G_SQ = TAN_I*TAN_I + TAN_R*TAN_R &
     &    + 2.*TAN_I*TAN_R*COS(DPHI)

!     ** hot spot
      G = SQRT(G_SQ)

!     ** HG phase function
      F = (1. - G_HG*G_HG)/&
     &     (1+G_HG*G_HG+2.*G_HG*COS_ALPHA)**1.5


!     ** BRDF semiempirical function
      BRDF = RHO0 &
     &      * (MU_I*MU_R*(MU_I+MU_R))**(KAPPA-1.)&
     &      * F&
     &      * (1. + ((1.-H0)/(1.+G)))

      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! +--------------------------------------------------------------------
      SUBROUTINE BRDF_ROSSLI(MU_I, MU_R, DPHI,&
     &                       K_ISO, K_VOL, K_GEO,&
     &                       ALPHA0,&
     &                       BRDF)

! +--------------------------------------------------------------------
! Version 3: Ross-Li BRDF
!   Input:
!
!   MU_I:    absolute cosine of incident polar angle (positive)
!   MU_R:    absolute cosine of reflected polar angle (positive)
!   DPHI:  relative azimuth to incident vector; (pi - dphi), sun-view re
!          azimuth sun located at phi = 180, while incident solar beam l
!          at phi = 0
!   K_ISO:   BRDF parameter, isotropic scattering kernel
!   K_VOL:   BRDF parameter, volume scattering kernel
!   K_GEO:   BRDF parameter, geometry scattering kernel
!   ALPHA0:  BRDF parameter, control hot spot (back scattering direction
!
!   Output:
!   BRDF:  Ross-Li BRDF
!
! +--------------------------------------------------------------------
      IMPLICIT NONE
      REAL MU_I, MU_R, DPHI
      REAL F_GEO, F_VOL
      REAL K_ISO, K_GEO, K_VOL
      REAL RATIO_HB, RATIO_BR
      REAL BRDF
      REAL PI
      REAL COS_ALPHA, SIN_ALPHA
      REAL COS_ALPHA1
      REAL ALPHA
      REAL SIN_I, SIN_R, TAN_I, TAN_R
      REAL SIN_I1, SIN_R1, COS_I1, COS_R1, TAN_I1, TAN_R1
      REAL G_SQ, COS_T, T       
      REAL C, ALPHA0
! +--------------------------------------------------------------------

!      PRINT *, MU_I, MU_R, DPHI,
!     &        K_ISO, K_GEO, K_VOL,
!     &        THETA0
!      PRINT *,

      RATIO_HB = 2.
      RATIO_BR = 1.
      PI       = 2.*ASIN(1.)

      SIN_I = SQRT(1. - MU_I*MU_I)
      SIN_R = SQRT(1. - MU_R*MU_R)
      TAN_I = SIN_I/MU_I
      TAN_R = SIN_R/MU_R

      COS_ALPHA = MU_I*MU_R - SIN_I*SIN_R&
     & *COS(DPHI)
      SIN_ALPHA = SQRT(1. - COS_ALPHA*COS_ALPHA)
      ALPHA = ACOS(COS_ALPHA)

!     ** Compute KERNEL RossThick
      C     = 1. + 1./(1.+ALPHA/ALPHA0)
      F_VOL = 4./(3.*PI) * (1./(MU_I+MU_R))&
     &       * ((PI/2. - ALPHA)*COS_ALPHA+SIN_ALPHA)*C - 1./3.

!      K1 = ((PI/2. - ALPHA)*COS_ALPHA + SIN_ALPHA)
!     &       /(MU_I + MU_R) - PI/4.


!     ** Compute KERNEL LSR
      TAN_I1 = RATIO_BR * TAN_I
      TAN_R1 = RATIO_BR * TAN_R
      SIN_I1 = TAN_I1/SQRT(1.+ TAN_I1*TAN_I1)
      SIN_R1 = TAN_R1/SQRT(1.+ TAN_R1*TAN_R1)
      COS_I1 = 1./SQRT(1.+ TAN_I1*TAN_I1)
      COS_R1 = 1./SQRT(1.+ TAN_R1*TAN_R1)

      COS_ALPHA1 = COS_I1*COS_R1 - SIN_I1*SIN_R1&
     &            *COS(DPHI)

      G_SQ = TAN_I1*TAN_I1 + TAN_R1*TAN_R1 &
     &      + 2.*TAN_I1*TAN_R1*COS(DPHI)

!      M = 1./COS_I1 + 1./COS_R1

      COS_T = RATIO_HB *(COS_I1*COS_R1)/(COS_I1+COS_R1)&
     &       *SQRT(G_SQ + (TAN_I1*TAN_R1*SIN(DPHI))**2)
  
      IF(COS_T .LE. 1. .AND. COS_T .GE. -1.) THEN
        T = ACOS(COS_T)
      ELSE
        T = 0.
      ENDIF

      F_GEO = (COS_I1+COS_R1)/(PI*COS_I1*COS_R1)*(T-SIN(T)*COS(T)-PI)   &
     &       + (1.+ COS_ALPHA1)/(2.*COS_I1*COS_R1)

!     Compute BRDF

!      PRINT *, RATIO_HB, D_SQ, 
!     &    TAN_I1*TAN_R1*SIN(DPHI),
!     &    M, COS_T

!      BRDF = K1
      BRDF = K_ISO + K_GEO*F_GEO + K_VOL*F_VOL

      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! +--------------------------------------------------------------------
      SUBROUTINE OCEABRDF2&
     &       ( DO_SHADOW, &
     &         REFRAC_INDEX, WS,&
     &         MU_I, MU_R, DPHI,&
     &         BRDF)

! +--------------------------------------------------------------------
! Version 3: 1D Gaussian Rough Ocean BRDF
!   Input:
!
!   mu_i:         absolute cosine of incident polar angle (positive)
!   mu_r:         absolute cosine of reflected polar angle (positive)
!   dphi:         relative azimuth (radians) 
!   do_shadow:    BRDF parameter, open/close shadow effect 
!   refrac_index: BRDF parameter, refractive index of boundary media (wa
!   ws:           BRDF parameter, wind speed (m/s)
!
!   Output:
!
!   brdf:         1D Gaussian Rough Ocean BRDF
!          
! +--------------------------------------------------------------------
      LOGICAL  DO_SHADOW
      REAL     REFRAC_INDEX, WS
      REAL     SIN_I, SIN_R, MU_I, MU_R, DPHI, BRDF
      REAL     COS_THETA, SIGMA_SQ, MU_N_SQ, P
      REAL     N_I, N_T, COS_LI, COS_LT, SIN_LI, SIN_LT
      REAL     R_S, R_P, R
      REAL     SHADOW
      REAL     PI

      PI = 2.*ASIN(1.)

!     ** Cox Munk slope distribution
      SIN_I = SQRT(1. - MU_I*MU_I)
      SIN_R = SQRT(1. - MU_R*MU_R)

      COS_THETA = -MU_I*MU_R + SIN_I*SIN_R*COS(DPHI)
      MU_N_SQ   = (MU_I + MU_R)*(MU_I + MU_R)/(2.*(1.-COS_THETA))   

      SIGMA_SQ  = 0.003 + 0.00512*WS

      P = 1./(PI*SIGMA_SQ) * EXP( -(1-MU_N_SQ)/(SIGMA_SQ*MU_N_SQ) )

!     ** Fresnel reflectance

      N_I = 1.0
      N_T = REFRAC_INDEX

      SIN_LI = SQRT( 1.-0.5*(1.-COS_THETA) ) 
      COS_LI = SQRT( 0.5*(1.-COS_THETA) ) 
      SIN_LT = N_I*SIN_LI/N_T
      COS_LT = SQRT(1. - SIN_LT*SIN_LT)

      R_S = (N_I*COS_LI-N_T*COS_LT)/(N_I*COS_LI+N_T*COS_LT)
      R_P = (N_T*COS_LI-N_I*COS_LT)/(N_I*COS_LT+N_T*COS_LI)

      R = 0.5*(R_S*R_S + R_P*R_P)

!     ** Rough surface BRDF
      BRDF = (P*R)/(4.*MU_I*MU_R*MU_N_SQ*MU_N_SQ)

!     Shadowing effect (see Tsang, Kong, Shin, Theory of Microwave Remot
!     Sensing, Wiley-Interscience, 1985) 
      IF(DO_SHADOW) THEN
        SHADOW = 1./( SHADOW_ETA(MU_I, SIGMA_SQ, PI) &
     &          + SHADOW_ETA(MU_R, SIGMA_SQ, PI) + 1. )
        BRDF = BRDF*SHADOW
      ENDIF

      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! +--------------------------------------------------------------------
      REAL FUNCTION SHADOW_ETA(COS_THETA, SIGMA_SQ, PI)
! +--------------------------------------------------------------------
! Version 3: shadow effect function
!            called by OCEABRDF2
!   Input:
!
!   COS_THETA     absolute cosine of incident/reflected polar angle (pos
!   SIGMA_SQ      slope variance 
!   PI            3.141592653... constant
!
!   Output:
!
!   SHADOW_ETA:   shadow function
! +--------------------------------------------------------------------
      REAL COS_THETA, SIN_THETA
      REAL MU, SIGMA_SQ, PI
      REAL TERM1, TERM2
      real, parameter :: max_exp=-log(epsilon(max_exp))

      SIN_THETA = SQRT(1.-COS_THETA*COS_THETA)
      if(SIN_THETA.lt.epsilon(SIN_THETA)) then
        SHADOW_ETA = 0
      else
        MU = COS_THETA/SIN_THETA

        TERM1 = SQRT(SIGMA_SQ/PI)/MU*EXP(-min(max_exp,MU*MU/(SIGMA_SQ)))
        TERM2 = ERFC( MU/SQRT(SIGMA_SQ) )

        SHADOW_ETA = 0.5*(TERM1 - TERM2)
      endif
      END
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end module
