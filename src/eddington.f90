module m_eddington
      use m_data_parameters, only: ireals,iintegers,mpiint,imp_comm,zero,one,pi,ireal128
      use m_helper_functions, only: approx,delta_scale_optprop
      implicit none

      private
      public :: eddington_coeff_rb,eddington_coeff_fab,eddington_coeff_cosmo

      integer(mpiint) :: VALUE_ERROR=-5

!      logical,parameter :: ldelta_scale=.True.

      contains

pure subroutine eddington_coeff_rb (dtau_in,omega_0_in,g_in,mu_0,a11,a12,a13,a23,a33)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0
          real(ireals),intent(out) :: a11,a12,a13,a23,a33

          real(ireals)            :: dtau,g,omega_0

          real(ireal128) ::  alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6;
          real(ireal128) ::  b_mmu_0, lambda, A, exp1, term1, bscr;
          real(ireal128) ::  den1, mu_0_inv;

          real(ireals),parameter ::  eps = sqrt(epsilon(omega_0))

          dtau   = max(( epsilon(dtau)   ), dtau_in)
          g      = max(( epsilon(g)      ), g_in)
          omega_0= max(( epsilon(omega_0)), omega_0_in)

          omega_0 = min(omega_0, 1.0_ireals - eps)
          if ( approx( omega_0 * g , 1.0_ireals ) ) omega_0 = omega_0 * (1.0_ireals - eps);

          call delta_scale_optprop( dtau, omega_0, g  )
!          print *,'eddington_coeff_rb dtau_in',dtau_in,'dtau',dtau

          mu_0_inv = 1._ireal128/mu_0;
    
          b_mmu_0 = 0.5_ireal128 - 0.75_ireal128 * g * mu_0;

          bscr = 0.5_ireal128 - 0.375_ireal128 * g;
          alpha_1 = 2._ireal128 * ( 1._ireal128 - omega_0 * ( 1._ireal128 - bscr ) ) - 0.25_ireal128;
          alpha_2 = 2._ireal128 * omega_0 * bscr - 0.25_ireal128;

          lambda = sqrt ( alpha_1 * alpha_1 - alpha_2 * alpha_2 );

          if ( lambda * dtau .gt. 100._ireal128 ) then
            a11 = zero;
            a12 = ( alpha_1 - lambda ) / alpha_2;
          else
            exp1  = exp( lambda * dtau );
            term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;

            A = 1.0_ireal128 / ( term1 - 1._ireal128 / term1 );

            a11 = A * 2.0_ireal128 * lambda / alpha_2;
            a12 = A * ( exp1 - 1._ireal128 / exp1 );
          endif

          den1 = 1._ireal128 / ( mu_0_inv * mu_0_inv - lambda * lambda )

          alpha_3 = - omega_0 * b_mmu_0;
          alpha_4 = omega_0 + alpha_3;
          alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1;
          alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1;

          a33     = exp ( - dtau  * mu_0_inv );   

          a13 = + alpha_5 * ( 1.0_ireal128 - a33 * a11 ) - alpha_6 * a12;
          a23 = - alpha_5 * a33 * a12 + alpha_6 * ( a33 - a11 );

          a13 = a13 * mu_0_inv !Fabian: Roberts coefficients a13 expect S to be
          a23 = a23 * mu_0_inv !        irradiance on tilted plane... we use irradiance on z-plane

!          a12 = max( zero, a12 )!min(one, )
!          a13 = max( zero, a13 )!min(one, )
!          a23 = max( zero, a23 )!min(one, )



      end subroutine

      subroutine eddington_coeff_cosmo (dtau_in,omega_0_in,g_in,mu_0_in,a11,a12,a13,a23,a33)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0_in
          real(ireals),intent(out) :: a11,a12,a13,a23,a33

          real(ireals)            :: dtau,g,omega_0,mu_0

!          real(ireals) ::  alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6;
          real(ireals) ::  bscr,ze,zeps,zm,ze1mwf,mu_0_inv,zg1,zg2,b_mmu_0,zmu0if
          real(ireals) ::  zod1,zod2,zod3,zod4,zod5

          real(ireals),parameter ::  eps = sqrt(epsilon(omega_0))
          REAL    (KIND=ireals   ), PARAMETER ::  &
              zargli  = 80.0     , &  ! argument limit for EXP  
              ztsec   = 1.0E-35  , &  ! (=exp(-zargli) avoids ALOG(0.0)
              zepres  = 1.0E-7  ! , &  ! for resonance case avoidance 
!              zangfa  = 1.648721271   ! exp(0.5)


          dtau   = max(( epsilon(dtau)   ), dtau_in)
          g      = max(( epsilon(g)      ), g_in)
          omega_0= max(( epsilon(omega_0)), omega_0_in)
          mu_0   = max(eps,min(one-eps,mu_0_in))

          omega_0 = min(omega_0, 1.0_ireals - eps)
          if ( approx( omega_0 * g , 1.0_ireals ) ) omega_0 = omega_0 * (1.0_ireals - eps);

          call delta_scale_optprop( dtau, omega_0, g  )

          
          mu_0_inv =  1._ireals/ mu_0;
          b_mmu_0 = 0.5_ireals - 0.75_ireals * g * mu_0;
          bscr = 0.5_ireal128 - 0.375_ireal128 * g;

          zod1 = 2._ireals * ( 1._ireals - omega_0 * ( 1._ireals - bscr ) ) - 0.25_ireals;
          zod2 = 2._ireals * omega_0 * bscr - 0.25_ireals;

          zod3 = - omega_0 * b_mmu_0;
          zod4 = omega_0 + zod3;

          zod5 = dtau

          zeps=SQRT(zod1*zod1-zod2*zod2)
          IF (zeps.LT.zargli) THEN
            ze = EXP  (-zeps)
          ELSE
            ze = ztsec
          END IF
          zm = zod2/(zod1+zeps)

          a11=ze*(1.-(zm**2))*(1./(1.-(zm**2)*(ze**2)))
          a12=zm*(1.-(ze**2))*(1./(1.-(zm**2)*(ze**2)))

          ze1mwf = zeps / zod5
          zmu0if = ze1mwf + SIGN ( MAX(ABS(mu_0_inv-ze1mwf),zepres) &
              ,(mu_0_inv-ze1mwf) )

          zod3 = zod3 * zmu0if
          zod4 = zod4 * zmu0if
          zod5 = zod5 * zmu0if

          IF (zod5.LT.ZARGLI) THEN
            a33 = EXP  (-zod5)
          ELSE
            a33 = ztsec
          END IF

          zg1 = ( zod3*(zod5-zod1) -zod2*zod4) /(zod5*zod5 - zeps*zeps)
          zg2 =-( zod4*(zod5+zod1) +zod2*zod3) /(zod5*zod5 - zeps*zeps)
          a23 = zg2*(a33-a11) -zg1*a12*a33
          a13 = zg1*(1.-a11*a33) -zg2*a12

      end subroutine
      subroutine eddington_coeff_fab (dtau_in,omega_0_in,g_in,mu_0,a11,a12,a13,a23,a33,g1,g2)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0
          real(ireals),intent(out) :: a11,a12,a13,a23,a33,g1,g2

          real(ireals)            :: dtau,g,omega_0

          real(ireal128) ::  alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6;
          real(ireal128) ::  b_mmu_0, lambda, A, exp1, term1, bscr, term2, exp2, g0
          real(ireal128) ::  den1, mu_0_inv;

          real(ireals),parameter ::  eps = 10._ireals * epsilon(omega_0)

          ! Singularities -- dont use values before here
          dtau   = max(( epsilon(dtau)   ), dtau_in)
          g      = g_in
          omega_0= max(( epsilon(omega_0)), omega_0_in)

          omega_0 = min(omega_0, one-eps)
          if ( approx( omega_0 * g , 1.0_ireals ) ) omega_0 = omega_0 * (one-eps);
          ! Singularities -- dont use values before here

          mu_0_inv =  1._ireal128/ mu_0;
          b_mmu_0 = (0.5_ireal128 - 0.75_ireal128 * g * mu_0);

          bscr = (0.5_ireal128 - 0.375_ireal128 * g);
          alpha_1 = 2._ireal128 * ( 1._ireal128 - omega_0 * ( 1._ireal128 - bscr ) ) - 0.25_ireal128;
          alpha_2 = 2._ireal128 * omega_0 * bscr - 0.25_ireal128;

          lambda = sqrt ( alpha_1 * alpha_1 - alpha_2 * alpha_2 );

          if ( lambda * dtau .gt. 40._ireal128 ) then
            a11 = zero;
            a12 = ( alpha_1 - lambda ) / alpha_2;
          else
            exp1  = exp( lambda * dtau );
            exp2  = exp(-lambda * dtau );
            term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;
            term2 = alpha_2 / ( alpha_1 + lambda ) * exp2;

            A = min( huge(A), 1.0_ireal128 / ( term1 - term2 ) );
!            A = 1.0_ireal128 / ( term1 - term2 );

            a11 = A * 2.0_ireal128 * lambda / alpha_2;
            a12 = A * ( exp1 - exp2 );
          endif

          den1 = min( huge(den1), 1._ireal128 / ( mu_0_inv * mu_0_inv - lambda * lambda ))

          alpha_3 = - omega_0 * b_mmu_0;
          alpha_4 = omega_0 + alpha_3;
          alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1;
          alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1;

          a33     = exp ( - dtau  * mu_0_inv );   

          a13 = + alpha_5 * ( 1.0_ireal128 - a33 * a11 ) - alpha_6 * a12;
          a23 = - alpha_5 * a33 * a12 + alpha_6 * ( a33 - a11 );

          a13 = a13 * mu_0_inv !Fabian: Roberts coefficients a13 expect S to be
          a23 = a23 * mu_0_inv !        irradiance on tilted plane... we use irradiance on z-plane

          a11 = min(one, max( zero, a11 ) )
          a12 = min(one, max( zero, a12 ) )
          a13 = min(one, max( zero, a13 ) )
          a23 = min(one, max( zero, a23 ) )

          g0 = 2._ireals*(one-omega_0) ! this is alpha3/pi in zdunkowsky for thermal part
          g1 = g0 / (alpha_1-alpha_2)
          g2 = g0 / (lambda**2)
          g1 = min(one, max( zero, g1 ) )
          g2 = min(one, max( zero, g2 ) )

!          print *,'eddington called with',dtau_in,omega_0_in,g_in,'::',a11,a12,a13,a23,a33

!          if(dtau_in.gt.1) then
!            print *,'eddington :: ',dtau_in, omega_0_in,g_in, mu_0
!            print *,'eddington :: ',dtau, omega_0,g, '::',a11,a12,a13,a23,a33
!          endif

      end subroutine

    end module
