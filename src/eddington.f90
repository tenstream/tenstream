!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

module m_eddington
      use m_data_parameters, only: ireals,iintegers,zero,one,pi,ireal128
      use m_helper_functions, only: approx,delta_scale_optprop

#ifdef _XLF
  use ieee_arithmetic
#define isnan ieee_is_nan
#endif

      implicit none

      logical,parameter :: ldebug=.False.

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

          if ( lambda * dtau .gt. 100._ireals ) then
            a11 = zero;
            a12 = real(( alpha_1 - lambda ) / alpha_2, ireals)
          else
            exp1  = exp( lambda * dtau );
            term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;

            A = one / ( term1 - one / term1 );

            a11 = real(A * 2 * lambda / alpha_2, ireals)
            a12 = real(A * ( exp1 - one / exp1 ), ireals)
          endif

          den1 = one / ( mu_0_inv * mu_0_inv - lambda * lambda )

          alpha_3 = - omega_0 * b_mmu_0;
          alpha_4 = omega_0 + alpha_3;
          alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1;
          alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1;

          a33     = real(exp ( - dtau  * mu_0_inv ), ireals)

          a13 = real(+ alpha_5 * ( one - a33 * a11 ) - alpha_6 * a12, ireals)
          a23 = real(- alpha_5 * a33 * a12 + alpha_6 * ( a33 - a11 ), ireals)

          a13 = a13 / mu_0 !Fabian: Roberts coefficients a13 expect S to be
          a23 = a23 / mu_0 !        irradiance on tilted plane... we use irradiance on z-plane
      end subroutine
pure subroutine eddington_coeff_bm (dtau_in,omega_0_in,g_in,mu_0,c11,c12,c13,c23,c33)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0
          real(ireals),intent(out) :: c11,c12,c13,c23,c33

          real(ireal128) :: dtau,g,omega0
          real(ireal128) :: a11,a12,a13,a23,a33

          real(ireal128) ::  alpha1, alpha2, alpha3, alpha4, alpha5, alpha6;
          real(ireal128) ::  lambda, A
          real(ireal128) ::  denom, mu0, b

          real(ireal128),parameter ::  eps = sqrt(epsilon(omega0))

          dtau   = max(( epsilon(dtau)  ), real(dtau_in, ireal128))
          g      = max(( epsilon(g)     ), real(g_in, ireal128))
          omega0 = max(( epsilon(omega0)), real(omega_0_in, ireal128))

          omega0 = min(omega0, 1.0_ireal128 - eps)
          mu0 = real(mu_0, ireal128)

          alpha1= (1.0-omega0)+0.75*(1.0-omega0*g);
          alpha2=-(1.0-omega0)+0.75*(1.0-omega0*g);

          lambda=sqrt(alpha1*alpha1-alpha2*alpha2);

          A=1.0/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau));

          a11=A*2.0*lambda/alpha2;
          a12=A*(exp(lambda*dtau)-exp(-lambda*dtau));

          b=0.5-0.75*g*mu0;
          alpha3=-omega0*b;
          alpha4=omega0*(1-b);

          denom = (1.0/mu0/mu0-lambda*lambda);
          alpha5=((alpha1-1.0/mu0)*alpha3-alpha2*alpha4)/denom;
          alpha6=(alpha2*alpha3-(alpha1+1.0/mu0)*alpha4)/denom;

          a33=exp(-dtau/mu0);

          a13=alpha5*(1.0-(a11)*(a33))-alpha6*(a12);
          a23=-(a12)*alpha5*(a33)+alpha6*((a33)-(a11));


          a13 = a13 / mu_0 !Fabian: Roberts coefficients a13 expect S to be
          a23 = a23 / mu_0 !        irradiance on tilted plane... we use irradiance on z-plane

          c11 = real(a11, ireals)
          c12 = real(a12, ireals)
          c13 = real(a13, ireals)
          c23 = real(a23, ireals)
          c33 = real(a33, ireals)
      end subroutine

      subroutine eddington_coeff_cosmo (dtau_in,omega_0_in,g_in,mu_0_in,a11,a12,a13,a23,a33)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0_in
          real(ireals),intent(out) :: a11,a12,a13,a23,a33

          real(ireals)            :: dtau,g,omega_0,mu_0

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
          bscr = 0.5_ireals - 0.375_ireals * g;

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
      subroutine eddington_coeff_fab(dtau_in,omega_0_in,g_in,mu_0,a11,a12,a13,a23,a33,g1,g2)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0
          real(ireals),intent(out) :: a11,a12,a13,a23,a33,g1,g2

          real(ireals)            :: dtau,g,omega_0

          real(ireal128) ::  alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6;
          real(ireal128) ::  b_mmu_0, lambda, A, exp1, term1, bscr, term2, exp2, g0
          real(ireal128) ::  den1, mu_0_inv;

          real(ireals),parameter ::  eps = 10._ireals * epsilon(omega_0)

          ! Singularities -- dont use values before here
          dtau   = max( epsilon(dtau_in)   , dtau_in   )
          g      = max( epsilon(g_in)      , g_in      )
          omega_0= max( epsilon(omega_0_in), omega_0_in)

          omega_0 = min(omega_0, one-eps)
          if ( approx( omega_0 * g , one ) ) omega_0 = omega_0 * (one-eps);
          ! Singularities -- dont use values before here

          mu_0_inv =  1._ireal128/ max(mu_0,epsilon(mu_0));
          b_mmu_0 = (0.5_ireal128 - 0.75_ireal128 * g * mu_0);

          bscr = (0.5_ireal128 - 0.375_ireal128 * g);
          alpha_1 = 2._ireal128 * ( 1._ireal128 - omega_0 * ( 1._ireal128 - bscr ) ) - 0.25_ireal128;
          alpha_2 = 2._ireal128 * omega_0 * bscr - 0.25_ireal128;

          lambda = sqrt ( alpha_1 * alpha_1 - alpha_2 * alpha_2 );

          if ( lambda * dtau .gt. 40._ireal128 ) then
            a11 = zero;
            a12 = real( ( alpha_1 - lambda ) / alpha_2, ireals)
          else
            exp1  = exp( lambda * dtau );
            exp2  = exp(-lambda * dtau );
            term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;
            term2 = alpha_2 / ( alpha_1 + lambda ) * exp2;

            A = min( huge(A), 1.0_ireal128 / ( term1 - term2 ) );
!            A = 1.0_ireal128 / ( term1 - term2 );

            a11 = real(A * 2.0_ireal128 * lambda / alpha_2, ireals)
            a12 = real(A * ( exp1 - exp2 ), ireals)
          endif

          alpha_3 = - omega_0 * b_mmu_0;
          alpha_4 = omega_0 + alpha_3;

          if(mu_0.gt.epsilon(mu_0)) then
            den1 = min( huge(den1), 1._ireal128 / ( mu_0_inv * mu_0_inv - lambda * lambda ))
            alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1;
            alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1;

            a33     = real(exp ( - dtau  * mu_0_inv ), ireals)

            a13 = real(+ alpha_5 * ( 1.0_ireal128 - a33 * a11 ) - alpha_6 * a12, ireals)
            a23 = real(- alpha_5 * a33 * a12 + alpha_6 * ( a33 - a11 ), ireals)

            a13 = a13 / mu_0 !Fabian: Roberts coefficients a13 expect S to be
            a23 = a23 / mu_0 !        irradiance on tilted plane... we use irradiance on z-plane
          else
            a33=zero
            a13=zero
            a23=zero
          endif

          a11 = min(one, max( zero, a11 ) )
          a12 = min(one, max( zero, a12 ) )
          a13 = min(one, max( zero, a13 ) )
          a23 = min(one, max( zero, a23 ) )

          g0 = 2._ireals*(one-omega_0) ! this is alpha3/pi in zdunkowsky for thermal part
          g1 = real(g0 / (alpha_1-alpha_2), ireals)
          g2 = real(g0 / (lambda**2), ireals)
          g1 = min(one, max( zero, g1 ) )
          g2 = min(one, max( zero, g2 ) )

!          print *,'eddington called with',dtau_in,omega_0_in,g_in,'::',a11,a12,a13,a23,a33

!          if(dtau_in.gt.1) then
!            print *,'eddington :: ',dtau_in, omega_0_in,g_in, mu_0
!            print *,'eddington :: ',dtau, omega_0,g, '::',a11,a12,a13,a23,a33
!          endif

      end subroutine
     subroutine eddington_coeff_zdun(dtau_in,omega_0_in,g_in,mu_0,c11,c12,c13,c23,c33,g1,g2)
          real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0
          real(ireals),intent(out) :: c11,c12,c13,c23,c33,g1,g2

          real(ireal128) :: dtau,g,omega_0
          real(ireal128) :: a11,a12,a13,a23,a33

          real(ireal128) ::  alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6
          real(ireal128) ::  beta11,beta21,beta12,beta22,beta13,beta23
          real(ireal128) ::  gamma12,gamma22
          real(ireal128) ::  mubar, bbar, b_minus_mu0
          real(ireal128) ::  lambda, A, den, g0, e1,e2
          real(ireal128) ::  alpha1_p_lambda, alpha1_m_lambda
          real(ireal128) ::  bscr

          real(ireal128),parameter ::  eps_resonance=1e-8_ireal128
          real(ireal128),parameter :: max_exponential = log(huge(max_exponential)/1e2)

          ! Singularities -- dont use values before here
          dtau   = max( epsilon(dtau_in)   , dtau_in   )
          dtau   = min( 500._ireal128      , dtau      )
          g      = max( 1e-6_ireals        , g_in      )
          omega_0= max( epsilon(omega_0_in), omega_0_in)

          omega_0 = min(omega_0, one-eps_resonance)
          ! Singularities -- dont use values before here


          mubar = .5_ireal128
          bbar  = 3._ireal128/8._ireal128*(one-g)
          b_minus_mu0 = .5_ireal128 - .75_ireal128 * g *mu_0

          alpha_1 = ( one - omega_0*(one-bbar) ) / mubar
          alpha_2 = omega_0*bbar/mubar

          bscr = 0.5_ireal128 - 0.375_ireal128 * g;
          alpha_1 = 2._ireal128 * ( 1._ireal128 - omega_0 * ( 1._ireal128 - bscr ) ) - 0.25_ireal128;
          alpha_2 = 2._ireal128 * omega_0 * bscr - 0.25_ireal128;

          lambda = sqrt(alpha_1**2 - alpha_2** 2)

          e1 = exp( min(max_exponential, lambda*dtau))
          e2 = exp(-min(max_exponential, lambda*dtau))

          alpha1_m_lambda = max(epsilon(alpha_1), alpha_1-lambda )
          alpha1_p_lambda = max(epsilon(alpha_1), alpha_1+lambda )

          A = one / ( alpha_2/alpha1_m_lambda*e1 - alpha_2/alpha1_p_lambda * e2 )

          beta11  =  A * alpha_2/alpha1_m_lambda
          beta21  = -A * alpha_2/alpha1_p_lambda
          beta12  = -A * e2
          beta22  =  A * e1

!          gamma11 = one
          gamma12 = alpha_2/alpha1_p_lambda * e1
!          gamma21 = one
          gamma22 = alpha_2/alpha1_m_lambda * e2

          a11 = beta11 + beta21
          a12 = beta12 + beta22

          a11 = max(0._ireal128,  min(1._ireal128, a11) )
          a12 = max(0._ireal128,  min(1._ireal128, a12) )

          if(mu_0.gt.epsilon(mu_0)) then
            a33     = exp ( - dtau  / mu_0 )

            alpha_3 = -omega_0 * b_minus_mu0
            alpha_4 =  omega_0 * (one-b_minus_mu0)

            den = (one/mu_0)**2 - lambda**2
            if( abs(den).le.eps_resonance ) then ! den.ge.-epsilon(den) .and. den.le.epsilon(den)  ) then !avoid resonance case
              if(mu_0.gt..5_ireal128) then
                den = one/ (mu_0**2 - eps_resonance)  - lambda**2
              else
                den = one/ (mu_0**2 + eps_resonance)  - lambda**2
              endif
            endif

            alpha_5 = ( (alpha_1-one/mu_0)*alpha_3 - alpha_2*alpha_4 ) / den
            alpha_6 = ( alpha_2*alpha_3 - (alpha_1+one/mu_0)*alpha_4 ) / den

            beta13  = -beta11*alpha_5 * a33 - beta12*alpha_6
            beta23  = -beta21*alpha_5 * a33 - beta22*alpha_6

            a13 = beta13         + beta23         + alpha_5
            a23 = beta13*gamma12 + beta23*gamma22 + alpha_6*a33

            a13 = a13 / mu_0 !Fabian: Roberts coefficients a13 expect S to be
            a23 = a23 / mu_0 !        irradiance on tilted plane... we use irradiance on z-plane

            a13 = max(0._ireal128, a13)
            a23 = max(0._ireal128, a23)

          else

            a33=zero
            a13=zero
            a23=zero
          endif

          g0 = (one-omega_0)/mubar ! this is alpha3/pi in zdunkowsky for thermal part
          g1 = real(g0 / (alpha_1-alpha_2), ireals)
          g2 = real(g0 / lambda**2, ireals)

          c11 = real(a11, ireals)
          c12 = real(a12, ireals)
          c13 = real(a13, ireals)
          c23 = real(a23, ireals)
          c33 = real(a33, ireals)

          if(ldebug) then
            if(      any([c11,c12,c13,c23,c33].gt.one)        &
                .or. any([c11,c12,c13,c23,c33,g1,g2].lt.zero) &
                .or. any(isnan([c11,c12,c13,c23,c33,g1,g2]))  ) then
            print *,'eddington ',dtau_in,omega_0_in,g_in,'::',a11,a12,a13,a23,a33,g1,g2
            print *,'eddington ',dtau_in,omega_0_in,g_in,'::',c11,c12,c13,c23,c33,g1,g2
            print *,'eddington ',dtau,omega_0,g,mu_0,':1:',A,den,lambda,alpha_1,alpha_2,e1,e2,g0
            print *,'eddington ',dtau,omega_0,g,mu_0,':2:',alpha_3,alpha_4,den,(one/mu_0)**2 - lambda**2,alpha_5,alpha_6
            print *,'eddington ',dtau,omega_0,g,mu_0,':3:',beta11,beta21,beta12,beta22,beta13,beta23
            print *,'eddington ',dtau,omega_0,g,mu_0,':4:',gamma12,gamma22
            print *,'eddington ',dtau,omega_0,g,mu_0,':5:',beta13*gamma12,beta23*gamma22,alpha_6*a33
            print *,'eddington ',dtau,omega_0,g,mu_0,':6:',alpha1_m_lambda,alpha1_p_lambda
            call exit()
          endif

          if( real(c11+c12) .gt. one  .or. real(c13+c23) .gt. one ) then
            print *,'eddington enercons',dtau_in,omega_0_in,g_in,'::',c11+c12,c13+c23,one-c33
            print *,'eddington enercons',dtau_in,omega_0_in,g_in,'::',c11,c12,c13,c23,c33,g1,g2
            print *,'eddington enercons',dtau,omega_0,g,mu_0,':1:',A,lambda,alpha_1,alpha_2,e1,e2,g0
            print *,'eddington enercons',dtau,omega_0,g,mu_0,':2:',alpha_3,alpha_4,(one/mu_0)**2 - lambda**2,alpha_5,alpha_6
            print *,'eddington enercons',dtau,omega_0,g,mu_0,':3:',beta11,beta21,beta12,beta22,beta13,beta23
            print *,'eddington enercons',dtau,omega_0,g,mu_0,':4:',gamma12,gamma22
            print *,'eddington enercons',dtau,omega_0,g,mu_0,':5:',beta13*gamma12,beta23*gamma22,alpha_6*a33
            call exit()
          endif
        endif

      end subroutine

    end module
