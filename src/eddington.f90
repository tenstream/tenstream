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

private
public :: eddington_coeff_zdun

logical,parameter :: ldebug=.False.

      contains

        pure elemental subroutine eddington_coeff_zdun(dtau_in,omega_0_in,g_in,mu_0,c11,c12,c13,c23,c33,g1,g2)
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

          real(ireal128),parameter :: eps_resonance=1e-8_ireal128
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

          alpha1_m_lambda = alpha_1-lambda ! max(epsilon(alpha_1), alpha_1-lambda )
          alpha1_p_lambda = alpha_1+lambda ! max(epsilon(alpha_1), alpha_1+lambda )
          if(approx(alpha1_m_lambda, 0._ireal128)) alpha1_m_lambda = sign(sqrt(tiny(alpha_1)), alpha1_m_lambda)
          if(approx(alpha1_p_lambda, 0._ireal128)) alpha1_p_lambda = sign(sqrt(tiny(alpha_1)), alpha1_p_lambda)

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
            a33     = exp ( - min(max_exponential, dtau / mu_0 ))

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

          !if(ldebug) then
          !  if(      any([c11,c12,c13,c23,c33].gt.one)        &
          !    .or. any([c11,c12,c13,c23,c33,g1,g2].lt.zero) &
          !    .or. any(isnan([c11,c12,c13,c23,c33,g1,g2]))  ) then
          !    print *,'eddington ',dtau_in,omega_0_in,g_in,'::',a11,a12,a13,a23,a33,g1,g2
          !    print *,'eddington ',dtau_in,omega_0_in,g_in,'::',c11,c12,c13,c23,c33,g1,g2
          !    print *,'eddington ',dtau,omega_0,g,mu_0,':1:',A,den,lambda,alpha_1,alpha_2,e1,e2,g0
          !    print *,'eddington ',dtau,omega_0,g,mu_0,':2:',alpha_3,alpha_4,den,(one/mu_0)**2 - lambda**2,alpha_5,alpha_6
          !    print *,'eddington ',dtau,omega_0,g,mu_0,':3:',beta11,beta21,beta12,beta22,beta13,beta23
          !    print *,'eddington ',dtau,omega_0,g,mu_0,':4:',gamma12,gamma22
          !    print *,'eddington ',dtau,omega_0,g,mu_0,':5:',beta13*gamma12,beta23*gamma22,alpha_6*a33
          !    print *,'eddington ',dtau,omega_0,g,mu_0,':6:',alpha1_m_lambda,alpha1_p_lambda
          !    call exit()
          !  endif

          !  if( real(c11+c12) .gt. one  .or. real(c13+c23) .gt. one ) then
          !    print *,'eddington enercons',dtau_in,omega_0_in,g_in,'::',c11+c12,c13+c23,one-c33
          !    print *,'eddington enercons',dtau_in,omega_0_in,g_in,'::',c11,c12,c13,c23,c33,g1,g2
          !    print *,'eddington enercons',dtau,omega_0,g,mu_0,':1:',A,lambda,alpha_1,alpha_2,e1,e2,g0
          !    print *,'eddington enercons',dtau,omega_0,g,mu_0,':2:',alpha_3,alpha_4,(one/mu_0)**2 - lambda**2,alpha_5,alpha_6
          !    print *,'eddington enercons',dtau,omega_0,g,mu_0,':3:',beta11,beta21,beta12,beta22,beta13,beta23
          !    print *,'eddington enercons',dtau,omega_0,g,mu_0,':4:',gamma12,gamma22
          !    print *,'eddington enercons',dtau,omega_0,g,mu_0,':5:',beta13*gamma12,beta23*gamma22,alpha_6*a33
          !    call exit()
          !  endif
          !endif
        end subroutine

      end module
