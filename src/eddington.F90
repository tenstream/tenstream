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
  use m_data_parameters, only: ireals,iintegers,zero,one,pi,irealeddington
  use m_helper_functions, only: approx, delta_scale_optprop

#ifdef _XLF
  use ieee_arithmetic
#define isnan ieee_is_nan
#endif

implicit none

private
public :: eddington_coeff_zdun, eddington_coeff_bm, eddington_coeff_ec

logical,parameter :: ldebug=.True.

contains

  pure elemental subroutine eddington_coeff_zdun(dtau_in,omega_0_in,g_in,mu_0,c11,c12,c13,c23,c33)
    real(ireals),intent(in) :: dtau_in,g_in,omega_0_in,mu_0
    real(ireals),intent(out) :: c11,c12,c13,c23,c33

    real(irealeddington) :: dtau,g,omega_0
    real(irealeddington) :: a11,a12,a13,a23,a33

    real(irealeddington) ::  alpha_1, alpha_2, alpha_3, alpha_4, alpha_5, alpha_6
    real(irealeddington) ::  beta11,beta21,beta12,beta22,beta13,beta23
    real(irealeddington) ::  gamma12,gamma22
    real(irealeddington) ::  mubar, bbar, b_minus_mu0
    real(irealeddington) ::  lambda, A, den, g0, e1,e2
    real(irealeddington) ::  alpha1_p_lambda, alpha1_m_lambda
    real(irealeddington) ::  bscr

    real(irealeddington),parameter :: eps_resonance=1e-8_irealeddington
    real(irealeddington),parameter :: max_exponential = log(huge(max_exponential)/1e2)

    ! Singularities -- dont use values before here
    dtau   = max( epsilon(dtau_in)   , dtau_in   )
    dtau   = min( 500._irealeddington      , dtau      )
    g      = max( 1e-6_ireals        , g_in      )
    omega_0= max( epsilon(omega_0_in), omega_0_in)

    omega_0 = min(omega_0, one-eps_resonance)
    ! Singularities -- dont use values before here

    mubar = .5_irealeddington
    bbar  = 3._irealeddington/8._irealeddington*(one-g)
    b_minus_mu0 = .5_irealeddington - .75_irealeddington * g *mu_0

    alpha_1 = ( one - omega_0*(one-bbar) ) / mubar
    alpha_2 = omega_0*bbar/mubar

    bscr = 0.5_irealeddington - 0.375_irealeddington * g;
    alpha_1 = 2._irealeddington * ( 1._irealeddington - omega_0 * ( 1._irealeddington - bscr ) ) - 0.25_irealeddington;
    alpha_2 = 2._irealeddington * omega_0 * bscr - 0.25_irealeddington;

    lambda = sqrt(alpha_1**2 - alpha_2** 2)

    e1 = exp( min(max_exponential, lambda*dtau))
    e2 = exp(-min(max_exponential, lambda*dtau))

    alpha1_m_lambda = alpha_1-lambda ! max(epsilon(alpha_1), alpha_1-lambda )
    alpha1_p_lambda = alpha_1+lambda ! max(epsilon(alpha_1), alpha_1+lambda )
    if(approx(alpha1_m_lambda, 0._irealeddington)) alpha1_m_lambda = sign(epsilon(alpha_1), alpha1_m_lambda)
    if(approx(alpha1_p_lambda, 0._irealeddington)) alpha1_p_lambda = sign(epsilon(alpha_1), alpha1_p_lambda)

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

    a11 = max(0._irealeddington,  min(1._irealeddington, a11) )
    a12 = max(0._irealeddington,  min(1._irealeddington, a12) )

    if(mu_0.gt.epsilon(mu_0)) then
      a33     = exp ( - min(max_exponential, dtau / mu_0 ))

      alpha_3 = -omega_0 * b_minus_mu0
      alpha_4 =  omega_0 * (one-b_minus_mu0)

      den = (one/mu_0)**2 - lambda**2
      if( abs(den).le.eps_resonance ) then ! den.ge.-epsilon(den) .and. den.le.epsilon(den)  ) then !avoid resonance case
        if(mu_0.gt..5_irealeddington) then
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

      a13 = max(0._irealeddington, a13)
      a23 = max(0._irealeddington, a23)

    else

      a33=zero
      a13=zero
      a23=zero
    endif

    g0 = (one-omega_0)/mubar ! this is alpha3/pi in zdunkowsky for thermal part
    c11 = real(a11, ireals)
    c12 = real(a12, ireals)
    c13 = real(a13, ireals)
    c23 = real(a23, ireals)
    c33 = real(a33, ireals)
  end subroutine

  pure elemental subroutine eddington_coeff_bm(&
      & dtau, omega0, g, mu0, &
      & t, r, rdir, sdir, tdir)
    real(ireals),intent(in) :: dtau, g, omega0, mu0
    real(ireals),intent(out) :: t, r, rdir, sdir, tdir
    real(ireals) :: alpha1, alpha2, alpha3, alpha4, alpha5, alpha6
    real(ireals) :: a11, a12, a13, a23, a33
    real(ireals) :: lambda, b, A, denom

    alpha1= (one-omega0)+0.75_ireals*(one-omega0*g)
    alpha2=-(one-omega0)+0.75_ireals*(one-omega0*g)

    lambda=sqrt(alpha1*alpha1-alpha2*alpha2)

    A=one/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau))

    a11=A*2.0_ireals*lambda/alpha2
    a12=A*(exp(lambda*dtau)-exp(-lambda*dtau))

    b=0.5-0.75*g*mu0
    alpha3=-omega0*b
    alpha4=omega0*(one-b)
    denom = (one/mu0/mu0-lambda*lambda)
    alpha5=((alpha1-one/mu0)*alpha3-alpha2*alpha4)/denom
    alpha6=(alpha2*alpha3-(alpha1+one/mu0)*alpha4)/denom

    a33=exp(-dtau/mu0)

    a13=alpha5*(one-(a11)*(a33))-alpha6*(a12)
    a23=-(a12)*alpha5*(a33)+alpha6*((a33)-(a11))

    t    = a11
    r    = a12
    tdir = a33
    rdir = a13 / mu0;
    sdir = a23 / mu0;
  end subroutine

  pure elemental subroutine eddington_coeff_ec(&
      & dtau, w0, g, mu0, &
      & t, r, rdir, sdir, tdir)

    real(ireals),intent(in) :: dtau, g, w0, mu0
    real(ireals),intent(out) :: t, r, rdir, sdir, tdir
    real(irealeddington) :: f, g1, g2, g3
    real(irealeddington) :: g4, alpha1, alpha2, A, beta
    real(irealeddington) :: e0, e, e2
    real(irealeddington) :: k_mu0, k_g3, k_g4
    real(irealeddington) :: k_2_e, dtau_slant

    f = 0.75_irealeddington * g
    g1 = 2._irealeddington  - w0 * (1.25_irealeddington + f)
    g2 = w0 * (0.75_irealeddington - f)
    g3 = 0.5_irealeddington  - mu0 * f

    dtau_slant = real(max(dtau / max(sqrt(tiny(mu0)), mu0), 0._ireals), irealeddington)

    if(dtau_slant.gt.1e-6_ireals) then
      g4 = 1._irealeddington - g3
      alpha1 = g1*g4 + g2*g3
      alpha2 = g1*g3 + g2*g4

      A = sqrt(max((g1 - g2) * (g1 + g2), 1e-12_irealeddington))
      k_mu0 = A*mu0
      k_g3 = A*g3
      k_g4 = A*g4

      e0 = exp(-dtau_slant)
      tdir = real(e0, ireals)

      e = exp(-A*dtau)
      e2 = e*e
      k_2_e = 2 * A * e
      if (approx(k_mu0,1._irealeddington)) then
        k_mu0 = 1 - 10*epsilon(k_mu0)
      end if
      beta = 1 / (A + g1 + (A - g1)*e2)
      r = real(g2 * (1 - e2) * beta, ireals)
      t = real(k_2_e * beta, ireals)

      beta = w0 * beta / (1 - k_mu0*k_mu0)

      sdir = real(beta * ( k_2_e*(g4 + alpha1*mu0) - e0 &
        & * ( (1 + k_mu0) * (alpha1 + k_g4) &
        &    -(1 - k_mu0) * (alpha1 - k_g4) * e2) ), ireals)
      rdir = real(beta &
        &  * ( (1 - k_mu0) * (alpha2 + k_g3) &
        &     -(1 + k_mu0) * (alpha2 - k_g3)*e2 &
        &     -k_2_e*(g3 - alpha2*mu0)*e0), ireals)
    else
      t    = real(1._irealeddington - g1 * real(dtau, irealeddington), ireals)
      r    = real(g2 * real(dtau, irealeddington), ireals)
      sdir = real((1._irealeddington - g3) * real(w0 * dtau, irealeddington), ireals)
      rdir = real(g3 * real(w0 * dtau, irealeddington), ireals)
      tdir = real(1._irealeddington - dtau_slant, ireals)
    endif
  end subroutine
end module
