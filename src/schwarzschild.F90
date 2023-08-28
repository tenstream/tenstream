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

module m_schwarzschild

#ifdef _XLF
  use ieee_arithmetic
#define isnan ieee_is_nan
#endif

  use m_data_parameters, only: ireals, iintegers, zero, one, pi, EXP_MINVAL
  use m_helper_functions, only: get_arg, expm1
  implicit none

  private
  public schwarzschild, B_eff

contains

  subroutine B_eff(B_far, B_near, tau, B)
    real(ireals), intent(in) :: B_far, B_near, tau
    real(ireals), intent(out) :: B
    integer(iintegers) :: imu
    real(ireals) :: mu
    integer(iintegers), parameter :: Nmu = 2
    real(ireals), save :: legendre_wi(Nmu) = -1._ireals
    real(ireals), save :: legendre_pt(Nmu) = -1._ireals

    if (legendre_wi(1) .lt. 0._ireals) call dgauss(size(legendre_wi), legendre_pt, legendre_wi)

    B = 0
    do imu = 1, size(legendre_pt)
      mu = legendre_pt(imu)
      B = B + B_eff_mu(B_far, B_near, tau, mu) * mu * legendre_wi(imu)
    end do
    B = B * 2

  contains
    real(ireals) function B_eff_mu(B_far, B_near, tau, mu)
      real(ireals), intent(in) :: B_far, B_near, tau, mu
      real(ireals) :: tm1, dtau
      real(ireals), parameter :: eps = 1e-3_ireals
      dtau = tau / mu
      if (dtau .lt. eps) then
        B_eff_mu = (B_far + B_near)*.5_ireals
      else
        tm1 = expm1(-dtau)
        B_eff_mu = (-B_near + B_far * (tm1 + 1)) / (tm1) + ((B_far - B_near) * mu) / tau
      end if
    end function
  end subroutine

  subroutine schwarzschild_radiance(tau, B_near, B_far, L)
    real(ireals), intent(in) :: B_far, B_near, tau
    real(ireals), intent(inout) :: L
    real(ireals) :: tm1
    if (tau .gt. 1e-3) then
      tm1 = expm1(-tau)
      L = L * (tm1 + 1) + (B_far - B_near) - (B_near - (B_far - B_near) / tau) * tm1
    else
      L = (B_near + B_far)*.5_ireals * tau + L * (1._ireals - tau)
    end if
  end subroutine

  subroutine schwarzschild(Nmu, dtau, albedo, Edn, Eup, planck, opt_srfc_emission)
    integer(iintegers), intent(in) :: Nmu
    real(ireals), intent(in), dimension(:) :: dtau
    real(ireals), intent(in) :: albedo
    real(ireals), dimension(:), intent(out) :: Edn, Eup
    real(ireals), dimension(:), intent(in) :: planck
    ! planck value specific for surface. i.e. overrides last planck value
    real(ireals), intent(in), optional :: opt_srfc_emission

    integer(iintegers) :: imu, k, ke, ke1

    real(ireals) :: Lup, Ldn, Bsrfc
    real(ireals) :: dmu, mu
    real(ireals) :: legendre_wi(Nmu)
    real(ireals) :: legendre_pt(Nmu)
    logical, parameter :: use_legendre = .true.

    Edn = zero
    Eup = zero

    ke = size(dtau)
    ke1 = ke + 1

    Bsrfc = get_arg(planck(ke1), opt_srfc_emission)

    if (use_legendre) then
      call dgauss(size(legendre_wi), legendre_pt, legendre_wi)
      do imu = 1, size(legendre_pt)
        mu = legendre_pt(imu)

        ! zero incoming radiation at TOA
        Ldn = zero
        !Edn(1) = Edn(1) + Ldn * mu * legendre_wi(imu)

        do k = 1, ke
          call schwarzschild_radiance(dtau(k) / mu, planck(k), planck(k + 1), Ldn)
          Edn(k + 1) = Edn(k + 1) + Ldn * mu * legendre_wi(imu)
        end do
      end do

      do imu = 1, size(legendre_pt)
        mu = legendre_pt(imu)

        ! Boundary conditions at surface
        Lup = Bsrfc * (one - albedo) + albedo * Edn(ke1) * 2
        Eup(ke1) = Eup(ke1) + Lup * mu * legendre_wi(imu)

        do k = ke, 1, -1
          call schwarzschild_radiance(dtau(k) / mu, planck(k + 1), planck(k), Lup)
          Eup(k) = Eup(k) + Lup * mu * legendre_wi(imu)
        end do
      end do ! enddo mu

      Eup = Eup * 2 * pi
      Edn = Edn * 2 * pi

    else

      dmu = one / real(Nmu, ireals)
      ! Transmission coefficients
      do imu = 1, Nmu
        mu = (real(imu, ireals) - .5_ireals) * dmu

        ! zero incoming radiation at TOA
        Ldn = zero
        !Edn(1) = Edn(1) + Ldn * mu

        do k = 1, ke
          call schwarzschild_radiance(dtau(k) / mu, planck(k), planck(k + 1), Ldn)
          Edn(k + 1) = Edn(k + 1) + Ldn * mu
        end do
      end do

      do imu = 1, Nmu
        mu = (real(imu, ireals) - .5_ireals) * dmu

        ! Boundary conditions at surface
        Lup = Bsrfc * (one - albedo) + albedo * Edn(ke1) * 2 * dmu
        Eup(ke1) = Eup(ke1) + Lup * mu

        do k = ke, 1, -1
          call schwarzschild_radiance(dtau(k) / mu, planck(k + 1), planck(k), Lup)
          Eup(k) = Eup(k) + Lup * mu
        end do

      end do

      Eup = Eup * 2 * pi * dmu
      Edn = Edn * 2 * pi * dmu
    end if
  end subroutine

  subroutine dgauss(m, gmu, gwt)
    !                                                                         *
    !     compute weights and abscissae for ordinary gaussian quadrature      *
    !     (no weight function inside integral) on the interval (0,1)          *
    !                                                                         *
    !     input :    m                      order of quadrature rule          *
    !                                                                         *
    !     output :   gmu(i)  i = 1 to m     array of abscissae                *
    !                gwt(i)  i = 1 to m     array of weights                  *
    !                                                                         *
    !     reference:  davis, p.j. and p. rabinowitz,                          *
    !                 methods of numerical integration,                       *
    !                 academic press, new york, pp. 87, 1975.                 *
    !                                                                         *
    !     method:  compute the abscissae as roots of the legendre             *
    !              polynomial p-sub-m using a cubically convergent            *
    !              refinement of newton's method.  compute the                *
    !              weights from eq. 2.7.3.8 of davis/rabinowitz.  note        *
    !              that newton's method can very easily diverge; only a       *
    !              very good initial guess can guarantee convergence.         *
    !              the initial guess used here has never led to divergence    *
    !              even for m up to 1000.                                     *
    !                                                                         *
    !     accuracy:  at least 13 significant digits                           *
    !                                                                         *
    !     internal variables:                                                 *
    !                                                                         *
    !    iter      : number of newton method iterations                       *
    !    maxit     : maximum allowed iterations of newton method              *
    !    pm2,pm1,p : 3 successive legendre polynomials                        *
    !    ppr       : derivative of legendre polynomial                        *
    !    p2pri     : 2nd derivative of legendre polynomial                    *
    !    tol       : convergence criterion for legendre poly root iteration   *
    !    x,xi      : successive iterates in cubically-convergent version      *
    !                of newtons method (seeking roots of legendre poly.)      *
    !**************************************************************************

    !implicit double precision ( a-h, o-z )
    real(ireals), dimension(:) :: gmu, gwt
    integer :: iter, lim, m, np1, k, nn

    double precision :: cona, en, nnp1, p, pm1, pm2, ppr, p2pri, prod, &
      tmp, x, xi, t

    integer, parameter :: maxit = 1000
    double precision, parameter :: one = 1, two = 2, eps = 1.d-14

    if (m .lt. 1) then
      write (*, *) 'dgauss--bad value of m'
      stop
    end if

    if (m .eq. 1) then
      gmu(1) = 0.5
      gwt(1) = 1.0
      return
    end if

    en = m
    np1 = m + 1
    nnp1 = m * np1
    cona = dble(m - 1) / (8 * m**3)

    lim = m / 2
    do k = 1, lim

      !        initial guess for k-th root of legendre polynomial,
      !        from davis/rabinowitz (2.7.3.3a)

      t = (4 * k - 1) * pi / (4 * m + 2)
      x = cos(t + cona / tan(t))
      iter = 0

      !        upward recurrence for legendre polynomials

10    iter = iter + 1
      pm2 = one
      pm1 = x
      do nn = 2, m
        p = ((2 * nn - 1) * x * pm1 - (nn - 1) * pm2) / nn
        pm2 = pm1
        pm1 = p
      end do

      !        newton method

      tmp = one / (one - x**2)
      ppr = en * (pm2 - x * p) * tmp
      p2pri = (two * x * ppr - nnp1 * p) * tmp
      xi = x - (p / ppr) * (one + (p / ppr) * p2pri / (two * ppr))

      !        check for convergence

      if (abs(xi - x) .gt. eps) then
        if (iter .gt. maxit) then
          write (*, *) 'dgauss--max iteration count'
          stop
        end if
        x = xi
        go to 10
      end if

      !        iteration finished--calculate weights, abscissae for (-1,1)

      gmu(k) = -real(x, kind(gmu))
      gwt(k) = real(two / (tmp * (en * pm2)**2), kind(gmu))
      gmu(np1 - k) = -gmu(k)
      gwt(np1 - k) = gwt(k)

    end do

    !     set middle abscissa and weight for rules of odd order

    if (mod(m, 2) .ne. 0) then
      gmu(lim + 1) = 0.0
      prod = one
      do k = 3, m, 2
        prod = prod * k / (k - 1)
      end do
      gwt(lim + 1) = real(two / prod**2, kind(gmu))
    end if

    !     convert from (-1,1) to (0,1) and resort

    do k = 1, m
      gmu(k) = 0.5 * gmu(k) + 0.5
      gwt(k) = 0.5 * gwt(k)
    end do

    return
  end subroutine dgauss
end module
