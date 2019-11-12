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

use m_data_parameters, only: ireals,iintegers,zero,one,pi,EXP_MINVAL
use m_helper_functions, only: get_arg
implicit none

private
public schwarzschild, B_eff

    contains

      subroutine B_eff(B_far, B_near, tau, B)
        real(ireals), intent(in) :: B_far, B_near, tau
        real(ireals), intent(out) :: B
        integer(iintegers) :: imu
        real(ireals) :: mu
        integer(iintegers), parameter :: Nmu=5
        real(ireals), save :: legendre_wi(Nmu)=-1._ireals
        real(ireals), save :: legendre_pt(Nmu)=-1._ireals

        if(legendre_wi(1).lt.0._ireals) call dgauss(size(legendre_wi), legendre_pt, legendre_wi)

        B=0
        do imu=1,size(legendre_pt)
          mu = legendre_pt(imu)
          B = B + B_eff_mu(B_far, B_near, tau, mu) * mu * legendre_wi(imu)
        enddo
        B = B * 2

        contains
          real(ireals) function B_eff_mu(B_far, B_near, tau, mu)
            real(ireals), intent(in) :: B_far, B_near, tau, mu
            real(ireals) :: t
            if(tau/mu.lt.sqrt(epsilon(tau))) then
              B_eff_mu = (B_far+B_near)*.5_ireals
            else
              t = exp(-tau/mu)
              B_eff_mu = (-B_near + B_far * t)/(-1._ireals + t) + ((B_far - B_near) *mu)/tau
            endif
          end function
      end subroutine

      subroutine schwarzschild(Nmu, dtau, albedo, Edn, Eup, planck, opt_srfc_emission)
        integer(iintegers), intent(in) :: Nmu
        real(ireals),intent(in),dimension(:) :: dtau
        real(ireals),intent(in) :: albedo
        real(ireals),dimension(:),intent(out):: Edn,Eup
        real(ireals),dimension(:),intent(in) :: planck
        ! planck value specific for surface. i.e. overrides last planck value
        real(ireals), intent(in), optional :: opt_srfc_emission

        integer(iintegers) :: imu,k,ke,ke1

        real(ireals) :: T(size(dtau)) ! Transmission coefficients
        real(ireals) :: Lup, Ldn, B, Bsrfc
        real(ireals) :: dmu, mu
        real(ireals) :: legendre_wi(Nmu)
        real(ireals) :: legendre_pt(Nmu)
        logical, parameter :: use_legendre=.True.

        Edn=zero
        Eup=zero

        ke = size(dtau)
        ke1 = ke+1

        Bsrfc = get_arg(planck(ke1), opt_srfc_emission)

        if(use_legendre) then
          call dgauss(size(legendre_wi), legendre_pt, legendre_wi)
          do imu=1,size(legendre_pt)
            mu = legendre_pt(imu)

            T = exp(- dtau/mu)

            ! zero incoming radiation at TOA
            Ldn = zero
            Edn(1) = Edn(1) + Ldn*mu*legendre_wi(imu)

            do k=1,ke
              call B_eff(planck(k), planck(k+1), dtau(k), B)
              Ldn = Ldn * T(k) + B*(one-T(k))
              Edn(k+1) = Edn(k+1) + Ldn*mu*legendre_wi(imu)
            enddo

            ! Boundary conditions at surface
            Lup = Bsrfc * (one-albedo) + albedo*Ldn
            Eup(ke1) = Eup(ke1) + Lup*mu*legendre_wi(imu)

            do k=ke,1,-1
              call B_eff(planck(k+1), planck(k), dtau(k), B)
              Lup = Lup * T(k) + B*(one-T(k))
              Eup(k) = Eup(k) + Lup*mu*legendre_wi(imu)
            enddo
          enddo ! enddo mu

          Eup = Eup*2*pi
          Edn = Edn*2*pi

        else

          dmu = one/real(Nmu, ireals)
          ! Transmission coefficients
          do imu=1,Nmu
            mu = (real(imu, ireals)-.5_ireals)*dmu
            T = exp(- dtau/mu)

            ! zero incoming radiation at TOA
            Ldn = zero
            Edn(1) = Edn(1) + Ldn*mu

            do k=1,ke
              call B_eff(planck(k), planck(k+1), dtau(k), B)
              Ldn = Ldn * T(k) + B*(one-T(k))
              Edn(k+1) = Edn(k+1) + Ldn*mu
            enddo

            ! Boundary conditions at surface
            Lup = Bsrfc * (one-albedo) + albedo*Ldn
            Eup(ke1) = Eup(ke1) + Lup*mu

            do k=ke,1,-1
              call B_eff(planck(k+1), planck(k), dtau(k), B)
              Lup = Lup * T(k) + B*(one-T(k))
              Eup(k) = Eup(k) + Lup*mu
            enddo

          enddo

          Eup = Eup*2*pi*dmu
          Edn = Edn*2*pi*dmu
        endif
      end subroutine

      subroutine  dgauss( m, gmu, gwt )
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
        real(ireals), dimension(:) :: gmu,gwt
        integer :: iter, lim, m, np1, k, nn

        double precision :: cona, en, nnp1, p, pm1, pm2, ppr, p2pri, prod, &
          tmp, x, xi, t

        integer, parameter :: maxit = 1000
        double precision, parameter :: one=1, two=2, eps=1.D-14

        if ( m.lt.1 ) then
          write(*,*)'dgauss--bad value of m'
          stop
        endif

        if ( m.eq.1 )  then
          gmu( 1 ) = 0.5
          gwt( 1 ) = 1.0
          return
        end if

        en   = m
        np1  = m + 1
        nnp1 = m * np1
        cona = dble( m-1 ) / ( 8 * m**3 )

        lim  = m / 2
        do  k = 1, lim

          !        initial guess for k-th root of legendre polynomial,
          !        from davis/rabinowitz (2.7.3.3a)

          t = ( 4*k - 1 ) * pi / ( 4*m + 2 )
          x = cos ( t + cona / tan( t ) )
          iter = 0

          !        upward recurrence for legendre polynomials

10        iter = iter + 1
          pm2 = one
          pm1 = x
          do nn = 2, m
            p   = ( ( 2*nn - 1 ) * x * pm1 - ( nn-1 ) * pm2 ) / nn
            pm2 = pm1
            pm1 = p
          enddo

          !        newton method

          tmp   = one / ( one - x**2 )
          ppr   = en * ( pm2 - x * p ) * tmp
          p2pri = ( two * x * ppr - nnp1 * p ) * tmp
          xi    = x - ( p / ppr ) * ( one +( p / ppr ) * p2pri / ( two* ppr ) )

          !        check for convergence

          if ( abs(xi-x) .gt. eps) then
            if( iter.gt.maxit ) then
              write(*,*)'dgauss--max iteration count'
              stop
            endif
            x = xi
            go to 10
          endif

          !        iteration finished--calculate weights, abscissae for (-1,1)

          gmu( k ) = - real(x, kind(gmu))
          gwt( k ) = real(two / ( tmp * ( en * pm2 )**2 ), kind(gmu))
          gmu( np1 - k ) = - gmu( k )
          gwt( np1 - k ) =   gwt( k )

        enddo

        !     set middle abscissa and weight for rules of odd order

        if ( mod( m,2 ) .ne. 0 )  then
          gmu( lim + 1 ) = 0.0
          prod = one
          do  k = 3, m, 2
            prod = prod * k / ( k-1 )
          enddo
          gwt( lim + 1 ) = real(two / prod**2, kind(gmu))
        endif

        !     convert from (-1,1) to (0,1) and resort

        do k = 1, m
          gmu(k) = 0.5 * gmu( k ) + 0.5
          gwt(k) = 0.5 * gwt( k )
        enddo

        return
      end subroutine dgauss
    end module
