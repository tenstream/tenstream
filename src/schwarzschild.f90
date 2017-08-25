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

      use m_data_parameters, only: ireals,iintegers,zero,one,pi
      implicit none

      private
      public schwarzschild

    contains

      subroutine schwarzschild(dtau, albedo, Edn, Eup, planck)
        real(ireals),intent(in),dimension(:) :: dtau
        real(ireals),intent(in) :: albedo
        real(ireals),dimension(:),intent(out):: Edn,Eup
        real(ireals),dimension(:),intent(in) :: planck

        integer(iintegers) :: imu,k,ke,ke1

        integer(iintegers),parameter :: Nmu = 20
        real(ireals),parameter       :: dmu = one/Nmu


        real(ireals) :: T(size(dtau)) ! Transmission coefficients
        real(ireals) :: Lup, Ldn
        real(ireals) :: mu

        Edn=zero
        Eup=zero

        ke = size(dtau)
        ke1 = ke+1

        ! Transmission coefficients
        do imu=1,Nmu
          mu = (imu-.5_ireals)*dmu
          T = exp(- dtau/mu)

          ! Boundary conditions at surface
          Lup = planck(ke1) * (one-albedo)
          Eup(ke1) = Eup(ke1) + Lup*mu

          ! zero incoming radiation at TOA
          Ldn = zero
          Edn(1) = Edn(1) + Ldn*mu

          do k=ke,1,-1
            Lup = Lup * T(k) + planck(k)*(one-T(k))
            Eup(k) = Eup(k) + Lup*mu
          enddo
          do k=1,ke
            Ldn = Ldn * T(k) + planck(k)*(one-T(k))
            Edn(k+1) = Edn(k+1) + Ldn*mu
          enddo

        enddo

        Eup = Eup*2*pi*dmu
        Edn = Edn*2*pi*dmu

        end subroutine


end module
