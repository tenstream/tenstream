!----------------------------------------------------------------------------
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
! Copyright (C) 2010-2015  Carolin Klinger, <carolin.klinger@physik.lmu.de>
!----------------------------------------------------------------------------
! Neighbouring Column Approximation
! RTE Solver for the thermal spectral range, calculation of heating rates
! Klinger and Mayer, 2015, The Neighbouring Column Approximation - A fast approach for the calculation of 3D thermal heating rates, JQSRT
! carolin.klinger@physik.lmu.de
!----------------------------------------------------------------------------

module m_ts_nca
  use m_data_parameters, only : ireals, iintegers, zero
  implicit none

  private
  public :: ts_nca

contains
  subroutine ts_nca (dx, dy, dz, B, kabs_3d, Edn, Eup, hr)

    real(ireals), intent(in) :: dx, dy
    real(ireals), intent(in), dimension(:,:,:) :: dz, kabs_3d   ! dimensions including ghost values (Nlevel,Nx,Ny)
    real(ireals), intent(in), dimension(:,:,:) :: B,Edn,Eup     ! dimensions including ghost values (Nlevel,Nx,Ny)
    real(ireals), intent(out),dimension(:,:,:) :: hr            ! dimensions including ghost values (Nlevel,Nx,Ny)

    hr = dz + kabs_3d + B + Edn + Eup * (dx+dy) ! remove compiler warnings

    hr = -9999999
    stop 'NCA not freely available.... please consult Carolin Klinger for an implementation.'

  end subroutine ts_nca
end module m_ts_nca
