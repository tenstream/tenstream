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

module m_f2c_twomax
  use iso_c_binding
  use iso_fortran_env, only: real64, int32
  use m_data_parameters, only: mpiint
  use m_helper_functions, only: CHKERR

  implicit none

  interface
    integer(c_int) function twostream_maxrandF90( &
      dtau_c, omega0_c, g_c, &
      dtau_f, omega0_f, g_f, &
      cf, nlev, S0, mu0, Ag, &
      Bg, B, delta, flagSolar, flagThermal, &
      Edir, Edn, Eup) bind(c, name='twostream_maxrand')
      use iso_c_binding
      integer(c_int), value :: nlev
      real(c_double) :: dtau_c(1:Nlev - 1)
      real(c_double) :: omega0_c(1:Nlev - 1)
      real(c_double) :: g_c(1:Nlev - 1)
      real(c_double) :: dtau_f(1:Nlev - 1)
      real(c_double) :: omega0_f(1:Nlev - 1)
      real(c_double) :: g_f(1:Nlev - 1)
      real(c_double) :: cf(1:Nlev - 1)
      real(c_double), value :: S0, mu0, Ag, Bg
      real(c_double) :: B(1:Nlev)
      integer(c_int), value :: delta, flagSolar, flagThermal
      real(c_double) :: Edir(1:Nlev)
      real(c_double) :: Edn(1:Nlev)
      real(c_double) :: Eup(1:Nlev)
    end function
  end interface
end module
