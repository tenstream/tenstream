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

module f2c_tenstream

      use iso_c_binding

      use m_data_parameters, only : iintegers, ireals, default_str_len, init_mpi_data_parameters

      use m_tenstr_rrtmg, only : tenstream_rrtmg, destroy_tenstream_rrtmg

#include "petsc/finclude/petsc.h"
      use petsc
      implicit none

      private
      public :: f2c_tenstream_rrtmg, f2c_destroy_tenstream_rrtmg

contains

  subroutine f2c_tenstream_rrtmg(comm, Nz, Nx, Ny, dx, dy, &
      phi0, theta0, albedo_thermal, albedo_solar,          &
      c_atm_filename, c_lthermal, c_lsolar,                &
      Nz_merged, cptr_edir, cptr_edn, cptr_eup, cptr_abso, &
      d_plev, d_tlev, d_lwc, d_reliq, d_iwc, d_reice,      &
      nprocx, nxproc, nprocy, nyproc) bind(C)

    integer(c_int), value :: comm
    integer(c_int), intent(in) :: Nx, Ny, Nz                       ! local subdomain size
    real(c_double), intent(in) :: dx, dy                           ! horizontal grid spacing in [m]
    real(c_double), intent(in) :: phi0, theta0                     ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(c_double), intent(in) :: albedo_thermal, albedo_solar     ! broadband ground albedo for solar and thermal spectrum

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(kind=c_char,len=1), intent(in) :: c_atm_filename(*)

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    integer(c_int), intent(in) :: c_lthermal, c_lsolar

    integer(c_int), intent(out) :: Nz_merged                       ! will determine the number of layers of the result
    type(c_ptr), intent(out) :: cptr_edir, cptr_edn, cptr_eup      ! fluxes edir, edn, eup have shape(Nz_merged+1, Nx, Ny)
    type(c_ptr), intent(out) :: cptr_abso                          ! abso just (Nz_merged, Nx, Ny)

    real(c_double), dimension(Nz+1, Nx, Ny), intent(in) :: d_plev  ! pressure on layer interfaces    [hPa]
    real(c_double), dimension(Nz+1, Nx, Ny), intent(in) :: d_tlev  ! Temperature on layer interfaces [K]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_lwc   ! liq water content               [g/kg]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_reliq ! effective radius                [micron]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_iwc   ! ice water content               [g/kg]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_reice ! ice effective radius            [micron]

    integer(c_int), intent(in) :: nprocx, nprocy                  ! number of processors in x and y
    integer(c_int), intent(in) :: nxproc(nprocx), nyproc(nprocy)  ! local size of subdomain along x and y


    ! From here on local variables
    real(ireals), allocatable, target, save, dimension(:,:,:) :: edir, edn, eup, abso

    character(default_str_len) :: atm_filename
    logical :: lthermal, lsolar

    call init_mpi_data_parameters(comm)

    if(kind(dx).ne.ireals) then
      print *,'ERROR: real datatypes in C_wrapper dont fit :: wrapper uses ',kind(dx), 'but Tenstream expects',ireals
      print *,'You have two possibilities: compile petsc with different datatype or copy files here'
      stop 'wrong datatypes'
    endif
    if(kind(Nx).ne.iintegers) then
      print *,'ERROR: integer datatypes in C_wrapper dont fit :: wrapper uses ',kind(dx), 'but Tenstream expects',ireals
      print *,'You have two possibilities: compile petsc with different datatype or copy files here'
      stop 'wrong datatypes'
    endif

    atm_filename = c_to_f_string(c_atm_filename)

    lthermal = c_int_2_logical(c_lthermal)
    lsolar   = c_int_2_logical(c_lsolar)

    call tenstream_rrtmg(comm, dx, dy, phi0, theta0,  &
      albedo_thermal, albedo_solar, atm_filename,     &
      lthermal, lsolar,                               &
      edir,edn,eup,abso,                              &
      d_plev, d_tlev, d_lwc=d_lwc, d_reliq=d_reliq,   &
      d_iwc=d_iwc, d_reice=d_reice,                   &
      nxproc=nxproc, nyproc=nyproc)

    cptr_edir = c_loc(edir)
    cptr_edn  = c_loc(edn )
    cptr_eup  = c_loc(eup )
    cptr_abso = c_loc(abso)
    Nz_merged = size(abso,dim=1)
  end subroutine

  subroutine f2c_destroy_tenstream_rrtmg(c_lfinalizepetsc) bind(C) ! Tidy up the solver
    integer(c_int), intent(in) :: c_lfinalizepetsc  ! determines if we drop the Petsc Environment. If you called petsc initialize in the C program, say False, i.e. 0
    call destroy_tenstream_rrtmg(c_int_2_logical(c_lfinalizepetsc))
  end subroutine

  function c_to_f_string(s) result(str)
    use iso_c_binding
    character(kind=c_char,len=1), intent(in) :: s(*)
    character(len=:), allocatable :: str
    integer i, nchars
    i = 1
    do
      if (s(i) == c_null_char) exit
      i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: str)
    str = transfer(s(1:nchars), str)
  end function c_to_f_string

  logical function c_int_2_logical(i)
    integer(c_int) :: i
    if(i.eq.0) then
      c_int_2_logical = .False.
    else
      c_int_2_logical = .True.
    endif
  end function

end module
