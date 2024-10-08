!-------------------------------------------------------------------------
! This file is part of the TenStream solver.
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

module m_mie_tables
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  use m_data_parameters, only: &
    & default_str_len, &
    & iintegers, &
    & irealLUT, &
    & ireals, &
    & mpiint, &
    & share_dir

  use m_helper_functions, only: &
    & CHKERR, &
    & deallocate_allocatable, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_bcast, &
    & toStr

  use m_netcdfIO, only: ncload, get_global_attribute

  use m_tenstream_interpolation, only: interp_2d
  use m_search, only: find_real_location

  implicit none

  private
  public :: &
    & destroy_mie_table, &
    & mie_optprop, &
    & mie_tables_init, &
    & t_mie_table

  logical, parameter :: ldebug = .true.

  type t_mie_table
    real(irealLUT), allocatable :: wvl(:)    ! in [mu]
    real(irealLUT), allocatable :: reff(:)   ! in [mu]
    real(irealLUT), allocatable :: qext(:, :) ! dim [reff, wvl], extinction coeff in km^-1 / (g/m^3)
    real(irealLUT), allocatable :: w0(:, :)   ! dim [reff, wvl], single scatter albedo
    real(irealLUT), allocatable :: g(:, :)    ! dim [reff, wvl], asymmetry parameter
  end type

  interface mie_optprop
    module procedure mie_optprop_general
    module procedure mie_optprop_index
  end interface

contains
  subroutine load_data(fname, mie_table, ierr, lverbose)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    type(t_mie_table), allocatable, intent(inout) :: mie_table
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2), creator_name
    logical :: lecckd_mie_file

    ierr = 0

    if (allocated(mie_table)) return
    allocate (mie_table)

    call get_global_attribute(fname, 'creator_name', creator_name, ierr)
    lecckd_mie_file = ierr .eq. 0 .and. trim(creator_name) .eq. "Robin J. Hogan"
    if (lecckd_mie_file) then
      groups(1) = trim(fname)
      groups(2) = 'wavelength'; call ncload(groups, mie_table%wvl, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'effective_radius'; call ncload(groups, mie_table%reff, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'mass_extinction_coefficient'; call ncload(groups, mie_table%qext, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'single_scattering_albedo'; call ncload(groups, mie_table%w0, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'asymmetry_factor'; call ncload(groups, mie_table%g, ierr, lverbose); call CHKERR(ierr)
      mie_table%wvl = mie_table%wvl * 1e6_ireallut
      mie_table%reff = mie_table%reff * 1e6_ireallut
      mie_table%qext = transpose(mie_table%qext)
      mie_table%w0 = transpose(mie_table%w0)
      mie_table%g = transpose(mie_table%g)
    else
      groups(1) = trim(fname)
      groups(2) = 'wvl'; call ncload(groups, mie_table%wvl, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'reff'; call ncload(groups, mie_table%reff, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'ext'; call ncload(groups, mie_table%qext, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'ssa'; call ncload(groups, mie_table%w0, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'g'; call ncload(groups, mie_table%g, ierr, lverbose); call CHKERR(ierr)
    end if

    if (get_arg(.false., lverbose)) then
      print *, 'Loaded Mie table:'//new_line('')// &
        & ' wvl  ('//toStr(shape(mie_table%wvl))//')'//toStr(mie_table%wvl)//new_line('')// &
        & ' reff ('//toStr(shape(mie_table%reff))//')'//toStr(mie_table%reff)//new_line('')// &
        & ' qext ('//toStr(shape(mie_table%qext))//')'//toStr(mie_table%qext)//new_line('')// &
        & ' w0   ('//toStr(shape(mie_table%w0))//')'//toStr(mie_table%w0)//new_line('')// &
        & ' g    ('//toStr(shape(mie_table%g))//')'//toStr(mie_table%g)
    end if
  end subroutine

  subroutine distribute_table(comm, mie_table, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_mie_table), allocatable, intent(inout) :: mie_table
    integer(mpiint), intent(out) :: ierr

    ierr = 0
    if (.not. allocated(mie_table)) allocate (mie_table)
    call imp_bcast(comm, mie_table%wvl, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%reff, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%qext, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%w0, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%g, 0_mpiint, ierr); call CHKERR(ierr)
  end subroutine

  subroutine mie_tables_init(comm, table, ierr, lverbose, path_table, prefix)
    integer(mpiint), intent(in) :: comm
    type(t_mie_table), allocatable, intent(inout) :: table
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose
    character(len=*), intent(in), optional :: path_table, prefix

    integer(mpiint) :: myid

    character(len=default_str_len), parameter :: default_path = share_dir//"mie_droplet_scattering.nc"

    character(len=default_str_len) :: path

    logical :: lset, lexists

    if (allocated(table)) return

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      path = trim(get_arg(default_path, path_table))

      call get_petsc_opt(get_arg('', prefix), '-mie_wc', path, lset, ierr); call CHKERR(ierr)

      inquire (file=trim(path), exist=lexists)
      if (.not. lexists) then
        call CHKERR(1_mpiint, "File at mie_table path "//toStr(path)// &
          & " does not exist"//new_line('')// &
          & " please make sure the file is at this location"// &
          & " or specify a correct path with option"//new_line('')// &
          & "   -mie_wc <path>")
      end if

      if (.not. allocated(table)) then
        call load_data(path, table, ierr, lverbose); call CHKERR(ierr)
      end if
    end if
    call distribute_table(comm, table, ierr); call CHKERR(ierr)
  end subroutine

  subroutine mie_optprop_general(table, wvl, reff, qext, w0, g, ierr)
    type(t_mie_table), intent(in) :: table
    real(ireals), intent(in) :: wvl, reff
    real(ireals), intent(out) :: qext, w0, g
    integer(mpiint), intent(out) :: ierr

    real(irealLUT) :: pt(2)
    real(irealLUT) :: rqext, rw0, rg

    ierr = 0

    pt(1) = find_real_location(table%reff, real(reff, irealLUT))
    pt(2) = find_real_location(table%wvl, real(wvl, irealLUT))
    call interp_2d(pt, table%qext, rqext)
    call interp_2d(pt, table%w0, rw0)
    call interp_2d(pt, table%g, rg)
    qext = real(rqext, irealLUT)
    w0 = real(rw0, irealLUT)
    g = real(rg, irealLUT)
  end subroutine

  subroutine mie_optprop_index(table, iwvl, reff, qext, w0, g, ierr)
    type(t_mie_table), intent(in) :: table
    integer(iintegers), intent(in) :: iwvl
    real(ireals), intent(in) :: reff
    real(ireals), intent(out) :: qext, w0, g
    integer(mpiint), intent(out) :: ierr

    real(irealLUT) :: pt(2)
    real(irealLUT) :: rqext, rw0, rg
    call CHKERR(1_mpiint, "optimized lookups for repwvl are not yet implemented."//&
      & " Note: one option would be to split solar and thermal and handle each case on its own."// &
      & " or we have them in one file but somehow make it easy for repwvl indices to query the mie tables."//&
      & " Have to think about it... For the moment, just use the general one")

    ierr = 0

    pt(1) = find_real_location(table%reff, real(reff, irealLUT))
    pt(2) = real(iwvl, kind(pt))
    call interp_2d(pt, table%qext, rqext)
    call interp_2d(pt, table%w0, rw0)
    call interp_2d(pt, table%g, rg)
    qext = real(rqext, irealLUT)
    w0 = real(rw0, irealLUT)
    g = real(rg, irealLUT)
  end subroutine

  subroutine destroy_mie_table(table, ierr)
    type(t_mie_table), allocatable, intent(inout) :: table
    integer(mpiint), intent(out) :: ierr

    ierr = 0
    if (.not. allocated(table)) return
    call deallocate_allocatable(table%wvl)
    call deallocate_allocatable(table%reff)
    call deallocate_allocatable(table%qext)
    call deallocate_allocatable(table%w0)
    call deallocate_allocatable(table%g)
    deallocate (table)
  end subroutine
end module
