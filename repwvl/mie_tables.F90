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
    & mpiint

  use m_helper_functions, only: &
    & CHKERR, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_bcast, &
    & toStr

  use m_netcdfIO, only: ncload

  use m_tenstream_interpolation, only: interp_2d
  use m_search, only: find_real_location

  implicit none

  private
  public :: mie_tables_init, mie_water_table, mie_ice_table, t_mie_table, mie_optprop

  logical, parameter :: ldebug = .true.

  type t_mie_table
    real(irealLUT), allocatable :: wvl(:)    ! in [mu]
    real(irealLUT), allocatable :: reff(:)   ! in [mu]
    real(irealLUT), allocatable :: qext(:, :) ! dim [reff, wvl], extinction coeff in km^-1 / (g/m^3)
    real(irealLUT), allocatable :: w0(:, :)   ! dim [reff, wvl], single scatter albedo
    real(irealLUT), allocatable :: g(:, :)    ! dim [reff, wvl], asymmetry parameter
  end type

  type(t_mie_table), allocatable :: mie_water_table
  type(t_mie_table), allocatable :: mie_ice_table
contains
  subroutine load_data(fname, mie_table, ierr, lverbose)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    type(t_mie_table), allocatable, intent(inout) :: mie_table
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2)

    ierr = 0

    if (allocated(mie_table)) return
    allocate (mie_table)

    groups(1) = trim(fname)
    groups(2) = 'wvl'; call ncload(groups, mie_table%wvl, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'reff'; call ncload(groups, mie_table%reff, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'ext'; call ncload(groups, mie_table%qext, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'ssa'; call ncload(groups, mie_table%w0, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'g'; call ncload(groups, mie_table%g, ierr, lverbose); call CHKERR(ierr)

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

    if (.not. allocated(mie_table)) allocate (mie_table)
    call imp_bcast(comm, mie_table%wvl, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%reff, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%qext, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%w0, ierr); call CHKERR(ierr)
    call imp_bcast(comm, mie_table%g, ierr); call CHKERR(ierr)
  end subroutine

  subroutine mie_tables_init(comm, ierr, lverbose, path_water_table, path_ice_table)
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose
    character(len=*), intent(in), optional :: path_water_table, path_ice_table

    integer(mpiint) :: myid

    character(len=default_str_len), parameter :: default_water_path = "mie.wc.table.nc"
    character(len=default_str_len), parameter :: default_ice_path = "mie.ic.table.nc"

    character(len=default_str_len) :: water_path
    character(len=default_str_len) :: ice_path

    logical :: lset, lexists

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    water_path = trim(get_arg(default_water_path, path_water_table))
    ice_path = trim(get_arg(default_ice_path, path_ice_table))

    call get_petsc_opt('', '-mie_wc', water_path, lset, ierr); call CHKERR(ierr)
    call get_petsc_opt('', '-mie_ic', ice_path, lset, ierr); call CHKERR(ierr)

    inquire (file=trim(water_path), exist=lexists)
    if (.not. lexists) then
      call CHKERR(1_mpiint, "File at mie_table path "//toStr(water_path)// &
        & " does not exist"//new_line('')// &
        & " please make sure the file is at this location"// &
        & " or specify a correct path with option"//new_line('')// &
        & "   -mie_wc <path>")
    end if

    inquire (file=trim(ice_path), exist=lexists)
    if (.not. lexists) then
      call CHKERR(1_mpiint, "File at mie_table path "//toStr(ice_path)// &
        & " does not exist"//new_line('')// &
        & " please make sure the file is at this location"// &
        & " or specify a correct path with option"//new_line('')// &
        & "   -mie_ic <path>")
    end if

    if (myid .eq. 0) then
      if (.not. allocated(mie_water_table)) then
        call load_data(water_path, mie_water_table, ierr, lverbose); call CHKERR(ierr)
      end if
      if (.not. allocated(mie_ice_table)) then
        call load_data(ice_path, mie_ice_table, ierr, lverbose); call CHKERR(ierr)
      end if
    end if
    call distribute_table(comm, mie_water_table, ierr); call CHKERR(ierr)
    call distribute_table(comm, mie_ice_table, ierr); call CHKERR(ierr)
  end subroutine

  subroutine mie_optprop(table, wvl, reff, qext, w0, g, ierr)
    type(t_mie_table), intent(in) :: table
    real(ireals), intent(in) :: wvl, reff
    real(ireals), intent(out) :: qext, w0, g
    integer(mpiint), intent(out) :: ierr

    real(irealLUT) :: pt(2)

    ierr = 0

    pt(1) = find_real_location(table%reff, reff)
    pt(2) = find_real_location(table%wvl, wvl)
    call interp_2d(pt, table%qext, qext)
    call interp_2d(pt, table%w0, w0)
    call interp_2d(pt, table%g, g)
  end subroutine
end module
