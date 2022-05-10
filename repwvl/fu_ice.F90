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

module m_fu_ice
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
  public :: fu_ice_init, fu_ice_data_solar, fu_ice_data_thermal, t_fu96_ice_data, t_fu98_ice_data, fu_ice_optprop

  logical, parameter :: ldebug = .true.

  type t_fu96_ice_data
    real(irealLUT), allocatable :: wvl(:)     ! in [mu]
    real(irealLUT), allocatable :: ext(:, :)  ! dim [2, band], Extinction coefficient b = IWC * (a0 + a1/D), eq. 3.9a
    real(irealLUT), allocatable :: w0(:, :)   ! dim [4, band], Coalbedo  1-w  =  b0  +  b1 * D  +  b2 * D**2  +b3 * D**3, eq. 3.9b
    real(irealLUT), allocatable :: g(:, :)    ! dim [4, band], Asymmetry factor g  =  c0  +  c1 * D  +  c2 * D**2  +c3 * D**3, eq. 3.9c
  end type

  type t_fu98_ice_data
    real(irealLUT), allocatable :: wvl(:)     ! in [mu]
    real(irealLUT), allocatable :: ext(:, :)  ! dim [3, band], Extinction coefficient b = IWC * (a0 + a1/D + a2/D**2)
    real(irealLUT), allocatable :: abso(:, :) ! dim [4, band], Absorption coefficient b = IWC/D * (b0 + b1*D + b2*D**2 + b3*D**3)
    real(irealLUT), allocatable :: g(:, :)    ! dim [4, band], Asymmetry parameter g = c0 +  c1 * D  + c2 * D**2  +  c3 * D**3
  end type

  type(t_fu96_ice_data), allocatable :: fu_ice_data_solar
  type(t_fu98_ice_data), allocatable :: fu_ice_data_thermal

  interface fu_ice_optprop
    module procedure fu_ice_optprop_solar
    module procedure fu_ice_optprop_thermal
  end interface

  real(ireals), parameter :: MaxAsymmetryFactor = 1.0_ireals - 10.0_ireals * epsilon(1.0_ireals)
  real(ireals), parameter :: MaxEffectiveRadius = 100.0e-6_ireals ! [metres]

contains

  subroutine load_data96(fname, data96, ierr, lverbose)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    type(t_fu96_ice_data), allocatable, intent(inout) :: data96
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2)

    ierr = 0

    if (allocated(data96)) return
    allocate (data96)

    groups(1) = trim(fname)
    groups(2) = 'band_wvl'; call ncload(groups, data96%wvl, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'repwvl_fu96_ext'; call ncload(groups, data96%ext, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'repwvl_fu96_ssa'; call ncload(groups, data96%w0, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'repwvl_fu96_asy'; call ncload(groups, data96%g, ierr, lverbose); call CHKERR(ierr)

    if (get_arg(.false., lverbose)) then
      print *, 'Loaded Fu96 Ice data:'//new_line('')// &
        & ' wvl  ('//toStr(shape(data96%wvl))//')'//toStr(data96%wvl)//new_line('')// &
        & ' qext ('//toStr(shape(data96%ext))//')'//toStr(data96%ext)//new_line('')// &
        & ' w0   ('//toStr(shape(data96%w0))//')'//toStr(data96%w0)//new_line('')// &
        & ' g    ('//toStr(shape(data96%g))//')'//toStr(data96%g)
    end if
  end subroutine

  subroutine load_data98(fname, data98, ierr, lverbose)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    type(t_fu98_ice_data), allocatable, intent(inout) :: data98
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2)

    ierr = 0

    if (allocated(data98)) return
    allocate (data98)

    groups(1) = trim(fname)
    groups(2) = 'fu98_wvl'; call ncload(groups, data98%wvl, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'repwvl_fu98_ext'; call ncload(groups, data98%ext, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'repwvl_fu98_abs'; call ncload(groups, data98%abso, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'repwvl_fu98_asy'; call ncload(groups, data98%g, ierr, lverbose); call CHKERR(ierr)

    if (get_arg(.false., lverbose)) then
      print *, 'Loaded Fu98 Ice data:'//new_line('')// &
        & ' wvl  ('//toStr(shape(data98%wvl))//')'//toStr(data98%wvl)//new_line('')// &
        & ' qext ('//toStr(shape(data98%ext))//')'//toStr(data98%ext)//new_line('')// &
        & ' w0   ('//toStr(shape(data98%abso))//')'//toStr(data98%abso)//new_line('')// &
        & ' g    ('//toStr(shape(data98%g))//')'//toStr(data98%g)
    end if
  end subroutine

  subroutine distribute_table96(comm, data96, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_fu96_ice_data), allocatable, intent(inout) :: data96
    integer(mpiint), intent(out) :: ierr

    ierr = 0
    if (.not. allocated(data96)) allocate (data96)
    call imp_bcast(comm, data96%wvl, ierr); call CHKERR(ierr)
    call imp_bcast(comm, data96%ext, ierr); call CHKERR(ierr)
    call imp_bcast(comm, data96%w0, ierr); call CHKERR(ierr)
    call imp_bcast(comm, data96%g, ierr); call CHKERR(ierr)
  end subroutine

  subroutine distribute_table98(comm, data98, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_fu98_ice_data), allocatable, intent(inout) :: data98
    integer(mpiint), intent(out) :: ierr

    ierr = 0
    if (.not. allocated(data98)) allocate (data98)
    call imp_bcast(comm, data98%wvl, ierr); call CHKERR(ierr)
    call imp_bcast(comm, data98%ext, ierr); call CHKERR(ierr)
    call imp_bcast(comm, data98%abso, ierr); call CHKERR(ierr)
    call imp_bcast(comm, data98%g, ierr); call CHKERR(ierr)
  end subroutine

  subroutine fu_ice_init(comm, ierr, lverbose, path_solar, path_thermal)
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose
    character(len=*), intent(in), optional :: path_solar, path_thermal

    integer(mpiint) :: myid

    character(len=default_str_len), parameter :: default_path_solar = "fu.ice.repwvl_solar.nc"
    character(len=default_str_len), parameter :: default_path_thermal = "fu.ice.repwvl_thermal.nc"

    character(len=default_str_len) :: psol, pth

    logical :: lset, lexists

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    psol = trim(get_arg(default_path_solar, path_solar))
    pth = trim(get_arg(default_path_thermal, path_thermal))

    call get_petsc_opt('', '-fu_ice_solar', psol, lset, ierr); call CHKERR(ierr)
    call get_petsc_opt('', '-fu_ice_thermal', pth, lset, ierr); call CHKERR(ierr)

    inquire (file=trim(psol), exist=lexists)
    if (.not. lexists) then
      call CHKERR(1_mpiint, "File at fu ice solar path "//toStr(psol)// &
        & " does not exist"//new_line('')// &
        & " please make sure the file is at this location"// &
        & " or specify a correct path with option"//new_line('')// &
        & "   -fu_ice_solar <path>")
    end if

    inquire (file=trim(pth), exist=lexists)
    if (.not. lexists) then
      call CHKERR(1_mpiint, "File at fu ice thermal path "//toStr(pth)// &
        & " does not exist"//new_line('')// &
        & " please make sure the file is at this location"// &
        & " or specify a correct path with option"//new_line('')// &
        & "   -fu_ice_thermal <path>")
    end if

    if (myid .eq. 0) then
      if (.not. allocated(fu_ice_data_solar)) then
        call load_data96(psol, fu_ice_data_solar, ierr, lverbose); call CHKERR(ierr)
      end if
      if (.not. allocated(fu_ice_data_thermal)) then
        call load_data98(pth, fu_ice_data_thermal, ierr, lverbose); call CHKERR(ierr)
      end if
    end if
    call distribute_table96(comm, fu_ice_data_solar, ierr); call CHKERR(ierr)
    call distribute_table98(comm, fu_ice_data_thermal, ierr); call CHKERR(ierr)
  end subroutine

  subroutine fu_ice_optprop_solar(table, iwvl, reff, qext, w0, g, ierr)
    type(t_fu96_ice_data), intent(in) :: table
    integer(iintegers) :: iwvl ! index of repwvl bands
    real(ireals), intent(in) :: reff
    real(ireals), intent(out) :: qext, w0, g
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: de_um

    ! Convert to effective diameter using the relationship in the IFS (as in EcRad)
    de_um = min(reff, MaxEffectiveRadius) * (1.0e6_ireals / 0.64952_ireals)
    !iwp_gm_2  = ice_wp * 1000.0_jprb

    qext = table%ext(1, iwvl) + table%ext(2, iwvl) / de_um

    w0 = 1.0_ireals - (&
      & table%w0(1, iwvl) &
      & + de_um * (table%w0(2, iwvl) &
      &   + de_um * (table%w0(3, iwvl) &
      &     + de_um * table%w0(4, iwvl))))

    g = min( &
      & table%g(1, iwvl) &
      & + de_um * (table%g(2, iwvl) &
      &   + de_um * (table%g(3, iwvl) &
      &     + de_um * table%g(4, iwvl))) &
      & , MaxAsymmetryFactor)

    ierr = 0
  end subroutine

  subroutine fu_ice_optprop_thermal(table, iwvl, reff, qext, w0, g, ierr)
    type(t_fu98_ice_data), intent(in) :: table
    integer(iintegers) :: iwvl ! index of repwvl bands
    real(ireals), intent(in) :: reff
    real(ireals), intent(out) :: qext, w0, g
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: de_um, inv_de_um

    ! Convert to effective diameter using the relationship in the IFS (as in EcRad)
    de_um = min(reff, MaxEffectiveRadius) * (1.0e6_ireals / 0.64952_ireals)
    inv_de_um = 1._ireals / de_um

    qext = table%ext(1, iwvl) &
      & + inv_de_um * (table%ext(2, iwvl) &
      &   + inv_de_um * table%ext(3, iwvl))

    w0 = 1.0_ireals - (&
      & inv_de_um * (table%abso(1, iwvl) &
      & + de_um * (table%abso(2, iwvl) &
      &   + de_um * (table%abso(3, iwvl) &
      &     + de_um * table%abso(4, iwvl)))))

    g = min( &
      & table%g(1, iwvl) &
      & + de_um * (table%g(2, iwvl) &
      &   + de_um * (table%g(3, iwvl) &
      &     + de_um * table%g(4, iwvl))) &
      & , MaxAsymmetryFactor)

    ierr = 0
  end subroutine
end module
