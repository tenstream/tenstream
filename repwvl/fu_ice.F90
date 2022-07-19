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
    & mpiint, &
    & share_dir

  use m_helper_functions, only: &
    & CHKERR, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_bcast, &
    & toStr

  use m_netcdfIO, only: ncload, get_global_attribute

  use m_tenstream_interpolation, only: interp_2d
  use m_search, only: find_real_location
  use m_repwvl_base, only: t_repwvl_data

  implicit none

  private
  public :: &
    & fu_ice_data_solar,   &
    & fu_ice_data_thermal, &
    & fu_ice_init,         &
    & fu_ice_optprop,      &
    & t_fu96_ice_data,     &
    & t_fu98_ice_data

  type t_fu96_ice_data
    logical :: is_repwvl                      ! flag to tell if its a general fu file or specifically for a repwvl case,
    ! this determines if we need to do binary searches for wavelengts
    ! or just use the indices
    real(irealLUT), allocatable :: wvl(:)     ! in [mu]
    real(irealLUT), allocatable :: ext(:, :)  ! dim [2, band], Extinction coefficient b = IWC * (a0 + a1/D), eq. 3.9a
    real(irealLUT), allocatable :: w0(:, :)   ! dim [4, band], Coalbedo  1-w  =  b0  +  b1 * D  +  b2 * D**2  +b3 * D**3, eq. 3.9b
    real(irealLUT), allocatable :: g(:, :)    ! dim [4, band], Asymmetry factor g  =  c0  +  c1 * D  +  c2 * D**2  +c3 * D**3, eq. 3.9c
  end type

  type t_fu98_ice_data
    logical :: is_repwvl                      ! flag to tell if its a general fu file or specifically for a repwvl case,
    ! this determines if we need to do binary searches for wavelengts
    ! or just use the indices
    real(irealLUT), allocatable :: wvl(:)     ! in [mu]
    real(irealLUT), allocatable :: ext(:, :)  ! dim [3, band], Extinction coefficient b = IWC * (a0 + a1/D + a2/D**2)
    real(irealLUT), allocatable :: abso(:, :) ! dim [4, band], Absorption coefficient b = IWC/D * (b0 + b1*D + b2*D**2 + b3*D**3)
    real(irealLUT), allocatable :: g(:, :)    ! dim [4, band], Asymmetry parameter g = c0 +  c1 * D  + c2 * D**2  +  c3 * D**3
  end type

  type(t_fu96_ice_data), allocatable :: fu_ice_data_solar
  type(t_fu98_ice_data), allocatable :: fu_ice_data_thermal

  interface fu_ice_optprop
    module procedure fu_ice_optprop_solar_index
    module procedure fu_ice_optprop_solar_general
    module procedure fu_ice_optprop_thermal_index
    module procedure fu_ice_optprop_thermal_general
  end interface

  real(ireals), parameter :: MaxAsymmetryFactor = 1.0_ireals - 10.0_ireals * epsilon(1.0_ireals)
  real(ireals), parameter :: MaxEffectiveRadius = 100.0e-6_ireals ! [metres]

contains

  subroutine load_data96(fname, data96, ierr, lverbose)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    type(t_fu96_ice_data), allocatable, intent(inout) :: data96
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2), attr
    integer(mpiint) :: have_attr

    ierr = 0

    if (allocated(data96)) return
    allocate (data96)

    groups(1) = trim(fname)
    groups(2) = 'fu96_wvl'; call ncload(groups, data96%wvl, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'fu96_ext'; call ncload(groups, data96%ext, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'fu96_ssa'; call ncload(groups, data96%w0, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'fu96_asy'; call ncload(groups, data96%g, ierr, lverbose); call CHKERR(ierr)

    attr = 'undef'
    call get_global_attribute(fname, 'repwvl_file_96', attr, have_attr)
    data96%is_repwvl = have_attr .eq. 0

    if (get_arg(.false., lverbose)) then
      print *, 'Loaded Fu96 Ice data '//new_line('')// &
        & ' wvl  ('//toStr(shape(data96%wvl))//')'//toStr(data96%wvl)//new_line('')// &
        & ' qext ('//toStr(shape(data96%ext))//')'//toStr(data96%ext)//new_line('')// &
        & ' w0   ('//toStr(shape(data96%w0))//')'//toStr(data96%w0)//new_line('')// &
        & ' g    ('//toStr(shape(data96%g))//')'//toStr(data96%g)//new_line('')// &
        & ' is_repwvl?'//toStr(data96%is_repwvl)//' <'//trim(attr)//'>'
    end if
  end subroutine

  subroutine load_data98(fname, data98, ierr, lverbose)
    character(len=*), intent(in) :: fname
    integer(mpiint), intent(out) :: ierr
    type(t_fu98_ice_data), allocatable, intent(inout) :: data98
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2), attr
    integer(mpiint) :: have_attr

    ierr = 0

    if (allocated(data98)) return
    allocate (data98)

    groups(1) = trim(fname)
    groups(2) = 'fu98_wvl'; call ncload(groups, data98%wvl, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'fu98_ext'; call ncload(groups, data98%ext, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'fu98_abs'; call ncload(groups, data98%abso, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'fu98_asy'; call ncload(groups, data98%g, ierr, lverbose); call CHKERR(ierr)

    attr = 'undef'
    call get_global_attribute(fname, 'repwvl_file_98', attr, have_attr)
    data98%is_repwvl = have_attr .eq. 0

    if (get_arg(.false., lverbose)) then
      print *, 'Loaded Fu98 Ice data with '//new_line('')// &
        & ' wvl  ('//toStr(shape(data98%wvl))//')'//toStr(data98%wvl)//new_line('')// &
        & ' qext ('//toStr(shape(data98%ext))//')'//toStr(data98%ext)//new_line('')// &
        & ' w0   ('//toStr(shape(data98%abso))//')'//toStr(data98%abso)//new_line('')// &
        & ' g    ('//toStr(shape(data98%g))//')'//toStr(data98%g)//new_line('')// &
        & ' is_repwvl?'//toStr(data98%is_repwvl)//' <'//trim(attr)//'>'
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
    call imp_bcast(comm, data96%is_repwvl, ierr); call CHKERR(ierr)
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
    call imp_bcast(comm, data98%is_repwvl, ierr); call CHKERR(ierr)
  end subroutine

  subroutine fu_ice_init(  &
      & comm,              &
      & ierr,              &
      & lverbose,          &
      & path)
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose
    character(len=*), intent(in), optional :: path

    integer(mpiint) :: myid

    character(len=default_str_len), parameter :: default_path = share_dir//"fu.ice.general.nc"

    character(len=default_str_len) :: prepwvl

    logical :: lset, lexists

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      prepwvl = trim(get_arg(default_path, path))

      call get_petsc_opt('', '-fu_ice', prepwvl, lset, ierr); call CHKERR(ierr)

      inquire (file=trim(prepwvl), exist=lexists)
      if (.not. lexists) then
        call CHKERR(1_mpiint, "File at fu ice path "//toStr(prepwvl)// &
          & " does not exist"//new_line('')// &
          & " please make sure the file is at this location"// &
          & " or specify a correct path with option"//new_line('')// &
          & "   -fu_ice <path>")
      end if

      if (.not. allocated(fu_ice_data_solar)) then
        call load_data96(prepwvl, fu_ice_data_solar, ierr, lverbose); call CHKERR(ierr)
      end if
      if (.not. allocated(fu_ice_data_thermal)) then
        call load_data98(prepwvl, fu_ice_data_thermal, ierr, lverbose); call CHKERR(ierr)
      end if
    end if

    call distribute_table96(comm, fu_ice_data_solar, ierr); call CHKERR(ierr)
    call distribute_table98(comm, fu_ice_data_thermal, ierr); call CHKERR(ierr)
  end subroutine

  subroutine fu_ice_optprop_solar_general(table, wvl, reff, qext, w0, g, ierr)
    type(t_fu96_ice_data), intent(in) :: table
    real(ireals), intent(in) :: wvl ! wvl in [mu]
    real(ireals), intent(in) :: reff
    real(ireals), intent(out) :: qext, w0, g
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: loc, wgt
    real(ireals) :: qext0, w00, g0
    real(ireals) :: qext1, w01, g1

    if (table%is_repwvl) &
      & call CHKERR(1_mpiint, 'tried to call fu_ice_optprop_solar with wvl value '// &
      & 'but loaded table is of type repwvl, use wvl index instead or load general table')
    ierr = 0
    loc = find_real_location(table%wvl, real(wvl, irealLUT))
    call fu_ice_optprop(table, int(floor(loc), iintegers), reff, qext0, w00, g0, ierr); call CHKERR(ierr)
    call fu_ice_optprop(table, int(ceiling(loc), iintegers), reff, qext1, w01, g1, ierr); call CHKERR(ierr)
    wgt = loc - floor(loc)
    qext = wgt * qext1 + (1._ireals - wgt) * qext0
    w0 = wgt * w01 + (1._ireals - wgt) * w00
    g = wgt * g1 + (1._ireals - wgt) * g0
  end subroutine

  subroutine fu_ice_optprop_solar_index(table, iwvl, reff, qext, w0, g, ierr)
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

  subroutine fu_ice_optprop_thermal_general(table, wvl, reff, qext, w0, g, ierr)
    type(t_fu98_ice_data), intent(in) :: table
    real(ireals), intent(in) :: wvl ! wvl in [mu]
    real(ireals), intent(in) :: reff
    real(ireals), intent(out) :: qext, w0, g
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: loc

    if (table%is_repwvl) &
      & call CHKERR(1_mpiint, 'tried to call fu_ice_optprop_thermal with wvl value '// &
      & 'but loaded table is of type repwvl, use wvl index instead or load general table')
    ierr = 0
    loc = find_real_location(table%wvl, real(wvl, irealLUT))
    call fu_ice_optprop(table, int(floor(loc), iintegers), reff, qext, w0, g, ierr); call CHKERR(ierr)
  end subroutine

  subroutine fu_ice_optprop_thermal_index(table, iwvl, reff, qext, w0, g, ierr)
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
