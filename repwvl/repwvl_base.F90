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

!> \page Routines to call tenstream with optical properties from a representative wavelength approach

module m_repwvl_base
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & default_str_len, &
    & iintegers, &
    & ireals, &
    & mpiint, &
    & share_dir

  use m_helper_functions, only: &
    & CHKERR, &
    & get_petsc_opt, &
    & get_arg, &
    & imp_bcast, &
    & toStr

  use m_netcdfIO, only: ncload

  implicit none

  private
  public :: &
    & repwvl_init, &
    & repwvl_log_events, &
    & setup_log_events, &
    & t_repwvl_data

  type t_repwvl_data
    ! dims xsec(pressure, wvl, tracer, temperature)
    !
    ! tracer as in ARTS:
    !   default(H2O_continuum, H2O, CO2, O3, N2O, CO, CH4, O2, HNO3, N2)
    !
    real(ireals), allocatable :: xsec(:, :, :, :) ! dim (Np_ref, Nwvl, Ntracer, Nt_pert)

    real(ireals), allocatable :: wvls(:)        ! dim (Nwvl)
    real(ireals), allocatable :: wgts(:)        ! dim (Nwvl)
    real(ireals), allocatable :: p_ref(:)       ! dim (Np_ref)
    real(ireals), allocatable :: t_ref(:)       ! dim (Nt_ref)
    real(ireals), allocatable :: t_pert(:)      ! dim (Nt_pert)
    real(ireals), allocatable :: vmrs_ref(:, :) ! dim (Np_ref, Ntracers)
    real(ireals), allocatable :: crs_o3(:, :)  ! additional tracer cross sections dim (3, Nwvl)
    real(ireals), allocatable :: crs_n2o(:, :) ! additional tracer cross sections dim (3, Nwvl)

    integer(iintegers) :: Nwvl, Ntracer
  end type

  type t_repwvl_log_events
    PetscLogStage :: stage_repwvl_solar
    PetscLogStage :: stage_repwvl_thermal
    PetscLogEvent :: repwvl_optprop
    PetscLogEvent :: repwvl_optprop_dtau
    PetscLogEvent :: repwvl_optprop_rayleigh
    PetscLogEvent :: repwvl_optprop_mie
    PetscLogEvent :: repwvl_optprop_fu_ice
  end type
  type(t_repwvl_log_events) :: repwvl_log_events

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .true.
#else
  logical, parameter :: ldebug = .true.
#endif

contains

  subroutine load_data(fname, repwvl_data, ierr, lverbose)
    character(len=*), intent(in) :: fname
    type(t_repwvl_data), allocatable, intent(inout) :: repwvl_data
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose

    integer(mpiint) :: ierr_ignore
    character(len=default_str_len) :: groups(2)

    if (allocated(repwvl_data)) return
    allocate (repwvl_data)

    groups(1) = trim(fname)
    groups(2) = 'xsec'; call ncload(groups, repwvl_data%xsec, ierr, lverbose); call CHKERR(ierr)
    repwvl_data%Nwvl = size(repwvl_data%xsec, dim=2)
    repwvl_data%Ntracer = size(repwvl_data%xsec, dim=3)

    groups(2) = 'wvl'; call ncload(groups, repwvl_data%wvls, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'wgts'; call ncload(groups, repwvl_data%wgts, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'p_ref'; call ncload(groups, repwvl_data%p_ref, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 't_ref'; call ncload(groups, repwvl_data%t_ref, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 't_pert'; call ncload(groups, repwvl_data%t_pert, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'vmrs_ref'; call ncload(groups, repwvl_data%vmrs_ref, ierr, lverbose); call CHKERR(ierr)

    groups(2) = 'crs_o3'; call ncload(groups, repwvl_data%crs_o3, ierr_ignore, lverbose)
    groups(2) = 'crs_n2o'; call ncload(groups, repwvl_data%crs_n2o, ierr_ignore, lverbose)

    if (get_arg(.false., lverbose)) then
      print *, 'xsec shape', shape(repwvl_data%xsec)
      print *, 'wvls', repwvl_data%wvls, 'shape', shape(repwvl_data%wvls)
      print *, 'wgts', repwvl_data%wgts, 'shape', shape(repwvl_data%wgts)
      print *, 'p_ref', repwvl_data%p_ref, 'shape', shape(repwvl_data%p_ref)
      print *, 't_ref', repwvl_data%t_ref, 'shape', shape(repwvl_data%t_ref)
      print *, 't_pert', repwvl_data%t_pert, 'shape', shape(repwvl_data%t_pert)
      print *, 'vmrs_ref', repwvl_data%vmrs_ref, 'shape', shape(repwvl_data%vmrs_ref)
      if (allocated(repwvl_data%crs_o3)) then
        print *, 'crs_o3', repwvl_data%crs_o3, 'shape', shape(repwvl_data%crs_o3)
      end if
      if (allocated(repwvl_data%crs_n2o)) then
        print *, 'crs_n2o', repwvl_data%crs_n2o, 'shape', shape(repwvl_data%crs_n2o)
      end if
    end if
  end subroutine

  subroutine repwvl_init(comm, repwvl_data_solar, repwvl_data_thermal, ierr, fname_repwvl_solar, fname_repwvl_thermal, lverbose)
    integer(mpiint), intent(in) :: comm
    type(t_repwvl_data), allocatable, intent(inout) :: repwvl_data_solar, repwvl_data_thermal
    integer(mpiint), intent(out) :: ierr
    character(len=*), intent(in), optional :: fname_repwvl_solar, fname_repwvl_thermal
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: basepath, fname_thermal, fname_solar
    character(len=*), parameter :: fname_thermal_default = share_dir//"repwvl_thermal.lut"
    character(len=*), parameter :: fname_solar_default = share_dir//"repwvl_solar.lut"
    logical :: lset, lexists
    integer(mpiint) :: myid

    ierr = 0

    call mpi_comm_rank(comm, myid, ierr)

    call setup_log_events(repwvl_log_events)

    if (myid .eq. 0) then
      fname_thermal = get_arg(fname_thermal_default, fname_repwvl_thermal)
      fname_solar = get_arg(fname_solar_default, fname_repwvl_solar)

      basepath = ''
      call get_petsc_opt('', '-repwvl_data', basepath, lset, ierr); call CHKERR(ierr)
      if (lset) then
        fname_thermal = trim(basepath)//'/'//trim(fname_thermal)
        fname_solar = trim(basepath)//'/'//trim(fname_solar)
      end if

      call get_petsc_opt('', '-repwvl_data_thermal', fname_thermal, lset, ierr); call CHKERR(ierr)
      call get_petsc_opt('', '-repwvl_data_solar', fname_solar, lset, ierr); call CHKERR(ierr)

      inquire (file=trim(fname_thermal), exist=lexists)
      if (.not. lexists) then
        call CHKERR(1_mpiint, "File at repwvl thermal path "//toStr(fname_thermal)// &
          & " does not exist"//new_line('')// &
          & " please make sure the file is at this location"// &
          & " you may set a directory with option "//new_line('')// &
          & "   -repwvl_data  <path with file repwvl_thermal.lut>"//new_line('')// &
          & " or specify a correct path with option"//new_line('')// &
          & "   -repwvl_data_thermal <path>")
      end if

      inquire (file=trim(fname_solar), exist=lexists)
      if (.not. lexists) then
        call CHKERR(1_mpiint, "File at repwvl solar path "//toStr(fname_solar)// &
          & " does not exist"//new_line('')// &
          & " please make sure the file is at this location"// &
          & " you may set a directory with option "//new_line('')// &
          & "   -repwvl_data  <path with file repwvl_solar.lut>"//new_line('')// &
          & " or specify a correct path with option"//new_line('')// &
          & "   -repwvl_data_solar <path>")
      end if

      if (get_arg(.false., lverbose)) print *, 'Reading representative wavelength data from '//trim(fname_thermal)
      call load_data(fname_thermal, repwvl_data_thermal, ierr, lverbose=lverbose); call CHKERR(ierr)

      if (get_arg(.false., lverbose)) print *, 'Reading representative wavelength data from '//trim(fname_solar)
      call load_data(fname_solar, repwvl_data_solar, ierr, lverbose=lverbose); call CHKERR(ierr)
    end if

    call distribute_repwvl_table(comm, repwvl_data_thermal, ierr); call CHKERR(ierr)
    call distribute_repwvl_table(comm, repwvl_data_solar, ierr); call CHKERR(ierr)

  contains

    subroutine distribute_repwvl_table(comm, table, ierr)
      integer(mpiint), intent(in) :: comm
      type(t_repwvl_data), allocatable, intent(inout) :: table
      integer(mpiint), intent(out) :: ierr
      logical :: lhave_crs_o3, lhave_crs_n2o
      integer(mpiint), parameter :: sendid = 0

      ierr = 0
      if (.not. allocated(table)) allocate (table)
      call imp_bcast(comm, table%xsec, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%wvls, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%wgts, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%p_ref, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%t_ref, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%t_pert, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%vmrs_ref, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%Nwvl, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%Ntracer, sendid, ierr); call CHKERR(ierr)

      lhave_crs_o3 = allocated(table%crs_o3)
      call imp_bcast(comm, lhave_crs_o3, sendid, ierr); call CHKERR(ierr)
      if (lhave_crs_o3) then
        call imp_bcast(comm, table%crs_o3, sendid, ierr); call CHKERR(ierr)
      end if
      lhave_crs_n2o = allocated(table%crs_n2o)
      call imp_bcast(comm, lhave_crs_n2o, sendid, ierr); call CHKERR(ierr)
      if (lhave_crs_n2o) then
        call imp_bcast(comm, table%crs_n2o, sendid, ierr); call CHKERR(ierr)
      end if
    end subroutine
  end subroutine

  subroutine setup_log_events(logs, solvername)
    type(t_repwvl_log_events), intent(inout) :: logs
    character(len=*), optional :: solvername
    character(len=default_str_len) :: s
    PetscClassId :: cid
    integer(mpiint) :: ierr

    s = get_arg('tenstr_repwvl.', solvername)
    call PetscClassIdRegister(trim(s), cid, ierr); call CHKERR(ierr)

    call setup_stage(trim(s)//'repwvl_solar', logs%stage_repwvl_solar)
    call setup_stage(trim(s)//'repwvl_thermal', logs%stage_repwvl_thermal)

    call PetscLogEventRegister(trim(s)//'repwvl_optprop', cid, logs%repwvl_optprop, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'repwvl_optprop_dtau', cid, logs%repwvl_optprop_dtau, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'repwvl_optprop_rayleigh', cid, logs%repwvl_optprop_rayleigh, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'repwvl_optprop_mie', cid, logs%repwvl_optprop_mie, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'repwvl_optprop_fu_ice', cid, logs%repwvl_optprop_fu_ice, ierr); call CHKERR(ierr)

  contains
    subroutine setup_stage(stagename, logstage)
      character(len=*), intent(in) :: stagename
      PetscLogStage, intent(inout) :: logstage
      call PetscLogStageGetId(stagename, logstage, ierr); call CHKERR(ierr)
      if (logstage .lt. 0_iintegers) then
        call PetscLogStageRegister(stagename, logstage, ierr); call CHKERR(ierr)
      end if
    end subroutine
  end subroutine
end module
