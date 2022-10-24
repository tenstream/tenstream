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

!> \page Routines to call tenstream with optical properties from ECRAD

module m_ecckd_base
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & default_str_len, &
    & iintegers, &
    & ireals, &
    & irealLUT, &
    & mpiint, &
    & share_dir

  use m_helper_functions, only: &
    & CHKERR, &
    & CHKWARN, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_bcast, &
    & split, &
    & toStr

  use m_dyn_atm_to_rrtmg, only: &
    & t_tenstr_atm

  use m_netcdfIO, only: ncload, get_global_attribute

  use m_mie_tables, only: t_mie_table, mie_optprop
  use m_fu_ice, only: t_fu_muskatel_ice_data, fu_ice_optprop

  implicit none

  private
  public :: &
    & ecckd_init, &
    & ecckd_log_events, &
    & IConcDependenceLinear, &
    & IConcDependenceLUT, &
    & IConcDependenceNone, &
    & IConcDependenceRelativeLinear, &
    & setup_log_events, &
    & t_ecckd_atm_gas, &
    & t_ecckd_data, &
    & t_ecckd_mie_table

! Concentration dependence of individual gases
  enum, bind(c)
    enumerator :: IConcDependenceNone = 0, &
      &        IConcDependenceLinear, &
      &        IConcDependenceLUT, &
      &        IConcDependenceRelativeLinear
  end enum

  type t_ecckd_atm_gas
    character(len=default_str_len) :: id
    integer(iintegers) :: conc_dependence_code
    real(ireals) :: reference_mole_fraction
    real(ireals), pointer :: mole_fraction1(:)
    real(ireals), allocatable :: log_mole_fraction1(:)
    !real(ireals), pointer :: mole_fraction2(:, :)
    real(ireals), pointer :: molar_absorption_coeff3(:, :, :)
    real(ireals), pointer :: molar_absorption_coeff4(:, :, :, :)
    real(ireals), pointer :: vmr(:, :) ! dim(k,icol) pointer into a gas, usually from a t_tenstr_atm
    real(ireals) :: vmr_const ! const vmr of a gas
  end type

  type t_ecckd_mie_table
    real(irealLUT), allocatable :: reff(:)   ! in [mu]
    real(irealLUT), allocatable :: qext(:, :) ! dim(reff, gpt)
    real(irealLUT), allocatable :: w0(:, :) ! dim(reff, gpt)
    real(irealLUT), allocatable :: g(:, :) ! dim(reff, gpt)
  end type

  type t_ecckd_fu_ice_table
    ! for units look at t_fu_muskatel_ice_data
    real(irealLUT), allocatable :: reff(:)   ! in [mu]
    real(irealLUT), allocatable :: qext(:, :) ! dim(reff, gpt)
    real(irealLUT), allocatable :: w0(:, :) ! dim(reff, gpt)
    real(irealLUT), allocatable :: g(:, :) ! dim(reff, gpt)
  end type

  type t_ecckd_data
    ! assume that gases are:
    ! :constituent_id = "composite h2o o3 co2 ch4 n2o cfc11 cfc12" ;
    ! :composite_constituent_id = "o2 n2 n2o ch4" ;
    integer(iintegers) :: n_gases ! Number of gases treated
    character(len=default_str_len) :: constituent_id, composite_constituent_id
    integer(iintegers) :: n_g_pnt ! Number of g-points
    real(ireals), allocatable :: temperature(:, :)     ! dim(pressure, temperature) [K]
    real(ireals), allocatable :: pressure(:)           ! dim(pressure) [Pa]
    real(ireals), allocatable :: log_pressure(:)       ! dim(pressure) [Pa]
    real(ireals), allocatable :: temperature_planck(:) ! dim(temperature_planck) [K] Temperature for Planck function look-up table
    real(ireals), allocatable :: planck_function(:, :) ! dim(g_point, temperature_planck) [Wm-2] Planck function look-up table
    real(ireals), allocatable :: solar_irradiance(:) ! dim(g_point) [Wm-2] Solar irradiance across each g point
    real(ireals), allocatable :: rayleigh_molar_scattering_coeff(:) ! dim(g_point) [m2 mol-1] ayleigh molar scattering coefficient in each gpt
    real(ireals), allocatable :: wavenumber1(:)        ! dim(wavenumber) [cm-1] Lower wavenumber bound of spectral interval
    real(ireals), allocatable :: wavenumber2(:)        ! dim(wavenumber) [cm-1] Upper wavenumber bound of spectral interval
    real(ireals), allocatable :: gpoint_fraction(:, :) ! dim(wavenumber, g_point) Fraction of spectrum contributing to each g-point
    real(ireals), allocatable :: wavenumber1_band(:)   ! dim(band) [cm-1] Lower wavenumber bound of band
    real(ireals), allocatable :: wavenumber2_band(:)   ! dim(band) [cm-1] Upper wavenumber bound of band
    integer(iintegers), allocatable :: band_number(:)  ! dim(g_point) Band number of each g point

    integer(iintegers) :: composite_conc_dependence_code ! COMPOSITE concentration dependence code
    real(ireals), allocatable :: composite_molar_absorption_coeff(:, :, :) ! dim(g_point, pressure, temperature) [m2 mol-1] Molar absorption coefficient of background gases
    real(ireals), allocatable :: composite_mole_fraction(:, :) ! dim(pressure, composite_gas) [1] Mole fractions of the gases that make up COMPOSITE

    integer(iintegers) :: h2o_conc_dependence_code ! H2O concentration dependence code
    real(ireals), allocatable :: h2o_mole_fraction(:) ! dim(h2o_mole_fraction) [1] H2O mole fraction for look-up table
    real(ireals), allocatable :: h2o_molar_absorption_coeff(:, :, :, :) ! dim(g_point, pressure, temperature, h2o_mole_fraction) [m2 mol-1] Molar absorption coefficient of H2O

    integer(iintegers) :: o3_conc_dependence_code ! o3 concentration dependence code
    real(ireals), allocatable :: o3_molar_absorption_coeff(:, :, :) ! dim(g_point, pressure, temperature) [m2 mol-1] Molar absorption coefficient of o3

    integer(iintegers) :: co2_conc_dependence_code ! co2 concentration dependence code
    real(ireals), allocatable :: co2_molar_absorption_coeff(:, :, :) ! dim(g_point, pressure, temperature) [m2 mol-1] Molar absorption coefficient of co2

    integer(iintegers) :: ch4_conc_dependence_code ! ch4 concentration dependence code
    real(ireals) :: ch4_reference_mole_fraction ! Reference mole fraction of CH4
    real(ireals), allocatable :: ch4_molar_absorption_coeff(:, :, :) ! dim(g_point, pressurem temperature) [m2 mol-1] Molar absorption coefficient of ch4

    integer(iintegers) :: n2o_conc_dependence_code ! n2o concentration dependence code
    real(ireals) :: n2o_reference_mole_fraction ! Reference mole fraction of n2o
    real(ireals), allocatable :: n2o_molar_absorption_coeff(:, :, :) ! dim(g_point, pressure, temperature) [m2 mol-1] Molar absorption coefficient of n2o

    integer(iintegers) :: cfc11_conc_dependence_code ! cfc11 concentration dependence code
    real(ireals), allocatable :: cfc11_molar_absorption_coeff(:, :, :) ! dim(g_point, pressure, temperature) [m2 mol-1] Molar absorption coefficient of cfc11

    integer(iintegers) :: cfc12_conc_dependence_code ! cfc12 concentration dependence code
    real(ireals), allocatable :: cfc12_molar_absorption_coeff(:, :, :) ! dim(g_point, pressure, temperature) [m2 mol-1] Molar absorption coefficient of cfc12

    ! following is populated to have a mapping between atm%gas_vmrs and ecckd_data entries
    type(t_ecckd_atm_gas), allocatable :: gases(:) ! dim(n_gases)

    type(t_ecckd_mie_table), allocatable :: mie_table
    type(t_ecckd_fu_ice_table), allocatable :: fu_ice_table
  end type

  type t_ecckd_log_events
    PetscLogStage :: stage_ecckd_solar
    PetscLogStage :: stage_ecckd_thermal
    PetscLogEvent :: ecckd_optprop
    PetscLogEvent :: ecckd_optprop_dtau
    PetscLogEvent :: ecckd_optprop_rayleigh
    PetscLogEvent :: ecckd_optprop_mie
    PetscLogEvent :: ecckd_optprop_fu_ice
  end type
  type(t_ecckd_log_events) :: ecckd_log_events

  ! constant trace gas volume mixing ratios
  real(ireals), parameter :: CFC11 = 0e-9_ireals
  real(ireals), parameter :: CFC12 = 0e-9_ireals

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .false.
#else
  logical, parameter :: ldebug = .true.
#endif

contains

  subroutine load_data(fname, lsw, llw, ecckd_data, ierr, lverbose)
    character(len=*), intent(in) :: fname
    logical, intent(in) :: lsw, llw
    type(t_ecckd_data), allocatable, intent(inout) :: ecckd_data
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2)

    if (allocated(ecckd_data)) return
    allocate (ecckd_data)

    call get_global_attribute(fname, 'constituent_id', ecckd_data%constituent_id, ierr); call CHKERR(ierr)
    call get_global_attribute(fname, 'composite_constituent_id', ecckd_data%composite_constituent_id, ierr); call CHKERR(ierr)
    groups(1) = trim(fname)
    groups(2) = 'n_gases         '; call ncload(groups, ecckd_data%n_gases, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'temperature     '; call ncload(groups, ecckd_data%temperature, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'pressure        '; call ncload(groups, ecckd_data%pressure, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'wavenumber1     '; call ncload(groups, ecckd_data%wavenumber1, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'wavenumber2     '; call ncload(groups, ecckd_data%wavenumber2, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'gpoint_fraction '; call ncload(groups, ecckd_data%gpoint_fraction, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'wavenumber1_band'; call ncload(groups, ecckd_data%wavenumber1_band, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'wavenumber2_band'; call ncload(groups, ecckd_data%wavenumber2_band, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'band_number     '; call ncload(groups, ecckd_data%band_number, ierr, lverbose); call CHKERR(ierr)
    ecckd_data%n_g_pnt = size(ecckd_data%band_number)

    groups(2) = 'composite_conc_dependence_code'
    call ncload(groups, ecckd_data%composite_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'composite_molar_absorption_coeff'
    call ncload(groups, ecckd_data%composite_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'composite_mole_fraction'
    call ncload(groups, ecckd_data%composite_mole_fraction, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'h2o_conc_dependence_code'
    call ncload(groups, ecckd_data%h2o_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'h2o_mole_fraction'
    call ncload(groups, ecckd_data%h2o_mole_fraction, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'h2o_molar_absorption_coeff'
    call ncload(groups, ecckd_data%h2o_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)

    groups(2) = 'o3_conc_dependence_code'
    call ncload(groups, ecckd_data%o3_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'o3_molar_absorption_coeff'
    call ncload(groups, ecckd_data%o3_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)

    groups(2) = 'co2_conc_dependence_code'
    call ncload(groups, ecckd_data%co2_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'co2_molar_absorption_coeff'
    call ncload(groups, ecckd_data%co2_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)

    groups(2) = 'ch4_conc_dependence_code'
    call ncload(groups, ecckd_data%ch4_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'ch4_reference_mole_fraction'
    call ncload(groups, ecckd_data%ch4_reference_mole_fraction, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'ch4_molar_absorption_coeff'
    call ncload(groups, ecckd_data%ch4_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)

    groups(2) = 'n2o_conc_dependence_code'
    call ncload(groups, ecckd_data%n2o_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'n2o_reference_mole_fraction'
    call ncload(groups, ecckd_data%n2o_reference_mole_fraction, ierr, lverbose); call CHKERR(ierr)
    groups(2) = 'n2o_molar_absorption_coeff'
    call ncload(groups, ecckd_data%n2o_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)

    if (lsw) then
      groups(2) = 'solar_irradiance'
      call ncload(groups, ecckd_data%solar_irradiance, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'rayleigh_molar_scattering_coeff'
      call ncload(groups, ecckd_data%rayleigh_molar_scattering_coeff, ierr, lverbose); call CHKERR(ierr)

    end if

    if (llw) then
      groups(2) = 'temperature_planck'; call ncload(groups, ecckd_data%temperature_planck, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'planck_function'; call ncload(groups, ecckd_data%planck_function, ierr, lverbose); call CHKERR(ierr)

      groups(2) = 'cfc11_conc_dependence_code'
      call ncload(groups, ecckd_data%cfc11_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'cfc11_molar_absorption_coeff'
      call ncload(groups, ecckd_data%cfc11_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)

      groups(2) = 'cfc12_conc_dependence_code'
      call ncload(groups, ecckd_data%cfc12_conc_dependence_code, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'cfc12_molar_absorption_coeff'
      call ncload(groups, ecckd_data%cfc12_molar_absorption_coeff, ierr, lverbose); call CHKERR(ierr)
    end if

    if (get_arg(.false., lverbose)) then
      print *, 'n_gases                         ', ecckd_data%n_gases
      print *, 'gases                           ', trim(ecckd_data%constituent_id)
      print *, 'composite gases                 ', trim(ecckd_data%composite_constituent_id)
      print *, 'temperature                     ', shape(ecckd_data%temperature)
      print *, 'pressure                        ', shape(ecckd_data%pressure)
      print *, 'wavenumber1                     ', shape(ecckd_data%wavenumber1)
      print *, 'wavenumber2                     ', shape(ecckd_data%wavenumber2)
      print *, 'gpoint_fraction                 ', shape(ecckd_data%gpoint_fraction)
      print *, 'wavenumber1_band                ', shape(ecckd_data%wavenumber1_band)
      print *, 'wavenumber2_band                ', shape(ecckd_data%wavenumber2_band)
      print *, 'band_number                     ', shape(ecckd_data%band_number)
      print *, 'composite_conc_dependence_code  ', shape(ecckd_data%composite_conc_dependence_code)
      print *, 'composite_molar_absorption_coeff', shape(ecckd_data%composite_molar_absorption_coeff)
      print *, 'composite_mole_fraction         ', shape(ecckd_data%composite_mole_fraction)
      print *, 'h2o_conc_dependence_code        ', shape(ecckd_data%h2o_conc_dependence_code)
      print *, 'h2o_mole_fraction               ', shape(ecckd_data%h2o_mole_fraction)
      print *, 'h2o_molar_absorption_coeff      ', shape(ecckd_data%h2o_molar_absorption_coeff)
      print *, 'o3_conc_dependence_code         ', shape(ecckd_data%o3_conc_dependence_code)
      print *, 'o3_molar_absorption_coeff       ', shape(ecckd_data%o3_molar_absorption_coeff)
      print *, 'co2_conc_dependence_code        ', shape(ecckd_data%co2_conc_dependence_code)
      print *, 'co2_molar_absorption_coeff      ', shape(ecckd_data%co2_molar_absorption_coeff)
      print *, 'ch4_conc_dependence_code        ', shape(ecckd_data%ch4_conc_dependence_code)
      print *, 'ch4_reference_mole_fraction     ', shape(ecckd_data%ch4_reference_mole_fraction)
      print *, 'ch4_molar_absorption_coeff      ', shape(ecckd_data%ch4_molar_absorption_coeff)

      print *, 'n2o_conc_dependence_code      ', shape(ecckd_data%n2o_conc_dependence_code)
      print *, 'n2o_reference_mole_fraction   ', shape(ecckd_data%n2o_reference_mole_fraction)
      print *, 'n2o_molar_absorption_coeff    ', shape(ecckd_data%n2o_molar_absorption_coeff)

      if (lsw) then
      end if

      if (llw) then
        print *, 'temperature_planck            ', shape(ecckd_data%temperature_planck)
        print *, 'planck_function               ', shape(ecckd_data%planck_function)
        print *, 'cfc11_conc_dependence_code    ', shape(ecckd_data%cfc11_conc_dependence_code)
        print *, 'cfc11_molar_absorption_coeff  ', shape(ecckd_data%cfc11_molar_absorption_coeff)
        print *, 'cfc12_conc_dependence_code    ', shape(ecckd_data%cfc12_conc_dependence_code)
        print *, 'cfc12_molar_absorption_coeff  ', shape(ecckd_data%cfc12_molar_absorption_coeff)
      end if

    end if
  end subroutine

  subroutine ecckd_init(comm, atm, general_mie_table, general_fu_ice_table, ecckd_data_solar, ecckd_data_thermal, ierr, &
      & fname_ecckd_solar, fname_ecckd_thermal, lverbose)
    integer(mpiint), intent(in) :: comm
    type(t_tenstr_atm), intent(in), target :: atm
    type(t_mie_table), intent(in) :: general_mie_table
    type(t_fu_muskatel_ice_data), intent(in) :: general_fu_ice_table
    type(t_ecckd_data), allocatable, intent(inout) :: ecckd_data_solar, ecckd_data_thermal
    integer(mpiint), intent(out) :: ierr
    character(len=*), intent(in), optional :: fname_ecckd_solar, fname_ecckd_thermal
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: basepath, fname_thermal, fname_solar
    character(len=*), parameter :: fname_thermal_default = share_dir//"ecckd-1.0_lw_climate_fsck-16_ckd-definition.nc"
    character(len=*), parameter :: fname_solar_default = share_dir//"ecckd-1.0_sw_climate_rgb-16_ckd-definition.nc"
    logical :: lset, lexists
    integer(mpiint) :: myid

    ierr = 0

    call mpi_comm_rank(comm, myid, ierr)

    call setup_log_events(ecckd_log_events)

    if (myid .eq. 0) then
      fname_thermal = get_arg(fname_thermal_default, fname_ecckd_thermal)
      fname_solar = get_arg(fname_solar_default, fname_ecckd_solar)

      basepath = ''
      call get_petsc_opt('', '-ecckd_data', basepath, lset, ierr); call CHKERR(ierr)
      if (lset) then
        fname_thermal = trim(basepath)//'/'//trim(fname_thermal)
        fname_solar = trim(basepath)//'/'//trim(fname_solar)
      end if

      call get_petsc_opt('', '-ecckd_data_thermal', fname_thermal, lset, ierr); call CHKERR(ierr)
      call get_petsc_opt('', '-ecckd_data_solar', fname_solar, lset, ierr); call CHKERR(ierr)

      inquire (file=trim(fname_thermal), exist=lexists)
      if (.not. lexists) then
        call CHKERR(1_mpiint, "File at ecckd thermal path "//toStr(fname_thermal)// &
          & " does not exist"//new_line('')// &
          & " please make sure the file is at this location"// &
          & " you may set a directory with option "//new_line('')// &
          & "   -ecckd_data  <path with file ecckd_thermal.lut>"//new_line('')// &
          & " or specify a correct path with option"//new_line('')// &
          & "   -ecckd_data_thermal <path>")
      end if

      inquire (file=trim(fname_solar), exist=lexists)
      if (.not. lexists) then
        call CHKERR(1_mpiint, "File at ecckd solar path "//toStr(fname_solar)// &
          & " does not exist"//new_line('')// &
          & " please make sure the file is at this location"// &
          & " you may set a directory with option "//new_line('')// &
          & "   -ecckd_data  <path with file ecckd_solar.lut>"//new_line('')// &
          & " or specify a correct path with option"//new_line('')// &
          & "   -ecckd_data_solar <path>")
      end if

      if (get_arg(.false., lverbose)) print *, 'Reading ecckd thermal data from '//trim(fname_thermal)
      call load_data(fname_thermal, lsw=.false., llw=.true., ecckd_data=ecckd_data_thermal, ierr=ierr, lverbose=lverbose)
      call CHKERR(ierr)

      if (get_arg(.false., lverbose)) print *, 'Reading ecckd solar data from '//trim(fname_solar)
      call load_data(fname_solar, lsw=.true., llw=.false., ecckd_data=ecckd_data_solar, ierr=ierr, lverbose=lverbose)
      call CHKERR(ierr)
    end if

    call distribute_ecckd_table(comm, ecckd_data_thermal, ierr); call CHKERR(ierr)
    call distribute_ecckd_table(comm, ecckd_data_solar, ierr); call CHKERR(ierr)

    call populate_gas_info(atm, ecckd_data_thermal, ierr); call CHKERR(ierr)
    call populate_gas_info(atm, ecckd_data_solar, ierr); call CHKERR(ierr)

    call init_mie_table(general_mie_table, ecckd_data_thermal, ierr); call CHKERR(ierr)
    call init_mie_table(general_mie_table, ecckd_data_solar, ierr); call CHKERR(ierr)

    call init_fu_ice_table(general_fu_ice_table, ecckd_data_thermal, ierr); call CHKERR(ierr)
    call init_fu_ice_table(general_fu_ice_table, ecckd_data_solar, ierr); call CHKERR(ierr)
  contains

    subroutine distribute_ecckd_table(comm, table, ierr)
      integer(mpiint), intent(in) :: comm
      type(t_ecckd_data), allocatable, intent(inout) :: table
      integer(mpiint), intent(out) :: ierr
      logical :: lhave_thermal, lhave_solar
      integer(mpiint), parameter :: sendid = 0

      ierr = 0
      if (.not. allocated(table)) allocate (table)
      call imp_bcast(comm, table%constituent_id, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%composite_constituent_id, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%n_gases, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%n_g_pnt, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%temperature, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%pressure, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%wavenumber1, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%wavenumber2, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%gpoint_fraction, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%wavenumber1_band, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%wavenumber2_band, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%band_number, sendid, ierr); call CHKERR(ierr)

      call imp_bcast(comm, table%composite_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%composite_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%composite_mole_fraction, sendid, ierr); call CHKERR(ierr)

      call imp_bcast(comm, table%h2o_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%h2o_mole_fraction, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%h2o_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)

      call imp_bcast(comm, table%o3_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%o3_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)

      call imp_bcast(comm, table%co2_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%co2_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)

      call imp_bcast(comm, table%ch4_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%ch4_reference_mole_fraction, sendid, ierr); call CHKERR(ierr)
      call imp_bcast(comm, table%ch4_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)

      lhave_solar = allocated(table%solar_irradiance)
      call imp_bcast(comm, lhave_solar, sendid, ierr); call CHKERR(ierr)
      if (lhave_solar) then
        call imp_bcast(comm, table%solar_irradiance, sendid, ierr); call CHKERR(ierr)
        call imp_bcast(comm, table%rayleigh_molar_scattering_coeff, sendid, ierr); call CHKERR(ierr)

        call imp_bcast(comm, table%n2o_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
        call imp_bcast(comm, table%n2o_reference_mole_fraction, sendid, ierr); call CHKERR(ierr)
        call imp_bcast(comm, table%n2o_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)
      end if
      lhave_thermal = allocated(table%temperature_planck)
      call imp_bcast(comm, lhave_thermal, sendid, ierr); call CHKERR(ierr)
      if (lhave_thermal) then
        call imp_bcast(comm, table%temperature_planck, sendid, ierr); call CHKERR(ierr)
        call imp_bcast(comm, table%planck_function, sendid, ierr); call CHKERR(ierr)

        call imp_bcast(comm, table%cfc11_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
        call imp_bcast(comm, table%cfc11_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)

        call imp_bcast(comm, table%cfc12_conc_dependence_code, sendid, ierr); call CHKERR(ierr)
        call imp_bcast(comm, table%cfc12_molar_absorption_coeff, sendid, ierr); call CHKERR(ierr)
      end if
    end subroutine
  end subroutine

  subroutine init_fu_ice_table(general_fu_ice_table, ecckd_data, ierr)
    type(t_fu_muskatel_ice_data), intent(in) :: general_fu_ice_table
    type(t_ecckd_data), target, intent(inout) :: ecckd_data
    integer(mpiint), intent(out) :: ierr

    integer(iintegers) :: ireff, igpt, iwvnr
    real(ireals) :: reff, wgt, wvl_lo, wvl_hi, wvl, gpt_qext, gpt_w0, gpt_g, qext, w0, g
    ierr = 0

    if (allocated(ecckd_data%fu_ice_table)) return
    allocate (ecckd_data%fu_ice_table)

    allocate (ecckd_data%fu_ice_table%reff(size(general_fu_ice_table%reff)))
    ecckd_data%fu_ice_table%reff(:) = general_fu_ice_table%reff(:)

    allocate (ecckd_data%fu_ice_table%qext(size(general_fu_ice_table%reff), ecckd_data%n_g_pnt))
    allocate (ecckd_data%fu_ice_table%w0(size(general_fu_ice_table%reff), ecckd_data%n_g_pnt))
    allocate (ecckd_data%fu_ice_table%g(size(general_fu_ice_table%reff), ecckd_data%n_g_pnt))

    do ireff = 1, size(ecckd_data%fu_ice_table%reff)
      reff = ecckd_data%fu_ice_table%reff(ireff)

      do igpt = 1, size(ecckd_data%gpoint_fraction, dim=2)

        gpt_qext = 0
        gpt_w0 = 0
        gpt_g = 0

        do iwvnr = 1, size(ecckd_data%gpoint_fraction, dim=1)
          wgt = ecckd_data%gpoint_fraction(iwvnr, igpt)
          if (wgt .gt. 0) then
            wvl_lo = 1e7_ireals / ecckd_data%wavenumber2(iwvnr)
            wvl_hi = 1e7_ireals / ecckd_data%wavenumber1(iwvnr)
            wvl = (wvl_lo + wvl_hi)*.5

            call fu_ice_optprop(&
              & general_fu_ice_table, &
              & wvl * 1e-3_ireals, &
              & reff, &
              & qext, w0, g, ierr); call CHKERR(ierr)

            gpt_qext = gpt_qext + wgt * qext
            gpt_w0 = gpt_w0 + wgt * w0
            gpt_g = gpt_g + wgt * g
          end if
        end do

        ecckd_data%fu_ice_table%qext(ireff, igpt) = real(gpt_qext, irealLUT)
        ecckd_data%fu_ice_table%w0(ireff, igpt) = real(gpt_w0, irealLUT)
        ecckd_data%fu_ice_table%g(ireff, igpt) = real(gpt_g, irealLUT)
      end do
    end do

  end subroutine

  subroutine init_mie_table(general_mie_table, ecckd_data, ierr)
    type(t_mie_table), intent(in) :: general_mie_table
    type(t_ecckd_data), target, intent(inout) :: ecckd_data
    integer(mpiint), intent(out) :: ierr

    integer(iintegers) :: ireff, igpt, iwvnr
    real(ireals) :: reff, wgt, wvl_lo, wvl_hi, wvl, gpt_qext, gpt_w0, gpt_g, qext, w0, g
    ierr = 0

    if (allocated(ecckd_data%mie_table)) return
    allocate (ecckd_data%mie_table)

    allocate (ecckd_data%mie_table%reff(size(general_mie_table%reff)))
    ecckd_data%mie_table%reff(:) = general_mie_table%reff(:)

    allocate (ecckd_data%mie_table%qext(size(general_mie_table%reff), ecckd_data%n_g_pnt))
    allocate (ecckd_data%mie_table%w0(size(general_mie_table%reff), ecckd_data%n_g_pnt))
    allocate (ecckd_data%mie_table%g(size(general_mie_table%reff), ecckd_data%n_g_pnt))

    do ireff = 1, size(ecckd_data%mie_table%reff)
      reff = ecckd_data%mie_table%reff(ireff)

      do igpt = 1, size(ecckd_data%gpoint_fraction, dim=2)

        gpt_qext = 0
        gpt_w0 = 0
        gpt_g = 0

        do iwvnr = 1, size(ecckd_data%gpoint_fraction, dim=1)
          wgt = ecckd_data%gpoint_fraction(iwvnr, igpt)
          if (wgt .gt. 0) then
            wvl_lo = 1e7_ireals / ecckd_data%wavenumber2(iwvnr)
            wvl_hi = 1e7_ireals / ecckd_data%wavenumber1(iwvnr)
            wvl = (wvl_lo + wvl_hi)*.5

            call mie_optprop(&
              & general_mie_table, &
              & wvl * 1e-3_ireals, &
              & reff, &
              & qext, w0, g, ierr); call CHKERR(ierr)

            gpt_qext = gpt_qext + wgt * qext
            gpt_w0 = gpt_w0 + wgt * w0
            gpt_g = gpt_g + wgt * g
          end if
        end do

        ecckd_data%mie_table%qext(ireff, igpt) = real(gpt_qext, irealLUT)
        ecckd_data%mie_table%w0(ireff, igpt) = real(gpt_w0, irealLUT)
        ecckd_data%mie_table%g(ireff, igpt) = real(gpt_g, irealLUT)
      end do
    end do

  end subroutine

  subroutine populate_gas_info(atm, ecckd_data, ierr)
    type(t_tenstr_atm), target, intent(in) :: atm
    type(t_ecckd_data), target, intent(inout) :: ecckd_data
    integer(mpiint), intent(out) :: ierr
    character(len=default_str_len), allocatable :: gas_strings(:), composite_gas_strings(:)
    integer(iintegers) :: i, j
    ierr = 0

    allocate (ecckd_data%log_pressure(size(ecckd_data%pressure)))
    ecckd_data%log_pressure = log(ecckd_data%pressure)

    call split(ecckd_data%constituent_id, gas_strings, ' ', ierr); call CHKERR(ierr)
    allocate (ecckd_data%gases(size(gas_strings)))

    do i = 1, size(gas_strings)
      print *, 'gas string', i, trim(gas_strings(i))
      associate (gas => ecckd_data%gases(i))
        gas%id = trim(gas_strings(i))
        gas%reference_mole_fraction = -1
        nullify (gas%mole_fraction1)
        !nullify (gas%mole_fraction2)
        nullify (gas%molar_absorption_coeff3)
        nullify (gas%molar_absorption_coeff4)
        nullify (gas%vmr)
        gas%vmr_const = -1

        select case (trim(gas_strings(i)))
        case ('composite')
          call split(ecckd_data%composite_constituent_id, composite_gas_strings, ' ', ierr); call CHKERR(ierr)
          do j = 1, size(composite_gas_strings)
            print *, 'composite gas string', j, trim(composite_gas_strings(j))
          end do
          gas%conc_dependence_code = ecckd_data%composite_conc_dependence_code
          !gas%mole_fraction2 => ecckd_data%composite_mole_fraction
          gas%molar_absorption_coeff3 => ecckd_data%composite_molar_absorption_coeff

        case ('h2o')
          gas%conc_dependence_code = ecckd_data%h2o_conc_dependence_code
          gas%mole_fraction1 => ecckd_data%h2o_mole_fraction
          allocate (gas%log_mole_fraction1(size(gas%mole_fraction1)))
          gas%log_mole_fraction1 = log(gas%mole_fraction1)
          gas%molar_absorption_coeff4 => ecckd_data%h2o_molar_absorption_coeff
          gas%vmr => atm%h2o_lay

        case ('o3')
          gas%conc_dependence_code = ecckd_data%o3_conc_dependence_code
          gas%molar_absorption_coeff3 => ecckd_data%o3_molar_absorption_coeff
          gas%vmr => atm%o3_lay

        case ('co2')
          gas%conc_dependence_code = ecckd_data%co2_conc_dependence_code
          gas%molar_absorption_coeff3 => ecckd_data%co2_molar_absorption_coeff
          gas%vmr => atm%co2_lay

        case ('ch4')
          gas%conc_dependence_code = ecckd_data%ch4_conc_dependence_code
          gas%molar_absorption_coeff3 => ecckd_data%ch4_molar_absorption_coeff
          gas%reference_mole_fraction = ecckd_data%ch4_reference_mole_fraction
          gas%vmr => atm%ch4_lay

        case ('n2o')
          gas%conc_dependence_code = ecckd_data%n2o_conc_dependence_code
          gas%molar_absorption_coeff3 => ecckd_data%n2o_molar_absorption_coeff
          gas%reference_mole_fraction = ecckd_data%n2o_reference_mole_fraction
          gas%vmr => atm%n2o_lay

        case ('cfc11')
          gas%conc_dependence_code = ecckd_data%cfc11_conc_dependence_code
          gas%molar_absorption_coeff3 => ecckd_data%cfc11_molar_absorption_coeff
          gas%vmr_const = CFC11

        case ('cfc12')
          gas%conc_dependence_code = ecckd_data%cfc12_conc_dependence_code
          gas%molar_absorption_coeff3 => ecckd_data%cfc12_molar_absorption_coeff
          gas%vmr_const = CFC12

        case default
          call CHKERR(1_mpiint, 'dont know how to handle gas string '//trim(gas_strings(i)))
        end select

        select case (gas%conc_dependence_code)
        case (IConcDependenceNone, IConcDependenceLinear, IConcDependenceLUT, IConcDependenceRelativeLinear)
          continue
        case default
          call CHKWARN(int(gas%conc_dependence_code, mpiint), 'Gas concentration dependence code '// &
            & toStr(gas%conc_dependence_code)// &
            & ' is unknown. The gas '//trim(gas_strings(i))//' will be ignored')
        end select
      end associate
    end do
  end subroutine

  subroutine setup_log_events(logs, solvername)
    type(t_ecckd_log_events), intent(inout) :: logs
    character(len=*), optional :: solvername
    character(len=default_str_len) :: s
    PetscClassId :: cid
    integer(mpiint) :: ierr

    s = get_arg('tenstr_ecckd.', solvername)
    call PetscClassIdRegister(trim(s), cid, ierr); call CHKERR(ierr)

    call setup_stage(trim(s)//'ecckd_solar', logs%stage_ecckd_solar)
    call setup_stage(trim(s)//'ecckd_thermal', logs%stage_ecckd_thermal)

    call PetscLogEventRegister(trim(s)//'ecckd_optprop', cid, logs%ecckd_optprop, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'ecckd_optprop_dtau', cid, logs%ecckd_optprop_dtau, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'ecckd_optprop_rayleigh', cid, logs%ecckd_optprop_rayleigh, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'ecckd_optprop_mie', cid, logs%ecckd_optprop_mie, ierr); call CHKERR(ierr)
    call PetscLogEventRegister(trim(s)//'ecckd_optprop_fu_ice', cid, logs%ecckd_optprop_fu_ice, ierr); call CHKERR(ierr)

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
