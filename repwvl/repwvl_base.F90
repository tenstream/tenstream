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

  use m_data_parameters, only: &
    & AVOGADRO, &
    & default_str_len, &
    & EARTHACCEL, &
    & iintegers, &
    & ireals, &
    & K_BOLTZMANN, &
    & MOLMASSAIR, &
    & mpiint, &
    & R_DRY_AIR

  use m_helper_functions, only: &
    & CHKERR, &
    & get_petsc_opt, &
    & get_arg, &
    & toStr

  use m_search, only: find_real_location

  use m_netcdfIO, only: ncload

  implicit none

  private
  public :: t_repwvl_data, repwvl_init, repwvl_dtau

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

  logical, parameter :: ldebug = .true.

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

  subroutine repwvl_init(repwvl_data_solar, repwvl_data_thermal, ierr, fname_repwvl_solar, fname_repwvl_thermal, lverbose)
    type(t_repwvl_data), allocatable, intent(inout) :: repwvl_data_solar, repwvl_data_thermal
    integer(mpiint), intent(out) :: ierr
    character(len=*), intent(in), optional :: fname_repwvl_solar, fname_repwvl_thermal
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: basepath, fname_thermal, fname_solar
    logical :: lset, lexists

    ierr = 0

    fname_thermal = get_arg('repwvl_thermal.lut', fname_repwvl_thermal)
    fname_solar = get_arg('repwvl_solar.lut', fname_repwvl_solar)

    basepath = ''
    call get_petsc_opt('', '-repwvl_data', basepath, lset, ierr); call CHKERR(ierr)
    fname_thermal = trim(basepath)//trim(fname_thermal)
    fname_solar = trim(basepath)//trim(fname_solar)

    call get_petsc_opt('', '-repwvl_thermal_lut', fname_thermal, lset, ierr); call CHKERR(ierr)
    call get_petsc_opt('', '-repwvl_solar_lut', fname_solar, lset, ierr); call CHKERR(ierr)

    inquire (file=trim(fname_thermal), exist=lexists)
    if (.not. lexists) then
      call CHKERR(1_mpiint, "File at repwvl thermal path "//toStr(fname_thermal)// &
        & " does not exist"//new_line('')// &
        & " please make sure the file is at this location"// &
        & " you may set a directory with option "//new_line('')// &
        & "   -repwvl_data  <path with file repwvl_thermal.lut>"//new_line('')// &
        & " or specify a correct path with option"//new_line('')// &
        & "   -repwvl_thermal_lut <path>")
    end if

    inquire (file=trim(fname_solar), exist=lexists)
    if (.not. lexists) then
      call CHKERR(1_mpiint, "File at repwvl solar path "//toStr(fname_solar)// &
        & " does not exist"//new_line('')// &
        & " please make sure the file is at this location"// &
        & " you may set a directory with option "//new_line('')// &
        & "   -repwvl_data  <path with file repwvl_solar.lut>"//new_line('')// &
        & " or specify a correct path with option"//new_line('')// &
        & "   -repwvl_solar_lut <path>")
    end if

    if (get_arg(.false., lverbose)) print *, 'Reading representative wavelength data from '//trim(fname_thermal)
    call load_data(fname_thermal, repwvl_data_thermal, ierr, lverbose=lverbose); call CHKERR(ierr)

    if (get_arg(.false., lverbose)) print *, 'Reading representative wavelength data from '//trim(fname_solar)
    call load_data(fname_solar, repwvl_data_solar, ierr, lverbose=lverbose); call CHKERR(ierr)
  end subroutine

  subroutine repwvl_dtau(    &
      & repwvl_data,         & !> data tables
      & iwvl,                & !> wavelenght index
      & P,                   & !> mean pressure of layer [Pa]
      & dP,                  & !> delta pressure of layer [Pa]
      & T,                   & !> mean Temperature of layer [K]
      & VMRS,                & !> vol mix ratios as in the repwvl data, e.g. (H2O, H2O, CO2, O3, N2O, CO, CH4, O2, HNO3, N2)
      & dtau,                & !> optical thickness of layer
      & ierr)

    type(t_repwvl_data), intent(in) :: repwvl_data
    integer(iintegers), intent(in) :: iwvl
    real(ireals), intent(in) :: P, dP, T
    real(ireals), intent(in) :: VMRS(:)
    real(ireals), intent(out) :: dtau
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: numDens, wP, wT
    real(ireals) :: wgt_P, wgt_T, xsecP0, xsecP1
    integer(iintegers) :: ip0, ip1, iT0, iT1, spec
    ierr = 0

    if (any([T, P, dP] .lt. 0)) then
      ierr = 1
      call CHKERR(ierr, "Found bad input, one of variables T,p,dP < 0:"//new_line('')// &
        & " T = "//toStr(T)//new_line('')// &
        & " p = "//toStr(P)//new_line('')// &
        & " dP= "//toStr(dP)//new_line('')// &
        & "")
    end if

    if (any(VMRS .lt. 0) .or. any(VMRS .gt. 1)) then
      ierr = 1
      call CHKERR(ierr, "Found bad input, one of VMRS < 0 or > 1:"//new_line('')// &
        & " VMRS = "//toStr(VMRS)//new_line('')// &
        & "")
    end if

    numDens = dP * AVOGADRO / MOLMASSAIR / EARTHACCEL

    wP = find_real_location(repwvl_data%p_ref, P)
    ip0 = int(floor(wP), iintegers)
    ip1 = ip0 + 1

    wT = find_real_location(repwvl_data%t_pert + repwvl_data%t_ref(ip0), T)
    iT0 = int(floor(wT), iintegers)
    iT1 = iT0 + 1

    wgt_P = wP - real(ip0, ireals)
    wgt_T = wT - real(iT0, ireals)

    dtau = 0

    ! H2O continuum is quadratic
    spec = 1
    xsecP0 = (repwvl_data%xsec(ip0, iwvl, spec, iT1) * wgt_T &
           & + repwvl_data%xsec(ip0, iwvl, spec, iT0) * (1._ireals - wgt_T)) &
           & / repwvl_data%vmrs_ref(ip0, spec)
    xsecP1 = (repwvl_data%xsec(ip1, iwvl, spec, iT1) * wgt_T &
           & + repwvl_data%xsec(ip1, iwvl, spec, iT0) * (1._ireals - wgt_T)) &
           & / repwvl_data%vmrs_ref(ip1, spec)
    dtau = dtau + (xsecP1 * wgt_P + xsecP0 * (1._ireals - wgt_P)) * numDens * (VMRS(spec) * VMRS(spec))

    do spec = 2, size(VMRS)
      xsecP0 = (repwvl_data%xsec(ip0, iwvl, spec, iT1) * wgt_T &
             & + repwvl_data%xsec(ip0, iwvl, spec, iT0) * (1._ireals - wgt_T))
      xsecP1 = (repwvl_data%xsec(ip1, iwvl, spec, iT1) * wgt_T &
             & + repwvl_data%xsec(ip1, iwvl, spec, iT0) * (1._ireals - wgt_T))
      dtau = dtau + (xsecP1 * wgt_P + xsecP0 * (1._ireals - wgt_P)) * numDens * VMRS(spec)
    end do

    if (allocated(repwvl_data%crs_o3)) &
      & call additional_tracer_from_bremen_data(repwvl_data%crs_o3(:, iwvl), VMRS(4), dtau)
    if (allocated(repwvl_data%crs_n2o)) &
      & call additional_tracer_from_bremen_data(repwvl_data%crs_n2o(:, iwvl), VMRS(5), dtau)

    if (dtau .lt. 0) then
      ierr = 2
      call CHKERR(ierr, "Resulted in negative dtau ?! dtau = "//toStr(dtau))
    end if

  contains
    subroutine additional_tracer_from_bremen_data(coeff, vmr, tau)
      real(ireals), intent(in) :: coeff(:), vmr
      real(ireals), intent(inout) :: tau

      real(ireals), parameter :: T0 = 273.15
      real(ireals) :: rho, dz, N, sigma, dtau

      rho = P / (R_DRY_AIR * T)
      dz = dP / (rho * EARTHACCEL)
      N = P / (K_BOLTZMANN * T) * 1e-4_ireals * dz
      sigma = max(0._ireals, (coeff(1) + coeff(2) * (T - T0) + coeff(3) * (T - T0)**2) * 1e-20_ireals)
      dtau = vmr * N * sigma
      tau = tau + dtau
    end subroutine
  end subroutine
end module
