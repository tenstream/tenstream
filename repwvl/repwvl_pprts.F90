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

module m_repwvl_pprts
#include "petsc/finclude/petsc.h"
  use petsc

  use m_helper_functions, only: &
    & approx, &
    & CHKERR, &
    & cross_3d, &
    & deg2rad, &
    & get_arg, &
    & get_petsc_opt, &
    & gradient, &
    & imp_allreduce_max, &
    & imp_allreduce_mean, &
    & imp_allreduce_min, &
    & imp_bcast, &
    & ind_1d_to_nd, &
    & meanvec, &
    & mpi_logical_all_same, &
    & read_ascii_file_2d, &
    & reverse, &
    & spherical_2_cartesian, &
    & toStr

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & iintegers, ireals, mpiint, &
    & zero, one, default_str_len, &
    & i0, i1, &
    & AVOGADRO, &
    & EARTHACCEL, &
    & MOLMASSAIR

  use m_tenstream_options, only: read_commandline_options

  use m_pprts_base, only: t_solver, destroy_pprts
  use m_pprts, only: &
    & init_pprts, &
    & pprts_get_result, &
    & pprts_get_result_toZero, &
    & set_angles, &
    & set_optical_properties, &
    & solve_pprts

  use m_dyn_atm_to_rrtmg, only: &
    & planck, &
    & plkint, &
    & print_tenstr_atm, &
    & t_tenstr_atm, &
    & vert_integral_coeff

  use m_buildings, only: &
    & t_pprts_buildings, &
    & PPRTS_BOT_FACE

  use m_repwvl_base, only: repwvl_init, t_repwvl_data, repwvl_dtau
  use m_mie_tables, only: mie_tables_init, mie_water_table, mie_optprop
  use m_fu_ice, only: fu_ice_init, fu_ice_optprop, fu_ice_data_solar, fu_ice_data_thermal
  use m_rayleigh, only: rayleigh

  implicit none

  private
  public :: repwvl_pprts, repwvl_pprts_destroy

  logical, parameter :: ldebug = .true.

  type(t_repwvl_data), allocatable :: repwvl_data_solar, repwvl_data_thermal

  real(ireals), parameter :: CO = 1e-9_ireals
  real(ireals), parameter :: HNO3 = 1e-9_ireals
  real(ireals), parameter :: N2 = 0.7808_ireals

contains

  subroutine repwvl_pprts(comm, solver, atm, ie, je,   &
      & dx, dy, sundir,                               &
      & albedo_thermal, albedo_solar,                 &
      & lthermal, lsolar,                             &
      & edir, edn, eup, abso,                         &
      & nxproc, nyproc, icollapse,                    &
      & opt_time, solar_albedo_2d, thermal_albedo_2d, &
      & opt_solar_constant,                           &
      & opt_buildings_solar, opt_buildings_thermal,   &
      & opt_tau_solar,                                &
      & opt_w0_solar,                                 &
      & opt_g_solar,                                  &
      & opt_tau_thermal,                              &
      & lonly_initialize)

    integer(mpiint), intent(in) :: comm ! MPI Communicator

    class(t_solver), intent(inout) :: solver                       ! solver type (e.g. t_solver_8_10)
    type(t_tenstr_atm), intent(in) :: atm                          ! contains info on atmospheric constituents
    integer(iintegers), intent(in) :: ie, je                       ! local domain size in x and y direction
    real(ireals), intent(in) :: dx, dy                       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: sundir(:)                    ! cartesian sun direction, pointing away from the sun, dim(3)
    real(ireals), intent(in) :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! nxproc dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    ! nyproc dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    ! if not present, we let petsc decide how to decompose the fields(probably does not fit the decomposition of a host model)
    integer(iintegers), intent(in), optional :: nxproc(:), nyproc(:)

    integer(iintegers), intent(in), optional :: icollapse ! experimental, dont use it if you dont know what you are doing.

    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions
    ! and compute new solutions only after threshold estimate is exceeded.
    ! If solar_albedo_2d is present, we use a 2D surface albedo
    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:, :), thermal_albedo_2d(:, :), opt_solar_constant

    ! buildings information, setup broadband thermal and solar albedo on faces inside the domain, not just on the surface
    ! see definition for details on how to set it up
    type(t_pprts_buildings), intent(inout), optional :: opt_buildings_solar
    type(t_pprts_buildings), intent(inout), optional :: opt_buildings_thermal

    ! optional optical properties with dims(nlyr(srfc to TOA), local_nx, local_ny, nr-g-points)
    ! used e.g. for aerosol or vegetation, you can provide only tau or (tau and w0) defaults to w0=0 and g=0
    ! note that first dimension can also be smaller than nlyr, we will fill it up from the ground,
    ! i.e. if you provide only two layers, the two lowermost layers near the surface will be filled with addition optprops
    real(ireals), intent(in), optional, dimension(:, :, :, :) :: opt_tau_solar, opt_w0_solar, opt_g_solar
    real(ireals), intent(in), optional, dimension(:, :, :, :) :: opt_tau_thermal

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :), intent(inout) :: edir, edn, eup  ! [nlyr+1, local_nx, local_ny ]
    real(ireals), allocatable, dimension(:, :, :), intent(inout) :: abso            ! [nlyr  , local_nx, local_ny ]

    ! if only_initialize we dont compute any radiation, merely setup the grid structures
    logical, intent(in), optional :: lonly_initialize

    ! ---------- end of API ----------------

    ! Counters
    integer(iintegers) :: ke, ke1

    ! for debug purposes, can output variables into netcdf files
    !character(default_str_len) :: output_path(2) ! [ filename, varname ]
    !logical :: lfile_exists

    integer(mpiint) :: myid, ierr
    logical :: lskip_thermal, lskip_solar, lprint_atm, lflg

    integer(iintegers) :: pprts_icollapse

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (ldebug) then ! make sure that all ranks give the same option for lsolar and lthermal
      if (.not. mpi_logical_all_same(comm, lsolar)) &
        call CHKERR(1_mpiint, 'all ranks have to give the same value for lsolar')
      if (.not. mpi_logical_all_same(comm, lthermal)) &
        call CHKERR(1_mpiint, 'all ranks have to give the same value for lthermal')
    end if

    if (.not. solver%linitialized) call read_commandline_options(comm) ! so that tenstream.options file are read in

    pprts_icollapse = get_arg(i1, icollapse)
    call get_petsc_opt(PETSC_NULL_CHARACTER, &
                       "-pprts_collapse", pprts_icollapse, lflg, ierr); call CHKERR(ierr)
    if (pprts_icollapse .eq. -1) then
      if (ldebug .and. myid .eq. 0) print *, 'Collapsing background atmosphere', atm%atm_ke
      pprts_icollapse = atm%atm_ke ! collapse the complete background atmosphere
    end if
    if (pprts_icollapse .ne. i1) then
      if (ldebug .and. myid .eq. 0) print *, 'Collapsing atmosphere', pprts_icollapse
    end if

    ke1 = ubound(atm%plev, 1)
    ke = ubound(atm%tlay, 1)

    if (.not. solver%linitialized) then
      call repwvl_init(repwvl_data_solar, repwvl_data_thermal, ierr); call CHKERR(ierr)
      call mie_tables_init(comm, ierr, lverbose=.false.); call CHKERR(ierr)
      call fu_ice_init(comm, ierr, lverbose=.true.); call CHKERR(ierr)
      call init_pprts_repwvl(comm, solver, &
                             dx, dy, atm%dz, &
                             sundir, &
                             ie, je, ke, &
                             nxproc, nyproc, pprts_icollapse)
    end if

    ! Allocate space for results -- for integrated values...
    if (.not. allocated(edn)) allocate (edn(solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
    if (.not. allocated(eup)) allocate (eup(solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
    if (.not. allocated(abso)) allocate (abso(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    edn = 0
    eup = 0
    abso = 0

    lprint_atm = ldebug
    call get_petsc_opt(PETSC_NULL_CHARACTER, &
                       "-repwvl_pprts_atm_view", lprint_atm, lflg, ierr); call CHKERR(ierr)
    if (lprint_atm .and. myid .eq. 0) then
      call print_tenstr_atm(atm)
    end if

    if (get_arg(.false., lonly_initialize)) return

    lskip_thermal = .false.
    call get_petsc_opt(PETSC_NULL_CHARACTER, &
                       "-skip_thermal", lskip_thermal, lflg, ierr); call CHKERR(ierr)
    if (lthermal .and. .not. lskip_thermal) then
      call compute_thermal(                    &
        & repwvl_data_thermal,                 &
        & solver,                              &
        & atm,                                 &
        & albedo_thermal,                      &
        & edn, eup, abso,                      &
        & ierr,                                &
        & opt_time=opt_time,                   &
        & thermal_albedo_2d=thermal_albedo_2d, &
        & opt_buildings=opt_buildings_thermal, &
        & opt_tau=opt_tau_thermal              &
        & )
    end if

    if (get_arg(.false., lskip_thermal) .and. get_arg(.false., lskip_solar)) then !DEBUG
      call CHKERR(1_mpiint, 'DEBUG this block just to get rid of unused warnings')
      ke = 0
      ke1 = 1
      print *, albedo_solar + albedo_thermal + dx + dy + sundir(1) + edir + icollapse + ie + je + ke1, nxproc + nyproc
      print *, opt_g_solar, opt_solar_constant, opt_tau_solar, opt_tau_thermal
      print *, opt_time, opt_w0_solar, pprts_icollapse, solar_albedo_2d, thermal_albedo_2d
      print *, opt_buildings_solar%iface
      print *, opt_buildings_thermal%iface
    end if

    if (lsolar .and. .not. allocated(edir)) allocate (edir(solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
    if (allocated(edir)) edir = zero

    if (lsolar) then
      lskip_solar = .false.
      call get_petsc_opt(PETSC_NULL_CHARACTER, &
                         "-skip_solar", lskip_solar, lflg, ierr); call CHKERR(ierr)
      if (.not. lskip_solar) then
        call compute_solar(                          &
          & repwvl_data_solar,                       &
          & solver,                                  &
          & atm,                                     &
          & sundir, albedo_solar,                    &
          & edir, edn, eup, abso,                    &
          & ierr,                                    &
          & opt_time=opt_time,            &
          & solar_albedo_2d=solar_albedo_2d,     &
          & opt_solar_constant=opt_solar_constant,  &
          & opt_buildings=opt_buildings_solar, &
          & opt_tau=opt_tau_solar,       &
          & opt_w0=opt_w0_solar,        &
          & opt_g=opt_g_solar          &
          & )
      end if
    end if
  end subroutine

  subroutine compute_thermal( &
      & repwvl_data_thermal,  &
      & solver,               &
      & atm,                  &
      & albedo,               &
      & edn,                  &
      & eup,                  &
      & abso,                 &
      & ierr,                 &
      & opt_time,             &
      & thermal_albedo_2d,    &
      & opt_buildings,        &
      & opt_tau)

    type(t_repwvl_data), intent(in) :: repwvl_data_thermal
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in), target :: atm
    real(ireals), intent(in) :: albedo
    real(ireals), intent(inout), dimension(:, :, :) :: edn, eup, abso
    integer(mpiint), intent(out) :: ierr

    real(ireals), optional, intent(in) :: opt_time, thermal_albedo_2d(:, :)
    type(t_pprts_buildings), intent(inout), optional :: opt_buildings
    real(ireals), intent(in), optional, dimension(:, :, :, :) :: opt_tau

    integer(iintegers) :: iwvl, i, j, k, icol

    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, kg ! [nlyr, local_nx, local_ny]
    real(ireals), allocatable :: Blev(:, :, :)

    real(ireals), dimension(:, :, :), allocatable :: spec_edn, spec_eup, spec_abso

    ierr = 0

    if (present(opt_buildings)) call CHKERR(1_mpiint, 'opt_buildings not yet implemented for repwvl_thermal')
    if (present(opt_tau)) call CHKERR(1_mpiint, 'opt_tau not yet implemented for repwvl_thermal')
    if (present(opt_time)) call CHKERR(1_mpiint, 'opt_time not yet implemented for repwvl_thermal')
    if (present(thermal_albedo_2d)) call CHKERR(1_mpiint, 'thermal_albedo_2d not yet implemented for repwvl_thermal')

    allocate (kabs(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    allocate (ksca(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    allocate (kg(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    allocate (Blev(solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))

    do iwvl = 1, size(repwvl_data_thermal%wvls)
      if (ldebug) then
        print *, 'Computing wavelengths '//toStr(iwvl)//' / '//toStr(size(repwvl_data_thermal%wvls))//&
          & ' -- '//toStr(100._ireals * real(iwvl, ireals) / real(size(repwvl_data_thermal%wvls), ireals))//' %'// &
          & ' ('//toStr(repwvl_data_thermal%wvls(iwvl))//' nm,  wgt='//toStr(repwvl_data_thermal%wgts(iwvl))//')'
      end if

      call repwvl_optprop(repwvl_data_thermal, atm, .false., iwvl, kabs, ksca, kg, ierr); call CHKERR(ierr)

      do j = 1, solver%C_one%ym
        do i = 1, solver%C_one%xm
          icol = i + (j - 1_iintegers) * solver%C_one%xm
          do k = 1, solver%C_one1%zm
            Blev(solver%C_one1%zm + 1 - k, i, j) = repwvl_data_thermal%wgts(iwvl) &
                                                  & * planck(repwvl_data_thermal%wvls(iwvl) * 1e-9_ireals, atm%tlev(k, icol)) &
                                                  & * 1e-9_ireals
          end do
        end do
      end do

      call set_optical_properties( &
        & solver,                  &
        & albedo,                  &
        & kabs,                    &
        & ksca,                    &
        & kg,                      &
        & planck=Blev)

      call solve_pprts(solver, lthermal=.true., lsolar=.false., edirTOA=-1._ireals)

      call pprts_get_result(solver, spec_edn, spec_eup, spec_abso)

      edn = edn + spec_edn
      eup = eup + spec_eup
      abso = abso + spec_abso
    end do !iwvl
  end subroutine

  subroutine compute_solar( &
      & repwvl_data_solar,  &
      & solver,             &
      & atm,                &
      & sundir,             &
      & albedo,             &
      & edir,               &
      & edn,                &
      & eup,                &
      & abso,               &
      & ierr,               &
      & opt_time,           &
      & solar_albedo_2d,    &
      & opt_solar_constant, &
      & opt_buildings,      &
      & opt_tau,            &
      & opt_w0,             &
      & opt_g)

    type(t_repwvl_data), intent(in) :: repwvl_data_solar
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in), target :: atm

    real(ireals), intent(in) :: albedo
    real(ireals), intent(in) :: sundir(3)

    real(ireals), intent(inout), dimension(:, :, :) :: edir, edn, eup, abso
    integer(mpiint), intent(out) :: ierr

    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:, :)
    real(ireals), intent(in), optional :: opt_solar_constant
    type(t_pprts_buildings), intent(inout), optional :: opt_buildings
    real(ireals), intent(in), optional, dimension(:, :, :, :) :: opt_tau, opt_w0, opt_g

    real(ireals) :: edirTOA

    integer(iintegers) :: iwvl

    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, kg      ! [nlyr, local_nx, local_ny]

    real(ireals), dimension(:, :, :), allocatable :: spec_edir, spec_edn, spec_eup, spec_abso

    ierr = 0

    if (present(opt_buildings)) call CHKERR(1_mpiint, 'opt_buildings not yet implemented for repwvl_solar')
    if (present(opt_tau)) call CHKERR(1_mpiint, 'opt_tau not yet implemented for repwvl_solar')
    if (present(opt_w0)) call CHKERR(1_mpiint, 'opt_w0 not yet implemented for repwvl_solar')
    if (present(opt_g)) call CHKERR(1_mpiint, 'opt_g not yet implemented for repwvl_solar')
    if (present(opt_time)) call CHKERR(1_mpiint, 'opt_time not yet implemented for repwvl_solar')
    if (present(solar_albedo_2d)) call CHKERR(1_mpiint, 'solar_albedo_2d not yet implemented for repwvl_solar')
    if (present(opt_solar_constant)) call CHKERR(1_mpiint, 'opt_solar_constant not yet implemented for repwvl_solar')

    allocate (kabs(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    allocate (ksca(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    allocate (kg(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

    call set_angles(solver, sundir)

    do iwvl = 1, size(repwvl_data_solar%wvls)
      if (ldebug) then
        print *, 'Computing wavelengths '//toStr(iwvl)//' / '//toStr(size(repwvl_data_solar%wvls))//&
          & ' -- '//toStr(100._ireals * real(iwvl, ireals) / real(size(repwvl_data_solar%wvls), ireals))//' %'// &
          & ' ('//toStr(repwvl_data_solar%wvls(iwvl))//' nm,  wgt='//toStr(repwvl_data_solar%wgts(iwvl))//')'
      end if

      call repwvl_optprop(repwvl_data_solar, atm, .true., iwvl, kabs, ksca, kg, ierr); call CHKERR(ierr)

      !call add_optional_optprop(tau, w0, g, opt_tau, opt_w0, opt_g)

      edirTOA = repwvl_data_solar%wgts(iwvl)

      call set_optical_properties( &
        & solver,                  &
        & albedo,                  &
        & kabs,                    &
        & ksca,                    &
        & kg)
      call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=edirTOA)

      call pprts_get_result(solver, spec_edn, spec_eup, spec_abso, spec_edir)

      edir = edir + spec_edir
      edn = edn + spec_edn
      eup = eup + spec_eup
      abso = abso + spec_abso
    end do !iwvl

  end subroutine compute_solar

  subroutine repwvl_optprop(repwvl_data, atm, lsolar, iwvl, kabs, ksca, kg, ierr)
    type(t_repwvl_data), intent(in) :: repwvl_data
    type(t_tenstr_atm), intent(in), target :: atm
    logical, intent(in) :: lsolar
    integer(iintegers) :: iwvl
    real(ireals), dimension(:, :, :) :: kabs, ksca, kg
    integer(mpiint), intent(out) :: ierr

    integer(iintegers) :: i, j, k, icol
    real(ireals) :: VMRS(repwvl_data%Ntracer)
    real(ireals) :: tabs, tsca, g
    real(ireals) :: P, dP, dtau, rayleigh_xsec, N, lwc_vmr, qext_cld, w0_cld, g_cld, iwp
    ierr = 0

    do j = lbound(kabs, 3), ubound(kabs, 3)
      do i = lbound(kabs, 2), ubound(kabs, 2)
        icol = i + (j - 1) * size(kabs, dim=2)
        do k = lbound(kabs, 1), ubound(kabs, 1)

          P = (atm%plev(k, icol) + atm%plev(k + 1, icol))*.5_ireals * 1e2_ireals
          dP = (atm%plev(k, icol) - atm%plev(k + 1, icol)) * 1e2_ireals

          VMRS(:) = [ &
            & atm%h2o_lay(k, icol), &
            & atm%h2o_lay(k, icol), &
            & atm%co2_lay(k, icol), &
            & atm%o3_lay(k, icol), &
            & atm%n2o_lay(k, icol), &
            & CO, &
            & atm%ch4_lay(k, icol), &
            & atm%o2_lay(k, icol), &
            & HNO3, &
            & N2]

          call repwvl_dtau(&
            & repwvl_data, &
            & iwvl, &
            & P, &
            & dP, &
            & atm%tlay(k, icol), &
            & VMRS, &
            & dtau, &
            & ierr); call CHKERR(ierr)

          tabs = dtau / atm%dz(k, icol)
          if (ldebug) then
            if (tabs .lt. 0) call CHKERR(1_mpiint, 'kabs from repwvl negative!'//toStr(tabs))
          end if

          call rayleigh(&
            & repwvl_data%wvls(iwvl) * 1e-3_ireals, &
            & atm%co2_lay(k, icol), &
            & rayleigh_xsec, &
            & ierr); call CHKERR(ierr)
          if (ldebug) then
            if (rayleigh_xsec .lt. 0) call CHKERR(1_mpiint, 'rayleigh xsec negative!'//toStr(rayleigh_xsec))
          end if

          N = dP * AVOGADRO / EARTHACCEL / MOLMASSAIR
          tsca = N * rayleigh_xsec * 1e-4 / atm%dz(k, icol) ! [1e-4 from cm2 to m2]
          if (ldebug) then
            if (tsca .lt. 0) call CHKERR(1_mpiint, 'rayleigh scattering coeff negative!'//toStr(tsca))
          end if

          ! rayleigh has symmetric asymmetry parameter
          g = 0

          if (atm%lwc(k, icol) > 0) then
            call mie_optprop(&
              & mie_water_table, &
              & repwvl_data%wvls(iwvl) * 1e-3_ireals, &
              & atm%reliq(k, icol), &
              & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)

            lwc_vmr = atm%lwc(k, icol) * dP / (EARTHACCEL * atm%dz(k, icol)) ! have lwc in [ g / kg ], lwc_vmr in [ g / m3 ]
            qext_cld = qext_cld * 1e-3 * lwc_vmr ! from [km^-1 / (g / m^3)] to [1/m]

            g = (g * tsca + g_cld * qext_cld) / (tsca + qext_cld)
            tabs = tabs + qext_cld * max(0._ireals, (1._ireals - w0_cld))
            tsca = tsca + qext_cld * w0_cld
          end if

          if (atm%iwc(k, icol) > 0) then
            if (lsolar) then
              call fu_ice_optprop(&
                & fu_ice_data_solar, &
                & iwvl, &
                & atm%reice(k, icol), &
                & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
            else
              call fu_ice_optprop(&
                & fu_ice_data_thermal, &
                & iwvl, &
                & atm%reice(k, icol), &
                & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
            end if

            iwp = atm%iwc(k, icol) * dP / (EARTHACCEL * atm%dz(k, icol)) ! have iwc in [ g / kg ], iwp in [ g / m3 ]
            qext_cld = qext_cld * iwp ! from [m^-1 / (g / m^3)] to [1/m]

            g = (g * tsca + g_cld * qext_cld) / (tsca + qext_cld)
            tabs = tabs + qext_cld * max(0._ireals, (1._ireals - w0_cld))
            tsca = tsca + qext_cld * w0_cld
          end if

          kabs(size(kabs, dim=1) + 1 - k, i, j) = tabs
          ksca(size(ksca, dim=1) + 1 - k, i, j) = tsca
          kg(size(kg, dim=1) + 1 - k, i, j) = g

        end do
      end do
    end do
  end subroutine

  subroutine init_pprts_repwvl(comm, solver, dx, dy, dz, &
                               sundir, &
                               xm, ym, zm, &
                               nxproc, nyproc, &
                               pprts_icollapse)

    integer(mpiint), intent(in) :: comm

    real(ireals), intent(in) :: dx, dy, sundir(3)
    real(ireals), intent(in) :: dz(:, :) ! bot to top, e.g. from atm%dz
    integer(iintegers), intent(in) :: xm, ym, zm
    class(t_solver), intent(inout) :: solver

    ! arrays containing xm and ym for all nodes :: dim[x-ranks, y-ranks]
    integer(iintegers), intent(in), optional :: nxproc(:), nyproc(:), pprts_icollapse

    integer(iintegers) :: i, j, icol

    ! vertical thickness in [m]
    real(ireals), allocatable :: dz_t2b(:, :, :) ! dz (t2b := top 2 bottom)

    if (present(nxproc) .neqv. present(nyproc)) then
      print *, 'Wrong call to init_tenstream_rrtm_lw --    &
            & in order to work, we need both arrays for &
            & the domain decomposition, call with nxproc AND nyproc'
      call CHKERR(1_mpiint, 'init_tenstream_rrtm_lw -- missing arguments nxproc or nyproc')
    end if

    allocate (dz_t2b(zm, xm, ym))
    do j = 1_iintegers, ym
      do i = 1_iintegers, xm
        icol = i + (j - 1_iintegers) * xm
        dz_t2b(:, i, j) = reverse(dz(:, icol))
      end do
    end do

    if (present(nxproc) .and. present(nyproc)) then
      call init_pprts(comm, zm, xm, ym, dx, dy, sundir, solver, nxproc=nxproc, nyproc=nyproc, dz3d=dz_t2b, &
                      collapseindex=pprts_icollapse)
    else ! we let petsc decide where to put stuff
      call init_pprts(comm, zm, xm, ym, dx, dy, sundir, solver, dz3d=dz_t2b, &
                      collapseindex=pprts_icollapse)
    end if
  end subroutine

  subroutine repwvl_pprts_destroy(solver, lfinalizepetsc)
    class(t_solver) :: solver
    logical, intent(in) :: lfinalizepetsc

    call destroy_pprts(solver, lfinalizepetsc=lfinalizepetsc)
  end subroutine
end module
