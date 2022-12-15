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

!> \page Routines to call tenstream with optical properties from a ecrad

module m_ecckd_pprts
#include "petsc/finclude/petsc.h"
  use petsc

  use m_helper_functions, only: &
    & CHKERR, &
    & get_arg, &
    & get_petsc_opt, &
    & mpi_logical_all_same, &
    & reverse, &
    & toStr

  use m_data_parameters, only: &
    & iintegers, ireals, mpiint, &
    & zero, one, default_str_len, &
    & i1, &
    & share_dir

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
    & print_tenstr_atm, &
    & t_tenstr_atm

  use m_buildings, only: &
    & clone_buildings, &
    & destroy_buildings, &
    & t_pprts_buildings

  use m_ecckd_base, only: ecckd_init, t_ecckd_data, ecckd_log_events
  use m_ecckd_optprop, only: ecckd_optprop, ecckd_planck
  use m_mie_tables, only: mie_tables_init, t_mie_table, destroy_mie_table
  use m_fu_ice, only: fu_ice_init, fu_muskatel_ice_data

  use m_pprts_rrtmg, only: smooth_surface_fluxes, slope_correction_fluxes

  implicit none

  private
  public :: ecckd_pprts, ecckd_pprts_destroy

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .false.
#else
  logical, parameter :: ldebug = .true.
#endif

  type(t_ecckd_data), allocatable :: ecckd_data_solar, ecckd_data_thermal
  type(t_mie_table), allocatable :: ecckd_general_mie_table

contains

  subroutine ecckd_pprts(comm, solver, atm, ie, je,   &
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

    ! optional optical properties with dims(nlyr(srfc to TOA), local_nx, local_ny)
    ! used e.g. for aerosol or vegetation, you can provide only tau or (tau and w0) defaults to w0=0 and g=0
    ! note that first dimension can also be smaller than nlyr, we will fill it up from the ground,
    ! i.e. if you provide only two layers, the two lowermost layers near the surface will be filled with addition optprops
    real(ireals), intent(in), optional, dimension(:, :, :) :: opt_tau_solar, opt_w0_solar, opt_g_solar
    real(ireals), intent(in), optional, dimension(:, :, :) :: opt_tau_thermal

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

      call mie_tables_init(comm, ecckd_general_mie_table, ierr, lverbose=.false., &
        & path_table=share_dir//"mie_droplet_scattering.nc"); call CHKERR(ierr)

      call fu_ice_init(comm, ierr, lverbose=.false.); call CHKERR(ierr)

      call ecckd_init(        &
        & comm,               &
        & atm,                &
        & ecckd_general_mie_table, &
        & fu_muskatel_ice_data, &
        & ecckd_data_solar,   &
        & ecckd_data_thermal, &
        & ierr,                &
        & lverbose=.false.); call CHKERR(ierr)

      call init_pprts_ecckd(comm, solver, &
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
                       "-ecckd_pprts_atm_view", lprint_atm, lflg, ierr); call CHKERR(ierr)
    if (lprint_atm .and. myid .eq. 0) then
      call print_tenstr_atm(atm)
    end if

    if (get_arg(.false., lonly_initialize)) return

    lskip_thermal = .false.
    call get_petsc_opt(PETSC_NULL_CHARACTER, &
                       "-skip_thermal", lskip_thermal, lflg, ierr); call CHKERR(ierr)
    if (lthermal .and. .not. lskip_thermal) then
      call PetscLogStagePush(ecckd_log_events%stage_ecckd_thermal, ierr); call CHKERR(ierr)
      call compute_thermal(                    &
        & comm,                                &
        & ecckd_data_thermal,                  &
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
      call PetscLogStagePop(ierr); call CHKERR(ierr) ! pop log_events%stage_ecckd_thermal
    end if

    if (lsolar .and. .not. allocated(edir)) allocate (edir(solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
    if (allocated(edir)) edir = zero

    if (lsolar) then
      lskip_solar = .false.
      call get_petsc_opt(PETSC_NULL_CHARACTER, &
                         "-skip_solar", lskip_solar, lflg, ierr); call CHKERR(ierr)
      if (.not. lskip_solar) then
        call PetscLogStagePush(ecckd_log_events%stage_ecckd_solar, ierr); call CHKERR(ierr)
        call compute_solar(                          &
          & comm,                                    &
          & ecckd_data_solar,                        &
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
        call PetscLogStagePop(ierr); call CHKERR(ierr) ! pop log_events%stage_ecckd_solar
      end if
    end if

    call smooth_surface_fluxes(solver, edn, eup)
    if (lsolar) call slope_correction_fluxes(solver, edir)
  end subroutine

  subroutine compute_thermal( &
      & comm,                 &
      & ecckd_data_thermal,  &
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

    integer(mpiint), intent(in) :: comm
    type(t_ecckd_data), intent(in) :: ecckd_data_thermal
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in), target :: atm
    real(ireals), intent(in) :: albedo
    real(ireals), intent(inout), dimension(:, :, :) :: edn, eup, abso
    integer(mpiint), intent(out) :: ierr

    real(ireals), optional, intent(in) :: opt_time, thermal_albedo_2d(:, :)
    type(t_pprts_buildings), intent(inout), optional :: opt_buildings
    real(ireals), intent(in), optional, dimension(:, :, :) :: opt_tau

    integer(iintegers) :: ig, i, j, k, icol, ke, ke1
    integer(mpiint) :: myid

    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, kg ! [nlyr, local_nx, local_ny]
    real(ireals), allocatable :: Blev(:, :, :)
    real(ireals), allocatable, dimension(:, :) :: Bsrfc             ! [local_nx, local_ny]

    real(ireals), dimension(:, :, :), allocatable :: spec_edn, spec_eup, spec_abso

    type(t_pprts_buildings), allocatable :: spec_buildings

    logical :: lflg
    integer(iintegers) :: argcnt, spectral_gpts(2)

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr)

    ke1 = ubound(atm%plev, 1)
    ke = ubound(atm%tlay, 1)

    allocate (kabs(ke, solver%C_one%xm, solver%C_one%ym))
    allocate (ksca(ke, solver%C_one%xm, solver%C_one%ym))
    allocate (kg(ke, solver%C_one%xm, solver%C_one%ym))
    allocate (Blev(ke1, solver%C_one1%xm, solver%C_one1%ym))
    allocate (Bsrfc(solver%C_one1%xm, solver%C_one1%ym))

    if (present(opt_buildings)) then

      call clone_buildings(opt_buildings, spec_buildings, l_copy_data=.false., ierr=ierr); call CHKERR(ierr)
      if (.not. allocated(spec_buildings%albedo)) allocate (spec_buildings%albedo(size(opt_buildings%albedo)))
      if (.not. allocated(opt_buildings%incoming)) allocate (opt_buildings%incoming(size(opt_buildings%iface)))
      if (.not. allocated(opt_buildings%outgoing)) allocate (opt_buildings%outgoing(size(opt_buildings%iface)))
      spec_buildings%albedo = opt_buildings%albedo
      opt_buildings%incoming = zero
      opt_buildings%outgoing = zero

      if (.not. (allocated(opt_buildings%temp))) &
        & call CHKERR(1_mpiint, 'Thermal computation but opt_buildings%temp is not allocated')
      if ((allocated(opt_buildings%planck))) &
        & call CHKERR(1_mpiint, 'Thermal computation but opt_buildings%planck is allocated... '//&
                              & 'you should only provide temperatures. I`ll take care of planck emissison')

      !allocate (spec_buildings%temp(size(opt_buildings%temp)))
      allocate (spec_buildings%planck(size(opt_buildings%temp)))
    end if

    spectral_gpts = [integer(iintegers) :: 1, ecckd_data_thermal%n_g_pnt]
    argcnt = size(spectral_gpts)
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-ecckd_bands", spectral_gpts, argcnt, lflg, ierr); call CHKERR(ierr)
    if (lflg) call CHKERR(int(argcnt - 2_iintegers, mpiint), "must provide 2 values for ecckd_bands, comma separated, no spaces")
    spectral_gpts = min(max(spectral_gpts, 1), ecckd_data_thermal%n_g_pnt)

    do ig = spectral_gpts(1), spectral_gpts(2)
      call spectral_integration_progress_report(comm, ecckd_data_thermal, ig, ierr); call CHKERR(ierr)
      call PetscLogEventBegin(ecckd_log_events%ecckd_optprop, ierr); call CHKERR(ierr)

      do j = 1, solver%C_one%ym
        do i = 1, solver%C_one%xm
          icol = i + (j - 1_iintegers) * solver%C_one%xm

          do k = 1, ke
            call ecckd_optprop(&
              & ecckd_data_thermal, atm, &
              & .false., k, icol, ig, &
              & kabs(size(kabs, dim=1) + 1 - k, i, j), &
              & ksca(size(ksca, dim=1) + 1 - k, i, j), &
              & kg(size(kg, dim=1) + 1 - k, i, j), &
              & ierr); call CHKERR(ierr)
          end do

          do k = 1, ke1
            call ecckd_planck( &
              & ecckd_data_thermal, &
              & ig, &
              & atm%tlev(k, icol), &
              & Blev(ke1 + 1 - k, i, j), &
              & ierr); call CHKERR(ierr)
          end do
          if (allocated(atm%tskin)) then
            call ecckd_planck( &
              & ecckd_data_thermal, &
              & ig, &
              & atm%tskin(icol), &
              & Bsrfc(i, j), &
              & ierr); call CHKERR(ierr)
          else
            Bsrfc(i, j) = Blev(ke1, i, j)
          end if
        end do
      end do

      if (present(opt_buildings)) then
        ! Use spec_buildings%temperature array to hold raw planck values ...
        do k = 1, size(opt_buildings%temp)
          call ecckd_planck( &
            & ecckd_data_thermal, &
            & ig, &
            & opt_buildings%temp(k), &
            & spec_buildings%planck(k), &
            & ierr); call CHKERR(ierr)
        end do
      end if

      call add_optional_optprop(solver%atm%dz, kabs=kabs, opt_tau=opt_tau)

      call PetscLogEventEnd(ecckd_log_events%ecckd_optprop, ierr); call CHKERR(ierr)

      call set_optical_properties( &
        & solver,                  &
        & albedo,                  &
        & kabs,                    &
        & ksca,                    &
        & kg,                      &
        & planck=Blev,             &
        & planck_srfc=Bsrfc,       &
        & albedo_2d=thermal_albedo_2d)

      call solve_pprts(               &
        & solver,                     &
        & lthermal=.true.,            &
        & lsolar=.false.,             &
        & edirTOA=-1._ireals,         &
        & opt_solution_uid=-ig,     &
        & opt_solution_time=opt_time, &
        & opt_buildings=spec_buildings)

      if (present(opt_buildings)) then
        call pprts_get_result(           &
          & solver,                      &
          & spec_edn,                    &
          & spec_eup,                    &
          & spec_abso,                   &
          & opt_solution_uid=-ig,      &
          & opt_buildings=spec_buildings)
        opt_buildings%incoming = opt_buildings%incoming + spec_buildings%incoming
        opt_buildings%outgoing = opt_buildings%outgoing + spec_buildings%outgoing
      else
        call pprts_get_result(    &
          & solver,               &
          & spec_edn,             &
          & spec_eup,             &
          & spec_abso,            &
          & opt_solution_uid=-ig)
      end if

      edn = edn + spec_edn
      eup = eup + spec_eup
      abso = abso + spec_abso

    end do !ig

    if (present(opt_buildings)) then
      call destroy_buildings(spec_buildings, ierr)
    end if
  end subroutine

  subroutine compute_solar( &
      & comm,               &
      & ecckd_data_solar,  &
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

    integer(mpiint), intent(in) :: comm
    type(t_ecckd_data), intent(in) :: ecckd_data_solar
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in), target :: atm

    real(ireals), intent(in) :: albedo
    real(ireals), intent(in) :: sundir(3)

    real(ireals), intent(inout), dimension(:, :, :) :: edir, edn, eup, abso
    integer(mpiint), intent(out) :: ierr

    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:, :)
    real(ireals), intent(in), optional :: opt_solar_constant
    type(t_pprts_buildings), intent(inout), optional :: opt_buildings
    real(ireals), intent(in), optional, dimension(:, :, :) :: opt_tau, opt_w0, opt_g

    real(ireals) :: edirTOA

    integer(iintegers) :: i, j, icol, k, ig, ke
    integer(mpiint) :: myid

    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, kg      ! [nlyr, local_nx, local_ny]

    real(ireals), dimension(:, :, :), allocatable :: spec_edir, spec_edn, spec_eup, spec_abso

    type(t_pprts_buildings), allocatable :: spec_buildings

    logical :: lflg
    integer(iintegers) :: argcnt, spectral_gpts(2)

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr)

    ke = ubound(atm%tlay, 1)

    allocate (kabs(ke, solver%C_one%xm, solver%C_one%ym))
    allocate (ksca(ke, solver%C_one%xm, solver%C_one%ym))
    allocate (kg(ke, solver%C_one%xm, solver%C_one%ym))

    if (present(opt_buildings)) then
      call clone_buildings(opt_buildings, spec_buildings, l_copy_data=.false., ierr=ierr); call CHKERR(ierr)
      if (.not. allocated(spec_buildings%albedo)) allocate (spec_buildings%albedo(size(opt_buildings%albedo)))
      if (.not. allocated(opt_buildings%edir)) allocate (opt_buildings%edir(size(opt_buildings%iface)))
      if (.not. allocated(opt_buildings%incoming)) allocate (opt_buildings%incoming(size(opt_buildings%iface)))
      if (.not. allocated(opt_buildings%outgoing)) allocate (opt_buildings%outgoing(size(opt_buildings%iface)))
      spec_buildings%albedo = opt_buildings%albedo
      opt_buildings%edir = zero
      opt_buildings%incoming = zero
      opt_buildings%outgoing = zero
    end if

    call set_angles(solver, sundir)

    spectral_gpts = [integer(iintegers) :: 1, ecckd_data_solar%n_g_pnt]
    argcnt = size(spectral_gpts)
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-ecckd_bands", spectral_gpts, argcnt, lflg, ierr); call CHKERR(ierr)
    if (lflg) call CHKERR(int(argcnt - 2_iintegers, mpiint), "must provide 2 values for ecckd_bands, comma separated, no spaces")
    spectral_gpts = min(max(spectral_gpts, 1), ecckd_data_solar%n_g_pnt)

    do ig = spectral_gpts(1), spectral_gpts(2)
      call spectral_integration_progress_report(comm, ecckd_data_solar, ig, ierr); call CHKERR(ierr)
      call PetscLogEventBegin(ecckd_log_events%ecckd_optprop, ierr); call CHKERR(ierr)

      do j = 1, solver%C_one%ym
        do i = 1, solver%C_one%xm
          icol = i + (j - 1_iintegers) * solver%C_one%xm

          do k = 1, ke
            call ecckd_optprop(&
              & ecckd_data_solar, atm, &
              & .true., k, icol, ig, &
              & kabs(size(kabs, dim=1) + 1 - k, i, j), &
              & ksca(size(ksca, dim=1) + 1 - k, i, j), &
              & kg(size(kg, dim=1) + 1 - k, i, j), &
              & ierr); call CHKERR(ierr)
          end do
        end do
      end do

      call add_optional_optprop(solver%atm%dz, kabs, ksca, kg, opt_tau, opt_w0, opt_g)
      call PetscLogEventEnd(ecckd_log_events%ecckd_optprop, ierr); call CHKERR(ierr)

      edirTOA = ecckd_data_solar%solar_irradiance(ig)
      if (present(opt_solar_constant)) then
        edirTOA = ecckd_data_solar%solar_irradiance(ig) / sum(ecckd_data_solar%solar_irradiance(:)) * opt_solar_constant
      end if

      call set_optical_properties( &
        & solver,                  &
        & albedo,                  &
        & kabs,                    &
        & ksca,                    &
        & kg,                      &
        & albedo_2d=solar_albedo_2d)

      call solve_pprts(          &
        & solver,                &
        & lthermal=.false.,      &
        & lsolar=.true.,         &
        & edirTOA=1._ireals,     &
        & opt_solution_uid=ig, &
        & opt_solution_time=opt_time,&
        & opt_buildings=spec_buildings)

      if (present(opt_buildings)) then
        call pprts_get_result(          &
          & solver,                     &
          & spec_edn,                   &
          & spec_eup,                   &
          & spec_abso,                  &
          & spec_edir,                  &
          & opt_solution_uid=ig,      &
          & opt_buildings=spec_buildings)

        opt_buildings%edir = opt_buildings%edir + spec_buildings%edir
        opt_buildings%incoming = opt_buildings%incoming + spec_buildings%incoming
        opt_buildings%outgoing = opt_buildings%outgoing + spec_buildings%outgoing

      else

        call pprts_get_result(          &
          & solver,                     &
          & spec_edn,                   &
          & spec_eup,                   &
          & spec_abso,                  &
          & spec_edir,                  &
          & opt_solution_uid=ig)

      end if

      edir = edir + spec_edir * edirTOA
      edn = edn + spec_edn * edirTOA
      eup = eup + spec_eup * edirTOA
      abso = abso + spec_abso * edirTOA
    end do !ig

    if (present(opt_buildings)) then
      call destroy_buildings(spec_buildings, ierr)
    end if
  end subroutine compute_solar

  subroutine init_pprts_ecckd(comm, solver, dx, dy, dz, &
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

  subroutine ecckd_pprts_destroy(solver, lfinalizepetsc, ierr)
    class(t_solver) :: solver
    logical, intent(in) :: lfinalizepetsc
    integer(mpiint), intent(out) :: ierr
    ierr = 0

    if (allocated(ecckd_data_solar)) deallocate (ecckd_data_solar)
    if (allocated(ecckd_data_thermal)) deallocate (ecckd_data_thermal)
    call destroy_mie_table(ecckd_general_mie_table, ierr); call CHKERR(ierr)
    call destroy_pprts(solver, lfinalizepetsc=lfinalizepetsc)
  end subroutine

  subroutine spectral_integration_progress_report(comm, ecckd_data, ig, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_ecckd_data), intent(in) :: ecckd_data
    integer(iintegers), intent(in) :: ig
    integer(mpiint), intent(out) :: ierr
    integer(iintegers) :: iband
    real(ireals) :: wvl_lo, wvl_hi
    integer(mpiint) :: myid

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) then
      iband = ecckd_data%band_number(ig) + i1

      wvl_lo = 1e7_ireals / ecckd_data%wavenumber2_band(iband)
      wvl_hi = 1e7_ireals / ecckd_data%wavenumber1_band(iband)

      print *, 'Computing wavelengths '//toStr(ig)//' / '//toStr(ecckd_data%n_g_pnt)//&
        & ' -- '//toStr(100._ireals * real(ig, ireals) / real(ecckd_data%n_g_pnt, ireals))//' %'// &
        & ' ('//toStr(wvl_lo)//' nm  - '//toStr(wvl_hi)//' nm)'
    end if

  end subroutine

  subroutine add_optional_optprop(dz, kabs, ksca, g, opt_tau, opt_w0, opt_g)
    real(ireals), intent(in), dimension(:, :, :) :: dz
    real(ireals), intent(inout), optional, dimension(:, :, :) :: kabs, ksca, g
    real(ireals), intent(in), optional, dimension(:, :, :) :: opt_tau, opt_w0, opt_g

    real(ireals) :: tausca_old, tausca_new, tau, w0
    integer(iintegers) :: k_r, k, i, j

    call check_shape('tau', kabs, opt_tau)
    call check_shape('w0', ksca, opt_w0)
    call check_shape('g', g, opt_g)

    if (present(opt_g) .and. .not. present(opt_w0)) call CHKERR(1_mpiint, 'if opt_g is provided, need also opt_w0')
    if ((present(opt_w0) .or. present(opt_g)) .and. .not. present(opt_tau)) &
      & call CHKERR(1_mpiint, 'if opt_g or opt_w0 is provided, need also opt_tau')

    if (present(opt_g)) then
      if (.not. all([present(kabs), present(ksca), present(g)])) &
        & call CHKERR(1_mpiint, 'if opt_g is provided, need to have kabs, ksca, and g in call to add_optional_optprop')
      do j = 1, size(opt_tau, 3)
        do i = 1, size(opt_tau, 2)
          do k = 1, size(opt_tau, 1)
            k_r = size(kabs, 1) - k + 1
            tausca_old = ksca(k_r, i, j) * dz(k_r, i, j)
            tausca_new = opt_tau(k, i, j) * opt_w0(k, i, j)
            g(k_r, i, j) = (g(k_r, i, j) * tausca_old + opt_g(k, i, j) * tausca_new) / (tausca_old + tausca_new)
            tau = (kabs(k_r, i, j) + ksca(k_r, i, j)) * dz(k_r, i, j) + opt_tau(k, i, j)
            w0 = (tausca_old + tausca_new) / tau
            kabs(k_r, i, j) = tau * (1.-w0) / dz(k_r, i, j)
            ksca(k_r, i, j) = tau * w0 / dz(k_r, i, j)
          end do
        end do
      end do
    elseif (present(opt_w0)) then
      if (.not. all([present(kabs), present(ksca)])) &
        & call CHKERR(1_mpiint, 'if opt_w0 is provided, need to have kabs, ksca in call to add_optional_optprop')
      do j = 1, size(opt_tau, 3)
        do i = 1, size(opt_tau, 2)
          do k = 1, size(opt_tau, 1)
            k_r = size(kabs, 1) - k + 1
            tausca_old = ksca(k_r, i, j) * dz(k_r, i, j)
            tausca_new = opt_tau(k, i, j) * opt_w0(k, i, j)
            g(k_r, i, j) = g(k_r, i, j) * tausca_old / (tausca_old + tausca_new)
            tau = (kabs(k_r, i, j) + ksca(k_r, i, j)) * dz(k_r, i, j) + opt_tau(k, i, j)
            w0 = (tausca_old + tausca_new) / tau
            kabs(k_r, i, j) = tau * (1.-w0) / dz(k_r, i, j)
            ksca(k_r, i, j) = tau * w0 / dz(k_r, i, j)
          end do
        end do
      end do
    elseif (present(opt_tau)) then
      if (any([present(ksca), present(g)])) &
        & call CHKERR(1_mpiint, 'if only opt_tau is provided, only kabs can be in call to add_optional_optprop')
      do j = 1, size(opt_tau, 3)
        do i = 1, size(opt_tau, 2)
          do k = 1, size(opt_tau, 1)
            k_r = size(kabs, 1) - k + 1
            kabs(k_r, i, j) = kabs(k_r, i, j) + opt_tau(k, i, j) / dz(k_r, i, j)
          end do
        end do
      end do
    end if

  contains
    subroutine check_shape(varname, var, opt_var)
      character(len=*), intent(in) :: varname
      real(ireals), intent(in), optional, dimension(:, :, :) :: var
      real(ireals), intent(in), optional, dimension(:, :, :) :: opt_var
      integer :: i

      if (.not. present(opt_var)) return

      if (present(var) .neqv. present(opt_var)) &
        & call CHKERR(-1_mpiint, trim(varname)//' / opt_'//trim(varname)//&
          & ' :: cannot have one argument, need both or none '// &
          & '('//toStr(present(var))//'/'//toStr(present(opt_var))//')')

      if (size(opt_var, 1) .gt. size(var, 1)) then
        call CHKERR(1_mpiint, trim(varname)//" :: first dimension (nlyr) of opt_var "// &
          & " cannot be larger than the target array. "//new_line('')// &
          & " shape(var)    "//toStr(shape(var))//new_line('')// &
          & " shape(opt_var)"//toStr(shape(opt_var)))
      end if
      do i = 2, 3
        if (size(opt_var, i) .gt. size(var, i)) then
          call CHKERR(1_mpiint, trim(varname)//" :: dimension "//toStr(i)//" of opt_var "// &
            & " has to be same as the target array. "//new_line('')// &
            & " shape(var)    "//toStr(shape(var))//new_line('')// &
            & " shape(opt_var)"//toStr(shape(opt_var)))
        end if
      end do
    end subroutine
  end subroutine
end module
