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

module m_repwvl_plexrt
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & iintegers, ireals, mpiint, &
    & zero, one, default_str_len

  use m_helper_functions, only: &
    & CHKERR, &
    & get_arg, &
    & get_petsc_opt, &
    & is_inrange, &
    & toStr

  use m_tenstream_options, only: read_commandline_options

  use m_dyn_atm_to_rrtmg, only: &
    & planck, &
    & print_tenstr_atm, &
    & t_tenstr_atm

  use m_repwvl_base, only: repwvl_init, t_repwvl_data, repwvl_log_events
  use m_repwvl_optprop, only: repwvl_optprop, check_fu_table_consistency
  use m_mie_tables, only: mie_tables_init, t_mie_table, destroy_mie_table
  use m_fu_ice, only: fu_ice_init

  use m_plex_grid, only: TOAFACE
  use m_plex_rt_base, only: t_plex_solver
  use m_plex_rt, only: run_plex_rt_solver, destroy_plexrt_solver, plexrt_get_result

  use m_icon_plex_utils, only: Nz_Ncol_vec_to_celldm1, Nz_Ncol_vec_to_horizface1_dm

  implicit none

  private
  public :: repwvl_plexrt, repwvl_plexrt_destroy

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .true.
#else
  logical, parameter :: ldebug = .true.
#endif

  type(t_repwvl_data), allocatable :: repwvl_data_solar, repwvl_data_thermal
  type(t_mie_table), allocatable :: repwvl_mie_table

  real(ireals), parameter :: CO = 1e-9_ireals
  real(ireals), parameter :: HNO3 = 1e-9_ireals
  real(ireals), parameter :: N2 = 0.78102_ireals

contains

  subroutine repwvl_plexrt(solver, atm, sundir, &
                           albedo_thermal, albedo_solar, &
                           lthermal, lsolar, &
                           edir, edn, eup, abso, &
                           opt_time, solar_albedo_2d, thermal_albedo_2d, &
                           opt_solar_constant)

    class(t_plex_solver), allocatable, intent(inout) :: solver ! solver type -- includes a dmplex and more
    type(t_tenstr_atm), intent(in) :: atm                      ! atmosphere construct, which describes physical properties on layers and levels up till TOA
    real(ireals), intent(in) :: sundir(:)                ! cartesian direction of sun rays
    real(ireals), intent(in) :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions and compute new solutions only after threshold estimate is exceeded.
    ! If solar_albedo_2d is present, we use a 2D surface albedo, dim=(ncol,)
    ! TODO: introduce solar diffuse albedo.
    ! Or hack something together ala Ritter/Geleyn -> see:
    ! icon-lem/src/atm_phy_nwp/mo_nwp_rrtm_interface.f90 l.1371
    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:), thermal_albedo_2d(:), opt_solar_constant

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :), intent(inout) :: edir, edn, eup, abso          ! [nlyr(+1), ncol]

    ! ---------- end of API ----------------

    type(tIS) :: toa_ids

    integer(iintegers) :: Ncol, ke, ke1

    integer(mpiint) :: myid, comm, ierr
    logical :: lskip_thermal, lskip_solar, lflg

    if (present(opt_time)) call CHKERR(1_mpiint, 'not yet implemented.... opt_time')
    if (present(solar_albedo_2d)) call CHKERR(1_mpiint, 'not yet implemented.... solar_albedo_2d')
    if (present(thermal_albedo_2d)) call CHKERR(1_mpiint, 'not yet implemented.... thermal_albedo_2d')
    if (present(opt_solar_constant)) call CHKERR(1_mpiint, 'not yet implemented.... opt_solar_constant')

    if (.not. allocated(solver)) call CHKERR(1_mpiint, 'solver has to be setup beforehand')
    if (.not. allocated(solver%plex)) call CHKERR(1_mpiint, 'Solver has to have a ready to go Plexgrid')

    call PetscObjectGetComm(solver%plex%dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    ke1 = solver%plex%Nlay + 1
    call CHKERR(int(ke1 - size(atm%plev, dim=1), mpiint), 'Vertical Size of atm and plex solver dont match')
    ke = ke1 - 1

    call DMGetStratumIS(solver%plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
    call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)

    if (ldebug .and. myid .eq. 0) then
      call print_tenstr_atm(atm, Ncol)
      print *, 'm_plexrt_rrtmg sundir:', sundir, 'albedo th,sol', albedo_thermal, albedo_solar, 'lth/lsol', lthermal, lsolar
      if (present(opt_time)) print *, 'time', opt_time
    end if

    if (ldebug) then
      if (present(solar_albedo_2d)) then
        if (any(.not. is_inrange(solar_albedo_2d, zero, one))) &
          call CHKERR(1_mpiint, 'Bad solar albedo value min: '//toStr(minval(solar_albedo_2d))// &
                      ' max: '//toStr(maxval(solar_albedo_2d)))
      end if
      if (present(thermal_albedo_2d)) then
        if (any(.not. is_inrange(thermal_albedo_2d, zero, one))) &
          call CHKERR(1_mpiint, 'Bad thermal albedo value min: '//toStr(minval(thermal_albedo_2d))// &
                      ' max: '//toStr(maxval(thermal_albedo_2d)))
      end if
    end if

    if (.not. allocated(repwvl_data_solar)) then
      call repwvl_init(        &
        & comm,                &
        & repwvl_data_solar,   &
        & repwvl_data_thermal, &
        & ierr,                &
        & lverbose=.false.); call CHKERR(ierr)

      call mie_tables_init(comm, repwvl_mie_table, ierr, lverbose=.false.); call CHKERR(ierr)

      call fu_ice_init(comm, ierr, lverbose=.true.); call CHKERR(ierr)
      call check_fu_table_consistency(repwvl_data_solar, repwvl_data_thermal)
    end if

    if (.not. allocated(edn)) allocate (edn(ke1, Ncol))
    if (.not. allocated(eup)) allocate (eup(ke1, Ncol))
    if (.not. allocated(abso)) allocate (abso(ke, Ncol))
    edn = zero
    eup = zero
    abso = zero

    ! make sure that optprop vecs are allocated
    call allocate_optprop_vec(solver%plex%cell1_dm, solver%kabs)
    call allocate_optprop_vec(solver%plex%cell1_dm, solver%ksca)
    call allocate_optprop_vec(solver%plex%cell1_dm, solver%g)
    call allocate_optprop_vec(solver%plex%srfc_boundary_dm, solver%albedo)

    lskip_thermal = .false.
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-skip_thermal", lskip_thermal, lflg, ierr); call CHKERR(ierr)
    if (lthermal .and. .not. lskip_thermal) then
      call allocate_optprop_vec(solver%plex%horizface1_dm, solver%plck)

      call PetscLogStagePush(repwvl_log_events%stage_repwvl_thermal, ierr); call CHKERR(ierr)
      call compute_thermal(comm, solver, &
                           & repwvl_data_thermal, &
                           & repwvl_mie_table, atm, &
                           & Ncol, ke1, &
                           & sundir, albedo_thermal, &
                           & edn, eup, abso, &
                           & ierr, &
                           & opt_time=opt_time, &
                           & thermal_albedo_2d=thermal_albedo_2d); call CHKERR(ierr)
      call PetscLogStagePop(ierr); call CHKERR(ierr) ! pop solver%logs%stage_repwvl_thermal
    end if

    if (lsolar .and. .not. allocated(edir)) allocate (edir(ke1, Ncol))
    if (allocated(edir)) edir = zero ! if user gave edir, make sure its a current value

    if (lsolar) then

      lskip_solar = .false.
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-skip_solar", lskip_solar, lflg, ierr); call CHKERR(ierr)
      if (.not. lskip_solar) then
        call PetscLogStagePush(repwvl_log_events%stage_repwvl_solar, ierr); call CHKERR(ierr)
        call compute_solar(comm, solver, &
                           & repwvl_data_solar, &
                           & repwvl_mie_table, atm, &
                           & Ncol, ke1, &
                           & sundir, albedo_solar, &
                           & edir, edn, eup, abso, &
                           & ierr, &
                           & opt_time=opt_time, &
                           & solar_albedo_2d=solar_albedo_2d, &
                           & opt_solar_constant=opt_solar_constant); call CHKERR(ierr)
        call PetscLogStagePop(ierr); call CHKERR(ierr) ! pop solver%logs%stage_repwvl_solar
      end if
    end if
  end subroutine

  subroutine compute_thermal(comm, solver, &
      & repwvl_data, mie_table, atm, &
      & Ncol, ke1, sundir, albedo, &
      & edn, eup, abso, &
      & ierr, opt_time, thermal_albedo_2d)

    integer(mpiint), intent(in) :: comm
    class(t_plex_solver), allocatable, intent(inout) :: solver
    type(t_repwvl_data), intent(in) :: repwvl_data
    type(t_mie_table), intent(in) :: mie_table
    type(t_tenstr_atm), intent(in), target :: atm
    integer(iintegers), intent(in) :: Ncol, ke1

    real(ireals), intent(in) :: sundir(:)
    real(ireals), intent(in) :: albedo

    real(ireals), intent(inout), dimension(:, :) :: edn, eup, abso
    integer(mpiint), intent(out) :: ierr

    real(ireals), optional, intent(in) :: opt_time, thermal_albedo_2d(:)

    real(ireals), dimension(:, :), allocatable :: spec_edn, spec_eup, spec_abso
    real(ireals), dimension(:, :), allocatable :: Blev, kabs, ksca, kg      ! [nlyr, ncol]
    integer(iintegers) :: spectral_bands(2), argcnt
    integer(iintegers) :: i, k, ke, iwvl
    logical :: lflg
    integer(mpiint) :: myid

    real(ireals), pointer :: xalbedo(:)

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    ke = ke1 - 1

    allocate (spec_edn(ke1, Ncol))
    allocate (spec_eup(ke1, Ncol))
    allocate (spec_abso(ke, Ncol))

    allocate (Blev(ke1, Ncol))
    allocate (kabs(ke, Ncol))
    allocate (ksca(ke, Ncol))
    allocate (kg(ke, Ncol))

    spectral_bands = [integer(iintegers) :: 1, size(repwvl_data%wvls)]
    argcnt = size(spectral_bands)
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-repwvl_bands", spectral_bands, argcnt, lflg, ierr); call CHKERR(ierr)
    if (lflg) call CHKERR(int(argcnt - 2_iintegers, mpiint), "must provide 2 values for repwvl_bands, comma separated, no spaces")
    spectral_bands = min(max(spectral_bands, 1), size(repwvl_data%wvls))

    do iwvl = spectral_bands(1), spectral_bands(2)
      if (myid .eq. 0 .and. (lflg .or. ldebug)) then
        print *, 'Computing wavelengths '//toStr(iwvl)//' / '//toStr(size(repwvl_data%wvls))//&
          & ' -- '//toStr(100._ireals * real(iwvl, ireals) / real(size(repwvl_data%wvls), ireals))//' %'// &
          & ' ('//toStr(repwvl_data%wvls(iwvl))//' nm,  wgt='//toStr(repwvl_data%wgts(iwvl))//')'
      end if

      call PetscLogEventBegin(repwvl_log_events%repwvl_optprop, ierr); call CHKERR(ierr)
      do i = 1, Ncol
        do k = 1, ke
          call repwvl_optprop(&
            & repwvl_data, atm, mie_table, &
            & .false., k, i, iwvl, &
            & kabs(size(kabs, dim=1) + 1 - k, i), &
            & ksca(size(ksca, dim=1) + 1 - k, i), &
            & kg(size(kg, dim=1) + 1 - k, i), &
            & ierr); call CHKERR(ierr)
        end do
        do k = 1, ke1
          Blev(ke1 + 1 - k, i) = repwvl_data%wgts(iwvl) &
                                & * planck(repwvl_data%wvls(iwvl) * 1e-9_ireals, atm%tlev(k, i)) &
                                & * 1e-9_ireals
        end do
      end do
      if (allocated(atm%tskin)) then
        do i = 1, Ncol
          Blev(ke1, i) = repwvl_data%wgts(iwvl) &
                        & * planck(repwvl_data%wvls(iwvl) * 1e-9_ireals, atm%tskin(i)) &
                        & * 1e-9_ireals
        end do
      end if
      call PetscLogEventEnd(repwvl_log_events%repwvl_optprop, ierr); call CHKERR(ierr)

      call VecGetArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
      if (present(thermal_albedo_2d)) then
        xalbedo = thermal_albedo_2d
      else
        xalbedo = albedo
      end if
      call VecRestoreArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)

      call Nz_Ncol_vec_to_horizface1_dm(solver%plex, Blev, solver%plck)
      call Nz_Ncol_vec_to_celldm1(solver%plex, kabs, solver%kabs)
      call Nz_Ncol_vec_to_celldm1(solver%plex, ksca, solver%ksca)
      call Nz_Ncol_vec_to_celldm1(solver%plex, kg, solver%g)

      call run_plex_rt_solver(&
        & solver, &
        & lthermal=.true., lsolar=.false., &
        & sundir=sundir, &
        & opt_solution_uid=-iwvl, opt_solution_time=opt_time)

      call plexrt_get_result(&
        & solver, &
        & spec_edn, spec_eup, spec_abso, &
        & opt_solution_uid=-iwvl)

      edn = edn + spec_edn
      eup = eup + spec_eup
      abso = abso + spec_abso
    end do !iwvl
  end subroutine

  subroutine compute_solar(comm, solver, &
      & repwvl_data, mie_table, atm, &
      & Ncol, ke1, sundir, albedo, &
      & edir, edn, eup, abso, &
      & ierr, &
      & opt_time, solar_albedo_2d, &
      & opt_solar_constant)

    integer(mpiint), intent(in) :: comm
    class(t_plex_solver), allocatable, intent(inout) :: solver
    type(t_repwvl_data), intent(in) :: repwvl_data
    type(t_mie_table), intent(in) :: mie_table
    type(t_tenstr_atm), intent(in), target :: atm
    integer(iintegers), intent(in) :: Ncol, ke1

    real(ireals), intent(in) :: albedo
    real(ireals), intent(in) :: sundir(:)

    real(ireals), intent(inout), dimension(:, :) :: edir, edn, eup, abso
    integer(mpiint), intent(out) :: ierr

    real(ireals), intent(in), optional :: opt_time, solar_albedo_2d(:)
    real(ireals), intent(in), optional :: opt_solar_constant

    real(ireals), dimension(:, :), allocatable :: spec_edir, spec_edn, spec_eup, spec_abso
    real(ireals), dimension(:, :), allocatable :: kabs, ksca, kg      ! [nlyr, ncol]
    integer(iintegers) :: spectral_bands(2), argcnt
    integer(iintegers) :: i, k, ke, iwvl
    logical :: lflg
    integer(mpiint) :: myid
    real(ireals) :: edirTOA, rescaled_sundir(3)

    real(ireals), pointer :: xalbedo(:)

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    ke = ke1 - 1

    allocate (spec_edir(ke1, Ncol))
    allocate (spec_edn(ke1, Ncol))
    allocate (spec_eup(ke1, Ncol))
    allocate (spec_abso(ke, Ncol))

    allocate (kabs(ke, Ncol))
    allocate (ksca(ke, Ncol))
    allocate (kg(ke, Ncol))

    spectral_bands = [integer(iintegers) :: 1, size(repwvl_data%wvls)]
    argcnt = size(spectral_bands)
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-repwvl_bands", spectral_bands, argcnt, lflg, ierr); call CHKERR(ierr)
    if (lflg) call CHKERR(int(argcnt - 2_iintegers, mpiint), "must provide 2 values for repwvl_bands, comma separated, no spaces")
    spectral_bands = min(max(spectral_bands, 1), size(repwvl_data%wvls))

    do iwvl = spectral_bands(1), spectral_bands(2)
      if (myid .eq. 0 .and. (lflg .or. ldebug)) then
        print *, 'Computing wavelengths '//toStr(iwvl)//' / '//toStr(size(repwvl_data%wvls))//&
          & ' -- '//toStr(100._ireals * real(iwvl, ireals) / real(size(repwvl_data%wvls), ireals))//' %'// &
          & ' ('//toStr(repwvl_data%wvls(iwvl))//' nm,  wgt='//toStr(repwvl_data%wgts(iwvl))//')'
      end if

      edirTOA = repwvl_data%wgts(iwvl)
      if (present(opt_solar_constant)) then
        edirTOA = repwvl_data%wgts(iwvl) / sum(repwvl_data%wgts(:)) * opt_solar_constant
      end if
      rescaled_sundir = sundir / norm2(sundir) * edirTOA

      call PetscLogEventBegin(repwvl_log_events%repwvl_optprop, ierr); call CHKERR(ierr)
      do i = 1, Ncol
        do k = 1, ke
          call repwvl_optprop(&
            & repwvl_data, atm, mie_table, &
            & .true., k, i, iwvl, &
            & kabs(size(kabs, dim=1) + 1 - k, i), &
            & ksca(size(ksca, dim=1) + 1 - k, i), &
            & kg(size(kg, dim=1) + 1 - k, i), &
            & ierr); call CHKERR(ierr)
        end do
      end do
      call PetscLogEventEnd(repwvl_log_events%repwvl_optprop, ierr); call CHKERR(ierr)

      call VecGetArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
      if (present(solar_albedo_2d)) then
        xalbedo = solar_albedo_2d
      else
        xalbedo = albedo
      end if
      call VecRestoreArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)

      call Nz_Ncol_vec_to_celldm1(solver%plex, kabs, solver%kabs)
      call Nz_Ncol_vec_to_celldm1(solver%plex, ksca, solver%ksca)
      call Nz_Ncol_vec_to_celldm1(solver%plex, kg, solver%g)

      call run_plex_rt_solver(&
        & solver, &
        & lthermal=.false., lsolar=.true., &
        & sundir=rescaled_sundir, &
        & opt_solution_uid=iwvl, opt_solution_time=opt_time)

      call plexrt_get_result(&
        & solver, &
        & spec_edn, spec_eup, spec_abso, spec_edir, &
        & opt_solution_uid=iwvl)

      edir = edir + spec_edir
      edn = edn + spec_edn
      eup = eup + spec_eup
      abso = abso + spec_abso
    end do !iwvl
  end subroutine

  subroutine repwvl_plexrt_destroy(solver, lfinalizepetsc, ierr)
    class(t_plex_solver), allocatable, intent(inout) :: solver
    logical, intent(in) :: lfinalizepetsc
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    call destroy_plexrt_solver(solver, lfinalizepetsc=lfinalizepetsc)
  end subroutine

  subroutine allocate_optprop_vec(dm, vec, val)
    type(tDM), intent(in) :: dm
    type(tVec), allocatable, intent(inout) :: vec
    real(ireals), intent(in), optional :: val
    integer(mpiint) :: ierr

    if (.not. allocated(vec)) then
      allocate (vec)
      call DMCreateGlobalVector(dm, vec, ierr); call CHKERR(ierr)
    end if
    if (present(val)) then
      call VecSet(vec, val, ierr); call CHKERR(ierr)
    end if
  end subroutine
end module
