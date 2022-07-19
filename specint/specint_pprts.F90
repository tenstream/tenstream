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

module m_specint_pprts
#include "petsc/finclude/petsc.h"
  use petsc

  use m_helper_functions, only: &
    & CHKERR, &
    & get_arg, &
    & get_petsc_opt, &
    & toStr

  use m_data_parameters, only: &
    & iintegers, ireals, mpiint, &
    & default_str_len

  use m_pprts_base, only: t_solver

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm

  use m_buildings, only: t_pprts_buildings

  use m_pprts_rrtmg, only: pprts_rrtmg, destroy_pprts_rrtmg
  use m_repwvl_pprts, only: repwvl_pprts, repwvl_pprts_destroy

  use m_xdmf_export, only: &
    & xdmf_pprts_buildings, &
    & xdmf_pprts_srfc_flux

  use m_petsc_helpers, only: f90vectopetsc

  implicit none

  private
  public :: specint_pprts, specint_pprts_destroy

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .false.
#else
  logical, parameter :: ldebug = .true.
#endif

contains

  subroutine specint_pprts(specint, comm, solver, atm,&
      & ie, je,                                       &
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

    character(len=*), intent(in) :: specint                  ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm                      ! MPI Communicator
    class(t_solver), intent(inout) :: solver                 ! solver type (e.g. t_solver_8_10)
    type(t_tenstr_atm), intent(in) :: atm                    ! contains info on atmospheric constituents
    integer(iintegers), intent(in) :: ie, je                 ! local domain size in x and y direction
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

    integer(mpiint) :: ierr

    character(len=default_str_len) :: spec
    call select_specint(specint, spec, ierr)
    select case (trim(spec))
    case ("rrtmg")
      call pprts_rrtmg(comm, solver, atm, ie, je,       &
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
    case ("repwvl")
      call repwvl_pprts(comm, solver, atm, ie, je,      &
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
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'invalid specint mode <'//trim(spec)//'>')
    end select

    call dump_results()
  contains
    subroutine dump_variable(var, dm, dumpstring, varname)
      real(ireals), intent(in) :: var(:, :, :)
      type(tDM), intent(in) :: dm
      character(len=*), intent(in) :: dumpstring, varname
      character(len=default_str_len) :: vname
      logical :: lflg
      type(tVec) :: dumpvec

      call PetscOptionsHasName(PETSC_NULL_OPTIONS, solver%prefix, &
                               trim(dumpstring), lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        vname = trim(varname)
        if (present(opt_time)) vname = trim(vname)//'.t'//trim(adjustl(toStr(opt_time)))
        if (.not. lsolar .or. .not. lthermal) vname = trim(vname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
        call DMGetGlobalVector(dm, dumpvec, ierr); call CHKERR(ierr)
        call PetscObjectSetName(dumpvec, trim(vname), ierr); call CHKERR(ierr)
        call f90VecToPetsc(var, dm, dumpvec)

        call PetscObjectViewFromOptions(dumpvec, PETSC_NULL_VEC, &
                                        trim(dumpstring), ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(dm, dumpvec, ierr); call CHKERR(ierr)
      end if
    end subroutine

    subroutine dump_results()
      logical :: lflg
      character(len=default_str_len) :: fname
      integer(mpiint) :: ierr

      if (get_arg(.false., lonly_initialize)) return

      if (lsolar) then
        call dump_variable(edir, solver%C_one1%da, "-specint_dump_edir", "edir")
      end if
      call dump_variable(edn, solver%C_one1%da, "-specint_dump_edn", "edn")
      call dump_variable(eup, solver%C_one1%da, "-specint_dump_eup", "eup")
      call dump_variable(abso, solver%C_one%da, "-specint_dump_abso", "abso")

      fname = ''
      if (lsolar .and. present(opt_buildings_solar)) then
        call get_petsc_opt(solver%prefix, '-specint_xdmf_buildings_solar', fname, lflg, ierr); call CHKERR(ierr)
        if (lflg) then
          if (present(opt_time)) fname = trim(fname)//'.t'//trim(adjustl(toStr(opt_time)))
          if (.not. lsolar .or. .not. lthermal) fname = trim(fname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
          call xdmf_pprts_buildings(solver, opt_buildings_solar, fname, ierr, verbose=.true.); call CHKERR(ierr)
        end if
      end if

      if (lthermal .and. present(opt_buildings_thermal)) then
        call get_petsc_opt(solver%prefix, '-specint_xdmf_buildings_thermal', fname, lflg, ierr); call CHKERR(ierr)
        if (lflg) then
          if (present(opt_time)) fname = trim(fname)//'.t'//trim(adjustl(toStr(opt_time)))
          if (.not. lsolar .or. .not. lthermal) fname = trim(fname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
          call xdmf_pprts_buildings(solver, opt_buildings_thermal, fname, ierr, verbose=.true.); call CHKERR(ierr)
        end if
      end if

      call get_petsc_opt(solver%prefix, '-specint_xdmf', fname, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        if (present(opt_time)) fname = trim(fname)//'.t'//trim(adjustl(toStr(opt_time)))
        if (.not. lsolar .or. .not. lthermal) fname = trim(fname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
        if (lsolar) then
          call xdmf_pprts_srfc_flux(solver, fname, edn, eup, ierr, edir=edir, verbose=.true.); call CHKERR(ierr)
        else
          call xdmf_pprts_srfc_flux(solver, fname, edn, eup, ierr, verbose=.true.); call CHKERR(ierr)
        end if
      end if
    end subroutine
  end subroutine

  subroutine select_specint(specint, spec, ierr)
    character(len=*), intent(in) :: specint
    character(len=default_str_len), intent(out) :: spec
    integer(mpiint), intent(out) :: ierr

    logical :: lflg

    ierr = 0
    spec = trim(specint)
    call get_petsc_opt("", "-specint", spec, lflg, ierr); call CHKERR(ierr)

    select case (trim(spec))
    case ("rrtmg")
    case ("repwvl")
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'invalid specint mode <'//trim(spec)//'>'//new_line('')// &
        & "use one of:"//new_line('')// &
        & "  -specint rrtmg"//new_line('')// &
        & "  -specint repwvl"//new_line('')// &
        & "")
    end select
  end subroutine

  subroutine specint_pprts_destroy(specint, solver, lfinalizepetsc, ierr)
    character(len=*), intent(in) :: specint
    class(t_solver) :: solver
    logical, intent(in) :: lfinalizepetsc
    integer(mpiint), intent(out) :: ierr

    character(len=default_str_len) :: spec

    ierr = 0
    call select_specint(specint, spec, ierr); call CHKERR(ierr)

    select case (trim(spec))
    case ("rrtmg")
      call destroy_pprts_rrtmg(solver, lfinalizepetsc)
    case ("repwvl")
      call repwvl_pprts_destroy(solver, lfinalizepetsc, ierr)
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'invalid specint mode <'//trim(spec)//'>')
    end select
  end subroutine
end module
