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

!> \page Routines to call tenstream with optical properties from different spectral integration approaches

module m_specint_plexrt
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

  use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

  use m_plex_rt_base, only: t_plex_solver

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm

  use m_repwvl_plexrt, only: repwvl_plexrt, repwvl_plexrt_destroy

  implicit none

  private
  public :: specint_plexrt, specint_plexrt_destroy

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .false.
#else
  logical, parameter :: ldebug = .true.
#endif

contains

  subroutine specint_plexrt(specint,        &
      & solver, atm,                        &
      & sundir,                             &
      & albedo_thermal, albedo_solar,       &
      & lthermal, lsolar,                   &
      & edir, edn, eup, abso,               &
      & opt_time,                           &
      & solar_albedo_2d, thermal_albedo_2d, &
      & opt_solar_constant)

    character(len=*), intent(in) :: specint                  ! name of module to use for spectral integration
    class(t_plex_solver), allocatable, intent(inout) :: solver            ! solver type (e.g. t_solver_8_10)
    type(t_tenstr_atm), intent(in) :: atm                    ! contains info on atmospheric constituents
    real(ireals), intent(in) :: sundir(:)                    ! cartesian sun direction, pointing away from the sun, dim(3)
    real(ireals), intent(in) :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions
    ! and compute new solutions only after threshold estimate is exceeded.
    ! If solar_albedo_2d is present, we use a 2D surface albedo
    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:), thermal_albedo_2d(:), opt_solar_constant

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :), intent(inout) :: edir, edn, eup  ! [nlyr+1, ncol ]
    real(ireals), allocatable, dimension(:, :), intent(inout) :: abso            ! [nlyr  , ncol ]

    ! ---------- end of API ----------------

    integer(mpiint) :: ierr

    character(len=default_str_len) :: spec
    call select_specint(specint, spec, ierr)
    select case (trim(spec))
    case ("rrtmg")
      call plexrt_rrtmg(solver, atm, sundir, &
        & albedo_thermal, albedo_solar, &
        & lthermal, lsolar, &
        & edir, edn, eup, abso, &
        & opt_time, solar_albedo_2d, thermal_albedo_2d, &
        & opt_solar_constant)
    case ("repwvl")
      call repwvl_plexrt(solver, atm, sundir, &
        & albedo_thermal, albedo_solar, &
        & lthermal, lsolar, &
        & edir, edn, eup, abso, &
        & opt_time, solar_albedo_2d, thermal_albedo_2d, &
        & opt_solar_constant)
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'invalid specint mode <'//trim(spec)//'>')
    end select
  end subroutine

  subroutine specint_plexrt_destroy(specint, solver, lfinalizepetsc, ierr)
    character(len=*), intent(in) :: specint
    class(t_plex_solver), allocatable, intent(inout) :: solver
    logical, intent(in) :: lfinalizepetsc
    integer(mpiint), intent(out) :: ierr

    character(len=default_str_len) :: spec

    ierr = 0
    call select_specint(specint, spec, ierr); call CHKERR(ierr)

    select case (trim(spec))
    case ("rrtmg")
      call destroy_plexrt_rrtmg(solver, lfinalizepetsc)
    case ("repwvl")
      call repwvl_plexrt_destroy(solver, lfinalizepetsc, ierr)
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'invalid specint mode <'//trim(spec)//'>')
    end select
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
end module
