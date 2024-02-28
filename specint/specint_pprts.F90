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

module m_specint_pprts
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi, only: mpi_barrier

  use m_helper_functions, only: &
    & atol, &
    & CHKERR, &
    & domain_decompose_2d_petsc, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_bcast, &
    & toStr

  use m_data_parameters, only: &
    & iintegers, ireals, mpiint, &
    & default_str_len

  use m_pprts_base, only: t_solver

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm

  use m_buildings, only: t_pprts_buildings

  use m_pprts_rrtmg, only: pprts_rrtmg, destroy_pprts_rrtmg
  use m_repwvl_pprts, only: repwvl_pprts, repwvl_pprts_destroy
  use m_ecckd_pprts, only: ecckd_pprts, ecckd_pprts_destroy

  use m_xdmf_export, only: &
    & xdmf_pprts_buildings, &
    & xdmf_pprts_srfc_flux

  use m_petsc_helpers, only: f90vectopetsc

  use m_netcdfio, only: &
    & get_dim_info, &
    & get_global_attribute, &
    & ncload, &
    & ncwrite, &
    & set_global_attribute

  implicit none

  private
  public :: specint_pprts, specint_pprts_destroy, load_input_dump

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

    integer(mpiint) :: ierr

    character(len=default_str_len) :: spec

    call select_specint(specint, spec, ierr)
    select case (trim(spec))
    case ("ecckd")
      call ecckd_pprts(comm, solver, atm, ie, je,       &
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

    call dump_input()
    call dump_results()
  contains
    subroutine dump_input()
      logical :: lflg
      character(len=default_str_len) :: fname, groups(2), dimnames(3)
      integer :: Nlev, Nlay, local_shape(3), global_shape(3), startp(3)
      integer(mpiint) :: ierr

      call get_petsc_opt(solver%prefix, '-specint_dump_input', fname, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        if (present(opt_time)) fname = trim(fname)//'.t'//trim(adjustl(toStr(opt_time)))
        if (.not. lsolar .or. .not. lthermal) fname = trim(fname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
        fname = trim(fname)//".nc"

        groups(1) = trim(fname)
        if (solver%myid .eq. 0) then ! some variables are only written by one process
          print *, 'Dumping inputs to file:', trim(fname)
          call set_global_attribute(fname, 'dx', dx, ierr); call CHKERR(ierr)
          call set_global_attribute(fname, 'dy', dy, ierr); call CHKERR(ierr)
          call set_global_attribute(fname, 'albedo_thermal', albedo_thermal, ierr); call CHKERR(ierr)
          call set_global_attribute(fname, 'albedo_solar', albedo_solar, ierr); call CHKERR(ierr)
          call set_global_attribute(fname, 'lthermal', toStr(lthermal), ierr); call CHKERR(ierr)
          call set_global_attribute(fname, 'lsolar', toStr(lsolar), ierr); call CHKERR(ierr)
          groups(2) = "sundir"; call ncwrite(groups, sundir, ierr, dimnames=["xyz"]); call CHKERR(ierr)
          call set_global_attribute(fname, 'atm.atm_ke', atm%atm_ke, ierr); call CHKERR(ierr)
          call set_global_attribute(fname, 'atm.d_ke', atm%d_ke, ierr); call CHKERR(ierr)
          call set_global_attribute(fname, 'atm.d_ke1', atm%d_ke1, ierr); call CHKERR(ierr)

          if (present(opt_solar_constant)) then
            call set_global_attribute(fname, 'solar_constant', opt_solar_constant, ierr); call CHKERR(ierr)
          end if
          if (present(opt_time)) then
            call set_global_attribute(fname, 'time', opt_time, ierr); call CHKERR(ierr)
          end if
        end if
        call mpi_barrier(comm, ierr); call CHKERR(ierr) ! wait for 0-rank before opening the file with parallel netcdf

        Nlev = size(atm%plev, dim=1)
        Nlay = Nlev - 1
        dimnames(:) = [character(len=default_str_len) :: "zlev", "x", "y"]
        local_shape = [integer :: Nlev, ie, je]
        global_shape = [integer :: Nlev, solver%C_one1%glob_xm, solver%C_one1%glob_ym]
        startp = [integer :: 1, solver%C_one1%xs + 1, solver%C_one1%ys + 1]
        print *, solver%myid, 'nxproc', nxproc, 'nyproc', nyproc
        print *, solver%myid, 'local_shape', local_shape, 'global_shape', global_shape

        call dump_input_atm_var(comm, fname, 'atm.plev', atm%plev, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.tlev', atm%tlev, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.zlev', atm%zt, dimnames, local_shape, global_shape, startp)

        dimnames(:) = [character(len=default_str_len) :: "zlay", "x", "y"]
        local_shape = [integer :: Nlay, ie, je]
        global_shape = [integer :: Nlay, solver%C_one%glob_xm, solver%C_one%glob_ym]
        call dump_input_atm_var(comm, fname, 'atm.zm', atm%zm, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.dz', atm%dz, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.tlay', atm%tlay, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.h2o_lay', atm%h2o_lay, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.o3_lay', atm%o3_lay, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.co2_lay', atm%co2_lay, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.ch4_lay', atm%ch4_lay, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.n2o_lay', atm%n2o_lay, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.o2_lay', atm%o2_lay, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.lwc', atm%lwc, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.reliq', atm%reliq, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.iwc', atm%iwc, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.reice', atm%reice, dimnames, local_shape, global_shape, startp)
        call dump_input_atm_var(comm, fname, 'atm.cfrac', atm%cfrac, dimnames, local_shape, global_shape, startp)

        if (allocated(atm%opt_tau)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'atm.opt_tau'], &
            & arr=reshape(atm%opt_tau, [integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)]), &
            & ierr=ierr, &
            & arr_shape=[integer :: Nlay, solver%C_one%glob_xm, solver%C_one%glob_ym, size(atm%opt_tau, dim=3)], &
            & dimnames=[character(len=default_str_len) :: "zlay", "x", "y", "Nwvl"], &
            & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1, 1], &
            & countp=[integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)])
          call CHKERR(ierr)
        end if
        if (allocated(atm%opt_w0)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'atm.opt_w0'], &
            & arr=reshape(atm%opt_w0, [integer :: Nlay, ie, je, size(atm%opt_w0, dim=3)]), &
            & ierr=ierr, &
            & arr_shape=[integer :: Nlay, solver%C_one%glob_xm, solver%C_one%glob_ym, size(atm%opt_w0, dim=3)], &
            & dimnames=[character(len=default_str_len) :: "zlay", "x", "y", "Nwvl"], &
            & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1, 1], &
            & countp=[integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)])
          call CHKERR(ierr)
        end if
        if (allocated(atm%opt_g)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'atm.opt_g'], &
            & arr=reshape(atm%opt_g, [integer :: Nlay, ie, je, size(atm%opt_g, dim=3)]), &
            & ierr=ierr, &
            & arr_shape=[integer :: Nlay, solver%C_one%glob_xm, solver%C_one%glob_ym, size(atm%opt_g, dim=3)], &
            & dimnames=[character(len=default_str_len) :: "zlay", "x", "y", "Nwvl"], &
            & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1, 1], &
            & countp=[integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)])
          call CHKERR(ierr)
        end if
        if (allocated(atm%tskin)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'atm.tskin'], &
            & arr=reshape(atm%tskin, [integer :: ie, je]), &
            & ierr=ierr, &
            & arr_shape=[integer :: solver%C_one%glob_xm, solver%C_one%glob_ym], &
            & dimnames=[character(len=default_str_len) :: "x", "y"], &
            & startp=[integer :: solver%C_one%xs + 1, solver%C_one%ys + 1], &
            & countp=[integer :: ie, je])
          call CHKERR(ierr)
        end if

        if (present(solar_albedo_2d)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'solar_albedo_2d'], &
            & arr=solar_albedo_2d, &
            & ierr=ierr, &
            & arr_shape=global_shape(2:3), &
            & dimnames=dimnames(2:3), &
            & startp=[integer :: solver%C_one%xs + 1, solver%C_one%ys + 1], &
            & countp=[integer :: ie, je])
          call CHKERR(ierr)
        end if
        if (present(thermal_albedo_2d)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'thermal_albedo_2d'], &
            & arr=thermal_albedo_2d, &
            & ierr=ierr, &
            & arr_shape=global_shape(2:3), &
            & dimnames=dimnames(2:3), &
            & startp=[integer :: solver%C_one%xs + 1, solver%C_one%ys + 1], &
            & countp=[integer :: ie, je])
          call CHKERR(ierr)
        end if
        if (present(opt_tau_solar)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'opt_tau_solar'], &
            & arr=opt_tau_solar, &
            & ierr=ierr, &
            & arr_shape=global_shape, &
            & dimnames=dimnames, &
            & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1], &
            & countp=shape(opt_tau_solar))
          call CHKERR(ierr)
        end if
        if (present(opt_w0_solar)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'opt_w0_solar'], &
            & arr=opt_w0_solar, &
            & ierr=ierr, &
            & arr_shape=global_shape, &
            & dimnames=dimnames, &
            & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1], &
            & countp=shape(opt_w0_solar))
          call CHKERR(ierr)
        end if
        if (present(opt_g_solar)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'opt_g_solar'], &
            & arr=opt_g_solar, &
            & ierr=ierr, &
            & arr_shape=global_shape, &
            & dimnames=dimnames, &
            & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1], &
            & countp=shape(opt_g_solar))
          call CHKERR(ierr)
        end if
        if (present(opt_tau_thermal)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: fname, 'opt_tau_thermal'], &
            & arr=opt_tau_thermal, &
            & ierr=ierr, &
            & arr_shape=global_shape, &
            & dimnames=dimnames, &
            & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1], &
            & countp=shape(opt_tau_thermal))
          call CHKERR(ierr)
        end if

      end if
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    end subroutine

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

  subroutine dump_input_atm_var(comm, fname, varname, arr, dimnames, local_shape, global_shape, startp)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: fname, varname, dimnames(:)
    real(ireals), allocatable, intent(in) :: arr(:, :)
    integer, intent(in) :: local_shape(3), global_shape(3), startp(3)
    integer(mpiint) :: ierr
    if (allocated(arr)) then
      call ncwrite(&
        & comm, &
        & groups=[character(len=default_str_len) :: fname, varname], &
        & arr=reshape(arr, local_shape), &
        & ierr=ierr, &
        & arr_shape=global_shape, &
        & dimnames=dimnames, &
        & startp=startp, &
        & countp=local_shape)
      call CHKERR(ierr)
    end if
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
    case ("ecckd")
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'invalid specint mode <'//trim(spec)//'>'//new_line('')// &
        & "use one of:"//new_line('')// &
        & "  -specint rrtmg"//new_line('')// &
        & "  -specint repwvl"//new_line('')// &
        & "  -specint ecckd"//new_line('')// &
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
    case ("ecckd")
      call ecckd_pprts_destroy(solver, lfinalizepetsc, ierr); call CHKERR(ierr)
    case ("rrtmg")
      call destroy_pprts_rrtmg(solver, lfinalizepetsc)
    case ("repwvl")
      call repwvl_pprts_destroy(solver, lfinalizepetsc, ierr); call CHKERR(ierr)
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'invalid specint mode <'//trim(spec)//'>')
    end select
  end subroutine

  subroutine load_input_dump(&
      & comm, &
      & inpfile, &
      & Nx_local, Ny_local, &
      & nxproc, nyproc, &
      & dx, dy, &
      & sundir, &
      & albedo_thermal, albedo_solar, &
      & lsolar, lthermal, &
      & atm, &
      & opt_time, &
      & solar_albedo_2d, thermal_albedo_2d, &
      & opt_solar_constant, &
      & ierr)
    integer(mpiint), intent(in) :: comm
    character(len=default_str_len) :: inpfile
    integer(iintegers), intent(out) :: Nx_local, Ny_local
    integer(iintegers), allocatable, intent(out) :: nxproc(:), nyproc(:)
    real(ireals), intent(out) :: dx, dy
    real(ireals), allocatable, intent(out) :: sundir(:)
    real(ireals), intent(out) :: albedo_solar, albedo_thermal
    logical, intent(out) :: lsolar, lthermal
    type(t_tenstr_atm), intent(out) :: atm
    real(ireals), allocatable, intent(out) :: opt_time, solar_albedo_2d(:, :), thermal_albedo_2d(:, :), opt_solar_constant
    integer(mpiint), intent(out) :: ierr

    integer(mpiint) :: myid
    integer(iintegers) :: Nx_global, Ny_global, xStart, yStart
    integer(iintegers) :: Nlev, Nlay
    character(len=default_str_len) :: lthermal_str, lsolar_str
    logical :: lhave_opt_time, lhave_opt_solar_constant

    integer :: ostart(4), ocount(4)

    call mpi_comm_rank(comm, myid, ierr)

    if (myid .eq. 0) then
      call get_dim_info(inpfile, 'x', ierr, dimsize=Nx_global); call CHKERR(ierr)
      call get_dim_info(inpfile, 'y', ierr, dimsize=Ny_global); call CHKERR(ierr)
    end if
    call imp_bcast(comm, Nx_global, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Ny_global, 0_mpiint, ierr); call CHKERR(ierr)

    call domain_decompose_2d_petsc(comm, Nx_global, Ny_global, &
      & Nx_local, Ny_local, &
      & xStart, yStart, &
      & nxproc, nyproc, ierr)
    call CHKERR(ierr)

    if (myid .eq. 0) print *, 'Domain Decomposition will be', nxproc, 'and', nyproc, 'xs', xStart, 'ys', yStart

    if (myid .eq. 0) then
      call get_global_attribute(inpfile, 'dx', dx, ierr); call CHKERR(ierr)
      call get_global_attribute(inpfile, 'dy', dy, ierr); call CHKERR(ierr)
      call get_global_attribute(inpfile, 'albedo_thermal', albedo_thermal, ierr); call CHKERR(ierr)
      call get_global_attribute(inpfile, 'albedo_solar', albedo_solar, ierr); call CHKERR(ierr)
      call get_global_attribute(inpfile, 'lthermal', lthermal_str, ierr); call CHKERR(ierr)
      lthermal = atol(lthermal_str)
      call get_global_attribute(inpfile, 'lsolar', lsolar_str, ierr); call CHKERR(ierr)
      lsolar = atol(lsolar_str)

      allocate (opt_solar_constant)
      call get_global_attribute(inpfile, 'solar_constant', opt_solar_constant, ierr)
      if (ierr .ne. 0_mpiint) deallocate (opt_solar_constant)
      lhave_opt_solar_constant = allocated(opt_solar_constant)

      allocate (opt_time)
      call get_global_attribute(inpfile, 'time', opt_time, ierr)
      if (ierr .ne. 0_mpiint) deallocate (opt_time)
      lhave_opt_time = allocated(opt_time)

      call ncload([character(default_str_len) :: inpfile, 'sundir'], sundir, ierr); call CHKERR(ierr)
    end if
    call imp_bcast(comm, dx, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dy, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, albedo_thermal, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, albedo_solar, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, lthermal, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, lsolar, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, sundir, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, lhave_opt_solar_constant, 0_mpiint, ierr); call CHKERR(ierr)
    if (lhave_opt_solar_constant) then
      if (.not. allocated(opt_solar_constant)) allocate (opt_solar_constant)
      call imp_bcast(comm, opt_solar_constant, 0_mpiint, ierr); call CHKERR(ierr)
    end if
    call imp_bcast(comm, lhave_opt_time, 0_mpiint, ierr); call CHKERR(ierr)
    if (lhave_opt_time) then
      if (.not. allocated(opt_time)) allocate (opt_time)
      call imp_bcast(comm, opt_time, 0_mpiint, ierr); call CHKERR(ierr)
    end if

    print *, myid, 'dx', dx
    print *, myid, 'dy', dy
    print *, myid, 'albedo_thermal', albedo_thermal
    print *, myid, 'albedo_solar', albedo_solar
    print *, myid, 'lthermal', lthermal
    print *, myid, 'lsolar', lsolar
    print *, myid, 'sundir', sundir
    print *, myid, 'have_opt_solar_constant', allocated(opt_solar_constant)
    if (lhave_opt_solar_constant) print *, myid, 'opt_solar_constant', opt_solar_constant
    print *, myid, 'have_opt_time', allocated(opt_time)
    if (lhave_opt_time) print *, myid, 'opt_time', opt_time

    ostart(1:2) = [integer :: xStart + 1, yStart + 1]
    ocount(1:2) = [integer :: Nx_local, Ny_local]
    call ncload(comm, [character(default_str_len) :: inpfile, 'solar_albedo_2d'], solar_albedo_2d, ierr, &
      & ostart=ostart, ocount=ocount)
    print *, myid, 'have_solar_albedo_2d', allocated(solar_albedo_2d)
    if (allocated(solar_albedo_2d)) print *, myid, 'solar_albedo_2d', solar_albedo_2d
    call ncload(comm, [character(default_str_len) :: inpfile, 'thermal_albedo_2d'], thermal_albedo_2d, ierr, &
      & ostart=ostart, ocount=ocount)
    print *, myid, 'have_thermal_albedo_2d', allocated(thermal_albedo_2d)
    if (allocated(thermal_albedo_2d)) print *, myid, 'thermal_albedo_2d', thermal_albedo_2d

    ! populate atm
    allocate (atm%atm_ke)
    if (myid .eq. 0) then
      call get_global_attribute(inpfile, 'atm.atm_ke', atm%atm_ke, ierr); call CHKERR(ierr)
      call get_global_attribute(inpfile, 'atm.d_ke', atm%d_ke, ierr); call CHKERR(ierr)
      call get_global_attribute(inpfile, 'atm.d_ke1', atm%d_ke1, ierr); call CHKERR(ierr)
      call get_dim_info(inpfile, 'zlev', ierr, dimsize=Nlev); call CHKERR(ierr)
      call get_dim_info(inpfile, 'zlay', ierr, dimsize=Nlay); call CHKERR(ierr)
    end if
    call imp_bcast(comm, atm%atm_ke, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, atm%d_ke, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, atm%d_ke1, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Nlev, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Nlay, 0_mpiint, ierr); call CHKERR(ierr)
    print *, myid, 'Nlev', Nlev
    print *, myid, 'Nlay', Nlay

    ostart(1:3) = [integer :: 1, xStart + 1, yStart + 1]
    ocount(1:3) = [integer :: Nlev, Nx_local, Ny_local]
    call load_input_atm_var('atm.plev', atm%plev, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.tlev', atm%tlev, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.zlev', atm%zt, ierr); call CHKERR(ierr)

    ocount(1:3) = [integer :: Nlay, Nx_local, Ny_local]
    call load_input_atm_var('atm.zm', atm%zm, ierr)!; call CHKERR(ierr)
    call load_input_atm_var('atm.dz', atm%dz, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.tlay', atm%tlay, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.h2o_lay', atm%h2o_lay, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.o3_lay', atm%o3_lay, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.co2_lay', atm%co2_lay, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.ch4_lay', atm%ch4_lay, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.n2o_lay', atm%n2o_lay, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.o2_lay', atm%o2_lay, ierr); call CHKERR(ierr)
    call load_input_atm_var('atm.lwc', atm%lwc, ierr)!; call CHKERR(ierr)
    call load_input_atm_var('atm.reliq', atm%reliq, ierr)!; call CHKERR(ierr)
    call load_input_atm_var('atm.iwc', atm%iwc, ierr)!; call CHKERR(ierr)
    call load_input_atm_var('atm.reice', atm%reice, ierr)!; call CHKERR(ierr)
    call load_input_atm_var('atm.cfrac', atm%cfrac, ierr)!; call CHKERR(ierr)

    !TODO missing loading:
    ! atm.opt_tau
    ! atm.opt_w0
    ! atm.opt_g
    ! atm.opt_tskin
    ! atm.opt_tau_solar
    ! atm.opt_w0_solar
    ! atm.opt_g_solar
    ! atm.opt_tau_thermal

    ierr = 0
  contains
    subroutine load_input_atm_var(varname, arr, ierr)
      character(len=*), intent(in) :: varname
      real(ireals), allocatable, intent(out) :: arr(:, :)
      integer(mpiint), intent(out) :: ierr

      real(ireals), allocatable :: arr3d(:, :, :)
      call ncload(comm, [character(default_str_len) :: inpfile, varname], arr3d, ierr, ostart=ostart, ocount=ocount)
      if (allocated(arr3d)) then
        allocate (arr(size(arr3d, dim=1), size(arr3d, dim=2) * size(arr3d, dim=3)))
        arr = reshape(arr3d, [size(arr3d, dim=1), size(arr3d, dim=2) * size(arr3d, dim=3)])
        print *, 'read ', trim(varname)
      else
        print *, 'no input for ', trim(varname)
      end if
    end subroutine
  end subroutine
end module
