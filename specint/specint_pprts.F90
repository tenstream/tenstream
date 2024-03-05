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
  use mpi

  use m_helper_functions, only: &
    & atol, &
    & CHKERR, &
    & cumsum, &
    & domain_decompose_2d_petsc, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_allreduce_max, &
    & imp_allreduce_sum, &
    & imp_bcast, &
    & imp_scan_sum, &
    & ind_1d_to_nd, &
    & toStr

  use m_data_parameters, only: &
    & iintegers, ireals, mpiint, &
    & imp_iinteger, &
    & imp_ireals, &
    & default_str_len

  use m_pprts_base, only: t_solver

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm

  use m_buildings, only: &
    & check_buildings_consistency, &
    & faceidx_by_cell_plus_offset, &
    & init_buildings, &
    & t_pprts_buildings

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
      character(len=default_str_len) :: fname
      logical :: lflg
      integer(mpiint) :: ierr
      logical :: in_time_range
      integer(iintegers) :: i, timerange_nargs
      real(ireals), allocatable :: timerange(:)

      call get_petsc_opt(solver%prefix, '-specint_dump_input', &
        & fname, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        in_time_range = .true. ! default is to assume we want the dump

        if (present(opt_time)) then
          timerange_nargs = 2048
          allocate (timerange(timerange_nargs))
          call get_petsc_opt(solver%prefix, '-specint_dump_input_time_range', &
            & timerange, timerange_nargs, lflg, ierr); call CHKERR(ierr)
          if (lflg) then
            call CHKERR(modulo(timerange_nargs, 2_iintegers), "timerange number of arguments must be multiple of two, &
              & i.e. pairs of start/end. Found "//toStr(timerange_nargs)//" arguments: "//toStr(timerange(1:timerange_nargs)))
            in_time_range = .false.
            do i = 1, timerange_nargs - 1, 2
              if (timerange(i) .gt. timerange(i + 1)) &
                & call CHKERR(int(i, mpiint), 'time range arguments have to ascending but '// &
                & 'arg nr #'//toStr(i)//' ('//toStr(timerange(i))//') is larger than '// &
                & 'arg nr #'//toStr(i + 1)//' ('//toStr(timerange(i + 1))//')')
              if (opt_time .ge. timerange(i) .and.&
                & opt_time .le. timerange(i + 1)) in_time_range = .true.
              if (in_time_range) exit
            end do
          end if
        end if

        if (in_time_range) then
          call dump_input_internal(fname)
        else
          if (solver%myid .eq. 0) print *, 'Skipping specint_dump_input because not in time range'
        end if
      end if

    end subroutine

    subroutine dump_input_internal(fname_in)
      character(len=default_str_len), intent(in) :: fname_in
      character(len=default_str_len) :: fname, groups(2), dimnames(3)
      integer :: Nlev, Nlay, Nveg, local_shape(3), global_shape(3), startp(3)
      integer(mpiint) :: ierr
      logical, parameter :: lverbose = .false.

      fname = trim(fname_in)

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
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'atm.opt_tau'], &
          & arr=reshape(atm%opt_tau, [integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)]), &
          & ierr=ierr, &
          & arr_shape=[integer :: Nlay, solver%C_one%glob_xm, solver%C_one%glob_ym, size(atm%opt_tau, dim=3)], &
          & dimnames=[character(len=default_str_len) :: "zlay", "x", "y", "Nwvl"], &
          & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1, 1], &
          & countp=[integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)], &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if
      if (allocated(atm%opt_w0)) then
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'atm.opt_w0'], &
          & arr=reshape(atm%opt_w0, [integer :: Nlay, ie, je, size(atm%opt_w0, dim=3)]), &
          & ierr=ierr, &
          & arr_shape=[integer :: Nlay, solver%C_one%glob_xm, solver%C_one%glob_ym, size(atm%opt_w0, dim=3)], &
          & dimnames=[character(len=default_str_len) :: "zlay", "x", "y", "Nwvl"], &
          & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1, 1], &
          & countp=[integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)], &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if
      if (allocated(atm%opt_g)) then
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'atm.opt_g'], &
          & arr=reshape(atm%opt_g, [integer :: Nlay, ie, je, size(atm%opt_g, dim=3)]), &
          & ierr=ierr, &
          & arr_shape=[integer :: Nlay, solver%C_one%glob_xm, solver%C_one%glob_ym, size(atm%opt_g, dim=3)], &
          & dimnames=[character(len=default_str_len) :: "zlay", "x", "y", "Nwvl"], &
          & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1, 1], &
          & countp=[integer :: Nlay, ie, je, size(atm%opt_tau, dim=3)], &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if
      if (allocated(atm%tskin)) then
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'atm.tskin'], &
          & arr=reshape(atm%tskin, [integer :: ie, je]), &
          & ierr=ierr, &
          & arr_shape=[integer :: solver%C_one%glob_xm, solver%C_one%glob_ym], &
          & dimnames=[character(len=default_str_len) :: "x", "y"], &
          & startp=[integer :: solver%C_one%xs + 1, solver%C_one%ys + 1], &
          & countp=[integer :: ie, je], &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if

      if (present(solar_albedo_2d)) then
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'solar_albedo_2d'], &
          & arr=solar_albedo_2d, &
          & ierr=ierr, &
          & arr_shape=global_shape(2:3), &
          & dimnames=dimnames(2:3), &
          & startp=[integer :: solver%C_one%xs + 1, solver%C_one%ys + 1], &
          & countp=[integer :: ie, je], &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if
      if (present(thermal_albedo_2d)) then
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'thermal_albedo_2d'], &
          & arr=thermal_albedo_2d, &
          & ierr=ierr, &
          & arr_shape=global_shape(2:3), &
          & dimnames=dimnames(2:3), &
          & startp=[integer :: solver%C_one%xs + 1, solver%C_one%ys + 1], &
          & countp=[integer :: ie, je], &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if

      dimnames(:) = [character(len=default_str_len) :: "zlay_veg", "x", "y"]
      call imp_allreduce_max(comm, size(opt_tau_solar, dim=1, kind=iintegers), Nveg)
      global_shape = [integer :: Nveg, solver%C_one%glob_xm, solver%C_one%glob_ym]
      if (present(opt_tau_solar)) then
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'opt_tau_solar'], &
          & arr=opt_tau_solar, &
          & ierr=ierr, &
          & arr_shape=global_shape, &
          & dimnames=dimnames, &
          & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1], &
          & countp=shape(opt_tau_solar), &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if
      if (present(opt_w0_solar)) then
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'opt_w0_solar'], &
          & arr=opt_w0_solar, &
          & ierr=ierr, &
          & arr_shape=global_shape, &
          & dimnames=dimnames, &
          & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1], &
          & countp=shape(opt_w0_solar), &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if
      if (present(opt_g_solar)) then
        call ncwrite(&
          & comm=comm, &
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
          & comm=comm, &
          & groups=[character(len=default_str_len) :: fname, 'opt_tau_thermal'], &
          & arr=opt_tau_thermal, &
          & ierr=ierr, &
          & arr_shape=global_shape, &
          & dimnames=dimnames, &
          & startp=[integer :: 1, solver%C_one%xs + 1, solver%C_one%ys + 1], &
          & countp=shape(opt_tau_thermal), &
          & verbose=lverbose)
        call CHKERR(ierr)
      end if

      if (present(opt_buildings_solar)) then
        call dump_input_buildings(comm, fname, opt_buildings_solar, 'solar', &
          & solver%C_one%zm, solver%C_one%xs, solver%C_one%ys, &
          & ierr); call CHKERR(ierr)
      end if
      if (present(opt_buildings_thermal)) then
        call dump_input_buildings(comm, fname, opt_buildings_thermal, 'thermal', &
          & solver%C_one%zm, solver%C_one%xs, solver%C_one%ys, &
          & ierr); call CHKERR(ierr)
      end if

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
        & comm=comm, &
        & groups=[character(len=default_str_len) :: fname, varname], &
        & arr=reshape(arr, local_shape), &
        & ierr=ierr, &
        & arr_shape=global_shape, &
        & dimnames=dimnames, &
        & startp=startp, &
        & countp=local_shape, &
        & deflate_lvl=0, &
        & verbose=.false.)
      call CHKERR(ierr)
    end if
  end subroutine

  subroutine dump_input_buildings(comm, fname, buildings, prefix, Nlay, xs, ys, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: fname, prefix
    type(t_pprts_buildings), intent(in) :: buildings
    integer(iintegers), intent(in) :: Nlay, xs, ys
    integer(mpiint), intent(out) :: ierr

    integer(mpiint) :: myid
    integer(iintegers) :: Nlocal, Nglobal, bStart
    integer(iintegers) :: m, idx(4)
    integer(iintegers), allocatable :: bldg_idx(:, :)
    character(len=default_str_len) :: groups(4), dimnames(2)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    associate (B => buildings)

      if (.not. associated(B%iface)) call CHKERR(1_mpiint, 'buildings%iface not associated! ('//trim(prefix)//')')
      Nlocal = size(B%iface)

      !print *,myid,'have', Nlocal, 'building faces ('//trim(prefix)//')'
      call imp_scan_sum(comm, Nlocal, bStart, ierr); call CHKERR(ierr)
      bStart = 1 + bStart - Nlocal
      !print *,myid,'my_startidx', bStart, 'building faces ('//trim(prefix)//')'

      call imp_allreduce_sum(comm, Nlocal, Nglobal)
      !print *,myid,'global', Nglobal, 'building faces ('//trim(prefix)//')'

      ! convert from local to global indices
      allocate (bldg_idx(4, Nlocal))
      do m = 1, size(B%iface)
        call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
        idx(3:4) = idx(3:4) + [xs, ys]
        idx(2) = Nlay - idx(2) + 1 ! count from bottom
        bldg_idx(:, m) = idx
        print *, myid, bStart + m, 'bldg', bldg_idx(:, m)
      end do

      groups = [character(len=default_str_len) :: fname, 'buildings', prefix, 'idx']
      dimnames = [character(len=default_str_len) :: 'fkij', 'Nbldgs.'//prefix]
      call ncwrite(&
        & comm=comm, &
        & groups=groups, &
        & arr=bldg_idx, &
        & ierr=ierr, &
        & arr_shape=[integer :: size(bldg_idx, dim=1), Nglobal], &
        & dimnames=dimnames, &
        & startp=[integer :: 1, bStart], &
        & countp=shape(bldg_idx), &
        & deflate_lvl=0, &
        & verbose=.false.)
      call CHKERR(ierr)

      dimnames(1:1) = 'Nbldgs.'//trim(prefix)
      if (allocated(B%albedo)) then
        print *, 'B%albedo', B%albedo
        groups(4) = 'albedo'
        call ncwrite(&
          & comm=comm, &
          & groups=groups, &
          & arr=B%albedo, &
          & ierr=ierr, &
          & arr_shape=[integer :: Nglobal], &
          & dimnames=dimnames, &
          & startp=[integer :: bStart], &
          & countp=shape(B%albedo), &
          & deflate_lvl=0, &
          & verbose=.true.)
        call CHKERR(ierr)
      end if

      if (allocated(B%temp)) then
        groups(4) = 'temp'
        call ncwrite(&
          & comm=comm, &
          & groups=groups, &
          & arr=B%temp, &
          & ierr=ierr, &
          & arr_shape=[integer :: Nglobal], &
          & dimnames=dimnames, &
          & startp=[integer :: bStart], &
          & countp=shape(B%temp), &
          & deflate_lvl=0, &
          & verbose=.false.)
        call CHKERR(ierr)
      end if

      if (allocated(B%planck)) then
        groups(4) = 'planck'
        call ncwrite(&
          & comm=comm, &
          & groups=groups, &
          & arr=B%planck, &
          & ierr=ierr, &
          & arr_shape=[integer :: Nglobal], &
          & dimnames=dimnames, &
          & startp=[integer :: bStart], &
          & countp=shape(B%planck), &
          & deflate_lvl=0, &
          & verbose=.false.)
        call CHKERR(ierr)
      end if

    end associate
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
      & opt_buildings_solar, opt_buildings_thermal, &
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
    type(t_pprts_buildings), allocatable, intent(out) :: opt_buildings_solar, opt_buildings_thermal
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

    ostart(1:2) = [integer :: xStart + 1, yStart + 1]
    ocount(1:2) = [integer :: Nx_local, Ny_local]
    call ncload([character(default_str_len) :: inpfile, 'solar_albedo_2d'], solar_albedo_2d, ierr, &
      & comm=comm, ostart=ostart, ocount=ocount)
    call ncload([character(default_str_len) :: inpfile, 'thermal_albedo_2d'], thermal_albedo_2d, ierr, &
      & comm=comm, ostart=ostart, ocount=ocount)

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

    call load_buildings_info(&
      & comm, inpfile, 'solar', &
      & Nlev - 1, Nx_local, Ny_local, &
      & nxproc, nyproc, &
      & opt_buildings_solar, ierr)

    call load_buildings_info(&
      & comm, inpfile, 'thermal', &
      & Nlev - 1, Nx_local, Ny_local, &
      & nxproc, nyproc, &
      & opt_buildings_solar, ierr)

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
      call ncload([character(default_str_len) :: inpfile, varname], arr3d, ierr, comm=comm, ostart=ostart, ocount=ocount)
      if (allocated(arr3d)) then
        allocate (arr(size(arr3d, dim=1), size(arr3d, dim=2) * size(arr3d, dim=3)))
        arr = reshape(arr3d, [size(arr3d, dim=1), size(arr3d, dim=2) * size(arr3d, dim=3)])
        if (myid .eq. 0) print *, 'read ', trim(varname)
      else
        if (myid .eq. 0) print *, 'no input for ', trim(varname)
      end if
    end subroutine

    subroutine find_proc(nxproc_cum, nyproc_cum, glob_i, glob_j, rank, loc_i, loc_j, ierr)
      integer(iintegers), intent(in) :: nxproc_cum(:), nyproc_cum(:), glob_i, glob_j
      integer(iintegers), intent(out) :: rank, loc_i, loc_j
      integer(mpiint), intent(out) :: ierr
      integer(iintegers) :: rank_i, rank_j

      do rank_i = 0, size(nxproc_cum) - 1
        if (nxproc_cum(rank_i + 1) .ge. glob_i) exit
      end do
      do rank_j = 0, size(nyproc_cum) - 1
        if (nyproc_cum(rank_j + 1) .ge. glob_j) exit
      end do
      rank = (rank_j) * size(nxproc_cum) + (rank_i)

      ierr = 0
      if (rank_i .ge. size(nxproc_cum)) ierr = ierr + 1
      if (rank_j .ge. size(nyproc_cum)) ierr = ierr + 1
      if (rank .gt. size(nyproc_cum) * size(nxproc_cum) - 1) ierr = ierr + 1
      if (ierr .ne. 0) then
        print *, 'glob_i, glob_j', glob_i, glob_j, 'lives on rank', rank_i, rank_j, "->", rank
        call CHKERR(ierr)
      end if

      if (rank_i .gt. 0) then
        loc_i = glob_i - nxproc_cum(rank_i)
      else
        loc_i = glob_i
      end if
      if (rank_j .gt. 0) then
        loc_j = glob_j - nyproc_cum(rank_j)
      else
        loc_j = glob_j
      end if
      !print *, 'glob_i, glob_j', glob_i, glob_j, 'lives on rank', rank_i, rank_j, "->", rank, &
      !  & 'loc_i', loc_i, 'loc_j', loc_j
    end subroutine

    subroutine load_buildings_info(comm, inpfile, prefix, Nlay, Nx, Ny, nxproc, nyproc, buildings, ierr)
      integer(mpiint), intent(in) :: comm
      character(len=*), intent(in) :: inpfile, prefix
      integer(iintegers), intent(in) :: Nlay, Nx, Ny, nxproc(:), nyproc(:)
      type(t_pprts_buildings), allocatable, intent(out) :: buildings
      integer(mpiint), intent(out) :: ierr

      integer(mpiint) :: myid, numnodes
      integer(iintegers), allocatable :: orig_idx(:, :) ! from netcdf
      integer(iintegers), allocatable :: idx(:, :) ! sorted by ranks
      integer(iintegers), allocatable :: my_idx(:, :) ! local_values after send
      integer(iintegers), allocatable :: netcdf_idx_to_rank_sorted(:)
      integer(iintegers), allocatable :: nxproc_cum(:), nyproc_cum(:)
      integer(iintegers), allocatable :: ioff_per_rank(:)
      integer(mpiint), allocatable :: send_cnts(:), displs(:)

      integer(iintegers) :: Nfaces, i, j, m !, fidx(4)
      integer(mpiint) :: rank, my_face_cnt
      logical :: lhave_building_data

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      allocate (send_cnts(0:numnodes - 1))
      allocate (displs(0:numnodes - 1))

      if (myid .eq. 0) then
        call ncload([character(len=default_str_len) :: inpfile, 'buildings', prefix, 'idx'], orig_idx, ierr)
        lhave_building_data = ierr .eq. 0
      end if
      call imp_bcast(comm, lhave_building_data, 0_mpiint, ierr); call CHKERR(ierr)
      if (.not. lhave_building_data) return

      if (myid .eq. 0) then
        Nfaces = size(orig_idx, dim=2)

        ! cumulative sum of of local domain sizes:
        allocate (nxproc_cum(size(nxproc)))
        nxproc_cum = cumsum(nxproc)
        allocate (nyproc_cum(size(nyproc)))
        nyproc_cum = cumsum(nyproc)

        send_cnts(:) = 0
        do m = 1, Nfaces
          call find_proc( &
            & nxproc_cum, nyproc_cum, &
            & orig_idx(3, m), orig_idx(4, m), &
            & rank, i, j, ierr); call CHKERR(ierr)

          send_cnts(rank) = send_cnts(rank) + 1
        end do

        displs(0) = 0
        if (numnodes .gt. 1_mpiint) then
          displs(1:numnodes - 1) = cumsum(send_cnts(0:numnodes - 2))
        end if

        allocate (ioff_per_rank(0:numnodes - 1))
        ioff_per_rank = displs

        !do m=0,numnodes-1
        !  print *,'rank', m, 'facecnt', send_cnts(m), 'displs', displs(m), 'cum_facecnt_per_rank', cum_facecnt_per_rank(m)
        !enddo

        allocate (netcdf_idx_to_rank_sorted(Nfaces))
        allocate (idx(4, Nfaces))
        do m = 1, Nfaces
          call find_proc( &
            & nxproc_cum, nyproc_cum, &
            & orig_idx(3, m), orig_idx(4, m), &
            & rank, i, j, ierr); call CHKERR(ierr)

          ioff_per_rank(rank) = ioff_per_rank(rank) + 1
          netcdf_idx_to_rank_sorted(m) = ioff_per_rank(rank)
          idx(:, netcdf_idx_to_rank_sorted(m)) = [orig_idx(1, m), orig_idx(2, m), i, j]
        end do
        !do m = 1, Nfaces
        !  print *, m, netcdf_idx_to_rank_sorted(m), 'idx', idx(:, netcdf_idx_to_rank_sorted(m))
        !end do
        deallocate (orig_idx)
      else
        allocate (idx(4, 0))
        allocate (netcdf_idx_to_rank_sorted(0))
        send_cnts(:) = 0
        displs(:) = 0
      end if
      call mpi_scatter(send_cnts, 1_mpiint, imp_iinteger, Nfaces, 1_mpiint, imp_iinteger, 0_mpiint, comm, ierr)
      call CHKERR(ierr)
      my_face_cnt = int(Nfaces, mpiint)

      call init_buildings(buildings, &
        & [integer(iintegers) :: 6, Nlay, Nx, Ny], &
        & my_face_cnt, ierr); call CHKERR(ierr)

      allocate (my_idx(4, my_face_cnt))
      call MPI_Scatterv(&
        & idx, &
        & send_cnts * 4, &
        & displs * 4, &
        & imp_iinteger, &
        & my_idx, &
        & my_face_cnt * 4, &
        & imp_iinteger, &
        & 0_mpiint, &
        & comm, &
        & ierr)
      call CHKERR(ierr)
      deallocate (idx)

      do m = 1, size(my_idx, dim=2)
        associate (d => my_idx(1, m), k => Nlay - my_idx(2, m) + 1, i => my_idx(3, m), j => my_idx(4, m))
          buildings%iface(m) = faceidx_by_cell_plus_offset(buildings%da_offsets, k, i, j, d)
          !print *, myid, trim(prefix), 'now have building at ', my_idx(:, m), 'kij', k, i, j, 'iface', buildings%iface(m)
        end associate
      end do

      call check_buildings_consistency(buildings, Nlay, Nx, Ny, ierr); call CHKERR(ierr)

      call scatter_buildings_var(comm, inpfile, prefix, 'albedo', &
        & netcdf_idx_to_rank_sorted, send_cnts, displs, my_face_cnt, &
        & buildings%albedo, ierr); call CHKERR(ierr)

      call scatter_buildings_var(comm, inpfile, prefix, 'temp', &
        & netcdf_idx_to_rank_sorted, send_cnts, displs, my_face_cnt, &
        & buildings%temp, ierr); call CHKERR(ierr)

      call scatter_buildings_var(comm, inpfile, prefix, 'planck', &
        & netcdf_idx_to_rank_sorted, send_cnts, displs, my_face_cnt, &
        & buildings%planck, ierr); call CHKERR(ierr)

      !do m = 1, size(buildings%iface)
      !  call ind_1d_to_nd(buildings%da_offsets, buildings%iface(m), fidx)
      !  associate (d => fidx(1), k => fidx(2), i => fidx(3), j => fidx(4))
      !    if (allocated(buildings%temp)) then
      !      print *, myid, 'now have building at ', my_idx(:, m), 'iface', buildings%iface(m), &
      !        & 'Ag', buildings%albedo(m), 'temp', buildings%temp(m)
      !    else
      !      print *, myid, 'now have building at ', my_idx(:, m), 'iface', buildings%iface(m), &
      !        & 'Ag', buildings%albedo(m), 'idx', d, k, i, j
      !    end if
      !  end associate
      !end do

    end subroutine

    subroutine scatter_buildings_var(comm, inpfile, prefix, varname, &
        & netcdf_idx_to_rank_sorted, send_cnts, displs, my_face_cnt, data_arr, ierr)
      integer(mpiint), intent(in) :: comm
      character(len=*), intent(in) :: inpfile, prefix, varname
      integer(iintegers), intent(in) :: netcdf_idx_to_rank_sorted(:)
      integer(mpiint), intent(in) :: send_cnts(0:), displs(0:), my_face_cnt
      real(ireals), allocatable, intent(inout) :: data_arr(:)
      integer(mpiint), intent(out) :: ierr

      integer(iintegers) :: m
      real(ireals), allocatable :: netcdf_arr(:)
      real(ireals), allocatable :: rank_sorted_arr(:)
      integer(mpiint) :: myid
      logical :: lhave_building_data

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      if (myid .eq. 0) then
        call ncload([character(len=default_str_len) :: inpfile, 'buildings', prefix, varname], netcdf_arr, ierr)
        lhave_building_data = ierr .eq. 0
      end if
      call imp_bcast(comm, lhave_building_data, 0_mpiint, ierr); call CHKERR(ierr)
      if (.not. lhave_building_data) return

      if (myid .eq. 0) then
        allocate (rank_sorted_arr(size(netcdf_arr)))
        do m = 1, size(rank_sorted_arr)
          rank_sorted_arr(netcdf_idx_to_rank_sorted(m)) = netcdf_arr(m)
        end do
        deallocate (netcdf_arr)
      else
        allocate (rank_sorted_arr(0))
      end if

      if (.not. allocated(data_arr)) then
        allocate (data_arr(my_face_cnt))
      else
        call CHKERR(size(data_arr, kind=mpiint) - my_face_cnt, 'wrong size of buildings variable'//varname)
      end if
      call MPI_Scatterv(&
        & rank_sorted_arr, &
        & send_cnts, &
        & displs, &
        & imp_ireals, &
        & data_arr, &
        & my_face_cnt, &
        & imp_ireals, &
        & 0_mpiint, &
        & comm, &
        & ierr)
      call CHKERR(ierr)
    end subroutine

  end subroutine
end module
