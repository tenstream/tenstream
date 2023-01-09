!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
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

module m_pprts

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, irealLUT, &
    & init_mpi_data_parameters, mpiint,                      &
    & imp_ireals,                                            &
    & zero, one, nan, pi,                                    &
    & nil, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi,      &
    & default_str_len

  use m_helper_functions, only: &
    & angle_between_two_normed_vec, &
    & approx, &
    & cartesian_2_spherical, &
    & CHKERR, &
    & CHKWARN, &
    & cross_3d, &
    & cstr, &
    & deallocate_allocatable, &
    & deg2rad, &
    & delta_scale, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_allreduce_max, &
    & imp_allreduce_mean, &
    & imp_allreduce_min, &
    & imp_bcast, &
    & imp_min_mean_max, &
    & inc, &
    & ind_1d_to_nd, &
    & ind_nd_to_1d, &
    & is_inrange, &
    & meanval, &
    & mpi_logical_and, &
    & ndarray_offsets, &
    & normalize_vec, &
    & rad2deg, &
    & rel_approx, &
    & rotation_matrix_world_to_local_basis, &
    & spherical_2_cartesian, &
    & toStr, &
    & triangle_area_by_vertices, &
    & vec_proj_on_plane

  use m_schwarzschild, only: schwarzschild, B_eff

  use m_optprop_parameters, only: ldebug_optprop

  use m_optprop, only: &
    & t_optprop, &
    & t_optprop_cube, &
    & t_optprop_1_2, &
    & t_optprop_3_6, &
    & t_optprop_3_10, &
    & t_optprop_3_16, &
    & t_optprop_3_24, &
    & t_optprop_3_30, &
    & t_optprop_8_10, &
    & t_optprop_8_16, &
    & t_optprop_8_18, &
    & t_optprop_3_10_ann

  use m_eddington, only: eddington_coeff_ec, eddington_coeff_zdun
  use m_geometric_coeffs, only: dir2dir3_geometric_coeffs

  use m_tenstream_options, only: read_commandline_options, luse_eddington, twostr_ratio, &
                                 options_max_solution_err, options_max_solution_time, &
                                 lcalc_nca, lskip_thermal, lschwarzschild, ltopography

  use m_petsc_helpers, only: petscGlobalVecToZero, scatterZerotoPetscGlobal, &
                             petscGlobalVecToAll, &
                             petscVecToF90, f90VecToPetsc, getVecPointer, restoreVecPointer, hegedus_trick

  use m_mcdmda, only: solve_mcdmda

  use m_pprts_base, only: &
    & allocate_pprts_solver_from_commandline, &
    & atmk, &
    & compute_gradient, &
    & destroy_pprts, &
    & destroy_solution, &
    & determine_ksp_tolerances, &
    & get_solution_uid, &
    & interpolate_cell_values_to_vertices, &
    & prepare_solution, &
    & set_dmda_cell_coordinates, &
    & setup_incSolar, &
    & setup_log_events, &
    & solver_to_str, &
    & t_atmosphere, &
    & t_coord, &
    & t_dof, &
    & t_mat_permute_info, &
    & t_solver, &
    & t_solver_1_2, &
    & t_solver_2str, &
    & t_solver_3_10, &
    & t_solver_3_16, &
    & t_solver_3_24, &
    & t_solver_3_30, &
    & t_solver_3_6, &
    & t_solver_8_10, &
    & t_solver_8_16, &
    & t_solver_8_18, &
    & t_solver_disort, &
    & t_solver_log_events, &
    & t_solver_mcdmda, &
    & t_solver_rayli, &
    & t_state_container, &
    & t_suninfo

  use m_pprts_shell, only: &
    & op_mat_getdiagonal, &
    & op_mat_mult_ediff, &
    & op_mat_mult_edir, &
    & op_mat_sor_ediff, &
    & op_mat_sor_edir, &
    & setup_matshell

  use m_pprts_explicit, only: &
    & explicit_edir, &
    & explicit_ediff

  use m_buildings, only: t_pprts_buildings, &
    & PPRTS_TOP_FACE, &
    & PPRTS_BOT_FACE, &
    & PPRTS_LEFT_FACE, &
    & PPRTS_RIGHT_FACE, &
    & PPRTS_REAR_FACE, &
    & PPRTS_FRONT_FACE

  use m_xdmf_export, only: &
    & xdmf_pprts_buildings, &
    & xdmf_pprts_srfc_flux

  use m_pprts_external_solvers, only: twostream, schwarz, pprts_rayli_wrapper, disort

  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use m_netcdfio, only: ncwrite, set_attribute, nc_var_exists

  implicit none
  private

  public :: &
    & init_pprts, &
    & set_optical_properties, set_global_optical_properties, &
    & solve_pprts, set_angles, pprts_get_result, &
    & pprts_get_result_toZero, &
    & gather_all_toZero, &
    & gather_all_to_all, &
    & scale_flx

  logical, parameter :: ldebug = .false.
  logical, parameter :: lcyclic_bc = .true.
  logical, parameter :: lprealloc = .true.

  integer(iintegers), parameter :: minimal_dimension = 3 ! this is the minimum number of gridpoints in x or y direction

contains

  !> @brief Main routine to setup PPRTS solver
  !> @details This will setup the PETSc DMDA grid and set other grid information, needed for the pprts
  !> \n Nx, Ny Nz are either global domain size or have to be local sizes if present(nxproc,nyproc)
  !> \n where nxproc and nyproc then are the number of pixel per rank for all ranks -- i.e. sum(nxproc) != Nx_global
  subroutine init_pprts(icomm, Nz, Nx, Ny, dx, dy, sundir, solver, dz1d, dz3d, nxproc, nyproc, collapseindex, solvername)
    MPI_Comm, intent(in) :: icomm         !< @param MPI_Communicator for this solver
    integer(iintegers), intent(in) :: Nz            !< @param[in] Nz     Nz is the number of layers and Nz+1 would be the number of levels
    integer(iintegers), intent(in) :: Nx            !< @param[in] Nx     number of boxes in x-direction
    integer(iintegers), intent(in) :: Ny            !< @param[in] Ny     number of boxes in y-direction
    real(ireals), intent(in) :: dx            !< @param[in] dx     physical size of grid in [m]
    real(ireals), intent(in) :: dy            !< @param[in] dy     physical size of grid in [m]
    real(ireals), intent(in) :: sundir(:)     !< @param[in] cartesian sun direction (pointing away from the sun), dim(3)

    class(t_solver), intent(inout) :: solver         !< @param[inout] solver
    real(ireals), optional, intent(in) :: dz1d(:)        !< @param[in]    dz1d    if given, dz1d is used everywhere on the rank
    real(ireals), optional, intent(in) :: dz3d(:, :, :)    !< @param[in]    dz3d    if given, dz3d has to be local domain size, cannot have global shape
    integer(iintegers), optional, intent(in) :: nxproc(:)      !< @param[in]    nxproc  if given, Nx has to be the local size, dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    integer(iintegers), optional, intent(in) :: nyproc(:)      !< @param[in]    nyproc  if given, Ny has to be the local size, dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    integer(iintegers), optional, intent(in) :: collapseindex  !< @param[in]    collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them
    character(len=*), optional, intent(in) :: solvername     !< @param[in] primarily for logging purposes, name will be prefix to logging stages

    integer(iintegers) :: k, i, j
    logical :: lpetsc_is_initialized, lview, luse_ann, lflg

    integer(mpiint) :: ierr

    if (.not. solver%linitialized) then
      call init_mpi_data_parameters(icomm)
      solver%comm = icomm
      call mpi_comm_rank(solver%comm, solver%myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(solver%comm, solver%numnodes, ierr); call CHKERR(ierr)

      call read_commandline_options(solver%comm)

      luse_ann = .false.
      call get_petsc_opt(solver%prefix, "-pprts_use_ANN", luse_ann, lflg, ierr); call CHKERR(ierr)

      solver%difftop%area_divider = 1
      solver%diffside%area_divider = 1

      solver%dirtop%area_divider = 1
      solver%dirside%area_divider = 1

      select type (solver)
      class is (t_solver_2str)
        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(0))

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = [.true.]

        allocate (solver%dirside%is_inward(0))

      class is (t_solver_disort)
        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(0))

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = [.true.]

        allocate (solver%dirside%is_inward(0))

      class is (t_solver_mcdmda)
        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(4))
        solver%diffside%is_inward = [.false., .true.]

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = [.true.]

        allocate (solver%dirside%is_inward(1))
        solver%dirside%is_inward = [.true.]

      class is (t_solver_rayli)
        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(4))
        solver%diffside%is_inward = [.false., .true., .false., .true.]

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = [.true.]

        allocate (solver%dirside%is_inward(1))
        solver%dirside%is_inward = [.true.]

      class is (t_solver_1_2)
        allocate (t_optprop_1_2 :: solver%OPP)
        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(0))

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = [.true.]

        allocate (solver%dirside%is_inward(0))

      class is (t_solver_3_6)
        allocate (t_optprop_3_6 :: solver%OPP)

        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(2))
        solver%diffside%is_inward = [.false., .true.]

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = [.true.]

        allocate (solver%dirside%is_inward(1))
        solver%dirside%is_inward = [.true.]

      class is (t_solver_3_10)
        if (luse_ann) then
          allocate (t_optprop_3_10_ann :: solver%OPP)
        else
          allocate (t_optprop_3_10 :: solver%OPP)
        end if

        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(4))
        solver%diffside%is_inward = [.false., .true., .false., .true.]

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = [.true.]

        allocate (solver%dirside%is_inward(1))
        solver%dirside%is_inward = [.true.]

      class is (t_solver_8_10)
        allocate (t_optprop_8_10 :: solver%OPP)

        allocate (solver%difftop%is_inward(2))
        solver%difftop%is_inward = [.false., .true.]

        allocate (solver%diffside%is_inward(4))
        solver%diffside%is_inward = [.false., .true., .false., .true.]

        allocate (solver%dirtop%is_inward(4))
        solver%dirtop%is_inward = .true.
        solver%dirtop%area_divider = 4

        allocate (solver%dirside%is_inward(2))
        solver%dirside%is_inward = .true.
        solver%dirside%area_divider = 2

      class is (t_solver_3_16)
        allocate (t_optprop_3_16 :: solver%OPP)

        allocate (solver%difftop%is_inward(8), source= &
                  [.false., .true., .false., .true., .false., .true., .false., .true.])

        allocate (solver%diffside%is_inward(4))
        solver%diffside%is_inward = [.false., .true., .false., .true.]

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = .true.

        allocate (solver%dirside%is_inward(1))
        solver%dirside%is_inward = .true.

      class is (t_solver_3_24)
        allocate (t_optprop_3_24 :: solver%OPP)

        allocate (solver%difftop%is_inward(8), source= &
                  [.false., .true., .false., .true., .false., .true., .false., .true.])

        allocate (solver%diffside%is_inward(8), source= &
                  [.true., .false., .true., .false., .true., .false., .true., .false.])

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = .true.

        allocate (solver%dirside%is_inward(1))
        solver%dirside%is_inward = .true.

      class is (t_solver_3_30)
        allocate (t_optprop_3_30 :: solver%OPP)

        allocate (solver%difftop%is_inward(10), source= &
                  [.false., .true., .false., .true., .false., .true., .false., .true., .false., .true.])

        allocate (solver%diffside%is_inward(10), source= &
                  [.true., .false., .true., .false., .true., .false., .true., .false., .true., .false.])

        allocate (solver%dirtop%is_inward(1))
        solver%dirtop%is_inward = .true.

        allocate (solver%dirside%is_inward(1))
        solver%dirside%is_inward = .true.

      class is (t_solver_8_16)
        allocate (t_optprop_8_16 :: solver%OPP)

        allocate (solver%difftop%is_inward(8), source= &
                  [.false., .true., .false., .true., .false., .true., .false., .true.])

        allocate (solver%diffside%is_inward(4), source=[.false., .true., .false., .true.])

        allocate (solver%dirtop%is_inward(4), source=.true.)
        solver%dirtop%area_divider = 4

        allocate (solver%dirside%is_inward(2), source=.true.)
        solver%dirside%area_divider = 2

      class is (t_solver_8_18)
        allocate (t_optprop_8_18 :: solver%OPP)

        allocate (solver%difftop%is_inward(10), source= &
                  [.false., .true., .false., .true., .false., .true., .false., .true., .false., .true.])

        allocate (solver%diffside%is_inward(4), source=[.false., .true., .false., .true.])

        allocate (solver%dirtop%is_inward(4), source=.true.)
        solver%dirtop%area_divider = 4

        allocate (solver%dirside%is_inward(2), source=.true.)
        solver%dirside%area_divider = 2

      class default
        call CHKERR(1_mpiint, 'unexpected type for solver')
      end select

      solver%difftop%dof = size(solver%difftop%is_inward)
      solver%diffside%dof = size(solver%diffside%is_inward)
      solver%dirtop%dof = size(solver%dirtop%is_inward)
      solver%dirside%dof = size(solver%dirside%is_inward)

      solver%difftop%streams = solver%difftop%dof / 2
      solver%diffside%streams = solver%diffside%dof / 2
      solver%dirtop%streams = solver%dirtop%dof
      solver%dirside%streams = solver%dirside%dof

      allocate (solver%solutions(-1000:1000))
      solver%lenable_solutions_err_estimates = &
        options_max_solution_err .gt. zero .and. options_max_solution_time .gt. zero

      if (.not. approx(dx, dy)) &
        call CHKERR(1_mpiint, 'dx and dy currently have to be the same '//toStr(dx)//' vs '//toStr(dy))

      lview = ldebug
      call get_petsc_opt(solver%prefix, "-pprts_solver_view", lview, lflg, ierr); call CHKERR(ierr)

      if (lview .and. solver%myid .eq. 0) then
        print *, 'Solver dirtop:', solver%dirtop%is_inward, ':', solver%dirtop%dof, ':', solver%dirtop%area_divider
        print *, 'Solver dirside:', solver%dirside%is_inward, ':', solver%dirside%dof, ':', solver%dirside%area_divider
        print *, 'Solver difftop:', solver%difftop%is_inward, ':', solver%difftop%dof, ':', solver%difftop%area_divider
        print *, 'Solver diffside:', solver%diffside%is_inward, ':', solver%diffside%dof, ':', solver%diffside%area_divider
      end if

      call PetscInitialized(lpetsc_is_initialized, ierr); call CHKERR(ierr)
      if (.not. lpetsc_is_initialized) call PetscInitialize(PETSC_NULL_CHARACTER, ierr); call CHKERR(ierr)
#ifdef _XLF
      call PetscPopSignalHandler(ierr); call CHKERR(ierr) ! in case of xlf ibm compilers, remove petsc signal handler -- otherwise we dont get fancy signal traps from boundschecking or FPE's
#endif

      if (present(nxproc) .and. present(nyproc)) then
        if (ldebug .and. solver%myid .eq. 0) print *, 'nxproc', shape(nxproc), '::', nxproc
        if (ldebug .and. solver%myid .eq. 0) print *, 'nyproc', shape(nyproc), '::', nyproc
        call setup_grid(solver, Nz, Nx, Ny, nxproc, nyproc, collapseindex=collapseindex)
      else
        call setup_grid(solver, Nz, max(minimal_dimension, Nx), max(minimal_dimension, Ny), collapseindex=collapseindex)
      end if

      call setup_atm()

      call print_domain_geometry_summary(solver, opt_lview=ldebug)

      ! init work vectors
      select type (solver)
      class is (t_solver_2str)
        continue
      class is (t_solver_disort)
        continue
      class is (t_solver_rayli)
        continue
      class is (t_solver_mcdmda)
        continue
      class default
        call init_memory(solver%C_dir, solver%C_diff, solver%incSolar, solver%b)
      end select

      if (present(solvername)) then
        solver%solvername = trim(solver%solvername)//trim(solvername)
      else
        solver%solvername = ''
      end if
      ! init petsc logging facilities
      call setup_log_events(solver%logs, solver%solvername)

      solver%linitialized = .true.
    else
      print *, solver%myid, 'You tried to initialize already initialized PPRTS     &
        &solver. This should not be done. If you need to reinitialize the grids, &
        &call destroy_pprts() first.'
      ierr = 1; call CHKERR(ierr)
    end if

    !Todo: this is just here so that we do not break the API. User could also
    !       directly use the set_angle routines?
    call set_angles(solver, sundir)

    ! dont need LUT, we just compute Twostream anyway
    select type (solver)
    class is (t_solver_2str)
      return
    class is (t_solver_disort)
      return
    class is (t_solver_rayli)
      return
    class is (t_solver_mcdmda)
      return
    end select

    if (.not. luse_eddington) then
      allocate (t_optprop_1_2 :: solver%OPP1d)
      call solver%OPP1d%init(solver%comm)
    end if

    call solver%OPP%init(solver%comm)
  contains
    subroutine setup_atm()
      if (.not. allocated(solver%atm)) allocate (solver%atm)

      solver%atm%dx = dx
      solver%atm%dy = dy

      if (.not. allocated(solver%atm%dz)) then
        allocate (solver%atm%dz(solver%C_one_atm%zs:solver%C_one_atm%ze, &
                                solver%C_one_atm%xs:solver%C_one_atm%xe, &
                                solver%C_one_atm%ys:solver%C_one_atm%ye), stat=ierr)
      end if

      if (present(dz1d)) then
        do j = solver%C_one_atm%ys, solver%C_one_atm%ye
          do i = solver%C_one_atm%xs, solver%C_one_atm%xe
            solver%atm%dz(:, i, j) = dz1d
          end do
        end do
      end if
      if (present(dz3d)) then
        if (any(shape(dz3d) .ne. shape(solver%atm%dz))) then
          print *, 'Whoops I got a 3D dz definition but this does not correspond to the grid definition :: shapes: ', &
            shape(dz3d), ' vs. ', shape(solver%atm%dz)
          print *, 'please know that providing a 3D dz profile has to be of local shape, it can not have global size'
          call MPI_Abort(icomm, 1 * ierr, ierr)
        end if
        solver%atm%dz = dz3d

      end if
      if (.not. present(dz1d) .and. .not. present(dz3d)) then
        print *, 'have to give either dz1d or dz3d in routine call....'
        ierr = 1; call CHKERR(ierr)
      end if

      if (.not. allocated(solver%atm%hhl)) then
        allocate (solver%atm%hhl)
        call compute_vertical_height_levels(dz=solver%atm%dz, C_hhl=solver%C_one_atm1_box, vhhl=solver%atm%hhl, &
          & prefix=solver%prefix)
      end if

      if (.not. allocated(solver%atm%hgrad)) then
        call compute_gradient(            &
          & comm=solver%comm,          &
          & atm=solver%atm,           &
          & C_hhl=solver%C_one_atm1_box,&
          & vhhl=solver%atm%hhl,       &
          & C_grad=solver%C_two1,        &
          & vgrad=solver%atm%hgrad)
      end if

      solver%sun%luse_topography = present(dz3d) .and. ltopography  ! if the user supplies 3d height levels and has set the topography option

      call determine_vertex_heights()

      call determine_1d_layers()

    end subroutine

    subroutine determine_1d_layers()
      real(ireals), pointer :: hhl(:, :, :, :) => null(), hhl1d(:) => null()
      real(ireals) :: pprts_1d_height
      real(ireals) :: N1dlayers, N1dlayers_max
      logical :: lflg

      if (.not. allocated(solver%atm%l1d)) then
        allocate (solver%atm%l1d(solver%C_one_atm%zs:solver%C_one_atm%ze))
      end if
      solver%atm%l1d = .false.

      select type (solver)
      class is (t_solver_2str)
        solver%atm%l1d = .true.
      class is (t_solver_disort)
        solver%atm%l1d = .true.
      class default
        !TODO if we have a horiz. staggered grid, this may lead to the point where one 3d box has a outgoing sideward flux but the adjacent
        !1d box does not send anything back --> therefore huge absorption :( -- would need to introduce mirror boundary conditions for
        !sideward fluxes in 1d boxes
        associate (C => solver%C_one_atm, dz => solver%atm%dz)
          solver%atm%l1d(C%ze) = any((dz(C%ze, C%xs:C%xe, C%ys:C%ye) / solver%atm%dx) .gt. twostr_ratio)
          do k = C%ze - 1, C%zs, -1
            if (any((solver%atm%dz(k, C%xs:C%xe, C%ys:C%ye) / solver%atm%dx) .gt. twostr_ratio)) then
              solver%atm%l1d(C%zs:k) = .true.
              exit
            end if
          end do
        end associate

        ! set layer to 1D if a top height of a box is higher than argument
        call get_petsc_opt(solver%prefix, "-pprts_1d_height", pprts_1d_height, lflg, ierr); call CHKERR(ierr)
        if (lflg) then
          associate ( &
              & C => solver%C_one_atm, &
              & C_hhl => solver%C_one_atm1_box, &
              & vhhl => solver%atm%hhl)

            call getVecPointer(C_hhl%da, vhhl, hhl1d, hhl, readonly=.true.)
            do k = C%ze, C%zs, -1
              if (any(hhl(i0, k, :, :) .gt. pprts_1d_height)) then
                solver%atm%l1d(C%zs:k) = .true.
                exit
              end if
            end do
            call restoreVecPointer(C_hhl%da, vhhl, hhl1d, hhl, readonly=.true.)
          end associate
        end if
      end select

      if (present(collapseindex)) then
        solver%atm%lcollapse = collapseindex .gt. i1
        solver%atm%icollapse = collapseindex
        if (solver%atm%lcollapse) then
          ierr = count(.not. solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs)))
          call CHKWARN(ierr, 'Found non 1D layers in an area that will be collapsed.'// &
            & 'This will change the results! '//&
            & 'collapse index: '//toStr(collapseindex))
          solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs)) = .true. ! if need to be collapsed, they have to be 1D.
          if (ldebug) print *, 'Using icollapse:', collapseindex, solver%atm%lcollapse
        end if
      end if

      N1dlayers = count(solver%atm%l1d)
      call imp_allreduce_max(solver%comm, N1dlayers, N1dlayers_max)
      if (N1dlayers_max .gt. 0) then
        call CHKWARN(int(N1dlayers - N1dlayers_max, mpiint), &
          & 'Nr of 1D layers does not match on all ranks'// &
          & ' while the global max(N1D_layers) is'//toStr(N1dlayers_max)// &
          & ' here (rank='//toStr(solver%myid)//') we have N1D_layers='//tostr(N1dlayers))
        solver%atm%l1d(solver%C_one_atm%zs:solver%C_one_atm%zs + N1dlayers_max - 1) = .true.
      end if

    end subroutine

    subroutine determine_vertex_heights()
      type(tVec) :: gvec_vertex_hhl
      if (allocated(solver%atm%vert_heights)) &
        & call CHKERR(1_mpiint, 'solver%vert_heights already allocated? Did not expect that could happen')

      call DMGetGlobalVector(solver%Cvert_one_atm1%da, gvec_vertex_hhl, ierr); call CHKERR(ierr)

      call interpolate_cell_values_to_vertices(&
        & solver%C_one_atm1_box, solver%atm%hhl, &
        & solver%Cvert_one_atm1, gvec_vertex_hhl)

      allocate (solver%atm%vert_heights)

      call DMCreateLocalVector(solver%Cvert_one_atm1%da, &
        & solver%atm%vert_heights, ierr); call CHKERR(ierr)

      call VecSet(solver%atm%vert_heights, zero, ierr); call CHKERR(ierr)

      call DMGlobalToLocalBegin(solver%Cvert_one_atm1%da, gvec_vertex_hhl, &
        & ADD_VALUES, solver%atm%vert_heights, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd(solver%Cvert_one_atm1%da, gvec_vertex_hhl, &
        & ADD_VALUES, solver%atm%vert_heights, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(solver%Cvert_one_atm1%da, gvec_vertex_hhl, ierr); call CHKERR(ierr)
    end subroutine
  end subroutine

  !> @brief Print a summary of the domain shape
  subroutine print_domain_geometry_summary(solver, opt_lview)
    class(t_solver), intent(in) :: solver
    logical, intent(in), optional :: opt_lview

    logical :: lview, lflg
    real(ireals), dimension(3) :: mdz, mhhl
    integer(mpiint) :: myid, numnodes, ierr
    integer(iintegers) :: k
    real(ireals), pointer :: hhl(:, :, :, :) => null(), hhl1d(:) => null()

    lview = get_arg(.false., opt_lview)
    call get_petsc_opt(solver%prefix, "-pprts_view_geometry", lview, lflg, ierr); call CHKERR(ierr)
    if (.not. lview) return

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    if (myid .eq. 0) print *, ''
    if (myid .eq. 0) print *, ' * Cell Domain Corners'//new_line('')//&
      & cstr(' rank  ', 'blue')// &
      & cstr(' zs:ze       ', 'green')// &
      & cstr(' xs:xe       ', 'blue')// &
      & cstr(' ys:ye', 'green')

    do k = 0, numnodes - 1
      if (myid .eq. k) then
        associate (C => solver%C_one)
          print *, ' '//cstr(toStr(myid), 'blue')//&
            & '     '//cstr(toStr(C%zs), 'green')//' : '//cstr(toStr(C%ze), 'green')// &
            & '     '//cstr(toStr(C%xs), 'blue')//' : '//cstr(toStr(C%xe), 'blue')// &
            & '     '//cstr(toStr(C%ys), 'green')//' : '//cstr(toStr(C%ye), 'green')
        end associate
      end if
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end do
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    if (myid .eq. 0) print *, ''

    if (myid .eq. 0) print *, '*  k(layer)'// &
      & ' '//cstr('all(1d) / any(1d)', 'blue')//' '// &
      & ' '//cstr('dz(min/mean/max)', 'red')//' [m]                   '// &
      & ' '//cstr('height [km]', 'blue')
    associate ( &
        & atm => solver%atm, &
        & C_one_atm => solver%C_one_atm, &
        & C_one_atm1 => solver%C_one_atm1)

      call getVecPointer(C_one_atm1%da, atm%hhl, hhl1d, hhl, readonly=.true.)
      do k = C_one_atm%zs, C_one_atm%ze

        call imp_min_mean_max(solver%comm, atm%dz(k, :, :), mdz)
        call imp_min_mean_max(solver%comm, (hhl(i0, k, :, :) + hhl(i0, k + 1, :, :)) / 2, mhhl)

        if (myid .eq. 0) &
          & print *, k, &
          & cstr(toStr(atm%l1d(k)), 'blue'), '        ', &
          & cstr(toStr(mdz), 'red'), ' ', cstr(toStr(mhhl * 1e-3_ireals), 'blue')
      end do
      call restoreVecPointer(C_one_atm1%da, atm%hhl, hhl1d, hhl, readonly=.true.)

      if (myid .eq. 0) print *, ' * atm dx/dy '//toStr(atm%dx)//' , '//toStr(atm%dy)
    end associate
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
  end subroutine

  !> @brief Construct PETSC grid information for regular DMDA
  !> @details setup DMDA grid containers for direct, diffuse and absorption grid
  !>  \n and fill user context containers(t_coord) which include local as well as global array sizes
  !>  \n every mpi rank has to call this
  subroutine setup_grid(solver, Nz_in, Nx, Ny, nxproc, nyproc, collapseindex)
    class(t_solver), intent(inout) :: solver
    integer(iintegers), intent(in) :: Nz_in, Nx, Ny                  !< @param[in] local number of grid boxes -- in the vertical we have Nz boxes and Nz+1 levels
    integer(iintegers), intent(in), optional :: nxproc(:), nyproc(:) !< @param[in] size of local domains on each node
    integer(iintegers), optional, intent(in) :: collapseindex        !< @param[in] collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them

    integer(iintegers) :: Nz
    DMBoundaryType, parameter :: &
      & bp = DM_BOUNDARY_PERIODIC, &
      & bn = DM_BOUNDARY_NONE,     &
      & bm = DM_BOUNDARY_MIRROR!,   &
    !& bg=DM_BOUNDARY_GHOSTED

    DMBoundaryType :: boundaries(3)
    integer(iintegers), allocatable :: nxprocp1(:), nyprocp1(:) !< last entry one larger for vertices
    integer(mpiint) :: ierr

    if (lcyclic_bc) then
      boundaries = [bn, bp, bp]
    else
      boundaries = [bn, bm, bm]
    end if

    Nz = Nz_in
    if (present(collapseindex)) then
      if (collapseindex .gt. 1) then
        Nz = Nz_in - collapseindex + i1
      end if
    end if

    if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid,&
      & 'Setting up the DMDA grid for', Nz, Nx, Ny, 'using', solver%numnodes, 'nodes'

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA C_diff', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%C_diff, Nz + 1, Nx, Ny, boundaries, &
      & solver%difftop%dof + 2 * solver%diffside%dof, nxproc, nyproc, prefix=solver%prefix)

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA C_dir', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%C_dir, Nz + 1, Nx, Ny, boundaries, &
      & solver%dirtop%dof + 2 * solver%dirside%dof, nxproc, nyproc, prefix=solver%prefix)

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA C_one', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%C_one, Nz, Nx, Ny, boundaries, &
      & i1, nxproc, nyproc, prefix=solver%prefix)
    call setup_dmda(solver%comm, solver%C_one1, Nz + 1, Nx, Ny, boundaries, &
      & i1, nxproc, nyproc, prefix=solver%prefix)

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA C_two', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%C_two1, Nz + 1, Nx, Ny, boundaries, &
      & i2, nxproc, nyproc, prefix=solver%prefix)

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA atm', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%C_one_atm, Nz_in, Nx, Ny, boundaries, &
      & i1, nxproc, nyproc, prefix=solver%prefix)

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA atm1', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%C_one_atm1, Nz_in + 1, Nx, Ny, boundaries, &
      & i1, nxproc, nyproc, prefix=solver%prefix)

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA atm1_box', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%C_one_atm1_box, Nz_in + 1, Nx, Ny, boundaries, &
      & i1, nxproc, nyproc, DMDA_STENCIL_BOX, prefix=solver%prefix)

    if (ldebug) then
      if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA srfc_one', 'green')
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end if
    call setup_dmda(solver%comm, solver%Csrfc_one, i1, Nx, Ny, boundaries, &
      & i1, nxproc, nyproc, prefix=solver%prefix)

    if (present(nxproc)) then
      allocate (nxprocp1(size(nxproc)), nyprocp1(size(nyproc)))
      nxprocp1(:) = nxproc
      nyprocp1(:) = nyproc
      nxprocp1(size(nxprocp1)) = nxprocp1(size(nxprocp1)) + 1
      nyprocp1(size(nyprocp1)) = nyprocp1(size(nyprocp1)) + 1

      if (ldebug) then
        if (solver%myid .eq. 0) print *, solver%myid, cstr('Configuring DMDA vert_one_atm1', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      end if
      call setup_dmda(solver%comm, solver%Cvert_one_atm1, Nz_in + 1, Nx + 1, Ny + 1, boundaries, &
        & i1, nxprocp1, nyprocp1, stencil_type=DMDA_STENCIL_BOX, prefix=solver%prefix)
    else
      call setup_dmda(solver%comm, solver%Cvert_one_atm1, Nz_in + 1, Nx + 1, Ny + 1, boundaries, &
        & i1, stencil_type=DMDA_STENCIL_BOX, prefix=solver%prefix)
    end if

    if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid, 'DMDA grid ready'
  end subroutine

  subroutine setup_dmda(icomm, C, Nz, Nx, Ny, boundary, dof, nxproc, nyproc, stencil_type, prefix)
    integer(mpiint), intent(in) :: icomm
    type(t_coord), allocatable :: C
    integer(iintegers), intent(in) :: Nz, Nx, Ny, dof
    DMBoundaryType, intent(in) :: boundary(:)
    integer(iintegers), intent(in), optional :: nxproc(:), nyproc(:) ! size of local domains on each node
    DMDAStencilType, intent(in), optional :: stencil_type
    character(len=*), intent(in), optional :: prefix
    DMDAStencilType :: opt_stencil_type

    integer(iintegers), parameter :: stencil_size = 1
    integer(mpiint) :: myid, ierr
    call mpi_comm_rank(icomm, myid, ierr); call CHKERR(ierr)

    if (any([present(nxproc), present(nyproc)])&
      & .and. .not. all([present(nxproc), present(nyproc)])) &
      & call CHKERR(1_mpiint, 'have to have both, nxproc AND nyproc present or none')

    if (ldebug .and. present(nxproc) .and. present(nyproc)) &
      & print *, myid, 'setup_dmda nxproc', nxproc, 'nyproc', nyproc

    opt_stencil_type = get_arg(DMDA_STENCIL_STAR, stencil_type)

    allocate (C)

    C%comm = icomm

    C%dof = i1 * dof
    if (present(nxproc) .and. present(nyproc)) then
      call DMDACreate3d(C%comm, &
                        boundary(1), boundary(2), boundary(3), &
                        opt_stencil_type, &
                        Nz, sum(nxproc), sum(nyproc), &
                        i1, size(nxproc, kind=iintegers), size(nyproc, kind=iintegers), &
                        C%dof, stencil_size, &
                        [Nz], nxproc, nyproc, &
                        C%da, ierr)
      call CHKERR(ierr)
    else
      call DMDACreate3d(C%comm, &
                        boundary(1), boundary(2), boundary(3), &
                        opt_stencil_type, &
                        i1 * Nz, Nx, Ny, &
                        i1, PETSC_DECIDE, PETSC_DECIDE, &
                        C%dof, stencil_size, &
                        [Nz], PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
                        C%da, ierr); call CHKERR(ierr)
    end if

    if (present(prefix)) then
      call PetscObjectSetOptionsPrefix(C%da, prefix, ierr); call CHKERR(ierr)
    end if
    ! need this first setfromoptions call because of a bug which happens with intel compilers
    call DMSetFromOptions(C%da, ierr); call CHKERR(ierr)

    call DMSetMatType(C%da, MATAIJ, ierr); call CHKERR(ierr)
    if (lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE, ierr); call CHKERR(ierr)
    call DMSetFromOptions(C%da, ierr); call CHKERR(ierr)
    call DMSetup(C%da, ierr); call CHKERR(ierr)

    if (ldebug) call DMView(C%da, PETSC_VIEWER_STDOUT_WORLD, ierr); call CHKERR(ierr)
    call setup_coords(C)
  contains
    subroutine setup_coords(C)
      type(t_coord) :: C
      DMBoundaryType :: bx, by, bz
      DMDAStencilType :: st
      integer(iintegers) :: stencil_width, nproc_x, nproc_y, nproc_z, Ndof
      integer(mpiint) :: numnodes

      call DMDAGetInfo(C%da, C%dim, &
                       C%glob_zm, C%glob_xm, C%glob_ym, &
                       nproc_z, nproc_x, nproc_y, Ndof, stencil_width, &
                       bx, by, bz, st, ierr); call CHKERR(ierr)

      call DMDAGetCorners(C%da, C%zs, C%xs, C%ys, C%zm, C%xm, C%ym, ierr); call CHKERR(ierr)
      C%xe = C%xs + C%xm - 1
      C%ye = C%ys + C%ym - 1
      C%ze = C%zs + C%zm - 1

      call DMDAGetGhostCorners(C%da, C%gzs, C%gxs, C%gys, C%gzm, C%gxm, C%gym, ierr); call CHKERR(ierr)
      C%gxe = C%gxs + C%gxm - 1
      C%gye = C%gys + C%gym - 1
      C%gze = C%gzs + C%gzm - 1

      allocate (C%neighbors(0:3**C%dim - 1))
      call DMDAGetNeighbors(C%da, C%neighbors, ierr); call CHKERR(ierr)
      call mpi_comm_size(icomm, numnodes, ierr); call CHKERR(ierr)
      if (numnodes .gt. i1) then
        if (ldebug .and. (C%dim .eq. 3)) print *, 'PETSC id', myid, C%dim, 'Neighbors are', C%neighbors(10), &
          C%neighbors(4), &
          C%neighbors(16), &
          C%neighbors(22), &
          'while I am ', C%neighbors(13)

        if (ldebug .and. (C%dim .eq. 2)) print *, 'PETSC id', myid, C%dim, 'Neighbors are', C%neighbors(1), &
          C%neighbors(3), &
          C%neighbors(7), &
          C%neighbors(5), &
          'while I am ', C%neighbors(4)
      end if
      if (C%glob_xm .lt. i2) call CHKWARN(1_mpiint, 'Global domain is too small in x-direction (Nx='//toStr(C%glob_xm)// &
                                          '). Iterative solvers may not converge. You may want to stick to 1D solvers anyway.')
      if (C%glob_ym .lt. i2) call CHKWARN(1_mpiint, 'Global domain is too small in y-direction (Ny='//toStr(C%glob_ym)// &
                                          '). Iterative solvers may not converge. You may want to stick to 1D solvers anyway.')

      if (C%xm .lt. i1) call CHKERR(1_mpiint, 'Local domain is too small in x-direction (Nx='//toStr(C%xm)// &
                                    '). However, need at least 1')
      if (C%ym .lt. i1) call CHKERR(1_mpiint, 'Local domain is too small in y-direction (Ny='//toStr(C%ym)// &
                                    '). However, need at least 1')
    end subroutine
  end subroutine

  !> @brief Determine height levels by summing up the atm%dz with the assumption that TOA is at a constant value
  !>        or a max_height is given in the option database
  subroutine compute_vertical_height_levels(dz, C_hhl, vhhl, prefix)
    type(t_coord), intent(in) :: C_hhl
    real(ireals), intent(in) :: dz(:, :, :)
    type(tVec), intent(inout) :: vhhl
    character(len=*), intent(in) :: prefix

    type(tVec) :: g_hhl
    real(ireals), pointer :: hhl(:, :, :, :) => null(), hhl1d(:) => null()
    real(ireals) :: max_height, global_max_height

    integer(mpiint) :: comm, ierr
    integer(iintegers) :: k, i, j
    integer(iintegers) :: kk, ii, jj
    logical :: lflg

    max_height = zero

    call DMGetGlobalVector(C_hhl%da, g_hhl, ierr); call CHKERR(ierr)
    call VecSet(g_hhl, zero, ierr); call CHKERR(ierr)

    call getVecPointer(C_hhl%da, g_hhl, hhl1d, hhl)

    do j = C_hhl%ys, C_hhl%ye
      jj = i1 + j - C_hhl%ys
      do i = C_hhl%xs, C_hhl%xe
        ii = i1 + i - C_hhl%xs

        hhl(i0, C_hhl%zs, i, j) = zero
        do k = C_hhl%zs, C_hhl%ze - 1
          kk = i1 + k - C_hhl%zs
          hhl(i0, k + 1, i, j) = hhl(i0, k, i, j) - dz(kk, ii, jj)
          max_height = min(max_height, hhl(i0, k + 1, i, j))
        end do
      end do
    end do

    call get_petsc_opt(prefix, "-pprts_global_max_height", global_max_height, lflg, ierr); call CHKERR(ierr)
    if (.not. lflg) then
      call PetscObjectGetComm(C_hhl%da, comm, ierr); call CHKERR(ierr)
      call imp_allreduce_min(comm, max_height, global_max_height)
    end if
    hhl(i0, :, :, :) = hhl(i0, :, :, :) - global_max_height

    call restoreVecPointer(C_hhl%da, g_hhl, hhl1d, hhl)

    call DMCreateLocalVector(C_hhl%da, vhhl, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(C_hhl%da, g_hhl, INSERT_VALUES, vhhl, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd(C_hhl%da, g_hhl, INSERT_VALUES, vhhl, ierr); call CHKERR(ierr)
    call DMRestoreGlobalVector(C_hhl%da, g_hhl, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(vhhl, C_hhl%da, "-pprts_show_hhl", ierr); call CHKERR(ierr)
  end subroutine

  !> @brief initialize basic memory structs like PETSc vectors and matrices
  subroutine init_memory(C_dir, C_diff, incSolar, b)
    type(t_coord), intent(in) :: C_dir, C_diff
    type(tVec), intent(inout), allocatable :: b, incSolar
    integer(mpiint) :: ierr

    if (.not. allocated(incSolar)) allocate (incSolar)
    if (.not. allocated(b)) allocate (b)

    call DMCreateGlobalVector(C_dir%da, incSolar, ierr); call CHKERR(ierr)
    call DMCreateGlobalVector(C_diff%da, b, ierr); call CHKERR(ierr)

    call VecSet(incSolar, zero, ierr); call CHKERR(ierr)
    call VecSet(b, zero, ierr); call CHKERR(ierr)
  end subroutine

  subroutine set_angles(solver, sundir)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: sundir(:)      !< @param[in] cartesian sun direction
    integer(mpiint) :: ierr

    if (.not. solver%linitialized) then
      print *, solver%myid, 'You tried to set angles in the PPRTS solver.  &
        & This should be called right after init_pprts'
      ierr = 1; call CHKERR(ierr)
    end if

    call setup_suninfo(solver, sundir, solver%sun)
  end subroutine

  !> @brief set direction where sun stands
  !> @details save sun azimuth and zenith angle
  !>   \n sun azimuth is reduced to the range of [0,90] and the transmission of direct radiation is contributed for by a integer increment,
  !>   \n determining which neighbouring box is used in the horizontal direction
  subroutine setup_suninfo(solver, sundir, sun)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(in) :: sundir(:)
    type(t_suninfo), intent(inout) :: sun

    logical :: lflg
    integer(mpiint) :: ierr

    if (.not. allocated(solver%atm%hgrad)) call CHKERR(1_mpiint, 'atm%hgrad not initialized!')

    call cartesian_2_spherical(sundir, sun%phi, sun%theta, ierr)

    call get_petsc_opt(solver%prefix, "-pprts_force_zenith", sun%theta, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(solver%prefix, "-pprts_force_azimuth", sun%phi, lflg, ierr); call CHKERR(ierr)

    sun%sundir = spherical_2_cartesian(sun%phi, sun%theta)

    if (sun%theta .ge. 90._ireals) sun%theta = -one

    sun%costheta = max(cos(deg2rad(sun%theta)), zero)
    sun%sintheta = max(sin(deg2rad(sun%theta)), zero)

    ! use symmetry for direct beam: always use azimuth [0,90] an just reverse the order where we insert the coeffs
    sun%symmetry_phi = sym_rot_phi(sun%phi)

    if (sin(deg2rad(sun%phi)) .gt. zero) then! phi between 0 and 180 degreee
      sun%xinc = i0
    else
      sun%xinc = i1
    end if

    if (cos(deg2rad(sun%phi)) .lt. zero) then ! phi between 90 and 270 degree
      sun%yinc = i1
    else
      sun%yinc = i0
    end if

    call print_suninfo_summary(solver, ldebug)

  contains
    pure elemental function sym_rot_phi(phi)
      real(ireals) :: sym_rot_phi
      real(ireals), intent(in) :: phi
      ! ''swap'' phi axis down to the range of [0,180]
      sym_rot_phi = acos(cos(deg2rad(phi)))
      !print *,'1st phi swap',phi,' :: ',sym_rot_phi,'=',phi*pi/180,cos(phi*pi/180),acos(cos(phi*pi/180))
      ! and then mirror it onto range [0,90]
      sym_rot_phi = min(90._ireals, max(0._ireals, rad2deg(asin(sin(sym_rot_phi)))))
      !print *,'2nd phi swap',phi,' :: ',sym_rot_phi,'=',&
      ! sin(sym_rot_phi),asin(sin(sym_rot_phi)),asin(sin(sym_rot_phi)) /pi * 180,int(asin(sin(sym_rot_phi)) /pi * 180)
    end function
  end subroutine

  !> @brief Print a summary of the sun_info
  subroutine print_suninfo_summary(solver, opt_lview)
    class(t_solver), intent(in) :: solver
    logical, intent(in), optional :: opt_lview

    logical :: lview, lflg
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: k
    real(ireals), dimension(3) :: msundir

    lview = get_arg(.false., opt_lview)
    call get_petsc_opt(solver%prefix, "-pprts_view_suninfo", lview, lflg, ierr); call CHKERR(ierr)
    if (.not. lview) return

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
    associate (sun => solver%sun, C_one => solver%C_one)

      do k = 1, 3
        call imp_min_mean_max(solver%comm, sun%sundir(k), msundir)
        if (myid .eq. 0) then
          if (.not. all(approx(sun%sundir(k), msundir))) then
            print *, 'Sundir rank=0 :: ', sun%sundir
            ierr = int(k, mpiint)
            call CHKERR(ierr, 'sundir component '//toStr(k)//' does not match across ranks '//&
              & 'min mean max sundir('//toStr(k)//') = '//toStr(msundir))
          end if
        end if
      end do

      if (solver%myid .eq. 0) &
        & print *, ' * '//cstr('sundir ['//toStr(sun%sundir)//'] ', 'red'), &
        & ' xinc', sun%xinc, 'yinc', sun%yinc

      if (myid .eq. 0) print *, ' * '// &
        & ' '//cstr('theta '//toStr(sun%theta), 'red')//'                    '// &
        & ' '//cstr('phi '//toStr(sun%phi), 'blue')

    end associate
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
  end subroutine

  !> @brief setup topography information
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information and
  !>   \n tweak the local sun angles to bend the rays.

  !> @brief create PETSc matrix and inserts diagonal elements
  !> @details one important step for performance is to set preallocation of matrix structure.
  !>  \n  i.e. determine the number of local vs. remote number of entries in each row.
  !>  \n  DMDA actually provides a preallocation but this assumes that degrees of freedom on neighbouring boxes are fully connected to local ones.
  !>  \n  this does of course drastically overestimate non-zeros as we need only the streams that actually send radiation in the respective direction.
  !>  \n  at the moment preallocation routines determine nonzeros by manually checking bounadries --
  !>  \n  !todo we should really use some form of iterating through the entries as it is done in the matrix assembly routines and just flag the rows
  subroutine init_Matrix(solver, C, A, prealloc_subroutine)
    interface
      subroutine preallocation_sub(solver, C, d_nnz, o_nnz)
        import :: t_solver, t_coord, iintegers
        class(t_solver), intent(in) :: solver
        type(t_coord), intent(in) :: C
        integer(iintegers), allocatable :: d_nnz(:)
        integer(iintegers), allocatable :: o_nnz(:)
      end subroutine
    end interface

    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    type(tMat), allocatable, intent(inout) :: A
    procedure(preallocation_sub), optional :: prealloc_subroutine

    integer(iintegers), dimension(:), allocatable :: o_nnz, d_nnz
    integer(mpiint) :: numnodes
    integer(mpiint) :: ierr

    if (allocated(A)) return

    allocate (A)
    call DMCreateMatrix(C%da, A, ierr); call CHKERR(ierr)
    call PetscObjectSetOptionsPrefix(A, trim(solver%prefix), ierr); call CHKERR(ierr)

    call mpi_comm_size(C%comm, numnodes, ierr); call CHKERR(ierr)

    if (numnodes .gt. 1) then
      if (lprealloc .and. present(prealloc_subroutine)) then
        call prealloc_subroutine(solver, C, d_nnz, o_nnz)
        call MatMPIAIJSetPreallocation(A, C%dof + 1, d_nnz, C%dof, o_nnz, ierr); call CHKERR(ierr)
      else ! poor mans perallocation uses way more memory...
        call CHKERR(1_mpiint, 'init_Matrix::setPreallocation : poor mans preallocation should really not be used...')
        call MatMPIAIJSetPreallocation(A, C%dof + 1, PETSC_NULL_INTEGER, C%dof, PETSC_NULL_INTEGER, ierr); call CHKERR(ierr)
      end if
    else
      call MatSeqAIJSetPreallocation(A, C%dof + i1, PETSC_NULL_INTEGER, ierr); call CHKERR(ierr)
    end if

    ! If matrix is resetted, keep nonzero pattern and allow to non-zero allocations -- those should not be many
    call MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE, ierr); call CHKERR(ierr)

    ! pressure mesh  may wiggle a bit and change atm%l1d -- keep the nonzeros flexible
    call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); call CHKERR(ierr)

    ! call MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr) ;call CHKERR(ierr) ! dont throw away the zero -- this completely destroys preallocation performance

    call MatSetUp(A, ierr); call CHKERR(ierr)
  end subroutine

  subroutine mat_set_diagonal(A, vdiag)
    type(tMat) :: A
    real(ireals), intent(in), optional :: vdiag

    integer(iintegers) :: is, ie, irow
    real(ireals) :: v
    integer(mpiint) :: ierr

    v = get_arg(one, vdiag)

    call MatGetOwnershipRange(A, is, ie, ierr); call CHKERR(ierr)
    do irow = is, ie - 1
      call MatSetValue(A, irow, irow, v, INSERT_VALUES, ierr); call CHKERR(ierr)
    end do
  end subroutine

  subroutine setup_direct_preallocation(solver, C, d_nnz, o_nnz)
    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    integer(iintegers), allocatable :: d_nnz(:)
    integer(iintegers), allocatable :: o_nnz(:)
    type(tVec) :: v_o_nnz, v_d_nnz
    type(tVec) :: g_o_nnz, g_d_nnz
    real(ireals), pointer :: xo(:, :, :, :) => null(), xd(:, :, :, :) => null()
    real(ireals), pointer :: xo1d(:) => null(), xd1d(:) => null()

    integer(iintegers) :: vsize, i, j, k, isrc, idst, xinc, yinc, src, dst, icnt

    logical :: llocal_src, llocal_dst
    integer(mpiint) :: myid, ierr
    MatStencil :: row(4, C%dof), col(4, C%dof)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, myid, 'building direct o_nnz for mat with', C%dof, 'dof'
    call DMGetLocalVector(C%da, v_o_nnz, ierr); call CHKERR(ierr)
    call DMGetLocalVector(C%da, v_d_nnz, ierr); call CHKERR(ierr)

    call getVecPointer(C%da, v_o_nnz, xo1d, xo)
    call getVecPointer(C%da, v_d_nnz, xd1d, xd)

    xd = zero
    xo = zero

    icnt = -1
    do j = C%ys, C%ye
      do i = C%xs, C%xe
        do k = C%zs, C%ze - 1

          if (solver%atm%l1d(atmk(solver%atm, k))) then
            do idst = i0, solver%dirtop%dof - 1
              call inc(xd(idst, k + 1, i, j), one)
            end do
          else
            xinc = solver%sun%xinc
            yinc = solver%sun%yinc

            do idst = 1, solver%dirtop%dof
              dst = idst
              row(MatStencil_j, dst) = i
              row(MatStencil_k, dst) = j
              row(MatStencil_i, dst) = k + 1
              row(MatStencil_c, dst) = dst - i1 ! Define transmission towards the lower/upper lid
            end do

            do idst = 1, solver%dirside%dof
              dst = idst + solver%dirtop%dof
              row(MatStencil_j, dst) = i + xinc
              row(MatStencil_k, dst) = j
              row(MatStencil_i, dst) = k
              row(MatStencil_c, dst) = dst - i1 ! Define transmission towards the left/right lid
            end do

            do idst = 1, solver%dirside%dof
              dst = idst + solver%dirtop%dof + solver%dirside%dof
              row(MatStencil_j, dst) = i
              row(MatStencil_k, dst) = j + yinc
              row(MatStencil_i, dst) = k
              row(MatStencil_c, dst) = dst - i1 ! Define transmission towards the front/back lid
            end do

            do isrc = 1, solver%dirtop%dof
              src = isrc
              col(MatStencil_j, src) = i
              col(MatStencil_k, src) = j
              col(MatStencil_i, src) = k
              col(MatStencil_c, src) = src - i1 ! Define transmission towards the lower/upper lid
            end do

            do isrc = 1, solver%dirside%dof
              src = isrc + solver%dirtop%dof
              col(MatStencil_j, src) = i + 1 - xinc
              col(MatStencil_k, src) = j
              col(MatStencil_i, src) = k
              col(MatStencil_c, src) = src - i1 ! Define transmission towards the left/right lid
            end do

            do isrc = 1, solver%dirside%dof
              src = isrc + solver%dirtop%dof + solver%dirside%dof
              col(MatStencil_j, src) = i
              col(MatStencil_k, src) = j + 1 - yinc
              col(MatStencil_i, src) = k
              col(MatStencil_c, src) = src - i1 ! Define transmission towards the front/back lid
            end do

            do idst = 1, C%dof
              icnt = icnt + 1
              llocal_dst = .true.
              if (C%neighbors(10) .ne. myid &
                & .and. C%neighbors(10) .ge. i0 &
                & .and. row(MatStencil_j, idst) .lt. C%xs) &
                & llocal_dst = .false. ! have real neighbor west  and is not local entry
              if (C%neighbors(16) .ne. myid &
                & .and. C%neighbors(16) .ge. i0 &
                & .and. row(MatStencil_j, idst) .gt. C%xe) &
                & llocal_dst = .false. ! have real neighbor east  and is not local entry
              if (C%neighbors(4) .ne. myid &
                & .and. C%neighbors(4) .ge. i0 &
                & .and. row(MatStencil_k, idst) .lt. C%ys) &
                & llocal_dst = .false. ! have real neighbor south and is not local entry
              if (C%neighbors(22) .ne. myid &
                & .and. C%neighbors(22) .ge. i0 &
                & .and. row(MatStencil_k, idst) .gt. C%ye) &
                & llocal_dst = .false. ! have real neighbor north and is not local entry

              do isrc = 1, C%dof
                llocal_src = .true.

                if (C%neighbors(10) .ne. myid &
                  & .and. C%neighbors(10) .ge. i0 &
                  & .and. col(MatStencil_j, isrc) .lt. C%xs) &
                  & llocal_src = .false. ! have real neighbor west  and is not local entry
                if (C%neighbors(16) .ne. myid &
                  & .and. C%neighbors(16) .ge. i0 &
                  & .and. col(MatStencil_j, isrc) .gt. C%xe) &
                  & llocal_src = .false. ! have real neighbor east  and is not local entry
                if (C%neighbors(4) .ne. myid &
                  & .and. C%neighbors(4) .ge. i0 &
                  & .and. col(MatStencil_k, isrc) .lt. C%ys) &
                  & llocal_src = .false. ! have real neighbor south and is not local entry
                if (C%neighbors(22) .ne. myid &
                  & .and. C%neighbors(22) .ge. i0 &
                  & .and. col(MatStencil_k, isrc) .gt. C%ye) &
                  & llocal_src = .false. ! have real neighbor north and is not local entry

                !if(myid.eq.0) print *,myid,icnt,k,i,j,'::',idst,isrc,'::',llocal_dst,llocal_src
                if (llocal_dst .and. llocal_src) then
                  call inc(xd(row(4, idst), row(3, idst), row(2, idst), row(1, idst)), one)
                else
                  call inc(xo(row(4, idst), row(3, idst), row(2, idst), row(1, idst)), one)
                end if
              end do
            end do
          end if
        end do
      end do
    end do

    call restoreVecPointer(C%da, v_o_nnz, xo1d, xo)
    call restoreVecPointer(C%da, v_d_nnz, xd1d, xd)

    call DMGetGlobalVector(C%da, g_o_nnz, ierr); call CHKERR(ierr)
    call DMGetGlobalVector(C%da, g_d_nnz, ierr); call CHKERR(ierr)
    call VecSet(g_o_nnz, zero, ierr); call CHKERR(ierr)
    call VecSet(g_d_nnz, zero, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da, v_o_nnz, ADD_VALUES, g_o_nnz, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd(C%da, v_o_nnz, ADD_VALUES, g_o_nnz, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da, v_d_nnz, ADD_VALUES, g_d_nnz, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd(C%da, v_d_nnz, ADD_VALUES, g_d_nnz, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(C%da, v_o_nnz, ierr); call CHKERR(ierr)
    call DMRestoreLocalVector(C%da, v_d_nnz, ierr); call CHKERR(ierr)

    call getVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.true.)
    call getVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.true.)

    call VecGetLocalSize(g_d_nnz, vsize, ierr); call CHKERR(ierr)
    allocate (o_nnz(0:vsize - 1))
    allocate (d_nnz(0:vsize - 1))

    o_nnz = int(xo1d, kind=iintegers)
    d_nnz = int(xd1d, kind=iintegers) + i1  ! +1 for diagonal entries

    call restoreVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.true.)
    call restoreVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.true.)

    call DMRestoreGlobalVector(C%da, g_o_nnz, ierr); call CHKERR(ierr)
    call DMRestoreGlobalVector(C%da, g_d_nnz, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, myid, 'direct d_nnz, ', sum(d_nnz), 'o_nnz', sum(o_nnz), &
      'together:', sum(d_nnz) + sum(o_nnz), 'expected less than', vsize * (C%dof + 1)
  end subroutine

  subroutine setup_diffuse_preallocation(solver, C, d_nnz, o_nnz)
    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    integer(iintegers), allocatable :: d_nnz(:)
    integer(iintegers), allocatable :: o_nnz(:)
    type(tVec) :: v_o_nnz, v_d_nnz
    type(tVec) :: g_o_nnz, g_d_nnz
    real(ireals), pointer :: xo(:, :, :, :) => null(), xd(:, :, :, :) => null()
    real(ireals), pointer :: xo1d(:) => null(), xd1d(:) => null()

    integer(iintegers) :: vsize, i, j, k, isrc, idst, src, dst, icnt
    integer(iintegers) :: dst_id, src_id, idof
    integer(mpiint) :: myid, ierr
    MatStencil :: row(4, 0:C%dof - 1), col(4, 0:C%dof - 1)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, myid, 'building diffuse o_nnz for mat with', C%dof, 'dof'
    call DMGetLocalVector(C%da, v_o_nnz, ierr); call CHKERR(ierr)
    call DMGetLocalVector(C%da, v_d_nnz, ierr); call CHKERR(ierr)

    call getVecPointer(C%da, v_o_nnz, xo1d, xo)
    call getVecPointer(C%da, v_d_nnz, xd1d, xd)

    xd = zero
    xo = zero

    icnt = -1
    do j = C%ys, C%ye
      do i = C%xs, C%xe
        do k = C%zs, C%ze - 1

          if (solver%atm%l1d(atmk(solver%atm, k))) then
            do idof = 1, solver%difftop%dof
              if (solver%difftop%is_inward(idof)) then
                call inc(xd(idof - 1, k + 1, i, j), real(solver%difftop%dof, ireals))
              else
                call inc(xd(idof - 1, k, i, j), real(solver%difftop%dof, ireals))
              end if
            end do
          else

            ! Diffuse Coefficients Code Begin
            do idof = 1, solver%difftop%dof
              src = idof - 1
              if (solver%difftop%is_inward(idof)) then
                col(MatStencil_j, src) = i
                col(MatStencil_k, src) = j
                col(MatStencil_i, src) = k
                col(MatStencil_c, src) = src
              else
                col(MatStencil_j, src) = i
                col(MatStencil_k, src) = j
                col(MatStencil_i, src) = k + 1
                col(MatStencil_c, src) = src
              end if
            end do

            do idof = 1, solver%diffside%dof
              src = solver%difftop%dof + idof - 1
              if (solver%diffside%is_inward(idof)) then
                col(MatStencil_j, src) = i
                col(MatStencil_k, src) = j
                col(MatStencil_i, src) = k
                col(MatStencil_c, src) = src
              else
                col(MatStencil_j, src) = i + 1
                col(MatStencil_k, src) = j
                col(MatStencil_i, src) = k
                col(MatStencil_c, src) = src
              end if
            end do

            do idof = 1, solver%diffside%dof
              src = solver%difftop%dof + solver%diffside%dof + idof - 1
              if (solver%diffside%is_inward(idof)) then
                col(MatStencil_j, src) = i
                col(MatStencil_k, src) = j
                col(MatStencil_i, src) = k
                col(MatStencil_c, src) = src
              else
                col(MatStencil_j, src) = i
                col(MatStencil_k, src) = j + 1
                col(MatStencil_i, src) = k
                col(MatStencil_c, src) = src
              end if
            end do

            do idof = 1, solver%difftop%dof
              dst = idof - 1
              if (solver%difftop%is_inward(idof)) then
                row(MatStencil_j, dst) = i
                row(MatStencil_k, dst) = j
                row(MatStencil_i, dst) = k + 1
                row(MatStencil_c, dst) = dst
              else
                row(MatStencil_j, dst) = i
                row(MatStencil_k, dst) = j
                row(MatStencil_i, dst) = k
                row(MatStencil_c, dst) = dst
              end if
            end do

            do idof = 1, solver%diffside%dof
              dst = solver%difftop%dof + idof - 1
              if (solver%diffside%is_inward(idof)) then
                row(MatStencil_j, dst) = i + 1
                row(MatStencil_k, dst) = j
                row(MatStencil_i, dst) = k
                row(MatStencil_c, dst) = dst
              else
                row(MatStencil_j, dst) = i
                row(MatStencil_k, dst) = j
                row(MatStencil_i, dst) = k
                row(MatStencil_c, dst) = dst
              end if
            end do

            do idof = 1, solver%diffside%dof
              dst = solver%difftop%dof + solver%diffside%dof + idof - 1
              if (solver%diffside%is_inward(idof)) then
                row(MatStencil_j, dst) = i
                row(MatStencil_k, dst) = j + 1
                row(MatStencil_i, dst) = k
                row(MatStencil_c, dst) = dst
              else
                row(MatStencil_j, dst) = i
                row(MatStencil_k, dst) = j
                row(MatStencil_i, dst) = k
                row(MatStencil_c, dst) = dst
              end if
            end do
            ! Diffuse Coefficients Code End

            do idst = 0, C%dof - 1
              icnt = icnt + 1
              dst_id = myid
              if (C%neighbors(10) .ne. myid &
                & .and. C%neighbors(10) .ge. i0 &
                & .and. row(MatStencil_j, idst) .lt. C%xs) &
                & dst_id = C%neighbors(10) ! have real neighbor west  and is not local entry
              if (C%neighbors(16) .ne. myid &
                & .and. C%neighbors(16) .ge. i0 &
                & .and. row(MatStencil_j, idst) .gt. C%xe) &
                & dst_id = C%neighbors(16) ! have real neighbor east  and is not local entry
              if (C%neighbors(4) .ne. myid &
                & .and. C%neighbors(4) .ge. i0 &
                & .and. row(MatStencil_k, idst) .lt. C%ys) &
                & dst_id = C%neighbors(4) ! have real neighbor south and is not local entry
              if (C%neighbors(22) .ne. myid &
                & .and. C%neighbors(22) .ge. i0 &
                & .and. row(MatStencil_k, idst) .gt. C%ye) &
                & dst_id = C%neighbors(22) ! have real neighbor north and is not local entry

              do isrc = 0, C%dof - 1
                src_id = myid

                if (C%neighbors(10) .ne. myid &
                  & .and. C%neighbors(10) .ge. i0 &
                  & .and. col(MatStencil_j, isrc) .lt. C%xs) &
                  & src_id = C%neighbors(10) ! have real neighbor west  and is not local entry
                if (C%neighbors(16) .ne. myid &
                  & .and. C%neighbors(16) .ge. i0 &
                  & .and. col(MatStencil_j, isrc) .gt. C%xe) &
                  & src_id = C%neighbors(16) ! have real neighbor east  and is not local entry
                if (C%neighbors(4) .ne. myid &
                  & .and. C%neighbors(4) .ge. i0 &
                  & .and. col(MatStencil_k, isrc) .lt. C%ys) &
                  & src_id = C%neighbors(4) ! have real neighbor south and is not local entry
                if (C%neighbors(22) .ne. myid &
                  & .and. C%neighbors(22) .ge. i0 &
                  & .and. col(MatStencil_k, isrc) .gt. C%ye) &
                  & src_id = C%neighbors(22) ! have real neighbor north and is not local entry

                if (src_id .eq. dst_id) then
                  call inc(xd(row(4, idst), row(3, idst), row(2, idst), row(1, idst)), one)
                else
                  call inc(xo(row(4, idst), row(3, idst), row(2, idst), row(1, idst)), one)
                end if

              end do
            end do

          end if ! atm_1d / 3d?
        end do ! k

        ! Surface entries E_up
        do idof = 1, solver%difftop%dof
          src = idof - 1
          if (.not. solver%difftop%is_inward(idof)) then
            call inc(xd(src, C%ze, i, j), real(solver%difftop%streams, ireals))
          end if
        end do

      end do
    end do

    call restoreVecPointer(C%da, v_o_nnz, xo1d, xo)
    call restoreVecPointer(C%da, v_d_nnz, xd1d, xd)

    call DMGetGlobalVector(C%da, g_o_nnz, ierr); call CHKERR(ierr)
    call DMGetGlobalVector(C%da, g_d_nnz, ierr); call CHKERR(ierr)
    call VecSet(g_o_nnz, zero, ierr); call CHKERR(ierr)
    call VecSet(g_d_nnz, zero, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da, v_o_nnz, ADD_VALUES, g_o_nnz, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd(C%da, v_o_nnz, ADD_VALUES, g_o_nnz, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da, v_d_nnz, ADD_VALUES, g_d_nnz, ierr); call CHKERR(ierr)
    call DMLocalToGlobalEnd(C%da, v_d_nnz, ADD_VALUES, g_d_nnz, ierr); call CHKERR(ierr)

    call DMRestoreLocalVector(C%da, v_o_nnz, ierr); call CHKERR(ierr)
    call DMRestoreLocalVector(C%da, v_d_nnz, ierr); call CHKERR(ierr)

    call getVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.true.)
    call getVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.true.)

    call VecGetLocalSize(g_d_nnz, vsize, ierr); call CHKERR(ierr)
    allocate (o_nnz(0:vsize - 1))
    allocate (d_nnz(0:vsize - 1))

    o_nnz = int(xo1d, kind=iintegers)
    d_nnz = int(xd1d, kind=iintegers) + i1  ! +1 for diagonal entries

    call restoreVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.true.)
    call restoreVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.true.)

    call DMRestoreGlobalVector(C%da, g_o_nnz, ierr); call CHKERR(ierr)
    call DMRestoreGlobalVector(C%da, g_d_nnz, ierr); call CHKERR(ierr)

    if (myid .eq. 0 .and. ldebug) print *, myid, 'diffuse d_nnz, ', sum(d_nnz), 'o_nnz', sum(o_nnz), &
      'together:', sum(d_nnz) + sum(o_nnz), 'expected less than', vsize * (C%dof + 1)

  end subroutine

  subroutine mat_info(comm, A)
    MPI_Comm, intent(in) :: comm
    type(tMat) :: A
    double precision :: info(MAT_INFO_SIZE)
    double precision :: mal, nz_allocated, nz_used, nz_unneeded
    integer(iintegers) :: m, n
    integer(mpiint) :: myid, ierr

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call MatGetInfo(A, MAT_LOCAL, info, ierr); call CHKWARN(ierr)
    mal = info(MAT_INFO_MALLOCS)
    nz_allocated = info(MAT_INFO_NZ_ALLOCATED)
    nz_used = info(MAT_INFO_NZ_USED)
    nz_unneeded = info(MAT_INFO_NZ_UNNEEDED)

    call MatGetOwnershipRange(A, m, n, ierr)

    if (myid .eq. 0 .and. ldebug) print *, myid, 'mat_info :: MAT_INFO_MALLOCS', mal, 'MAT_INFO_NZ_ALLOCATED', nz_allocated
    if (myid .eq. 0 .and. ldebug) print *, myid, 'mat_info :: MAT_INFO_USED', nz_used, 'MAT_INFO_NZ_unneded', nz_unneeded
    if (myid .eq. 0 .and. ldebug) print *, myid, 'mat_info :: Ownership range', m, n
  end subroutine

  subroutine set_optical_properties(solver, albedo, &
                                    kabs, ksca, g, &
                                    planck, planck_srfc, &
                                    albedo_2d, &
                                    ldelta_scaling)
    class(t_solver) :: solver
    real(ireals), intent(in) :: albedo
    real(ireals), intent(in), dimension(:, :, :), optional :: kabs, ksca, g ! dimensions (Nz  , Nx, Ny)
    real(ireals), intent(in), dimension(:, :, :), optional :: planck                    ! dimensions (Nz+1, Nx, Ny) planck radiation on levels
    real(ireals), intent(in), dimension(:, :), optional :: planck_srfc               ! dimensions (Nx, Ny)
    real(ireals), intent(in), dimension(:, :), optional :: albedo_2d                 ! dimensions (Nx, Ny)
    logical, intent(in), optional :: ldelta_scaling ! determines if we should try to delta scale these optprops

    real(ireals) :: pprts_delta_scale_max_g
    integer(iintegers) :: k, i, j
    logical :: lpprts_delta_scale, lzdun, lflg
    real(ireals) :: pprts_set_absorption, pprts_set_scatter, pprts_set_asymmetry, pprts_set_albedo
    real(irealLUT) :: c1d_dir2dir(1), c1d_dir2diff(2), c1d_diff2diff(4)
    integer(mpiint) :: ierr

    call PetscLogEventBegin(solver%logs%set_optprop, ierr); call CHKERR(ierr)

    associate ( &
        & atm => solver%atm, &
        & C_one_atm => solver%C_one_atm, &
        & C_one_atm1 => solver%C_one_atm1, &
        & sun => solver%sun, &
        & C_one => solver%C_one)

      if (.not. allocated(atm%kabs)) &
        allocate (atm%kabs(C_one_atm%zs:C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if (.not. allocated(atm%ksca)) &
        allocate (atm%ksca(C_one_atm%zs:C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if (.not. allocated(atm%g)) &
        allocate (atm%g(C_one_atm%zs:C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))

      if (.not. allocated(atm%albedo)) allocate (atm%albedo(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      atm%albedo = albedo
      if (present(albedo_2d)) atm%albedo = albedo_2d

      if (present(kabs)) atm%kabs = kabs
      if (present(ksca)) atm%ksca = ksca
      if (present(g)) atm%g = g

      if (present(planck)) then
        if (.not. allocated(atm%planck)) &
          allocate (atm%planck(C_one_atm1%zs:C_one_atm1%ze, C_one_atm1%xs:C_one_atm1%xe, C_one_atm1%ys:C_one_atm1%ye))
        atm%planck = planck
      else
        if (allocated(atm%planck)) deallocate (atm%planck)
      end if

      if (present(planck_srfc)) then
        if (.not. allocated(atm%Bsrfc)) &
          allocate (atm%Bsrfc(C_one_atm1%xs:C_one_atm1%xe, C_one_atm1%ys:C_one_atm1%ye))
        atm%Bsrfc = planck_srfc
      else
        if (allocated(atm%Bsrfc)) deallocate (atm%Bsrfc)
      end if

      if (ldebug) then
        if (any([kabs, ksca, g] .lt. zero)) then
          print *, solver%myid, 'set_optical_properties :: found illegal value in local_optical properties! abort!'
          do k = C_one_atm%zs, C_one_atm%ze
            print *, solver%myid, k, 'kabs', kabs(k, :, :)
            print *, solver%myid, k, 'ksca', ksca(k, :, :)
          end do
          call CHKERR(1_mpiint, 'set_optical_properties :: found illegal value in local_optical properties! '// &
                      ' '//toStr(solver%myid)// &
                      ' kabs min '//toStr(minval(kabs))//' max '//toStr(maxval(kabs))// &
                      ' ksca min '//toStr(minval(ksca))//' max '//toStr(maxval(ksca))// &
                      ' g    min '//toStr(minval(g))//' max '//toStr(maxval(g)))
        end if
        if (any(isnan([kabs, ksca, g]))) then
          call CHKERR(1_mpiint, 'set_optical_properties :: found NaN value in optical properties!'// &
                      ' NaN in kabs? '//toStr(any(isnan(kabs)))// &
                      ' NaN in ksca? '//toStr(any(isnan(ksca)))// &
                      ' NaN in g   ? '//toStr(any(isnan(g))))
        end if
      end if
      if (ldebug) then
        if ((any([atm%kabs, atm%ksca, atm%g] .lt. zero)) .or. (any(isnan([atm%kabs, atm%ksca, atm%g])))) then
          call CHKERR(1_mpiint, 'set_optical_properties :: found illegal value in optical properties! '// &
                      ' '//toStr(solver%myid)// &
                      ' kabs min '//toStr(minval(atm%kabs))//' max '//toStr(maxval(atm%kabs))// &
                      ' ksca min '//toStr(minval(atm%ksca))//' max '//toStr(maxval(atm%ksca))// &
                      ' g    min '//toStr(minval(atm%g))//' max '//toStr(maxval(atm%g)))
        end if
      end if

      pprts_set_albedo = -1
      call get_petsc_opt(solver%prefix, "-pprts_set_albedo", pprts_set_albedo, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        atm%albedo = pprts_set_albedo
      end if

      pprts_set_absorption = -1
      call get_petsc_opt(solver%prefix, "-pprts_set_kabs", pprts_set_absorption, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        atm%kabs = pprts_set_absorption
      end if

      pprts_set_scatter = -1
      call get_petsc_opt(solver%prefix, "-pprts_set_ksca", pprts_set_scatter, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        atm%ksca = pprts_set_scatter
      end if

      pprts_set_asymmetry = -1
      call get_petsc_opt(solver%prefix, "-pprts_set_asymmetry", pprts_set_asymmetry, lflg, ierr); call CHKERR(ierr)
      if (lflg) then
        atm%g = pprts_set_asymmetry
      end if

      lpprts_delta_scale = get_arg(.true., ldelta_scaling)
      call get_petsc_opt(solver%prefix, "-pprts_delta_scale", lpprts_delta_scale, lflg, ierr); call CHKERR(ierr)

      pprts_delta_scale_max_g = .85_ireals - epsilon(pprts_delta_scale_max_g)
      call get_petsc_opt(solver%prefix, "-pprts_delta_scale_max_g", pprts_delta_scale_max_g, lflg, ierr); call CHKERR(ierr)

      if (lpprts_delta_scale) then
        call delta_scale(atm%kabs, atm%ksca, atm%g, max_g=pprts_delta_scale_max_g)
      else
        if (solver%myid .eq. 0 .and. lflg) print *, "Skipping Delta scaling of optprops"
        if (any(atm%g .ge. 0.85_ireals)) &
          call CHKWARN(1_mpiint, 'Skipping delta scaling but now we have values of '// &
                       'g > '//toStr(pprts_delta_scale_max_g)// &
                       ' (max='//toStr(maxval(atm%g))//')')
      end if

      call print_optical_properties_summary(solver, opt_lview=ldebug)
      call dump_optical_properties(solver, ierr); call CHKERR(ierr)

      if ((any([atm%kabs, atm%ksca, atm%g] .lt. zero)) .or. (any(isnan([atm%kabs, atm%ksca, atm%g])))) then
        call CHKERR(1_mpiint, 'set_optical_properties :: found illegal value in delta_scaled optical properties!'//new_line('')// &
        & 'min(atm%kabs) '//toStr(minval(atm%kabs))//' isnan? '//toStr(any(isnan(atm%kabs)))//new_line('')// &
        & 'min(atm%ksca) '//toStr(minval(atm%ksca))//' isnan? '//toStr(any(isnan(atm%ksca)))//new_line('')// &
        & 'min(atm%g   ) '//toStr(minval(atm%g))//' isnan? '//toStr(any(isnan(atm%g)))//new_line('')// &
        & '')
      end if

      select type (solver)
      class is (t_solver_2str)
        return ! twostream should not depend on eddington coeffs... it will have to calculate it on its own.
      class is (t_solver_disort)
        return
      class is (t_solver_rayli)
        return
      class is (t_solver_mcdmda)
        return
      end select

      ! allocate space for twostream coefficients
      associate ( &
          & zs => C_one_atm%zs, ze => C_one_atm%ze, &
          & xs => C_one_atm%xs, xe => C_one_atm%xe, &
          & ys => C_one_atm%ys, ye => C_one_atm%ye)
        if (.not. allocated(atm%a11)) allocate (atm%a11(zs:ze, xs:xe, ys:ye))
        if (.not. allocated(atm%a12)) allocate (atm%a12(zs:ze, xs:xe, ys:ye))
        if (.not. allocated(atm%a21)) allocate (atm%a21(zs:ze, xs:xe, ys:ye))
        if (.not. allocated(atm%a22)) allocate (atm%a22(zs:ze, xs:xe, ys:ye))
        if (.not. allocated(atm%a13)) allocate (atm%a13(zs:ze, xs:xe, ys:ye))
        if (.not. allocated(atm%a23)) allocate (atm%a23(zs:ze, xs:xe, ys:ye))
        if (.not. allocated(atm%a33)) allocate (atm%a33(zs:ze, xs:xe, ys:ye))
      end associate

      if (luse_eddington) then
        lzdun = .false.
        call get_petsc_opt(solver%prefix, "-pprts_eddington_zdun", lzdun, lflg, ierr); call CHKERR(ierr)
        !DIR$ IVDEP
        do j = C_one_atm%ys, C_one_atm%ye
          do i = C_one_atm%xs, C_one_atm%xe
            do k = C_one_atm%zs, C_one_atm%ze
              if (atm%l1d(k)) then
                if (lzdun) then
                  call eddington_coeff_zdun( &
                    atm%dz(k, i, j) * max(tiny(one), atm%kabs(k, i, j) + atm%ksca(k, i, j)), & ! tau
                    atm%ksca(k, i, j) / max(tiny(one), atm%kabs(k, i, j) + atm%ksca(k, i, j)), & ! w0
                    atm%g(k, i, j), &
                    sun%costheta, &
                    atm%a11(k, i, j), &
                    atm%a12(k, i, j), &
                    atm%a13(k, i, j), &
                    atm%a23(k, i, j), &
                    atm%a33(k, i, j))
                else
                  call eddington_coeff_ec( &
                    atm%dz(k, i, j) * max(tiny(one), atm%kabs(k, i, j) + atm%ksca(k, i, j)), & ! tau
                    atm%ksca(k, i, j) / max(tiny(one), atm%kabs(k, i, j) + atm%ksca(k, i, j)), & ! w0
                    atm%g(k, i, j), &
                    sun%costheta, &
                    atm%a11(k, i, j), &
                    atm%a12(k, i, j), &
                    atm%a13(k, i, j), &
                    atm%a23(k, i, j), &
                    atm%a33(k, i, j))
                end if
              else
                !TODO :: we should really not have this memeory accessible at all....
                !     :: the fix would be trivial at the moment, as long as all 1d layers start at same 'k',
                !     :: however if that is to change i.e. on staggered grid, we would need staggered array constructs....
                if (ldebug) then
                  atm%a11(k, i, j) = nil
                  atm%a12(k, i, j) = nil
                  atm%a13(k, i, j) = nil
                  atm%a23(k, i, j) = nil
                  atm%a33(k, i, j) = nil
                end if !ldebug
              end if !l1d
            end do !k
          end do !i
        end do !j
      else
        do k = C_one_atm%zs, C_one_atm%ze
          if (atm%l1d(k)) then
            do j = C_one_atm%ys, C_one_atm%ye
              do i = C_one_atm%xs, C_one_atm%xe
                if (is_inrange(sun%theta, zero, 90._ireals)) then
                  call get_coeff( &
                    & solver%OPP1d, &
                    & atm%kabs(atmk(atm, k), i, j), &
                    & atm%ksca(atmk(atm, k), i, j), &
                    & atm%g(atmk(atm, k), i, j), &
                    & atm%dz(atmk(atm, k), i, j), &
                    & atm%dz(atmk(atm, k), i, j), &
                    & .true., &
                    & c1d_dir2dir, &
                    & [0._ireallut, real(sun%theta, irealLUT)], &
                    & lswitch_east=sun%xinc .eq. 0, lswitch_north=sun%yinc .eq. 0 &
                    & )
                  atm%a33(k, i, j) = real(c1d_dir2dir(1), ireals)

                  call get_coeff( &
                    & solver%OPP1d, &
                    & atm%kabs(atmk(atm, k), i, j), &
                    & atm%ksca(atmk(atm, k), i, j), &
                    & atm%g(atmk(atm, k), i, j), &
                    & atm%dz(atmk(atm, k), i, j), &
                    & atm%dz(atmk(atm, k), i, j), &
                    & .false., &
                    & c1d_dir2diff, &
                    & [0._ireallut, real(sun%theta, irealLUT)], &
                    & lswitch_east=sun%xinc .eq. 0, lswitch_north=sun%yinc .eq. 0 &
                    & )
                  atm%a13(k, i, j) = real(c1d_dir2diff(1), ireals)
                  atm%a23(k, i, j) = real(c1d_dir2diff(2), ireals)
                end if

                call get_coeff( &
                  & solver%OPP1d, &
                  & atm%kabs(atmk(atm, k), i, j), &
                  & atm%ksca(atmk(atm, k), i, j), &
                  & atm%g(atmk(atm, k), i, j), &
                  & atm%dz(atmk(atm, k), i, j), &
                  & atm%dz(atmk(atm, k), i, j), &
                  & .false., &
                  & c1d_diff2diff &
                  & )
                atm%a11(k, i, j) = real(c1d_diff2diff(1), ireals)
                atm%a12(k, i, j) = real(c1d_diff2diff(2), ireals)
                atm%a21(k, i, j) = real(c1d_diff2diff(3), ireals)
                atm%a22(k, i, j) = real(c1d_diff2diff(4), ireals)
              end do !i
            end do !j
          end if !l1d
        end do !k
      end if

      ! Set symmetric up and down transport coefficients. If only for one
      ! layer, they would be indeed symmetric. If we collapse the atmosphere
      ! however, we have different transmission if we come from top or from
      ! bottom.
      atm%a21 = atm%a12
      atm%a22 = atm%a11
      call handle_atm_collapse()

    end associate
    call PetscLogEventEnd(solver%logs%set_optprop, ierr); call CHKERR(ierr)

  contains
    subroutine handle_atm_collapse()
      real(ireals), allocatable :: Eup(:), Edn(:)
      integer(iintegers) :: i, j, ak
      associate (atm => solver%atm, &
                 C_one => solver%C_one, &
                 C_one_atm => solver%C_one_atm)

        if (atm%lcollapse) then
          ak = atmk(atm, C_one%zs)
          if (present(planck)) then
            allocate (Edn(C_one_atm%zs:ak + 1), Eup(C_one_atm%zs:ak + 1))
            if (.not. allocated(atm%Bbot)) &
              allocate (atm%Bbot(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
            if (.not. allocated(atm%Btop)) &
              allocate (atm%Btop(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          end if
          do j = C_one_atm%ys, C_one_atm%ye
            do i = C_one_atm%xs, C_one_atm%xe
              if (present(planck)) then
                call adding( &
                  atm%a11(C_one_atm%zs:ak, i, j), &
                  atm%a12(C_one_atm%zs:ak, i, j), &
                  atm%a21(C_one_atm%zs:ak, i, j), &
                  atm%a22(C_one_atm%zs:ak, i, j), &
                  atm%a13(C_one_atm%zs:ak, i, j), &
                  atm%a23(C_one_atm%zs:ak, i, j), &
                  atm%a33(C_one_atm%zs:ak, i, j), &
                  atm%dz(C_one_atm%zs:ak, i, j) * atm%kabs(C_one_atm%zs:ak, i, j), &
                  atm%planck(C_one_atm%zs:ak + 1, i, j), &
                  Eup, Edn, atm%Btop(i, j), atm%Bbot(i, j))
              else
                call adding( &
                  atm%a11(C_one_atm%zs:ak, i, j), &
                  atm%a12(C_one_atm%zs:ak, i, j), &
                  atm%a21(C_one_atm%zs:ak, i, j), &
                  atm%a22(C_one_atm%zs:ak, i, j), &
                  atm%a13(C_one_atm%zs:ak, i, j), &
                  atm%a23(C_one_atm%zs:ak, i, j), &
                  atm%a33(C_one_atm%zs:ak, i, j))
              end if
            end do !i
          end do !j
        end if !lcollapse
      end associate
    end subroutine
    subroutine adding(a11, a12, a21, a22, a13, a23, a33, dtau, planck, Eup, Edn, Btop, Bbot)
      real(ireals), intent(inout), dimension(:) :: a11, a12, a21, a22, a13, a23, a33
      real(ireals), intent(in), dimension(:), optional :: dtau, planck
      real(ireals), intent(out), dimension(:), optional :: Eup, Edn
      real(ireals), intent(out), optional :: Btop, Bbot
      real(ireals) :: t, r, rdir, sdir, tdir

      integer(iintegers) :: N ! = size(a11)
      integer(iintegers) :: k
      real(ireals) :: rl, tl, Tbot, Ttop, Rbot, Rtop
      N = size(a11)

      t = a11(1)
      r = a12(1)

      tdir = a33(1)
      rdir = a13(1)
      sdir = a23(1)

      ! Reflectivity as seen from top
      do k = 2, N
        rl = r
        tl = t

        r = r + (a12(k) * t**2) / (one - r * a12(k))
        t = t * a11(k) / (one - rl * a12(k))

        sdir = (a11(k) * sdir + tdir * a13(k) * rl * a11(k)) / (one - rl * a12(k)) + tdir * a23(k)
        rdir = rdir + (tdir * a13(k) + sdir * a12(k)) * tl
        tdir = tdir * a33(k)
      end do

      Ttop = t
      Rtop = r

      a13(N) = rdir
      a23(N) = sdir
      a33(N) = tdir

      t = a22(N)
      r = a21(N)
      ! Reflectivity as seen from bottom
      do k = N - 1, 1, -1
        rl = r
        tl = t

        r = a12(k) + (r * a11(k)**2) / (one - r * a12(k))
        t = t * a11(k) / (one - rl * a21(k))
      end do

      Tbot = t
      Rbot = r

      a12(N) = Rtop
      a22(N) = Ttop
      a11(N) = Tbot
      a21(N) = Rbot

      a11(1:N - 1) = nil
      a12(1:N - 1) = nil
      a21(1:N - 1) = nil
      a22(1:N - 1) = nil

      a13(1:N - 1) = nil
      a23(1:N - 1) = nil
      a33(1:N - 1) = nil

      if (present(planck)) then
        call schwarzschild(2_iintegers, dtau, 0._ireals, Edn, Eup, planck, &
                           opt_srfc_emission=0._ireals)
        Bbot = Edn(N + 1) / pi
        Btop = Eup(1) / pi
      end if
    end subroutine
  end subroutine

  subroutine print_optical_properties_summary(solver, opt_lview)
    class(t_solver), intent(in) :: solver
    logical, intent(in), optional :: opt_lview

    type(tVec) :: valbedo
    logical :: lview, lflg
    real(ireals), dimension(3) :: mkabs, mksca, mg, malbedo, mplck
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: k

    lview = get_arg(.false., opt_lview)
    call get_petsc_opt(solver%prefix, "-pprts_view_optprop", lview, lflg, ierr); call CHKERR(ierr)
    if (.not. lview) return

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
    if (myid .eq. 0) print *, '*      k(layer)  '// &
      '     '//cstr('kabs(min/mean/max)', 'blue')//'                  '// &
      '     '//cstr('ksca              ', 'red')//'                  '// &
      '     '//cstr('g                 ', 'blue')//'                  '// &
      '     '//cstr('plck              ', 'red')
    associate (atm => solver%atm, C_one_atm => solver%C_one_atm)

      do k = C_one_atm%zs, C_one_atm%ze

        if (allocated(atm%kabs)) call imp_min_mean_max(solver%comm, atm%kabs(k, :, :), mkabs)
        if (allocated(atm%ksca)) call imp_min_mean_max(solver%comm, atm%ksca(k, :, :), mksca)
        if (allocated(atm%g)) call imp_min_mean_max(solver%comm, atm%g(k, :, :), mg)
        if (allocated(atm%planck)) then
          call imp_min_mean_max(solver%comm, atm%planck(k, :, :), mplck)
        else
          mplck = nan
        end if

        if (myid .eq. 0) &
          & print *, k, cstr(toStr(mkabs), 'blue'), cstr(toStr(mksca), 'red'), cstr(toStr(mg), 'blue'), cstr(toStr(mplck), 'red')

      end do

      if (allocated(atm%albedo)) then
        call imp_min_mean_max(solver%comm, atm%albedo(:, :), malbedo)
        if (myid .eq. 0) print *, ' * Albedo (min/mean,max)', malbedo

        call DMGetGlobalVector(solver%Csrfc_one%da, valbedo, ierr); call CHKERR(ierr)
        call f90VecToPetsc(atm%albedo, solver%Csrfc_one%da, valbedo)
        call PetscObjectViewFromOptions(valbedo, solver%Csrfc_one%da, "-pprts_view_albedo", ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(solver%Csrfc_one%da, valbedo, ierr); call CHKERR(ierr)
      end if

      if (myid .eq. 0) then
        print *, ' * Number of 1D layers: ', count(atm%l1d), size(atm%l1d), &
          & '(', (100._ireals * count(atm%l1d)) / size(atm%l1d), '%)'
      end if

    end associate
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
  end subroutine

  subroutine dump_optical_properties(solver, ierr)
    class(t_solver), intent(in) :: solver
    integer(mpiint), intent(out) :: ierr

    character(len=default_str_len) :: fname
    logical :: lflg

    character(len=default_str_len) :: dimnames(3)

    ierr = 0

    call get_petsc_opt(solver%prefix, '-pprts_dump_optprop', fname, lflg, ierr); call CHKERR(ierr)
    if (lflg) then
      if (len_trim(fname) .eq. 0) then
        ierr = 1
        call CHKERR(ierr, '-pprts_dump_optprop options needs file name as argument')
      end if
      dimnames(1) = 'nlay'
      dimnames(2) = 'nx'
      dimnames(3) = 'ny'
      call dump_var(solver%C_one_atm, solver%atm%kabs, 'kabs', 'm-1', dimnames, ierr); call CHKERR(ierr)
      call dump_var(solver%C_one_atm, solver%atm%ksca, 'ksca', 'm-1', dimnames, ierr); call CHKERR(ierr)
      call dump_var(solver%C_one_atm, solver%atm%g, 'g', '', dimnames, ierr); call CHKERR(ierr)

      if (allocated(solver%atm%planck)) then
        dimnames(1) = 'nlev'
        call dump_var(solver%C_one_atm1, solver%atm%planck, 'planck', 'W m-2', dimnames, ierr); call CHKERR(ierr)
      end if
    end if

  contains
    subroutine dump_var(C, var, varname, units, dimnames, ierr)
      type(t_coord), intent(in) :: C
      real(ireals), intent(in) :: var(:, :, :)
      character(len=*), intent(in) :: varname, units, dimnames(:)
      integer(mpiint), intent(out) :: ierr

      real(ireals), allocatable :: var0(:, :, :)
      character(len=default_str_len) :: groups(2)
      logical :: var_exists
      integer(iintegers) :: prefixid

      call gather_all_toZero(C, var, var0)

      if (solver%myid .eq. 0) then
        groups(1) = trim(fname)
        var_exists = .true.
        prefixid = 0
        do while (var_exists)
          groups(2) = trim(varname)//'.'//trim(toStr(prefixid))

          call nc_var_exists(groups, var_exists, ierr, verbose=.false.)
          prefixid = prefixid + 1
        end do
        call ncwrite(groups, var0, ierr, dimnames=dimnames, verbose=.false., &
          & deflate_lvl=5, &
          & chunksizes=[integer :: 1, min(64, size(var0, 2)), min(64, size(var0, 3))]); call CHKERR(ierr)
        call set_attribute(groups(1), groups(2), 'units', units, ierr); call CHKERR(ierr)
      end if
      ierr = 0
    end subroutine
  end subroutine

  subroutine set_global_optical_properties(solver, global_albedo, global_kabs, global_ksca, global_g, global_planck)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(inout), optional :: global_albedo
    real(ireals), intent(inout), dimension(:, :, :), allocatable, optional :: global_kabs, global_ksca, global_g
    real(ireals), intent(inout), dimension(:, :, :), allocatable, optional :: global_planck
    real(ireals), dimension(:, :, :), allocatable :: local_kabs, local_ksca, local_g
    real(ireals), dimension(:, :, :), allocatable :: local_planck
    real(ireals) :: local_albedo
    logical :: lhave_planck, lhave_kabs, lhave_ksca, lhave_g
    integer(mpiint) :: ierr

    call PetscLogEventBegin(solver%logs%set_optprop, ierr); call CHKERR(ierr)

    if (.not. solver%linitialized) then
      call CHKERR(1_mpiint, 'You tried to set global optical properties but pprts environment seems not to be initialized....'// &
                  'please call init first!')
    end if

    lhave_kabs = present(global_kabs); call imp_bcast(solver%comm, lhave_kabs, 0_mpiint, ierr); call CHKERR(ierr)
    lhave_ksca = present(global_ksca); call imp_bcast(solver%comm, lhave_ksca, 0_mpiint, ierr); call CHKERR(ierr)
    lhave_g = present(global_g); call imp_bcast(solver%comm, lhave_g, 0_mpiint, ierr); call CHKERR(ierr)
    lhave_planck = present(global_planck); call imp_bcast(solver%comm, lhave_planck, 0_mpiint, ierr); call CHKERR(ierr)

    ! Make sure that our domain has at least 3 entries in each dimension.... otherwise violates boundary conditions
    if (solver%myid .eq. 0) then
      if (present(global_albedo)) local_albedo = global_albedo
      if (lhave_kabs) call extend_arr(global_kabs)
      if (lhave_ksca) call extend_arr(global_ksca)
      if (lhave_g) call extend_arr(global_g)
      if (lhave_planck) call extend_arr(global_planck)
    end if

    ! Scatter global optical properties to MPI nodes
    call local_optprop()
    ! Now global_fields are local to mpi subdomain.
    call PetscLogEventEnd(solver%logs%set_optprop, ierr); call CHKERR(ierr)

    if (lhave_planck) then
      call set_optical_properties(solver, local_albedo, local_kabs, local_ksca, local_g, local_planck)
    else
      call set_optical_properties(solver, local_albedo, local_kabs, local_ksca, local_g)
    end if
  contains
    subroutine local_optprop()
      type(tVec) :: local_vec

      if (solver%myid .eq. 0 .and. ldebug .and. lhave_kabs) &
        print *, solver%myid, 'copying optprop: global to local :: shape kabs', shape(global_kabs), &
        'xstart/end', solver%C_one_atm%xs, solver%C_one_atm%xe, &
        'ystart/end', solver%C_one_atm%ys, solver%C_one_atm%ye

      call imp_bcast(solver%comm, local_albedo, 0_mpiint, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(solver%C_one_atm%da, local_vec, ierr); call CHKERR(ierr)

      if (lhave_kabs) then
        call scatterZerotoPetscGlobal(global_kabs, solver%C_one_atm%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm%da, local_kabs)
      end if

      if (lhave_ksca) then
        call scatterZerotoPetscGlobal(global_ksca, solver%C_one_atm%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm%da, local_ksca)
      end if

      if (lhave_g) then
        call scatterZerotoPetscGlobal(global_g, solver%C_one_atm%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm%da, local_g)
      end if

      call DMRestoreGlobalVector(solver%C_one_atm%da, local_vec, ierr); call CHKERR(ierr)

      if (lhave_planck) then
        call DMGetGlobalVector(solver%C_one_atm1%da, local_vec, ierr); call CHKERR(ierr)
        call scatterZerotoPetscGlobal(global_planck, solver%C_one_atm1%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm1%da, local_planck)
        call DMRestoreGlobalVector(solver%C_one_atm1%da, local_vec, ierr); call CHKERR(ierr)
      end if
    end subroutine
    subroutine extend_arr(arr)
      real(ireals), intent(inout), allocatable :: arr(:, :, :)
      real(ireals), allocatable :: tmp(:, :)
      integer(iintegers) :: dims(3), i

      if (.not. allocated(arr)) print *, solver%myid, 'ERROR in SUBROUTINE extend_arr :: Cannot extend non allocated array!'

      dims = shape(arr)
      if (dims(2) .eq. 1) then
        allocate (tmp(dims(1), dims(3)), SOURCE=arr(:, 1, :))
        deallocate (arr)
        allocate (arr(dims(1), minimal_dimension, dims(3)))
        do i = 1, minimal_dimension
          arr(:, i, :) = tmp
        end do
        deallocate (tmp)
      end if

      dims = shape(arr)
      if (dims(3) .eq. 1) then
        allocate (tmp(dims(1), dims(2)), SOURCE=arr(:, :, 1))
        deallocate (arr)
        allocate (arr(dims(1), dims(2), minimal_dimension))
        do i = 1, minimal_dimension
          arr(:, :, i) = tmp
        end do
        deallocate (tmp)
      end if
      if (any(shape(arr) .lt. minimal_dimension)) &
        call CHKERR(1_mpiint, 'set_optprop -> extend_arr :: dimension is smaller than we support.. please think of something here')
    end subroutine

  end subroutine

  recursive subroutine solve_pprts(solver, lthermal, lsolar, edirTOA, opt_solution_uid, opt_solution_time, opt_buildings)
    class(t_solver), intent(inout) :: solver
    logical, intent(in) :: lthermal, lsolar
    real(ireals), intent(in) :: edirTOA
    integer(iintegers), optional, intent(in) :: opt_solution_uid
    real(ireals), optional, intent(in) :: opt_solution_time
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    integer(iintegers) :: uid, last_uid
    logical :: derived_lsolar, luse_rayli, lrayli_snapshot
    logical :: linitial_guess_from_last_uid, linitial_guess_from_2str, linitial_guess_from_sc, lflg
    integer(mpiint) :: ierr

    if (.not. allocated(solver%atm)) call CHKERR(1_mpiint, 'atmosphere is not allocated?!')
    if (.not. allocated(solver%atm%kabs)) &
      call CHKERR(1_mpiint, 'atmosphere%kabs is not allocated! - maybe you need to call set_optical_properties() first')
    if (.not. allocated(solver%atm%ksca)) &
      call CHKERR(1_mpiint, 'atmosphere%ksca is not allocated! - maybe you need to call set_optical_properties() first')

    uid = get_solution_uid(solver%solutions, opt_solution_uid)

    associate (solution => solver%solutions(uid))

      if (lthermal .and. (.not. allocated(solver%atm%planck))) &
        & call CHKERR(1_mpiint, 'you asked to compute thermal radiation but '// &
        & 'did not provide planck emissions in set_optical_properties')

      derived_lsolar = mpi_logical_and(solver%comm, &
        & lsolar .and. &
        & edirTOA .gt. zero .and. &
        & solver%sun%theta .ge. zero)

      if (ldebug) print *, 'uid', uid, 'lsolar', derived_lsolar, &
        & 'edirTOA', edirTOA, 'any_theta in [0..90]?', solver%sun%theta .ge. zero

      if (derived_lsolar .and. lthermal) then
        print *, 'uid', uid, 'lthermal', lthermal, 'lsolar', lsolar, derived_lsolar, &
          & 'edirTOA', edirTOA, 'any_theta in [0..90]?', solver%sun%theta .ge. zero
        call CHKERR(1_mpiint, &
          & 'Somehow ended up with a request to compute solar and thermal radiation in one call.'//new_line('')// &
          & '       This may currently not work or at least has to be tested further...'//new_line('')// &
          & '       I recommend you call it one after another')
      end if

      if (.not. solution%lset) then
        call prepare_solution(solver%C_dir%da, solver%C_diff%da, solver%C_one%da, &
                              lsolar=derived_lsolar, lthermal=lthermal, solution=solution, uid=uid)

        call PetscLogEventBegin(solver%logs%setup_initial_guess, ierr); call CHKERR(ierr)
        linitial_guess_from_last_uid = .true.
        call get_petsc_opt('', "-initial_guess_from_last_uid", &
          & linitial_guess_from_last_uid, lflg, ierr); call CHKERR(ierr)
        call get_petsc_opt(solver%prefix, "-initial_guess_from_last_uid", &
          & linitial_guess_from_last_uid, lflg, ierr); call CHKERR(ierr)
        if (linitial_guess_from_last_uid) then
          last_uid = get_solution_uid(solver%solutions, uid - 1)
          if (solver%solutions(last_uid)%lset) then
            if (solver%solutions(last_uid)%lsolar_rad) then
              call VecCopy(solver%solutions(last_uid)%edir, solution%edir, ierr); call CHKERR(ierr)
              solution%lWm2_dir = solver%solutions(last_uid)%lWm2_dir
            end if
            call VecCopy(solver%solutions(last_uid)%ediff, solution%ediff, ierr); call CHKERR(ierr)
            solution%lWm2_diff = solver%solutions(last_uid)%lWm2_diff
          end if
        end if

        linitial_guess_from_2str = .false.
        call get_petsc_opt('', "-initial_guess_from_2str", linitial_guess_from_2str, lflg, ierr); call CHKERR(ierr)
        call get_petsc_opt(solver%prefix, "-initial_guess_from_2str", linitial_guess_from_2str, lflg, ierr); call CHKERR(ierr)
        if (linitial_guess_from_2str) then
          call PetscLogEventBegin(solver%logs%solve_twostream, ierr); call CHKERR(ierr)
          call twostream(solver, edirTOA, solution, opt_buildings)
          call PetscLogEventEnd(solver%logs%solve_twostream, ierr); call CHKERR(ierr)
        end if

        linitial_guess_from_sc = .false.
        call get_petsc_opt('', "-initial_guess_from_sc", linitial_guess_from_sc, lflg, ierr); call CHKERR(ierr)
        call get_petsc_opt(solver%prefix, "-initial_guess_from_sc", linitial_guess_from_sc, lflg, ierr); call CHKERR(ierr)
        if (linitial_guess_from_sc) then
          call initial_guess_from_single_column_comp(solver, edirTOA, solution, ierr, opt_buildings); call CHKERR(ierr)
        end if
        call PetscLogEventEnd(solver%logs%setup_initial_guess, ierr); call CHKERR(ierr)

      else
        if (solution%lsolar_rad .neqv. derived_lsolar) then
          call destroy_solution(solution)
          call prepare_solution(solver%C_dir%da, solver%C_diff%da, solver%C_one%da, &
                                lsolar=derived_lsolar, lthermal=lthermal, solution=solution, uid=uid)
        end if
      end if

      if ((solution%lthermal_rad .eqv. .false.) .and. (solution%lsolar_rad .eqv. .false.)) then ! nothing else to do
        if (allocated(solution%edir)) then
          call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
        end if
        if (allocated(solution%ediff)) then
          call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
        end if
        solution%lchanged = .true.
        goto 99 ! quick exit
      end if

      ! --------- Skip Thermal Computation (-lskip_thermal) --
      if (lskip_thermal .and. (solution%lsolar_rad .eqv. .false.)) then ! nothing else to do
        if (ldebug .and. solver%myid .eq. 0) print *, 'skipping thermal calculation -- returning zero flux'
        call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
        solution%lchanged = .true.
        goto 99 ! quick exit
      end if

      ! --------- Calculate Radiative Transfer with RayLi ------------
      luse_rayli = .false.
      select type (solver)
      class is (t_solver_rayli)
        luse_rayli = .true.
      end select

      lrayli_snapshot = .false.
      call PetscOptionsHasName(PETSC_NULL_OPTIONS, solver%prefix, &
                               "-rayli_snapshot", lrayli_snapshot, ierr); call CHKERR(ierr)

      if (luse_rayli .or. lrayli_snapshot) then
        call PetscLogEventBegin(solver%logs%solve_mcrts, ierr); call CHKERR(ierr)
        call pprts_rayli_wrapper(luse_rayli, lrayli_snapshot, solver, edirTOA, solution, ierr, opt_buildings); call CHKERR(ierr)
        call PetscLogEventEnd(solver%logs%solve_mcrts, ierr); call CHKERR(ierr)
        if (luse_rayli) goto 99
      end if

      select type (solver)
      class is (t_solver_2str)
        if (solution%lthermal_rad .and. lschwarzschild) then
          call PetscLogEventBegin(solver%logs%solve_schwarzschild, ierr); call CHKERR(ierr)
          call schwarz(solver, solution, opt_buildings)
          call PetscLogEventEnd(solver%logs%solve_schwarzschild, ierr); call CHKERR(ierr)
        else
          call PetscLogEventBegin(solver%logs%solve_twostream, ierr); call CHKERR(ierr)
          call twostream(solver, edirTOA, solution, opt_buildings)
          call PetscLogEventEnd(solver%logs%solve_twostream, ierr); call CHKERR(ierr)
        end if
        goto 99

      class is (t_solver_disort)
        call PetscLogEventBegin(solver%logs%solve_disort, ierr); call CHKERR(ierr)
        call disort(solver, edirTOA, solution, opt_buildings)
        call PetscLogEventEnd(solver%logs%solve_disort, ierr); call CHKERR(ierr)
        goto 99

      class is (t_solver_mcdmda)
        call PetscLogEventBegin(solver%logs%solve_mcrts, ierr); call CHKERR(ierr)
        call solve_mcdmda(solver, edirTOA, solution, ierr, opt_buildings); call CHKERR(ierr)
        call PetscLogEventEnd(solver%logs%solve_mcrts, ierr); call CHKERR(ierr)
        goto 99

      end select

      call pprts(solver, edirTOA, solution, opt_buildings)

99    continue ! this is the quick exit final call where we clean up before the end of the routine

      call restore_solution(solver, solution, opt_solution_time)

    end associate
  end subroutine

  !> @brief call the matrix assembly and petsc solve routines for pprts solvers
  subroutine pprts(solver, edirTOA, solution, opt_buildings)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), intent(inout) :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    logical :: lflg, lskip_diffuse_solve
    logical :: lexplicit_dir, lexplicit_diff
    real(ireals) :: b_norm, rtol, atol
    integer(mpiint) :: ierr

    character(len=default_str_len) :: prefix

    ! Populate transport coeffs
    if (solution%lsolar_rad) then
      call alloc_coeff_dir2dir(solver, solver%dir2dir, opt_buildings)
      call alloc_coeff_dir2diff(solver, solver%dir2diff)
    end if
    call alloc_coeff_diff2diff(solver, solver%diff2diff, opt_buildings)

    ! --------- scale from [W/m**2] to [W] -----------------
    call scale_flx(solver, solution, lWm2=.false.)

    ! ---------------------------- Edir  -------------------
    if (solution%lsolar_rad) then
      call PetscLogEventBegin(solver%logs%compute_Edir, ierr)
      prefix = "solar_dir_"
      if (len_trim(solver%prefix) .gt. 0) prefix = trim(solver%prefix)//prefix

      lexplicit_dir = .true.
      call get_petsc_opt(prefix, "-explicit", lexplicit_dir, lflg, ierr); call CHKERR(ierr)
      if (lexplicit_dir) then
        call explicit_edir(solver, prefix, edirTOA, solution, ierr); call CHKERR(ierr)
      else
        call edir(prefix)
      end if
      call PetscLogEventEnd(solver%logs%compute_Edir, ierr)
    end if

    ! ---------------------------- Source Term -------------
    call PetscLogEventBegin(solver%logs%compute_Ediff, ierr); call CHKERR(ierr)
    call setup_b(solver, solution, solver%b, opt_buildings)

    if (solution%lsolar_rad) then
      call determine_ksp_tolerances(solver%C_diff, solver%atm%l1d, rtol, atol, solver%ksp_solar_diff)
    else
      call determine_ksp_tolerances(solver%C_diff, solver%atm%l1d, rtol, atol, solver%ksp_thermal_diff)
    end if

    ! ---------------------------- Ediff -------------------
    call VecNorm(solver%b, NORM_1, b_norm, ierr); call CHKERR(ierr)
    lskip_diffuse_solve = b_norm .lt. atol
    call get_petsc_opt(solver%prefix, "-skip_diffuse_solve", lskip_diffuse_solve, lflg, ierr); call CHKERR(ierr)

    if (lskip_diffuse_solve) then
      if (ldebug) &
        & call CHKWARN(1_mpiint, 'Skipping diffuse solver :: b_norm='//toStr(b_norm))
      call VecCopy(solver%b, solution%ediff, ierr); call CHKERR(ierr)
      solution%Niter_diff = i0
      solution%diff_ksp_residual_history = atol
    else

      if (solution%lsolar_rad) then
        prefix = "solar_diff_"
      else
        prefix = "thermal_diff_"
      end if
      if (len_trim(solver%prefix) .gt. 0) prefix = trim(solver%prefix)//prefix

      lexplicit_diff = .false.
      call get_petsc_opt(prefix, "-explicit", lexplicit_diff, lflg, ierr); call CHKERR(ierr)
      if (lexplicit_diff) then
        call explicit_ediff(solver, prefix, solver%b, solution, ierr); call CHKERR(ierr)
      else
        if (solution%lsolar_rad) then
          call ediff(solver%Mdiff, solver%Mdiff_perm, solver%ksp_solar_diff, prefix)
        else
          call ediff(solver%Mth, solver%Mth_perm, solver%ksp_thermal_diff, prefix)
        end if
      end if
    end if

    solution%lchanged = .true.
    solution%lWm2_diff = .false. !Tenstream solver returns fluxes as [W]

    call PetscLogEventEnd(solver%logs%compute_Ediff, ierr)

  contains

    subroutine edir(prefix)
      character(len=*), intent(in) :: prefix
      logical :: lmat_permute, lmat_permute_reuse, lshell

      call VecSet(solver%incSolar, zero, ierr); call CHKERR(ierr)
      call setup_incSolar(solver, edirTOA, solver%incSolar)

      lshell = .false.
      call get_petsc_opt(prefix, "-shell", lshell, lflg, ierr); call CHKERR(ierr)

      if (lshell) then
        call setup_matshell(solver, solver%C_dir, solver%Mdir, op_mat_mult_edir, op_mat_sor_edir, op_mat_getdiagonal)
      else
        call PetscLogEventBegin(solver%logs%setup_Mdir, ierr)
        call init_Matrix(solver, solver%C_dir, solver%Mdir, setup_direct_preallocation)
        call set_dir_coeff(solver, solver%sun, solver%Mdir, solver%C_dir)
        call PetscLogEventEnd(solver%logs%setup_Mdir, ierr)
        if (ldebug) call mat_info(solver%comm, solver%Mdir)
      end if

      lmat_permute = .false.
      call get_petsc_opt(prefix, "-mat_permute", lmat_permute, lflg, ierr); call CHKERR(ierr)
      if (lmat_permute) then ! prepare mat permutation
        call PetscLogEventBegin(solver%logs%permute_mat_gen_dir, ierr)
        call gen_mat_permutation( &
          & A=solver%Mdir, &
          & C=solver%C_dir, &
          & rev_x=solver%sun%xinc .eq. 0, & ! reverse if sun is east
          & rev_y=solver%sun%yinc .eq. 0, & ! reverse if sun is north, &
          & rev_z=.false., &
          & zlast=.true., &
          & switch_xy=abs(solver%sun%sundir(2)) .gt. abs(solver%sun%sundir(1)), &
          & perm_info=solver%perm_dir, &
          & prefix=solver%prefix)
        call PetscLogEventEnd(solver%logs%permute_mat_gen_dir, ierr)

        call PetscLogEventBegin(solver%logs%permute_mat_dir, ierr)
        call VecPermute(solver%incSolar, solver%perm_dir%is, PETSC_FALSE, ierr); call CHKERR(ierr)
        call VecPermute(solution%edir, solver%perm_dir%is, PETSC_FALSE, ierr); call CHKERR(ierr)
        if (.not. allocated(solver%Mdir_perm)) then
          allocate (solver%Mdir_perm)
          !call MatPermute(solver%Mdir, solver%perm_dir%is, solver%perm_dir%is, solver%Mdir_perm, ierr); call CHKERR(ierr)
          call MatCreateSubMatrix(solver%Mdir, solver%perm_dir%is, solver%perm_dir%is, &
            & MAT_INITIAL_MATRIX, solver%Mdir_perm, ierr); call CHKERR(ierr)
        else

          lmat_permute_reuse = .false.
          call get_petsc_opt(prefix, "-mat_permute_reuse", lmat_permute_reuse, lflg, ierr); call CHKERR(ierr)
          ! TODO: we should be able to reuse the mat, but results change, I could not figure out why.
          ! Instead we destroy the mat and build and initial one every time.
          ! This hits performance but is still better than not reordering...
          if (lmat_permute_reuse) then
            if (solver%myid .eq. 0) &
              & call CHKWARN(1_mpiint, 'mat_permute_reuse &
              & :: not sure if this always gives the correct results, please check carefully')
            call MatCreateSubMatrix(solver%Mdir, solver%perm_dir%is, solver%perm_dir%is, &
              & MAT_REUSE_MATRIX, solver%Mdir_perm, ierr); call CHKERR(ierr)
          else
            call MatDestroy(solver%Mdir_perm, ierr); call CHKERR(ierr); 
            !call MatPermute(solver%Mdir, solver%perm_dir%is, solver%perm_dir%is, solver%Mdir_perm, ierr); call CHKERR(ierr)
            call MatCreateSubMatrix(solver%Mdir, solver%perm_dir%is, solver%perm_dir%is, &
              & MAT_INITIAL_MATRIX, solver%Mdir_perm, ierr); call CHKERR(ierr)
          end if
          call KSPSetOperators(solver%ksp_solar_dir, solver%Mdir_perm, solver%Mdir_perm, ierr); call CHKERR(ierr)
        end if
        call PetscLogEventEnd(solver%logs%permute_mat_dir, ierr)

        call setup_ksp(solver, solver%ksp_solar_dir, solver%C_dir, solver%Mdir_perm, prefix=prefix)
      else
        call setup_ksp(solver, solver%ksp_solar_dir, solver%C_dir, solver%Mdir, prefix=prefix)
      end if

      call PetscLogEventBegin(solver%logs%solve_Mdir, ierr)
      call solve(solver, &
                 solver%ksp_solar_dir, &
                 solver%incSolar, &
                 solution%edir, &
                 solution%uid, &
                 solution%Niter_dir, &
                 solution%dir_ksp_residual_history)
      call PetscLogEventEnd(solver%logs%solve_Mdir, ierr)

      if (lmat_permute) then
        call PetscLogEventBegin(solver%logs%permute_mat_dir, ierr)
        call VecPermute(solver%incSolar, solver%perm_dir%is, PETSC_TRUE, ierr); call CHKERR(ierr)
        call VecPermute(solution%edir, solver%perm_dir%is, PETSC_TRUE, ierr); call CHKERR(ierr)
        call PetscLogEventEnd(solver%logs%permute_mat_dir, ierr)
      end if

      solution%lchanged = .true.
      solution%lWm2_dir = .false.
      call PetscObjectSetName(solution%edir, 'debug_edir', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_debug_edir", ierr); call CHKERR(ierr)
    end subroutine

    subroutine ediff(A, Aperm, ksp, prefix)
      type(tMat), allocatable, intent(inout) :: A, Aperm
      type(tKSP), allocatable, intent(inout) :: ksp
      character(len=*), intent(in) :: prefix

      logical :: lmat_permute, lmat_permute_reuse, lshell

      lshell = .false.
      call get_petsc_opt(prefix, "-shell", lshell, lflg, ierr); call CHKERR(ierr)

      if (lshell) then
        call setup_matshell(solver, solver%C_diff, A, op_mat_mult_ediff, op_mat_sor_ediff, op_mat_getdiagonal)
      else
        call init_Matrix(solver, solver%C_diff, A, setup_diffuse_preallocation)

        call PetscLogEventBegin(solver%logs%setup_Mdiff, ierr)
        call set_diff_coeff(solver, A, solver%C_diff)
        call PetscLogEventEnd(solver%logs%setup_Mdiff, ierr)
        if (ldebug) call mat_info(solver%comm, A)
      end if

      lmat_permute = .false.
      call get_petsc_opt(prefix, "-mat_permute", lmat_permute, lflg, ierr); call CHKERR(ierr)

      if (lmat_permute) then ! prepare mat permutation

        call PetscLogEventBegin(solver%logs%permute_mat_gen_diff, ierr)
        call gen_mat_permutation( &
          & A=A, &
          & C=solver%C_diff, &
          & rev_x=.false., &
          & rev_y=.false., &
          & rev_z=.false., &
          & zlast=.true., &
          & switch_xy=.false., &
          & perm_info=solver%perm_diff, &
          & prefix=solver%prefix)
        call PetscLogEventEnd(solver%logs%permute_mat_gen_diff, ierr)

        call PetscLogEventBegin(solver%logs%permute_mat_diff, ierr)

        call VecPermute(solver%b, solver%perm_diff%is, PETSC_FALSE, ierr); call CHKERR(ierr)
        call VecPermute(solution%ediff, solver%perm_diff%is, PETSC_FALSE, ierr); call CHKERR(ierr)

        if (.not. allocated(Aperm)) then
          allocate (Aperm)
          call MatCreateSubMatrix(A, solver%perm_diff%is, solver%perm_diff%is, &
            & MAT_INITIAL_MATRIX, Aperm, ierr); call CHKERR(ierr)
        else

          lmat_permute_reuse = .false.
          call get_petsc_opt(prefix, "-mat_permute_reuse", lmat_permute_reuse, lflg, ierr); call CHKERR(ierr)

          if (lmat_permute_reuse) then
            if (solver%myid .eq. 0) &
              & call CHKWARN(1_mpiint, 'mat_permute_reuse :: &
              & not sure if this always gives the correct results, please check carefully')
            call MatCreateSubMatrix(A, solver%perm_diff%is, solver%perm_diff%is, &
              & MAT_REUSE_MATRIX, Aperm, ierr); call CHKERR(ierr)
          else
            call MatDestroy(Aperm, ierr); call CHKERR(ierr)
            call MatCreateSubMatrix(A, solver%perm_diff%is, solver%perm_diff%is, &
              & MAT_INITIAL_MATRIX, Aperm, ierr); call CHKERR(ierr)
          end if
          call KSPSetOperators(ksp, Aperm, Aperm, ierr); call CHKERR(ierr)
        end if
        call PetscLogEventEnd(solver%logs%permute_mat_diff, ierr)

        call setup_ksp(solver, ksp, solver%C_diff, Aperm, prefix=prefix)

      else ! without permutation (normal case)

        call setup_ksp(solver, ksp, solver%C_diff, A, prefix=prefix)

      end if

      call PetscLogEventBegin(solver%logs%solve_Mdiff, ierr)
      call solve(solver, &
                 ksp, &
                 solver%b, &
                 solution%ediff, &
                 solution%uid, &
                 solution%Niter_diff, &
                 solution%diff_ksp_residual_history)
      call PetscLogEventEnd(solver%logs%solve_Mdiff, ierr)

      if (lmat_permute) then
        call PetscLogEventBegin(solver%logs%permute_mat_diff, ierr)
        call VecPermute(solver%b, solver%perm_diff%is, PETSC_TRUE, ierr); call CHKERR(ierr)
        call VecPermute(solution%ediff, solver%perm_diff%is, PETSC_TRUE, ierr); call CHKERR(ierr)
        call PetscLogEventEnd(solver%logs%permute_mat_diff, ierr)
      end if

    end subroutine
  end subroutine

  subroutine read_cmd_line_opts_get_coeffs( &
      & prefix, &
      & lgeometric_coeffs, &
      & ltop_bottom_faces_planar, &
      & ltop_bottom_planes_parallel &
      & )
    character(len=*), intent(in) :: prefix
    logical, intent(out) :: lgeometric_coeffs, ltop_bottom_faces_planar, ltop_bottom_planes_parallel
    logical :: lflg
    integer(mpiint) :: ierr

    lgeometric_coeffs = .false.
    call get_petsc_opt(prefix, "-pprts_geometric_coeffs", lgeometric_coeffs, lflg, ierr); call CHKERR(ierr)

    ltop_bottom_faces_planar = lgeometric_coeffs
    call get_petsc_opt(prefix, "-pprts_top_bottom_faces_planar", ltop_bottom_faces_planar, lflg, ierr); call CHKERR(ierr)

    ltop_bottom_planes_parallel = lgeometric_coeffs
    call get_petsc_opt(prefix, "-pprts_top_bottom_planes_parallel", ltop_bottom_planes_parallel, lflg, ierr); call CHKERR(ierr)
  end subroutine

  subroutine init_vertices( &
      & hs, &
      & dz, &
      & ltop_bottom_faces_planar, &
      & ltop_bottom_planes_parallel, &
      & vertices &
      & )
    logical, intent(in) :: ltop_bottom_faces_planar, ltop_bottom_planes_parallel
    real(ireals), intent(in) :: hs(:, :, :), dz
    real(ireals), intent(inout) :: vertices(:)

    vertices(3) = hs(2, 1, 1)
    vertices(6) = hs(2, 2, 1)
    vertices(9) = hs(2, 1, 2)
    vertices(12) = hs(2, 2, 2)
    vertices(15) = hs(1, 1, 1)
    vertices(18) = hs(1, 2, 1)
    vertices(21) = hs(1, 1, 2)
    vertices(24) = hs(1, 2, 2)

    if (ltop_bottom_faces_planar) then
      vertices(12) = vertices(9) + (vertices(6) - vertices(3))
      vertices(24) = vertices(21) + (vertices(18) - vertices(15))
    end if

    if (ltop_bottom_planes_parallel) then
      vertices(15:24:3) = vertices(3:12:3) + dz
    end if
  end subroutine

  !> @brief precompute transport coefficients
  subroutine alloc_coeff_dir2dir(solver, coeffs, opt_buildings)
    class(t_solver), intent(in) :: solver
    real(ireals), target, allocatable, intent(inout) :: coeffs(:, :, :, :)
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings
    real(irealLUT), allocatable :: v(:)
    integer(iintegers) :: src, k, i, j
    integer(mpiint) :: ierr

    real(ireals), pointer :: xhhl(:, :, :, :) => null(), xhhl1d(:) => null()
    real(ireals), allocatable :: vertices(:)
    real(ireals) :: norm
    real(ireals), pointer :: c(:, :)
    logical :: lgeometric_coeffs, ltop_bottom_faces_planar, ltop_bottom_planes_parallel

    associate ( &
        & atm => solver%atm, &
        & sun => solver%sun, &
        & C_dir => solver%C_dir)

      if (.not. allocated(coeffs)) &
        & allocate (coeffs(&
        & 1:C_dir%dof**2, &
        & C_dir%zs:C_dir%ze - 1, &
        & C_dir%xs:C_dir%xe, &
        & C_dir%ys:C_dir%ye))
      allocate (v(1:C_dir%dof**2))

      call PetscLogEventBegin(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)

      call read_cmd_line_opts_get_coeffs( &
        & solver%prefix, &
        & lgeometric_coeffs, &
        & ltop_bottom_faces_planar, &
        & ltop_bottom_planes_parallel &
        & )

      if (lgeometric_coeffs .and. C_dir%dof .ne. i3) &
        & call CHKERR(int(C_dir%dof, mpiint), 'geometric coeffs currently only implemented for 3 direct streams')
      call setup_default_unit_cube_geometry(atm%dx, atm%dy, -one, vertices)
      call getVecPointer(solver%Cvert_one_atm1%da, atm%vert_heights, xhhl1d, xhhl, readonly=.true.)

      do k = C_dir%zs, C_dir%ze - 1
        do j = C_dir%ys, C_dir%ye
          do i = C_dir%xs, C_dir%xe
            if (.not. atm%l1d(atmk(atm, k))) then
              call init_vertices( &
                & xhhl(i0, atmk(atm, k):atmk(atm, k + 1), i:i + 1, j:j + 1), &
                & atm%dz(atmk(atm, k), i, j), &
                & ltop_bottom_faces_planar, &
                & ltop_bottom_planes_parallel, &
                & vertices &
                & )

              if (lgeometric_coeffs) then
                vertices(3:24:3) = vertices(3:24:3) - minval(vertices(3:24:3))
                call dir2dir3_geometric_coeffs( &
                  & vertices, &
                  & sun%sundir, &
                  & atm%kabs(atmk(atm, k), i, j) + atm%ksca(atmk(atm, k), i, j), &
                  & coeffs(:, k, i, j) &
                  & )
              else
                call get_coeff( &
                  & solver%OPP, &
                  & atm%kabs(atmk(atm, k), i, j), &
                  & atm%ksca(atmk(atm, k), i, j), &
                  & atm%g(atmk(atm, k), i, j), &
                  & atm%dz(atmk(atm, k), i, j), &
                  & atm%dx, &
                  & .true., &
                  & v, &
                  & [real(sun%symmetry_phi, irealLUT), real(sun%theta, irealLUT)], &
                  & lswitch_east=sun%xinc .eq. 0, lswitch_north=sun%yinc .eq. 0, &
                  & opt_vertices=vertices &
                  & )
                coeffs(:, k, i, j) = real(v, ireals)
              end if

              if (ldebug_optprop) then
                c(1:C_dir%dof, 1:C_dir%dof) => coeffs(:, k, i, j) ! dim(src,dst)
                do src = 1, C_dir%dof
                  norm = sum(c(src, :))
                  if (norm .gt. one) then ! could renormalize
                    if (norm .gt. one + 100._ireals * sqrt(epsilon(norm))) then ! fatally off
                      print *, 'direct sum(dst==', src, ') gt one', norm
                      print *, 'direct coeff', norm, '::', c(src, :)
                      call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
                    end if ! fatally off
                    c(src, :) = c(src, :) / norm
                  end if ! could renormalize
                end do
              end if
            end if
          end do
        end do
      end do
      call restoreVecPointer(solver%Cvert_one_atm1%da, atm%vert_heights, xhhl1d, xhhl, readonly=.true.)

      call PetscLogEventEnd(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)
    end associate

    if (present(opt_buildings)) then
      call set_buildings_coeff()
    end if

  contains

    !> @brief   apply blocking of direct radiation from buildings
    !> @details Goal: set all src dof on a buildings face towards all dst dof to zero
    !> \n       albedo is not used in the dir2dir case, we only set blocking of radiation
    subroutine set_buildings_coeff()
      integer(iintegers) :: m, idx(4)

      associate ( &
          & C => solver%C_dir, &
          & B => opt_buildings)
        do m = 1, size(B%iface)
          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) - 1 + [C%zs, C%xs, C%ys]

          associate (k => idx(2), i => idx(3), j => idx(4))
            coeffs(:, k, i, j) = -zero
          end associate
        end do
      end associate
    end subroutine

  end subroutine

  subroutine alloc_coeff_dir2diff(solver, coeffs)
    class(t_solver), intent(in), target :: solver
    real(ireals), target, allocatable, intent(inout) :: coeffs(:, :, :, :)

    real(irealLUT), allocatable :: v(:), T_LUT(:)
    real(ireals), allocatable :: T_GOMTRC(:), S_GOMTRC(:), S_LUT(:)
    integer(iintegers) :: src, k, i, j
    integer(mpiint) :: ierr

    real(ireals), pointer :: xhhl(:, :, :, :) => null(), xhhl1d(:) => null()
    real(ireals), allocatable :: vertices(:)
    real(ireals) :: norm_diff, norm_dir, normref
    real(ireals) :: S_LUT_norm, T_LUT_norm, T_GOMTRC_norm
    real(ireals), pointer :: c(:, :), cdir2dir(:, :)
    real(ireals), parameter :: eps = one + sqrt(sqrt(epsilon(eps)))

    logical :: &
      & lbmc_online,                 &
      & lcheck_coeff_sums,           &
      & lconserve_lut_atm_abso,      &
      & lflg,                        &
      & lgeometric_coeffs,           &
      & ltop_bottom_faces_planar,    &
      & ltop_bottom_planes_parallel

    associate ( &
        & atm => solver%atm, &
        & sun => solver%sun, &
        & C_dir => solver%C_dir, &
        & C_diff => solver%C_diff)

      if (.not. allocated(coeffs)) &
        & allocate (coeffs(1:C_dir%dof * C_diff%dof, &
                         & C_dir%zs:C_dir%ze - 1,    &
                         & C_dir%xs:C_dir%xe,      &
                         & C_dir%ys:C_dir%ye))
      allocate (v(1:C_dir%dof * C_diff%dof))
      allocate (S_GOMTRC(1:C_dir%dof * C_diff%dof))
      allocate (T_LUT(1:C_dir%dof**2))
      allocate (S_LUT(1:C_dir%dof * C_diff%dof))
      allocate (T_GOMTRC(1:C_dir%dof**2))
      call PetscLogEventBegin(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)

      lbmc_online = .false.
      call get_petsc_opt(solver%prefix, "-bmc_online", lbmc_online, lflg, ierr); call CHKERR(ierr)

      lcheck_coeff_sums = .not. lbmc_online
      call get_petsc_opt(solver%prefix, "-pprts_check_coeff_sums", lcheck_coeff_sums, lflg, ierr); call CHKERR(ierr)

      call read_cmd_line_opts_get_coeffs( &
        & solver%prefix, &
        & lgeometric_coeffs, &
        & ltop_bottom_faces_planar, &
        & ltop_bottom_planes_parallel &
        & )

      lconserve_lut_atm_abso = lgeometric_coeffs
      call get_petsc_opt(solver%prefix, "-pprts_conserve_lut_atm_abso", lconserve_lut_atm_abso, lflg, ierr); call CHKERR(ierr)

      if (lgeometric_coeffs .and. lbmc_online .and. lconserve_lut_atm_abso) then
        lbmc_online = .false.
      end if
      call setup_default_unit_cube_geometry(atm%dx, atm%dy, -one, vertices)
      call getVecPointer(solver%Cvert_one_atm1%da, solver%atm%vert_heights, xhhl1d, xhhl, readonly=.true.)

      do k = C_dir%zs, C_dir%ze - 1
        do j = C_dir%ys, C_dir%ye
          do i = C_dir%xs, C_dir%xe
            if (.not. atm%l1d(atmk(atm, k))) then
              call init_vertices( &
                & xhhl(i0, atmk(atm, k):atmk(atm, k + 1), i:i + 1, j:j + 1), &
                & atm%dz(atmk(atm, k), i, j), &
                & ltop_bottom_faces_planar, &
                & ltop_bottom_planes_parallel, &
                & vertices &
                & )

              if (lbmc_online) then
                call get_coeff( &
                  & solver%OPP, &
                  & atm%kabs(atmk(atm, k), i, j), &
                  & atm%ksca(atmk(atm, k), i, j), &
                  & atm%g(atmk(atm, k), i, j), &
                  & atm%dz(atmk(atm, k), i, j), &
                  & atm%dx, &
                  & .false., &
                  & v, &
                  & [real(sun%symmetry_phi, irealLUT), real(sun%theta, irealLUT)], &
                  & lswitch_east=sun%xinc .eq. 0, lswitch_north=sun%yinc .eq. 0, &
                  & opt_vertices=vertices &
                  & )
              else
                call get_coeff( &
                  & solver%OPP, &
                  & atm%kabs(atmk(atm, k), i, j), &
                  & atm%ksca(atmk(atm, k), i, j), &
                  & atm%g(atmk(atm, k), i, j), &
                  & atm%dz(atmk(atm, k), i, j), &
                  & atm%dx, &
                  & .false., &
                  & v, &
                  & [real(sun%symmetry_phi, irealLUT), real(sun%theta, irealLUT)], &
                  & lswitch_east=sun%xinc .eq. 0, lswitch_north=sun%yinc .eq. 0 &
                  & )
              end if
              S_LUT = real(v, ireals)
              coeffs(:, k, i, j) = S_LUT

              if (lgeometric_coeffs .and. lconserve_lut_atm_abso) then
                call get_coeff( &
                  & solver%OPP, &
                  & atm%kabs(atmk(atm, k), i, j), &
                  & atm%ksca(atmk(atm, k), i, j), &
                  & atm%g(atmk(atm, k), i, j), &
                  & atm%dz(atmk(atm, k), i, j), &
                  & atm%dx, &
                  & .true., &
                  & T_LUT, &
                  & [real(sun%symmetry_phi, irealLUT), real(sun%theta, irealLUT)], &
                  & lswitch_east=sun%xinc .eq. 0, lswitch_north=sun%yinc .eq. 0 &
                  & )

                T_GOMTRC = solver%dir2dir(:, k, i, j)

                do src = 1, C_dir%dof
                  S_LUT_norm = sum(S_LUT(src:C_dir%dof * C_diff%dof:C_dir%dof))
                  T_LUT_norm = sum(T_LUT(src:C_dir%dof**2:C_dir%dof))
                  T_GOMTRC_norm = sum(T_GOMTRC(src:C_dir%dof**2:C_dir%dof))

                  if (S_LUT_norm .le. epsilon(S_LUT_norm)) then
                    S_GOMTRC(src:C_dir%dof * C_diff%dof:C_dir%dof) = zero
                  else
                    S_GOMTRC(src:C_dir%dof * C_diff%dof:C_dir%dof) = S_LUT(src:C_dir%dof * C_diff%dof:C_dir%dof) / S_LUT_norm
                  end if
                  S_GOMTRC(src:C_dir%dof * C_diff%dof:C_dir%dof) = S_GOMTRC(src:C_dir%dof * C_diff%dof:C_dir%dof) * &
                                                                   (one - (one - T_LUT_norm - S_LUT_norm) - T_GOMTRC_norm)
                end do

                coeffs(:, k, i, j) = S_GOMTRC
              end if ! lgeometric_coeffs .and. lconserve_lut_atm_abso

              if (ldebug_optprop) then
                c(1:C_dir%dof, 1:C_diff%dof) => coeffs(:, k, i, j) ! dim(src,dst)
                cdir2dir(1:C_dir%dof, 1:C_dir%dof) => solver%dir2dir(:, k, i, j)
                do src = 1, C_dir%dof
                  norm_diff = sum(c(src, :))
                  norm_dir = sum(cdir2dir(src, :))
                  if (lcheck_coeff_sums) then
                    normref = norm_dir + norm_diff
                    if (normref .gt. eps) &
                      call CHKERR(1, 'Failed since for src'//toStr(src)//new_line('A')//&
                        & ', norm(dir2dir(src)) + norm(dir2diff(src)) = '//new_line('A')//&
                        & toStr(norm_dir)//' + '//toStr(norm_diff)//&
                        & ' = '//toStr(normref)//' > '//toStr(eps)//new_line('')//&
                        & 'S = '//toStr(c(src, :))//new_line('A')//&
                        & 'T = '//toStr(cdir2dir(src, :))//new_line('A')//&
                        & '; box indize: k = '//toStr(k)//', i = '//toStr(i)//', j = '//toStr(j)//new_line('A')//&
                        & '; vertices: '//toStr(vertices)//new_line('A')//&
                        & '; c_abso = '//toStr(atm%kabs(atmk(solver%atm, k), i, j))//new_line('A')//&
                        & '; c_scat = '//toStr(atm%ksca(atmk(solver%atm, k), i, j))//&
                        & '; g = '//toStr(atm%g(atmk(solver%atm, k), i, j)))
                  end if
                  if (norm_diff .gt. one) then ! could renormalize
                    if (norm_diff .gt. eps) then ! fatally off
                      print *, 'dir2diff sum(dst==', src, ') gt one', norm_diff
                      print *, 'dir2diff coeff', norm_diff, '::', c(src, :)
                      call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
                    end if ! fatally off
                    c(src, :) = c(src, :) / norm_diff
                  end if ! could renormalize
                end do
              end if
            end if
          end do
        end do
      end do
      call restoreVecPointer(solver%Cvert_one_atm1%da, solver%atm%vert_heights, xhhl1d, xhhl, readonly=.true.)
      call PetscLogEventEnd(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)
    end associate
  end subroutine

  subroutine alloc_coeff_diff2diff(solver, coeffs, opt_buildings)
    class(t_solver), intent(in) :: solver
    real(ireals), target, allocatable, intent(inout) :: coeffs(:, :, :, :)
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings
    real(ireals), allocatable :: vertices(:)
    real(irealLUT), allocatable :: v(:)
    integer(iintegers) :: src, k, i, j
    integer(mpiint) :: ierr

    real(ireals) :: norm
    real(ireals), pointer :: c(:, :)
    real(ireals), pointer :: xhhl(:, :, :, :) => null(), xhhl1d(:) => null()

    logical :: lgeometric_coeffs, ltop_bottom_faces_planar, ltop_bottom_planes_parallel

    associate ( &
        & atm => solver%atm, &
        & C_diff => solver%C_diff)
      if (.not. allocated(coeffs)) &
        & allocate (coeffs(         &
          & 1:C_diff%dof**2,       &
          & C_diff%zs:C_diff%ze - 1, &
          & C_diff%xs:C_diff%xe,   &
          & C_diff%ys:C_diff%ye))
      allocate (v(1:C_diff%dof**2))
      call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)

      call read_cmd_line_opts_get_coeffs( &
        & solver%prefix, &
        & lgeometric_coeffs, &
        & ltop_bottom_faces_planar, &
        & ltop_bottom_planes_parallel &
        & )

      call setup_default_unit_cube_geometry(atm%dx, atm%dy, -one, vertices)
      call getVecPointer(solver%Cvert_one_atm1%da, atm%vert_heights, xhhl1d, xhhl, readonly=.true.)

      do k = C_diff%zs, C_diff%ze - 1
        do j = C_diff%ys, C_diff%ye
          do i = C_diff%xs, C_diff%xe
            if (.not. atm%l1d(atmk(atm, k))) then
              call init_vertices( &
                & xhhl(i0, atmk(atm, k):atmk(atm, k + 1), i:i + 1, j:j + 1), &
                & atm%dz(atmk(atm, k), i, j), &
                & ltop_bottom_faces_planar, &
                & ltop_bottom_planes_parallel, &
                & vertices &
                & )

              call get_coeff( &
                & solver%OPP, &
                & atm%kabs(atmk(atm, k), i, j), &
                & atm%ksca(atmk(atm, k), i, j), &
                & atm%g(atmk(atm, k), i, j), &
                & atm%dz(atmk(atm, k), i, j), &
                & atm%dx, &
                & .false., &
                & v, &
                & opt_vertices=vertices &
                & )
              coeffs(:, k, i, j) = real(v, ireals)

              if (ldebug_optprop) then
                c(1:C_diff%dof, 1:C_diff%dof) => coeffs(:, k, i, j) ! dim(src,dst)
                do src = 1, C_diff%dof
                  norm = sum(c(src, :))
                  if (norm .gt. one) then ! could renormalize
                    if (norm .gt. one + 100._ireals * sqrt(epsilon(norm))) then ! fatally off
                      print *, 'diffuse sum(dst==', src, ') gt one', norm
                      print *, 'diffuse coeff', norm, '::', c(src, :)
                      call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
                    end if ! fatally off
                    c(src, :) = c(src, :) / norm
                  end if ! could renormalize
                end do
              end if
            end if
          end do
        end do
      end do
      call restoreVecPointer(solver%Cvert_one_atm1%da, atm%vert_heights, xhhl1d, xhhl, readonly=.true.)
      call PetscLogEventEnd(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
    end associate

    if (present(opt_buildings)) then
      call set_buildings_coeff(coeffs)
    end if
  contains

    !> @brief   apply blocking of diffuse radiation from buildings and do lambertian reflections
    !> @details Goal: first set all src dof on a buildings face towards all dst dof to zero
    !> \n       note, this assumes that the building is the full cell,
    !> \n       i.e. the albedo is applied on the outside of the cell
    subroutine set_buildings_coeff(coeffs)
      real(ireals), target, intent(inout) :: coeffs(:, :, :, :)
      integer(iintegers) :: m, idx(4)
      integer(iintegers) :: isrc, idst, off
      real(ireals), pointer :: v(:, :) ! dim(src, dst)

      associate ( &
          & C => solver%C_diff, &
          & B => opt_buildings)
        do m = 1, size(B%iface)
          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)

          associate (k => idx(2), i => idx(3), j => idx(4))
            v(0:C%dof - 1, 0:C%dof - 1) => coeffs(1:C%dof**2, k, i, j)
            select case (idx(1))
            case (PPRTS_TOP_FACE)
              do idst = 0, solver%difftop%dof - 1
                if (.not. solver%difftop%is_inward(i1 + idst)) then ! eup
                  v(:, idst) = 0._ireals
                  do isrc = 0, solver%difftop%dof - 1
                    if (solver%difftop%is_inward(i1 + isrc)) then ! edn
                      v(isrc, idst) = real(B%albedo(m), ireals) / real(solver%difftop%streams, ireals)
                    end if
                  end do
                end if
              end do

            case (PPRTS_BOT_FACE)
              do idst = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + idst)) then ! edn
                  v(:, idst) = 0._ireals
                  do isrc = 0, solver%difftop%dof - 1
                    if (.not. solver%difftop%is_inward(i1 + isrc)) then ! eup
                      v(isrc, idst) = real(B%albedo(m), ireals) / real(solver%difftop%streams, ireals)
                    end if
                  end do
                end if
              end do

            case (PPRTS_LEFT_FACE)
              off = solver%difftop%dof
              do idst = 0, solver%diffside%dof - 1
                if (.not. solver%diffside%is_inward(i1 + idst)) then ! leftward
                  v(:, off + idst) = 0._ireals
                  do isrc = 0, solver%diffside%dof - 1
                    if (solver%diffside%is_inward(i1 + isrc)) then ! right ward
                      v(off + isrc, off + idst) = real(B%albedo(m), ireals) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end if
              end do

            case (PPRTS_RIGHT_FACE)
              off = solver%difftop%dof
              do idst = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idst)) then ! rightward
                  v(:, off + idst) = 0._ireals
                  do isrc = 0, solver%diffside%dof - 1
                    if (.not. solver%diffside%is_inward(i1 + isrc)) then ! leftward
                      v(off + isrc, off + idst) = real(B%albedo(m), ireals) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end if
              end do

            case (PPRTS_REAR_FACE)
              off = solver%difftop%dof + solver%diffside%dof
              do idst = 0, solver%diffside%dof - 1
                if (.not. solver%diffside%is_inward(i1 + idst)) then ! backward
                  v(:, off + idst) = 0._ireals
                  do isrc = 0, solver%diffside%dof - 1
                    if (solver%diffside%is_inward(i1 + isrc)) then ! forward
                      v(off + isrc, off + idst) = real(B%albedo(m), ireals) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end if
              end do

            case (PPRTS_FRONT_FACE)
              off = solver%difftop%dof + solver%diffside%dof
              do idst = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idst)) then ! forward
                  v(:, off + idst) = 0._ireals
                  do isrc = 0, solver%diffside%dof - 1
                    if (.not. solver%diffside%is_inward(i1 + isrc)) then ! backward
                      v(off + isrc, off + idst) = real(B%albedo(m), ireals) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end if
              end do
            end select
          end associate
        end do
      end associate
    end subroutine

  end subroutine

  !> @brief renormalize fluxes with the size of a face(sides or lid)
  subroutine scale_flx(solver, solution, lWm2)
    class(t_solver), intent(inout) :: solver
    type(t_state_container), intent(inout) :: solution   !< @param solution container with computed fluxes
    logical, intent(in) :: lWm2  !< @param determines direction of scaling, if true, scale to W/m**2
    integer(mpiint) :: ierr

    call PetscLogEventBegin(solver%logs%scale_flx, ierr); call CHKERR(ierr)
    if (solution%lsolar_rad) then
      if (.not. allocated(solver%dir_scalevec_Wm2_to_W)) then
        allocate (solver%dir_scalevec_Wm2_to_W)
        call VecDuplicate(solution%edir, solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        call gen_scale_dir_flx_vec(solver, solver%dir_scalevec_Wm2_to_W, solver%C_dir)
      end if
      if (.not. allocated(solver%dir_scalevec_W_to_Wm2)) then
        allocate (solver%dir_scalevec_W_to_Wm2)
        call VecDuplicate(solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        call VecSet(solver%dir_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
        call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
          solver%dir_scalevec_W_to_Wm2, &
          solver%dir_scalevec_W_to_Wm2, &
          solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      end if
      if (solution%lWm2_dir .neqv. lWm2) then
        if (lWm2) then
          call VecPointwiseMult(solution%edir, solution%edir, solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        else
          call VecPointwiseMult(solution%edir, solution%edir, solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        end if
        solution%lWm2_dir = lWm2
      end if
    end if

    if (.not. allocated(solver%diff_scalevec_Wm2_to_W)) then
      allocate (solver%diff_scalevec_Wm2_to_W)
      call VecDuplicate(solution%ediff, solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      call gen_scale_diff_flx_vec(solver, solver%diff_scalevec_Wm2_to_W, solver%C_diff)
    end if
    if (.not. allocated(solver%diff_scalevec_W_to_Wm2)) then
      allocate (solver%diff_scalevec_W_to_Wm2)
      call VecDuplicate(solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
      call VecSet(solver%diff_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
      call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
        solver%diff_scalevec_W_to_Wm2, &
        solver%diff_scalevec_W_to_Wm2, &
        solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
    end if

    if (solution%lWm2_diff .neqv. lWm2) then
      if (lWm2) then
        call VecPointwiseMult(solution%ediff, solution%ediff, solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
      else
        call VecPointwiseMult(solution%ediff, solution%ediff, solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      end if
      solution%lWm2_diff = lWm2
    end if
    call PetscLogEventEnd(solver%logs%scale_flx, ierr); call CHKERR(ierr)

  contains
    subroutine gen_scale_dir_flx_vec(solver, v, coord)
      class(t_solver) :: solver
      type(tVec) :: v
      type(t_coord) :: coord
      real(ireals), pointer :: xv(:, :, :, :) => null()
      real(ireals), pointer :: xv1d(:) => null()
      integer(iintegers) :: i, j, k, l, iside, ak
      real(ireals) :: Ax, Ay, Az, fac

      real(ireals), allocatable :: vertices(:)
      real(ireals), pointer :: xhhl(:, :, :, :) => null(), xhhl1d(:) => null()

      call getVecPointer(solver%Cvert_one_atm1%da, solver%atm%vert_heights, xhhl1d, xhhl, readonly=.true.)
      call setup_default_unit_cube_geometry(solver%atm%dx, solver%atm%dy, one, vertices)

      associate ( &
         & atm => solver%atm,    &
         & A => vertices(1:3), &
         & B => vertices(4:6), &
         & C => vertices(7:9), &
         & D => vertices(10:12), &
         & E => vertices(13:15), &
         & F => vertices(16:18), &
         & G => vertices(19:21))
        !& H => vertices(22:24)  )

        if (solver%myid .eq. 0 .and. ldebug) print *, 'rescaling direct fluxes', coord%zm, coord%xm, coord%ym
        call getVecPointer(coord%da, v, xv1d, xv)

        ! Scaling top faces
        do j = coord%ys, coord%ye
          do i = coord%xs, coord%xe
            do k = coord%zs, coord%ze
              ak = atmk(atm, k)
              A(3) = xhhl(i0, ak, i, j)
              B(3) = xhhl(i0, ak, i + 1, j)
              C(3) = xhhl(i0, ak, i, j + 1)
              D(3) = xhhl(i0, ak, i + 1, j + 1)

              Az = triangle_area_by_vertices(A, B, D) + triangle_area_by_vertices(A, D, C)
              fac = Az / real(solver%dirtop%area_divider, ireals)
              do iside = 0, solver%dirtop%dof - 1
                xv(iside, k, i, j) = fac
              end do
            end do
          end do
        end do

        ! Scaling side faces
        do j = coord%ys, coord%ye
          do i = coord%xs, coord%xe
            do k = coord%zs, coord%ze - 1
              ak = atmk(atm, k)

              ! First the faces in x-direction
              A(3) = xhhl(i0, ak + 1, i, j)
              B(3) = xhhl(i0, ak + 1, i + 1, j)
              C(3) = xhhl(i0, ak + 1, i, j + 1)
              E(3) = xhhl(i0, ak, i, j)
              F(3) = xhhl(i0, ak, i + 1, j)
              G(3) = xhhl(i0, ak, i, j + 1)

              !Ax = solver%atm%dy*solver%atm%dz(ak,i,j)
              Ax = triangle_area_by_vertices(A, B, F) + triangle_area_by_vertices(A, F, E)
              fac = Ax / real(solver%dirside%area_divider, ireals)
              do iside = 0, solver%dirside%dof - 1
                l = solver%dirtop%dof + iside
                xv(l, k, i, j) = fac
              end do

              ! Then the rest of the faces in y-direction
              !Ay = atm%dy*atm%dz(ak,i,j) / real(solver%dirside%area_divider, ireals)
              Ay = triangle_area_by_vertices(A, C, G) + triangle_area_by_vertices(A, G, E)
              fac = Ay / real(solver%dirside%area_divider, ireals)
              do iside = 0, solver%dirside%dof - 1
                l = solver%dirtop%dof + solver%dirside%dof + iside
                xv(l, k, i, j) = fac
              end do
            end do
            ! the side faces underneath the surface are always scaled by unity
            xv(solver%dirtop%dof:ubound(xv, dim=1), coord%ze, i, j) = one
          end do
        end do
        call restoreVecPointer(coord%da, v, xv1d, xv)
      end associate

      call restoreVecPointer(solver%Cvert_one_atm1%da, solver%atm%vert_heights, xhhl1d, xhhl, readonly=.true.)
    end subroutine
    subroutine gen_scale_diff_flx_vec(solver, v, C)
      class(t_solver) :: solver
      type(tVec) :: v
      type(t_coord) :: C
      real(ireals), pointer :: xv(:, :, :, :) => null()
      real(ireals), pointer :: xv1d(:) => null()

      integer(iintegers) :: iside, src, i, j, k, ak
      real(ireals) :: Az, Ax, Ay, fac

      if (solver%myid .eq. 0 .and. ldebug) print *, 'rescaling fluxes', C%zm, C%xm, C%ym
      call getVecPointer(C%da, v, xv1d, xv)

      if (C%dof .eq. i3 .or. C%dof .eq. i8) then
        print *, 'scale_flx_vec is just for diffuse radia'
      end if

      ! Scaling top faces
      Az = solver%atm%dx * solver%atm%dy / real(solver%difftop%area_divider, ireals)
      fac = Az

      do j = C%ys, C%ye
        do i = C%xs, C%xe
          do k = C%zs, C%ze
            do iside = 1, solver%difftop%dof
              src = iside - 1
              xv(src, k, i, j) = fac                  ! diffuse radiation
            end do
          end do
        end do
      end do

      ! Scaling side faces
      do j = C%ys, C%ye
        do i = C%xs, C%xe
          do k = C%zs, C%ze - 1
            ak = atmk(solver%atm, k)
            ! faces in x-direction
            Ax = solver%atm%dy * solver%atm%dz(ak, i, j) / real(solver%diffside%area_divider, ireals)
            fac = Ax

            do iside = 1, solver%diffside%dof
              src = solver%difftop%dof + iside - 1
              xv(src, k, i, j) = fac
            end do

            ! faces in y-direction
            Ay = solver%atm%dx * solver%atm%dz(ak, i, j) / real(solver%difftop%area_divider, ireals)
            fac = Ay

            do iside = 1, solver%diffside%dof
              src = solver%difftop%dof + solver%diffside%dof + iside - 1
              xv(src, k, i, j) = fac
            end do
          end do
          ! the side faces underneath the surface are always scaled by unity
          xv(solver%difftop%dof:ubound(xv, dim=1), k, i, j) = one
        end do
      end do
      call restoreVecPointer(C%da, v, xv1d, xv)

    end subroutine
  end subroutine

  !> @brief get solution to a final state, e.g. update absorption
  subroutine restore_solution(solver, solution, time)
    ! restore_solution:: if flux have changed, we need to update absorption, save the residual history
    class(t_solver) :: solver
    type(t_state_container) :: solution
    real(ireals), intent(in), optional :: time

    character(default_str_len) :: vecname
    real(ireals) :: inf_norm
    type(tVec) :: abso_old
    integer(mpiint) :: ierr

    if (.not. solution%lset) call CHKERR(1_mpiint, 'cant restore solution that was not initialized')
    if (.not. allocated(solution%abso)) call CHKERR(1_mpiint, 'cant restore solution that was not initialized')

    if (present(time)) then
      solution%time = eoshift(solution%time, shift=-1) !shift all values by 1 to the right
      solution%time(1) = time
      if (ldebug .and. solver%myid .eq. 0) print *, 'Setting solution time '//toStr(time)//'('//toStr(solution%time)//')'
    end if

    if (.not. solution%lchanged) return

    if (present(time) .and. solver%lenable_solutions_err_estimates) then ! Create working vec to determine difference between old and new absorption vec
      call DMGetGlobalVector(solver%C_one%da, abso_old, ierr); call CHKERR(ierr)
      call VecCopy(solution%abso, abso_old, ierr); call CHKERR(ierr)
    end if

    ! update absorption
    call calc_flx_div(solver, solution)

    if (ldebug .and. solver%myid .eq. 0) &
      print *, 'Saving Solution ', solution%uid

    ! make sure to bring the fluxes into [W/m**2]
    call scale_flx(solver, solution, lWm2=.true.)

    if (ldebug .and. solver%myid .eq. 0) &
      print *, 'Saving Solution done'
    solution%lchanged = .false.

    if (present(time) .and. solver%lenable_solutions_err_estimates) then ! Compute norm between old absorption and new one
      call VecAXPY(abso_old, -one, solution%abso, ierr); call CHKERR(ierr) ! overwrite abso_old with difference to new one
      call VecNorm(abso_old, NORM_INFINITY, inf_norm, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(solver%C_one%da, abso_old, ierr); call CHKERR(ierr)

      ! Save norm for later analysis
      solution%maxnorm = eoshift(solution%maxnorm, shift=-1) !shift all values by 1 to the right
      solution%maxnorm(1) = inf_norm

      if (ldebug .and. solver%myid .eq. 0) &
        print *, 'Updating error statistics for solution ', solution%uid, 'at time ', time, '::', solution%time(1), &
        ':: norm', inf_norm, '[W] :: hr_norm approx:', inf_norm * 86.1, '[K/d]'
    end if !present(time) .and. solver%lenable_solutions_err_estimates

    if (allocated(solution%edir)) then
      write (vecname, FMT='("edir",I0)') solution%uid
      call PetscObjectSetName(solution%edir, vecname, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_edir", ierr); call CHKERR(ierr)
    end if

    if (allocated(solution%ediff)) then
      write (vecname, FMT='("ediff",I0)') solution%uid
      call PetscObjectSetName(solution%ediff, vecname, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-show_ediff", ierr); call CHKERR(ierr)
    end if

    if (allocated(solution%abso)) then
      write (vecname, FMT='("abso",I0)') solution%uid
      call PetscObjectSetName(solution%abso, vecname, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, "-show_abso", ierr); call CHKERR(ierr)
    end if

    if (allocated(solver%b)) then
      write (vecname, FMT='("b",I0)') solution%uid
      call PetscObjectSetName(solver%b, vecname, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solver%b, PETSC_NULL_VEC, "-show_b", ierr); call CHKERR(ierr)
    end if

    if (allocated(solver%incSolar)) then
      write (vecname, FMT='("incSolar",I0)') solution%uid
      call PetscObjectSetName(solver%incSolar, vecname, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solver%incSolar, PETSC_NULL_VEC, "-show_incSolar", ierr); call CHKERR(ierr)
    end if
  end subroutine

  !> @brief Compute flux divergence, i.e. absorption
  !> @details from gauss divergence theorem, the divergence in the volume is the integral of the flux through the surface
  !> \n we therefore sum up the incoming and outgoing fluxes to compute the divergence
  subroutine calc_flx_div(solver, solution)
    class(t_solver), target :: solver
    type(t_state_container) :: solution

    real(ireals), pointer, dimension(:, :, :, :) :: xediff => null(), xedir => null(), xabso => null()
    real(ireals), pointer, dimension(:) :: xediff1d => null(), xedir1d => null(), xabso1d => null()

    integer(iintegers) :: isrc, src
    integer(iintegers) :: i, j, k, xinc, yinc
    type(tVec) :: ledir, lediff ! local copies of vectors, including ghosts
    real(ireals) :: Volume, Az

    logical :: by_coeff_divergence, ldirect_absorption_only, lflg
    integer(mpiint) :: ierr

    if (solver%myid .eq. 0 .and. ldebug) print *, 'Calculating flux divergence solar?', solution%lsolar_rad, 'NCA?', lcalc_nca

    if (allocated(solution%edir)) then
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_flxdiv_edir", ierr); call CHKERR(ierr)
    end if
    call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-show_flxdiv_ediff", ierr); call CHKERR(ierr)

    if ((solution%lsolar_rad .eqv. .false.) .and. lcalc_nca) then ! if we should calculate NCA (Klinger), we can just return afterwards
      call scale_flx(solver, solution, lWm2=.true.)
      call nca_wrapper(solver, solution%ediff, solution%abso)
      return
    end if

    call PetscLogEventBegin(solver%logs%compute_absorption, ierr); call CHKERR(ierr)

    if (.not. allocated(solver%abso_scalevec)) &
      & call gen_abso_scalevec()

    call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)

    ! make sure to bring the fluxes into [W] before the absorption calculation
    call scale_flx(solver, solution, lWm2=.false.)

    ldirect_absorption_only = .false.
    call get_petsc_opt('', "-direct_absorption_only", ldirect_absorption_only, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(solver%prefix, "-direct_absorption_only", ldirect_absorption_only, lflg, ierr); call CHKERR(ierr)
    if (ldirect_absorption_only) then
      call direct_absorption_only()
      goto 99
    end if

    ! if there are no 3D layers globally, we should skip the ghost value copying....
    !lhave_no_3d_layer = mpi_logical_and(solver%comm, all(atm%l1d.eqv..True.))
    select type (solver)
    class is (t_solver_2str)
      call compute_1D_absorption()
      goto 99
    end select

    by_coeff_divergence = .true.
    call get_petsc_opt(solver%prefix, "-absorption_by_coeff_divergence", by_coeff_divergence, lflg, ierr); call CHKERR(ierr)

    if (by_coeff_divergence) then
      call compute_absorption_by_coeff_divergence()
    else
      call compute_absorption_by_flx_divergence()
    end if

99  continue ! cleanup
    call VecPointwiseMult(solution%abso, solution%abso, solver%abso_scalevec, ierr); call CHKERR(ierr)

    call PetscLogEventEnd(solver%logs%compute_absorption, ierr); call CHKERR(ierr)
  contains

    subroutine direct_absorption_only()
      real(ireals) :: dtau
      real(irealLUT), target :: dir2dir(solver%C_dir%dof**2)
      real(irealLUT), pointer :: pdir2dir(:, :) ! dim(src,dst)
      integer(iintegers) :: ak
      associate ( &
        sun => solver%sun, &
        atm => solver%atm, &
        C_dir => solver%C_dir, &
        C_one => solver%C_one)

        pdir2dir(0:solver%C_dir%dof - 1, 0:solver%C_dir%dof - 1) => dir2dir

        if (solution%lsolar_rad) then
          ! Copy ghosted values for direct vec
          call DMGetLocalVector(C_dir%da, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C_dir%da, solution%edir, INSERT_VALUES, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd(C_dir%da, solution%edir, INSERT_VALUES, ledir, ierr); call CHKERR(ierr)
          call getVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.true.)
        end if

        call getVecPointer(C_one%da, solution%abso, xabso1d, xabso)

        if (solution%lsolar_rad) then
          do j = C_one%ys, C_one%ye
            do i = C_one%xs, C_one%xe
              do k = C_one%zs, C_one%ze
                ak = atmk(atm, k)

                if (atm%l1d(atmk(atm, k))) then ! one dimensional i.e. twostream
                  do isrc = 0, solver%dirtop%dof - 1
                    dtau = atm%kabs(ak, i, j) * atm%dz(ak, i, j) / sun%costheta
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(isrc, k, i, j) * (one - exp(-dtau))
                  end do

                else ! 3D-radiation
                  ! direct part of absorption
                  xinc = sun%xinc
                  yinc = sun%yinc

                  call get_coeff( &
                    & solver%OPP, &
                    & atm%kabs(ak, i, j), &
                    & zero, &
                    & atm%g(ak, i, j), &
                    & atm%dz(ak, i, j), &
                    & atm%dx, &
                    & .true., &
                    & dir2dir, &
                    & [real(irealLUT) :: sun%symmetry_phi, sun%theta], &
                    & lswitch_east=xinc .eq. 0, lswitch_north=yinc .eq. 0 &
                    & )

                  src = 0
                  do isrc = 0, solver%dirtop%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(src, k, i, j) * (one - sum(pdir2dir(src, :)))
                    src = src + 1
                  end do

                  do isrc = 0, solver%dirside%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(src, k, i + 1 - xinc, j) * (one - sum(pdir2dir(src, :)))
                    src = src + 1
                  end do

                  do isrc = 0, solver%dirside%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(src, k, i, j + i1 - yinc) * (one - sum(pdir2dir(src, :)))
                    src = src + 1
                  end do
                end if ! 1d/3D
              end do
            end do
          end do
        end if

        if (solution%lsolar_rad) then
          call restoreVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.true.)
          call DMRestoreLocalVector(C_dir%da, ledir, ierr); call CHKERR(ierr)
        end if

        call restoreVecPointer(C_one%da, solution%abso, xabso1d, xabso)
      end associate
    end subroutine

    subroutine compute_1D_absorption()
      associate ( &
        C_dir => solver%C_dir, &
        C_diff => solver%C_diff, &
        C_one => solver%C_one)
        if (ldebug .and. solver%myid .eq. 0) print *, 'lhave_no_3d_layer => will use 1D absorption computation'

        if (solution%lsolar_rad) call getVecPointer(C_dir%da, solution%edir, xedir1d, xedir)

        call getVecPointer(C_diff%da, solution%ediff, xediff1d, xediff, readonly=.true.)
        call getVecPointer(C_one%da, solution%abso, xabso1d, xabso)

        !! calculate absorption by flux divergence

        do j = C_one%ys, C_one%ye
          do i = C_one%xs, C_one%xe
            do k = C_one%zs, C_one%ze
              if (solution%lsolar_rad) then
                do isrc = i0, solver%dirtop%dof - 1
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(isrc, k, i, j) - xedir(isrc, k + i1, i, j)
                end do
              end if

              do isrc = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + isrc)) then
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(isrc, k, i, j) - xediff(isrc, k + 1, i, j))
                else
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(isrc, k + 1, i, j) - xediff(isrc, k, i, j))
                end if
              end do
            end do
          end do
        end do

        if (solution%lsolar_rad) call restoreVecPointer(C_dir%da, solution%edir, xedir1d, xedir)

        call restoreVecPointer(C_diff%da, solution%ediff, xediff1d, xediff, readonly=.true.)
        call restoreVecPointer(C_one%da, solution%abso, xabso1d, xabso)
      end associate
    end subroutine

    subroutine compute_absorption_by_coeff_divergence()
      real(ireals) :: cdiv
      real(ireals), pointer :: dir2dir(:, :) ! dim(src, dst)
      real(ireals), pointer :: dir2diff(:, :) ! dim(src, dst)
      real(ireals), pointer :: diff2diff(:, :) ! dim(src, dst)

      type(tVec) :: local_b
      real(ireals), pointer, dimension(:, :, :, :) :: xsrc => null()
      real(ireals), pointer, dimension(:) :: xsrc1d => null()

      integer(iintegers) :: idof, msrc

      call getVecPointer(solver%C_one%da, solution%abso, xabso1d, xabso)

      associate ( &
        sun => solver%sun, &
        atm => solver%atm, &
        C_dir => solver%C_dir, &
        C_diff => solver%C_diff, &
        C_one => solver%C_one)

        if (solution%lsolar_rad) then
          ! Copy ghosted values for direct vec
          call DMGetLocalVector(C_dir%da, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C_dir%da, solution%edir, INSERT_VALUES, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd(C_dir%da, solution%edir, INSERT_VALUES, ledir, ierr); call CHKERR(ierr)
          call getVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.true.)

          do j = C_one%ys, C_one%ye
            do i = C_one%xs, C_one%xe
              do k = C_one%zs, C_one%ze

                ! Divergence = Incoming - Outgoing

                if (atm%l1d(atmk(atm, k))) then ! one dimensional i.e. twostream

                  cdiv = max(zero, one - atm%a33(atmk(atm, k), i, j) - atm%a13(atmk(atm, k), i, j) - atm%a23(atmk(atm, k), i, j))
                  do isrc = 0, solver%dirtop%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(isrc, k, i, j) * cdiv
                  end do

                else
                  dir2dir(0:C_dir%dof - 1, 0:C_dir%dof - 1) => solver%dir2dir(1:C_dir%dof * C_dir%dof, k, i, j)
                  dir2diff(0:C_dir%dof - 1, 0:C_diff%dof - 1) => solver%dir2diff(1:C_dir%dof * C_diff%dof, k, i, j)

                  idof = 0
                  do isrc = 0, solver%dirtop%dof - 1
                    cdiv = one - sum(dir2dir(isrc, :)) - sum(dir2diff(isrc, :))
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(isrc, k, i, j) * cdiv
                    idof = idof + 1
                  end do
                  do isrc = 0, solver%dirside%dof - 1
                    cdiv = one - sum(dir2dir(idof, :)) - sum(dir2diff(idof, :))
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(idof, k, i + 1 - sun%xinc, j) * cdiv
                    idof = idof + 1
                  end do
                  do isrc = 0, solver%dirside%dof - 1
                    cdiv = one - sum(dir2dir(idof, :)) - sum(dir2diff(idof, :))
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + xedir(idof, k, i, j + 1 - sun%yinc) * cdiv
                    idof = idof + 1
                  end do

                end if

              end do
            end do
          end do

          call restoreVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.true.)
          call DMRestoreLocalVector(C_dir%da, ledir, ierr); call CHKERR(ierr)
        end if
      end associate

      associate ( &
        atm => solver%atm, &
        C_diff => solver%C_diff, &
        C_one => solver%C_one)

        ! Copy ghosted values for diffuse vec
        call DMGetLocalVector(C_diff%da, lediff, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(C_diff%da, solution%ediff, INSERT_VALUES, lediff, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd(C_diff%da, solution%ediff, INSERT_VALUES, lediff, ierr); call CHKERR(ierr)
        call getVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.true.)

        do j = C_one%ys, C_one%ye
          do i = C_one%xs, C_one%xe
            do k = C_one%zs, C_one%ze

              ! Divergence = Incoming - Outgoing

              if (atm%l1d(atmk(atm, k))) then ! one dimensional i.e. twostream

                cdiv = max(zero, one - atm%a11(atmk(atm, k), i, j) - atm%a12(atmk(atm, k), i, j))
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) + xediff(isrc, msrc, i, j) * cdiv
                end do

              else

                diff2diff(0:C_diff%dof - 1, 0:C_diff%dof - 1) => solver%diff2diff(1:C_diff%dof**2, k, i, j)

                idof = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  cdiv = one - sum(diff2diff(idof, :))
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) + xediff(idof, msrc, i, j) * cdiv
                  idof = idof + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  cdiv = one - sum(diff2diff(idof, :))
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) + xediff(idof, k, msrc, j) * cdiv
                  idof = idof + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  cdiv = one - sum(diff2diff(idof, :))
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) + xediff(idof, k, i, msrc) * cdiv
                  idof = idof + 1
                end do

              end if

            end do
          end do
        end do

        if (solution%lthermal_rad) then ! subtract local emissions
          call DMGetLocalVector(C_diff%da, local_b, ierr); call CHKERR(ierr)
          call VecSet(local_b, zero, ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C_diff%da, solver%b, ADD_VALUES, local_b, ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd(C_diff%da, solver%b, ADD_VALUES, local_b, ierr); call CHKERR(ierr)

          call getVecPointer(C_diff%da, local_b, xsrc1d, xsrc, readonly=.true.)
          do j = C_one%ys, C_one%ye
            do i = C_one%xs, C_one%xe
              do k = C_one%zs, C_one%ze
                idof = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k + 1, k, solver%difftop%is_inward(i1 + isrc))
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) - xsrc(idof, msrc, i, j)
                  idof = idof + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i + 1, i, solver%diffside%is_inward(i1 + isrc))
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) - xsrc(idof, k, msrc, j)
                  idof = idof + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j + 1, j, solver%diffside%is_inward(i1 + isrc))
                  xabso(i0, k, i, j) = xabso(i0, k, i, j) - xsrc(idof, k, i, msrc)
                  idof = idof + 1
                end do
              end do
            end do
          end do
          call restoreVecPointer(C_diff%da, local_b, xsrc1d, xsrc, readonly=.true.)
          call DMRestoreLocalVector(C_diff%da, local_b, ierr); call CHKERR(ierr)
        end if

        call restoreVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.true.)
        call DMRestoreLocalVector(C_diff%da, lediff, ierr); call CHKERR(ierr)
      end associate

      call restoreVecPointer(solver%C_one%da, solution%abso, xabso1d, xabso)
    end subroutine

    subroutine compute_absorption_by_flx_divergence()
      associate ( &
        atm => solver%atm, &
        C_dir => solver%C_dir, &
        C_diff => solver%C_diff, &
        C_one => solver%C_one)

        if (solution%lsolar_rad) then
          ! Copy ghosted values for direct vec
          call DMGetLocalVector(C_dir%da, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C_dir%da, solution%edir, INSERT_VALUES, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd(C_dir%da, solution%edir, INSERT_VALUES, ledir, ierr); call CHKERR(ierr)
          call getVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.true.)
        end if

        ! Copy ghosted values for diffuse vec
        call DMGetLocalVector(C_diff%da, lediff, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(C_diff%da, solution%ediff, INSERT_VALUES, lediff, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd(C_diff%da, solution%ediff, INSERT_VALUES, lediff, ierr); call CHKERR(ierr)
        call getVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.true.)

        call getVecPointer(C_one%da, solution%abso, xabso1d, xabso)

        ! calculate absorption by flux divergence
        !DIR$ IVDEP
        do j = C_one%ys, C_one%ye
          do i = C_one%xs, C_one%xe
            do k = C_one%zs, C_one%ze

              ! Divergence = Incoming - Outgoing

              if (atm%l1d(atmk(atm, k))) then ! one dimensional i.e. twostream
                if (solution%lsolar_rad) then
                  do isrc = 0, solver%dirtop%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xedir(isrc, k, i, j) - xedir(isrc, k + i1, i, j))
                  end do
                end if

                do isrc = 0, solver%difftop%dof - 1
                  if (solver%difftop%is_inward(i1 + isrc)) then
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(isrc, k, i, j) - xediff(isrc, k + 1, i, j))
                  else
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(isrc, k + 1, i, j) - xediff(isrc, k, i, j))
                  end if
                end do

              else ! 3D-radiation

                ! direct part of absorption
                if (solution%lsolar_rad) then
                  xinc = solver%sun%xinc
                  yinc = solver%sun%yinc

                  src = 0
                  do isrc = 0, solver%dirtop%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xedir(src, k, i, j) - xedir(src, k + i1, i, j))
                    src = src + 1
                  end do

                  do isrc = 0, solver%dirside%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xedir(src, k, i + 1 - xinc, j) - xedir(src, k, i + xinc, j))
                    src = src + 1
                  end do

                  do isrc = 0, solver%dirside%dof - 1
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xedir(src, k, i, j + i1 - yinc) - xedir(src, k, i, j + yinc))
                    src = src + 1
                  end do
                end if

                ! diffuse part of absorption
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  if (solver%difftop%is_inward(i1 + isrc)) then
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(src, k, i, j) - xediff(src, k + 1, i, j))
                  else
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(src, k + 1, i, j) - xediff(src, k, i, j))
                  end if
                  src = src + 1
                end do

                do isrc = 0, solver%diffside%dof - 1
                  if (solver%diffside%is_inward(i1 + isrc)) then
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(src, k, i, j) - xediff(src, k, i + 1, j))
                  else
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(src, k, i + 1, j) - xediff(src, k, i, j))
                  end if
                  src = src + 1
                end do

                do isrc = 0, solver%diffside%dof - 1
                  if (solver%diffside%is_inward(i1 + isrc)) then
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(src, k, i, j) - xediff(src, k, i, j + 1))
                  else
                    xabso(i0, k, i, j) = xabso(i0, k, i, j) + (xediff(src, k, i, j + 1) - xediff(src, k, i, j))
                  end if
                  src = src + 1
                end do

                if (ldebug) then
                  if (isnan(xabso(i0, k, i, j))) then
                    print *, 'nan in flxdiv', k, i, j, '::', xabso(i0, k, i, j)
                  end if
                end if
              end if ! 1d/3D
            end do
          end do
        end do

        if (solution%lsolar_rad) then
          call restoreVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.true.)
          call DMRestoreLocalVector(C_dir%da, ledir, ierr); call CHKERR(ierr)
        end if

        call restoreVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.true.)
        call DMRestoreLocalVector(C_diff%da, lediff, ierr); call CHKERR(ierr)

        call restoreVecPointer(C_one%da, solution%abso, xabso1d, xabso)
      end associate
    end subroutine

    subroutine gen_abso_scalevec()
      associate (atm => solver%atm, C_one => solver%C_one)

        allocate (solver%abso_scalevec)
        call VecDuplicate(solution%abso, solver%abso_scalevec, ierr); call CHKERR(ierr)
        call getVecPointer(C_one%da, solver%abso_scalevec, xabso1d, xabso)
        Az = atm%dx * atm%dy

        do j = C_one%ys, C_one%ye
          do i = C_one%xs, C_one%xe
            do k = C_one%zs, C_one%ze
              Volume = Az * atm%dz(atmk(atm, k), i, j)
              xabso(i0, k, i, j) = one / Volume
            end do
          end do
        end do

        ! here a special case for icollapse, take the dz of all layers above dynamical grid
        do j = C_one%ys, C_one%ye
          do i = C_one%xs, C_one%xe
            Volume = Az * sum(atm%dz(C_one%zs:atmk(atm, C_one%zs), i, j))
            xabso(i0, C_one%zs, i, j) = one / Volume
          end do
        end do

        call restoreVecPointer(C_one%da, solver%abso_scalevec, xabso1d, xabso)
      end associate
    end subroutine
  end subroutine

  !> @brief generate matrix col/row permutations
  subroutine gen_mat_permutation(A, C, rev_x, rev_y, rev_z, zlast, switch_xy, perm_info, prefix)
    type(tMat), intent(in) :: A
    type(t_coord), intent(in) :: C
    logical, intent(in) :: rev_x, rev_y, rev_z, zlast, switch_xy
    type(t_mat_permute_info), allocatable, intent(inout) :: perm_info
    character(len=*), intent(in) :: prefix

    logical :: lflg
    integer(iintegers), dimension(3) :: dd, dx, dy, dz ! start, end, increment for each dimension

    integer(iintegers), allocatable :: is_data(:)
    integer(iintegers) :: Astart, Aend, m, k, i, j, d, da_offsets(4)
    logical :: opt_rev_x, opt_rev_y, opt_rev_z, opt_zlast, opt_switch_xy

    integer(mpiint) :: comm, ierr

    dd = [i0, C%dof - i1, i1]
    dz = [i0, C%zm - i1, i1]
    dx = [i0, C%xm - i1, i1]
    dy = [i0, C%ym - i1, i1]

    opt_rev_y = rev_y
    call get_petsc_opt(prefix, "-mat_permute_rev_y", opt_rev_y, lflg, ierr); call CHKERR(ierr)
    if (opt_rev_y) dy = [dy(2), dy(1), -dy(3)]

    opt_rev_x = rev_x
    call get_petsc_opt(prefix, "-mat_permute_rev_x", opt_rev_x, lflg, ierr); call CHKERR(ierr)
    if (opt_rev_x) dx = [dx(2), dx(1), -dx(3)]

    opt_rev_z = rev_z
    call get_petsc_opt(prefix, "-mat_permute_rev_z", opt_rev_z, lflg, ierr); call CHKERR(ierr)
    if (opt_rev_z) dz = [dz(2), dz(1), -dz(3)]

    opt_switch_xy = switch_xy
    call get_petsc_opt(prefix, "-mat_permute_ij", opt_switch_xy, lflg, ierr); call CHKERR(ierr)

    opt_zlast = zlast
    call get_petsc_opt(prefix, "-mat_permute_z", opt_zlast, lflg, ierr); call CHKERR(ierr)

    if (.not. allocated(perm_info)) then
      allocate (perm_info)
    else
      if (all([ &
        & opt_rev_x .eqv. perm_info%rev_x,    &
        & opt_rev_y .eqv. perm_info%rev_y,    &
        & opt_rev_z .eqv. perm_info%rev_z,    &
        & opt_zlast .eqv. perm_info%zlast,    &
        & opt_switch_xy .eqv. perm_info%switch_xy &
        & ])) then ! up to date
        return
      else
        call ISDestroy(perm_info%is, ierr); call CHKERR(ierr)
      end if
    end if

    perm_info%rev_x = opt_rev_x
    perm_info%rev_y = opt_rev_y
    perm_info%rev_z = opt_rev_z
    perm_info%zlast = opt_zlast
    perm_info%switch_xy = opt_switch_xy

    call MatGetOwnershipRange(A, Astart, Aend, ierr); call CHKERR(ierr)
    call ndarray_offsets([C%dof, C%zm, C%xm, C%ym], da_offsets)

    allocate (is_data(Astart:Aend - 1))
    m = Astart
    if (perm_info%zlast) then
      if (perm_info%switch_xy) then
        do k = dz(1), dz(2), dz(3)
          do i = dx(1), dx(2), dx(3)
            do j = dy(1), dy(2), dy(3)
              do d = dd(1), dd(2), dd(3)
                is_data(m) = Astart + ind_nd_to_1d(da_offsets, [d, k, i, j], cstyle=.true.)
                m = m + 1
              end do
            end do
          end do
        end do
      else
        do k = dz(1), dz(2), dz(3)
          do j = dy(1), dy(2), dy(3)
            do i = dx(1), dx(2), dx(3)
              do d = dd(1), dd(2), dd(3)
                is_data(m) = Astart + ind_nd_to_1d(da_offsets, [d, k, i, j], cstyle=.true.)
                m = m + 1
              end do
            end do
          end do
        end do
      end if
    else
      if (perm_info%switch_xy) then
        do i = dx(1), dx(2), dx(3)
          do j = dy(1), dy(2), dy(3)
            do k = dz(1), dz(2), dz(3)
              do d = dd(1), dd(2), dd(3)
                is_data(m) = Astart + ind_nd_to_1d(da_offsets, [d, k, i, j], cstyle=.true.)
                m = m + 1
              end do
            end do
          end do
        end do
      else
        do j = dy(1), dy(2), dy(3)
          do i = dx(1), dx(2), dx(3)
            do k = dz(1), dz(2), dz(3)
              do d = dd(1), dd(2), dd(3)
                is_data(m) = Astart + ind_nd_to_1d(da_offsets, [d, k, i, j], cstyle=.true.)
                m = m + 1
              end do
            end do
          end do
        end do
      end if
    end if

    call PetscObjectGetComm(A, comm, ierr); call CHKERR(ierr)
    call ISCreateGeneral(comm, size(is_data, kind=iintegers), is_data, PETSC_COPY_VALUES, perm_info%is, ierr); call CHKERR(ierr)
  end subroutine

  !> @brief call PETSc Krylov Subspace Solver
  !> @details solve with ksp and save residual history of solver
  !> \n -- this may be handy later to decide next time if we have to calculate radiation again
  !> \n if we did not get convergence, we try again with standard GMRES and a resetted(zero) initial guess -- if that doesnt help, we got a problem!
  subroutine solve(solver, ksp, b, x, uid, iter, ksp_residual_history)
    class(t_solver) :: solver
    type(tKSP) :: ksp
    type(tVec) :: b
    type(tVec) :: x
    integer(iintegers) :: uid
    integer(iintegers), intent(out) :: iter
    real(ireals), intent(inout), optional :: ksp_residual_history(:)

    KSPConvergedReason :: reason

    character(len=default_str_len) :: prefix
    KSPType :: old_ksp_type

    logical :: lskip_ksp_solve, laccept_incomplete_solve, lflg
    integer(mpiint) :: ierr

    if (solver%myid .eq. 0 .and. ldebug) print *, 'Solving Matrix'

    if (present(ksp_residual_history)) then
      call KSPSetResidualHistory(ksp, ksp_residual_history, &
        & size(ksp_residual_history, kind=iintegers), PETSC_TRUE, ierr); call CHKERR(ierr)
    end if

    lskip_ksp_solve = .false.
    call get_petsc_opt(solver%prefix, "-skip_ksp_solve", lskip_ksp_solve, lflg, ierr); call CHKERR(ierr)
    if (lskip_ksp_solve) then
      call VecCopy(b, x, ierr); call CHKERR(ierr)
      return
    end if

    call hegedus_trick(ksp, b, x)
    call KSPSolve(ksp, b, x, ierr); call CHKERR(ierr)
    call KSPGetIterationNumber(ksp, iter, ierr); call CHKERR(ierr)
    call KSPGetConvergedReason(ksp, reason, ierr); call CHKERR(ierr)

    ! if(reason.eq.KSP_DIVERGED_ITS) then
    !   if(solver%myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
    !   return
    ! endif
    laccept_incomplete_solve = .false.
    call get_petsc_opt(solver%prefix, "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)
    if (laccept_incomplete_solve) return

    if (reason .le. 0) then
      call KSPGetOptionsPrefix(ksp, prefix, ierr); call CHKERR(ierr)
      if (solver%myid .eq. 0) &
        & call CHKWARN(int(reason, mpiint), trim(prefix)//' :: Resetted initial guess to zero '// &
        & 'and try again with gmres (uid='//toStr(uid)//')')
      call VecSet(x, zero, ierr); call CHKERR(ierr)
      call KSPGetType(ksp, old_ksp_type, ierr); call CHKERR(ierr)
      call KSPSetType(ksp, KSPGMRES, ierr); call CHKERR(ierr)
      call KSPSetUp(ksp, ierr); call CHKERR(ierr)
      call KSPSolve(ksp, b, x, ierr); call CHKERR(ierr)
      call KSPGetIterationNumber(ksp, iter, ierr); call CHKERR(ierr)
      call KSPGetConvergedReason(ksp, reason, ierr); call CHKERR(ierr)

      ! And return to normal solver...
      call KSPSetType(ksp, old_ksp_type, ierr); call CHKERR(ierr)
      call KSPSetFromOptions(ksp, ierr); call CHKERR(ierr)
      call KSPSetUp(ksp, ierr); call CHKERR(ierr)
      if (solver%myid .eq. 0 .and. ldebug) &
        & print *, solver%myid, 'Solver took', iter, 'iterations and converged', reason .gt. 0, 'because', reason, 'uid', uid
    end if

    if (reason .le. 0) then
      call CHKERR(int(reason, mpiint), &
        & '***** SOLVER did NOT converge :( -- (uid='//toStr(uid)//')'// &
        & 'if you know what you are doing, you can use the option -accept_incomplete_solve to continue')
    end if
  end subroutine

  !> @brief initialize PETSc Krylov Subspace Solver
  !> @details default KSP solver is a FGMRES with BJCAOBI // ILU(1)
  !> \n -- the default does however not scale well -- and we configure petsc solvers per commandline anyway
  !> \n -- see documentation for details on how to do so
  subroutine setup_ksp(solver, ksp, C, A, prefix)
    class(t_solver) :: solver
    type(tKSP), intent(inout), allocatable :: ksp
    type(t_coord), intent(inout) :: C
    type(tMat), intent(in) :: A
    character(len=*), intent(in), optional :: prefix

    integer(iintegers), parameter :: maxiter = 1000

    integer(mpiint) :: myid, numnodes

    real(ireals) :: rtol, atol

    type(tPC) :: prec
    logical :: prec_is_set
    integer(mpiint) :: ierr
    character(len=default_str_len) :: kspprefix

    if (allocated(ksp)) return

    call set_dmda_cell_coordinates(solver, solver%atm, C%da, ierr); call CHKERR(ierr)

    call mpi_comm_rank(C%comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(C%comm, numnodes, ierr); call CHKERR(ierr)

    allocate (ksp)
    call KSPCreate(C%comm, ksp, ierr); call CHKERR(ierr)
    if (present(prefix)) then
      call KSPAppendOptionsPrefix(ksp, trim(prefix), ierr); call CHKERR(ierr)
    end if
    call KSPGetOptionsPrefix(ksp, kspprefix, ierr); call CHKERR(ierr)

    call KSPSetType(ksp, KSPBCGS, ierr); call CHKERR(ierr)
    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr); call CHKERR(ierr)

    prec_is_set = .false.
    call PetscOptionsHasName(PETSC_NULL_OPTIONS, trim(kspprefix), '-pc_type', prec_is_set, ierr); call CHKERR(ierr)

    if (.not. prec_is_set) then
      !call CHKWARN(1_mpiint, 'no preconditioner setting found, applying defaults')
      call KSPGetPC(ksp, prec, ierr); call CHKERR(ierr)
      if (numnodes .le. 1_mpiint) then
        if (C%glob_xm .le. 3_iintegers .and. C%glob_ym .le. 3_iintegers) then ! small domains can use direct solves
          call KSPSetInitialGuessNonzero(ksp, PETSC_FALSE, ierr); call CHKERR(ierr)
          call KSPSetType(ksp, KSPPREONLY, ierr); call CHKERR(ierr)
          call PCSetType(prec, PCLU, ierr); call CHKERR(ierr)
          call PCFactorSetFill(prec, 16._ireals, ierr); call CHKERR(ierr)
        else
          call PCSetType(prec, PCILU, ierr); call CHKERR(ierr)
        end if
      else
        if (trim(prefix) .eq. 'solar_dir_') then
          call PCSetType(prec, PCSOR, ierr); call CHKERR(ierr)
          call PCSORSetSymmetric(prec, SOR_LOCAL_FORWARD_SWEEP, ierr); call CHKERR(ierr)
        else
          call PCSetType(prec, PCBJACOBI, ierr); call CHKERR(ierr)

          !call PCSetType(prec, PCSOR, ierr); call CHKERR(ierr)

          !call PCSetType(prec, PCASM, ierr); call CHKERR(ierr)
          !call PCASMSetOverlap(prec, i1, ierr); call CHKERR(ierr)
        end if
      end if
    end if

    call determine_ksp_tolerances(C, solver%atm%l1d, rtol, atol)
    call KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT_REAL, maxiter, ierr); call CHKERR(ierr)

    call KSPSetConvergenceTest(ksp, MyKSPConverged, 0, PETSC_NULL_FUNCTION, ierr); call CHKERR(ierr)

    call KSPSetOperators(ksp, A, A, ierr); call CHKERR(ierr)

    call KSPSetDM(ksp, C%da, ierr); call CHKERR(ierr)
    call KSPSetDMActive(ksp, PETSC_FALSE, ierr); call CHKERR(ierr)

    call KSPSetFromOptions(ksp, ierr); call CHKERR(ierr)
    call KSPSetUp(ksp, ierr); call CHKERR(ierr)

    if (.not. prec_is_set) then
      default_preconditioner_settings: block
        integer(iintegers) :: i, asm_N, asm_iter
        integer(iintegers) :: first_local
        type(tKSP), allocatable :: asm_ksps(:)
        type(tPC) :: subpc
        PCType :: pctype
        logical :: lflg

        call PCGetType(prec, pctype, ierr); call CHKERR(ierr)
        if (pctype .eq. PCASM) then

          asm_iter = 1
          call get_petsc_opt(solver%prefix, "-ts_ksp_iter", asm_iter, lflg, ierr); call CHKERR(ierr)

          call PCASMGetSubKSP(prec, asm_N, first_local, PETSC_NULL_KSP, ierr); call CHKERR(ierr)
          allocate (asm_ksps(asm_N))
          call PCASMGetSubKSP(prec, asm_N, first_local, asm_ksps, ierr); call CHKERR(ierr)
          do i = 1, asm_N
            call KSPSetType(asm_ksps(i), KSPRICHARDSON, ierr); call CHKERR(ierr)
            call KSPRichardsonSetSelfScale(asm_ksps(i), PETSC_TRUE, ierr); call CHKERR(ierr)
            call KSPSetNormType(asm_ksps(i), KSP_NORM_PRECONDITIONED, ierr); call CHKERR(ierr)
            call KSPSetTolerances(asm_ksps(i), &
              & PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
              & asm_iter, ierr); call CHKERR(ierr)
            call KSPGetPC(asm_ksps(i), subpc, ierr); call CHKERR(ierr)
            call PCSetType(subpc, PCSOR, ierr); call CHKERR(ierr)
            call KSPSetFromOptions(asm_ksps(i), ierr); call CHKERR(ierr)
          end do

        else if (pctype .eq. PCBJACOBI) then

          call PCBJacobiGetSubKSP(prec, asm_N, first_local, PETSC_NULL_KSP, ierr); call CHKERR(ierr)
          allocate (asm_ksps(asm_N))
          call PCBjacobiGetSubKSP(prec, asm_N, first_local, asm_ksps, ierr); call CHKERR(ierr)
          do i = 1, asm_N
            call KSPSetType(asm_ksps(i), KSPPREONLY, ierr); call CHKERR(ierr)
            call KSPGetPC(asm_ksps(i), subpc, ierr); call CHKERR(ierr)
            call PCSetType(subpc, PCILU, ierr); call CHKERR(ierr)
            call KSPSetFromOptions(asm_ksps(i), ierr); call CHKERR(ierr)
          end do

        end if
      end block default_preconditioner_settings
    end if

    if (myid .eq. 0 .and. ldebug) print *, 'Setup KSP done'
  end subroutine

  !> @brief override convergence tests -- the normal KSPConverged returns bad solution if no iterations are needed for convergence
  subroutine MyKSPConverged(ksp, n, rnorm, flag, dummy, ierr)
    ! Input Parameters:
    !    ksp   - iterative context
    !    n     - iteration number
    !    rnorm - 2-norm (preconditioned) residual value (may be estimated)
    !    dummy - optional user-defined monitor context (unused here)
    type(tKSP) :: ksp
    integer(mpiint) :: ierr
    integer(iintegers) :: n, dummy
    KSPConvergedReason :: flag
    real(ireals) :: rnorm

    real(ireals) :: rtol, atol, dtol
    real(ireals), save :: initial_rnorm
    integer(iintegers) :: maxits

    call KSPGetTolerances(ksp, rtol, atol, dtol, maxits, ierr)

    flag = 0
    if (n .eq. 0) then
      initial_rnorm = max(tiny(rnorm), rnorm)
      return
    end if

    if (rnorm / initial_rnorm .le. rtol) then
      flag = 2
      return
    end if

    if (rnorm .le. atol) then
      flag = 3
      return
    end if

    if (n .gt. maxits) then
      flag = -3
      return
    end if

    if (rnorm / initial_rnorm .ge. dtol) then
      flag = -4
      return
    end if

    if (isnan(rnorm)) then
      flag = -9
      return
    end if
    if (.false.) dummy = dummy + 1 ! stupid statement to remove unused variable warning
  end subroutine

  !> @brief build direct radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_dir_coeff(solver, sun, A, C)
    class(t_solver), intent(in) :: solver
    type(t_suninfo), intent(in) :: sun
    type(tMat), intent(inout) :: A
    type(t_coord), intent(in) :: C

    integer(iintegers) :: i, j, k
    integer(mpiint) :: ierr

    if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid, 'setup_direct_matrix ...'

    call MatZeroEntries(A, ierr); call CHKERR(ierr)
    call mat_set_diagonal(A)

    do j = C%ys, C%ye
      do i = C%xs, C%xe
        do k = C%zs, C%ze - 1

          if (solver%atm%l1d(atmk(solver%atm, k))) then
            call set_eddington_coeff(solver%atm, A, k, i, j)
          else
            call set_pprts_coeff(solver, C, A, k, i, j)
          end if

        end do
      end do
    end do

    if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid, 'setup_direct_matrix done'

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, "-show_Mdir", ierr); call CHKERR(ierr)
  contains

    subroutine set_pprts_coeff(solver, C, A, k, i, j)
      class(t_solver) :: solver
      type(t_coord), intent(in) :: C
      type(tMat), intent(inout) :: A
      integer(iintegers), intent(in) :: i, j, k

      MatStencil :: row(4, C%dof), col(4, C%dof)
      real(ireals) :: norm

      integer(iintegers) :: dst, src, xinc, yinc, isrc, idst

      xinc = sun%xinc
      yinc = sun%yinc

      do idst = 1, solver%dirtop%dof
        dst = idst
        row(MatStencil_j, dst) = i
        row(MatStencil_k, dst) = j
        row(MatStencil_i, dst) = k + 1
        row(MatStencil_c, dst) = dst - i1 ! Define transmission towards the lower/upper lid
      end do

      do idst = 1, solver%dirside%dof
        dst = idst + solver%dirtop%dof
        row(MatStencil_j, dst) = i + xinc
        row(MatStencil_k, dst) = j
        row(MatStencil_i, dst) = k
        row(MatStencil_c, dst) = dst - i1 ! Define transmission towards the left/right lid
      end do

      do idst = 1, solver%dirside%dof
        dst = idst + solver%dirtop%dof + solver%dirside%dof
        row(MatStencil_j, dst) = i
        row(MatStencil_k, dst) = j + yinc
        row(MatStencil_i, dst) = k
        row(MatStencil_c, dst) = dst - i1 ! Define transmission towards the front/back lid
      end do

      do isrc = 1, solver%dirtop%dof
        src = isrc
        col(MatStencil_j, src) = i
        col(MatStencil_k, src) = j
        col(MatStencil_i, src) = k
        col(MatStencil_c, src) = src - i1 ! Define transmission towards the lower/upper lid
      end do

      do isrc = 1, solver%dirside%dof
        src = isrc + solver%dirtop%dof
        col(MatStencil_j, src) = i + 1 - xinc
        col(MatStencil_k, src) = j
        col(MatStencil_i, src) = k
        col(MatStencil_c, src) = src - i1 ! Define transmission towards the left/right lid
      end do

      do isrc = 1, solver%dirside%dof
        src = isrc + solver%dirtop%dof + solver%dirside%dof
        col(MatStencil_j, src) = i
        col(MatStencil_k, src) = j + 1 - yinc
        col(MatStencil_i, src) = k
        col(MatStencil_c, src) = src - i1 ! Define transmission towards the front/back lid
      end do

      call MatSetValuesStencil(A, C%dof, row, C%dof, col, &
        & real(-solver%dir2dir(:, k, i, j), ireals), INSERT_VALUES, ierr); call CHKERR(ierr)

      if (ldebug) then
        do src = 1, C%dof
          norm = sum(solver%dir2dir(src:C%dof**2:C%dof, k, i, j))
          if (norm .gt. one + 10._ireals * epsilon(norm)) then ! could renormalize
            if (norm .gt. one + 10._ireals * sqrt(epsilon(norm))) then ! fatally off
              print *, 'direct sum(dst==', dst, ') gt one', norm
              print *, 'direct coeff', norm, '::', solver%dir2dir(src:C%dof**2:C%dof, k, i, j)
              call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
            end if ! fatally off
          end if ! could renormalize
        end do
      end if
    end subroutine

    subroutine set_eddington_coeff(atm, A, k, i, j)
      type(t_atmosphere), intent(in) :: atm
      type(tMat), intent(inout) :: A
      integer(iintegers), intent(in) :: i, j, k

      MatStencil :: row(4, 1), col(4, 1)
      real(ireals) :: v(1)
      integer(iintegers) :: src

      v = atm%a33(atmk(atm, k), i, j)

      col(MatStencil_j, i1) = i; col(MatStencil_k, i1) = j; col(MatStencil_i, i1) = k
      row(MatStencil_j, i1) = i; row(MatStencil_k, i1) = j; row(MatStencil_i, i1) = k + 1

      do src = i1, solver%dirtop%dof
        col(MatStencil_c, i1) = src - i1 ! Source may be the upper/lower lid:
        row(MatStencil_c, i1) = src - i1 ! Define transmission towards the lower/upper lid
        call MatSetValuesStencil(A, i1, row, i1, col, -v, INSERT_VALUES, ierr); call CHKERR(ierr)
      end do
    end subroutine

  end subroutine set_dir_coeff

  !> @brief setup source term for diffuse radiation
  !> @details this is either direct radiation scattered into one of the diffuse coeffs:
  !> \n direct source term is
  !> \n   direct radiation times the dir2diff coeffs
  !> \n or it may be that we have a source term due to thermal emission --
  !> \n   to determine emissivity of box, we use the forward transport coefficients backwards
  !> \n   a la: transmissivity $T = \sum(coeffs)$ and therefore emissivity $E = 1 - T$
  subroutine setup_b(solver, solution, b, opt_buildings)
    class(t_solver), target, intent(in) :: solver
    type(t_state_container), intent(in) :: solution
    type(tVec), intent(inout) :: b
    type(t_pprts_buildings), intent(in), optional :: opt_buildings

    type(tVec) :: local_b
    type(tVec) :: local_edir

    real(ireals), pointer, dimension(:, :, :, :) :: xsrc => null()
    real(ireals), pointer, dimension(:) :: xsrc1d => null()
    integer(mpiint) :: ierr

    associate (atm => solver%atm, &
               C_dir => solver%C_dir, &
               C_diff => solver%C_diff)

      if (solver%myid .eq. 0 .and. ldebug) print *, 'src Vector Assembly...'

      call DMGetLocalVector(C_diff%da, local_b, ierr); call CHKERR(ierr)
      call VecSet(local_b, zero, ierr); call CHKERR(ierr)

      call getVecPointer(C_diff%da, local_b, xsrc1d, xsrc)

      if (solution%lsolar_rad) then
        ! Copy ghosted values for direct vec
        call DMGetLocalVector(C_dir%da, local_edir, ierr); call CHKERR(ierr)
        call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(C_dir%da, solution%edir, ADD_VALUES, local_edir, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd(C_dir%da, solution%edir, ADD_VALUES, local_edir, ierr); call CHKERR(ierr)

        call set_solar_source(solver, local_edir)
      end if

      if (solution%lthermal_rad) &
        call set_thermal_source()

      if (present(opt_buildings)) then
        if (solution%lsolar_rad) then
          call set_buildings_reflection(solver, local_edir, opt_buildings)
        end if

        if (solution%lthermal_rad) then
          call set_buildings_emission(solver, opt_buildings, ierr); call CHKERR(ierr)
        end if
      end if

      if (solver%myid .eq. 0 .and. ldebug) print *, 'src Vector Assembly... setting coefficients ...done'

      if (solution%lsolar_rad) then
        call DMRestoreLocalVector(C_dir%da, local_edir, ierr); call CHKERR(ierr)
      end if

      call restoreVecPointer(C_diff%da, local_b, xsrc1d, xsrc)

      call VecSet(b, zero, ierr); call CHKERR(ierr) ! reset global Vec

      call DMLocalToGlobalBegin(C_diff%da, local_b, ADD_VALUES, b, ierr); call CHKERR(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
      call DMLocalToGlobalEnd(C_diff%da, local_b, ADD_VALUES, b, ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(C_diff%da, local_b, ierr); call CHKERR(ierr)

      if (solver%myid .eq. 0 .and. ldebug) print *, 'src Vector Assembly done'
    end associate
  contains
    subroutine set_thermal_source()

      real(ireals) :: Ax, Ay, Az, emis, b0, b1, btop, bbot, bfac, tauz
      integer(iintegers) :: k, i, j, src, iside, ak
      real(ireals), pointer :: diff2diff(:, :) ! dim(dst, src)

      associate (atm => solver%atm, &
                 C_diff => solver%C_diff)

        if (.not. allocated(atm%planck)) then
          ierr = 1
          call CHKERR(ierr, 'have thermal computation but planck data is not allocated')
        end if
        if (solver%myid .eq. 0 .and. ldebug) &
          print *, 'Assembly of SRC-Vector ... setting thermal source terms min/max planck', &
          minval(atm%planck), maxval(atm%planck)

        Az = atm%dx * atm%dy / real(solver%difftop%area_divider, ireals)

        do j = C_diff%ys, C_diff%ye
          do i = C_diff%xs, C_diff%xe
            do k = C_diff%zs, C_diff%ze - 1
              ak = atmk(atm, k)
              b0 = atm%planck(ak, i, j)
              b1 = atm%planck(ak + 1, i, j)
              tauz = atm%kabs(ak, i, j) * atm%dz(ak, i, j)

              if (atm%l1d(ak)) then

                bfac = pi * Az / real(solver%difftop%streams, ireals)

                if (atm%lcollapse .and. k .eq. i0) then
                  btop = atm%Btop(i, j) * bfac
                  bbot = atm%Bbot(i, j) * bfac

                else
                  emis = min(one, max(zero, one - atm%a11(ak, i, j) - atm%a12(ak, i, j)))

                  call B_eff(b1, b0, tauz, btop)
                  call B_eff(b0, b1, tauz, bbot)
                  btop = btop * bfac * emis
                  bbot = bbot * bfac * emis
                end if

                do src = 0, solver%difftop%dof - 1
                  if (solver%difftop%is_inward(i1 + src)) then !Edn
                    xsrc(src, k + 1, i, j) = xsrc(src, k + 1, i, j) + bbot
                  else !E_up
                    xsrc(src, k, i, j) = xsrc(src, k, i, j) + btop
                  end if
                end do

              else ! Tenstream source terms
                Ax = atm%dy * atm%dz(ak, i, j) / real(solver%diffside%area_divider, ireals)
                Ay = atm%dx * atm%dz(ak, i, j) / real(solver%diffside%area_divider, ireals)

                call B_eff(b1, b0, tauz, btop)
                call B_eff(b0, b1, tauz, bbot)

                diff2diff(0:C_diff%dof - 1, 0:C_diff%dof - 1) => solver%diff2diff(1:C_diff%dof**2, k, i, j)

                src = 0
                bfac = pi * Az / real(solver%difftop%streams, ireals)
                do iside = 1, solver%difftop%dof
                  emis = one - sum(diff2diff(src, :))
                  if (solver%difftop%is_inward(iside) .eqv. .false.) then ! outgoing means Eup
                    xsrc(src, k, i, j) = xsrc(src, k, i, j) + btop * bfac * emis
                  else
                    xsrc(src, k + 1, i, j) = xsrc(src, k + 1, i, j) + bbot * bfac * emis
                  end if
                  src = src + 1
                end do

                bfac = pi * Ax / real(solver%diffside%streams, ireals)
                do iside = 1, solver%diffside%dof
                  emis = one - sum(diff2diff(src, :))
                  if (iside .gt. solver%diffside%dof / 2) then ! upward streams
                    emis = btop * emis
                  else
                    emis = bbot * emis
                  end if
                  if (solver%diffside%is_inward(iside) .eqv. .false.) then ! outgoing means towards +x
                    xsrc(src, k, i, j) = xsrc(src, k, i, j) + emis * bfac
                  else
                    xsrc(src, k, i + 1, j) = xsrc(src, k, i + 1, j) + emis * bfac
                  end if
                  src = src + 1
                end do

                bfac = pi * Ay / real(solver%diffside%streams, ireals)
                do iside = 1, solver%diffside%dof
                  emis = one - sum(diff2diff(src, :))
                  if (iside .gt. solver%diffside%dof / 2) then ! upward streams
                    emis = btop * emis
                  else
                    emis = bbot * emis
                  end if
                  if (solver%diffside%is_inward(iside) .eqv. .false.) then ! outgoing means towards +y
                    xsrc(src, k, i, j) = xsrc(src, k, i, j) + emis * bfac
                  else
                    xsrc(src, k, i, j + 1) = xsrc(src, k, i, j + 1) + emis * bfac
                  end if
                  src = src + 1
                end do

              end if ! 1D or Tenstream?

            end do ! k
          end do
        end do

        ! Thermal emission at surface
        k = C_diff%ze
        if (allocated(atm%Bsrfc)) then
          do j = C_diff%ys, C_diff%ye
            do i = C_diff%xs, C_diff%xe
              do src = 0, solver%difftop%dof - 1
                if (.not. solver%difftop%is_inward(i1 + src)) then !Eup
                  xsrc(src, k, i, j) = xsrc(src, k, i, j) + atm%Bsrfc(i, j) &
                                       * Az * (one - atm%albedo(i, j)) * pi / real(solver%difftop%streams, ireals)
                end if
              end do
            end do
          end do
        else
          ak = atmk(atm, k)
          do j = C_diff%ys, C_diff%ye
            do i = C_diff%xs, C_diff%xe
              do src = 0, solver%difftop%dof - 1
                if (.not. solver%difftop%is_inward(i1 + src)) then !Eup
                  xsrc(src, k, i, j) = xsrc(src, k, i, j) + atm%planck(ak, i, j) &
                                       * Az * (one - atm%albedo(i, j)) * pi / real(solver%difftop%streams, ireals)
                end if
              end do
            end do
          end do
        end if
      end associate
    end subroutine

    subroutine set_solar_source(solver, local_edir)
      class(t_solver), target, intent(in) :: solver
      type(tVec), intent(in) :: local_edir

      real(ireals), pointer :: dir2diff(:)
      real(ireals) :: solrad
      integer(iintegers) :: idof, idofdst, idiff
      integer(iintegers) :: k, i, j, src, dst

      real(ireals), pointer, dimension(:, :, :, :) :: xedir => null()
      real(ireals), pointer, dimension(:) :: xedir1d => null()

      associate (atm => solver%atm, &
                 C_dir => solver%C_dir, &
                 C_diff => solver%C_diff, &
                 sun => solver%sun)

        call getVecPointer(C_dir%da, local_edir, xedir1d, xedir)

        if (solver%myid .eq. 0 .and. ldebug) print *, 'Assembly of SRC-Vector .. setting solar source', &
          sum(xedir(i0, C_dir%zs:C_dir%ze, C_dir%xs:C_dir%xe, C_dir%ys:C_dir%ye)) / &
          size(xedir(i0, C_dir%zs:C_dir%ze, C_dir%xs:C_dir%xe, C_dir%ys:C_dir%ye))

        do j = C_diff%ys, C_diff%ye
          do i = C_diff%xs, C_diff%xe
            do k = C_diff%zs, C_diff%ze - 1

              if (any(xedir(:, k, i, j) .gt. epsilon(one))) then
                if (atm%l1d(atmk(atm, k))) then

                  ! Only transport the 4 tiles from dir0 to the Eup and Edn
                  do src = 1, solver%dirtop%dof

                    do idiff = 1, solver%difftop%dof
                      if (solver%difftop%is_inward(idiff)) then
                        ! fetch all diffuse downward fluxes at k+1
                        xsrc(idiff - 1, k + 1, i, j) = xsrc(idiff - 1, k + 1, i, j) &
                                                     & + xedir(src - 1, k, i, j) * atm%a23(atmk(atm, k), i, j) &
                                                     & / real(solver%difftop%streams, ireals)
                      else
                        ! fetch all diffuse upward fluxes at k
                        xsrc(idiff - 1, k, i, j) = xsrc(idiff - 1, k, i, j) &
                                                 & + xedir(src - 1, k, i, j) * atm%a13(atmk(atm, k), i, j) &
                                                 & / real(solver%difftop%streams, ireals)
                      end if
                    end do
                  end do

                else ! Tenstream source terms

                  dir2diff => solver%dir2diff(:, k, i, j)
                  dst = 0
                  do idofdst = 1, solver%difftop%dof
                    src = 1
                    do idof = 1, solver%dirtop%dof
                      solrad = xedir(src - 1, k, i, j)

                      if (solver%difftop%is_inward(idofdst)) then
                        xsrc(dst, k + 1, i, j) = xsrc(dst, k + 1, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do

                    do idof = 1, solver%dirside%dof
                      solrad = xedir(src - 1, k, i + i1 - sun%xinc, j)
                      if (solver%difftop%is_inward(idofdst)) then
                        xsrc(dst, k + 1, i, j) = xsrc(dst, k + 1, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do

                    do idof = 1, solver%dirside%dof
                      solrad = xedir(src - 1, k, i, j + i1 - sun%yinc)
                      if (solver%difftop%is_inward(idofdst)) then
                        xsrc(dst, k + 1, i, j) = xsrc(dst, k + 1, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do
                    dst = dst + 1
                  end do

                  do idofdst = 1, solver%diffside%dof
                    src = 1
                    do idof = 1, solver%dirtop%dof
                      solrad = xedir(src - 1, k, i, j)
                      if (solver%diffside%is_inward(idofdst)) then
                        xsrc(dst, k, i + 1, j) = xsrc(dst, k, i + 1, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do

                    do idof = 1, solver%dirside%dof
                      solrad = xedir(src - 1, k, i + i1 - sun%xinc, j)
                      if (solver%diffside%is_inward(idofdst)) then
                        xsrc(dst, k, i + 1, j) = xsrc(dst, k, i + 1, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do

                    do idof = 1, solver%dirside%dof
                      solrad = xedir(src - 1, k, i, j + i1 - sun%yinc)
                      if (solver%diffside%is_inward(idofdst)) then
                        xsrc(dst, k, i + 1, j) = xsrc(dst, k, i + 1, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do
                    dst = dst + 1
                  end do

                  do idofdst = 1, solver%diffside%dof
                    src = 1
                    do idof = 1, solver%dirtop%dof
                      solrad = xedir(src - 1, k, i, j)
                      if (solver%diffside%is_inward(idofdst)) then
                        xsrc(dst, k, i, j + 1) = xsrc(dst, k, i, j + 1) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do

                    do idof = 1, solver%dirside%dof
                      solrad = xedir(src - 1, k, i + i1 - sun%xinc, j)
                      if (solver%diffside%is_inward(idofdst)) then
                        xsrc(dst, k, i, j + 1) = xsrc(dst, k, i, j + 1) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do

                    do idof = 1, solver%dirside%dof
                      solrad = xedir(src - 1, k, i, j + i1 - sun%yinc)
                      if (solver%diffside%is_inward(idofdst)) then
                        xsrc(dst, k, i, j + 1) = xsrc(dst, k, i, j + 1) + solrad * dir2diff((dst) * C_dir%dof + src)
                      else
                        xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + solrad * dir2diff((dst) * C_dir%dof + src)
                      end if
                      src = src + 1
                    end do
                    dst = dst + 1
                  end do

                  if (ldebug) then
                    do src = 1, C_dir%dof
                      if (sum(dir2diff(src:C_dir%dof * C_diff%dof:C_dir%dof)) .gt. one .or. &
                          sum(dir2diff(src:C_dir%dof * C_diff%dof:C_dir%dof)) .lt. zero) &
                        print *, 'DEBUG Found dir2diff gt one:', src, '::', &
                        sum(dir2diff(src:C_dir%dof * C_diff%dof:C_dir%dof)), &
                        ':', dir2diff(src:C_dir%dof * C_diff%dof:C_dir%dof), &
                        '   :::::::     ', dir2diff
                    end do
                  end if

                end if ! 1D or Tenstream?
              end if ! if solar

            end do
          end do
        end do

        ! Ground Albedo reflecting direct radiation, the diffuse part is considered by the solver(Matrix)
        k = C_diff%ze
        do j = C_diff%ys, C_diff%ye
          do i = C_diff%xs, C_diff%xe
            do dst = 0, solver%difftop%dof - 1
              if (.not. solver%difftop%is_inward(dst + 1)) then
                do src = 0, solver%dirtop%dof - 1
                  xsrc(dst, k, i, j) = xsrc(dst, k, i, j) + xedir(src, k, i, j) * atm%albedo(i, j) &
                                       / real(solver%difftop%streams, ireals)
                end do
              end if
            end do
          end do
        end do

        call restoreVecPointer(C_dir%da, local_edir, xedir1d, xedir)
      end associate
    end subroutine

    subroutine set_buildings_reflection(solver, local_edir, buildings)
      class(t_solver), intent(in) :: solver
      type(tVec), intent(in) :: local_edir
      type(t_pprts_buildings), intent(in) :: buildings

      real(ireals), pointer, dimension(:, :, :, :) :: xedir => null()
      real(ireals), pointer, dimension(:) :: xedir1d => null()

      integer(iintegers) :: m, idx(4), isrc, idiff, dof_offset

      associate (                           &
          & B => buildings,                &
          & C => solver%C_dir)

        if (solution%lsolar_rad) then
          call getVecPointer(C%da, local_edir, xedir1d, xedir)
          do m = 1, size(B%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) - 1 + [C%zs, C%xs, C%ys]

            associate (k => idx(2), i => idx(3), j => idx(4))

              select case (idx(1))
              case (PPRTS_TOP_FACE)
                do idiff = 0, solver%difftop%dof - 1
                  if (.not. solver%difftop%is_inward(i1 + idiff)) xsrc(idiff, k, i, j) = 0
                end do
                do isrc = 0, solver%dirtop%dof - 1
                  do idiff = 0, solver%difftop%dof - 1
                    if (.not. solver%difftop%is_inward(i1 + idiff)) then ! eup on upper face
                      xsrc(idiff, k, i, j) = xsrc(idiff, k, i, j) &
                        & + xedir(isrc, k, i, j) * B%albedo(m) / real(solver%difftop%streams, ireals)
                    end if
                  end do
                end do

              case (PPRTS_BOT_FACE)
                do idiff = 0, solver%difftop%dof - 1
                  if (solver%difftop%is_inward(i1 + idiff)) xsrc(idiff, k + 1, i, j) = 0
                end do
                do isrc = 0, solver%dirtop%dof - 1
                  do idiff = 0, solver%difftop%dof - 1
                    if (solver%difftop%is_inward(i1 + idiff)) then ! edn on lower face
                      xsrc(idiff, k + 1, i, j) = xsrc(idiff, k + 1, i, j) &
                        & + xedir(isrc, k + 1, i, j) * B%albedo(m) / real(solver%difftop%streams, ireals)
                    end if
                  end do
                end do

              case (PPRTS_LEFT_FACE)
                dof_offset = solver%difftop%dof
                do idiff = 0, solver%diffside%dof - 1
                  if (.not. solver%diffside%is_inward(i1 + idiff)) xsrc(dof_offset + idiff, k, i, j) = 0
                end do

                do isrc = solver%dirtop%dof, solver%dirtop%dof + solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof - 1
                    if (.not. solver%diffside%is_inward(i1 + idiff)) then ! e_left on left face
                      xsrc(dof_offset + idiff, k, i, j) = xsrc(dof_offset + idiff, k, i, j) &
                        & + xedir(isrc, k, i, j) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end do

              case (PPRTS_RIGHT_FACE)
                dof_offset = solver%difftop%dof
                do idiff = 0, solver%diffside%dof - 1
                  if (solver%diffside%is_inward(i1 + idiff)) xsrc(dof_offset + idiff, k, i + 1, j) = 0
                end do

                do isrc = solver%dirtop%dof, solver%dirtop%dof + solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof - 1
                    if (solver%diffside%is_inward(i1 + idiff)) then ! e_right on right face
                      xsrc(dof_offset + idiff, k, i + 1, j) = xsrc(dof_offset + idiff, k, i + 1, j) &
                        & + xedir(isrc, k, i + 1, j) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end do

              case (PPRTS_REAR_FACE)
                dof_offset = solver%difftop%dof + solver%diffside%dof
                do idiff = 0, solver%diffside%dof - 1
                  if (.not. solver%diffside%is_inward(i1 + idiff)) xsrc(dof_offset + idiff, k, i, j) = 0
                end do

                do isrc = solver%dirtop%dof + solver%dirside%dof, solver%dirtop%dof + 2 * solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof - 1
                    if (.not. solver%diffside%is_inward(i1 + idiff)) then ! e_backward
                      xsrc(dof_offset + idiff, k, i, j) = xsrc(dof_offset + idiff, k, i, j) &
                        & + xedir(isrc, k, i, j) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end do

              case (PPRTS_FRONT_FACE)
                dof_offset = solver%difftop%dof + solver%diffside%dof
                do idiff = 0, solver%diffside%dof - 1
                  if (solver%diffside%is_inward(i1 + idiff)) xsrc(dof_offset + idiff, k, i, j + 1) = 0
                end do

                do isrc = solver%dirtop%dof + solver%dirside%dof, solver%dirtop%dof + 2 * solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof - 1
                    if (solver%diffside%is_inward(i1 + idiff)) then ! e_forward
                      xsrc(dof_offset + idiff, k, i, j + 1) = xsrc(dof_offset + idiff, k, i, j + 1) &
                        & + xedir(isrc, k, i, j + 1) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    end if
                  end do
                end do

              case default
                call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1) + 1))
              end select
            end associate
          end do ! loop over building faces
          call restoreVecPointer(C%da, local_edir, xedir1d, xedir)
        end if ! is_solar
      end associate
    end subroutine

    subroutine set_buildings_emission(solver, buildings, ierr)
      class(t_solver), intent(in) :: solver
      type(t_pprts_buildings), intent(in) :: buildings

      integer(iintegers) :: m, idx(4), idiff, dof_offset, ak
      real(ireals) :: emis, Az, Ax, Ay
      integer(mpiint) :: ierr

      if (.not. allocated(buildings%planck)) then
        ierr = 1
        call CHKERR(ierr, 'tried to set building emissions but planck data is not allocated')
      end if

      associate (atm => solver%atm, B => buildings, C => solver%C_diff)

        Az = atm%dx * atm%dy / real(solver%difftop%area_divider, ireals)

        do m = 1, size(B%iface)

          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) - 1 + [C%zs, C%xs, C%ys]

          associate (k => idx(2), i => idx(3), j => idx(4))

            ak = atmk(atm, k)
            Ax = atm%dy * atm%dz(ak, i, j) / real(solver%diffside%area_divider, ireals)
            Ay = atm%dx * atm%dz(ak, i, j) / real(solver%diffside%area_divider, ireals)

            emis = pi * B%planck(m) * (one - B%albedo(m))

            select case (idx(1))
            case (PPRTS_TOP_FACE)
              do idiff = 0, solver%difftop%dof - 1
                if (.not. solver%difftop%is_inward(i1 + idiff)) then ! eup on upper face
                  xsrc(idiff, k, i, j) = Az * emis / real(solver%difftop%streams, ireals)
                end if
              end do

            case (PPRTS_BOT_FACE)
              do idiff = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + idiff)) then ! edn on lower face
                  xsrc(idiff, k + 1, i, j) = Az * emis / real(solver%difftop%streams, ireals)
                end if
              end do

            case (PPRTS_LEFT_FACE)
              dof_offset = solver%difftop%dof
              do idiff = 0, solver%diffside%dof - 1
                if (.not. solver%diffside%is_inward(i1 + idiff)) then ! e_left on left face
                  xsrc(dof_offset + idiff, k, i, j) = Ax * emis / real(solver%diffside%streams, ireals)
                end if
              end do

            case (PPRTS_RIGHT_FACE)
              dof_offset = solver%difftop%dof
              do idiff = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idiff)) then ! e_right on right face
                  xsrc(dof_offset + idiff, k, i + 1, j) = Ax * emis / real(solver%diffside%streams, ireals)
                end if
              end do

            case (PPRTS_REAR_FACE)
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idiff = 0, solver%diffside%dof - 1
                if (.not. solver%diffside%is_inward(i1 + idiff)) then ! e_backward
                  xsrc(dof_offset + idiff, k, i, j) = Ay * emis / real(solver%diffside%streams, ireals)
                end if
              end do

            case (PPRTS_FRONT_FACE)
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idiff = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idiff)) then ! e_forward
                  xsrc(dof_offset + idiff, k, i, j + 1) = Ay * emis / real(solver%diffside%streams, ireals)
                end if
              end do

            case default
              call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1) + 1))
            end select
          end associate
        end do
      end associate

    end subroutine

  end subroutine setup_b

  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(OPP, kabs, ksca, g, dz, dx, ldir, coeff, &
                       angles, lswitch_east, lswitch_north, opt_vertices)
    class(t_optprop_cube), intent(in) :: OPP
    real(ireals), intent(in) :: kabs, ksca, g, dz, dx
    logical, intent(in) :: ldir
    real(irealLUT), intent(out) :: coeff(:)

    real(irealLUT), intent(in), optional :: angles(2)
    logical, intent(in), optional :: lswitch_east, lswitch_north
    real(ireals), intent(in), optional :: opt_vertices(:)

    real(irealLUT), save :: coeff_cache(1000)
    logical :: lcurrent

    real(irealLUT) :: aspect_zx, tauz, w0
    integer(mpiint) :: ierr

    call check_cache(lcurrent)
    if (lcurrent) then
      coeff = coeff_cache(1:size(coeff))
      return
    end if

    aspect_zx = real(dz / dx, irealLUT)
    w0 = real(ksca / max(kabs + ksca, epsilon(kabs)), irealLUT)
    tauz = real((kabs + ksca) * dz, irealLUT)

    if (present(angles)) then
      aspect_zx = max(OPP%dev%dirconfig%dims(3)%vrange(1), aspect_zx)
      tauz = max(OPP%dev%dirconfig%dims(1)%vrange(1), &
           & min(OPP%dev%dirconfig%dims(1)%vrange(2), tauz))
      w0 = max(OPP%dev%dirconfig%dims(2)%vrange(1), &
           & min(OPP%dev%dirconfig%dims(2)%vrange(2), w0))
    else
      aspect_zx = max(OPP%dev%diffconfig%dims(3)%vrange(1), aspect_zx)
      tauz = max(OPP%dev%diffconfig%dims(1)%vrange(1), &
           & min(OPP%dev%diffconfig%dims(1)%vrange(2), tauz))
      w0 = max(OPP%dev%diffconfig%dims(2)%vrange(1), &
           & min(OPP%dev%diffconfig%dims(2)%vrange(2), w0))
    end if

    call OPP%get_coeff(tauz, w0, real(g, irealLUT), aspect_zx, ldir, coeff, ierr, &
                       angles, lswitch_east, lswitch_north, opt_vertices); call CHKERR(ierr)

    coeff_cache(1:size(coeff)) = coeff

  contains
    subroutine check_cache(lcurrent)
      logical, intent(out) :: lcurrent

      logical, save :: lpresent_angles = .false.
      logical, save :: c_ldir = .false., c_lswitch_east = .false., c_lswitch_north = .false.
      real(ireals), save :: c_kabs = -1, c_ksca = -1, c_g = -1, c_dz = -1, c_vertices(24) = -1
      real(irealLUT), save :: c_angles(2) = -1

      real(ireals), parameter :: cache_limit = 1e-3 ! rel change cache limit
      real(irealLUT), parameter :: cache_limit2 = cache_limit ! rel change cache limit

      logical, parameter :: lenable_cache = .true.

      if (.not. lenable_cache) then
        lcurrent = .false.
        return
      end if

      if (.not. rel_approx(c_kabs, kabs, cache_limit)) goto 99
      if (.not. rel_approx(c_ksca, ksca, cache_limit)) goto 99
      if (.not. rel_approx(c_g, g, cache_limit)) goto 99
      if (.not. rel_approx(c_dz, dz, cache_limit)) goto 99
      if (lpresent_angles .neqv. present(angles)) goto 99

      if (present(opt_vertices)) then
        if (any(.not. rel_approx(c_vertices, opt_vertices, cache_limit))) goto 99
      end if

      if (present(angles)) then
        if (any(.not. rel_approx(c_angles, angles, cache_limit2))) goto 99
      end if

      if (c_ldir .neqv. ldir) goto 99
      if (present(lswitch_east)) then
        if (c_lswitch_east .neqv. lswitch_east) goto 99
      end if
      if (present(lswitch_north)) then
        if (c_lswitch_north .neqv. lswitch_north) goto 99
      end if
      lcurrent = .true.
      !print *,'hit cache'//new_line(''),  &
      !  & 'kabs', c_kabs, kabs,new_line(''), &
      !  & 'ksca', c_ksca, ksca,new_line(''), &
      !  & 'g',    c_g   , g   ,new_line(''), &
      !  & 'dz',   c_dz  , dz  ,new_line(''), &
      !  & 'ldir', c_ldir, ldir
      !call CHKERR(1_mpiint, 'DEBUG')
      return

99    continue ! update sample pts and go on to compute them
      !print *,'missed cache'//new_line(''),  &
      !  & 'kabs', c_kabs, kabs,new_line(''), &
      !  & 'ksca', c_ksca, ksca,new_line(''), &
      !  & 'g',    c_g   , g   ,new_line(''), &
      !  & 'dz',   c_dz  , dz  ,new_line(''), &
      !  & 'ldir', c_ldir, ldir

      lcurrent = .false.
      c_kabs = kabs
      c_ksca = ksca
      c_g = g
      c_dz = dz
      c_ldir = ldir
      lpresent_angles = present(angles)
      if (present(angles)) c_angles = angles
      if (present(opt_vertices)) c_vertices = opt_vertices
      if (present(lswitch_east)) c_lswitch_east = lswitch_east
      if (present(lswitch_north)) c_lswitch_north = lswitch_north

    end subroutine
  end subroutine

  !> @brief nca wrapper to call NCA of Carolin Klinger
  !> @details This is supposed to work on a 1D solution which has to be calculated beforehand
  !> \n the wrapper copies fluxes and optical properties on one halo and then gives that to NCA
  !> \n the result is the 3D approximation of the absorption, considering neighbouring information
  subroutine nca_wrapper(solver, ediff, abso)
    use m_ts_nca, only: ts_nca
    class(t_solver) :: solver
    type(tVec) :: ediff, abso
    type(tVec) :: gnca ! global nca vector
    type(tVec) :: lnca ! local nca vector with ghost values -- in dimension 0 and 1 are fluxes followed by dz,planck,kabs

    real(ireals), pointer, dimension(:, :, :, :) :: xv => null()
    real(ireals), pointer, dimension(:) :: xv1d => null()
    real(ireals), pointer, dimension(:, :, :, :) :: xvlnca => null(), xvgnca => null()
    real(ireals), pointer, dimension(:) :: xvlnca1d => null(), xvgnca1d => null()
    real(ireals), pointer, dimension(:, :, :, :) :: xhr => null()
    real(ireals), pointer, dimension(:) :: xhr1d => null()
    integer(iintegers) :: k

    integer(iintegers), parameter :: E_up = i0, E_dn = i1, idz = i2, iplanck = i3, ikabs = i4, ihr = i5
    integer(mpiint) :: ierr

    associate (atm => solver%atm, &
               C_one => solver%C_one, &
               C_diff => solver%C_diff)

      !call CHKERR(1_mpiint, 'nca_wrapper not implemented')
      !print *, 'DEBUG: Stupid print statement to prevent unused compiler warnings', ediff, abso
      if (C_diff%dof .lt. i6) call CHKERR(1_mpiint, 'For NCA, need a solver with at least 6 diffuse streams to copy over some data')

      ! put additional values into a local ediff vec .. TODO: this is a rather dirty hack but is straightforward

      ! get ghost values for dz, planck, kabs and fluxes, ready to give it to NCA
      call DMGetGlobalVector(solver%C_diff%da, gnca, ierr); call CHKERR(ierr)

      call getVecPointer(solver%C_diff%da, gnca, xvgnca1d, xvgnca)
      xvgnca(idz, C_diff%zs:C_diff%ze - 1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = atm%dz
      xvgnca(iplanck, C_diff%zs:C_diff%ze, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = atm%planck
      xvgnca(ikabs, C_diff%zs:C_diff%ze - 1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = atm%kabs

      ! Copy Edn and Eup to local convenience vector
      call getVecPointer(C_diff%da, ediff, xv1d, xv)
      xvgnca(E_up, :, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv(E_up, :, :, :)
      xvgnca(E_dn, :, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv(E_dn, :, :, :)
      call restoreVecPointer(C_diff%da, ediff, xv1d, xv)

      call restoreVecPointer(C_diff%da, gnca, xvgnca1d, xvgnca)

      ! retrieve ghost values into l(ocal) nca vec
      call DMGetLocalVector(C_diff%da, lnca, ierr); call CHKERR(ierr)
      call VecSet(lnca, zero, ierr); call CHKERR(ierr)

      call DMGlobalToLocalBegin(C_diff%da, gnca, ADD_VALUES, lnca, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd(C_diff%da, gnca, ADD_VALUES, lnca, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(C_diff%da, gnca, ierr); call CHKERR(ierr)

      ! call NCA
      call getVecPointer(C_diff%da, lnca, xvlnca1d, xvlnca)

      call ts_nca(atm%dx, atm%dy, &
                  xvlnca(idz, :, :, :), &
                  xvlnca(iplanck, :, :, :), &
                  xvlnca(ikabs, :, :, :), &
                  xvlnca(E_dn, :, :, :), &
                  xvlnca(E_up, :, :, :), &
                  xvlnca(ihr, :, :, :))

      ! return absorption
      call getVecPointer(C_one%da, abso, xhr1d, xhr)

      do k = C_one%zs, C_one%ze
        xhr(i0, k, :, :) = xvlnca(ihr, k, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) / &
                           xvlnca(idz, k, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye)
      end do
      call restoreVecPointer(C_one%da, abso, xhr1d, xhr)

      !return convenience vector that holds optical properties
      call restoreVecPointer(C_diff%da, lnca, xvlnca1d, xvlnca)
      call DMRestoreLocalVector(C_diff%da, lnca, ierr); call CHKERR(ierr)
    end associate
  end subroutine

  !> @brief build diffuse radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_diff_coeff(solver, A, C)
    class(t_solver) :: solver
    type(tMat) :: A
    type(t_coord) :: C

    integer(iintegers) :: i, j, k
    integer(mpiint) :: ierr

    if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid, 'Setting coefficients for diffuse Light'

    call MatZeroEntries(A, ierr); call CHKERR(ierr)
    call mat_set_diagonal(A)

    do j = C%ys, C%ye
      do i = C%xs, C%xe
        do k = C%zs, C%ze - 1

          if (solver%atm%l1d(atmk(solver%atm, k))) then
            call set_eddington_coeff(solver%atm, A, k, i, j)
          else
            call set_pprts_coeff(solver, C, A, k, i, j, ierr); call CHKERR(ierr)
          end if

        end do
      end do
    end do

    call set_albedo_coeff(solver, C, A)

    if (solver%myid .eq. 0 .and. ldebug) print *, solver%myid, 'Final diffuse Matrix Assembly'
    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, "-show_Mdiff", ierr); call CHKERR(ierr)
  contains
    subroutine set_pprts_coeff(solver, C, A, k, i, j, ierr)
      class(t_solver) :: solver
      type(t_coord), intent(in) :: C
      type(tMat), intent(inout) :: A
      integer(iintegers), intent(in) :: k, i, j
      integer(mpiint), intent(out) :: ierr

      MatStencil :: row(4, 0:C%dof - 1), col(4, 0:C%dof - 1)
      real(ireals) :: norm

      integer(iintegers) :: dst, src, idof

      src = 0
      do idof = 1, solver%difftop%dof
        if (solver%difftop%is_inward(idof)) then
          col(MatStencil_j, src) = i
          col(MatStencil_k, src) = j
          col(MatStencil_i, src) = k
          col(MatStencil_c, src) = src
        else
          col(MatStencil_j, src) = i
          col(MatStencil_k, src) = j
          col(MatStencil_i, src) = k + 1
          col(MatStencil_c, src) = src
        end if
        src = src + 1
      end do

      do idof = 1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          col(MatStencil_j, src) = i
          col(MatStencil_k, src) = j
          col(MatStencil_i, src) = k
          col(MatStencil_c, src) = src
        else
          col(MatStencil_j, src) = i + 1
          col(MatStencil_k, src) = j
          col(MatStencil_i, src) = k
          col(MatStencil_c, src) = src
        end if
        src = src + 1
      end do

      do idof = 1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          col(MatStencil_j, src) = i
          col(MatStencil_k, src) = j
          col(MatStencil_i, src) = k
          col(MatStencil_c, src) = src
        else
          col(MatStencil_j, src) = i
          col(MatStencil_k, src) = j + 1
          col(MatStencil_i, src) = k
          col(MatStencil_c, src) = src
        end if
        src = src + 1
      end do

      dst = 0
      do idof = 1, solver%difftop%dof
        if (solver%difftop%is_inward(idof)) then
          row(MatStencil_j, dst) = i
          row(MatStencil_k, dst) = j
          row(MatStencil_i, dst) = k + 1
          row(MatStencil_c, dst) = dst
        else
          row(MatStencil_j, dst) = i
          row(MatStencil_k, dst) = j
          row(MatStencil_i, dst) = k
          row(MatStencil_c, dst) = dst
        end if
        dst = dst + 1
      end do

      do idof = 1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          row(MatStencil_j, dst) = i + 1
          row(MatStencil_k, dst) = j
          row(MatStencil_i, dst) = k
          row(MatStencil_c, dst) = dst
        else
          row(MatStencil_j, dst) = i
          row(MatStencil_k, dst) = j
          row(MatStencil_i, dst) = k
          row(MatStencil_c, dst) = dst
        end if
        dst = dst + 1
      end do

      do idof = 1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          row(MatStencil_j, dst) = i
          row(MatStencil_k, dst) = j + 1
          row(MatStencil_i, dst) = k
          row(MatStencil_c, dst) = dst
        else
          row(MatStencil_j, dst) = i
          row(MatStencil_k, dst) = j
          row(MatStencil_i, dst) = k
          row(MatStencil_c, dst) = dst
        end if
        dst = dst + 1
      end do

      call MatSetValuesStencil(A, C%dof, row, C%dof, col, &
        & real(-solver%diff2diff(:, k, i, j), ireals), INSERT_VALUES, ierr); call CHKERR(ierr)

      if (ldebug) then
        do src = 1, C%dof
          norm = sum(solver%diff2diff(src:C%dof**2:C%dof, k, i, j))
          if (norm .gt. one + 10._ireals * epsilon(norm)) then ! could renormalize
            if (norm .gt. one + 10._ireals * sqrt(epsilon(norm))) then ! fatally off
              print *, 'diffuse sum(src==', src, ') gt one',&
                & 'norm', norm, &
                & '=>', solver%diff2diff(src:C%dof**2:C%dof, k, i, j)
              print *, 'get_coeff', solver%atm%kabs(atmk(solver%atm, k), i, j), &
                & solver%atm%ksca(atmk(solver%atm, k), i, j), &
                & solver%atm%g(atmk(solver%atm, k), i, j), &
                & solver%atm%dz(atmk(solver%atm, k), i, j), &
                & .false., solver%atm%l1d(atmk(solver%atm, k)), &
                & '=> all coeff', solver%diff2diff(:, k, i, j)
              call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
            end if ! fatal
          end if ! could renormalize
        end do
      end if

    end subroutine

    subroutine set_eddington_coeff(atm, A, k, i, j)
      type(t_atmosphere) :: atm
      type(tMat), intent(inout) :: A
      integer(iintegers), intent(in) :: k, i, j
      integer(iintegers) :: src, dst, idof

      MatStencil :: row(4, 0:solver%difftop%dof - 1), col(4, 0:solver%difftop%dof - 1)
      real(ireals) :: v(solver%difftop%dof**2)

      do idof = 1, solver%difftop%dof
        src = idof - 1
        if (solver%difftop%is_inward(idof)) then
          col(MatStencil_j, src) = i
          col(MatStencil_k, src) = j
          col(MatStencil_i, src) = k
          col(MatStencil_c, src) = src
        else
          col(MatStencil_j, src) = i
          col(MatStencil_k, src) = j
          col(MatStencil_i, src) = k + 1
          col(MatStencil_c, src) = src
        end if
      end do

      do idof = 1, solver%difftop%dof
        dst = idof - 1
        if (solver%difftop%is_inward(idof)) then
          row(MatStencil_j, dst) = i
          row(MatStencil_k, dst) = j
          row(MatStencil_i, dst) = k + 1
          row(MatStencil_c, dst) = dst
        else
          row(MatStencil_j, dst) = i
          row(MatStencil_k, dst) = j
          row(MatStencil_i, dst) = k
          row(MatStencil_c, dst) = dst
        end if
      end do

      ! for each destination, find all transmission coeffs
      v(:) = zero
      do dst = 0, solver%difftop%dof - 1
        do src = 0, solver%difftop%dof - 1
          if (col(MatStencil_i, src) .eq. row(MatStencil_i, dst)) then ! for reflection, has to be the same k layers

            if (src .ne. inv_dof(dst)) cycle ! in 1D has to be the inverse stream
            v(i1 + dst * solver%difftop%dof + src) = atm%a12(atmk(atm, k), i, j) ! TODO: should be using a21/a22 for upward streams?
            !print *,i,j,k,'setting r ',toStr(src)//' (k='//toStr(col(MatStencil_i,src))//') ->', &
            !  toStr(dst)//' (k='//toStr(row(MatStencil_i,dst))//')'// &
            !  ':', i1 + dst*solver%difftop%dof + src, v(i1 + dst*solver%difftop%dof + src), &
            !  'invdof',src,dst,inv_dof(dst)
          else

            if (src .ne. dst) cycle ! in 1D has to be the same
            v(i1 + dst * solver%difftop%dof + src) = atm%a11(atmk(atm, k), i, j) ! TODO: should be using a21/a22 for upward streams?

            !print *,i,j,k,'setting t ',toStr(src)//' (k='//toStr(col(MatStencil_i,src))//') ->', &
            !  toStr(dst)//' (k='//toStr(row(MatStencil_i,dst))//') :', i1 + dst*solver%difftop%dof + src, v(i1 + dst*solver%difftop%dof + src)
          end if ! which k-lev
        end do
      end do

      call MatSetValuesStencil(A, solver%difftop%dof, row, solver%difftop%dof, col, -v, INSERT_VALUES, ierr); call CHKERR(ierr)
    end subroutine
    pure function inv_dof(dof) ! returns the dof that is the same stream but the opposite direction
      integer(iintegers), intent(in) :: dof
      integer(iintegers) :: inv_dof, inc
      if (solver%difftop%is_inward(1)) then ! starting with downward streams
        inc = 1
      else
        inc = -1
      end if
      if (solver%difftop%is_inward(i1 + dof)) then ! downward stream
        inv_dof = dof + inc
      else
        inv_dof = dof - inc
      end if
    end function

    !> @brief insert lower boundary condition, i.e. diffuse reflection of downward radiation
    subroutine set_albedo_coeff(solver, C, A)
      class(t_solver), intent(in) :: solver
      type(t_coord), intent(in) :: C
      type(tMat), intent(inout) :: A

      MatStencil :: row(4, 1), col(4, 1)
      integer(iintegers) :: i, j, src, dst

      ! Set surface albedo values
      col(MatStencil_i, i1) = C%ze
      row(MatStencil_i, i1) = C%ze

      do j = C%ys, C%ye
        row(MatStencil_k, i1) = j
        col(MatStencil_k, i1) = j

        do i = C%xs, C%xe
          row(MatStencil_j, i1) = i
          col(MatStencil_j, i1) = i

          do dst = 1, solver%difftop%dof
            if (.not. solver%difftop%is_inward(dst)) then
              do src = 1, solver%difftop%dof
                if (solver%difftop%is_inward(src)) then
                  col(MatStencil_c, i1) = src - 1
                  row(MatStencil_c, i1) = dst - 1
                  !print *,solver%myid, 'i '//toStr(i)//' j '//toStr(j), &
                  !  & ' Setting albedo for dst '//toStr(row)//' src '//toStr(col), &
                  !  & ':', -solver%atm%albedo(i,j) / real(solver%difftop%streams, ireals)
                  call MatSetValuesStencil(A, i1, row, i1, col, &
                                           [-solver%atm%albedo(i, j) / real(solver%difftop%streams, ireals)], &
                                           INSERT_VALUES, ierr); call CHKERR(ierr)
                end if
              end do
            end if
          end do

        end do
      end do
    end subroutine

  end subroutine

  subroutine pprts_get_result(solver, redn, reup, rabso, redir, opt_solution_uid, opt_buildings)
    class(t_solver) :: solver
    real(ireals), dimension(:, :, :), intent(inout), allocatable :: redn, reup, rabso
    real(ireals), dimension(:, :, :), intent(inout), allocatable, optional :: redir
    integer(iintegers), optional, intent(in) :: opt_solution_uid
    type(t_pprts_buildings), optional, intent(inout) :: opt_buildings

    integer(iintegers) :: uid, iside
    real(ireals), pointer :: x1d(:) => null(), x4d(:, :, :, :) => null()
    integer(mpiint) :: ierr

    uid = get_solution_uid(solver%solutions, opt_solution_uid)

    associate (solution => solver%solutions(uid))
      if (solution%lsolar_rad) then
        call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-pprts_show_edir", ierr); call CHKERR(ierr)
      end if
      call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-pprts_show_ediff", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, "-pprts_show_abso", ierr); call CHKERR(ierr)

      if (ldebug .and. solver%myid .eq. 0) print *, 'calling pprts_get_result', present(redir), 'for uid', uid

      if (solution%lchanged) &
        & call CHKERR(1_mpiint, 'tried to get results from unrestored solution -- call restore_solution first')

      if (present(opt_solution_uid)) then
        if (.not. solution%lWm2_diff) &
          & call CHKERR(1_mpiint, 'solution vecs for diffuse radiation are not in W/m2 ... this is not what I expected')
      end if

      call alloc_or_check_size(solver%C_diff, redn, 'edn')
      call alloc_or_check_size(solver%C_diff, reup, 'eup')
      call alloc_or_check_size(solver%C_one, rabso, 'abso')

      if (present(redir)) then
        call alloc_or_check_size(solver%C_dir, redir, 'edir')

        if (.not. solution%lsolar_rad) then
          if (ldebug) then
            call CHKWARN(1_mpiint, 'Hey, You called pprts_get_result for uid '//toStr(uid)// &
                         ' and provided an array for direct radiation.'// &
                         ' However in this particular band we haven`t computed direct radiation.'// &
                         ' I will return with edir=0 but are you sure this is what you intended?')
          end if
          redir = zero
        else
          if (.not. solution%lWm2_dir) &
            & call CHKERR(1_mpiint, 'tried to get result from a result vector(dir) which is not in [W/m2]')

          call getVecPointer(solver%C_dir%da, solution%edir, x1d, x4d)
          ! average of direct radiation of all fluxes through top faces
          redir = sum(x4d(0:solver%dirtop%dof - 1, :, :, :), dim=1) / real(solver%dirtop%area_divider, ireals)
          call restoreVecPointer(solver%C_dir%da, solution%edir, x1d, x4d)

          if (ldebug) then
            if (solver%myid .eq. 0) print *, 'Edir vertically first column', redir(:, lbound(redir, 2), lbound(redir, 3))
            if (any(redir .lt. -one)) then
              print *, 'Found direct radiation smaller than 0 in dir result... that should not happen', minval(redir)
              call CHKERR(1_mpiint)
            end if
          end if
        end if
      end if

      if (.not. solution%lWm2_diff) &
        & call CHKERR(1_mpiint, 'tried to get result from a result vector(diff) which is not in [W/m2]')

      redn = zero
      reup = zero
      call getVecPointer(solver%C_diff%da, solution%ediff, x1d, x4d)
      do iside = 1, solver%difftop%dof
        if (solver%difftop%is_inward(iside)) then
          redn = redn + x4d(iside - 1, :, :, :) / real(solver%difftop%area_divider, ireals)
        else
          reup = reup + x4d(iside - 1, :, :, :) / real(solver%difftop%area_divider, ireals)
        end if
      end do
      call restoreVecPointer(solver%C_diff%da, solution%ediff, x1d, x4d)

      if (solver%myid .eq. 0 .and. ldebug .and. present(redir)) &
        print *, 'mean surface Edir', meanval(redir(ubound(redir, 1), :, :))
      if (solver%myid .eq. 0 .and. ldebug) print *, 'mean surface Edn', meanval(redn(ubound(redn, 1), :, :))
      if (solver%myid .eq. 0 .and. ldebug) print *, 'mean surface Eup', meanval(reup(ubound(reup, 1), :, :))

      if (ldebug .and. solution%lsolar_rad) then
        if (any(redn .lt. -one)) then
          print *, 'Found radiation smaller than 0 in edn result... that should not happen', minval(redn)
          call exit(1)
        end if
        if (any(reup .lt. -one)) then
          print *, 'Found radiation smaller than 0 in eup result... that should not happen', minval(reup)
          call exit(1)
        end if
      end if

      call getVecPointer(solver%C_one%da, solution%abso, x1d, x4d, readonly=.true.)
      rabso = x4d(i0, :, :, :)
      call restoreVecPointer(solver%C_one%da, solution%abso, x1d, x4d, readonly=.true.)

      if (present(opt_buildings)) then
        call fill_buildings()
      end if

      call dump_to_xdmf()

      if (solver%myid .eq. 0 .and. ldebug) print *, 'get_result done'

    end associate

  contains

    subroutine dump_to_xdmf()
      character(len=default_str_len) :: fname
      logical :: lflg

      fname = ''
      if (present(opt_buildings)) then
        call get_petsc_opt(solver%prefix, '-pprts_xdmf_buildings', fname, lflg, ierr); call CHKERR(ierr)
        if (lflg) then
          fname = trim(fname)//'.'//toStr(uid)
          call xdmf_pprts_buildings(solver, opt_buildings, fname, ierr, verbose=.true.); call CHKERR(ierr)
        end if
      end if

      associate (solution => solver%solutions(uid))
        call get_petsc_opt(solver%prefix, '-pprts_xdmf', fname, lflg, ierr); call CHKERR(ierr)
        if (lflg) then
          fname = trim(fname)//toStr(uid)
          call xdmf_pprts_srfc_flux(solver, fname, redn, reup, ierr, edir=redir, verbose=.true.); call CHKERR(ierr)
        end if
      end associate
    end subroutine

    subroutine fill_buildings
      integer(iintegers) :: m, idx(4), dof_offset, idof
      type(tVec) :: ledir, lediff

      associate (                               &
          & solution => solver%solutions(uid), &
          & B => opt_buildings,                &
          & C => solver%C_dir)

        if (solution%lsolar_rad) then
          if (.not. allocated(B%edir)) allocate (B%edir(size(B%iface)))

          call DMGetLocalVector(C%da, ledir, ierr); call CHKERR(ierr)
          call VecSet(ledir, zero, ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C%da, solution%edir, ADD_VALUES, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd(C%da, solution%edir, ADD_VALUES, ledir, ierr); call CHKERR(ierr)

          call getVecPointer(C%da, ledir, x1d, x4d, readonly=.true.)
          ! average of direct radiation of all fluxes through top faces
          !redir = sum(x4d(0:solver%dirtop%dof-1, :, :, :), dim=1) / real(solver%dirtop%area_divider, ireals)

          do m = 1, size(B%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) - 1 + [C%zs, C%xs, C%ys]

            associate (k => idx(2), i => idx(3), j => idx(4))

              select case (idx(1))
              case (PPRTS_TOP_FACE)
                dof_offset = 0
                B%edir(m) = sum(x4d(dof_offset:dof_offset + solver%dirtop%dof - 1, k, i, j)) &
                  & / real(solver%dirtop%area_divider, ireals)
              case (PPRTS_BOT_FACE)
                dof_offset = 0
                B%edir(m) = sum(x4d(dof_offset:dof_offset + solver%dirtop%dof - 1, k + 1, i, j)) &
                  & / real(solver%dirtop%area_divider, ireals)
              case (PPRTS_LEFT_FACE)
                dof_offset = solver%dirtop%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset + solver%dirside%dof - 1, k, i, j)) &
                  & / real(solver%dirside%area_divider, ireals)
              case (PPRTS_RIGHT_FACE)
                dof_offset = solver%dirtop%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset + solver%dirside%dof - 1, k, i + 1, j)) &
                  & / real(solver%dirside%area_divider, ireals)
              case (PPRTS_REAR_FACE)
                dof_offset = solver%dirtop%dof + solver%dirside%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset + solver%dirside%dof - 1, k, i, j)) &
                  & / real(solver%dirside%area_divider, ireals)
              case (PPRTS_FRONT_FACE)
                dof_offset = solver%dirtop%dof + solver%dirside%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset + solver%dirside%dof - 1, k, i, j + 1)) &
                  & / real(solver%dirside%area_divider, ireals)
              case default
                call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1) + 1))
              end select
            end associate
          end do

          call restoreVecPointer(C%da, ledir, x1d, x4d, readonly=.true.)
          call DMRestoreLocalVector(C%da, ledir, ierr); call CHKERR(ierr)
        end if
      end associate

      associate (                               &
          & solution => solver%solutions(uid), &
          & B => opt_buildings,                &
          & C => solver%C_diff)
        if (.not. allocated(B%incoming)) allocate (B%incoming(size(B%iface)))
        if (.not. allocated(B%outgoing)) allocate (B%outgoing(size(B%iface)))

        call DMGetLocalVector(C%da, lediff, ierr); call CHKERR(ierr)
        call VecSet(lediff, zero, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(C%da, solution%ediff, ADD_VALUES, lediff, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd(C%da, solution%ediff, ADD_VALUES, lediff, ierr); call CHKERR(ierr)
        call getVecPointer(C%da, lediff, x1d, x4d, readonly=.true.)

        do m = 1, size(B%iface)
          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) - 1 + [C%zs, C%xs, C%ys]

          associate (k => idx(2), i => idx(3), j => idx(4))

            B%incoming(m) = zero
            B%outgoing(m) = zero

            select case (idx(1))

            case (PPRTS_TOP_FACE)
              do idof = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + idof)) then ! e_down
                  B%incoming(m) = B%incoming(m) + x4d(idof, k, i, j)
                else
                  B%outgoing(m) = B%outgoing(m) + x4d(idof, k, i, j)
                end if
              end do
              B%incoming(m) = B%incoming(m) / real(solver%difftop%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%difftop%area_divider, ireals)

            case (PPRTS_BOT_FACE)
              do idof = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + idof)) then ! e_down
                  B%outgoing(m) = B%outgoing(m) + x4d(idof, k + 1, i, j)
                else
                  B%incoming(m) = B%incoming(m) + x4d(idof, k + 1, i, j)
                end if
              end do
              B%incoming(m) = B%incoming(m) / real(solver%difftop%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%difftop%area_divider, ireals)

            case (PPRTS_LEFT_FACE)
              dof_offset = solver%difftop%dof
              do idof = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idof)) then ! e_right
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i, j)
                else
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i, j)
                end if
              end do
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case (PPRTS_RIGHT_FACE)
              dof_offset = solver%difftop%dof
              do idof = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idof)) then ! e_right
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i + 1, j)
                else
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i + 1, j)
                end if
              end do
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case (PPRTS_REAR_FACE)
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idof = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idof)) then ! e_forward
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i, j)
                else
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i, j)
                end if
              end do
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case (PPRTS_FRONT_FACE)
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idof = 0, solver%diffside%dof - 1
                if (solver%diffside%is_inward(i1 + idof)) then ! e_forward
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i, j + 1)
                else
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i, j + 1)
                end if
              end do
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case default
              call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1) + 1))
            end select
          end associate
        end do

        call restoreVecPointer(C%da, solution%ediff, x1d, x4d, readonly=.true.)
        call DMRestoreLocalVector(C%da, lediff, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    subroutine alloc_or_check_size(C, arr, varname)
      type(t_coord), intent(in) :: C
      real(ireals), intent(inout), allocatable :: arr(:, :, :)
      character(len=*) :: varname
      if (allocated(arr)) then
        if (.not. all(shape(arr) .eq. [C%zm, C%xm, C%ym])) then
          call CHKERR(1_mpiint, 'the shape of result array which you provided does not conform to output size. '// &
            & 'Either call with unallocated object or make sure it has the correct size! var='//trim(varname)// &
            & ' you provided arr with shape: '//toStr(shape(arr))//' but expect '//toStr([C%zm, C%xm, C%ym]))
        end if
      else
        allocate (arr(C%zm, C%xm, C%ym))
      end if
    end subroutine
  end subroutine

  subroutine pprts_get_result_toZero(solver, gedn, geup, gabso, gedir, opt_solution_uid)
    ! after solving equations -- retrieve the results for edir,edn,eup and absorption
    ! only zeroth node gets the results back.
    class(t_solver) :: solver

    real(ireals), intent(inout), dimension(:, :, :), allocatable :: gedn
    real(ireals), intent(inout), dimension(:, :, :), allocatable :: geup
    real(ireals), intent(inout), dimension(:, :, :), allocatable :: gabso
    real(ireals), intent(inout), dimension(:, :, :), allocatable, optional :: gedir
    integer(iintegers), optional, intent(in) :: opt_solution_uid

    real(ireals), allocatable, dimension(:, :, :) :: redir, redn, reup, rabso
    integer(mpiint) :: ierr

    call PetscLogEventBegin(solver%logs%scatter_to_Zero, ierr); call CHKERR(ierr)
    if (solver%myid .eq. 0) then
      call check_arr_size(solver%C_one1, gedn)
      call check_arr_size(solver%C_one1, geup)
      call check_arr_size(solver%C_one, gabso)
      if (present(gedir)) call check_arr_size(solver%C_one1, gedir)
    end if

    if (present(gedir)) then
      call pprts_get_result(solver, redn, reup, rabso, redir=redir, opt_solution_uid=opt_solution_uid)
      call gather_all_toZero(solver%C_one1, redir, gedir)
    else
      call pprts_get_result(solver, redn, reup, rabso, opt_solution_uid=opt_solution_uid)
    end if

    call gather_all_toZero(solver%C_one1, redn, gedn)
    call gather_all_toZero(solver%C_one1, reup, geup)
    call gather_all_toZero(solver%C_one, rabso, gabso)

    if (solver%myid .eq. 0 .and. ldebug) then
      print *, 'Retrieving results:'
      if (present(gedir)) print *, sum(gedir) / size(gedir)
      print *, sum(gedn) / size(gedn)
      print *, sum(geup) / size(geup)
      print *, sum(gabso) / size(gabso)
    end if
    call PetscLogEventEnd(solver%logs%scatter_to_Zero, ierr); call CHKERR(ierr)

  contains
    subroutine check_arr_size(C, inp)
      type(t_coord), intent(in) :: C
      real(ireals), intent(in), allocatable :: inp(:, :, :)
      if (allocated(inp)) then
        if (size(inp) .ne. C%dof * C%glob_zm * C%glob_xm * C%glob_ym) then
          print *, 'pprts_get_result_toZero was called with an already allocated output array but it has the wrong size.'
          print *, 'while I could just re-allocate it, this may not be just what you intended. Please do that yourself.'
          print *, 'Size of global Dimensions of simulation:', &
            & C%dof * C%glob_zm * C%glob_xm * C%glob_ym, &
            & '.vs. your input:', size(inp)
          call CHKERR(1_mpiint, 'pprts_get_result_toZero :: should not be called with already allocated array with wrong size')
        end if
      end if
    end subroutine
  end subroutine

  subroutine gather_all_toZero(C, inp, outp)
    type(t_coord), intent(in) :: C
    real(ireals), intent(in) :: inp(:, :, :) ! local array from get_result
    real(ireals), intent(inout), allocatable :: outp(:, :, :) ! global sized array on rank 0

    type(tVec) :: vec, lvec_on_zero

    integer(mpiint) :: myid, ierr
    call mpi_comm_rank(C%comm, myid, ierr)

    if (ldebug) then
      print *, myid, 'exchange_var', allocated(outp)
      print *, myid, 'exchange_var shape', shape(inp)
    end if

    call DMGetGlobalVector(C%da, vec, ierr); call CHKERR(ierr)
    call f90VecToPetsc(inp, C%da, vec)
    call petscGlobalVecToZero(vec, C%da, lvec_on_zero)
    call DMRestoreGlobalVector(C%da, vec, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      call petscVecToF90(lvec_on_zero, C%da, outp, only_on_rank0=.true.)
    end if
    call VecDestroy(lvec_on_zero, ierr); call CHKERR(ierr)
  end subroutine

  subroutine gather_all_to_all(C, inp, outp)
    type(t_coord), intent(in) :: C
    real(ireals), intent(in), allocatable :: inp(:, :, :) ! local array from get_result
    real(ireals), intent(inout), allocatable :: outp(:, :, :) ! global sized array on rank 0

    type(tVec) :: vec, lvec

    integer(mpiint) :: myid, ierr
    call mpi_comm_rank(C%comm, myid, ierr)

    if (ldebug) then
      print *, myid, 'exchange_var', allocated(inp), allocated(outp)
      print *, myid, 'exchange_var shape', shape(inp)
    end if

    call DMGetGlobalVector(C%da, vec, ierr); call CHKERR(ierr)
    call f90VecToPetsc(inp, C%da, vec)
    call petscGlobalVecToAll(vec, C%da, lvec)
    call DMRestoreGlobalVector(C%da, vec, ierr); call CHKERR(ierr)

    call petscVecToF90(lvec, C%da, outp)
    call VecDestroy(lvec, ierr); call CHKERR(ierr)
  end subroutine

  !> @brief single column tenstream solve from horizontally averaged optical properties
  recursive subroutine initial_guess_from_single_column_comp(solver, edirTOA, solution, ierr, opt_buildings)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), intent(inout) :: solution
    integer(mpiint), intent(out) :: ierr
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    class(t_solver), allocatable :: solver1d
    character(len=default_str_len) :: prefix

    integer(iintegers) :: k, idof
    integer(iintegers), parameter :: solution_uid = 0

    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g, planck
    real(ireals), allocatable, dimension(:, :) :: planck_srfc
    real(ireals), allocatable, dimension(:) :: dz

    real(ireals), pointer :: x1d(:) => null(), x4d(:, :, :, :) => null()
    real(ireals), pointer :: y1d(:) => null(), y4d(:, :, :, :) => null()

    ierr = 0
    if (solver%C_one%glob_xm .le. 5_iintegers .and. solver%C_one%glob_ym .le. 5_iintegers) return ! if domains are really small we dont do single column initializations

    if (len_trim(solver%solvername) .gt. 0) then
      prefix = trim(solver%solvername)//'_initial_guess_'
    else
      prefix = 'initial_guess_'
    end if

    call allocate_pprts_solver_from_commandline(solver1d, solver_to_str(solver), ierr, trim(prefix)); call CHKERR(ierr)

    allocate (dz(solver%C_one_atm%zs:solver%C_one_atm%ze))
    do k = solver%C_one_atm%zs, solver%C_one_atm%ze
      dz(k) = meanval(solver%atm%dz(k, :, :))
    end do

    call init_pprts(MPI_COMM_SELF, &
      & solver%C_one_atm%zm, &
      & i1, &
      & i1, &
      & solver%atm%dx, &
      & solver%atm%dy, &
      & solver%sun%sundir, &
      & solver1d, &
      & dz1d=dz, &
      & collapseindex=solver%atm%icollapse, &
      & solvername=prefix)

    associate ( &
        & Catm => solver1d%C_one_atm, &
        & Catm1 => solver1d%C_one_atm1 &
        )
      allocate (kabs(Catm%zs:Catm%ze, Catm%xs:Catm%xe, Catm%ys:Catm%ye))
      allocate (ksca(Catm%zs:Catm%ze, Catm%xs:Catm%xe, Catm%ys:Catm%ye))
      allocate (g(Catm%zs:Catm%ze, Catm%xs:Catm%xe, Catm%ys:Catm%ye))
      do k = Catm%zs, Catm%ze
        kabs(k, :, :) = meanval(solver%atm%kabs(k, :, :))
        ksca(k, :, :) = meanval(solver%atm%ksca(k, :, :))
        g(k, :, :) = meanval(solver%atm%g(k, :, :))
      end do
      if (.not. solution%lsolar_rad) then
        allocate (planck(Catm1%zs:Catm1%ze, Catm1%xs:Catm1%xe, Catm1%ys:Catm1%ye))
        allocate (planck_srfc(Catm1%xs:Catm1%xe, Catm1%ys:Catm1%ye))
        do k = Catm1%zs, Catm1%ze
          planck(k, :, :) = meanval(solver%atm%planck(k, :, :))
        end do
        planck_srfc(:, :) = meanval(solver%atm%Bsrfc)
      end if
    end associate

    if (solution%lsolar_rad) then
      call set_optical_properties(solver1d, meanval(solver%atm%albedo), kabs, ksca, g)
      call set_angles(solver1d, solver%sun%sundir)
    else
      call set_optical_properties(solver1d, meanval(solver%atm%albedo), kabs, ksca, g, &
        & planck=planck, planck_srfc=planck_srfc)
    end if

    call solve_pprts(solver1d, &
      & lthermal=.not. solution%lsolar_rad, &
      & lsolar=solution%lsolar_rad, &
      & edirTOA=edirTOA, &
      & opt_solution_uid=solution_uid)

    associate (solution1d => solver1d%solutions(solution_uid))
      if (solution%lsolar_rad) then
        call getVecPointer(solver1d%C_dir%da, solution1d%edir, x1d, x4d)
        call getVecPointer(solver%C_dir%da, solution%edir, y1d, y4d)

        do k = solver%C_dir%zs, solver%C_dir%ze
          do idof = 0, solver%C_dir%dof - 1
            y4d(idof, k, :, :) = meanval(x4d(idof, k, :, :))
          end do
        end do

        call restoreVecPointer(solver1d%C_dir%da, solution1d%edir, y1d, y4d)
        call restoreVecPointer(solver%C_dir%da, solution%edir, x1d, x4d)
      end if

      call getVecPointer(solver1d%C_diff%da, solution1d%ediff, x1d, x4d)
      call getVecPointer(solver%C_diff%da, solution%ediff, y1d, y4d)

      do k = solver%C_diff%zs, solver%C_diff%ze
        do idof = 0, solver%C_diff%dof - 1
          y4d(idof, k, :, :) = meanval(x4d(idof, k, :, :))
        end do
      end do

      call restoreVecPointer(solver%C_diff%da, solution%ediff, y1d, y4d)
      call restoreVecPointer(solver1d%C_diff%da, solution1d%ediff, x1d, x4d)

      solution%lWm2_dir = solution1d%lWm2_dir
      solution%lWm2_diff = solution1d%lWm2_diff
    end associate

    call destroy_pprts(solver1d)
  end subroutine
end module

