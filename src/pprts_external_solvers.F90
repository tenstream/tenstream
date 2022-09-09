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

module m_pprts_external_solvers
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, imp_iinteger, &
    & i0, i1, i2, i3, i4, zero, one

  use m_helper_functions, only: &
    & CHKERR, CHKWARN, &
    & get_petsc_opt, &
    & imp_allreduce_sum, &
    & ind_1d_to_nd, &
    & ind_nd_to_1d, &
    & meanval, &
    & ndarray_offsets, &
    & spherical_2_cartesian, &
    & toStr

  use m_petsc_helpers, only: &
    & f90VecToPetsc, &
    & gen_shared_scatter_ctx, &
    & gen_shared_subcomm, &
    & getVecPointer, &
    & restoreVecPointer

  use m_pprts_base, only: &
    & atmk, &
    & destroy_solution, &
    & interpolate_cell_values_to_vertices, &
    & prepare_solution, &
    & set_dmda_cell_coordinates, &
    & t_solver, &
    & t_state_container

  use m_schwarzschild, only: schwarzschild

  use m_twostream, only: delta_eddington_twostream, adding_delta_eddington_twostream

  use m_icon_plex_utils, only: create_2d_regular_plex, dmplex_2D_to_3D

  use m_plex_grid, only: &
    & destroy_plexgrid, &
    & print_dmplex, &
    & setup_abso_dmplex, &
    & setup_ediff_dmplex, &
    & setup_edir_dmplex, &
    & setup_plexgrid, &
    & t_plexgrid

  use m_plex2rayli, only: rayli_wrapper

  use m_pprts2plex, only: pprts_buildings_to_plex, find_face_idx_by_orientation, pprts_cell_to_plex_cell_idx

  use m_tenstr_disort, only: default_flx_computation

  use m_buildings, only: &
    & t_pprts_buildings, t_plex_buildings, &
    & clone_buildings, &
    & destroy_buildings, &
    & check_buildings_consistency, &
    & PPRTS_TOP_FACE, &
    & PPRTS_BOT_FACE, &
    & PPRTS_LEFT_FACE, &
    & PPRTS_REAR_FACE

  implicit none
  private
  public :: disort, twostream, schwarz, pprts_rayli_wrapper, destroy_rayli_info

  logical, parameter :: ldebug = .false.

  type t_rayli_info_buildings
    type(t_pprts_buildings), allocatable :: subcomm_buildings
    type(t_plex_buildings), allocatable :: plex_buildings
    integer(iintegers) :: Nglob ! global number of building entries (0 if not a subcomm master)
    type(tVecScatter) :: ctx_albedo
  end type

  type t_rayli_info
    type(t_state_container) :: plex_solution
    type(t_plexgrid), allocatable :: plex
    integer(mpiint) :: subcomm
    integer(iintegers) :: num_subcomm_masters
    type(tVecScatter) :: ctx_hhl
    type(tVecScatter) :: ctx_albedo
    type(tVecScatter) :: ctx_optprop
    type(tVecScatter), allocatable :: ctx_planck
    type(tVecScatter), allocatable :: ctx_edir
    type(tVecScatter) :: ctx_ediff
    type(tVecScatter) :: ctx_abso
    type(tVec) :: albedo, kabs, ksca, g, planck
    type(tVec), allocatable :: planck_srfc
    type(t_rayli_info_buildings), allocatable :: buildings_info
  end type

  type(t_rayli_info), allocatable :: rayli_info

contains

  !> @brief wrapper for the rayli montecarlo solver
  !> @details solve the radiative transfer equation with the wedge bindings to  rayli
  ! tasks:
  ! * copy pprts dm to Zero (and all shared mem root-ranks)
  ! * implement pprts to wedge interface
  ! * average results over shared mem comm
  ! * distribute results
  subroutine pprts_rayli_wrapper(lcall_solver, lcall_snap, solver, edirTOA, solution, ierr, opt_buildings)
    logical, intent(in) :: lcall_solver, lcall_snap
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), intent(inout) :: solution
    integer(mpiint), intent(out) :: ierr
    type(t_pprts_buildings), intent(in), optional :: opt_buildings

    integer(mpiint) :: myid, numnodes
    integer(mpiint) :: submyid, subnumnodes

    real(ireals) :: sundir(3)

    ierr = 0

    if (all([lcall_solver, lcall_snap] .eqv. .false.)) return

    sundir = spherical_2_cartesian(solver%sun%phi, solver%sun%theta) &
            & * edirTOA

    call init_pprts_rayli_wrapper(solver, solution, rayli_info, opt_buildings=opt_buildings)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)
    call mpi_comm_rank(rayli_info%subcomm, submyid, ierr); call CHKERR(ierr)
    call mpi_comm_size(rayli_info%subcomm, subnumnodes, ierr); call CHKERR(ierr)

    call prepare_input()
    call call_solver(rayli_info%plex_solution)
    call transfer_result(rayli_info%plex_solution, solution)

  contains
    subroutine prepare_input()
      type(tVec) :: glob_albedo, glob_kabs, glob_ksca, glob_g, glob_B, glob_Bsrfc
      character(len=*), parameter :: log_event_name = "pprts_rayli_prepare_input"

      PetscClassId :: cid
      PetscLogEvent :: log_event

      call PetscClassIdRegister("pprts_rayli", cid, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(log_event_name, cid, log_event, ierr); call CHKERR(ierr)
      call PetscLogEventBegin(log_event, ierr); call CHKERR(ierr)

      associate ( &
          & atm => solver%atm,            &
          & Cs => solver%Csrfc_one,      &
          & Ca1 => solver%C_one_atm1,     &
          & Ca => solver%C_one_atm,      &
          & Cv => solver%Cvert_one_atm1, &
          & ri => rayli_info)

        call DMGetGlobalVector(Cs%da, glob_albedo, ierr); call CHKERR(ierr)
        call f90VecToPetsc(atm%albedo, Cs%da, glob_albedo)
        call VecScatterBegin(ri%ctx_albedo, glob_albedo, ri%albedo, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call VecScatterEnd(ri%ctx_albedo, glob_albedo, ri%albedo, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(Cs%da, glob_albedo, ierr); call CHKERR(ierr)

        if (allocated(atm%Bsrfc)) then
          if (.not. allocated(ri%planck_srfc)) call CHKERR(1_mpiint, 'ri%planck_srfc not allocated. Developer error!')
          call DMGetGlobalVector(Cs%da, glob_Bsrfc, ierr); call CHKERR(ierr)
          call f90VecToPetsc(atm%Bsrfc, Cs%da, glob_Bsrfc)
          call VecScatterBegin(ri%ctx_albedo, glob_Bsrfc, ri%planck_srfc, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call VecScatterEnd(ri%ctx_albedo, glob_Bsrfc, ri%planck_srfc, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call DMRestoreGlobalVector(Cs%da, glob_Bsrfc, ierr); call CHKERR(ierr)
        end if

        if (present(opt_buildings)) then
          ! Create a pprts buildings object on each subcomm
          call local_buildings_to_subcomm_buildings(&
            & opt_buildings, &
            & ri%buildings_info, &
            & ri%buildings_info%subcomm_buildings, &
            & ierr); call CHKERR(ierr)

          ! Then build a plex_buildings object out of that one
          if (submyid .eq. 0) then
            call pprts_buildings_to_plex(solver, &
              & ri%plex, &
              & ri%buildings_info%subcomm_buildings, &
              & ri%buildings_info%plex_buildings, ierr); call CHKERR(ierr)
          end if
        else
          if (allocated(ri%buildings_info)) then
            if (allocated(ri%buildings_info%plex_buildings)) deallocate (ri%buildings_info%plex_buildings)
            if (allocated(ri%buildings_info%subcomm_buildings)) deallocate (ri%buildings_info%subcomm_buildings)
          end if
        end if

        call DMGetGlobalVector(Ca%da, glob_kabs, ierr); call CHKERR(ierr)
        call f90VecToPetsc(atm%kabs, Ca%da, glob_kabs)
        call VecScatterBegin(ri%ctx_optprop, glob_kabs, ri%kabs, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call VecScatterEnd(ri%ctx_optprop, glob_kabs, ri%kabs, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(Ca%da, glob_kabs, ierr); call CHKERR(ierr)

        call DMGetGlobalVector(Ca%da, glob_ksca, ierr); call CHKERR(ierr)
        call f90VecToPetsc(atm%ksca, Ca%da, glob_ksca)
        call VecScatterBegin(ri%ctx_optprop, glob_ksca, ri%ksca, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call VecScatterEnd(ri%ctx_optprop, glob_ksca, ri%ksca, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(Ca%da, glob_ksca, ierr); call CHKERR(ierr)

        call DMGetGlobalVector(Ca%da, glob_g, ierr); call CHKERR(ierr)
        call f90VecToPetsc(atm%g, Ca%da, glob_g)
        call VecScatterBegin(ri%ctx_optprop, glob_g, ri%g, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call VecScatterEnd(ri%ctx_optprop, glob_g, ri%g, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(Ca%da, glob_g, ierr); call CHKERR(ierr)

        if (solution%lthermal_rad) then
          call DMGetGlobalVector(Ca1%da, glob_B, ierr); call CHKERR(ierr)
          call f90VecToPetsc(atm%planck, Ca1%da, glob_B)
          call VecScatterBegin(ri%ctx_planck, glob_B, ri%planck, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call VecScatterEnd(ri%ctx_planck, glob_B, ri%planck, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call DMRestoreGlobalVector(Ca1%da, glob_B, ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(ri%planck, PETSC_NULL_VEC, '-show_rayli_planck', ierr); call CHKERR(ierr)
        end if

        call PetscObjectViewFromOptions(ri%kabs, PETSC_NULL_VEC, '-show_rayli_kabs', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(ri%ksca, PETSC_NULL_VEC, '-show_rayli_ksca', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(ri%g, PETSC_NULL_VEC, '-show_rayli_g', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(ri%albedo, PETSC_NULL_VEC, '-show_rayli_albedo', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(ri%planck_srfc, PETSC_NULL_VEC, '-show_rayli_planck_srfc', ierr); call CHKERR(ierr)
      end associate
      call PetscLogEventEnd(log_event, ierr); call CHKERR(ierr)
    end subroutine

    subroutine local_buildings_to_subcomm_buildings(localB, B_info, subB, ierr)
      type(t_pprts_buildings), intent(in) :: localB
      type(t_rayli_info_buildings), intent(in) :: B_info
      type(t_pprts_buildings), allocatable, intent(inout) :: subB
      type(tVec) :: vlocal, vsub
      real(ireals), allocatable :: glob_idx(:), sub_idx(:)
      integer(iintegers) :: m, global_da_offsets(4), idx(4), nlocal
      integer(mpiint), intent(out) :: ierr
      ierr = 0

      associate (atm => solver%atm, C1 => solver%C_one)
        if (.not. allocated(subB)) then
          call clone_buildings(localB, subB, .false., ierr); call CHKERR(ierr)
          call ndarray_offsets([6_iintegers, C1%glob_zm, C1%glob_xm, C1%glob_ym], subB%da_offsets)
          allocate (subB%iface(B_info%Nglob))
          allocate (subB%albedo(B_info%Nglob))
          allocate (subB%planck(B_info%Nglob))
        end if

        nlocal = size(localB%iface)

        ! Setup Albedo
        call VecCreateMPIWithArray(solver%comm, i1, &
          & nlocal, PETSC_DECIDE, localB%albedo, &
          & vlocal, ierr); call CHKERR(ierr)

        call VecCreateSeqWithArray(PETSC_COMM_SELF, i1, B_info%Nglob, subB%albedo, vsub, ierr); call CHKERR(ierr)
        call VecScatterBegin(B_info%ctx_albedo, vlocal, vsub, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call VecScatterEnd(B_info%ctx_albedo, vlocal, vsub, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

        call PetscObjectViewFromOptions(vsub, PETSC_NULL_VEC, '-show_rayli_buildings_albedo', ierr); call CHKERR(ierr)
        call VecDestroy(vlocal, ierr); call CHKERR(ierr)
        call VecDestroy(vsub, ierr); call CHKERR(ierr)

        if (allocated(localB%planck)) then ! Setup Planck
          call VecCreateMPIWithArray(solver%comm, i1, &
            & nlocal, PETSC_DECIDE, localB%planck, &
            & vlocal, ierr); call CHKERR(ierr)

          call VecCreateSeqWithArray(PETSC_COMM_SELF, i1, B_info%Nglob, subB%planck, vsub, ierr); call CHKERR(ierr)
          call VecScatterBegin(B_info%ctx_albedo, vlocal, vsub, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call VecScatterEnd(B_info%ctx_albedo, vlocal, vsub, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

          call PetscObjectViewFromOptions(vsub, PETSC_NULL_VEC, '-show_rayli_buildings_planck', ierr); call CHKERR(ierr)
          call VecDestroy(vlocal, ierr); call CHKERR(ierr)
          call VecDestroy(vsub, ierr); call CHKERR(ierr)
        end if

        ! Setup buildings faces
        ! i.e. find out what the face idx of buildings would be with global offsets
        call ndarray_offsets( &
          & [6_iintegers, C1%glob_zm, C1%glob_xm, C1%glob_ym], &
          & global_da_offsets)

        allocate (glob_idx(nlocal), sub_idx(B_info%Nglob))
        do m = 1, nlocal
          call ind_1d_to_nd(localB%da_offsets, localB%iface(m), idx)

          associate (d => idx(1), k => idx(2), i => idx(3), j => idx(4))
            glob_idx(m) = real(ind_nd_to_1d(global_da_offsets, [d, C1%zs + k, C1%xs + i, C1%ys + j]), ireals)
          end associate

        end do
        call VecCreateMPIWithArray(solver%comm, i1, &
          & nlocal, PETSC_DECIDE, glob_idx, &
          & vlocal, ierr); call CHKERR(ierr)

        call VecCreateSeqWithArray(PETSC_COMM_SELF, i1, B_info%Nglob, sub_idx, vsub, ierr); call CHKERR(ierr)
        call VecScatterBegin(B_info%ctx_albedo, vlocal, vsub, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call VecScatterEnd(B_info%ctx_albedo, vlocal, vsub, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        subB%iface(:) = int(sub_idx(:), kind=iintegers)

        call PetscObjectViewFromOptions(vsub, PETSC_NULL_VEC, '-show_rayli_buildings_iface', ierr); call CHKERR(ierr)
        call VecDestroy(vlocal, ierr); call CHKERR(ierr)
        call VecDestroy(vsub, ierr); call CHKERR(ierr)

        call check_buildings_consistency(subB, &
          & C1%glob_zm, C1%glob_xm, C1%glob_ym, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    subroutine call_solver(plex_solution)
      type(t_state_container), intent(inout) :: plex_solution
      logical :: gotmsg
      integer(iintegers) :: idummy
      integer(mpiint) :: isub, status(MPI_STATUS_SIZE), imp_request
      integer(mpiint), parameter :: FINALIZEMSG = 1
      integer(mpiint) :: run_rank = 0

      integer(iintegers) :: Nphotons
      real(ireals) :: Nphotons_r
      logical :: lflg

      character(len=*), parameter :: log_event_name = "pprts_rayli_call_solver"
      PetscClassId :: cid
      PetscLogEvent :: log_event

      call PetscClassIdRegister("pprts_rayli", cid, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(log_event_name, cid, log_event, ierr); call CHKERR(ierr)
      call PetscLogEventBegin(log_event, ierr); call CHKERR(ierr)

      associate (ri => rayli_info)
        call VecGetSize(ri%albedo, Nphotons, ierr); call CHKERR(ierr)
        nphotons_r = real(Nphotons * 10, ireals)
        call get_petsc_opt(solver%prefix, &
                           "-pprts_rayli_photons", nphotons_r, lflg, ierr); call CHKERR(ierr)

        Nphotons_r = Nphotons_r / real(ri%num_subcomm_masters, ireals)

        if (submyid .eq. run_rank) then
          if (plex_solution%lthermal_rad) then
            if (present(opt_buildings)) then
              call rayli_wrapper(lcall_solver, lcall_snap, &
                & ri%plex, ri%kabs, ri%ksca, ri%g, ri%albedo, &
                & plex_solution, plck=ri%planck, plck_srfc=ri%planck_srfc, &
                & nr_photons=Nphotons_r, petsc_log=solver%logs%rayli_tracing, &
                & opt_buildings=ri%buildings_info%plex_buildings, &
                & opt_Nthreads=int(subnumnodes, iintegers))
            else
              call rayli_wrapper(&
                & lcall_solver, lcall_snap, &
                & ri%plex, ri%kabs, ri%ksca, ri%g, ri%albedo, &
                & plex_solution, plck=ri%planck, plck_srfc=ri%planck_srfc, &
                & nr_photons=Nphotons_r, petsc_log=solver%logs%rayli_tracing)
            end if

          else
            if (present(opt_buildings)) then
              call rayli_wrapper(lcall_solver, lcall_snap, &
                & ri%plex, ri%kabs, ri%ksca, ri%g, ri%albedo, &
                & plex_solution, sundir=sundir, &
                & nr_photons=Nphotons_r, petsc_log=solver%logs%rayli_tracing, &
                & opt_buildings=ri%buildings_info%plex_buildings, &
                & opt_Nthreads=int(subnumnodes, iintegers))
            else
              call rayli_wrapper(lcall_solver, lcall_snap, &
                & ri%plex, ri%kabs, ri%ksca, ri%g, ri%albedo, &
                & plex_solution, sundir=sundir, &
                & nr_photons=Nphotons_r, petsc_log=solver%logs%rayli_tracing)
            end if

            call PetscObjectViewFromOptions(plex_solution%edir, PETSC_NULL_VEC, '-show_plex_rayli_edir', ierr); call CHKERR(ierr)
          end if

          call PetscObjectViewFromOptions(plex_solution%ediff, PETSC_NULL_VEC, '-show_plex_rayli_ediff', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(plex_solution%abso, PETSC_NULL_VEC, '-show_plex_rayli_abso', ierr); call CHKERR(ierr)

          do isub = 0, subnumnodes - 1 ! send finalize msg to all others to stop waiting
            if (isub .ne. run_rank) then
              call mpi_isend(-i1, 1_mpiint, imp_iinteger, isub, FINALIZEMSG, &
                & ri%subcomm, imp_request, ierr); call CHKERR(ierr)
            end if
          end do
        else
          lazy_wait: do ! prevent eager MPI polling while the one rank performs rayli computations
            call mpi_iprobe(MPI_ANY_SOURCE, FINALIZEMSG, ri%subcomm, gotmsg, status, ierr); call CHKERR(ierr)
            if (gotmsg) then
              call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, &
                & ri%subcomm, status, ierr); call CHKERR(ierr)
              exit lazy_wait
            end if
            call sleep(1)
          end do lazy_wait
        end if
      end associate
      call PetscLogEventEnd(log_event, ierr); call CHKERR(ierr)
    end subroutine

    subroutine transfer_result(plex_solution, solution)
      type(t_state_container), intent(in) :: plex_solution
      type(t_state_container), intent(inout) :: solution

      real(ireals), pointer :: x(:, :, :, :) => null(), x1d(:) => null()
      real(ireals) :: fac

      type(tVec) :: abso_in_W
      type(tPetscSection) :: abso_section, geomSection
      real(ireals), pointer :: xabso(:), geoms(:)
      integer(iintegers) :: icell, cStart, cEnd, voff, geom_offset

      character(len=*), parameter :: log_event_name = "pprts_rayli_transfer_result"
      PetscClassId :: cid
      PetscLogEvent :: log_event

      call PetscClassIdRegister("pprts_rayli", cid, ierr); call CHKERR(ierr)
      call PetscLogEventRegister(log_event_name, cid, log_event, ierr); call CHKERR(ierr)
      call PetscLogEventBegin(log_event, ierr); call CHKERR(ierr)

      fac = one / real(rayli_info%num_subcomm_masters, ireals)

      if (allocated(solution%edir)) then
        call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
      end if
      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
      call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)

      if (allocated(solution%edir)) then
        call VecScatterBegin(rayli_info%ctx_edir, plex_solution%edir, solution%edir, &
                             ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
      end if
      call VecScatterBegin(rayli_info%ctx_ediff, plex_solution%ediff, solution%ediff, &
                           ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
      if (allocated(solution%edir)) then
        call VecScatterEnd(rayli_info%ctx_edir, plex_solution%edir, solution%edir, &
                           ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
      end if
      call VecScatterEnd(rayli_info%ctx_ediff, plex_solution%ediff, solution%ediff, &
                         ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

      ! need absorption as W, not W/m3 to handle collapsed area correctly
      associate (plex => rayli_info%plex)
        call VecDuplicate(plex_solution%abso, abso_in_W, ierr); call CHKERR(ierr)
        call VecCopy(plex_solution%abso, abso_in_W, ierr); call CHKERR(ierr)
        call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
        call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)
        call VecGetArrayF90(abso_in_W, xabso, ierr); call CHKERR(ierr)
        call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
        do icell = cStart, cEnd - 1
          call PetscSectionGetFieldOffset(geomSection, icell, i3, geom_offset, ierr); call CHKERR(ierr) ! cell dz
          call PetscSectionGetOffset(abso_section, icell, voff, ierr); call CHKERR(ierr)
          xabso(i1 + voff) = xabso(i1 + voff) * geoms(i1 + geom_offset)
        end do
        call VecRestoreArrayF90(abso_in_W, xabso, ierr); call CHKERR(ierr)
        call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
      end associate

      call VecScatterBegin(rayli_info%ctx_abso, abso_in_W, solution%abso, &
                           ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

      call VecScatterEnd(rayli_info%ctx_abso, abso_in_W, solution%abso, &
                         ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

      call VecDestroy(abso_in_W, ierr); call CHKERR(ierr)

      if (allocated(solution%edir)) then
        call getVecPointer(solver%C_dir%da, solution%edir, x1d, x)
        x(i0, :, :, :) = x(i0, :, :, :) * fac / 2._ireals
        x(i1:i2, :, :, :) = x(i1:i2, :, :, :) * fac
        call restoreVecPointer(solver%C_dir%da, solution%edir, x1d, x)
      end if
      call VecScale(solution%ediff, fac / 2._ireals, ierr); call CHKERR(ierr)
      call VecScale(solution%abso, fac / 2._ireals, ierr); call CHKERR(ierr)

      call getVecPointer(solver%C_one%da, solution%abso, x1d, x)
      x(i0, solver%C_one%zs + 1:, :, :) = x(i0, solver%C_one%zs + 1:, :, :) &
        & / solver%atm%dz(atmk(solver%atm, solver%C_one_atm%zs + 1):solver%C_one_atm%ze, :, :)

      ! the top value may have been collapsed, if this is the case take weighted mean with dz
      x(i0, solver%C_one%zs, :, :) = x(i0, solver%C_one%zs, :, :) &
        & / sum(solver%atm%dz(:atmk(solver%atm, solver%C_one_atm%zs), :, :), dim=1)
      call restoreVecPointer(solver%C_one%da, solution%abso, x1d, x)

      if (allocated(solution%edir)) then
        call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, '-show_rayli_edir', ierr); call CHKERR(ierr)
      end if
      call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, '-show_rayli_ediff', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-show_rayli_abso', ierr); call CHKERR(ierr)

      !Rayli solver returns fluxes as [W]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is up to date (to prevent absoprtion computations)
      solution%lchanged = .false.
      call PetscLogEventEnd(log_event, ierr); call CHKERR(ierr)
    end subroutine
  end subroutine

  subroutine destroy_rayli_info()
    integer(mpiint) :: ierr
    if (.not. allocated(rayli_info)) return

    associate (ri => rayli_info)
      call destroy_solution(ri%plex_solution)
      call VecScatterDestroy(ri%ctx_hhl, ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_albedo, ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_optprop, ierr); call CHKERR(ierr)
      if (allocated(ri%ctx_planck)) then
        call VecScatterDestroy(ri%ctx_planck, ierr); call CHKERR(ierr)
        deallocate (ri%ctx_planck)
      end if
      if (allocated(ri%ctx_edir)) then
        call VecScatterDestroy(ri%ctx_edir, ierr); call CHKERR(ierr)
        deallocate (ri%ctx_edir)
      end if
      call VecScatterDestroy(ri%ctx_ediff, ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_abso, ierr); call CHKERR(ierr)
      call VecDestroy(ri%albedo, ierr); call CHKERR(ierr)
      call VecDestroy(ri%kabs, ierr); call CHKERR(ierr)
      call VecDestroy(ri%ksca, ierr); call CHKERR(ierr)
      call VecDestroy(ri%g, ierr); call CHKERR(ierr)
      call VecDestroy(ri%planck, ierr); call CHKERR(ierr)
      if (allocated(ri%planck_srfc)) then
        call VecDestroy(ri%planck_srfc, ierr); call CHKERR(ierr)
        deallocate (ri%planck_srfc)
      end if
      if (allocated(ri%buildings_info)) then
        if (allocated(ri%buildings_info%subcomm_buildings)) then
          call destroy_buildings(ri%buildings_info%subcomm_buildings, ierr); call CHKERR(ierr)
        end if
        if (allocated(ri%buildings_info%plex_buildings)) then
          call destroy_buildings(ri%buildings_info%plex_buildings, ierr); call CHKERR(ierr)
        end if
      end if
      if (allocated(ri%plex)) then
        call destroy_plexgrid(ri%plex)
      end if
      call mpi_comm_free(ri%subcomm, ierr); call CHKERR(ierr)
    end associate

    deallocate (rayli_info)
  end subroutine

  !> @brief setup mpi scatters to move pprts domains into plexrt meshes on each local subcomm root
  subroutine init_pprts_rayli_wrapper(solver, solution, rayli_info, opt_buildings)
    class(t_solver), intent(in) :: solver
    type(t_state_container), intent(in) :: solution
    type(t_rayli_info), allocatable, intent(inout) :: rayli_info
    type(t_pprts_buildings), intent(in), optional :: opt_buildings
    integer(iintegers) :: i_am_master
    integer(mpiint) :: submyid, ierr

    if (.not. allocated(rayli_info)) then
      allocate (rayli_info)

      ! Scatter height info to subcomm[0]
      call gen_shared_subcomm(solver%comm, rayli_info%subcomm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(rayli_info%subcomm, submyid, ierr); call CHKERR(ierr)

      if (submyid .eq. 0) then
        i_am_master = 1
      else
        i_am_master = 0
      end if
      call imp_allreduce_sum(solver%comm, i_am_master, rayli_info%num_subcomm_masters)

      call setup_mesh()

      if (submyid .eq. 0) then
        call prepare_solution(&
          & rayli_info%plex%edir_dm, &
          & rayli_info%plex%ediff_dm, &
          & rayli_info%plex%abso_dm, &
          & lsolar=.true., &
          & lthermal=.true., &
          & solution=rayli_info%plex_solution)
      end if

      ! Setup albedo scatter context
      call setup_srfc_opp_sct_ctx()

      ! Setup result scatter contexts
      call setup_ediff_scatter_context()
      call setup_abso_scatter_context()
    end if

    if (solution%lsolar_rad) call setup_edir_scatter_context()
    if (solution%lthermal_rad) call setup_planck_scatter_context()

    call setup_buildings_scatter_context(opt_buildings)

    rayli_info%plex_solution%lsolar_rad = solution%lsolar_rad
    rayli_info%plex_solution%lthermal_rad = solution%lthermal_rad
    rayli_info%plex_solution%uid = solution%uid

  contains

    subroutine setup_mesh()
      integer(iintegers) :: Nhhl
      type(tVec) :: vertex_hhl, hhl
      real(ireals), pointer :: xhhl(:, :, :, :) => null(), xhhl1d(:) => null()
      type(tPetscSection) :: coord_section
      type(tVec) :: coordinates
      real(ireals), pointer :: coords(:)
      integer(iintegers) :: i, vStart, vEnd, ivert, coord_offset

      AO :: dmda_ao
      type(tIS) :: is_in, is_out

      type(tDM) :: dm2d, dm2d_dist, dm3d
      integer(iintegers), allocatable :: zindex(:)

      call DMGetGlobalVector(solver%Cvert_one_atm1%da, vertex_hhl, ierr); call CHKERR(ierr)

      call interpolate_cell_values_to_vertices(&
        & solver%C_one_atm1_box, solver%atm%hhl, &
        & solver%Cvert_one_atm1, vertex_hhl)

      call PetscObjectViewFromOptions(vertex_hhl, PETSC_NULL_VEC, '-show_rayli_hhl', ierr); call CHKERR(ierr)

      if (submyid .eq. 0) then
        call VecGetSize(vertex_hhl, Nhhl, ierr); call CHKERR(ierr)
      else
        Nhhl = 0
      end if

      call VecCreateSeq(PETSC_COMM_SELF, Nhhl, hhl, ierr); call CHKERR(ierr)

      call ISCreateStride(PETSC_COMM_SELF, Nhhl, 0_iintegers, 1_iintegers, is_in, ierr); call CHKERR(ierr)
      call ISCreateStride(PETSC_COMM_SELF, Nhhl, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)

      call DMDAGetAO(solver%Cvert_one_atm1%da, dmda_ao, ierr); call CHKERR(ierr)
      call AOApplicationToPetscIS(dmda_ao, is_in, ierr); call CHKERR(ierr)

      call gen_shared_scatter_ctx(vertex_hhl, hhl, rayli_info%ctx_hhl, ierr, is_in, is_out); call CHKERR(ierr)
      call VecScatterBegin(rayli_info%ctx_hhl, vertex_hhl, hhl, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
      call VecScatterEnd(rayli_info%ctx_hhl, vertex_hhl, hhl, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(solver%Cvert_one_atm1%da, vertex_hhl, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(hhl, PETSC_NULL_VEC, '-show_rayli_hhl', ierr); call CHKERR(ierr)

      ! Setup Plex
      associate ( &
          & atm => solver%atm, &
          & Ca1 => solver%C_one_atm1, &
          & Ca => solver%C_one_atm, &
          & Cv => solver%Cvert_one_atm1)

        if (submyid .eq. 0) then
          call create_2d_regular_plex(PETSC_COMM_SELF, Ca%glob_xm + 1, Ca%glob_ym + 1, &
            & dm2d, dm2d_dist, opt_dx=atm%dx, opt_dy=atm%dy, lverbose=.true.)

          call VecGetArrayReadF90(hhl, xhhl1d, ierr); call CHKERR(ierr)
          xhhl(i1:i1, i1:Cv%glob_zm, i1:Cv%glob_xm, i1:Cv%glob_ym) => xhhl1d

          call dmplex_2D_to_3D(dm2d, Ca1%glob_zm, xhhl(i1, :, i1, i1), [zero, zero, -huge(zero) / 100], dm3d, zindex)

          !set height lvls on vertices
          call DMGetCoordinateSection(dm3d, coord_section, ierr); call CHKERR(ierr)
          call DMGetCoordinatesLocal(dm3d, coordinates, ierr); call CHKERR(ierr)
          call DMPlexGetDepthStratum(dm3d, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices
          call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

          do i = 1, size(xhhl1d)
            ivert = vStart + i - i1
            call PetscSectionGetOffset(coord_section, ivert, coord_offset, ierr); call CHKERR(ierr)
            coords(i1 + coord_offset + i2) = xhhl1d(i)
          end do
          call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

          call setup_plexgrid(dm2d, dm3d, Ca%glob_zm, zindex, rayli_info%plex)

          call PetscObjectViewFromOptions(dm3d, PETSC_NULL_DM, '-show_rayli_dm3d', ierr); call CHKERR(ierr)
          nullify (xhhl)
          call VecRestoreArrayReadF90(hhl, xhhl1d, ierr); call CHKERR(ierr)

          call setup_edir_dmplex(rayli_info%plex, rayli_info%plex%dm, i1, i0, i1, rayli_info%plex%horizface1_dm)
          call setup_edir_dmplex(rayli_info%plex, rayli_info%plex%dm, solver%dirtop%dof, i1, i1, rayli_info%plex%edir_dm)
          call setup_ediff_dmplex(rayli_info%plex, rayli_info%plex%dm, i1, i1, i2, rayli_info%plex%ediff_dm)
          call setup_abso_dmplex(rayli_info%plex%dm, rayli_info%plex%abso_dm)

        else ! subcomm slave

          allocate (rayli_info%plex_solution%edir)
          call VecCreateSeq(PETSC_COMM_SELF, i0, rayli_info%plex_solution%edir, ierr); call CHKERR(ierr)
          allocate (rayli_info%plex_solution%ediff)
          call VecCreateSeq(PETSC_COMM_SELF, i0, rayli_info%plex_solution%ediff, ierr); call CHKERR(ierr)
          allocate (rayli_info%plex_solution%abso)
          call VecCreateSeq(PETSC_COMM_SELF, i0, rayli_info%plex_solution%abso, ierr); call CHKERR(ierr)

        end if ! is subcomm master
      end associate
      call VecDestroy(hhl, ierr); call CHKERR(ierr)
    end subroutine

    subroutine setup_buildings_scatter_context(opt_buildings)
      type(t_pprts_buildings), intent(in), optional :: opt_buildings
      type(tVec) :: loc_albedo ! a petsc vec version of buildings%albedo
      type(tVec) :: albedo ! tmp vec to generate scatter context
      integer(iintegers) :: nlocal, iStart
      type(tIS) :: is_in
      type(tIS) :: is_out
      integer(mpiint) :: numnodes, ierr

      if (.not. present(opt_buildings)) return

      if (allocated(rayli_info%buildings_info)) return
      allocate (rayli_info%buildings_info)

      call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)
      associate (B => opt_buildings, BI => rayli_info%buildings_info)
        nlocal = size(B%iface)
        call VecCreateMPIWithArray(solver%comm, i1, nlocal, PETSC_DECIDE, B%albedo, &
          & loc_albedo, ierr); call CHKERR(ierr)
        call VecGetOwnershipRange(loc_albedo, iStart, PETSC_NULL_INTEGER, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(loc_albedo, PETSC_NULL_VEC, &
          & '-show_rayli_buildings_loc_vec', ierr); call CHKERR(ierr)

        if (submyid .eq. 0) then
          call VecGetSize(loc_albedo, BI%Nglob, ierr); call CHKERR(ierr)
        else
          BI%Nglob = 0
        end if

        call VecCreateSeq(PETSC_COMM_SELF, BI%Nglob, albedo, ierr); call CHKERR(ierr)
        call ISCreateStride(PETSC_COMM_SELF, BI%Nglob, 0_iintegers, 1_iintegers, is_in, ierr); call CHKERR(ierr)
        call ISCreateStride(PETSC_COMM_SELF, BI%Nglob, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)
        call gen_shared_scatter_ctx(loc_albedo, albedo, BI%ctx_albedo, ierr, &
          & is_in, is_out); call CHKERR(ierr)
        call ISDestroy(is_out, ierr); call CHKERR(ierr)
        call ISDestroy(is_in, ierr); call CHKERR(ierr)

        call VecDestroy(loc_albedo, ierr); call CHKERR(ierr)
        call VecDestroy(albedo, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    subroutine setup_srfc_opp_sct_ctx()
      integer(iintegers) :: Nalbedo, Noptprop, i, k
      AO :: dmda_ao
      type(tIS) :: is_in, is_out
      integer(iintegers), allocatable :: is_data(:)
      type(tVec) :: vec_albedo, vec_optprop

      call DMGetGlobalVector(solver%Csrfc_one%da, vec_albedo, ierr); call CHKERR(ierr)
      if (submyid .eq. 0) then
        call VecGetSize(vec_albedo, Nalbedo, ierr); call CHKERR(ierr)
      else
        Nalbedo = i0
      end if
      call VecCreateSeq(PETSC_COMM_SELF, Nalbedo * i2, rayli_info%albedo, ierr); call CHKERR(ierr)
      if (allocated(solver%atm%Bsrfc)) then
        if (.not. allocated(rayli_info%planck_srfc)) allocate (rayli_info%planck_srfc)
        call VecCreateSeq(PETSC_COMM_SELF, Nalbedo * i2, rayli_info%planck_srfc, ierr); call CHKERR(ierr)
      end if
      allocate (is_data(Nalbedo * i2), source=(/(i / i2, i=i0, Nalbedo * i2 - i1)/))
      call ISCreateGeneral(PETSC_COMM_SELF, Nalbedo * i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
      call ISCreateStride(PETSC_COMM_SELF, Nalbedo * i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)

      call DMDAGetAO(solver%Csrfc_one%da, dmda_ao, ierr); call CHKERR(ierr)
      call AOApplicationToPetscIS(dmda_ao, is_in, ierr); call CHKERR(ierr)

      call gen_shared_scatter_ctx(vec_albedo, rayli_info%albedo, rayli_info%ctx_albedo, ierr, &
        & is_in, is_out); call CHKERR(ierr)
      call ISDestroy(is_out, ierr); call CHKERR(ierr)
      call ISDestroy(is_in, ierr); call CHKERR(ierr)
      deallocate (is_data)
      call DMRestoreGlobalVector(solver%Csrfc_one%da, vec_albedo, ierr); call CHKERR(ierr)

      ! Setup optprop_scatter_ctx
      call DMGetGlobalVector(solver%C_one_atm%da, vec_optprop, ierr); call CHKERR(ierr)
      if (submyid .eq. 0) then
        call VecGetSize(vec_optprop, Noptprop, ierr); call CHKERR(ierr)
      else
        Noptprop = 0
      end if
      call VecCreateSeq(PETSC_COMM_SELF, Noptprop * i2, rayli_info%kabs, ierr); call CHKERR(ierr)
      call VecCreateSeq(PETSC_COMM_SELF, Noptprop * i2, rayli_info%ksca, ierr); call CHKERR(ierr)
      call VecCreateSeq(PETSC_COMM_SELF, Noptprop * i2, rayli_info%g, ierr); call CHKERR(ierr)
      allocate (is_data(i0:Noptprop * i2 - i1))
      if (submyid .eq. 0) then
        do i = i0, solver%C_one_atm%glob_xm * solver%C_one_atm%glob_ym - i1
          do k = i0, solver%C_one_atm%glob_zm - i1
            is_data((i * i2) * solver%C_one_atm%glob_zm + k) = i * solver%C_one_atm%glob_zm + k
            is_data((i * i2 + 1) * solver%C_one_atm%glob_zm + k) = i * solver%C_one_atm%glob_zm + k
          end do
        end do
      end if
      call ISCreateGeneral(PETSC_COMM_SELF, Noptprop * i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
      call ISCreateStride(PETSC_COMM_SELF, Noptprop * i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)

      call DMDAGetAO(solver%C_one_atm%da, dmda_ao, ierr); call CHKERR(ierr)
      call AOApplicationToPetscIS(dmda_ao, is_in, ierr); call CHKERR(ierr)

      call gen_shared_scatter_ctx(vec_optprop, rayli_info%kabs, rayli_info%ctx_optprop, ierr, &
        & is_in, is_out); call CHKERR(ierr)
      call ISDestroy(is_out, ierr); call CHKERR(ierr)
      call ISDestroy(is_in, ierr); call CHKERR(ierr)
      deallocate (is_data)
      call DMRestoreGlobalVector(solver%C_one_atm%da, vec_optprop, ierr); call CHKERR(ierr)
    end subroutine

    subroutine setup_edir_scatter_context()
      integer(iintegers) :: i1d, k, i, j, l, m, iface, numDof, voff, ak
      integer(iintegers) :: pprts_offsets(4), plex_cells(2)

      AO :: dmda_ao
      type(tPetscSection) :: edirsection
      type(tIS) :: is_in, is_out
      integer(iintegers) :: is_data_size
      integer(iintegers), allocatable :: is_data_in(:)
      integer(iintegers), allocatable :: is_data_out(:)

      associate (ri => rayli_info, Cdir => solver%C_dir, C1 => solver%C_one, Ca => solver%C_one_atm)

        if (allocated(ri%ctx_edir)) return
        allocate (ri%ctx_edir)

        if (submyid .eq. 0) then
          is_data_size = Cdir%glob_zm * Cdir%glob_xm * Cdir%glob_ym * i2 & ! all horizontal fluxes twice
            & + C1%glob_zm * C1%glob_xm * C1%glob_ym * i2 ! plus fluxes on vertical faces once for x and y
        else
          is_data_size = 0
        end if
        allocate (is_data_in(0:is_data_size - 1))
        allocate (is_data_out(0:is_data_size - 1))
        is_data_in(:) = -1
        is_data_out(:) = -1

        if (submyid .eq. 0) then
          ! First do it for all the atmosphere dofs
          call ndarray_offsets(      &
            &[Cdir%dof,      &
            & Cdir%glob_zm,  &
            & Cdir%glob_xm,  &
            & Cdir%glob_ym], &
            & pprts_offsets)

          call DMGetSection(ri%plex%edir_dm, edirsection, ierr); call CHKERR(ierr)

          l = 0
          do j = 0, Cdir%glob_ym - 1
            do i = 0, Cdir%glob_xm - 1
              do k = 0, Cdir%glob_zm - 2
                if (k .eq. 0) then
                  ak = 0
                else
                  ak = atmk(solver%atm, k)
                end if
                call pprts_cell_to_plex_cell_idx(Ca, [ak, i, j], ri%plex, plex_cells, ierr); call CHKERR(ierr)

                ! on top faces:
                i1d = ind_nd_to_1d(pprts_offsets, [i0, k, i, j], cstyle=.true.)
                do m = 1, 2
                  iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_TOP_FACE)
                  if (iface .ge. 0) then
                    call PetscSectionGetOffset(edirsection, iface, voff, ierr); call CHKERR(ierr)
                    is_data_in(l) = i1d
                    is_data_out(l) = voff
                    l = l + 1
                  end if
                end do

                ! on x-side faces:
                i1d = ind_nd_to_1d(pprts_offsets, [i1, k, i, j], cstyle=.true.)
                do m = 1, 2
                  iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_LEFT_FACE)
                  if (iface .ge. 0) then
                    call PetscSectionGetDof(edirSection, iface, numDof, ierr); call CHKERR(ierr)
                    if (numDof .gt. 0) then
                      call PetscSectionGetOffset(edirsection, iface, voff, ierr); call CHKERR(ierr)
                      is_data_in(l) = i1d
                      is_data_out(l) = voff
                      l = l + 1
                    end if
                  end if
                end do

                ! on y-side faces:
                i1d = ind_nd_to_1d(pprts_offsets, [i2, k, i, j], cstyle=.true.)
                do m = 1, 2
                  iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_REAR_FACE)
                  if (iface .ge. 0) then
                    call PetscSectionGetDof(edirSection, iface, numDof, ierr); call CHKERR(ierr)
                    if (numDof .gt. 0) then
                      call PetscSectionGetOffset(edirsection, iface, voff, ierr); call CHKERR(ierr)
                      is_data_in(l) = i1d
                      is_data_out(l) = voff
                      l = l + 1
                    end if
                  end if
                end do

              end do

              ! and at the surface
              k = Cdir%glob_zm - 1
              ak = atmk(solver%atm, k)
              call pprts_cell_to_plex_cell_idx(Ca, [ak - 1, i, j], ri%plex, plex_cells, ierr); call CHKERR(ierr)

              ! on top faces:
              i1d = ind_nd_to_1d(pprts_offsets, [i0, k, i, j], cstyle=.true.)

              do m = 1, 2
                iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_BOT_FACE)
                if (iface .ge. 0) then
                  call PetscSectionGetOffset(edirsection, iface, voff, ierr); call CHKERR(ierr)
                  is_data_in(l) = i1d
                  is_data_out(l) = voff
                  l = l + 1
                end if
              end do
            end do
          end do
          is_data_size = l
        end if

        call ISCreateGeneral(PETSC_COMM_SELF, is_data_size, is_data_in, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
        call ISCreateGeneral(PETSC_COMM_SELF, is_data_size, is_data_out, PETSC_USE_POINTER, is_out, ierr); call CHKERR(ierr)

        call DMDAGetAO(Cdir%da, dmda_ao, ierr); call CHKERR(ierr)
        call AOApplicationToPetscIS(dmda_ao, is_in, ierr); call CHKERR(ierr)

        call PetscObjectSetName(is_in, "rayli_dir_iss_pprts_idx", ierr); call CHKERR(ierr)
        call PetscObjectSetName(is_out, "rayli_dir_iss_plex_idx", ierr); call CHKERR(ierr)

        call PetscObjectViewFromOptions(is_in, PETSC_NULL_IS, '-show_rayli_dir_iss', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(is_out, PETSC_NULL_IS, '-show_rayli_dir_iss', ierr); call CHKERR(ierr)

        call gen_shared_scatter_ctx(solution%edir, rayli_info%plex_solution%edir, rayli_info%ctx_edir, ierr, &
          & is_in, is_out); call CHKERR(ierr)

        call ISDestroy(is_out, ierr); call CHKERR(ierr)
        call ISDestroy(is_in, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    subroutine setup_ediff_scatter_context()
      integer(iintegers) :: i, j, k, voff, idof, i1d, iface, l, m, ak, numDof
      integer(iintegers) :: pprts_offsets(4), plex_cells(2)

      AO :: dmda_ao
      type(tPetscSection) :: ediffsection
      type(tIS) :: is_in, is_out
      integer(iintegers) :: is_data_size
      integer(iintegers), allocatable :: is_data_in(:)
      integer(iintegers), allocatable :: is_data_out(:)

      associate (ri => rayli_info, Cdiff => solver%C_diff, C1 => solver%C_one, Ca => solver%C_one_atm)
        if (submyid .eq. 0) then
          is_data_size = Cdiff%glob_zm * Cdiff%glob_xm * Cdiff%glob_ym * solver%difftop%dof * i2 & ! all horizontal fluxes
            & + C1%glob_zm * C1%glob_xm * C1%glob_ym * solver%diffside%dof * i2 ! plus fluxes on vertical faces once for x and y
        else
          is_data_size = 0
        end if

        allocate (is_data_in(0:is_data_size - 1))
        allocate (is_data_out(0:is_data_size - 1))
        is_data_in = -1
        is_data_out = -1
        if (submyid .eq. 0) then
          ! First do it for all the main dofs
          call ndarray_offsets(      &
            &[Cdiff%dof,      &
            & Cdiff%glob_zm,  &
            & Cdiff%glob_xm,  &
            & Cdiff%glob_ym], &
            & pprts_offsets)

          call DMGetSection(ri%plex%ediff_dm, ediffsection, ierr); call CHKERR(ierr)

          l = 0
          do j = 0, Cdiff%glob_ym - 1
            do i = 0, Cdiff%glob_xm - 1
              do k = 0, Cdiff%glob_zm - 2
                if (k .eq. 0) then
                  ak = 0
                else
                  ak = atmk(solver%atm, k)
                end if
                call pprts_cell_to_plex_cell_idx(Ca, [ak, i, j], ri%plex, plex_cells, ierr); call CHKERR(ierr)

                ! on top faces:
                do idof = 0, solver%difftop%dof - 1
                  i1d = ind_nd_to_1d(pprts_offsets, [idof, k, i, j], cstyle=.true.) ! i1d points to dof on pprts topface
                  do m = 1, 2
                    iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_TOP_FACE)
                    if (iface .ge. 0) then ! have a top face
                      call PetscSectionGetFieldOffset(ediffsection, iface, i0, voff, ierr); call CHKERR(ierr)
                      if (solver%difftop%is_inward(i1 + idof)) then ! in pprts: edn
                        is_data_in(l) = i1d
                        is_data_out(l) = voff
                      else
                        is_data_in(l) = i1d
                        is_data_out(l) = voff + 1
                      end if
                      l = l + 1
                    end if
                  end do
                end do

                !! on x-side faces:
                do idof = 0, solver%diffside%dof - 1
                  ! i1d points to dof on pprts sideface
                  i1d = ind_nd_to_1d(pprts_offsets, [solver%difftop%dof + idof, k, i, j], cstyle=.true.)
                  do m = 1, 2
                    iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_LEFT_FACE)
                    if (iface .ge. 0) then ! have a left face
                      call PetscSectionGetFieldOffset(ediffsection, iface, i1, voff, ierr); call CHKERR(ierr)
                      call PetscSectionGetDof(ediffSection, iface, numDof, ierr); call CHKERR(ierr)
                      if (numDof .gt. 0) then
                        if (solver%diffside%is_inward(i1 + idof)) then ! in pprts: e_right
                          is_data_in(l) = i1d
                          is_data_out(l) = voff
                        else
                          is_data_in(l) = i1d
                          is_data_out(l) = voff + 1
                        end if
                        l = l + 1
                      end if
                    end if
                  end do
                end do

                !! on y-side faces:
                do idof = 0, solver%diffside%dof - 1
                  ! i1d points to dof on pprts sideface
                  i1d = ind_nd_to_1d(pprts_offsets, &
                    & [solver%difftop%dof + solver%diffside%dof + idof, k, i, j], cstyle=.true.)
                  do m = 1, 2
                    iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_REAR_FACE)
                    if (iface .ge. 0) then ! have a rear face
                      call PetscSectionGetFieldOffset(ediffsection, iface, i1, voff, ierr); call CHKERR(ierr)
                      call PetscSectionGetDof(ediffSection, iface, numDof, ierr); call CHKERR(ierr)
                      if (numDof .gt. 0) then
                        if (solver%diffside%is_inward(i1 + idof)) then ! in pprts: e_forward
                          is_data_in(l) = i1d
                          is_data_out(l) = voff
                        else
                          is_data_in(l) = i1d
                          is_data_out(l) = voff + 1
                        end if
                        l = l + 1
                      end if
                    end if
                  end do
                end do

              end do ! k

              ! and at the surface
              k = Cdiff%glob_zm - 1
              ak = atmk(solver%atm, k)
              call pprts_cell_to_plex_cell_idx(Ca, [ak - 1, i, j], ri%plex, plex_cells, ierr); call CHKERR(ierr)

              do idof = 0, solver%difftop%dof - 1
                ! on top faces:
                i1d = ind_nd_to_1d(pprts_offsets, [idof, k, i, j], cstyle=.true.)
                do m = 1, 2
                  iface = find_face_idx_by_orientation(ri%plex, plex_cells(m), PPRTS_BOT_FACE)
                  if (iface .ge. 0) then
                    call PetscSectionGetFieldOffset(ediffsection, iface, i0, voff, ierr); call CHKERR(ierr)
                    is_data_in(l) = i1d
                    if (solver%difftop%is_inward(i1 + idof)) then ! in pprts: edn
                      is_data_out(l) = voff + 1
                    else
                      is_data_out(l) = voff
                    end if
                    l = l + 1
                  end if
                end do
              end do
            end do
          end do
          is_data_size = l
        end if

        call ISCreateGeneral(PETSC_COMM_SELF, is_data_size, is_data_in, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
        call ISCreateGeneral(PETSC_COMM_SELF, is_data_size, is_data_out, PETSC_USE_POINTER, is_out, ierr); call CHKERR(ierr)

        call DMDAGetAO(Cdiff%da, dmda_ao, ierr); call CHKERR(ierr)
        call AOApplicationToPetscIS(dmda_ao, is_in, ierr); call CHKERR(ierr)

        call PetscObjectSetName(is_in, "rayli_diff_iss_pprts_idx", ierr); call CHKERR(ierr)
        call PetscObjectSetName(is_out, "rayli_diff_iss_plex_idx", ierr); call CHKERR(ierr)

        call PetscObjectViewFromOptions(is_in, PETSC_NULL_IS, '-show_rayli_diff_iss', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(is_out, PETSC_NULL_IS, '-show_rayli_diff_iss', ierr); call CHKERR(ierr)
        call gen_shared_scatter_ctx(solution%ediff, rayli_info%plex_solution%ediff, rayli_info%ctx_ediff, ierr, &
          & is_in, is_out); call CHKERR(ierr)
        call ISDestroy(is_out, ierr); call CHKERR(ierr)
        call ISDestroy(is_in, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    subroutine setup_planck_scatter_context()
      integer(iintegers) :: Nplanck, i, k
      AO :: dmda_ao
      type(tIS) :: is_in, is_out
      integer(iintegers), allocatable :: is_data(:)
      type(tVec) :: vec_planck ! tmp vec to create scatter ctx

      associate ( &
          & ri => rayli_info, &
          & Ca1 => solver%C_one_atm1 &
          )

        if (allocated(ri%ctx_planck)) return
        allocate (ri%ctx_planck)

        call DMGetGlobalVector(Ca1%da, vec_planck, ierr); call CHKERR(ierr)
        if (submyid .eq. 0) then
          call VecGetSize(vec_planck, Nplanck, ierr); call CHKERR(ierr)
        else
          Nplanck = 0
        end if
        call VecCreateSeq(PETSC_COMM_SELF, Nplanck * i2, ri%planck, ierr); call CHKERR(ierr)
        allocate (is_data(i0:Nplanck * i2 - i1))
        if (submyid .eq. 0) then
          do i = i0, Ca1%glob_xm * Ca1%glob_ym - i1
            do k = i0, Ca1%glob_zm - i1
              is_data((i * i2) * Ca1%glob_zm + k) = i * Ca1%glob_zm + k
              is_data((i * i2 + 1) * Ca1%glob_zm + k) = i * Ca1%glob_zm + k
            end do
          end do
        end if
        call ISCreateGeneral(PETSC_COMM_SELF, Nplanck * i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
        call ISCreateStride(PETSC_COMM_SELF, Nplanck * i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)

        call DMDAGetAO(Ca1%da, dmda_ao, ierr); call CHKERR(ierr)
        call AOApplicationToPetscIS(dmda_ao, is_in, ierr); call CHKERR(ierr)

        call gen_shared_scatter_ctx(vec_planck, ri%planck, ri%ctx_planck, ierr, is_in, is_out); call CHKERR(ierr)
        call ISDestroy(is_out, ierr); call CHKERR(ierr)
        call ISDestroy(is_in, ierr); call CHKERR(ierr)
        deallocate (is_data)
        call DMRestoreGlobalVector(Ca1%da, vec_planck, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    subroutine setup_abso_scatter_context()
      integer(iintegers) :: is_data_size
      integer(iintegers) :: pprts_offsets(4), plex_cells(2)
      integer(iintegers) :: i, j, k, l, voff, i1d, m, ak

      AO :: dmda_ao
      type(tPetscSection) :: cellSection
      type(tIS) :: is_in, is_out
      integer(iintegers), allocatable :: is_data_in(:)
      integer(iintegers), allocatable :: is_data_out(:)

      associate (ri => rayli_info, C1 => solver%C_one, Ca => solver%C_one_atm)
        if (submyid .eq. 0) then
          is_data_size = Ca%glob_zm * C1%glob_xm * C1%glob_ym * i2
        else
          is_data_size = 0
        end if

        allocate (is_data_in(0:is_data_size - 1))
        allocate (is_data_out(0:is_data_size - 1))
        is_data_in = -1
        is_data_out = -1

        if (submyid .eq. 0) then
          ! First do it for all the main dofs
          call ndarray_offsets(&
            &[i1,          &
            & C1%glob_zm,  &
            & C1%glob_xm,  &
            & C1%glob_ym], &
            & pprts_offsets)

          if (.not. allocated(ri%plex%cell1_dm)) call CHKERR(1_mpiint, 'ri%plex%cell1_dm not allocated but needed')
          call DMGetSection(ri%plex%cell1_dm, cellSection, ierr); call CHKERR(ierr)

          l = 0
          do j = 0, C1%glob_ym - 1
            do i = 0, C1%glob_xm - 1
              ! First build sum of upper atmosphere and put them in the 0th layer, i.e. sum up the collapsed values
              do ak = 0, solver%atm%icollapse - 1
                k = 0
                call pprts_cell_to_plex_cell_idx(Ca, [ak, i, j], ri%plex, plex_cells, ierr); call CHKERR(ierr)
                i1d = ind_nd_to_1d(pprts_offsets, [i0, k, i, j], cstyle=.true.) ! i1d points to dof on pprts topface
                do m = 1, 2
                  call PetscSectionGetOffset(cellSection, plex_cells(m), voff, ierr); call CHKERR(ierr)
                  is_data_in(l) = i1d
                  is_data_out(l) = voff
                  l = l + 1
                end do
              end do
              ! then define the normal layers
              do k = i1, C1%glob_zm - 1
                ak = atmk(solver%atm, k)
                call pprts_cell_to_plex_cell_idx(Ca, [ak, i, j], ri%plex, plex_cells, ierr); call CHKERR(ierr)
                i1d = ind_nd_to_1d(pprts_offsets, [i0, k, i, j], cstyle=.true.) ! i1d points to dof on pprts cell
                do m = 1, 2
                  call PetscSectionGetOffset(cellSection, plex_cells(m), voff, ierr); call CHKERR(ierr)
                  is_data_in(l) = i1d
                  is_data_out(l) = voff
                  l = l + 1
                end do
              end do
            end do
          end do
          is_data_size = l
        end if

        call ISCreateGeneral(PETSC_COMM_SELF, is_data_size, is_data_in, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
        call ISCreateGeneral(PETSC_COMM_SELF, is_data_size, is_data_out, PETSC_USE_POINTER, is_out, ierr); call CHKERR(ierr)

        call DMDAGetAO(C1%da, dmda_ao, ierr); call CHKERR(ierr)
        call AOApplicationToPetscIS(dmda_ao, is_in, ierr); call CHKERR(ierr)

        call PetscObjectSetName(is_in, "rayli_abso_iss_pprts_idx", ierr); call CHKERR(ierr)
        call PetscObjectSetName(is_out, "rayli_abso_iss_plex_idx", ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(is_in, PETSC_NULL_IS, '-show_rayli_abso_iss', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(is_out, PETSC_NULL_IS, '-show_rayli_abso_iss', ierr); call CHKERR(ierr)

        call gen_shared_scatter_ctx(solution%abso, rayli_info%plex_solution%abso, rayli_info%ctx_abso, ierr, &
          & is_in, is_out); call CHKERR(ierr)
        call ISDestroy(is_out, ierr); call CHKERR(ierr)
        call ISDestroy(is_in, ierr); call CHKERR(ierr)
      end associate
    end subroutine

  end subroutine

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine twostream(solver, edirTOA, solution, opt_buildings)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container) :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    real(ireals), pointer, dimension(:, :, :, :) :: xv_dir => null(), xv_diff => null(), xv_abso => null()
    real(ireals), pointer, dimension(:) :: xv_dir1d => null(), xv_diff1d => null(), xv_abso1d => null()
    integer(iintegers) :: i, j, k, src

    real(ireals), allocatable :: dtau(:), kext(:), w0(:), g(:), S(:), Edn(:), Eup(:)
    real(ireals) :: mu0, incSolar, fac, Ag, Bsrfc
    integer(mpiint) :: ierr

    integer(iintegers), allocatable :: buildings_mink(:, :)
    integer(iintegers), allocatable :: buildings_face(:, :)
    integer(iintegers) :: m, bk, idx(4)

    associate (atm => solver%atm, &
               C_diff => solver%C_diff, &
               C_dir => solver%C_dir, &
               C_one => solver%C_one, &
               C_one_atm => solver%C_one_atm, &
               C_one_atm1 => solver%C_one_atm1)

      if (present(opt_buildings)) then
        associate (B => opt_buildings)
          ! on each pixel, find the top most face
          allocate (buildings_mink(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          allocate (buildings_face(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          buildings_mink(:, :) = C_dir%ze
          buildings_face(:, :) = -1
          do m = 1, size(opt_buildings%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) - 1 + [C_dir%zs, C_dir%xs, C_dir%ys]

            associate (d => idx(1), k => idx(2), i => idx(3), j => idx(4))
              if (d .eq. PPRTS_TOP_FACE) then
                if (k .lt. buildings_mink(i, j)) then
                  buildings_mink(i, j) = k
                  buildings_face(i, j) = m
                end if
              end if
            end associate
          end do
        end associate
      end if

      if (solution%lsolar_rad) then
        call PetscObjectSetName(solution%edir, 'twostream_edir_vec uid='//toStr(solution%uid), ierr); call CHKERR(ierr)
        call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
      end if

      call PetscObjectSetName(solution%ediff, 'twostream_ediff_vec uid='//toStr(solution%uid), ierr); call CHKERR(ierr)
      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

      allocate (dtau(C_one_atm%zs:C_one_atm%ze))
      allocate (kext(C_one_atm%zs:C_one_atm%ze))
      allocate (w0(C_one_atm%zs:C_one_atm%ze))
      allocate (g(C_one_atm%zs:C_one_atm%ze))

      if (solution%lsolar_rad) then
        call getVecPointer(C_dir%da, solution%edir, xv_dir1d, xv_dir)
        mu0 = real(solver%sun%costheta)
      else
        mu0 = 0
      end if
      incSolar = edirTOA * mu0

      call getVecPointer(C_diff%da, solution%ediff, xv_diff1d, xv_diff)
      call getVecPointer(C_one%da, solution%abso, xv_abso1d, xv_abso)

      allocate (S(C_one_atm1%zs:C_one_atm1%ze))
      allocate (Eup(C_one_atm1%zs:C_one_atm1%ze))
      allocate (Edn(C_one_atm1%zs:C_one_atm1%ze))

      do j = C_one_atm%ys, C_one_atm%ye
        do i = C_one_atm%xs, C_one_atm%xe

          kext = atm%kabs(:, i, j) + atm%ksca(:, i, j)
          dtau = atm%dz(:, i, j) * kext
          w0 = atm%ksca(:, i, j) / max(kext, epsilon(kext))
          g = atm%g(:, i, j)

          Ag = atm%albedo(i, j)

          if (allocated(atm%planck)) then
            if (allocated(atm%Bsrfc)) then
              Bsrfc = atm%Bsrfc(i, j)
            else
              Bsrfc = atm%planck(C_one_atm1%ze, i, j)
            end if
          end if

          if (present(opt_buildings)) then
            m = buildings_face(i, j)
            if (m .gt. i0) then
              Ag = opt_buildings%albedo(m)
              bk = atmk(atm, buildings_mink(i, j) - 1)
              if (allocated(atm%planck)) then
                Bsrfc = opt_buildings%planck(m)
              end if
            else
              bk = C_one_atm%ze
            end if

            S(:) = zero
            Edn(:) = zero
            Eup(:) = zero

            if (allocated(atm%planck)) then
              call delta_eddington_twostream(&
                & dtau(C_one_atm%zs:bk), &
                & w0(C_one_atm%zs:bk), &
                & g(C_one_atm%zs:bk), &
                & mu0, incSolar, Ag, &
                & S(C_one_atm%zs:bk + 1), &
                & Edn(C_one_atm%zs:bk + 1), &
                & Eup(C_one_atm%zs:bk + 1), &
                & planck=atm%planck(C_one_atm%zs:bk + 1, i, j), &
                & planck_srfc=Bsrfc)
            else

              call delta_eddington_twostream(&
                & dtau(C_one_atm%zs:bk), &
                & w0(C_one_atm%zs:bk), &
                & g(C_one_atm%zs:bk), &
                & mu0, incSolar, Ag, &
                & S(C_one_atm%zs:bk + 1), &
                & Edn(C_one_atm%zs:bk + 1), &
                & Eup(C_one_atm%zs:bk + 1))
            end if

!            print *,'twostream i,j',i,j,'bk',bk, m
!            if(i.eq.2.and.j.eq.2) then
!              do bk=lbound(Edn,1),ubound(Edn,1)
!                print *,i,j,'k',bk,'edn',Edn(bk),'eup',Eup(bk)
!              enddo
!            endif

          else ! not buildings

            if (allocated(atm%planck)) then
              call delta_eddington_twostream(dtau, w0, g, &
                                             mu0, incSolar, Ag, &
                                             S, Edn, Eup, &
                                             planck=atm%planck(:, i, j), &
                                             planck_srfc=Bsrfc)
            else
              call adding_delta_eddington_twostream(dtau, w0, g, mu0, incSolar, atm%albedo(i, j), S, Edn, Eup)

              !call delta_eddington_twostream(&
              !  & dtau, w0, g,&
              !  & mu0, incSolar, Ag, &
              !  & S, Edn, Eup)
            end if
          end if

          if (solution%lsolar_rad) then
            fac = real(solver%dirtop%area_divider, ireals) / real(solver%dirtop%streams, ireals)
            do src = i0, solver%dirtop%dof - 1
              xv_dir(src, C_dir%zs + 1:C_dir%ze, i, j) = S(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_dir(src, C_dir%zs, i, j) = S(C_one_atm1%zs) * fac
            end do
          end if

          fac = real(solver%difftop%area_divider, ireals) / real(solver%difftop%streams, ireals)
          do src = 1, solver%difftop%dof
            if (solver%difftop%is_inward(src)) then
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = Edn(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = Edn(C_one_atm1%zs) * fac
            else
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = Eup(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = Eup(C_one_atm1%zs) * fac
            end if
          end do

          xv_abso(i0, :, i, j) = &
            & +Edn(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & - Edn(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) &
            & - Eup(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & + Eup(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)

          if (solution%lsolar_rad) then
            xv_abso(i0, :, i, j) = xv_abso(i0, :, i, j) &
              & + S(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
              & - S(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)
          end if
          do k = C_one%zs, C_one%ze
            xv_abso(i0, k, i, j) = xv_abso(i0, k, i, j) / atm%dz(atmk(atm, k), i, j)
          end do
        end do
      end do

      if (solution%lsolar_rad) &
        call restoreVecPointer(C_dir%da, solution%edir, xv_dir1d, xv_dir)
      call restoreVecPointer(C_diff%da, solution%ediff, xv_diff1d, xv_diff)
      call restoreVecPointer(C_one%da, solution%abso, xv_abso1d, xv_abso)

      !Twostream solver returns fluxes as [W]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is not up to date
      solution%lchanged = .false.

      deallocate (S)
      deallocate (Edn)
      deallocate (Eup)

    end associate
  end subroutine

  !> @brief wrapper for the disort solver
  subroutine disort(solver, edirTOA, solution, opt_buildings)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container) :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    real(ireals), pointer, dimension(:, :, :, :) :: xv_dir => null(), xv_diff => null()
    real(ireals), pointer, dimension(:) :: xv_dir1d => null(), xv_diff1d => null()
    integer(iintegers) :: i, j, src, nstreams

    real, allocatable :: &
      & Bfrac(:), &
      & Blev(:), &
      & tlev(:), &
      & kext(:), &
      & dtau(:), &
      & w0(:), &
      & g(:), &
      & FLDIR(:),&
      & FLDN(:), &
      & FLUP(:), &
      & DFDT(:), &
      & UAVG(:)
    real :: mu0, Ag, Bskin
    real(ireals) :: fac
    integer(mpiint) :: ierr
    logical :: lflg

    if (present(opt_buildings)) call CHKERR(1_mpiint, "buildings not implemented for pprts_disort")

    associate (atm => solver%atm, &
               C_diff => solver%C_diff, &
               C_dir => solver%C_dir, &
               C_one_atm => solver%C_one_atm, &
               C_one_atm1 => solver%C_one_atm1)

      nstreams = 16
      call get_petsc_opt(solver%prefix, &
                         "-disort_streams", nstreams, lflg, ierr); call CHKERR(ierr)

      if (solution%lsolar_rad) then
        call PetscObjectSetName(solution%edir, 'twostream_edir_vec uid='//toStr(solution%uid), ierr); call CHKERR(ierr)
        call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
      end if

      call PetscObjectSetName(solution%ediff, 'twostream_ediff_vec uid='//toStr(solution%uid), ierr); call CHKERR(ierr)
      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

      allocate (dtau(C_one_atm%zs:C_one_atm%ze))
      allocate (kext(C_one_atm%zs:C_one_atm%ze))
      allocate (w0(C_one_atm%zs:C_one_atm%ze))
      allocate (g(C_one_atm%zs:C_one_atm%ze))
      allocate (Bfrac(C_one_atm%zs:C_one_atm%ze))
      allocate (tlev(C_one_atm1%zs:C_one_atm1%ze))
      Bfrac = 1
      tlev = 300

      if (solution%lsolar_rad) then
        call getVecPointer(C_dir%da, solution%edir, xv_dir1d, xv_dir)
        mu0 = real(solver%sun%costheta)
      else
        mu0 = 0
        allocate (Blev(C_one_atm%zs:C_one_atm%ze))
      end if

      call getVecPointer(C_diff%da, solution%ediff, xv_diff1d, xv_diff)

      allocate (FLDIR(C_one_atm1%zs:C_one_atm1%ze))
      allocate (FLDN(C_one_atm1%zs:C_one_atm1%ze))
      allocate (FLUP(C_one_atm1%zs:C_one_atm1%ze))
      allocate (DFDT(C_one_atm1%zs:C_one_atm1%ze))
      allocate (UAVG(C_one_atm1%zs:C_one_atm1%ze))

      do j = C_one_atm%ys, C_one_atm%ye
        do i = C_one_atm%xs, C_one_atm%xe

          kext = real(atm%kabs(:, i, j) + atm%ksca(:, i, j))
          dtau = real(atm%dz(:, i, j) * kext)
          w0 = real(atm%ksca(:, i, j) / max(kext, epsilon(kext)))
          g = real(atm%g(:, i, j))
          Ag = real(atm%albedo(i, j))
          if (.not. solution%lsolar_rad) then
            Blev = real(atm%planck(:, i, j))
            if (allocated(atm%Bsrfc)) then
              Bskin = real(atm%Bsrfc(i, j))
            else
              Bskin = real(atm%planck(ubound(atm%planck, 1), i, j))
            end if
          end if

          call default_flx_computation(       &
            & mu0,                            &
            & real(max(0._ireals, edirTOA)),  &
            & Ag,                             &
            & 300.,                           & ! tskin (ignored because we provide planck values directly)
            & .not. solution%lsolar_rad,      & ! lthermal
            & [0., 1.],                       & ! wavenumbers (ignored because we provide planck values directly)
            & Bfrac,                          &
            & dtau,                           &
            & w0,                             &
            & g,                              &
            & tlev,                           & ! (ignored because we provide planck values directly)
            & FLDIR, FLDN, FLUP, DFDT, UAVG,  &
            & int(nstreams), lverbose=.false.,&
            & Blev=Blev, Bskin=Bskin)

          if (solution%lsolar_rad) then
            fac = real(solver%dirtop%area_divider, ireals) / real(solver%dirtop%streams, ireals)
            do src = i0, solver%dirtop%dof - 1
              xv_dir(src, C_dir%zs + 1:C_dir%ze, i, j) = FLDIR(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_dir(src, C_dir%zs, i, j) = FLDIR(C_one_atm1%zs) * fac
            end do
          end if

          fac = real(solver%difftop%area_divider, ireals) / real(solver%difftop%streams, ireals)
          do src = 1, solver%difftop%dof
            if (solver%difftop%is_inward(src)) then
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = FLDN(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = FLDN(C_one_atm1%zs) * fac
            else
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = FLUP(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = FLUP(C_one_atm1%zs) * fac
            end if
          end do
        end do
      end do

      if (solution%lsolar_rad) &
        call restoreVecPointer(C_dir%da, solution%edir, xv_dir1d, xv_dir)
      call restoreVecPointer(C_diff%da, solution%ediff, xv_diff1d, xv_diff)

      !Disort returns fluxes as [W]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is not up to date
      solution%lchanged = .true.

    end associate
  end subroutine

  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine schwarz(solver, solution, opt_buildings)
    class(t_solver) :: solver
    type(t_state_container) :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    real(ireals), pointer, dimension(:, :, :, :) :: xv_diff => null(), xv_abso => null()
    real(ireals), pointer, dimension(:) :: xv_diff1d => null(), xv_abso1d => null()
    integer(iintegers) :: i, j, k, idof
    integer(iintegers) :: Nmu, ak
    logical :: lflg

    real(ireals), allocatable :: dtau(:), Edn(:), Eup(:)
    real(ireals) :: Bsrfc, Ag
    integer(mpiint) :: ierr

    integer(iintegers), allocatable :: buildings_mink(:, :)
    integer(iintegers), allocatable :: buildings_face(:, :)
    integer(iintegers) :: m, bk, idx(4)

    associate ( &
      atm => solver%atm, &
      C_diff => solver%C_diff, &
      C_one => solver%C_one, &
      C_one1 => solver%C_one1, &
      C_one_atm => solver%C_one_atm, &
      C_one_atm1 => solver%C_one_atm1)

      if (present(opt_buildings)) then
        associate (B => opt_buildings)
          ! on each pixel, find the top most face
          allocate (buildings_mink(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          allocate (buildings_face(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          buildings_mink(:, :) = C_diff%ze
          buildings_face(:, :) = -1
          do m = 1, size(B%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) - 1 + [C_diff%zs, C_diff%xs, C_diff%ys]

            associate (d => idx(1), k => idx(2), i => idx(3), j => idx(4))
              if (d .eq. PPRTS_TOP_FACE) then
                if (k .lt. buildings_mink(i, j)) then
                  buildings_mink(i, j) = k
                  buildings_face(i, j) = m
                end if
              end if
            end associate
          end do
        end associate
      end if

      if (solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling schwarschild solver for solar calculation -- stopping!')
      if (.not. allocated(atm%planck)) &
        & call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

      allocate (dtau(C_one_atm%zs:C_one_atm%ze))

      call getVecPointer(C_one%da, solution%abso, xv_abso1d, xv_abso)
      call getVecPointer(C_diff%da, solution%ediff, xv_diff1d, xv_diff)

      allocate (Eup(C_one_atm1%zs:C_one_atm1%ze))
      allocate (Edn(C_one_atm1%zs:C_one_atm1%ze))

      Nmu = 4
      call get_petsc_opt(solver%prefix, &
                         "-schwarzschild_Nmu", Nmu, lflg, ierr); call CHKERR(ierr)

      if (solver%myid .eq. 0 .and. ldebug) print *, ' CALCULATING schwarzschild ::'

      do j = C_diff%ys, C_diff%ye
        do i = C_diff%xs, C_diff%xe

          dtau = atm%dz(:, i, j) * atm%kabs(:, i, j)
          Ag = atm%albedo(i, j)

          if (allocated(atm%Bsrfc)) then
            Bsrfc = atm%Bsrfc(i, j)
          else
            Bsrfc = atm%planck(C_one_atm1%ze, i, j)
          end if

          if (present(opt_buildings)) then
            m = buildings_face(i, j)
            if (m .gt. i0) then
              Ag = opt_buildings%albedo(m)
              bk = atmk(atm, buildings_mink(i, j) - 1)
              if (allocated(atm%planck)) then
                Bsrfc = opt_buildings%planck(m)
              end if
            else
              bk = C_one_atm%ze
            end if
            call schwarzschild(          &
              & Nmu,                     &
              & dtau(C_one_atm%zs:bk),   &
              & Ag,                      &
              & Edn(C_one_atm%zs:bk + 1),  &
              & Eup(C_one_atm%zs:bk + 1),  &
              & atm%planck(C_one_atm%zs:bk + 1, i, j), &
              & opt_srfc_emission=Bsrfc)
          else ! not buildings
            call schwarzschild(          &
              & Nmu,                     &
              & dtau,                    &
              & Ag,                      &
              & Edn, Eup,                &
              & atm%planck(:, i, j),     &
              & opt_srfc_emission=Bsrfc)
          end if

          ! icollapse needs special case for TOA flx's
          do idof = 0, solver%difftop%dof - 1
            if (solver%difftop%is_inward(i1 + idof)) then ! Edn
              xv_diff(idof, C_diff%zs, i, j) = Edn(0) / real(solver%difftop%streams, ireals)
            else ! Eup
              xv_diff(idof, C_diff%zs, i, j) = Eup(0) / real(solver%difftop%streams, ireals)
            end if
          end do

          ! rest of the atmosphere
          do k = C_diff%zs + 1, C_diff%ze
            ak = atmk(atm, k)
            do idof = 0, solver%difftop%dof - 1
              if (solver%difftop%is_inward(i1 + idof)) then ! Edn
                xv_diff(idof, k, i, j) = Edn(ak) / real(solver%difftop%streams, ireals)
              else ! Eup
                xv_diff(idof, k, i, j) = Eup(ak) / real(solver%difftop%streams, ireals)
              end if
            end do
          end do

          xv_abso(i0, :, i, j) = &
            & +Edn(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & - Edn(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) &
            & - Eup(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & + Eup(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)

          do k = C_one%zs, C_one%ze
            xv_abso(i0, k, i, j) = xv_abso(i0, k, i, j) / atm%dz(atmk(atm, k), i, j)
          end do

        end do
      end do

      call restoreVecPointer(C_diff%da, solution%ediff, xv_diff1d, xv_diff)
      call restoreVecPointer(C_one%da, solution%abso, xv_abso1d, xv_abso)

      !Schwarzschild solver returns fluxes as [W/m^2]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is up to date
      solution%lchanged = .false.

      deallocate (Edn)
      deallocate (Eup)

    end associate
  end subroutine
end module
