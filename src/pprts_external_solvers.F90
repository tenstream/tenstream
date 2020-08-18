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

  use m_data_parameters, only : ireals, iintegers, mpiint, imp_iinteger, &
    & i0, i1, i2, i3, i4, zero, one
  use m_helper_functions, only : CHKERR, CHKWARN, itoa, ndarray_offsets, ind_nd_to_1d, &
    meanval, spherical_2_cartesian, imp_allreduce_sum
  use m_petsc_helpers, only : getVecPointer, restoreVecPointer, &
    & f90VecToPetsc, &
    & gen_shared_subcomm, gen_shared_scatter_ctx

  use m_pprts_base, only : t_solver, &
    & t_state_container, prepare_solution, destroy_solution, &
    & atmk, set_dmda_cell_coordinates, &
    & interpolate_cell_values_to_vertices
  use m_schwarzschild, only: schwarzschild, B_eff
  use m_twostream, only: delta_eddington_twostream, adding_delta_eddington_twostream

  use m_icon_plex_utils, only: create_2d_regular_plex, dmplex_2D_to_3D
  use m_plex_grid, only: t_plexgrid, destroy_plexgrid, &
    & setup_plexgrid, setup_edir_dmplex, setup_ediff_dmplex, setup_abso_dmplex
  use m_plex2rayli, only: rayli_wrapper

  implicit none
  private
  public :: twostream, schwarz, pprts_rayli_wrapper, destroy_rayli_info

  logical,parameter :: ldebug=.True.

  type t_rayli_info
    type(t_state_container) :: plex_solution
    type(t_plexgrid), allocatable :: plex
    integer(mpiint) :: subcomm
    integer(iintegers) :: num_subcomm_masters
    type(tVecScatter) :: ctx_hhl
    type(tVecScatter) :: ctx_albedo
    type(tVecScatter) :: ctx_optprop
    type(tVecScatter) :: ctx_edir
    type(tVecScatter) :: ctx_ediff
    type(tVecScatter) :: ctx_abso
    type(tVec) :: albedo, kabs, ksca, g
  end type

  type(t_rayli_info), allocatable :: rayli_info

contains

  subroutine destroy_rayli_info()
    integer(mpiint) :: ierr
    if(.not.allocated(rayli_info)) return
    associate(ri => rayli_info)
      call destroy_solution(ri%plex_solution)
      call destroy_plexgrid(ri%plex)
      call mpi_comm_free(ri%subcomm, ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_hhl    , ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_albedo , ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_optprop, ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_edir   , ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_ediff  , ierr); call CHKERR(ierr)
      call VecScatterDestroy(ri%ctx_abso   , ierr); call CHKERR(ierr)
      call VecDestroy(ri%albedo, ierr); call CHKERR(ierr)
      call VecDestroy(ri%kabs  , ierr); call CHKERR(ierr)
      call VecDestroy(ri%ksca  , ierr); call CHKERR(ierr)
      call VecDestroy(ri%g     , ierr); call CHKERR(ierr)
    end associate

    deallocate(rayli_info)
  end subroutine

  subroutine init_pprts_rayli_wrapper(solver, solution, rayli_info)
    class(t_solver), intent(in) :: solver
    type(t_state_container),intent(in) :: solution
    type(t_rayli_info), allocatable, intent(inout) :: rayli_info
    integer(mpiint) :: submyid, ierr

    integer(iintegers) :: Nhhl, Nalbedo, Noptprop, Nedir, Nediff, Nabso
    type(tVec) :: vertex_hhl, hhl
    type(tVec) :: vec_albedo, vec_optprop
    real(ireals), pointer :: xhhl(:,:,:,:) => null(), xhhl1d(:) => null()
    type(tPetscSection) :: coord_section
    type(tVec) :: coordinates
    real(ireals), pointer :: coords(:)
    integer(iintegers) :: i, k, kk, vStart, vEnd, ivert, coord_offset, idof, voff, off_plex(2)

    type(tDM) :: dm2d, dm2d_dist, dm3d
    integer(iintegers), allocatable :: zindex(:)

    type(tIS) :: is_in, is_out
    integer(iintegers), allocatable :: is_data(:)

    if(allocated(rayli_info)) return
    allocate(rayli_info)

    ! Scatter height info to subcomm[0]
    call gen_shared_subcomm(solver%comm, rayli_info%subcomm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(rayli_info%subcomm,submyid,ierr); call CHKERR(ierr)

    if(submyid.eq.0) then
      k=1
    else
      k=0
    endif
    call imp_allreduce_sum(solver%comm, k, rayli_info%num_subcomm_masters)

    call DMGetGlobalVector(solver%Cvert_one_atm1%da, vertex_hhl, ierr); call CHKERR(ierr)

    call interpolate_cell_values_to_vertices(&
      & solver%C_one_atm1_box, solver%atm%hhl, &
      & solver%Cvert_one_atm1, vertex_hhl)

    call PetscObjectViewFromOptions(vertex_hhl, PETSC_NULL_VEC, '-show_rayli_hhl', ierr); call CHKERR(ierr)

    if(submyid.eq.0) then
      call VecGetSize(vertex_hhl, Nhhl, ierr); call CHKERR(ierr)
    else
      Nhhl = 0
    endif

    call VecCreateSeq(PETSC_COMM_SELF, Nhhl, hhl, ierr); call CHKERR(ierr)
    call gen_shared_scatter_ctx(vertex_hhl, hhl, rayli_info%ctx_hhl, ierr); call CHKERR(ierr)
    call VecScatterBegin(rayli_info%ctx_hhl, vertex_hhl, hhl, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
    call VecScatterEnd  (rayli_info%ctx_hhl, vertex_hhl, hhl, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

    call DMRestoreGlobalVector(solver%Cvert_one_atm1%da, vertex_hhl, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(hhl, PETSC_NULL_VEC, '-show_rayli_hhl', ierr); call CHKERR(ierr)

    associate( &
        & atm => solver%atm, &
        & Ca1 => solver%C_one_atm1, &
        & Ca  => solver%C_one_atm, &
        & Cv  => solver%Cvert_one_atm1 )

      if(submyid.eq.0) then
        call create_2d_regular_plex(PETSC_COMM_SELF, Ca%glob_xm+1, Ca%glob_ym+1, &
          & dm2d, dm2d_dist, opt_dx=atm%dx, opt_dy=atm%dy, lverbose=.True.)

        call VecGetArrayReadF90(hhl, xhhl1d, ierr); call CHKERR(ierr)
        xhhl(i1:i1, i1:Cv%glob_zm, i1:Cv%glob_xm, i1:Cv%glob_ym) => xhhl1d

        call dmplex_2D_to_3D(dm2d, solver%C_one_atm1%glob_zm, xhhl(i1, :, i1, i1), dm3d, zindex, lpolar_coords=.False.)

        !set height lvls on vertices
        call DMGetCoordinateSection(dm3d, coord_section, ierr); call CHKERR(ierr)
        call DMGetCoordinatesLocal(dm3d, coordinates, ierr); call CHKERR(ierr)
        call DMPlexGetDepthStratum(dm3d, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices
        call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)


        do i = 1, size(xhhl1d)
          ivert = vStart + i - i1
          call PetscSectionGetOffset(coord_section, ivert, coord_offset, ierr); call CHKERR(ierr)
          coords(i1+coord_offset+i2) = xhhl1d(i)
        enddo
        call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

        call setup_plexgrid(dm2d, dm3d, solver%C_one_atm%glob_zm, zindex, rayli_info%plex, xhhl(i1, :, i1, i1))

        call PetscObjectViewFromOptions(dm3d, PETSC_NULL_DM, '-show_rayli_dm3d', ierr); call CHKERR(ierr)
        nullify(xhhl)
        call VecRestoreArrayReadF90(hhl, xhhl1d, ierr); call CHKERR(ierr)

        call setup_edir_dmplex (rayli_info%plex, rayli_info%plex%dm, solver%dirtop%dof, i0, i1, rayli_info%plex%edir_dm)
        call setup_ediff_dmplex(rayli_info%plex, rayli_info%plex%dm, solver%difftop%dof/2, i0, i2, rayli_info%plex%ediff_dm)
        call setup_abso_dmplex (rayli_info%plex%dm, rayli_info%plex%abso_dm)

        call prepare_solution(rayli_info%plex%edir_dm, rayli_info%plex%ediff_dm, rayli_info%plex%abso_dm, &
          lsolar=solution%lsolar_rad, solution=rayli_info%plex_solution)

      else
        allocate(rayli_info%plex_solution%edir)
        call VecCreateSeq(PETSC_COMM_SELF, i0, rayli_info%plex_solution%edir, ierr); call CHKERR(ierr)
        allocate(rayli_info%plex_solution%ediff)
        call VecCreateSeq(PETSC_COMM_SELF, i0, rayli_info%plex_solution%ediff, ierr); call CHKERR(ierr)
        allocate(rayli_info%plex_solution%abso)
        call VecCreateSeq(PETSC_COMM_SELF, i0, rayli_info%plex_solution%abso, ierr); call CHKERR(ierr)
      endif ! is subcomm master
    end associate
    call VecDestroy(hhl, ierr); call CHKERR(ierr)

    ! Setup albedo scatter context
    call DMGetGlobalVector(solver%Csrfc_one%da, vec_albedo, ierr); call CHKERR(ierr)
    if(submyid.eq.0) then
      call VecGetSize(vec_albedo, Nalbedo, ierr); call CHKERR(ierr)
    else
      Nalbedo = i0
    endif
    call VecCreateSeq(PETSC_COMM_SELF, Nalbedo*i2, rayli_info%albedo, ierr); call CHKERR(ierr)
    allocate(is_data(Nalbedo*i2), source=(/(i/i2, i=i0,Nalbedo*i2-i1)/))
    call ISCreateGeneral(PETSC_COMM_SELF, Nalbedo*i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
    call ISCreateStride(PETSC_COMM_SELF, Nalbedo*i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)
    call gen_shared_scatter_ctx(vec_albedo, rayli_info%albedo, rayli_info%ctx_albedo, ierr, &
      & is_in, is_out); call CHKERR(ierr)
    call ISDestroy(is_out, ierr); call CHKERR(ierr)
    call ISDestroy(is_in, ierr); call CHKERR(ierr)
    deallocate(is_data)
    call DMRestoreGlobalVector(solver%Csrfc_one%da, vec_albedo, ierr); call CHKERR(ierr)


    ! Setup optprop_scatter_ctx
    call DMGetGlobalVector(solver%C_one_atm%da, vec_optprop, ierr); call CHKERR(ierr)
    if(submyid.eq.0) then
      call VecGetSize(vec_optprop, Noptprop, ierr); call CHKERR(ierr)
    else
      Noptprop = 0
    endif
    call VecCreateSeq(PETSC_COMM_SELF, Noptprop*i2, rayli_info%kabs, ierr); call CHKERR(ierr)
    call VecCreateSeq(PETSC_COMM_SELF, Noptprop*i2, rayli_info%ksca, ierr); call CHKERR(ierr)
    call VecCreateSeq(PETSC_COMM_SELF, Noptprop*i2, rayli_info%g   , ierr); call CHKERR(ierr)
    allocate(is_data(i0:Noptprop*i2-i1))
    if(submyid.eq.0) then
      do i = i0, solver%C_one_atm%glob_xm*solver%C_one_atm%glob_ym-i1
        do k = i0, solver%C_one_atm%glob_zm-i1
          is_data((i*i2  )*solver%C_one_atm%glob_zm+k) = i*solver%C_one_atm%glob_zm+k
          is_data((i*i2+1)*solver%C_one_atm%glob_zm+k) = i*solver%C_one_atm%glob_zm+k
        enddo
      enddo
    endif
    call ISCreateGeneral(PETSC_COMM_SELF, Noptprop*i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
    call ISCreateStride(PETSC_COMM_SELF, Noptprop*i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)
    call gen_shared_scatter_ctx(vec_optprop, rayli_info%kabs, rayli_info%ctx_optprop, ierr, &
      & is_in, is_out); call CHKERR(ierr)
    call ISDestroy(is_out, ierr); call CHKERR(ierr)
    call ISDestroy(is_in, ierr); call CHKERR(ierr)
    deallocate(is_data)
    call DMRestoreGlobalVector(solver%C_one_atm%da, vec_optprop, ierr); call CHKERR(ierr)

    ! Setup result scatter contexts
    if(submyid.eq.0) then
      Nedir = solver%C_dir%glob_xm*solver%C_dir%glob_ym*solver%C_dir%glob_zm
    else
      Nedir = 0
    endif
    allocate(is_data(i0:Nedir*i2-i1))
    if(submyid.eq.0) then
      do i = i0, solver%C_dir%glob_xm*solver%C_dir%glob_ym-i1
        do k = i0, solver%C_dir%glob_zm-i1
          kk = atmk(solver%atm,k)
          voff = i*solver%C_dir%glob_zm*solver%C_dir%dof + k*solver%C_dir%dof
          is_data((i*i2  )*solver%C_one_atm1%glob_zm+kk) = voff
          is_data((i*i2+1)*solver%C_one_atm1%glob_zm+kk) = voff
        enddo
      enddo
    endif
    call ISCreateGeneral(PETSC_COMM_SELF, Nedir*i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
    call ISCreateStride(PETSC_COMM_SELF, Nedir*i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)
    call gen_shared_scatter_ctx(solution%edir, rayli_info%plex_solution%edir, rayli_info%ctx_edir, ierr, &
      & is_in, is_out); call CHKERR(ierr)
    call PetscObjectViewFromOptions(is_in , PETSC_NULL_IS, '-show_rayli_iss', ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(is_out, PETSC_NULL_IS, '-show_rayli_iss', ierr); call CHKERR(ierr)
    call ISDestroy(is_out, ierr); call CHKERR(ierr)
    call ISDestroy(is_in, ierr); call CHKERR(ierr)
    deallocate(is_data)

    if(submyid.eq.0) then
      Nediff = solver%C_diff%glob_xm*solver%C_diff%glob_ym*solver%C_diff%glob_zm*i2
    else
      Nediff = 0
    endif
    allocate(is_data(i0:Nediff*i2-i1))
    if(submyid.eq.0) then
      do i = i0, solver%C_diff%glob_xm*solver%C_diff%glob_ym-i1
        do k = i0, solver%C_diff%glob_zm-i2
          kk = atmk(solver%atm,k)
          voff = i*solver%C_diff%glob_zm*solver%C_diff%dof + k*solver%C_diff%dof

          off_plex(1) = (i*i2   ) * solver%C_one_atm1%glob_zm * i2  + kk   * i2
          off_plex(2) = (i*i2+i1) * solver%C_one_atm1%glob_zm * i2  + kk   * i2

          do idof = i0, i1
            if(solver%difftop%is_inward(i1+idof)) then ! edn
              is_data(off_plex(1)) = voff+idof
              is_data(off_plex(2)) = voff+idof
            else
              is_data(off_plex(1)+i1) = voff+idof
              is_data(off_plex(2)+i1) = voff+idof
            endif
          enddo
        enddo
        ! at the surface reversed inout dof (see details in plex_grid)
        kk = atmk(solver%atm,k)
        voff = i*solver%C_diff%glob_zm*solver%C_diff%dof + k*solver%C_diff%dof

        off_plex(1) = (i*i2   ) * solver%C_one_atm1%glob_zm * i2  + kk   * i2
        off_plex(2) = (i*i2+i1) * solver%C_one_atm1%glob_zm * i2  + kk   * i2

        do idof = i0, i1
          if(solver%difftop%is_inward(i1+idof)) then ! edn
            is_data(off_plex(1)+i1) = voff+idof
            is_data(off_plex(2)+i1) = voff+idof
          else
            is_data(off_plex(1)) = voff+idof
            is_data(off_plex(2)) = voff+idof
          endif
        enddo
      enddo
    endif
    call ISCreateGeneral(PETSC_COMM_SELF, Nediff*i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
    call ISCreateStride(PETSC_COMM_SELF, Nediff*i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(is_in , PETSC_NULL_IS, '-show_rayli_iss', ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(is_out, PETSC_NULL_IS, '-show_rayli_iss', ierr); call CHKERR(ierr)
    call gen_shared_scatter_ctx(solution%ediff, rayli_info%plex_solution%ediff, rayli_info%ctx_ediff, ierr, &
      & is_in, is_out); call CHKERR(ierr)
    call ISDestroy(is_out, ierr); call CHKERR(ierr)
    call ISDestroy(is_in, ierr); call CHKERR(ierr)
    deallocate(is_data)

    if(submyid.eq.0) then
      Nabso = solver%C_one%glob_xm*solver%C_one%glob_ym*solver%C_one%glob_zm
    else
      Nabso = 0
    endif
    allocate(is_data(i0:Nabso*i2-i1))
    if(submyid.eq.0) then
      do i = i0, solver%C_one%glob_xm*solver%C_one%glob_ym-i1
        do k = i0, solver%C_one%glob_zm-i1
          kk = atmk(solver%atm,k)
          voff = i*solver%C_one%glob_zm + k
          is_data((i*i2  )*solver%C_one_atm%glob_zm+kk) = voff
          is_data((i*i2+1)*solver%C_one_atm%glob_zm+kk) = voff
        enddo
      enddo
    endif
    call ISCreateGeneral(PETSC_COMM_SELF, Nabso*i2, is_data, PETSC_USE_POINTER, is_in, ierr); call CHKERR(ierr)
    call ISCreateStride(PETSC_COMM_SELF, Nabso*i2, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)
    call gen_shared_scatter_ctx(solution%abso, rayli_info%plex_solution%abso, rayli_info%ctx_abso, ierr, &
      & is_in, is_out); call CHKERR(ierr)
    call PetscObjectViewFromOptions(is_in , PETSC_NULL_IS, '-show_rayli_iss', ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(is_out, PETSC_NULL_IS, '-show_rayli_iss', ierr); call CHKERR(ierr)
    call ISDestroy(is_out, ierr); call CHKERR(ierr)
    call ISDestroy(is_in, ierr); call CHKERR(ierr)
    deallocate(is_data)
  end subroutine

  !> @brief wrapper for the rayli montecarlo solver
  !> @details solve the radiative transfer equation with the wedge bindings to  rayli
    ! TODO: tasks:
    ! * copy pprts dm to Zero (and all shared mem root-ranks)
    ! * implement pprts to wedge interface
    ! * average results over shared mem comm
    ! * distribute results
  subroutine pprts_rayli_wrapper(lcall_solver, lcall_snap, solver, edirTOA, solution)
    logical, intent(in) :: lcall_solver, lcall_snap
    class(t_solver), intent(inout)        :: solver
    real(ireals),intent(in)               :: edirTOA
    type(t_state_container),intent(inout) :: solution

    integer(mpiint) :: ierr
    integer(mpiint) :: myid, numnodes
    integer(mpiint) :: submyid, subnumnodes

    real(ireals) :: sundir(3)

    if(all([lcall_solver, lcall_snap].eqv..False.)) return

    if(solution%lsolar_rad.eqv..False.) &
      & call CHKERR(1_mpiint, "Cannot use Rayli for thermal computations")

    sundir = spherical_2_cartesian(meanval(solver%sun%phi), meanval(solver%sun%theta)) &
      & * edirTOA

    call init_pprts_rayli_wrapper(solver, solution, rayli_info)

    call mpi_comm_rank(solver%comm,myid,ierr); call CHKERR(ierr)
    call mpi_comm_size(solver%comm,numnodes,ierr); call CHKERR(ierr)
    call mpi_comm_rank(rayli_info%subcomm,submyid,ierr); call CHKERR(ierr)
    call mpi_comm_size(rayli_info%subcomm,subnumnodes,ierr); call CHKERR(ierr)

    call prepare_input()
    call call_solver(rayli_info%plex_solution)
    call transfer_result(solution)

    contains
      subroutine prepare_input()
        type(tVec) :: glob_vec, nat_vec
        character(len=*), parameter :: log_event_name="pprts_rayli_prepare_input"
        PetscClassId :: cid
        PetscLogEvent :: log_event

        call PetscClassIdRegister("pprts_rayli", cid, ierr); call CHKERR(ierr)
        call PetscLogEventRegister(log_event_name, cid, log_event, ierr); call CHKERR(ierr)
        call PetscLogEventBegin(log_event, ierr); call CHKERR(ierr)

        associate( &
            & atm => solver%atm,            &
            & Cs  => solver%Csrfc_one,      &
            & Ca1 => solver%C_one_atm1,     &
            & Ca  => solver%C_one_atm,      &
            & Cv  => solver%Cvert_one_atm1, &
            & ri  => rayli_info             )

          call DMDACreateNaturalVector(Cs%da, nat_vec, ierr); call CHKERR(ierr)
          call DMGetGlobalVector(Cs%da, glob_vec, ierr); call CHKERR(ierr)
          call f90VecToPetsc(atm%albedo, Cs%da, glob_vec)
          call DMDAGlobalToNaturalBegin(Cs%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call DMDAGlobalToNaturalEnd  (Cs%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call VecScatterBegin(ri%ctx_albedo, nat_vec, ri%albedo, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call VecScatterEnd  (ri%ctx_albedo, nat_vec, ri%albedo, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call DMRestoreGlobalVector(Cs%da, glob_vec, ierr); call CHKERR(ierr)
          call VecDestroy(nat_vec, ierr); call CHKERR(ierr)


          call DMDACreateNaturalVector(Ca%da, nat_vec, ierr); call CHKERR(ierr)
          call DMGetGlobalVector(Ca%da, glob_vec, ierr); call CHKERR(ierr)

          call f90VecToPetsc(atm%kabs, Ca%da, glob_vec)
          call DMDAGlobalToNaturalBegin(Ca%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call DMDAGlobalToNaturalEnd  (Ca%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call VecScatterBegin(ri%ctx_optprop, nat_vec, ri%kabs, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call VecScatterEnd  (ri%ctx_optprop, nat_vec, ri%kabs, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

          call f90VecToPetsc(atm%ksca, Ca%da, glob_vec)
          call DMDAGlobalToNaturalBegin(Ca%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call DMDAGlobalToNaturalEnd  (Ca%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call VecScatterBegin(ri%ctx_optprop, nat_vec, ri%ksca, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call VecScatterEnd  (ri%ctx_optprop, nat_vec, ri%ksca, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

          call f90VecToPetsc(atm%g, Ca%da, glob_vec)
          call DMDAGlobalToNaturalBegin(Ca%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call DMDAGlobalToNaturalEnd  (Ca%da, glob_vec, INSERT_VALUES, nat_vec, ierr); call CHKERR(ierr)
          call VecScatterBegin(ri%ctx_optprop, nat_vec, ri%g, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
          call VecScatterEnd  (ri%ctx_optprop, nat_vec, ri%g, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

          call DMRestoreGlobalVector(Ca%da, glob_vec, ierr); call CHKERR(ierr)
          call VecDestroy(nat_vec, ierr); call CHKERR(ierr)

          call PetscObjectViewFromOptions(rayli_info%kabs  , PETSC_NULL_VEC, '-show_rayli_kabs', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(rayli_info%ksca  , PETSC_NULL_VEC, '-show_rayli_ksca', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(rayli_info%g     , PETSC_NULL_VEC, '-show_rayli_g', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(rayli_info%albedo, PETSC_NULL_VEC, '-show_rayli_albedo', ierr); call CHKERR(ierr)
        end associate
        call PetscLogEventEnd(log_event, ierr); call CHKERR(ierr)
      end subroutine

      subroutine call_solver(plex_solution)
        type(t_state_container),intent(inout) :: plex_solution
        logical :: gotmsg
        integer(iintegers) :: idummy
        integer(mpiint) :: isub, status(MPI_STATUS_SIZE), mpi_request
        integer(mpiint), parameter :: FINALIZEMSG=1
        integer(mpiint) :: run_rank=0

        character(len=*), parameter :: log_event_name="pprts_rayli_call_solver"
        PetscClassId :: cid
        PetscLogEvent :: log_event

        call PetscClassIdRegister("pprts_rayli", cid, ierr); call CHKERR(ierr)
        call PetscLogEventRegister(log_event_name, cid, log_event, ierr); call CHKERR(ierr)
        call PetscLogEventBegin(log_event, ierr); call CHKERR(ierr)

        if(submyid.eq.run_rank) then
          plex_solution%uid = solution%uid
          call rayli_wrapper(lcall_solver, lcall_snap, &
            rayli_info%plex, rayli_info%kabs, rayli_info%ksca, rayli_info%g, rayli_info%albedo, &
            & sundir, plex_solution, petsc_log=solver%logs%rayli_tracing)

          call PetscObjectViewFromOptions(plex_solution%edir , PETSC_NULL_VEC, '-show_plex_rayli_edir', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(plex_solution%ediff, PETSC_NULL_VEC, '-show_plex_rayli_ediff', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(plex_solution%abso , PETSC_NULL_VEC, '-show_plex_rayli_abso', ierr); call CHKERR(ierr)

          do isub=0,subnumnodes-1 ! send finalize msg to all others to stop waiting
            if(isub.ne.run_rank) then
              call mpi_isend(-i1, 1_mpiint, imp_iinteger, isub, FINALIZEMSG, &
                & rayli_info%subcomm, mpi_request, ierr); call CHKERR(ierr)
            endif
          enddo
        else
          lazy_wait: do ! prevent eager MPI polling while the one rank performs rayli computations
            call mpi_iprobe(MPI_ANY_SOURCE, FINALIZEMSG, rayli_info%subcomm, gotmsg, status, ierr); call CHKERR(ierr)
            if(gotmsg) then
              call mpi_recv(idummy, 1_mpiint, imp_iinteger, status(MPI_SOURCE), FINALIZEMSG, &
                & rayli_info%subcomm, status, ierr); call CHKERR(ierr)
              exit lazy_wait
            endif
            call sleep(1)
          enddo lazy_wait
        endif
        call PetscLogEventEnd(log_event, ierr); call CHKERR(ierr)
      end subroutine

      subroutine transfer_result(solution)
        type(t_state_container),intent(inout) :: solution

        real(ireals) :: fac
        type(tVec) :: nat_vec

        character(len=*), parameter :: log_event_name="pprts_rayli_transfer_result"
        PetscClassId :: cid
        PetscLogEvent :: log_event

        call PetscClassIdRegister("pprts_rayli", cid, ierr); call CHKERR(ierr)
        call PetscLogEventRegister(log_event_name, cid, log_event, ierr); call CHKERR(ierr)
        call PetscLogEventBegin(log_event, ierr); call CHKERR(ierr)

        associate( &
            & psol => rayli_info%plex_solution, &
            & ri   => rayli_info                )
          fac = one / 2._ireals / real(ri%num_subcomm_masters, ireals)

          if(allocated(solution%edir)) then
            call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
          endif
          call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
          call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)

          if(allocated(solution%edir)) then
            call DMDACreateNaturalVector(solver%C_dir%da, nat_vec, ierr); call CHKERR(ierr)
            call VecScatterBegin(ri%ctx_edir, psol%edir, nat_vec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
            call VecScatterEnd  (ri%ctx_edir, psol%edir, nat_vec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
            call DMDANaturalToGlobalBegin(solver%C_dir%da, nat_vec, INSERT_VALUES, solution%edir, ierr); call CHKERR(ierr)
            call DMDANaturalToGlobalEnd  (solver%C_dir%da, nat_vec, INSERT_VALUES, solution%edir, ierr); call CHKERR(ierr)
            call VecDestroy(nat_vec, ierr); call CHKERR(ierr)
          endif

          call DMDACreateNaturalVector(solver%C_diff%da, nat_vec, ierr); call CHKERR(ierr)
          call VecScatterBegin(ri%ctx_ediff, psol%ediff, nat_vec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
          call VecScatterEnd  (ri%ctx_ediff, psol%ediff, nat_vec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
          call DMDANaturalToGlobalBegin(solver%C_diff%da, nat_vec, INSERT_VALUES, solution%ediff, ierr); call CHKERR(ierr)
          call DMDANaturalToGlobalEnd  (solver%C_diff%da, nat_vec, INSERT_VALUES, solution%ediff, ierr); call CHKERR(ierr)
          call VecDestroy(nat_vec, ierr); call CHKERR(ierr)

          call DMDACreateNaturalVector(solver%C_one%da, nat_vec, ierr); call CHKERR(ierr)
          call VecScatterBegin(ri%ctx_abso , psol%abso , nat_vec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
          call VecScatterEnd  (ri%ctx_abso , psol%abso , nat_vec, ADD_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
          call DMDANaturalToGlobalBegin(solver%C_one%da, nat_vec, INSERT_VALUES, solution%abso, ierr); call CHKERR(ierr)
          call DMDANaturalToGlobalEnd  (solver%C_one%da, nat_vec, INSERT_VALUES, solution%abso, ierr); call CHKERR(ierr)
          call VecDestroy(nat_vec, ierr); call CHKERR(ierr)
        end associate

        if(allocated(solution%edir)) then
          call VecScale(solution%edir, fac, ierr); call CHKERR(ierr)
        endif
        call VecScale(solution%ediff, fac, ierr); call CHKERR(ierr)
        call VecScale(solution%abso, fac, ierr); call CHKERR(ierr)

        if(allocated(solution%edir)) then
          call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, '-show_rayli_edir', ierr); call CHKERR(ierr)
        endif
        call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, '-show_rayli_ediff', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-show_rayli_abso', ierr); call CHKERR(ierr)

        !Rayli solver returns fluxes as [W]
        solution%lWm2_dir  = .True.
        solution%lWm2_diff = .True.
        ! and mark solution that it is up to date (to prevent absoprtion computations)
        solution%lchanged  = .False.
        call PetscLogEventEnd(log_event, ierr); call CHKERR(ierr)
      end subroutine
  end subroutine

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine twostream(solver, edirTOA, solution)
    class(t_solver), intent(inout) :: solver
    real(ireals),intent(in)       :: edirTOA
    type(t_state_container)       :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: xv_dir=>null(),xv_diff=>null()
    real(ireals),pointer,dimension(:) :: xv_dir1d=>null(),xv_diff1d=>null()
    integer(iintegers) :: i,j,src

    real(ireals),allocatable :: dtau(:),kext(:),w0(:),g(:),S(:),Edn(:),Eup(:)
    real(ireals) :: mu0,incSolar,fac
    integer(mpiint) :: ierr

    associate(atm         => solver%atm, &
        C_diff      => solver%C_diff, &
        C_dir       => solver%C_dir, &
        C_one_atm   => solver%C_one_atm, &
        C_one_atm1  => solver%C_one_atm1)

      if(solution%lsolar_rad) then
        call PetscObjectSetName(solution%edir,'twostream_edir_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
        call VecSet(solution%edir ,zero,ierr); call CHKERR(ierr)
      endif

      call PetscObjectSetName(solution%ediff,'twostream_ediff_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
      call VecSet(solution%ediff,zero,ierr); call CHKERR(ierr)

      allocate( dtau(C_one_atm%zm) )
      allocate( kext(C_one_atm%zm) )
      allocate(   w0(C_one_atm%zm) )
      allocate(    g(C_one_atm%zm) )

      if(solution%lsolar_rad) &
        call getVecPointer(solution%edir  ,C_dir%da  ,xv_dir1d , xv_dir)
      call getVecPointer(solution%ediff ,C_diff%da ,xv_diff1d, xv_diff)

      allocate( S  (C_one_atm1%zs:C_one_atm1%ze) )
      allocate( Eup(C_one_atm1%zs:C_one_atm1%ze) )
      allocate( Edn(C_one_atm1%zs:C_one_atm1%ze) )

      do j=C_one_atm%ys,C_one_atm%ye
        do i=C_one_atm%xs,C_one_atm%xe

          mu0 = solver%sun%costheta(C_one_atm1%zs,i,j)
          incSolar = edirTOA* mu0

          kext = atm%kabs(:,i,j) + atm%ksca(:,i,j)
          dtau = atm%dz(:,i,j)* kext
          w0   = atm%ksca(:,i,j) / max(kext, epsilon(kext))
          g    = atm%g(:,i,j)

          if(allocated(atm%planck) ) then
            if(allocated(atm%Bsrfc)) then
              call delta_eddington_twostream(dtau, w0, g, &
                mu0, incSolar, atm%albedo(i,j), &
                S, Edn, Eup, &
                planck=atm%planck(:,i,j), &
                planck_srfc=atm%Bsrfc(i,j))
            else
              call delta_eddington_twostream(dtau, w0, g, &
                mu0, incSolar, atm%albedo(i,j), &
                S, Edn, Eup, &
                planck=atm%planck(:,i,j) )
            endif
          else
            !
            ! call adding_delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup )
            !
            !TODO investigate if this one is really ok...
            ! I recently had valgrind errors in VecNorm after calling this:
            ! make -j ex_pprts_rrtm_lw_sw &&
            ! mpirun -np 1 -wdir ../examples/rrtm_lw_sw/ valgrind $(pwd)/bin/ex_pprts_rrtm_lw_sw -twostr_only
            call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup )
          endif

          if(solution%lsolar_rad) then
            fac = real(solver%dirtop%area_divider, ireals) / real(solver%dirtop%streams, ireals)
            do src=i0,solver%dirtop%dof-1
              xv_dir(src,C_dir%zs+1:C_dir%ze,i,j) = S(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
              xv_dir(src,C_dir%zs           ,i,j) = S(C_one_atm1%zs) * fac
            enddo
          endif

          fac = real(solver%difftop%area_divider, ireals) / real(solver%difftop%streams, ireals)
          do src = 1, solver%difftop%dof
            if(solver%difftop%is_inward(src)) then
              xv_diff(src-1,C_diff%zs+1:C_diff%ze,i,j) = Edn(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
              xv_diff(src-1,C_diff%zs            ,i,j) = Edn(C_one_atm1%zs) * fac
            else
              xv_diff(src-1,C_diff%zs+1:C_diff%ze,i,j) = Eup(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
              xv_diff(src-1,C_diff%zs            ,i,j) = Eup(C_one_atm1%zs) * fac
            endif
          enddo
        enddo
      enddo

      if(solution%lsolar_rad) &
        call restoreVecPointer(solution%edir, xv_dir1d, xv_dir  )
      call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

      !Twostream solver returns fluxes as [W]
      solution%lWm2_dir  = .True.
      solution%lWm2_diff = .True.
      ! and mark solution that it is not up to date
      solution%lchanged  = .True.

      deallocate(S)
      deallocate(Edn)
      deallocate(Eup)

    end associate
  end subroutine

  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine schwarz(solver, solution)
    class(t_solver)         :: solver
    type(t_state_container) :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: xv_diff=>null()
    real(ireals),pointer,dimension(:)       :: xv_diff1d=>null()
    integer(iintegers) :: i,j,k,idof
    integer(iintegers) :: Nmu, ak
    logical :: lflg

    real(ireals),allocatable :: dtau(:),Edn(:),Eup(:)
    integer(mpiint) :: ierr


    associate( &
        C_diff => solver%C_diff, &
        atm    => solver%atm,    &
        C_one  => solver%C_one,  &
        C_one1 => solver%C_one1)

      if(solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling schwarschild solver for solar calculation -- stopping!')
      if(.not.allocated(atm%planck)) call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

      allocate(dtau(size(atm%dz,dim=1)))

      call getVecPointer(solution%ediff, C_diff%da, xv_diff1d, xv_diff)

      allocate( Eup(0:size(atm%dz,dim=1)) )
      allocate( Edn(0:size(atm%dz,dim=1)) )

      Nmu = 10
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
        "-schwarzschild_Nmu" , Nmu, lflg , ierr) ;call CHKERR(ierr)

      if(solver%myid.eq.0 .and. ldebug) print *,' CALCULATING schwarzschild ::'

      do j = C_diff%ys, C_diff%ye
        do i = C_diff%xs, C_diff%xe

          dtau = atm%dz(:, i, j) * atm%kabs(:, i, j)

          if(allocated(atm%Bsrfc)) then
            call schwarzschild(Nmu, dtau, atm%albedo(i,j), Edn, Eup, &
              atm%planck(:, i, j), opt_srfc_emission=atm%Bsrfc(i,j))
          else
            call schwarzschild(Nmu, dtau, atm%albedo(i,j), Edn, Eup, &
              atm%planck(:, i, j))
          endif

          ! icollapse needs special case for TOA flx's
          do idof = 0, solver%difftop%dof-1
            if (solver%difftop%is_inward(i1+idof)) then ! Edn
              xv_diff(idof,C_diff%zs,i,j) = Edn(0)
            else ! Eup
              xv_diff(idof,C_diff%zs,i,j) = Eup(0)
            endif
          enddo

          ! rest of the atmosphere
          do k=C_diff%zs+1,C_diff%ze
            ak = atmk(atm,k)
            do idof = 0, solver%difftop%dof-1
              if (solver%difftop%is_inward(i1+idof)) then ! Edn
                xv_diff(idof,k,i,j) = Edn(ak)
              else ! Eup
                xv_diff(idof,k,i,j) = Eup(ak)
              endif
            enddo
          enddo
        enddo
      enddo
      xv_diff = xv_diff / real(solver%difftop%streams, ireals)

      call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

      !Schwarzschild solver returns fluxes as [W/m^2]
      solution%lWm2_dir  = .True.
      solution%lWm2_diff = .True.
      ! and mark solution that it is not up to date
      solution%lchanged         = .True.

      deallocate(Edn)
      deallocate(Eup)

    end associate
  end subroutine
end module
