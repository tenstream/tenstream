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

module m_mcdmda
  use iso_fortran_env, only: int64, output_unit
  use iso_c_binding, only: c_backspace

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, ireal_dp, &
                               mpiint, imp_iinteger, imp_int8, &
                               zero, i0, i1, i2, i3, i4, i5, i6, pi

  use m_helper_functions, only: &
    & CHKERR, &
    & CHKWARN, &
    & cstr, &
    & deg2rad, &
    & expm1, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_allreduce_sum, &
    & ind_1d_to_nd, &
    & ind_nd_to_1d, &
    & ndarray_offsets, &
    & rotate_angle_x, &
    & rotate_angle_y, &
    & spherical_2_cartesian, &
    & toStr

  use m_boxmc, only: t_photon, print_photon, scatter_photon, roulette, R, &
                     tau, distance, update_photon_loc, absorb_photon, &
                     t_boxmc, t_boxmc_1_2, t_boxmc_3_6, &
                     imp_t_photon

  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use m_intersection, only: pnt_in_cube

  use m_pprts_base, only: &
    & atmk, &
    & t_atmosphere, &
    & t_coord, &
    & t_solver, &
    & t_solver_1_2, &
    & t_solver_3_6, &
    & t_solver_mcdmda, &
    & t_state_container

  use m_petsc_helpers, only: getVecPointer, restoreVecPointer

  use m_buildings, only: t_pprts_buildings

  use m_boxmc_geometry, only: &
    & PPRTS_TOP_FACE, &
    & PPRTS_BOT_FACE, &
    & PPRTS_LEFT_FACE, &
    & PPRTS_RIGHT_FACE, &
    & PPRTS_FRONT_FACE, &
    & PPRTS_REAR_FACE

  use m_linked_list_int64, only: t_list_int64, t_node

  implicit none

  ! queue status indices
  integer(iintegers), parameter :: &
    PQ_SELF = 1, &
    PQ_NORTH = 2, &
    PQ_EAST = 3, &
    PQ_SOUTH = 4, &
    PQ_WEST = 5

  character(len=16), parameter :: id2name(5) = [ &
    & character(len=16) :: &
    & 'PQ_SELF', &
    & 'PQ_NORTH', &
    & 'PQ_EAST', &
    & 'PQ_SOUTH', &
    & 'PQ_WEST' &
    & ]

  type :: t_distributed_photon
    type(t_photon) :: p
    integer(mpiint) :: request
  end type

  type :: t_photon_queue
    type(t_distributed_photon), allocatable :: photons(:)
    type(t_list_int64) :: ready ! linked list for read_to_go photon indices
    type(t_list_int64) :: empty ! linked_list of empty slots in this queue
    type(t_list_int64) :: sending ! linked_list of sending slots in this queue
    integer(mpiint) :: owner ! owner is the owning rank, i.e. myid or the neighbor id
    integer(iintegers) :: queue_index ! is the STATUS integer, i.e. one of PQ_SELF, PQ_NORTH etc.
  end type

  logical, parameter :: ldebug = .false.
  logical, parameter :: ldebug_tracing = .false.
  logical :: lfull3d = .false. ! option if 1d layers are treated 1d or really 3d (heavy toll on performance for high resolution)

  real(ireal_dp), parameter :: loceps = 0 !sqrt(epsilon(loceps))

  real(ireals), parameter :: blocking_waittime = 5 ! sec
contains

  subroutine solve_mcdmda(solver, edirTOA, solution, ierr, opt_buildings)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container) :: solution
    integer(mpiint), intent(out) :: ierr
    type(t_pprts_buildings), intent(in), optional :: opt_buildings

    integer(mpiint) :: myid

    type(t_photon_queue) :: pqueues(5) ! [own, north, east, south, west]

    integer(int64) :: Nqueuesize, Nbatchsize
    real(ireals) :: rNqueuesize, rNbatchsize
    integer(int64) :: Nphotons_global
    integer(int64) :: locally_started_photons, globally_started_photons, globally_killed_photons
    integer(int64) :: Nphotons_local
    integer(int64) :: started_photons
    integer(int64) :: killed_photons
    integer(int64) :: ip, kp
    integer(iintegers) :: iter, percent_printed, last_percent_printed
    integer(mpiint) :: started_request, killed_request, stat(mpi_status_size)
    logical :: lcomm_finished, lflg, lfirst_print, lfinish_border_photons_first
    real(ireals) :: photon_weight

    class(t_boxmc), allocatable :: bmc
    type(t_boxmc_1_2) :: bmc_1_2 ! 1_2 tracer needed for 1D layers

    real(ireal_dp), dimension(:, :, :, :), allocatable :: edir, ediff, abso
    integer(iintegers), dimension(:, :, :, :), allocatable :: Nediff, buildings_idx

    ierr = 0

    call determine_Nphotons(solver, Nphotons_local, ierr); call CHKERR(ierr)
    call mpi_allreduce(Nphotons_local, Nphotons_global, 1_mpiint, imp_int8, &
                       MPI_SUM, solver%comm, ierr); call CHKERR(ierr)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    call get_petsc_opt(solver%prefix, "-mcdmda_full3d", lfull3d, lflg, ierr); call CHKERR(ierr)

    Nbatchsize = 10000
    call get_petsc_opt(solver%prefix, "-mcdmda_batch_size", rNbatchsize, lflg, ierr); call CHKERR(ierr)
    if (lflg) Nbatchsize = int(rNbatchsize, kind(Nbatchsize))

    Nqueuesize = int(Nphotons_local, int64)
    call get_petsc_opt(solver%prefix, "-mcdmda_queue_size", rNqueuesize, lflg, ierr); call CHKERR(ierr)
    if (lflg) Nqueuesize = int(rNqueuesize, kind(Nqueuesize))

    lfinish_border_photons_first = .false.
    call get_petsc_opt(solver%prefix, "-mcdmda_finish_border_photons_first", lfinish_border_photons_first, lflg, ierr)
    call CHKERR(ierr)

    associate (C => solver%C_one_atm, C1 => solver%C_one_atm1, Cdir => solver%C_dir, Cdiff => solver%C_diff)
      if (ldebug) then
        print *, myid, 'Edir TOA', edirTOA, ': Nphotons_local', Nphotons_local, 'Nphotons_global', Nphotons_global
        print *, myid, 'Domain start:', C%zs, C%xs, C%ys
        print *, myid, 'Domain end  :', C%ze, C%xe, C%ye
        print *, myid, 'Domain Size :', C%zm, C%xm, C%ym
        print *, myid, 'Global Size :', C%glob_zm, C%glob_xm, C%glob_ym

        print *, myid, 'my neighs NESW', &
          C%neighbors(22), C%neighbors(16), C%neighbors(4), C%neighbors(10)
      end if

      if (solution%lsolar_rad) then
        allocate (edir(0:Cdir%dof - 1, C1%zs:C1%ze, C%gxs:C%gxe, C%gys:C%gye), source=0._ireal_dp)
        photon_weight = edirTOA * real(solver%sun%costheta, ireals) &
                       & * solver%C_one_atm%xm * solver%C_one_atm%ym / real(Nphotons_local, ireals)
      else
        photon_weight = 0
      end if

      allocate ( &
        & ediff(0:Cdiff%dof - 1, C1%zs:C1%ze, C%gxs:C%gxe, C%gys:C%gye), &
        & abso(0:C%dof - 1, C%zs:C%ze, C%xs:C%xe, C%ys:C%ye), &
        & source=0._ireal_dp)
      allocate ( &
        & Nediff(0:Cdiff%dof - 1, C1%zs:C1%ze, C%gxs:C%gxe, C%gys:C%gye), &
        & source=0_iintegers)

      call fill_buildings_idx(solver, opt_buildings, buildings_idx)

      call setup_photon_queue(pqueues(PQ_SELF), Nphotons_local, myid, PQ_SELF)
      call setup_photon_queue(pqueues(PQ_NORTH), Nqueuesize, C%neighbors(22), PQ_NORTH)
      call setup_photon_queue(pqueues(PQ_EAST), Nqueuesize, C%neighbors(16), PQ_EAST)
      call setup_photon_queue(pqueues(PQ_SOUTH), Nqueuesize, C%neighbors(4), PQ_SOUTH)
      call setup_photon_queue(pqueues(PQ_WEST), Nqueuesize, C%neighbors(10), PQ_WEST)

      if (C%neighbors(22) .lt. zero) call CHKERR(C%neighbors(22), 'Bad neighbor id for PQ_NORTH')
      if (C%neighbors(16) .lt. zero) call CHKERR(C%neighbors(16), 'Bad neighbor id for PQ_EAST')
      if (C%neighbors(4) .lt. zero) call CHKERR(C%neighbors(4), 'Bad neighbor id for PQ_SOUTH')
      if (C%neighbors(10) .lt. zero) call CHKERR(C%neighbors(10), 'Bad neighbor id for PQ_WEST')
    end associate
    select type (solver)
    class is (t_solver_1_2)
      allocate (t_boxmc_1_2 :: bmc)

    class is (t_solver_3_6)
      allocate (t_boxmc_3_6 :: bmc)

    class is (t_solver_mcdmda)
      allocate (t_boxmc_3_6 :: bmc)

    class default
      call CHKERR(1_mpiint, 'initialize bmc for mcdmda: unexpected type for solver object for DMDA computations!'// &
                  '-- call with -solver 1_2 or -solver 3_6')

    end select

    if (ldebug_tracing) then
      call bmc%init(MPI_COMM_SELF, rngseed=9, luse_random_seed=.false.)
      call bmc_1_2%init(MPI_COMM_SELF, rngseed=9, luse_random_seed=.false.)
    else
      call bmc%init(MPI_COMM_SELF, luse_random_seed=.true.)
      call bmc_1_2%init(MPI_COMM_SELF, luse_random_seed=.true.)
    end if

    ! Initialize the locally owned photons
    call prepare_locally_owned_photons(solver, bmc, solution%lsolar_rad, pqueues(PQ_SELF), Nphotons_local, weight=photon_weight)

    killed_photons = 0
    call mpi_iallreduce(killed_photons, globally_killed_photons, 1_mpiint, imp_int8, &
                        MPI_SUM, solver%comm, killed_request, ierr); call CHKERR(ierr)

    if (ldebug) then
      lfirst_print = .true.
      locally_started_photons = 0
      call mpi_iallreduce(locally_started_photons, globally_started_photons, 1_mpiint, imp_int8, &
                          MPI_SUM, solver%comm, started_request, ierr); call CHKERR(ierr)
    end if

    last_percent_printed = -1
    iter = 0
    do
      iter = iter + 1

      call run_photon_queue( &
        & solver, bmc, bmc_1_2, &
        & pqueues, PQ_SELF, &
        & edir, ediff, Nediff, abso, &
        & started_photons=ip, &
        & killed_photons=kp, &
        & limit_number_photons=Nbatchsize, &
        & opt_buildings=opt_buildings, &
        & buildings_idx=buildings_idx)

      locally_started_photons = locally_started_photons + ip

      started_photons = ip; killed_photons = killed_photons + kp
      if (ldebug) print *, 'SELF', 'started_photons', started_photons, 'killed_photons', killed_photons

      remote_photons: do ! Run remote photons until they are done, only then go on, this hopefully keeps the queue small
        started_photons = 0
        call exchange_photons(solver, pqueues)

        call run_photon_queue(&
          & solver, bmc, bmc_1_2, &
          & pqueues, PQ_NORTH, &
          & edir, ediff, Nediff, abso, &
          & started_photons=ip, &
          & killed_photons=kp, &
          & opt_buildings=opt_buildings, &
          & buildings_idx=buildings_idx)

        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'NORTH', 'started_photons', ip, 'killed_photons', kp

        call exchange_photons(solver, pqueues)

        call run_photon_queue( &
          & solver, bmc, bmc_1_2, &
          & pqueues, PQ_EAST, &
          & edir, ediff, Nediff, abso, &
          & started_photons=ip, &
          & killed_photons=kp, &
          & opt_buildings=opt_buildings, &
          & buildings_idx=buildings_idx)

        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'EAST', 'started_photons', ip, 'killed_photons', kp

        call exchange_photons(solver, pqueues)

        call run_photon_queue( &
          & solver, bmc, bmc_1_2, &
          & pqueues, PQ_SOUTH, &
          & edir, ediff, Nediff, abso, &
          & started_photons=ip, &
          & killed_photons=kp, &
          & opt_buildings=opt_buildings, &
          & buildings_idx=buildings_idx)

        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'SOUTH', 'started_photons', ip, 'killed_photons', kp

        call exchange_photons(solver, pqueues)

        call run_photon_queue( &
          & solver, bmc, bmc_1_2, &
          & pqueues, PQ_WEST, &
          & edir, ediff, Nediff, abso, &
          & started_photons=ip, &
          & killed_photons=kp, &
          & opt_buildings=opt_buildings, &
          & buildings_idx=buildings_idx)
        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'WEST', 'started_photons', ip, 'killed_photons', kp

        call exchange_photons(solver, pqueues)

        if (.not. lfinish_border_photons_first) then
          exit remote_photons
        elseif (started_photons .eq. 0) then
          exit remote_photons
        end if
      end do remote_photons

      if (ldebug) then
        if (lfirst_print .or. globally_started_photons .ne. Nphotons_global) then
          call mpi_test(started_request, lcomm_finished, stat, ierr); call CHKERR(ierr)
          if (lcomm_finished) then
            if (globally_started_photons .eq. Nphotons_global) lfirst_print = .false.
            if (myid .eq. 0) print *, '#it ', iter, &
              'Globally started photons', globally_started_photons, '/', Nphotons_global, &
              '('//toStr(100 * real(globally_started_photons) / real(Nphotons_global))//' % )'
            call mpi_iallreduce(locally_started_photons, globally_started_photons, 1_mpiint, imp_int8, &
                                MPI_SUM, solver%comm, started_request, ierr); call CHKERR(ierr)
          end if
        end if
      end if

      call mpi_test(killed_request, lcomm_finished, stat, ierr); call CHKERR(ierr)
      if (lcomm_finished) then
        if (myid .eq. 0) then
          percent_printed = int(100 * real(globally_killed_photons) / real(Nphotons_global), kind(percent_printed))
          if (ldebug .or. percent_printed .ne. last_percent_printed) then
            !print *, iter, 'Globally killed photons', globally_killed_photons, '/', Nphotons_global, &
            !  '('//toStr(percent_printed)//' % )'
            write (*, fmt="(A)", advance='no') repeat(c_backspace, 100)//'#it '//toStr(iter)// &
              & ' Globally killed photons '//toStr(globally_killed_photons)//' / '//toStr(Nphotons_global)// &
              & '( '//toStr(percent_printed)//' % )'
            call flush (output_unit)
            last_percent_printed = percent_printed
          end if
        end if
        if (ldebug) call debug_output_queues()
        if (globally_killed_photons .eq. Nphotons_global) then
          if (myid .eq. 0) print *, ''
          exit
        end if
        ! if we reach here this means there is still work todo, setup a new allreduce
        call mpi_iallreduce(killed_photons, globally_killed_photons, 1_mpiint, imp_int8, &
                            MPI_SUM, solver%comm, killed_request, ierr); call CHKERR(ierr)
      end if

    end do

    call get_result()

    ! cleanup
    call bmc%destroy(ierr); call CHKERR(ierr)
    call bmc_1_2%destroy(ierr); call CHKERR(ierr)
    do ip = 1, size(pqueues)
      call photon_queue_destroy(pqueues(ip), ierr); call CHKERR(ierr)
    end do
  contains

    subroutine fill_buildings_idx(solver, opt_buildings, buildings_idx)
      class(t_solver), intent(in) :: solver
      type(t_pprts_buildings), intent(in), optional :: opt_buildings
      integer(iintegers), dimension(:, :, :, :), allocatable :: buildings_idx
      integer(iintegers) :: m, idx(4)

      if (present(opt_buildings)) then
        associate (C => solver%C_one_atm, C_diff => solver%C_diff, atm => solver%atm)
          allocate (buildings_idx(6, C%zs:C%ze, C%xs:C%xe, C%ys:C%ye))
          buildings_idx = -1
          associate (B => opt_buildings)
            do m = 1, size(B%iface)
              call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
              idx(2:4) = idx(2:4) - 1 + [C_diff%zs, C_diff%xs, C_diff%ys]
              associate (k => idx(2), i => idx(3), j => idx(4))
                select case (idx(1))
                case (PPRTS_TOP_FACE)
                  buildings_idx(PPRTS_TOP_FACE, atmk(atm, k), i, j) = m

                case (PPRTS_BOT_FACE)
                  buildings_idx(PPRTS_BOT_FACE, atmk(atm, k), i, j) = m

                case (PPRTS_LEFT_FACE)
                  buildings_idx(PPRTS_LEFT_FACE, atmk(atm, k), i, j) = m

                case (PPRTS_RIGHT_FACE)
                  buildings_idx(PPRTS_RIGHT_FACE, atmk(atm, k), i, j) = m

                case (PPRTS_REAR_FACE)
                  buildings_idx(PPRTS_REAR_FACE, atmk(atm, k), i, j) = m

                case (PPRTS_FRONT_FACE)
                  buildings_idx(PPRTS_FRONT_FACE, atmk(atm, k), i, j) = m

                end select
              end associate
            end do
          end associate
        end associate
      end if
    end subroutine

    subroutine get_result()
      type(tVec) :: ediff_local, edir_local
      real(ireals), pointer, dimension(:, :, :, :) :: xv_dir => null(), xv_diff => null(), xv_abso => null()
      real(ireals), pointer, dimension(:) :: xv_dir1d => null(), xv_diff1d => null(), xv_abso1d => null()
      integer(iintegers) :: Nundersampled

      associate (&
          & atm => solver%atm,&
          & C_one => solver%C_one, &
          & C_dir => solver%C_dir, &
          & C_diff => solver%C_diff, &
          & C_one_atm => solver%C_one_atm,&
          & C_one_atm1 => solver%C_one_atm1)

        if (solution%lsolar_rad) then
          ! handle edir
          call DMGetLocalVector(C_dir%da, edir_local, ierr); call CHKERR(ierr)
          call VecSet(edir_local, 0._ireals, ierr); call CHKERR(ierr)

          call getVecPointer(C_dir%da, edir_local, xv_dir1d, xv_dir)

          xv_dir(:, C_dir%zs + 1:, :, :) = real(&
            & edir(:, atmk(atm, C_one_atm1%zs + 1):C_one_atm1%ze, :, :), &
            & kind(xv_dir))
          xv_dir(:, C_dir%zs, :, :) = real(edir(:, C_one_atm1%zs, :, :), kind(xv_dir))

          call restoreVecPointer(C_dir%da, edir_local, xv_dir1d, xv_dir)

          call VecSet(solution%edir, 0._ireals, ierr); call CHKERR(ierr)
          call DMLocalToGlobalBegin(C_dir%da, edir_local, ADD_VALUES, solution%edir, ierr); call CHKERR(ierr)
          call DMLocalToGlobalEnd(C_dir%da, edir_local, ADD_VALUES, solution%edir, ierr); call CHKERR(ierr)
          call DMRestoreLocalVector(C_dir%da, edir_local, ierr); call CHKERR(ierr)

          call PetscObjectSetName(solution%edir, 'edir', ierr); call CHKERR(ierr)
          call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, '-mcdmda_show_edir', ierr); call CHKERR(ierr)
        end if

        ! handle normalization of ediff
        call DMGetLocalVector(C_diff%da, ediff_local, ierr); call CHKERR(ierr)
        call VecSet(ediff_local, zero, ierr); call CHKERR(ierr)

        call getVecPointer(C_diff%da, ediff_local, xv_diff1d, xv_diff)

        if (.not. solution%lsolar_rad) then
          Nundersampled = count(Nediff .gt. 0_iintegers .and. Nediff .le. 10)
          call CHKWARN(int(Nundersampled, mpiint), 'Found '//toStr(Nundersampled)//' fluxes that were certainly reacheable,'// &
            & 'but were sampled crudely. This may lead to large biases in the results. Please consider using more photons!')

          ediff = ediff * real(pi, kind(ediff)) / real(max(1_iintegers, Nediff), kind(ediff))
        end if

        xv_diff(:, C_diff%zs + 1:, :, :) = real(&
          & ediff(:, atmk(atm, C_one_atm1%zs + 1):C_one_atm1%ze, :, :), &
          & kind(xv_diff))
        xv_diff(:, C_diff%zs, :, :) = real(ediff(:, C_one_atm1%zs, :, :), kind(xv_diff))

        call restoreVecPointer(C_diff%da, ediff_local, xv_diff1d, xv_diff)

        call VecSet(solution%ediff, 0._ireals, ierr); call CHKERR(ierr)
        call DMLocalToGlobalBegin(C_diff%da, ediff_local, ADD_VALUES, solution%ediff, ierr); call CHKERR(ierr)
        call DMLocalToGlobalEnd(C_diff%da, ediff_local, ADD_VALUES, solution%ediff, ierr); call CHKERR(ierr)
        call DMRestoreLocalVector(C_diff%da, ediff_local, ierr); call CHKERR(ierr)

        call PetscObjectSetName(solution%ediff, 'ediff', ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, '-mcdmda_show_ediff', ierr); call CHKERR(ierr)

        ! absorption
        call VecSet(solution%abso, 0._ireals, ierr); call CHKERR(ierr)
        call getVecPointer(C_one%da, solution%abso, xv_abso1d, xv_abso)

        if (.not. solution%lsolar_rad) then
          abso = abso * pi * C_one%glob_xm * C_one%glob_ym / Nphotons_global
        end if

        xv_abso(i0, :, :, :) = real(&
          & abso(i0, atmk(atm, C_one_atm%zs):C_one_atm%ze, :, :) &
          & / atm%dz(atmk(atm, C_one_atm%zs):C_one_atm%ze, :, :), &
          & kind(xv_abso))

        ! eventually collapsed entry at the top
        xv_abso(i0, C_one%zs, :, :) = real( &
          & sum(abso(i0, :atmk(atm, C_one_atm%zs), :, :), dim=1) &
          & / sum(atm%dz(:atmk(atm, C_one_atm%zs), :, :), dim=1), kind(xv_abso))

        call restoreVecPointer(C_one%da, solution%abso, xv_abso1d, xv_abso)
      end associate

      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-mcdmda_show_abso', ierr); call CHKERR(ierr)

      !Rayli solver returns fluxes as [W]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is up to date (to prevent absoprtion computations)
      solution%lchanged = .false.
    end subroutine

    subroutine debug_output_queues()
      integer(mpiint) :: myid, numnodes, ierr, i
      call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)
      do i = 0, numnodes - 1
        if (i .eq. myid) then
          print *, '---------- rank ', myid
          call print_pqueue(pqueues(PQ_SELF))
          call print_pqueue(pqueues(PQ_NORTH))
          call print_pqueue(pqueues(PQ_EAST))
          call print_pqueue(pqueues(PQ_SOUTH))
          call print_pqueue(pqueues(PQ_WEST))
        end if
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      end do
    end subroutine
  end subroutine

  subroutine determine_Nphotons(solver, Nphotons_local, ierr)
    class(t_solver), intent(in) :: solver
    integer(int64), intent(out) :: Nphotons_local
    integer(mpiint), intent(out) :: ierr
    integer(int64) :: mcdmda_photons_per_pixel ! has to be constant over all TOA faces
    real(ireals) :: rmcdmda_photons_per_pixel ! has to be constant over all TOA faces
    real(ireals) :: rN
    logical :: lflg
    integer(mpiint) :: numnodes

    ierr = 0

    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    mcdmda_photons_per_pixel = 1000
    call get_petsc_opt(solver%prefix, "-mcdmda_photons_per_px", &
                       rmcdmda_photons_per_pixel, lflg, ierr); call CHKERR(ierr)
    if (lflg) mcdmda_photons_per_pixel = int(rmcdmda_photons_per_pixel, int64)
    if (ldebug) print *, 'mcdmda_photons_per_pixel', mcdmda_photons_per_pixel

    call get_petsc_opt(solver%prefix, "-mcdmda_photons", rN, lflg, ierr); call CHKERR(ierr)
    if (lflg) then
      mcdmda_photons_per_pixel = 1_int64 + int((rN / solver%C_one_atm%glob_xm) / solver%C_one_atm%glob_ym, kind(Nphotons_local))
      if (ldebug) print *, 'override by -mcdmda_photons', rN, '-> mcdmda_photons_per_pixel ', mcdmda_photons_per_pixel
    end if

    Nphotons_local = int(solver%C_one_atm%xm, int64) * int(solver%C_one_atm%ym, int64) * mcdmda_photons_per_pixel
    if (ldebug) print *, 'Nphotons_local', Nphotons_local
  end subroutine

  subroutine photon_queue_destroy(pq, ierr)
    type(t_photon_queue), intent(inout) :: pq
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    if (allocated(pq%photons)) deallocate (pq%photons)
    call pq%ready%finalize()
    call pq%empty%finalize()
    call pq%sending%finalize()
  end subroutine

  subroutine run_photon_queue(solver, bmc, bmc_1_2, &
      & pqueues, ipq, &
      & edir, ediff, Nediff, abso, &
      & started_photons, killed_photons, &
      & limit_number_photons, &
      & opt_buildings, buildings_idx)
    class(t_solver), intent(in) :: solver
    class(t_boxmc), intent(in) :: bmc
    class(t_boxmc_1_2), intent(in) :: bmc_1_2
    type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
    integer(iintegers), intent(in) :: ipq
    real(ireal_dp), allocatable, dimension(:, :, :, :), intent(inout) :: edir, ediff, abso
    integer(iintegers), allocatable, dimension(:, :, :, :), intent(inout) :: Nediff
    integer(int64), intent(out) :: started_photons
    integer(int64), intent(out) :: killed_photons
    integer(int64), optional, intent(in) :: limit_number_photons
    type(t_pprts_buildings), intent(in), optional :: opt_buildings
    integer(iintegers), allocatable, dimension(:, :, :, :), intent(in) :: buildings_idx
    integer(int64) :: Nphotmax
    logical :: lkilled_photon

    if (ldebug) then
      print *, 'run photon queue ', id2name(pqueues(ipq)%queue_index)
      call print_pqueue(pqueues(ipq))
    end if

    Nphotmax = get_arg(huge(Nphotmax), limit_number_photons)
    Nphotmax = min(size(pqueues(ipq)%photons, kind=int64), Nphotmax)

    killed_photons = 0
    started_photons = 0
    call pqueues(ipq)%ready%for_each(run_ready)

  contains
    subroutine run_ready(idx, node, iphoton)
      integer(int64), intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(int64), intent(inout) :: iphoton
      integer(mpiint) :: ierr

      if (started_photons .ge. Nphotmax) return

      call run_photon(&
        & solver, &
        & bmc, bmc_1_2, &
        & pqueues, &
        & ipq, &
        & iphoton, &
        & edir, &
        & ediff, &
        & Nediff, &
        & abso, &
        & lkilled_photon, &
        & opt_buildings=opt_buildings, &
        & buildings_idx=buildings_idx)

      started_photons = started_photons + 1
      if (lkilled_photon) killed_photons = killed_photons + 1

      call pqueues(ipq)%ready%del_node(node, ierr); call CHKERR(ierr)
      return
      print *, 'remove unused var warning', idx
    end subroutine
  end subroutine

  subroutine run_photon(solver, bmc, bmc_1_2, pqueues, ipq, iphoton, &
      & edir, ediff, Nediff, abso, &
      & lkilled_photon, &
      & opt_buildings, buildings_idx)
    class(t_solver), intent(in) :: solver
    class(t_boxmc), intent(in) :: bmc
    class(t_boxmc_1_2), intent(in) :: bmc_1_2
    type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
    integer(iintegers), intent(in) :: ipq
    integer(int64), intent(in) :: iphoton
    real(ireal_dp), allocatable, dimension(:, :, :, :), intent(inout) :: edir, ediff, abso
    integer(iintegers), allocatable, dimension(:, :, :, :), intent(inout) :: Nediff
    logical, intent(out) :: lkilled_photon
    type(t_pprts_buildings), intent(in), optional :: opt_buildings
    integer(iintegers), allocatable, dimension(:, :, :, :), intent(in) :: buildings_idx

    logical :: lexit_cell, lthermal, lexit_domain
    integer(mpiint) :: myid, ierr

    lkilled_photon = .false.

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    associate (p => pqueues(ipq)%photons(iphoton)%p)

      if (ldebug_tracing) print *, myid, cstr('Start of run_photon :: QUEUE:', 'pink'), id2name(ipq), 'iphoton', iphoton
      if (ldebug_tracing) call print_photon(p)
      call check_if_photon_is_in_domain(solver%C_one_atm, p)

      lthermal = allocated(solver%atm%planck)

      if (p%src_side .eq. PPRTS_TOP_FACE) then ! started at the top of the box, lets increment TOA downward flux
        if (p%direct) then
          p%side = PPRTS_TOP_FACE
          call update_flx(solver, p, p%k, p%i, p%j, edir, ediff, Nediff)
        end if
      end if

      lexit_domain = .false.
      lexit_cell = .true. ! mark es exited because it obviously came from some neighbouring box

      move: do while (.not. lexit_domain) ! this loop will move the photon to the edge of the subdomain
        if (.not. lthermal) then
          !call roulette(p)
          if (.not. p%alive) then
            lkilled_photon = .true.
            exit move
          end if
        end if
        if (ldebug_tracing) then
          print *, 'start of move', p%k, p%i, p%j
          call print_photon(p)
        end if

        call check_if_photon_is_in_domain(solver%C_one_atm, p)

        if (lexit_cell) then
          ! if photon exited cell before, this is the reentrant state
          abso(i0, p%k, p%i, p%j) = abso(i0, p%k, p%i, p%j) + real(p%weight, ireals)
        end if

        if (present(opt_buildings)) then
          call building_interaction(lexit_cell)
        else
          lexit_cell = .false.
        end if

        if (.not. lexit_cell) call move_inside_cell(lexit_cell)

        ! Define new cellindex and update local position
        if (.not. lexit_cell) then
          call scatter_photon_in_cell()
        else
          abso(i0, p%k, p%i, p%j) = abso(i0, p%k, p%i, p%j) - real(p%weight, ireals)
          call exit_cell(lexit_domain)
        end if

      end do move
      call pqueues(ipq)%empty%add(iphoton)
    end associate

  contains

    subroutine scatter_photon_in_cell()
      associate (p => pqueues(ipq)%photons(iphoton)%p)
        call scatter_photon(p, real(solver%atm%g(p%k, p%i, p%j), ireal_dp))
        p%tau_travel = tau(R())
        if (ldebug_tracing) print *, myid, cstr('********************* SCATTERING', 'peach'), p%k, p%i, p%j
      end associate
    end subroutine

    subroutine move_inside_cell(lexit_cell)
      logical, intent(out) :: lexit_cell

      real(ireal_dp) :: kabs, ksca, pathlen
      real(ireal_dp) :: Btop, Bbot, B1, B2, dz, tauabs, tm1
      real(ireal_dp), allocatable :: vertices(:)

      associate (p => pqueues(ipq)%photons(iphoton)%p)

        dz = solver%atm%dz(p%k, p%i, p%j)
        call setup_default_unit_cube_geometry( &
          & real(solver%atm%dx, ireal_dp), &
          & real(solver%atm%dy, ireal_dp), &
          & dz, vertices)

        if (.not. pnt_in_cube(vertices, p%loc)) then
          ! on warped grids we encounter local photon coordinates that dont match the uniform cube
          ! in this case, one option would be to use global coordinates and search for the correct vertical index
          ! however, this brings along the need for epsilons which we try hard to get rid of when tracking cell ids
          ! a dirty hack here just moves the photon inside the cube vertically. On strongly warped grids, this introduces errors.
          ! However, I currently don't have a good idea how to do it. Use RayLi if you need this to be handled.
          p%loc(3) = min(max(p%loc(3), 0._ireal_dp), dz - 1e-3_ireal_dp)
        end if

        kabs = solver%atm%kabs(p%k, p%i, p%j)
        ksca = solver%atm%ksca(p%k, p%i, p%j)
        if (lthermal) then
          Btop = solver%atm%planck(p%k, p%i, p%j)
          Bbot = solver%atm%planck(p%k + 1, p%i, p%j)
          B1 = (Btop * p%loc(3) + Bbot * (dz - p%loc(3))) / dz ! planck at start of the ray
        end if

        call move_photon(solver%atm, bmc, bmc_1_2, vertices, ksca, p, pathlen, lexit_cell)

        if (lthermal) then
          B2 = (Btop * p%loc(3) + Bbot * (dz - p%loc(3))) / dz ! planck at end of the ray
          tauabs = kabs * pathlen
          if (tauabs > 1e-10_ireal_dp) then
            tm1 = expm1(-tauabs)
            p%weight = p%weight * (tm1 + 1._ireal_dp) + (B2 - B1) - (B1 - (B2 - B1) / tauabs) * tm1
          else
            p%weight = p%weight * (1._ireal_dp - tauabs) + (B1 + B2)*.5_ireal_dp * tauabs
          end if
          !print *,p%k, p%i, p%j,'Btop/bot', Btop, Bbot, 'B1,2', B1, B2, 'weight', p%weight, 'kabs', kabs, 'ksca', ksca
        else ! lsolar
          call absorb_photon(p, pathlen, kabs)
        end if
      end associate
    end subroutine

    subroutine building_interaction(lexit_cell)
      logical, intent(out) :: lexit_cell
      integer(iintegers) :: bidx
      real(ireal_dp) :: mu, phi

      lexit_cell = .false.

      associate (p => pqueues(ipq)%photons(iphoton)%p)
        ! p%src_side - side where it is coming from, i.e. where it starts,
        ! .e.g if flying down, i.e. entering from top, it will be PPRTS_TOP_FACE
        bidx = buildings_idx(p%src_side, p%k, p%i, p%j)
        if (bidx .gt. 0) then

          if (ldebug_tracing) print *, 'hit building @ '//toStr([p%src_side, p%k, p%i, p%j])
          p%direct = .false.
          p%scattercnt = 0

          p%weight = p%weight * opt_buildings%albedo(bidx)
          if (lthermal) then
            p%weight = p%weight + (1._ireals - opt_buildings%albedo(bidx)) * opt_buildings%planck(bidx)
          end if

          ! send right back to where it came from:
          p%side = p%src_side

          mu = sqrt(R())
          phi = deg2rad(R() * 360)
          p%dir = [sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu]

          select case (p%src_side)
          case (PPRTS_TOP_FACE)
            continue

          case (PPRTS_BOT_FACE)
            p%dir = rotate_angle_y(p%dir, 180._ireal_dp)

          case (PPRTS_LEFT_FACE)
            p%dir = rotate_angle_y(p%dir, 90._ireal_dp)

          case (PPRTS_RIGHT_FACE)
            p%dir = rotate_angle_y(p%dir, 270._ireal_dp)

          case (PPRTS_REAR_FACE)
            p%dir = rotate_angle_x(p%dir, 270._ireal_dp)

          case (PPRTS_FRONT_FACE)
            p%dir = rotate_angle_x(p%dir, 90._ireal_dp)

          case default
            call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(p%src_side))
          end select

          lexit_cell = .true.
        end if
      end associate
    end subroutine

    subroutine exit_cell(lexit_domain)
      logical, intent(out) :: lexit_domain
      real(ireal_dp) :: mu, phi

      lexit_domain = .false.

      associate (p => pqueues(ipq)%photons(iphoton)%p)

        call update_flx(solver, p, p%k, p%i, p%j, edir, ediff, Nediff)

        ! Move photon to new cell
        if (ldebug_tracing) print *, myid, cstr('* MOVE Photon to new cell', 'green')

        ! reset scatter count in new cell, remove scatter info from last box,
        ! i.e. don't allow intersections with source side immediately
        p%scattercnt = 0

        select case (p%side)

        case (PPRTS_TOP_FACE) ! exit on top

          if (p%k .eq. solver%C_one_atm%zs) then ! outgoing at TOA
            if (ldebug_tracing) print *, myid, '********************* Exit TOA', p%k, p%i, p%j
            lkilled_photon = .true.
            lexit_domain = .true.
          end if

          p%loc(3) = zero + loceps
          p%k = p%k - 1
          p%src_side = PPRTS_BOT_FACE

        case (PPRTS_BOT_FACE) ! exit cell at bottom

          if (p%k .eq. solver%C_one_atm%ze) then ! hit the surface, need reflection
            if (ldebug_tracing) print *, myid, '********************* Before Reflection', p%k, p%i, p%j

            p%weight = p%weight * solver%atm%albedo(p%i, p%j)
            if (lthermal) then
              if (allocated(solver%atm%Bsrfc)) then
                p%weight = p%weight + (1._ireals - solver%atm%albedo(p%i, p%j)) * solver%atm%Bsrfc(p%i, p%j)
              else
                p%weight = p%weight + (1._ireals - solver%atm%albedo(p%i, p%j)) * solver%atm%planck(p%k + 1, p%i, p%j)
              end if
            end if
            p%direct = .false.
            !p%scattercnt = p%scattercnt + 1

            mu = sqrt(R())
            phi = deg2rad(R() * 360)
            p%dir = [sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu]
            p%loc(3) = zero + loceps

            p%side = PPRTS_TOP_FACE
            p%src_side = PPRTS_BOT_FACE
            call update_flx(solver, p, p%k + 1, p%i, p%j, edir, ediff, Nediff)
            if (ldebug_tracing) print *, myid, cstr('********************* After  Reflection', 'aqua'), p%k, p%i, p%j
            lexit_domain = .false.

          else

            p%loc(3) = solver%atm%dz(p%k + 1, p%i, p%j) - loceps
            p%k = p%k + 1
            p%src_side = PPRTS_TOP_FACE

          end if

        case (PPRTS_LEFT_FACE)
          p%loc(1) = solver%atm%dx - loceps
          p%i = p%i - 1
          p%src_side = PPRTS_RIGHT_FACE
          if (p%i .eq. solver%C_one_atm%xs - 1) then
            if (ldebug_tracing) print *, myid, cstr('* Sending to WEST', 'blue'), pqueues(PQ_WEST)%owner, p%k, p%i, p%j
            call send_photon_to_neighbor(solver, solver%C_one_atm, p, pqueues(PQ_WEST))
            lexit_domain = .true.
          end if
        case (PPRTS_RIGHT_FACE)
          p%loc(1) = zero + loceps
          p%i = p%i + 1
          p%src_side = PPRTS_LEFT_FACE
          if (p%i .eq. solver%C_one_atm%xe + 1) then
            if (ldebug_tracing) print *, myid, cstr('* Sending to EAST', 'blue'), pqueues(PQ_EAST)%owner, p%k, p%i, p%j
            call send_photon_to_neighbor(solver, solver%C_one_atm, p, pqueues(PQ_EAST))
            lexit_domain = .true.
          end if
        case (PPRTS_REAR_FACE)
          p%loc(2) = solver%atm%dy - loceps
          p%j = p%j - 1
          p%src_side = PPRTS_FRONT_FACE
          if (p%j .eq. solver%C_one_atm%ys - 1) then
            if (ldebug_tracing) print *, myid, cstr('* Sending to SOUTH', 'blue'), pqueues(PQ_SOUTH)%owner, p%k, p%i, p%j
            call send_photon_to_neighbor(solver, solver%C_one_atm, p, pqueues(PQ_SOUTH))
            lexit_domain = .true.
          end if
        case (PPRTS_FRONT_FACE)
          p%loc(2) = zero + loceps
          p%j = p%j + 1
          p%src_side = PPRTS_REAR_FACE
          if (p%j .eq. solver%C_one_atm%ye + 1) then
            if (ldebug_tracing) print *, myid, cstr('* Sending to NORTH', 'blue'), pqueues(PQ_NORTH)%owner, p%k, p%i, p%j
            call send_photon_to_neighbor(solver, solver%C_one_atm, p, pqueues(PQ_NORTH))
            lexit_domain = .true.
          end if
        end select
      end associate
    end subroutine
  end subroutine

  subroutine update_flx(solver, p, k, i, j, edir, ediff, Nediff)
    class(t_solver), intent(in) :: solver
    type(t_photon), intent(in) :: p
    integer(iintegers), intent(in) :: k, i, j ! layer/box indices
    real(ireal_dp), allocatable, dimension(:, :, :, :), intent(inout) :: edir, ediff
    integer(iintegers), allocatable, dimension(:, :, :, :), intent(inout) :: Nediff
    integer(iintegers) :: off, dof

    if (ldebug_tracing) print *, 'Update Flux', k, i, j, p%direct, p%side

    select case (p%side)
    case (PPRTS_TOP_FACE)
      if (k .lt. lbound(ediff, 2) .or. k .gt. ubound(ediff, 2)) &
        & call CHKERR(1_mpiint, 'invalid k '//toStr(k)//' '// &
                       & toStr(lbound(ediff, 2))//'/'//toStr(ubound(ediff, 2)))
    case (PPRTS_BOT_FACE)
      if (k .lt. lbound(ediff, 2) .or. k .gt. ubound(ediff, 2) - 1) &
        & call CHKERR(1_mpiint, 'invalid k '//toStr(k)//' '// &
                      & toStr(lbound(ediff, 2))//'/'//toStr(ubound(ediff, 2) - 1))
    case (3:6)
      continue
    case default
      call print_photon(p)
      call CHKERR(1_mpiint, 'hmpf .. didnt expect a p%side gt 6, have '//toStr(p%side))
    end select

    if (i .lt. lbound(ediff, 3) .or. i .gt. ubound(ediff, 3)) &
      & call CHKERR(1_mpiint, 'invalid i '//toStr(i)//' '// &
                    & '['//toStr(lbound(ediff, 3))//','//toStr(ubound(ediff, 3))//']')
    if (j .lt. lbound(ediff, 4) .or. j .gt. ubound(ediff, 4)) &
      & call CHKERR(1_mpiint, 'invalid j '//toStr(j)//' '// &
                    & '['//toStr(lbound(ediff, 4))//','//toStr(ubound(ediff, 4))//']')

    if (p%direct) then
      select case (p%side)
      case (PPRTS_TOP_FACE)
        edir(i0, k, i, j) = edir(i0, k, i, j) + real(p%weight, kind(edir))
      case (PPRTS_BOT_FACE)
        edir(i0, k + 1, i, j) = edir(i0, k + 1, i, j) + real(p%weight, kind(edir))

      case (PPRTS_LEFT_FACE)
        off = solver%dirtop%dof
        do dof = 0, solver%dirside%dof - 1
          edir(off + dof, k, i, j) = edir(off + dof, k, i, j) + real(p%weight, kind(edir))
        end do
      case (PPRTS_RIGHT_FACE)
        off = solver%dirtop%dof
        do dof = 0, solver%dirside%dof - 1
          edir(off + dof, k, i + 1, j) = edir(off + dof, k, i + 1, j) + real(p%weight, kind(edir))
        end do

      case (PPRTS_REAR_FACE) ! rear
        off = solver%dirtop%dof + solver%dirside%dof
        do dof = 0, solver%dirside%dof - 1
          edir(off + dof, k, i, j) = edir(off + dof, k, i, j) + real(p%weight, kind(edir))
        end do
      case (PPRTS_FRONT_FACE) ! front
        off = solver%dirtop%dof + solver%dirside%dof
        do dof = 0, solver%dirside%dof - 1
          edir(off + dof, k, i, j + 1) = edir(off + dof, k, i, j + 1) + real(p%weight, kind(edir))
        end do
      case default
        call print_photon(p)
        call CHKERR(1_mpiint, 'hmpf .. didnt expect p%side '//toStr(p%side))
      end select

    else

      select case (p%side)
      case (PPRTS_TOP_FACE)
        do dof = 0, solver%difftop%dof - 1
          if (.not. solver%difftop%is_inward(i1 + dof)) then !Eup
            ediff(dof, k, i, j) = ediff(dof, k, i, j) + real(p%weight, kind(ediff))
            Nediff(dof, k, i, j) = Nediff(dof, k, i, j) + 1_iintegers
          end if
        end do
      case (PPRTS_BOT_FACE)
        do dof = 0, solver%difftop%dof - 1
          if (solver%difftop%is_inward(i1 + dof)) then !Edn
            ediff(dof, k + 1, i, j) = ediff(dof, k + 1, i, j) + real(p%weight, kind(ediff))
            Nediff(dof, k + 1, i, j) = Nediff(dof, k + 1, i, j) + 1_iintegers
          end if
        end do

      case (PPRTS_LEFT_FACE)
        off = solver%difftop%dof
        do dof = 0, solver%diffside%dof - 1
          if (.not. solver%diffside%is_inward(i1 + dof)) then !Eleft
            ediff(off + dof, k, i, j) = ediff(off + dof, k, i, j) + real(p%weight, kind(ediff))
            Nediff(off + dof, k, i, j) = Nediff(off + dof, k, i, j) + 1_iintegers
          end if
        end do
      case (PPRTS_RIGHT_FACE)
        off = solver%difftop%dof
        do dof = 0, solver%diffside%dof - 1
          if (solver%diffside%is_inward(i1 + dof)) then !Eright
            ediff(off + dof, k, i + 1, j) = ediff(off + dof, k, i + 1, j) + real(p%weight, kind(ediff))
            Nediff(off + dof, k, i + 1, j) = Nediff(off + dof, k, i + 1, j) + 1_iintegers
          end if
        end do

      case (PPRTS_REAR_FACE) ! rear
        off = solver%difftop%dof + solver%diffside%dof
        do dof = 0, solver%diffside%dof - 1
          if (.not. solver%diffside%is_inward(i1 + dof)) then !Erear
            ediff(off + dof, k, i, j) = ediff(off + dof, k, i, j) + real(p%weight, kind(ediff))
            Nediff(off + dof, k, i, j) = Nediff(off + dof, k, i, j) + 1_iintegers
          end if
        end do
      case (PPRTS_FRONT_FACE) ! front
        off = solver%difftop%dof + solver%diffside%dof
        do dof = 0, solver%diffside%dof - 1
          if (solver%diffside%is_inward(i1 + dof)) then !Eforward
            ediff(off + dof, k, i, j + 1) = ediff(off + dof, k, i, j + 1) + real(p%weight, kind(ediff))
            Nediff(off + dof, k, i, j + 1) = Nediff(off + dof, k, i, j + 1) + 1_iintegers
          end if
        end do

      case default
        call print_photon(p)
        call CHKERR(1_mpiint, 'hmpf .. didnt expect this p%side, have'//toStr(p%side))
      end select
    end if
  end subroutine

  subroutine prepare_locally_owned_photons(solver, bmc, lsolar, pqueue, Nphotons, weight)
    class(t_solver), intent(in) :: solver
    class(t_boxmc), intent(in) :: bmc
    logical, intent(in) :: lsolar
    type(t_photon_queue), intent(inout) :: pqueue
    integer(int64), intent(in) :: Nphotons
    real(ireals), intent(in) :: weight

    type(t_photon) :: p
    real(ireals) :: phi0, theta0
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: i, j
    integer(int64) :: ip, Nphotons_per_pixel, l
    real(ireals), allocatable :: vertices(:)
    real(ireals) :: initial_dir(3)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    phi0 = solver%sun%phi
    theta0 = solver%sun%theta
    initial_dir = spherical_2_cartesian(phi0, theta0)

    Nphotons_per_pixel = max(1_int64, Nphotons / int(solver%C_one_atm%xm * solver%C_one_atm%ym, kind(Nphotons)))
    if (modulo(Nphotons, Nphotons_per_pixel) .ne. 0) &
      & call CHKERR(1_mpiint, 'Nphotons '//toStr(Nphotons)//' not divisible by Nphotons_per_pixel '//toStr(Nphotons_per_pixel))

    do l = 1, Nphotons_per_pixel
      do i = solver%C_one_atm%xs, solver%C_one_atm%xe
        do j = solver%C_one_atm%ys, solver%C_one_atm%ye

          call setup_default_unit_cube_geometry(solver%atm%dx, solver%atm%dy, &
                                                solver%atm%dz(i0, i, j), vertices)

          if (lsolar) then
            call bmc%init_dir_photon(p, i1, .true., real(initial_dir, ireal_dp), real(vertices, ireal_dp), ierr)
          else
            call bmc%init_diff_photon(p, i2, real(vertices, ireal_dp), ierr)
          end if
          p%i = i
          p%j = j
          p%k = i0
          p%src_side = PPRTS_TOP_FACE
          p%weight = weight
          p%tau_travel = tau(R())
          if (ldebug_tracing) then
            call antialiased_photon_start(Nphotons_per_pixel, l, p%loc(1), p%loc(2))
          else
            p%loc(1) = R()
            p%loc(2) = R()
          end if
          p%loc(1) = p%loc(1) * solver%atm%dx
          p%loc(2) = p%loc(2) * solver%atm%dy

          call pqueue_add_photon(pqueue, p, 0_mpiint, ip, ierr); call CHKERR(ierr)
          call pqueue%ready%add(ip)

          if (ldebug) then
            print *, 'Prepared Photon', i, j, ': iphoton', ip, Nphotons
            call print_photon(pqueue%photons(ip)%p)
          end if
        end do
      end do
    end do
  end subroutine

  subroutine antialiased_photon_start(Nmax, ip, x, y, tilt_angle)
    integer(int64), intent(in) :: Nmax, ip ! total number of photons per pixel, and the ip'th point for that
    real(ireal_dp), intent(in), optional :: tilt_angle ! angle by which the grid is rotated in [rad], default: 30deg
    real(ireal_dp), intent(out) :: x, y ! gives x and y position for a regularly sampled tilted grid

    real(ireal_dp) :: tilt_grid, rot_x, rot_y
    integer(int64) :: Nx, Ny ! number of pixels in x and y direction
    integer(int64) :: i, j

    tilt_grid = get_arg(deg2rad(26.6_ireal_dp), tilt_angle)

    if (Nmax .eq. 1_int64) then
      x = .5_ireals
      y = .5_ireals
      return
    end if

    Nx = int(sqrt(real(Nmax, ireal_dp)), int64)
    Ny = Nx

    if (ip .gt. Nx * Ny) then
      x = R()
      y = R()
      return
    end if

    j = (ip - 1) / int(Nx, kind(ip)) ! zero based offsets
    i = (ip - 1) - j * int(Nx, kind(ip))

    x = real(i + 1, ireals) / real(Nx + 1, ireals) * sqrt(5._ireals) / 2._ireals
    y = real(j + 1, ireals) / real(Ny + 1, ireals) * sqrt(5._ireals) / 2._ireals

    ! Translate to center (0,0)
    x = x - .5_ireals
    y = y - .5_ireals

    rot_x = x * cos(tilt_grid) - y * sin(tilt_grid) + .5_ireals
    rot_y = x * sin(tilt_grid) + y * cos(tilt_grid) + .5_ireals

    x = modulo(rot_x, 1._ireal_dp)
    y = modulo(rot_y, 1._ireal_dp)

    !print *,'plot(', x,',',y,',"o") #',ip
  end subroutine

  subroutine find_empty_entry_in_pqueue(pqueue, emptyid, ierr)
    type(t_photon_queue), intent(inout) :: pqueue
    integer(int64), intent(out) :: emptyid
    integer(mpiint), intent(out) :: ierr
    call pqueue%empty%pop(emptyid, ierr)
    if (ierr .ne. 0) then
      call print_pqueue(pqueue)
    end if
  end subroutine

  subroutine print_pqueue(pq, list_queue)
    type(t_photon_queue), intent(in) :: pq
    logical, intent(in), optional :: list_queue
    if (get_arg(.false., list_queue)) then
      print *, 'PQUEUE:', pq%owner, '::', pq%queue_index, 'name ', cstr(trim(id2name(pq%queue_index)), 'blue')
      print *, 'empty:'
      call pq%empty%view()
      print *, 'ready:'
      call pq%ready%view()
      print *, 'sending:'
      call pq%sending%view()
    else
      print *, pq%owner, cstr(id2name(pq%queue_index), 'blue'), &
        & '  '//'empty', pq%empty%len(), &
        & '  '//cstr('ready '//toStr(pq%ready%len()), 'green'), &
        & '  '//cstr('send  '//toStr(pq%sending%len()), 'peach')
    end if
  end subroutine

  subroutine pqueue_add_photon(pqueue, p, request, ind, ierr)
    type(t_photon_queue), intent(inout) :: pqueue
    type(t_photon), intent(in) :: p
    integer(mpiint), intent(in) :: request
    integer(int64), intent(out) :: ind
    integer(mpiint), intent(out) :: ierr

    call find_empty_entry_in_pqueue(pqueue, ind, ierr)
    if (ierr .ne. 0) then
      call finalize_msgs_blocking(pqueue)
      call find_empty_entry_in_pqueue(pqueue, ind, ierr)
      call CHKERR(ierr, 'Could not find an empty slot in neighbor queue ', pqueue%queue_index)
    end if

    ! Put photon to send queue to neighbor
    pqueue%photons(ind)%p = p
    pqueue%photons(ind)%request = request
  end subroutine

  subroutine setup_photon_queue(pq, N, owner, queue_index)
    type(t_photon_queue), intent(inout) :: pq
    integer(int64), intent(in) :: N
    integer(mpiint), intent(in) :: owner
    integer(iintegers), intent(in) :: queue_index
    integer(int64) :: i

    if (allocated(pq%photons)) then
      call CHKERR(1_mpiint, 'photon queue already allocated')
    else
      allocate (pq%photons(N))
    end if
    pq%owner = owner
    pq%queue_index = queue_index
    do i = 1, N
      call pq%empty%add(i)
    end do
  end subroutine

  subroutine send_photon_to_neighbor(solver, C, p_in, pqueue)
    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    type(t_photon), intent(in) :: p_in
    type(t_photon_queue), intent(inout) :: pqueue
    integer(mpiint) :: myid, tag, ierr
    integer(int64) :: iphoton

    logical, parameter :: lcyclic_boundary = .true.

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    ! Put photon to send queue to neighbor
    call pqueue_add_photon(pqueue, p_in, 0_mpiint, iphoton, ierr); call CHKERR(ierr)
    call pqueue%sending%add(iphoton)

    associate (p => pqueue%photons(iphoton)%p)

      if (ldebug) then
        select case (pqueue%queue_index)
        case (PQ_SELF)
          call CHKERR(1_mpiint, 'should not happen that a photon is sent to ourselves?')
        case (PQ_WEST)
          !print *,'Sending Photon WEST'
          if (p%side .ne. PPRTS_LEFT_FACE) call CHKERR(1_mpiint, 'I would assume that side should be 3 when sending WEST')
          if (p%src_side .ne. PPRTS_RIGHT_FACE) call CHKERR(1_mpiint, 'I would assume that src_side should be 4 when sending WEST')
        case (PQ_EAST)
          !print *,'Sending Photon EAST'
          if (p%side .ne. PPRTS_RIGHT_FACE) call CHKERR(1_mpiint, 'I would assume that side should be 4 when sending EAST')
          if (p%src_side .ne. PPRTS_LEFT_FACE) call CHKERR(1_mpiint, 'I would assume that src_side should be 3 when sending EAST')
        case (PQ_NORTH)
          !print *,'Sending Photon NORTH'
          if (p%side .ne. PPRTS_FRONT_FACE) call CHKERR(1_mpiint, 'I would assume that side should be 6 when sending NORTH')
          if (p%src_side .ne. PPRTS_REAR_FACE) call CHKERR(1_mpiint, 'I would assume that src_side should be 5 when sending NORTH')
        case (PQ_SOUTH)
          !print *,'Sending Photon SOUTH'
          if (p%side .ne. PPRTS_REAR_FACE) call CHKERR(1_mpiint, 'I would assume that side should be 5 when sending SOUTH')
          if (p%src_side .ne. PPRTS_FRONT_FACE) call CHKERR(1_mpiint, 'I would assume that src_side should be 6 when sending SOUTH')
        end select
      end if

      if (lcyclic_boundary) then
        if (p%i .eq. -1) p%i = C%glob_xm - 1
        if (p%j .eq. -1) p%j = C%glob_ym - 1
        if (p%i .eq. C%glob_xm) p%i = 0
        if (p%j .eq. C%glob_ym) p%j = 0
      else
        call CHKERR(1_mpiint, 'NON-cyclic boundary conditions are not implemented yet!')
      end if

      ! asynchronous SEND starts here
      tag = int(pqueue%queue_index, kind(tag))
      call mpi_isend(p, 1_mpiint, imp_t_photon, &
                     pqueue%owner, tag, solver%comm, pqueue%photons(iphoton)%request, ierr); call CHKERR(ierr, 'mpi isend failed')
      !call mpi_send(p, 1_mpiint, imp_t_photon, &
      !  pqueue%owner, tag, solver%comm, ierr); call CHKERR(ierr, 'mpi isend failed')

      if (ldebug) then
        print *, 'Sending the following photon to rank', pqueue%owner, ':tag', tag
        call print_photon(p)
      end if

    end associate
  end subroutine

  subroutine finalize_msgs_non_blocking(pqueue)
    type(t_photon_queue), intent(inout) :: pqueue

    call pqueue%sending%for_each(check_sending)

  contains

    subroutine check_sending(idx, node, iphoton)
      integer(int64), intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(int64), intent(inout) :: iphoton

      integer(mpiint) :: stat(mpi_status_size), ierr
      logical :: lcomm_finished

      call mpi_test(pqueue%photons(iphoton)%request, lcomm_finished, stat, ierr); call CHKERR(ierr)
      if (lcomm_finished) then
        call pqueue%empty%add(iphoton)
        call pqueue%sending%del_node(node, ierr); call CHKERR(ierr)
        if (ldebug) then
          print *, 'Finalized Sending a photon:', iphoton
        end if
      end if
      return
      print *, 'unused var warning', idx
    end subroutine
  end subroutine

  subroutine finalize_msgs_blocking(pqueue)
    type(t_photon_queue), intent(inout) :: pqueue

    integer(iintegers) :: cnt_finished_msgs

    real(ireal_dp) :: tstart, t
    integer(iintegers), save :: count_warnings = 0
    integer(iintegers), parameter :: max_warnings = 10

    call cpu_time(tstart)

    cnt_finished_msgs = 0
    do
      call pqueue%sending%for_each(check_sending)
      call cpu_time(t)
      count_warnings = count_warnings + 1
      if (count_warnings .lt. max_warnings) then
        call CHKWARN(1_mpiint, id2name(pqueue%queue_index)// &
          & ': waited for '//toStr(t - tstart)//' s.'//new_line('')// &
          & ' But this is bad for performance!'//new_line('')// &
        & ' Maybe increasing the queue size helps at the cost of more memory.'//new_line('')// &
          & ' -mcdmda_queue_size <int>'//new_line('')// &
          & ' or try to reduce the batch_size with '//new_line('')// &
          & ' -mcdmda_batch_size <int>'//new_line('')// &
          & ' or try to wait for all boundary photons to be finished emitting new photons locally'//new_line('')// &
          & ' -mcdmda_finish_border_photons_first')
      elseif (count_warnings .eq. max_warnings) then
        call CHKWARN(1_mpiint, 'waiting warning has been issued '//toStr(max_warnings)//' times... now suppressing it.')
      end if
      if (cnt_finished_msgs .ne. 0) return
      if (t - tstart .gt. blocking_waittime) then
        call CHKERR(1_mpiint, 'waited for '//toStr(t - tstart)//' s but the queues havent freed up... ')
      end if
    end do
  contains
    subroutine check_sending(idx, node, iphoton)
      integer(int64), intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(int64), intent(inout) :: iphoton

      integer(mpiint) :: stat(mpi_status_size), ierr
      logical :: lcomm_finished

      call mpi_test(pqueue%photons(iphoton)%request, lcomm_finished, stat, ierr); call CHKERR(ierr)
      if (lcomm_finished) then
        call pqueue%empty%add(iphoton)
        call pqueue%sending%del_node(node, ierr); call CHKERR(ierr)
        cnt_finished_msgs = cnt_finished_msgs + 1
        if (ldebug) then
          print *, 'Finalized Sending a photon:', iphoton
          call print_photon(pqueue%photons(iphoton)%p)
        end if
      end if
      return
      print *, 'unused var warning', idx
    end subroutine
  end subroutine

  subroutine exchange_photons(solver, pqueues)
    class(t_solver), intent(in) :: solver
    type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
    integer(mpiint) :: myid, mpi_status(mpi_status_size), ierr, tag
    integer(int64) :: ipq, iphoton
    logical :: lgot_msg

    call mpi_comm_rank(solver%comm, myid, ierr)

    do ipq = 1, size(pqueues)
      call finalize_msgs_non_blocking(pqueues(ipq))
    end do

    ! receive all messages that we can get
    do
      call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, solver%comm, lgot_msg, mpi_status, ierr); call CHKERR(ierr)
      if (lgot_msg) then
        tag = mpi_status(MPI_TAG)
        select case (tag)
        case (PQ_WEST) ! was sent towards WEST, i.e. it arrives EAST
          ipq = PQ_EAST
        case (PQ_EAST)
          ipq = PQ_WEST
        case (PQ_SOUTH)
          ipq = PQ_NORTH
        case (PQ_NORTH)
          ipq = PQ_SOUTH
        case default
          call CHKERR(1_mpiint, 'received unexpected message with tag '//toStr(tag))
        end select
        if (mpi_status(MPI_SOURCE) .ne. pqueues(ipq)%owner) call CHKERR(1_mpiint, 'Something unexpected happened')

        call find_empty_entry_in_pqueue(pqueues(ipq), iphoton, ierr)
        if (ierr .ne. 0) then
          call finalize_msgs_blocking(pqueues(ipq))
          call find_empty_entry_in_pqueue(pqueues(ipq), iphoton, ierr)
          call CHKERR(ierr, 'no space in queue to receive a msg')
        end if

        call mpi_recv(pqueues(ipq)%photons(iphoton)%p, 1_mpiint, imp_t_photon, &
                      mpi_status(MPI_SOURCE), mpi_status(MPI_TAG), solver%comm, mpi_status, ierr); call CHKERR(ierr)

        call pqueues(ipq)%ready%add(iphoton)

        if (ldebug) then
          print *, myid, 'Got Message from rank:', mpi_status(MPI_SOURCE), 'Receiving ipq ', id2name(ipq)
          call print_photon(pqueues(ipq)%photons(iphoton)%p)
        end if
      else
        exit
      end if
    end do

  end subroutine

  subroutine move_photon(atm, bmc, bmc_1_2, vertices, ksca, p, pathlen, lexit_cell)
    type(t_atmosphere), intent(in) :: atm
    class(t_boxmc) :: bmc
    class(t_boxmc_1_2), intent(in) :: bmc_1_2
    real(ireal_dp), intent(in) :: vertices(:), ksca
    type(t_photon), intent(inout) :: p
    real(ireal_dp) :: pathlen, old_loc(3)
    logical, intent(out) :: lexit_cell

    real(ireal_dp) :: dist, intersec_dist

    if ((.not. lfull3d) .and. atm%l1d(p%k)) then
      old_loc = p%loc
      call bmc_1_2%intersect_distance(vertices, p, intersec_dist)
      dist = distance(p%tau_travel, ksca)
      pathlen = min(intersec_dist, dist)
      call update_photon_loc(p, pathlen, ksca)
      p%loc(1:2) = old_loc(1:2)
    else
      call bmc%intersect_distance(vertices, p, intersec_dist)
      dist = distance(p%tau_travel, ksca)
      pathlen = min(intersec_dist, dist)
      call update_photon_loc(p, pathlen, ksca)
    end if

    if (intersec_dist .le. dist) then
      lexit_cell = .true.
    else
      lexit_cell = .false.
    end if
  end subroutine move_photon

  subroutine check_if_photon_is_in_domain(C, p)
    type(t_coord), intent(in) :: C
    type(t_photon), intent(in) :: p
    if (p%k .lt. C%zs .or. p%k .gt. C%ze) &
      call CHKERR(1_mpiint, 'Wrong index(dim1) '//toStr(p%k)//' not in ('//toStr(C%zs)//'/'//toStr(C%ze)//')')
    if (p%i .lt. C%xs .or. p%i .gt. C%xe) &
      call CHKERR(1_mpiint, 'Wrong index(dim2) '//toStr(p%i)//' not in ('//toStr(C%xs)//'/'//toStr(C%xe)//')')
    if (p%j .lt. C%ys .or. p%j .gt. C%ye) &
      call CHKERR(1_mpiint, 'Wrong index(dim3) '//toStr(p%j)//' not in ('//toStr(C%ys)//'/'//toStr(C%ye)//')')
  end subroutine
end module
