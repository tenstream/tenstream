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
  use iso_fortran_env, only: int64

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, ireal_dp, &
                               mpiint, imp_iinteger, imp_int8, &
                               zero, i0, i1, i2, i3, i4, i5, i6

  use m_helper_functions, only: &
    & CHKERR, &
    & cstr, &
    & deg2rad, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_allreduce_sum, &
    & ind_1d_to_nd, &
    & ind_nd_to_1d, &
    & ndarray_offsets, &
    & spherical_2_cartesian, &
    & toStr

  use m_boxmc, only: t_photon, print_photon, scatter_photon, roulette, R, &
                     tau, distance, update_photon_loc, &
                     t_boxmc, t_boxmc_1_2, t_boxmc_3_6, &
                     imp_t_photon

  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use m_pprts_base, only: t_solver, t_solver_1_2, t_solver_3_6, t_solver_mcdmda, &
                          t_state_container, t_coord

  use m_petsc_helpers, only: getVecPointer, restoreVecPointer

  use m_buildings, only: &
    & t_pprts_buildings!, &
  !& t_plex_buildings, &
  !& clone_buildings, &
  !& destroy_buildings, &
  !& check_buildings_consistency, &
  !& PPRTS_TOP_FACE, &
  !& PPRTS_BOT_FACE, &
  !& PPRTS_LEFT_FACE, &
  !& PPRTS_REAR_FACE

  use m_linked_list_mpiint, only: t_list_mpiint, t_node

  implicit none

  ! queue status indices
  integer(mpiint), parameter :: &
    PQ_SELF = 1, &
    PQ_NORTH = 2, &
    PQ_EAST = 3, &
    PQ_SOUTH = 4, &
    PQ_WEST = 5

  ! states a photon can have
  integer(mpiint), parameter :: &
    PQ_NULL = 6, &
    PQ_READY_TO_RUN = 7, &
    PQ_RUNNING = 8, &
    PQ_SENDING = 9

  character(len=16), parameter :: id2name(9) = [ &
    & character(len=16) :: &
    & 'PQ_SELF', &
    & 'PQ_NORTH', &
    & 'PQ_EAST', &
    & 'PQ_SOUTH', &
    & 'PQ_WEST', &
    & 'PQ_NULL', &
    & 'PQ_READY_TO_RUN', &
    & 'PQ_RUNNING', &
    & 'PQ_SENDING' &
    & ]

  type :: t_distributed_photon
    type(t_photon) :: p
    integer(mpiint) :: request
  end type

  type :: t_photon_queue
    type(t_distributed_photon), allocatable :: photons(:)
    integer(mpiint) :: current, current_size
    type(t_list_mpiint) :: readyslots ! linked list for read_to_go photon indices
    type(t_list_mpiint) :: emptyslots ! linked_list of empty slots in this queue
    type(t_list_mpiint) :: sendingslots ! linked_list of sending slots in this queue
    integer(mpiint) :: owner ! owner is the owning rank, i.e. myid or the neighbor id
    integer(mpiint) :: queue_index ! is the STATUS integer, i.e. one of PQ_SELF, PQ_NORTH etc.
  end type

  logical, parameter :: ldebug = .false.

  real(ireal_dp), parameter :: loceps = 0 !sqrt(epsilon(loceps))
  integer(iintegers), parameter :: E_up = 0, E_dn = 1

contains

  subroutine determine_Nphotons(solver, Nphotons_local, ierr)
    class(t_solver), intent(in) :: solver
    integer(iintegers), intent(out) :: Nphotons_local
    integer(mpiint), intent(out) :: ierr
    integer(iintegers) :: mcdmda_photons_per_pixel ! has to be constant over all TOA faces
    real(ireals) :: rN
    logical :: lflg
    integer(mpiint) :: numnodes

    ierr = 0

    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    mcdmda_photons_per_pixel = 1000
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-mcdmda_photons_per_px", &
                       mcdmda_photons_per_pixel, lflg, ierr); call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-mcdmda_photons", rN, lflg, ierr); call CHKERR(ierr)
    if (lflg) then
      mcdmda_photons_per_pixel = int(rN, kind(Nphotons_local)) / (numnodes * solver%C_one%xm * solver%C_one%ym)
      mcdmda_photons_per_pixel = max(1_iintegers, mcdmda_photons_per_pixel)
    end if

    Nphotons_local = solver%C_one%xm * solver%C_one%ym * mcdmda_photons_per_pixel

  end subroutine

  subroutine solve_mcdmda(solver, edirTOA, solution, ierr, opt_buildings)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container) :: solution
    integer(mpiint), intent(out) :: ierr
    type(t_pprts_buildings), intent(in), optional :: opt_buildings

    integer(mpiint) :: myid

    type(t_photon_queue) :: pqueues(5) ! [own, north, east, south, west]

    integer(iintegers), parameter :: Nqueuesize=100000
    integer(iintegers) :: Nphotons_global, globally_killed_photons
    integer(iintegers) :: Nphotons_local, photon_limit
    integer(iintegers) :: started_photons
    integer(iintegers) :: killed_photons
    integer(iintegers) :: ip, kp
    integer(iintegers) :: iter
    integer(mpiint) :: killed_request, stat(mpi_status_size)
    logical :: lcomm_finished

    class(t_boxmc), allocatable :: bmc

    ierr = 0

    call determine_Nphotons(solver, Nphotons_local, ierr); call CHKERR(ierr)
    call mpi_allreduce(Nphotons_local, Nphotons_global, 1_mpiint, imp_iinteger, &
                       MPI_SUM, solver%comm, ierr); call CHKERR(ierr)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    if (ldebug) then
      print *, myid, 'Edir TOA', edirTOA, ': Nphotons_local', Nphotons_local, 'Nphotons_global', Nphotons_global
      print *, myid, 'Domain start:', solver%C_one%zs, solver%C_one%xs, solver%C_one%ys
      print *, myid, 'Domain end  :', solver%C_one%ze, solver%C_one%xe, solver%C_one%ye
      print *, myid, 'Domain Size :', solver%C_one%zm, solver%C_one%xm, solver%C_one%ym
      print *, myid, 'Global Size :', solver%C_one%glob_zm, solver%C_one%glob_xm, solver%C_one%glob_ym

      print *, myid, 'my neighs NESW', &
        solver%C_one%neighbors(22), solver%C_one%neighbors(16), solver%C_one%neighbors(4), solver%C_one%neighbors(10)
    end if

    call setup_photon_queue(pqueues(PQ_SELF) , Nphotons_local, myid, PQ_SELF)
    call setup_photon_queue(pqueues(PQ_NORTH), Nqueuesize, solver%C_one%neighbors(22), PQ_NORTH)
    call setup_photon_queue(pqueues(PQ_EAST) , Nqueuesize, solver%C_one%neighbors(16), PQ_EAST)
    call setup_photon_queue(pqueues(PQ_SOUTH), Nqueuesize, solver%C_one%neighbors(4), PQ_SOUTH)
    call setup_photon_queue(pqueues(PQ_WEST) , Nqueuesize, solver%C_one%neighbors(10), PQ_WEST)

    if (solver%C_one%neighbors(22) .lt. zero) call CHKERR(solver%C_one%neighbors(22), 'Bad neighbor id for PQ_NORTH')
    if (solver%C_one%neighbors(16) .lt. zero) call CHKERR(solver%C_one%neighbors(16), 'Bad neighbor id for PQ_EAST')
    if (solver%C_one%neighbors(4) .lt. zero) call CHKERR(solver%C_one%neighbors(4), 'Bad neighbor id for PQ_SOUTH')
    if (solver%C_one%neighbors(10) .lt. zero) call CHKERR(solver%C_one%neighbors(10), 'Bad neighbor id for PQ_WEST')

    call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
    call VecSet(solution%abso, zero, ierr); call CHKERR(ierr)

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

    if (ldebug) then
      call bmc%init(MPI_COMM_WORLD, rngseed=9, luse_random_seed=.false.)
    else
      call bmc%init(MPI_COMM_WORLD, rngseed=9, luse_random_seed=.false.)
    end if

    ! Initialize the locally owned photons
    call prepare_locally_owned_photons(solver, bmc, pqueues(PQ_SELF), Nphotons_local)

    pqueues(PQ_SELF)%photons(:)%p%weight = edirTOA * solver%C_one%xm * solver%C_one%ym / real(Nphotons_local, ireals)

    killed_photons = 0
    call mpi_iallreduce(killed_photons, globally_killed_photons, 1_mpiint, imp_iinteger, &
                        MPI_SUM, solver%comm, killed_request, ierr); call CHKERR(ierr)

    iter = 0
    do
      iter = iter + 1

      ! Start local photons in batches
      photon_limit = max(1_iintegers, size(pqueues(PQ_NORTH)%photons, kind=iintegers) / 2_iintegers)

      call run_photon_queue(solver, bmc, solution, pqueues, PQ_SELF, &
                            started_photons=ip, killed_photons=kp, limit_number_photons=photon_limit)
      started_photons = ip; killed_photons = killed_photons + kp
      if (ldebug) print *, 'SELF', 'started_photons', started_photons, 'killed_photons', killed_photons

      remote_photons: do ! Run remote photons until they are done, only then go on, this hopefully keeps the queue small
        started_photons = 0
        call exchange_photons(solver, pqueues)

        call run_photon_queue(solver, bmc, solution, pqueues, PQ_NORTH, started_photons=ip, killed_photons=kp)
        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'NORTH', 'started_photons', ip, 'killed_photons', kp

        call run_photon_queue(solver, bmc, solution, pqueues, PQ_EAST, started_photons=ip, killed_photons=kp)
        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'EAST', 'started_photons', ip, 'killed_photons', kp

        call run_photon_queue(solver, bmc, solution, pqueues, PQ_SOUTH, started_photons=ip, killed_photons=kp)
        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'SOUTH', 'started_photons', ip, 'killed_photons', kp

        call run_photon_queue(solver, bmc, solution, pqueues, PQ_WEST, started_photons=ip, killed_photons=kp)
        started_photons = started_photons + ip; killed_photons = killed_photons + kp
        if (ldebug) print *, 'WEST', 'started_photons', ip, 'killed_photons', kp

        if (started_photons .eq. 0) then
          if (ldebug) print *, 'exit remote_photons branch'
          exit remote_photons
        end if
      end do remote_photons

      call mpi_test(killed_request, lcomm_finished, stat, ierr); call CHKERR(ierr)
      if (lcomm_finished) then
        if (myid .eq. 0) print *, iter, &
          'Globally killed photons', globally_killed_photons, '/', Nphotons_global, &
          '('//toStr(100 * real(globally_killed_photons) / real(Nphotons_global))//' % )'
        if (ldebug) then
          print *, myid, 'Globally killed photons', globally_killed_photons, '/', Nphotons_global
          call print_pqueue(pqueues(PQ_NORTH))
          call print_pqueue(pqueues(PQ_EAST))
          call print_pqueue(pqueues(PQ_SOUTH))
          call print_pqueue(pqueues(PQ_WEST))
        end if
        if (globally_killed_photons .eq. Nphotons_global) then
          exit
        end if
        ! if we reach here this means there is still work todo, setup a new allreduce
        call mpi_iallreduce(killed_photons, globally_killed_photons, 1_mpiint, imp_iinteger, &
                            MPI_SUM, solver%comm, killed_request, ierr); call CHKERR(ierr)
      end if

    end do

    !call VecScale(solution%edir, one/Nphotons, ierr); call CHKERR(ierr)
    !call VecScale(solution%ediff,one/Nphotons, ierr); call CHKERR(ierr)
    !call VecScale(solution%abso, one/Nphotons, ierr); call CHKERR(ierr)

    call PetscObjectSetName(solution%edir, 'edir', ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, '-mcdmda_show_edir', ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, '-mcdmda_show_ediff', ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-mcdmda_show_abso', ierr); call CHKERR(ierr)

    ! cleanup
    call bmc%destroy(ierr); call CHKERR(ierr)
    do ip = 1, size(pqueues)
      call photon_queue_destroy(pqueues(ip), ierr); call CHKERR(ierr)
    enddo
  end subroutine

  subroutine photon_queue_destroy(pq, ierr)
    type(t_photon_queue), intent(inout) :: pq
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    if(allocated(pq%photons)) deallocate(pq%photons)
    call pq%readyslots%finalize()
    call pq%emptyslots%finalize()
    call pq%sendingslots%finalize()
  end subroutine

  subroutine run_photon_queue(solver, bmc, solution, pqueues, ipq, started_photons, killed_photons, limit_number_photons)
    class(t_solver), intent(in) :: solver
    class(t_boxmc), intent(in) :: bmc
    type(t_state_container), intent(in) :: solution
    type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
    integer(mpiint), intent(in) :: ipq
    integer(iintegers), intent(out) :: started_photons
    integer(iintegers), intent(out) :: killed_photons
    integer(iintegers) :: Nphotmax
    integer(iintegers), optional :: limit_number_photons
    real(ireals), pointer, dimension(:, :, :, :) :: xv_dir => null(), xv_diff => null(), xv_abso => null()
    real(ireals), pointer, dimension(:) :: xv_dir1d => null(), xv_diff1d => null(), xv_abso1d => null()
    logical :: lkilled_photon

    if (ldebug) then
      print *, 'run photon queue ', id2name(pqueues(ipq)%queue_index)
      call print_pqueue(pqueues(ipq))
    end if

    Nphotmax = get_arg(huge(Nphotmax), limit_number_photons)
    Nphotmax = min(size(pqueues(ipq)%photons, kind=iintegers), Nphotmax)

    call getVecPointer(solver%C_dir%da, solution%edir, xv_dir1d, xv_dir)
    call getVecPointer(solver%C_diff%da, solution%ediff, xv_diff1d, xv_diff)
    call getVecPointer(solver%C_one%da, solution%abso, xv_abso1d, xv_abso)

    killed_photons = 0
    started_photons = 0
    call pqueues(ipq)%readyslots%for_each(run_ready_slots)

    call restoreVecPointer(solver%C_dir%da, solution%edir, xv_dir1d, xv_dir)
    call restoreVecPointer(solver%C_diff%da, solution%ediff, xv_diff1d, xv_diff)
    call restoreVecPointer(solver%C_one%da, solution%abso, xv_abso1d, xv_abso)
    contains
      subroutine run_ready_slots(idx, node, iphoton)
        integer, intent(in) :: idx
        type(t_node), pointer, intent(inout) :: node
        integer(mpiint), intent(inout) :: iphoton
        integer(mpiint) :: ierr

        if (started_photons .ge. Nphotmax) return

        call run_photon(&
          & solver, &
          & bmc, &
          & pqueues, &
          & ipq, &
          & iphoton, &
          & xv_dir, &
          & xv_diff, &
          & xv_abso, &
          & lkilled_photon)

        started_photons = started_photons + 1
        if (lkilled_photon) killed_photons = killed_photons + 1

        call pqueues(ipq)%readyslots%del_node(node, ierr); call CHKERR(ierr)
        return
        print *,'remove unused var warning', idx
      end subroutine
  end subroutine

  subroutine run_photon(solver, bmc, pqueues, ipq, iphoton, xv_dir, xv_diff, xv_abso, lkilled_photon)
    class(t_solver), intent(in) :: solver
    class(t_boxmc), intent(in) :: bmc
    type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
    integer(mpiint), intent(in) :: ipq
    integer(iintegers), intent(in) :: iphoton
    real(ireals), pointer, dimension(:, :, :, :), intent(in) :: xv_dir, xv_diff, xv_abso
    logical, intent(out) :: lkilled_photon

    real(ireal_dp) :: kabs, ksca, g, mu, phi
    real(ireals), allocatable :: vertices(:)
    logical :: lexit_cell
    integer(mpiint) :: myid, ierr

    lkilled_photon = .false.

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    associate (p => pqueues(ipq)%photons(iphoton)%p)

      if (ldebug) print *, myid, cstr('Start of run_photon :: QUEUE:', 'pink'), id2name(ipq), 'iphoton', iphoton
      if (ldebug) call print_photon(p)
      call check_if_photon_is_in_domain(solver%C_one, p)

      if (p%src_side .eq. i1) then ! started at the top of the box, lets increment TOA downward flux
        if (p%direct) then
          p%side = 1
          call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
        else
          call CHKERR(1_mpiint, 'shouldnt happen?')
        end if
      end if

      move: do ! this loop will move the photon to the edge of the subdomain
        call roulette(p)
        if (.not. p%alive) then
          lkilled_photon = .true.
          exit move
        end if

        call check_if_photon_is_in_domain(solver%C_one, p)

        kabs = solver%atm%kabs(p%k, p%i, p%j)
        ksca = solver%atm%ksca(p%k, p%i, p%j)

        call setup_default_unit_cube_geometry(solver%atm%dx, solver%atm%dy, &
                                              solver%atm%dz(p%k, p%i, p%j), vertices)

        xv_abso(i0, p%k, p%i, p%j) = xv_abso(i0, p%k, p%i, p%j) + real(p%weight, ireals)

        call move_photon(bmc, real(vertices, ireal_dp), kabs, ksca, p, lexit_cell)

        xv_abso(i0, p%k, p%i, p%j) = xv_abso(i0, p%k, p%i, p%j) - real(p%weight, ireals)

        !print *,'start of move', k, i, j
        !call print_photon(p)

        ! Define new cellindex and update local position

        if (.not. lexit_cell) then
          g = solver%atm%g(p%k, p%i, p%j)
          call scatter_photon(p, g)
          p%tau_travel = tau(R())
          if (ldebug) print *, myid, cstr('********************* SCATTERING', 'peach'), p%k, p%i, p%j
        else ! lexit_cell
          ! Determine actions on boundaries
          select case (p%side)
          case (1)
            if (p%k .eq. solver%C_one%zs) then ! outgoing at TOA
              call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
              if (ldebug) print *, myid, '********************* Exit TOA', p%k, p%i, p%j
              lkilled_photon = .true.
              exit move
            end if
          case (2)
            if (p%k .eq. solver%C_one%ze) then ! hit the surface, need reflection
              call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
              if (ldebug) print *, myid, '********************* Before Reflection', p%k, p%i, p%j

              p%weight = p%weight * solver%atm%albedo(p%i, p%j)
              p%direct = .false.
              p%scattercnt = p%scattercnt + 1

              mu = sqrt(R())
              phi = deg2rad(R() * 360)
              p%dir = (/sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu/)
              p%loc(3) = zero + loceps

              p%side = i1
              p%src_side = i2
              call update_flx(p, p%k + 1, p%i, p%j, xv_dir, xv_diff)
              if (ldebug) print *, myid, cstr('********************* After  Reflection', 'aqua'), p%k, p%i, p%j
              cycle move
            end if
          end select

          ! Move photon to new box
          if (ldebug) print *, myid, cstr('* MOVE Photon', 'green')
          select case (p%side)
          case (1)
            p%loc(3) = zero + loceps
            call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
            p%k = p%k - 1
            p%src_side = 2
          case (2)
            p%loc(3) = solver%atm%dz(p%k + 1, p%i, p%j) - loceps
            call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
            p%k = p%k + 1
            p%src_side = 1
          case (3)
            p%loc(1) = solver%atm%dx - loceps
            p%i = p%i - 1
            p%src_side = 4
            if (p%i .eq. solver%C_one%xs - 1) then
              if (ldebug) print *, myid, cstr('* Sending to WEST', 'blue'), pqueues(PQ_WEST)%owner, p%k, p%i, p%j
              call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_WEST))
              exit move
            end if
          case (4)
            p%loc(1) = zero + loceps
            p%i = p%i + 1
            p%src_side = 3
            if (p%i .eq. solver%C_one%xe + 1) then
              if (ldebug) print *, myid, cstr('* Sending to EAST', 'blue'), pqueues(PQ_EAST)%owner, p%k, p%i, p%j
              call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_EAST))
              exit move
            end if
          case (5)
            p%loc(2) = solver%atm%dy - loceps
            p%j = p%j - 1
            p%src_side = 6
            if (p%j .eq. solver%C_one%ys - 1) then
              if (ldebug) print *, myid, cstr('* Sending to SOUTH', 'blue'), pqueues(PQ_SOUTH)%owner, p%k, p%i, p%j
              call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_SOUTH))
              exit move
            end if
          case (6)
            p%loc(2) = zero + loceps
            p%j = p%j + 1
            p%src_side = 5
            if (p%j .eq. solver%C_one%ye + 1) then
              if (ldebug) print *, myid, cstr('* Sending to NORTH', 'blue'), pqueues(PQ_NORTH)%owner, p%k, p%i, p%j
              call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_NORTH))
              exit move
            end if
          end select

        end if
      end do move
      call pqueues(ipq)%emptyslots%add(iphoton)
    end associate
  end subroutine

  subroutine update_flx(p, k, i, j, xv_dir, xv_diff)
    type(t_photon), intent(in) :: p
    integer(iintegers), intent(in) :: k, i, j ! layer/box indices
    real(ireals), pointer, dimension(:, :, :, :) :: xv_dir, xv_diff

    if (ldebug) print *, 'Update Flux', k, i, j, p%direct, p%side

    select case (p%side)
    case (1)
      if (k .lt. lbound(xv_diff, 2) .or. k .gt. ubound(xv_diff, 2)) &
        & call CHKERR(1_mpiint, 'invalid k '//toStr(k)//' '// &
                       & toStr(lbound(xv_diff, 2))//'/'//toStr(ubound(xv_diff, 2)))
    case (2)
      if (k .lt. lbound(xv_diff, 2) .or. k .gt. ubound(xv_diff, 2) - 1) &
        & call CHKERR(1_mpiint, 'invalid k '//toStr(k)//' '// &
                      & toStr(lbound(xv_diff, 2))//'/'//toStr(ubound(xv_diff, 2) - 1))
    case (3:6)
      continue
    case default
      call CHKERR(1_mpiint, 'hmpf .. didnt expect a p%side gt 6'//toStr(p%side))
    end select

    if (i .lt. lbound(xv_diff, 3) .or. i .gt. ubound(xv_diff, 3)) &
      & call CHKERR(1_mpiint, 'invalid i '//toStr(i)//' '// &
                    & toStr(lbound(xv_diff, 3))//'/'//toStr(ubound(xv_diff, 3)))
    if (j .lt. lbound(xv_diff, 4) .or. j .gt. ubound(xv_diff, 4)) &
      & call CHKERR(1_mpiint, 'invalid j '//toStr(j)//' '// &
                    & toStr(lbound(xv_diff, 4))//'/'//toStr(ubound(xv_diff, 4)))

    if (p%direct) then
      select case (p%side)
      case (1)
        xv_dir(i0, k, i, j) = xv_dir(i0, k, i, j) + real(p%weight, ireals)
      case (2)
        xv_dir(i0, k + 1, i, j) = xv_dir(i0, k + 1, i, j) + real(p%weight, ireals)
      case default
        !call print_photon(p)
        call CHKERR(1_mpiint, 'hmpf .. didnt expect a p%side gt 2'//toStr(p%side))
      end select
    else
      select case (p%side)
      case (1) ! Eup
        xv_diff(E_up, k, i, j) = xv_diff(E_up, k, i, j) + real(p%weight, ireals)
      case (2) ! Edn
        xv_diff(E_dn, k + 1, i, j) = xv_diff(E_dn, k + 1, i, j) + real(p%weight, ireals)
      case default
        !call print_photon(p)
        call CHKERR(1_mpiint, 'hmpf .. didnt expect a p%side gt 2'//toStr(p%side))
      end select
    end if
  end subroutine

  subroutine prepare_locally_owned_photons(solver, bmc, pqueue, Nphotons)
    class(t_solver), intent(in) :: solver
    class(t_boxmc), intent(in) :: bmc
    type(t_photon_queue), intent(inout) :: pqueue
    integer(iintegers), intent(in) :: Nphotons

    type(t_photon) :: p
    real(ireals) :: phi0, theta0
    integer(mpiint) :: ip, myid, ierr
    integer(iintegers) :: i, j
    integer(iintegers) :: Nphotons_per_pixel, iphoton, l
    real(ireals), allocatable :: vertices(:)
    real(ireals) :: initial_dir(3)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    phi0 = solver%sun%phi
    theta0 = solver%sun%theta
    initial_dir = spherical_2_cartesian(phi0, theta0)

    Nphotons_per_pixel = max(1_iintegers, 1_iintegers + Nphotons / int(solver%C_one%xm * solver%C_one%ym, kind(Nphotons)))

    iphoton = 0
    do l = 1, Nphotons_per_pixel
      do i = solver%C_one%xs, solver%C_one%xe
        do j = solver%C_one%ys, solver%C_one%ye
          iphoton = iphoton + 1
          !print *, 'Preparing Photon', i, j, ': iphoton', iphoton, Nphotons
          !call print_pqueue(pqueue)
          if (iphoton .gt. Nphotons) then
            return
          end if

          call setup_default_unit_cube_geometry(solver%atm%dx, solver%atm%dy, &
                                                solver%atm%dz(i0, i, j), vertices)

          call bmc%init_dir_photon(p, i1, .true., real(initial_dir, ireal_dp), real(vertices, ireal_dp), ierr)
          p%i = i
          p%j = j
          p%k = i0
          p%src_side = i1
          p%tau_travel = tau(R())
          !p%loc(1) = solver%atm%dx/2
          !p%loc(2) = solver%atm%dy/2
          call antialiased_photon_start(Nphotons_per_pixel, l, p%loc(1), p%loc(2))
          p%loc(1) = p%loc(1) * solver%atm%dx
          p%loc(2) = p%loc(2) * solver%atm%dy

          call pqueue_add_photon(pqueue, p, 0_mpiint, ip)
          call pqueue%readyslots%add(ip)

          pqueue%photons(ip)%p%cellid = myid * 100 + ip

          if (ldebug) then
            print *, 'Prepared Photon', i, j, ': iphoton', ip, iphoton, Nphotons
            call print_photon(pqueue%photons(ip)%p)
          end if
        end do
      end do
    end do
  end subroutine

  subroutine antialiased_photon_start(Nmax, ip, x, y, tilt_angle)
    integer(iintegers), intent(in) :: Nmax, ip ! total number of photons per pixel, and the ip'th point for that
    real(ireal_dp), intent(in), optional :: tilt_angle ! angle by which the grid is rotated in [rad], default: 30deg
    real(ireal_dp), intent(out) :: x, y ! gives x and y position for a regularly sampled tilted grid

    real(ireal_dp) :: tilt_grid, rot_x, rot_y
    integer(iintegers) :: Nx, Ny ! number of pixels in x and y direction
    integer(iintegers) :: i, j

    tilt_grid = get_arg(deg2rad(26.6_ireal_dp), tilt_angle)

    if (Nmax .eq. i1) then
      x = .5_ireals
      y = .5_ireals
      return
    end if

    Nx = int(sqrt(real(Nmax, ireal_dp)), iintegers)
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
    integer(mpiint), intent(out) :: emptyid
    integer(mpiint), intent(out) :: ierr
    call pqueue%emptyslots%pop(emptyid, ierr)
  end subroutine

  subroutine print_pqueue(pq)
    type(t_photon_queue), intent(in) :: pq
    print *, 'PQUEUE:', pq%owner, '::', pq%queue_index, 'name ', cstr(trim(id2name(pq%queue_index)), 'blue')
    print *, 'current', pq%current, 'current_size', pq%current_size
    print *, 'readyslots:'
    call pq%readyslots%view()
    print *, 'emptyslots:'
    call pq%emptyslots%view()
    print *, 'sendingslots:'
    call pq%sendingslots%view()
  end subroutine

  subroutine pqueue_add_photon(pqueue, p, request, ind)
    type(t_photon_queue), intent(inout) :: pqueue
    type(t_photon), intent(in) :: p
    integer(mpiint), intent(in) :: request
    integer(mpiint), intent(out) :: ind
    integer(mpiint) :: ierr

    call find_empty_entry_in_pqueue(pqueue, ind, ierr)
    if (ierr .ne. 0) then
      call finalize_msgs_blocking(pqueue)
      call find_empty_entry_in_pqueue(pqueue, ind, ierr)
      call CHKERR(ierr, 'Could not find an empty slot in neighbor queue '//toStr(pqueue%queue_index))
    end if

    pqueue%current = ind
    ! Put photon to send queue to neighbor
    pqueue%photons(ind)%p = p
    pqueue%photons(ind)%request = request
  end subroutine

  subroutine setup_photon_queue(pq, N, owner, queue_index)
    type(t_photon_queue), intent(inout) :: pq
    integer(iintegers), intent(in) :: N
    integer(mpiint), intent(in) :: owner, queue_index
    integer(mpiint) :: i

    if (allocated(pq%photons)) call CHKERR(1_mpiint, 'photon queue already allocated')
    allocate (pq%photons(N))
    pq%current = 0
    pq%owner = owner
    pq%queue_index = queue_index
    pq%current_size = int(N, mpiint)
    do i = 1, pq%current_size
      call pq%emptyslots%add(i)
    enddo
    !pq%readyslots = 0
  end subroutine

  subroutine send_photon_to_neighbor(solver, C, p_in, pqueue)
    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    type(t_photon), intent(in) :: p_in
    type(t_photon_queue), intent(inout) :: pqueue
    integer(mpiint) :: myid, tag, ierr, iphoton

    logical, parameter :: lcyclic_boundary = .true.

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    ! Put photon to send queue to neighbor
    call pqueue_add_photon(pqueue, p_in, 0_mpiint, iphoton)
    call pqueue%sendingslots%add(iphoton)

    associate (p => pqueue%photons(iphoton)%p)

      if(ldebug) then
        select case (pqueue%queue_index)
        case (PQ_SELF)
          call CHKERR(1_mpiint, 'should not happen that a photon is sent to ourselves?')
        case (PQ_WEST)
          !print *,'Sending Photon WEST'
          if (p%side .ne. i3) call CHKERR(1_mpiint, 'I would assume that side should be 3 when sending WEST')
          if (p%src_side .ne. i4) call CHKERR(1_mpiint, 'I would assume that src_side should be 4 when sending WEST')
        case (PQ_EAST)
          !print *,'Sending Photon EAST'
          if (p%side .ne. i4) call CHKERR(1_mpiint, 'I would assume that side should be 4 when sending EAST')
          if (p%src_side .ne. i3) call CHKERR(1_mpiint, 'I would assume that src_side should be 3 when sending EAST')
        case (PQ_NORTH)
          !print *,'Sending Photon NORTH'
          if (p%side .ne. i6) call CHKERR(1_mpiint, 'I would assume that side should be 6 when sending NORTH')
          if (p%src_side .ne. i5) call CHKERR(1_mpiint, 'I would assume that src_side should be 5 when sending NORTH')
        case (PQ_SOUTH)
          !print *,'Sending Photon SOUTH'
          if (p%side .ne. i5) call CHKERR(1_mpiint, 'I would assume that side should be 5 when sending SOUTH')
          if (p%src_side .ne. i6) call CHKERR(1_mpiint, 'I would assume that src_side should be 6 when sending SOUTH')
        end select
      endif

      if (lcyclic_boundary) then
        if (p%i .eq. -1) p%i = C%glob_xm - 1
        if (p%j .eq. -1) p%j = C%glob_ym - 1
        if (p%i .eq. C%glob_xm) p%i = 0
        if (p%j .eq. C%glob_ym) p%j = 0
      else
        call CHKERR(1_mpiint, 'NON-cyclic boundary conditions are not implemented yet!')
      end if

      ! asynchronous SEND starts here
      tag = pqueue%queue_index
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

    call pqueue%sendingslots%for_each(check_sending_slots)

  contains

    subroutine check_sending_slots(idx, node, iphoton)
      integer, intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(mpiint), intent(inout) :: iphoton

      integer(mpiint) :: stat(mpi_status_size), ierr
      logical :: lcomm_finished

      call mpi_test(pqueue%photons(iphoton)%request, lcomm_finished, stat, ierr); call CHKERR(ierr)
      if (lcomm_finished) then
        call pqueue%emptyslots%add(iphoton)
        call pqueue%sendingslots%del_node(node, ierr); call CHKERR(ierr)
        if (ldebug) then
          print *, 'Finalized Sending a photon:', iphoton
          call print_photon(pqueue%photons(iphoton)%p)
        end if
      end if
      return
      print *,'unused var warning', idx
    end subroutine
  end subroutine

  subroutine finalize_msgs_blocking(pqueue)
    type(t_photon_queue), intent(inout) :: pqueue

    integer(mpiint) :: iter, cnt_finished_msgs

    cnt_finished_msgs = 0
    do iter = 1, 1000
      call pqueue%sendingslots%for_each(check_sending_slots)
      if (cnt_finished_msgs .ne. 0) return
    end do
    call CHKERR(1_mpiint, 'waited too long but nothing happened... exiting...')
  contains
    subroutine check_sending_slots(idx, node, iphoton)
      integer, intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(mpiint), intent(inout) :: iphoton

      integer(mpiint) :: stat(mpi_status_size), ierr
      logical :: lcomm_finished

      call mpi_test(pqueue%photons(iphoton)%request, lcomm_finished, stat, ierr); call CHKERR(ierr)
      if (lcomm_finished) then
        call pqueue%emptyslots%add(iphoton)
        call pqueue%sendingslots%del_node(node, ierr); call CHKERR(ierr)
        cnt_finished_msgs = cnt_finished_msgs + 1
        if (ldebug) then
          print *, 'Finalized Sending a photon:', iphoton
          call print_photon(pqueue%photons(iphoton)%p)
        end if
      end if
      return
      print *,'unused var warning', idx
    end subroutine
  end subroutine

  subroutine exchange_photons(solver, pqueues)
    class(t_solver), intent(in) :: solver
    type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
    integer(mpiint) :: myid, ipq, mpi_status(mpi_status_size), ierr, tag
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
        end select
        if (mpi_status(MPI_SOURCE) .ne. pqueues(ipq)%owner) call CHKERR(1_mpiint, 'Something unexpected happened')

        call find_empty_entry_in_pqueue(pqueues(ipq), pqueues(ipq)%current, ierr)
        if (ierr .ne. 0) then
          call finalize_msgs_blocking(pqueues(ipq))
          call find_empty_entry_in_pqueue(pqueues(ipq), pqueues(ipq)%current, ierr)
          call CHKERR(ierr, 'no space in queue to receive a msg')
        end if

        call mpi_recv(pqueues(ipq)%photons(pqueues(ipq)%current)%p, 1_mpiint, imp_t_photon, &
                      mpi_status(MPI_SOURCE), mpi_status(MPI_TAG), solver%comm, mpi_status, ierr); call CHKERR(ierr)

        call pqueues(ipq)%readyslots%add(pqueues(ipq)%current)

        if (ldebug) then
          print *, myid, 'Got Message from rank:', mpi_status(MPI_SOURCE), 'Receiving ipq ', id2name(ipq)
          call print_photon(pqueues(ipq)%photons(pqueues(ipq)%current)%p)
          !call check_if_cellid_is_in_domain(solver%C_one, pqueues(ipq)%photons(pqueues(ipq)%current)%p%cellid)
        end if
      else
        exit
      end if
    end do

  end subroutine

  subroutine move_photon(bmc, vertices, kabs, ksca, p, lexit_cell)
    class(t_boxmc) :: bmc
    real(ireal_dp), intent(in) :: vertices(:), kabs, ksca
    type(t_photon), intent(inout) :: p
    logical, intent(out) :: lexit_cell

    real(ireal_dp) :: dist, intersec_dist

    call bmc%intersect_distance(vertices, p, intersec_dist)

    dist = distance(p%tau_travel, ksca)

    if (intersec_dist .le. dist) then
      call update_photon_loc(p, intersec_dist, kabs, ksca)
      lexit_cell = .true.
    else
      call update_photon_loc(p, dist, kabs, ksca)
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
