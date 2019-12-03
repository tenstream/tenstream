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

module m_mcrts_dmda

#include "petsc/finclude/petsc.h"
  use petsc
  use m_tenstream_options, only: mcrts_photons_per_pixel

  use m_data_parameters, only : ireals, iintegers, ireal_dp,     &
    init_mpi_data_parameters, mpiint, imp_iinteger,              &
    zero, one, nil, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi, &
    default_str_len

  use m_helper_functions, only: CHKERR, spherical_2_cartesian, &
    ind_nd_to_1d, ind_1d_to_nd, ndarray_offsets, imp_allreduce_sum, &
    get_arg, itoa, ftoa, cstr

  use m_helper_functions_dp, only : deg2rad

  use m_boxmc, only: t_photon, print_photon, scatter_photon, roulette, R, &
    tau, distance, update_photon_loc, &
    t_boxmc, t_boxmc_1_2, t_boxmc_3_6, &
    imp_t_photon

  use m_boxmc_geometry, only: setup_default_unit_cube_geometry

  use m_pprts_base, only: t_solver, t_solver_1_2, t_solver_3_6, &
    t_state_container, t_coord

  use m_petsc_helpers, only: getVecPointer, restoreVecPointer

  implicit none

  integer(mpiint), parameter :: &
    PQ_SELF     =  1, &
    PQ_NORTH    =  2, &
    PQ_EAST     =  3, &
    PQ_SOUTH    =  4, &
    PQ_WEST     =  5, &
    PQ_NULL     = 10, &
    PQ_READY_TO_RUN  = 12, &
    PQ_RUNNING  = 13, &
    PQ_SENDING  = 14

  type :: t_distributed_photon
    type(t_photon) :: p
    integer(mpiint) :: pstatus, request
  end type

  type :: t_photon_queue
    type(t_distributed_photon), allocatable :: photons(:)
    integer(mpiint) :: current, current_size
    integer(mpiint) :: emptyslots ! count of empty slots in this queue
    integer(mpiint) :: owner ! owner is the owning rank, i.e. myid or the neighbor id
    integer(mpiint) :: queue_index ! is the STATUS integer, i.e. one of PQ_SELF, PQ_NORTH etc.
  end type

  !logical, parameter :: ldebug=.True.
  logical, parameter :: ldebug=.False.

  real(ireal_dp), parameter :: loceps = zero !sqrt(epsilon(loceps))
  integer(iintegers), parameter :: E_up=0, E_dn=1

contains
  subroutine solve_mcrts(solver, edirTOA, solution)
    class(t_solver),intent(in) :: solver
    real(ireals),intent(in) :: edirTOA
    type(t_state_container) :: solution

    integer(mpiint) :: myid, numnodes, ierr

    type(t_photon_queue) :: pqueues(5) ! [own, north, east, south, west]

    integer(iintegers) :: Nphotons, photon_limit
    integer(iintegers) :: ip, kp, iter
    integer(iintegers) :: started_photons, globally_started_photons
    integer(iintegers) :: killed_photons, globally_killed_photons
    integer(mpiint) :: killed_request, stat(mpi_status_size)
    logical :: lcomm_finished


    class(t_boxmc), allocatable :: bmc

    globally_started_photons = solver%C_one%glob_xm * solver%C_one%glob_ym * mcrts_photons_per_pixel
    Nphotons = solver%C_one%xm * solver%C_one%ym * mcrts_photons_per_pixel

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    if(ldebug) then
      print *,myid,'Edir TOA', edirTOA, ':', Nphotons
      print *,myid,'Domain start:', solver%C_one%zs , solver%C_one%xs,  solver%C_one%ys
      print *,myid,'Domain end  :', solver%C_one%ze , solver%C_one%xe,  solver%C_one%ye
      print *,myid,'Domain Size :', solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym
      print *,myid,'Global Size :', solver%C_one%glob_zm , solver%C_one%glob_xm,  solver%C_one%glob_ym

      print *, myid, 'my neighs NESW', &
        solver%C_one%neighbors(22), solver%C_one%neighbors(16), solver%C_one%neighbors( 4), solver%C_one%neighbors(10)
    endif

    !if(numnodes.ne.1) call CHKERR(1_mpiint, 'cannot run mcrts in parallel, would need to implement message passing and for that, need global unique indices')

    call setup_photon_queue(pqueues(PQ_SELF), Nphotons, myid, PQ_SELF)
    call setup_photon_queue(pqueues(PQ_NORTH),Nphotons*10, solver%C_one%neighbors(22), PQ_NORTH)
    call setup_photon_queue(pqueues(PQ_EAST), Nphotons*10, solver%C_one%neighbors(16), PQ_EAST)
    call setup_photon_queue(pqueues(PQ_SOUTH),Nphotons*10, solver%C_one%neighbors( 4), PQ_SOUTH)
    call setup_photon_queue(pqueues(PQ_WEST), Nphotons*10, solver%C_one%neighbors(10), PQ_WEST)

    if(solver%C_one%neighbors(22).lt.zero) call CHKERR(solver%C_one%neighbors(22), 'Bad neighbor id for PQ_NORTH')
    if(solver%C_one%neighbors(16).lt.zero) call CHKERR(solver%C_one%neighbors(16), 'Bad neighbor id for PQ_EAST')
    if(solver%C_one%neighbors( 4).lt.zero) call CHKERR(solver%C_one%neighbors( 4), 'Bad neighbor id for PQ_SOUTH')
    if(solver%C_one%neighbors(10).lt.zero) call CHKERR(solver%C_one%neighbors(10), 'Bad neighbor id for PQ_WEST')

    call VecSet(solution%edir, zero,ierr); call CHKERR(ierr)
    call VecSet(solution%ediff,zero,ierr); call CHKERR(ierr)
    call VecSet(solution%abso, zero,ierr); call CHKERR(ierr)

    select type (solver)
    class is (t_solver_1_2)
      allocate(t_boxmc_1_2::bmc)

    class is (t_solver_3_6)
      allocate(t_boxmc_3_6::bmc)

    class default
      call CHKERR(1_mpiint, 'initialize bmc for mcrts: unexpected type for solver object for DMDA computations!'// &
      '-- call with -solver 1_2 or -solver 3_6')

  end select

  if(ldebug) then
    call bmc%init(MPI_COMM_WORLD, rngseed=9, luse_random_seed=.False.)
  else
    call bmc%init(MPI_COMM_WORLD, rngseed=9, luse_random_seed=.False.)
  endif

  ! Initialize the locally owned photons
  call prepare_locally_owned_photons(solver, bmc, pqueues(PQ_SELF), mcrts_photons_per_pixel)

  !pqueues(PQ_SELF)%photons(:)%p%weight = edirTOA*solver%atm%dx*solver%atm%dy/mcrts_photons_per_pixel
  pqueues(PQ_SELF)%photons(:)%p%weight = edirTOA/real(mcrts_photons_per_pixel, ireals)

  killed_photons = 0
  call mpi_iallreduce(killed_photons, globally_killed_photons, 1_mpiint, imp_iinteger, &
    MPI_SUM, solver%comm, killed_request, ierr); call CHKERR(ierr)

  iter = 0
  do
    iter = iter + 1

    ! Start local photons in batches
    photon_limit = size(pqueues(PQ_NORTH)%photons, kind=iintegers)/100_iintegers
    photon_limit = max(mcrts_photons_per_pixel, photon_limit)

    call run_photon_queue(solver, bmc, solution, pqueues, PQ_SELF, &
      started_photons=ip, killed_photons=kp, limit_number_photons=photon_limit)
    started_photons = ip; killed_photons = killed_photons + kp

    remote_photons: do ! Run remote photons until they are done, only then go on, this hopefully keeps the queue small
      started_photons = 0
      call exchange_photons(solver, pqueues)

      call run_photon_queue(solver, bmc, solution, pqueues, PQ_NORTH, started_photons=ip, killed_photons=kp)
      started_photons = started_photons + ip; killed_photons = killed_photons + kp

      call run_photon_queue(solver, bmc, solution, pqueues, PQ_EAST,  started_photons=ip, killed_photons=kp)
      started_photons = started_photons + ip; killed_photons = killed_photons + kp

      call run_photon_queue(solver, bmc, solution, pqueues, PQ_SOUTH, started_photons=ip, killed_photons=kp)
      started_photons = started_photons + ip; killed_photons = killed_photons + kp

      call run_photon_queue(solver, bmc, solution, pqueues, PQ_WEST,  started_photons=ip, killed_photons=kp)
      started_photons = started_photons + ip; killed_photons = killed_photons + kp

      if(started_photons.eq.0) exit remote_photons
    enddo remote_photons

    call mpi_test(killed_request, lcomm_finished, stat, ierr); call CHKERR(ierr)
    if(lcomm_finished) then
      if(myid.eq.0) print *, iter, &
        'Globally killed photons', globally_killed_photons, '/', globally_started_photons, &
        '('//ftoa(100 * real(globally_killed_photons)/real(globally_started_photons))//' % )'
      if(ldebug) then
        print *, myid, 'Globally killed photons', globally_killed_photons,'/',globally_started_photons
        call print_pqueue(pqueues(PQ_NORTH))
        call print_pqueue(pqueues(PQ_EAST))
        call print_pqueue(pqueues(PQ_SOUTH))
        call print_pqueue(pqueues(PQ_WEST))
      endif
      if(globally_killed_photons.eq.globally_started_photons) then
        exit
      endif
      ! if we reach here this means there is still work todo, setup a new allreduce
      call mpi_iallreduce(killed_photons, globally_killed_photons, 1_mpiint, imp_iinteger, &
        MPI_SUM, solver%comm, killed_request, ierr); call CHKERR(ierr)
    endif

  enddo

  !call VecScale(solution%edir, one/Nphotons, ierr); call CHKERR(ierr)
  !call VecScale(solution%ediff,one/Nphotons, ierr); call CHKERR(ierr)
  !call VecScale(solution%abso, one/Nphotons, ierr); call CHKERR(ierr)

  call PetscObjectSetName(solution%edir,'edir',ierr) ; call CHKERR(ierr)
  call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, '-mcrts_show_edir', ierr); call CHKERR(ierr)
  call PetscObjectViewFromOptions(solution%ediff,PETSC_NULL_VEC, '-mcrts_show_ediff', ierr); call CHKERR(ierr)
  call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, '-mcrts_show_abso', ierr); call CHKERR(ierr)
end subroutine

subroutine run_photon_queue(solver, bmc, solution, pqueues, ipq, started_photons, killed_photons, limit_number_photons)
  class(t_solver),intent(in) :: solver
  class(t_boxmc), intent(in) :: bmc
  type(t_state_container), intent(in) :: solution
  type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
  integer(mpiint), intent(in) :: ipq
  integer(iintegers), intent(out) :: started_photons
  integer(iintegers), intent(out) :: killed_photons
  integer(iintegers) :: iphoton, Nphotmax
  integer(iintegers), optional :: limit_number_photons
  real(ireals),pointer,dimension(:,:,:,:) :: xv_dir=>null(),xv_diff=>null(),xv_abso=>null()
  real(ireals),pointer,dimension(:) :: xv_dir1d=>null(),xv_diff1d=>null(),xv_abso1d=>null()
  logical :: lkilled_photon

  Nphotmax = get_arg(huge(Nphotmax), limit_number_photons)
  Nphotmax = min(size(pqueues(ipq)%photons, kind=iintegers), Nphotmax)

  call getVecPointer(solution%edir , solver%C_dir%da , xv_dir1d , xv_dir)
  call getVecPointer(solution%ediff, solver%C_diff%da, xv_diff1d, xv_diff)
  call getVecPointer(solution%abso , solver%C_one%da , xv_abso1d, xv_abso)

  killed_photons = 0
  started_photons = 0
  do iphoton=1,size(pqueues(ipq)%photons)
    if(pqueues(ipq)%photons(iphoton)%pstatus .eq. PQ_READY_TO_RUN) then
      call run_photon(solver, bmc, pqueues, ipq, iphoton, xv_dir, xv_diff, xv_abso, lkilled_photon)
      if(lkilled_photon) killed_photons = killed_photons + 1
      started_photons = started_photons + 1
      if(started_photons.ge.Nphotmax) exit
    endif
  enddo

  call restoreVecPointer(solution%edir , xv_dir1d , xv_dir)
  call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff)
  call restoreVecPointer(solution%abso , xv_abso1d, xv_abso)
end subroutine

subroutine run_photon(solver, bmc, pqueues, ipq, iphoton, xv_dir, xv_diff, xv_abso, lkilled_photon)
  class(t_solver),intent(in) :: solver
  class(t_boxmc), intent(in) :: bmc
  type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
  integer(mpiint), intent(in) :: ipq
  integer(iintegers), intent(in) :: iphoton
  real(ireals),pointer,dimension(:,:,:,:),intent(in) :: xv_dir, xv_diff, xv_abso
  logical, intent(out) :: lkilled_photon

  !integer(iintegers) :: k,i,j
  real(ireal_dp) :: kabs, ksca, g, mu, phi
  real(ireals), allocatable :: vertices(:)
  logical :: lexit_cell
  integer(mpiint) :: myid, ierr

  lkilled_photon=.False.

  call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

  if(pqueues(ipq)%photons(iphoton)%pstatus.ne.PQ_READY_TO_RUN) &
    call CHKERR(1_mpiint, 'We cannot run photons that are not in READY_TO_RUN state')

  call pqueue_set_status(pqueues(ipq), iphoton, PQ_RUNNING)

  associate(p => pqueues(ipq)%photons(iphoton)%p)

    if(ldebug) print *,myid,cstr('Start of run_photon :: QUEUE:','pink'), ipq, 'iphoton', iphoton
    if(ldebug) call print_photon(p)
    call check_if_photon_is_in_domain(solver%C_one, p)

    if(p%src_side.eq.i1) then ! started at the top of the box, lets increment TOA downward flux
      if(p%direct) then
        p%side = 1
        call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
      else
        call CHKERR(1_mpiint, 'shouldnt happen?')
      endif
    endif

    move: do ! this loop will move the photon to the edge of the subdomain
      call roulette(p)
      if(.not. p%alive) then
        lkilled_photon=.True.
        exit move
      endif

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

      if(.not.lexit_cell) then
        g = solver%atm%g(p%k,p%i,p%j)
        call scatter_photon(p, g)
        p%tau_travel = tau(R())
        if(ldebug) print *,myid,cstr('********************* SCATTERING','peach'),p%k,p%i,p%j
      else ! lexit_cell
        ! Determine actions on boundaries
        select case(p%side)
        case(1)
          if(p%k.eq.solver%C_one%zs) then ! outgoing at TOA
            call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
            if(ldebug) print *,myid,'********************* Exit TOA',p%k,p%i,p%j
            lkilled_photon=.True.
            exit move
          endif
        case(2)
          if(p%k.eq.solver%C_one%ze) then ! hit the surface, need reflection
            call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
            if(ldebug) print *,myid,'********************* Before Reflection',p%k,p%i,p%j

            p%weight = p%weight * solver%atm%albedo(p%i,p%j)
            p%direct=.False.
            p%scattercnt = p%scattercnt +1

            mu = sqrt(R())
            phi = deg2rad( R()*360 )
            p%dir = (/sin(phi)*sin(acos(mu)) , cos(phi)*sin(acos(mu)) , mu  /)
            p%loc(3) = zero + loceps

            p%side = i1
            p%src_side = i2
            call update_flx(p, p%k+1, p%i, p%j, xv_dir, xv_diff)
            if(ldebug) print *,myid,cstr('********************* After  Reflection','aqua'),p%k,p%i,p%j
            cycle move
          endif
        end select

        ! Move photon to new box
        if(ldebug) print *,myid,cstr('********************* MOVE Photon','green')
        select case(p%side)
        case(1)
          p%loc(3) = zero + loceps
          call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
          p%k = p%k-1
          p%src_side = 2
        case(2)
          p%loc(3) = solver%atm%dz(p%k+1, p%i, p%j) - loceps
          call update_flx(p, p%k, p%i, p%j, xv_dir, xv_diff)
          p%k = p%k+1
          p%src_side = 1
        case(3)
          p%loc(1) = solver%atm%dx - loceps
          p%i = p%i-1
          p%src_side = 4
          if(p%i.eq.solver%C_one%xs-1) then
            if(ldebug) print *,myid,cstr('********************* Sending to WEST','blue'), pqueues(PQ_WEST)%owner, p%k,p%i,p%j
            call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_WEST))
            exit move
          endif
        case(4)
          p%loc(1) = zero + loceps
          p%i = p%i+1
          p%src_side = 3
          if(p%i.eq.solver%C_one%xe+1) then
            if(ldebug) print *,myid,cstr('********************* Sending to EAST','blue'), pqueues(PQ_EAST)%owner, p%k,p%i,p%j
            call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_EAST))
            exit move
          endif
        case(5)
          p%loc(2) = solver%atm%dy - loceps
          p%j = p%j-1
          p%src_side = 6
          if(p%j.eq.solver%C_one%ys-1) then
            if(ldebug) print *,myid,cstr('********************* Sending to SOUTH','blue'), pqueues(PQ_SOUTH)%owner, p%k,p%i,p%j
            call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_SOUTH))
            exit move
          endif
        case(6)
          p%loc(2) = zero + loceps
          p%j = p%j+1
          p%src_side = 5
          if(p%j.eq.solver%C_one%ye+1) then
            if(ldebug) print *,myid,cstr('********************* Sending to NORTH','blue'), pqueues(PQ_NORTH)%owner, p%k,p%i,p%j
            call send_photon_to_neighbor(solver, solver%C_one, p, pqueues(PQ_NORTH))
            exit move
          endif
        end select

      endif
    enddo move
    call pqueue_set_status(pqueues(ipq), iphoton, PQ_NULL)
  end associate
end subroutine

subroutine update_flx(p, k, i, j, xv_dir, xv_diff)
  type(t_photon), intent(in) :: p
  integer(iintegers), intent(in) :: k, i, j ! layer/box indices
  real(ireals),pointer,dimension(:,:,:,:) :: xv_dir,xv_diff

  if(ldebug) print *,'Update Flux', k, i, j, p%direct, p%side

  select case(p%side)
  case (1)
    if(k.lt.lbound(xv_diff,2) .or. k.gt.ubound(xv_diff,2)) call CHKERR(1_mpiint, 'invalid k '//itoa(k)//' '//&
      itoa(lbound(xv_diff,2))//'/'//itoa(ubound(xv_diff,2)))
  case (2)
    if(k.lt.lbound(xv_diff,2) .or. k.gt.ubound(xv_diff,2)-1) call CHKERR(1_mpiint, 'invalid k '//itoa(k)//' '//&
      itoa(lbound(xv_diff,2))//'/'//itoa(ubound(xv_diff,2)-1))
  case (3:6)
    continue
  case default
    call CHKERR(1_mpiint, 'hmpf .. didnt expect a p%side gt 6'//itoa(p%side))
  end select

  if(i.lt.lbound(xv_diff,3) .or. i.gt.ubound(xv_diff,3)) call CHKERR(1_mpiint, 'invalid i '//itoa(i)//' '//&
    itoa(lbound(xv_diff,3))//'/'//itoa(ubound(xv_diff,3)))
  if(j.lt.lbound(xv_diff,4) .or. j.gt.ubound(xv_diff,4)) call CHKERR(1_mpiint, 'invalid j '//itoa(j)//' '//&
    itoa(lbound(xv_diff,4))//'/'//itoa(ubound(xv_diff,4)))

  if(p%direct) then
    select case(p%side)
    case (1)
      xv_dir(i0, k, i, j) = xv_dir(i0, k, i, j) + real(p%weight, ireals)
    case (2)
      xv_dir(i0, k+1, i, j) = xv_dir(i0, k+1, i, j) + real(p%weight, ireals)
    case default
      !call print_photon(p)
      call CHKERR(1_mpiint, 'hmpf .. didnt expect a p%side gt 2'//itoa(p%side))
    end select
  else
    select case(p%side)
    case (1) ! Eup
      xv_diff(E_up, k, i, j) = xv_diff(E_up, k, i, j) + real(p%weight, ireals)
    case (2) ! Edn
      xv_diff(E_dn, k+1, i, j) = xv_diff(E_dn, k+1, i, j) + real(p%weight, ireals)
    case default
      !call print_photon(p)
      call CHKERR(1_mpiint, 'hmpf .. didnt expect a p%side gt 2'//itoa(p%side))
    end select
  endif
end subroutine

subroutine prepare_locally_owned_photons(solver, bmc, pqueue, Nphotons_per_pixel)
  class(t_solver),intent(in) :: solver
  class(t_boxmc), intent(in) :: bmc
  type(t_photon_queue), intent(inout) :: pqueue
  integer(iintegers),intent(in) :: Nphotons_per_pixel

  type(t_photon) :: p
  real(ireals) :: phi0, theta0
  integer(mpiint) :: ip, myid, ierr
  integer(iintegers) :: i, j, l
  real(ireals), allocatable :: vertices(:)
  real(ireals) :: initial_dir(3)

  call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

  do i = solver%C_one%xs, solver%C_one%xe
    do j = solver%C_one%ys, solver%C_one%ye
      call setup_default_unit_cube_geometry(solver%atm%dx, solver%atm%dy, &
        solver%atm%dz(i0, i, j), vertices)

      phi0   = solver%sun%phi  (i0, i, j)
      theta0 = solver%sun%theta(i0, i, j)
      initial_dir = spherical_2_cartesian(phi0, theta0)

      do l = 1,Nphotons_per_pixel
        call bmc%init_dir_photon(p, i1, .True., real(initial_dir, ireal_dp), real(vertices, ireal_dp), ierr)
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

        call pqueue_add_photon(pqueue, p, PQ_READY_TO_RUN, 0_mpiint, ip)
        pqueue%photons(ip)%p%cellid = myid*100 + ip

        if(ldebug) then
          print *,'Prepared Photon', i, j, ': iphoton', ip
          call print_photon(pqueue%photons(ip)%p)
        endif
      enddo
    enddo
  enddo
end subroutine

subroutine antialiased_photon_start(Nmax, ip, x, y, tilt_angle)
  integer(iintegers), intent(in) :: Nmax, ip ! total number of photons per pixel, and the ip'th point for that
  real(ireal_dp), intent(in), optional :: tilt_angle ! angle by which the grid is rotated in [rad], default: 30deg
  real(ireal_dp), intent(out) :: x, y ! gives x and y position for a regularly sampled tilted grid

  real(ireal_dp) :: tilt_grid, rot_x, rot_y
  integer(iintegers) :: Nx, Ny ! number of pixels in x and y direction
  integer(iintegers) :: i, j

  tilt_grid = get_arg(deg2rad(26.6_ireal_dp), tilt_angle)

  if(Nmax.eq.i1) then
    x = .5_ireals
    y = .5_ireals
    return
  endif

  Nx = int(sqrt(real(Nmax, ireal_dp)), iintegers)
  Ny = Nx

  if(ip.gt.Nx*Ny) then
    x = R()
    y = R()
    return
  endif

  j = (ip-1) / Nx ! zero based offsets
  i = (ip-1) - j*Nx

  x = real(i+1, ireals)/real(Nx+1, ireals) * sqrt(5._ireals)/2._ireals
  y = real(j+1, ireals)/real(Ny+1, ireals) * sqrt(5._ireals)/2._ireals

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
  integer(mpiint) :: i
  if(pqueue%emptyslots.eq.0) then
    pqueue%current_size = pqueue%current_size +1
    pqueue%emptyslots = 1
  endif
  if(pqueue%current_size.gt.size(pqueue%photons)) call CHKERR(1_mpiint, 'cant enlarge pqueue%currentsize, reached max')
  ierr = 0
  do i=pqueue%current+1, pqueue%current_size
    if(pqueue%photons(i)%pstatus.eq.PQ_NULL) then
      emptyid = i
      return
    endif
  enddo
  do i = 1, pqueue%current
    if(pqueue%photons(i)%pstatus.eq.PQ_NULL) then
      emptyid = i
      return
    endif
  enddo
  emptyid = -1
  ierr = 1
end subroutine

subroutine print_pqueue(pq)
  type(t_photon_queue), intent(in) :: pq
  !integer(iintegers) :: i
  print *,'PQUEUE:', pq%owner, '::', pq%queue_index
  print *,'READY:  ', count(pq%photons(:)%pstatus.eq.PQ_READY_TO_RUN)
  print *,'RUNNING:', count(pq%photons(:)%pstatus.eq.PQ_RUNNING)
  print *,'SENDING:', count(pq%photons(:)%pstatus.eq.PQ_SENDING)
  print *,'PQ_NULL:', count(pq%photons(:)%pstatus.eq.PQ_NULL)
  !do i=1,size(pq%photons)
  !  select case(pq%photons(i)%pstatus)
  !  case (PQ_NULL)
  !    !print *,i,':',pq%photons(i)%pstatus
  !  case (PQ_READY_TO_RUN )
  !    print *,i,': READY_TO_RUN'
  !  case (PQ_RUNNING)
  !    print *,i,': RUNNING'
  !  case (PQ_SENDING )
  !    print *,i,': SENDING'
  !  end select
  !enddo
end subroutine

subroutine pqueue_add_photon(pqueue, p, pstatus, request, ind)
  type(t_photon_queue), intent(inout) :: pqueue
  type(t_photon), intent(in) :: p
  integer(mpiint), intent(in) :: pstatus, request
  integer(mpiint), intent(out) :: ind
  integer(mpiint) :: ierr

  call find_empty_entry_in_pqueue(pqueue, ind, ierr)
  if(ierr.ne.0) then
    call finalize_msgs(pqueue, lwait=.True.)
    call find_empty_entry_in_pqueue(pqueue, ind, ierr)
    call CHKERR(ierr, 'Could not find an empty slot in neighbor queue '//itoa(pqueue%queue_index))
  endif

  pqueue%current = ind
  ! Put photon to send queue to neighbor
  pqueue%photons(ind)%p = p
  pqueue%photons(ind)%pstatus = pstatus
  pqueue%photons(ind)%request = request
  pqueue%emptyslots = pqueue%emptyslots -1
end subroutine

subroutine pqueue_set_status(pqueue, ind, pstatus)
  type(t_photon_queue), intent(inout) :: pqueue
  integer(iintegers), intent(in) :: ind
  integer(mpiint), intent(in) :: pstatus
  if(pqueue%photons(ind)%pstatus.eq.PQ_NULL) pqueue%emptyslots = pqueue%emptyslots -1
  pqueue%photons(ind)%pstatus = pstatus
  if(pstatus.eq.PQ_NULL) pqueue%emptyslots = pqueue%emptyslots +1
end subroutine

subroutine setup_photon_queue(pq, N, owner, queue_index)
  type(t_photon_queue), intent(inout) :: pq
  integer(iintegers), intent(in) :: N
  integer(mpiint), intent(in) :: owner, queue_index

  if(allocated(pq%photons)) call CHKERR(1_mpiint, 'photon queue already allocated')
  allocate(pq%photons(N))
  pq%photons(:)%pstatus = PQ_NULL
  pq%current = 0
  pq%owner = owner
  pq%queue_index = queue_index
  pq%current_size = int(N, mpiint)
  pq%emptyslots = pq%current_size
end subroutine

subroutine send_photon_to_neighbor(solver, C, p_in, pqueue)
  class(t_solver),intent(in) :: solver
  type(t_coord), intent(in) :: C
  type(t_photon),intent(in) :: p_in
  type(t_photon_queue), intent(inout) :: pqueue
  integer(mpiint) :: myid, tag, ierr, iphoton
  !integer(iintegers) :: i,j,k

  logical, parameter :: lcyclic_boundary = .True.

  call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

  ! Put photon to send queue to neighbor
  call pqueue_add_photon(pqueue, p_in, PQ_SENDING, 0_mpiint, iphoton)
  associate(p=>pqueue%photons(iphoton)%p)

    select case(pqueue%queue_index)
    case (PQ_SELF)
      call CHKERR(1_mpiint, 'should not happen that a photon is sent to ourselves?')
      print *, C%xm
    case (PQ_WEST)
      !print *,'Sending Photon WEST'
      if(p%side.ne.i3) call CHKERR(1_mpiint, 'I would assume that side should be 3 when sending WEST')
      if(p%src_side.ne.i4) call CHKERR(1_mpiint, 'I would assume that src_side should be 4 when sending WEST')
    case (PQ_EAST)
      !print *,'Sending Photon EAST'
      if(p%side.ne.i4) call CHKERR(1_mpiint, 'I would assume that side should be 4 when sending EAST')
      if(p%src_side.ne.i3) call CHKERR(1_mpiint, 'I would assume that src_side should be 3 when sending EAST')
    case (PQ_NORTH)
      !print *,'Sending Photon NORTH'
      if(p%side.ne.i6) call CHKERR(1_mpiint, 'I would assume that side should be 6 when sending NORTH')
      if(p%src_side.ne.i5) call CHKERR(1_mpiint, 'I would assume that src_side should be 5 when sending NORTH')
    case (PQ_SOUTH)
      !print *,'Sending Photon SOUTH'
      if(p%side.ne.i5) call CHKERR(1_mpiint, 'I would assume that side should be 5 when sending SOUTH')
      if(p%src_side.ne.i6) call CHKERR(1_mpiint, 'I would assume that src_side should be 6 when sending SOUTH')
    end select

    if(lcyclic_boundary) then
      if(p%i.eq.-1) p%i = C%glob_xm-1
      if(p%j.eq.-1) p%j = C%glob_ym-1
      if(p%i.eq.C%glob_xm) p%i = 0
      if(p%j.eq.C%glob_ym) p%j = 0
    else
      call CHKERR(1_mpiint, 'NON-cyclic boundary conditions are not implemented yet!')
    endif

    ! asynchronous SEND starts here
    tag = pqueue%queue_index
    call mpi_isend(p, 1_mpiint, imp_t_photon, &
      pqueue%owner, tag, solver%comm, pqueue%photons(iphoton)%request, ierr); call CHKERR(ierr, 'mpi isend failed')
    !call mpi_send(p, 1_mpiint, imp_t_photon, &
    !  pqueue%owner, tag, solver%comm, ierr); call CHKERR(ierr, 'mpi isend failed')

    call pqueue_set_status(pqueue, int(iphoton, iintegers), PQ_SENDING)

    if(ldebug) then
      print *,'Sending the following photon to', pqueue%owner,':tag', tag
      call print_photon(p)
    endif

  end associate
end subroutine

subroutine finalize_msgs(pqueue, lwait)
  type(t_photon_queue), intent(inout) :: pqueue
  logical, intent(in) :: lwait

  integer(iintegers) :: i
  integer(mpiint) :: iter, stat(mpi_status_size), ierr
  integer(mpiint) :: cnt_finished_msgs
  logical :: lcomm_finished

  cnt_finished_msgs = 0
  do iter = 1, 1000
    ! See if we can close messages already
    do i = 1, size(pqueue%photons)
      if(pqueue%photons(i)%pstatus .eq. PQ_SENDING) then
        call mpi_test(pqueue%photons(i)%request, lcomm_finished, stat, ierr); call CHKERR(ierr)
        if(lcomm_finished) then
          call pqueue_set_status(pqueue, i, PQ_NULL)
          cnt_finished_msgs = cnt_finished_msgs +1
          if(ldebug) then
            print *,'Finalized Sending a photon:', i
            call print_photon(pqueue%photons(i)%p)
          endif
        endif
      endif
    enddo
    if(.not.lwait) return
    if(cnt_finished_msgs.ne.0) return
  enddo
  call CHKERR(1_mpiint, 'waited too long but nothing happened... exiting...')
end subroutine

subroutine exchange_photons(solver, pqueues)
  class(t_solver),intent(in) :: solver
  type(t_photon_queue), intent(inout) :: pqueues(:) ! [own, north, east, south, west]
  integer(mpiint) :: myid, ipq, mpi_status(mpi_status_size), ierr, tag
  logical :: lgot_msg

  call mpi_comm_rank(solver%comm, myid, ierr)

  do ipq = 1, size(pqueues)
    call finalize_msgs(pqueues(ipq), lwait=.False.)
  enddo

  ! receive all messages that we can get
  do
    call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, solver%comm, lgot_msg, mpi_status, ierr); call CHKERR(ierr)
    if(lgot_msg) then
      tag = mpi_status(MPI_TAG)
      select case(tag)
      case (PQ_WEST) ! was sent towards WEST, i.e. it arrives EAST
        ipq = PQ_EAST
      case (PQ_EAST)
        ipq = PQ_WEST
      case (PQ_SOUTH)
        ipq = PQ_NORTH
      case (PQ_NORTH)
        ipq = PQ_SOUTH
      end select
      if(mpi_status(MPI_SOURCE).ne.pqueues(ipq)%owner) call CHKERR(1_mpiint, 'Something unexpected happened')

      call find_empty_entry_in_pqueue(pqueues(ipq), pqueues(ipq)%current, ierr)
      if(ierr.ne.0) then
        call finalize_msgs(pqueues(ipq), lwait=.True.)
        call find_empty_entry_in_pqueue(pqueues(ipq), pqueues(ipq)%current, ierr)
        call CHKERR(ierr, 'no space in queue to receive a msg')
      endif

      call mpi_recv(pqueues(ipq)%photons(pqueues(ipq)%current)%p, 1_mpiint, imp_t_photon, &
        mpi_status(MPI_SOURCE), mpi_status(MPI_TAG), solver%comm, mpi_status, ierr); call CHKERR(ierr)

      call pqueue_set_status(pqueues(ipq), int(pqueues(ipq)%current, iintegers), PQ_READY_TO_RUN)

      if(ldebug) then
        print *,myid, 'Got Message from:', mpi_status(MPI_SOURCE), 'Receiving ipq', ipq
        call print_photon(pqueues(ipq)%photons(pqueues(ipq)%current)%p)
        !call check_if_cellid_is_in_domain(solver%C_one, pqueues(ipq)%photons(pqueues(ipq)%current)%p%cellid)
      endif
    else
      exit
    endif
  enddo

end subroutine

subroutine move_photon(bmc, vertices, kabs, ksca, p, lexit_cell)
  class(t_boxmc) :: bmc
  real(ireal_dp), intent(in) :: vertices(:), kabs, ksca
  type(t_photon),intent(inout) :: p
  logical,intent(out) :: lexit_cell

  real(ireal_dp) :: dist,intersec_dist

  call bmc%intersect_distance(vertices, p, intersec_dist)

  dist = distance(p%tau_travel, ksca)

  if(intersec_dist .le. dist) then
    call update_photon_loc(p, intersec_dist, kabs, ksca)
    lexit_cell = .True.
  else
    call update_photon_loc(p, dist, kabs, ksca)
    lexit_cell = .False.
  endif
end subroutine move_photon

!function cellindex_from_local_kij(C, k, i, j, lcheck_local_bounds) result(cellid)
!  type(t_coord), intent(in) :: C
!  integer(iintegers), intent(in) :: k, i, j
!  logical, intent(in), optional :: lcheck_local_bounds
!  logical :: lcheck_bounds
!  integer(iintegers) :: cellid
!  integer(iintegers) :: dmda_offsets(3)
!
!  lcheck_bounds = get_arg(.True., lcheck_local_bounds)
!  if(lcheck_bounds) then
!    if(k.le.C%zs .or. k.gt.C%ze+1) call CHKERR(1_mpiint, 'Wrong index in dimension 1 : '//itoa(k)//' not in ('//itoa(C%zs)//'/'//itoa(C%ze)//')')
!    if(i.le.C%xs .or. i.gt.C%xe+1) call CHKERR(1_mpiint, 'Wrong index in dimension 2 : '//itoa(k)//' not in ('//itoa(C%xs)//'/'//itoa(C%xe)//')')
!    if(j.le.C%ys .or. j.gt.C%ye+1) call CHKERR(1_mpiint, 'Wrong index in dimension 3 : '//itoa(k)//' not in ('//itoa(C%ys)//'/'//itoa(C%ye)//')')
!  endif
!
!  call ndarray_offsets([C%glob_zm, C%glob_xm, C%glob_ym], dmda_offsets)
!  cellid = ind_nd_to_1d(dmda_offsets, [k, i, j])
!  if(ldebug) print *,'local kij',k,i,j,'=>', cellid
!end function
!
!subroutine kij_from_cellid(C, cellid, k, i, j)
!  type(t_coord), intent(in) :: C
!  integer(iintegers), intent(in) :: cellid
!  integer(iintegers), intent(out) :: k, i, j
!  integer(iintegers) :: dmda_offsets(3)
!  integer(iintegers) :: cellindices3d(3)
!
!  call ndarray_offsets([C%glob_zm, C%glob_xm, C%glob_ym], dmda_offsets)
!  call ind_1d_to_nd(dmda_offsets, cellid, cellindices3d)
!  k = cellindices3d(1)
!  i = cellindices3d(2)
!  j = cellindices3d(3)
!
!  if(ldebug) then
!    call check_if_cellid_is_in_domain(C, cellid)
!  endif
!end subroutine
!
subroutine check_if_photon_is_in_domain(C, p)
  type(t_coord), intent(in) :: C
  type(t_photon),intent(in) :: p
  if(p%k.lt.C%zs .or. p%k.gt.C%ze) &
    call CHKERR(1_mpiint, 'Wrong index(dim1) '//itoa(p%k)//' not in ('//itoa(C%zs)//'/'//itoa(C%ze)//')')
  if(p%i.lt.C%xs .or. p%i.gt.C%xe) &
    call CHKERR(1_mpiint, 'Wrong index(dim2) '//itoa(p%i)//' not in ('//itoa(C%xs)//'/'//itoa(C%xe)//')')
  if(p%j.lt.C%ys .or. p%j.gt.C%ye) &
    call CHKERR(1_mpiint, 'Wrong index(dim3) '//itoa(p%j)//' not in ('//itoa(C%ys)//'/'//itoa(C%ye)//')')
end subroutine
end module
