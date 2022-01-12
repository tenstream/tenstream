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

module m_pprts_explicit

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : &
    & i0, i1, &
    & iintegers, mpiint, &
    & ireals, imp_ireals, &
    & zero, one

  use m_helper_functions, only : &
    & approx, &
    & CHKERR, &
    & imp_min_mean_max, &
    & mpi_logical_and, &
    & toStr

  use m_pprts_base, only : &
    & atmk, &
    & determine_ksp_tolerances, &
    & setup_incSolar, &
    & t_solver, &
    & t_state_container

  use m_petsc_helpers, only: &
    getVecPointer, restoreVecPointer

  implicit none

  private
  public :: &
    & exchange_diffuse_boundary, &
    & exchange_direct_boundary, &
    & explicit_ediff, &
    & explicit_ediff_sor_sweep, &
    & explicit_edir, &
    & explicit_edir_forward_sweep

  logical,parameter :: ldebug=.False.

contains

  !> @brief explicit loop to compute direct radiation
  subroutine explicit_edir(solver, prefix, edirTOA, solution, ierr)
  class(t_solver), target, intent(in) :: solver
    character(len=*), intent(in) :: prefix
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), intent(inout) :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: x0=>null(), xg=>null()
    real(ireals),pointer,dimension(:) :: x01d=>null(), xg1d=>null()
    type(tVec) :: v0, b, lb

    integer(iintegers), dimension(3) :: dx, dy ! start, end, increment for each dimension
    integer(iintegers) :: iter, maxiter

    real(ireals), allocatable :: residual(:)
    real(ireals) :: residual_mmm(3), rel_residual, atol, rtol

    real(ireals) :: ignore_max_it! Ignore max iter setting if time is less

    logical :: lsun_north, lsun_east, lpermute, lskip_residual, lmonitor_residual, lconverged, lflg, lflg2
    logical :: laccept_incomplete_solve, lconverged_reason

    integer(mpiint) :: ierr

    ierr = 0

    call PetscLogEventBegin(solver%logs%compute_Edir, ierr)
    associate( &
        & vedir => solution%edir, &
        & sun => solver%sun, &
        & atm => solver%atm, &
        & C   => solver%C_dir)


      maxiter=1000
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, prefix, &
        "-ksp_max_it", maxiter, lflg, ierr) ;call CHKERR(ierr)
      
      ignore_max_it=huge(ignore_max_it)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, prefix, "-ksp_ignore_max_it", &
        ignore_max_it, lflg , ierr) ;call CHKERR(ierr)
      if (solution%time(1).lt.ignore_max_it) then
         maxiter=1001
      endif
      allocate(residual(maxiter))

      lskip_residual = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
        "-ksp_skip_residual", lskip_residual, lflg , ierr) ;call CHKERR(ierr)

      lmonitor_residual = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
        & "-ksp_monitor", lmonitor_residual, lflg , ierr) ;call CHKERR(ierr)
      if(.not.lmonitor_residual) then
        call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
          & "-ksp_monitor_true_residual", lmonitor_residual, lflg , ierr) ;call CHKERR(ierr)
      endif

      lconverged_reason = lmonitor_residual
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
        & "-ksp_converged_reason", lconverged_reason, lflg , ierr) ;call CHKERR(ierr)

      laccept_incomplete_solve = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)

      call determine_ksp_tolerances(C, atm%l1d, rtol, atol)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, prefix, "-ksp_atol", &
        atol, lflg , ierr) ;call CHKERR(ierr)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, prefix, "-ksp_rtol", &
        rtol, lflg2, ierr) ;call CHKERR(ierr)
      if(lmonitor_residual) then
        if(lflg.or.lflg2) then
          if(solver%myid.eq.0) print *,'pprts_explicit_edir setting rtol/atol', rtol, atol
        endif
      endif

      lsun_north = sun%yinc.eq.i0
      lsun_east  = sun%xinc.eq.i0

      dx = [C%xs, C%xe, i1]
      dy = [C%ys, C%ye, i1]

      lpermute = .True.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
        "-explicit_edir_permute", lpermute, lflg , ierr) ;call CHKERR(ierr)
      if(lpermute) then
        if(lsun_east) dx = [dx(2), dx(1), -dx(3)]
        if(lsun_north) dy = [dy(2), dy(1), -dy(3)]
        if(solver%myid.eq.0.and.ldebug) &
          & print *,'Using permutations for pprts_explicit_edir x-iterator:', dx, 'y-iterator:', dy
      endif

      call DMGetGlobalVector(C%da, b, ierr); call CHKERR(ierr)
      call VecSet(b, zero, ierr); call CHKERR(ierr)
      call setup_incSolar(solver, edirTOA, b)

      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call DMGlobalToLocal(C%da, b, INSERT_VALUES, lb, ierr) ;call CHKERR(ierr)

      call DMGetLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMGlobalToLocal(C%da, vedir, INSERT_VALUES, v0, ierr) ;call CHKERR(ierr)

      do iter = 1, maxiter

        call explicit_edir_forward_sweep(solver, solver%dir2dir, dx, dy, lb, v0)

        call exchange_direct_boundary(solver, lsun_north, lsun_east, v0, ierr); call CHKERR(ierr)

        ! Residual computations
        if(.not.lskip_residual) then
          call getVecPointer(C%da, vedir, xg1d, xg)
          call getVecPointer(C%da, v0, x01d, x0, readonly=.True.)
          residual(iter) = max(tiny(one), norm2(xg - x0(:,:,C%xs:C%xe,C%ys:C%ye)))
          xg = x0(:,:,C%xs:C%xe,C%ys:C%ye)
          call restoreVecPointer(C%da, vedir, xg1d, xg)
          call restoreVecPointer(C%da, v0, x01d, x0, readonly=.True.)
          if(residual(1).le.sqrt(tiny(residual))) then
            rel_residual = 0
          else
            rel_residual = residual(iter)/residual(1)
          endif
          if(solver%myid.eq.0.and.lmonitor_residual) then
            call imp_min_mean_max(solver%comm, residual(iter), residual_mmm)
            residual(iter) = residual_mmm(2)
            if(residual(1).le.sqrt(tiny(residual))) then
              rel_residual = 0
            else
              rel_residual = residual(iter)/residual(1)
            endif
            print *,trim(prefix)//" iter "//toStr(iter)//' residual (min/mean/max)', residual_mmm, &
              & 'rel res', rel_residual
          endif
          solution%dir_ksp_residual_history(min(size(solution%dir_ksp_residual_history, kind=iintegers), iter)) = residual(iter)

          lconverged = mpi_logical_and(solver%comm, residual(iter).lt.atol.or.rel_residual.lt.rtol)
          if(lconverged) then
            if(solver%myid.eq.0.and.lconverged_reason) &
              & print *,trim(prefix)//' solve converged after', iter, 'iterations'
            exit
          endif
        endif

        if(iter.eq.maxiter) then
          if(.not.laccept_incomplete_solve) then
            call CHKERR(int(iter,mpiint), trim(prefix)//" did not converge")
          endif
        endif
      enddo ! iter

      ! update solution vec
      call getVecPointer(C%da, vedir, xg1d, xg)
      call getVecPointer(C%da, v0, x01d, x0, readonly=.True.)
      xg = x0(:,:,C%xs:C%xe,C%ys:C%ye)
      call restoreVecPointer(C%da, vedir, xg1d, xg)
      call restoreVecPointer(C%da, v0, x01d, x0, readonly=.True.)

      call DMRestoreLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call DMRestoreGlobalVector(C%da, b, ierr); call CHKERR(ierr)

      call PetscObjectSetName(vedir,'debug_edir',ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(vedir, PETSC_NULL_VEC, "-show_debug_edir", ierr); call CHKERR(ierr)
    end associate

    solution%lchanged=.True.
    solution%lWm2_dir=.False.
    call PetscLogEventEnd(solver%logs%compute_Edir, ierr)

  end subroutine

  subroutine exchange_direct_boundary(solver, lsun_north, lsun_east, x, ierr)
  class(t_solver), intent(in) :: solver
    logical, intent(in) :: lsun_north, lsun_east
    type(tVec), intent(inout) :: x
    integer(mpiint), intent(out) :: ierr

    real(ireals),pointer :: x0(:,:,:,:)=>null(), x01d(:)=>null()

    integer(mpiint), parameter :: tag_x=1, tag_y=2
    integer(mpiint) :: neigh_s, neigh_r, requests(4), statuses(mpi_status_size,4)
    integer(iintegers) :: dofstart, dofend

    real(ireals), allocatable :: mpi_send_bfr_x(:,:,:), mpi_send_bfr_y(:,:,:)
    real(ireals), allocatable :: mpi_recv_bfr_x(:,:,:), mpi_recv_bfr_y(:,:,:)

    associate( &
        & C   => solver%C_dir)

      allocate(mpi_send_bfr_x(solver%dirside%dof, C%zm, C%ys:C%ye), mpi_recv_bfr_x(solver%dirside%dof, C%zm, C%ys:C%ye))
      allocate(mpi_send_bfr_y(solver%dirside%dof, C%zm, C%xs:C%xe), mpi_recv_bfr_y(solver%dirside%dof, C%zm, C%xs:C%xe))

      call getVecPointer(C%da, x, x01d, x0)

      ! Boundary exchanges
      ! x direction scatters
      dofstart = solver%dirtop%dof
      dofend   = -1 + solver%dirtop%dof + solver%dirside%dof
      if(lsun_east) then
        mpi_send_bfr_x = x0(dofstart:dofend, :, C%xs, C%ys:C%ye)
        neigh_s = int(C%neighbors(10), mpiint) ! neigh west
        neigh_r = int(C%neighbors(16), mpiint) ! neigh east
      else
        mpi_send_bfr_x = x0(dofstart:dofend, :, C%xe+1, C%ys:C%ye)
        neigh_s = int(C%neighbors(16), mpiint) ! neigh east
        neigh_r = int(C%neighbors(10), mpiint) ! neigh west
      endif
      call MPI_Irecv(mpi_recv_bfr_x, size(mpi_recv_bfr_x, kind=mpiint), &
        & imp_ireals, neigh_r, tag_x, solver%comm, requests(1), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_x, size(mpi_send_bfr_x, kind=mpiint), &
        & imp_ireals, neigh_s, tag_x, solver%comm, requests(2), ierr); call CHKERR(ierr)

      ! y direction scatters
      dofstart = solver%dirtop%dof + solver%dirside%dof
      dofend   = -1 + solver%dirtop%dof + solver%dirside%dof * 2
      if(lsun_north) then
        mpi_send_bfr_y = x0(dofstart:dofend, :, C%xs:C%xe, C%ys)
        neigh_s = int(C%neighbors( 4), mpiint) ! neigh south
        neigh_r = int(C%neighbors(22), mpiint) ! neigh north
      else
        mpi_send_bfr_y = x0(dofstart:dofend, :, C%xs:C%xe, C%ye+1)
        neigh_s = int(C%neighbors(22), mpiint) ! neigh north
        neigh_r = int(C%neighbors( 4), mpiint) ! neigh south
      endif
      call MPI_Irecv(mpi_recv_bfr_y, size(mpi_recv_bfr_y, kind=mpiint), &
        & imp_ireals, neigh_r, tag_y, solver%comm, requests(3), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_y, size(mpi_send_bfr_y, kind=mpiint), &
        & imp_ireals, neigh_s, tag_y, solver%comm, requests(4), ierr); call CHKERR(ierr)

      call MPI_Waitall(4_mpiint, requests, statuses, ierr); call CHKERR(ierr)

      dofstart = solver%dirtop%dof
      dofend   = -1 + solver%dirtop%dof + solver%dirside%dof
      if(lsun_east) then
        x0(dofstart:dofend, :, C%xe+1, C%ys:C%ye) = mpi_recv_bfr_x
      else
        x0(dofstart:dofend, :, C%xs, C%ys:C%ye) = mpi_recv_bfr_x
      endif

      dofstart = solver%dirtop%dof + solver%dirside%dof
      dofend   = -1 + solver%dirtop%dof + solver%dirside%dof * 2
      if(lsun_north) then
        x0(dofstart:dofend, :, C%xs:C%xe, C%ye+1) = mpi_recv_bfr_y
      else
        x0(dofstart:dofend, :, C%xs:C%xe, C%ys) = mpi_recv_bfr_y
      endif
      call restoreVecPointer(C%da, x, x01d, x0)
    end associate
    ierr = 0
  end subroutine

  subroutine explicit_edir_forward_sweep(solver, coeffs, dx, dy, b, x)
  class(t_solver), intent(in) :: solver
    real(ireals), target, intent(in) :: coeffs(:,:,:,:)
    integer(iintegers), dimension(3), intent(in) :: dx, dy ! start, end, increment for each dimension
    type(tVec), intent(in) :: b, x

    real(ireals),pointer :: x0(:,:,:,:)=>null(), x01d(:)=>null()
    real(ireals),pointer :: xb(:,:,:,:)=>null(), xb1d(:)=>null()

    integer(iintegers) :: i, j, k
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:,:) ! dim(src, dst)

    associate( &
        & atm => solver%atm, &
        & C   => solver%C_dir, &
        & xinc => solver%sun%xinc, &
        & yinc => solver%sun%yinc )

      call getVecPointer(C%da, x, x01d, x0)
      call getVecPointer(C%da, b, xb1d, xb, readonly=.True.)
      x0(0:solver%dirtop%dof-1, C%zs, :, :) = xb(0:solver%dirtop%dof-1, C%zs, :, :)
      call restoreVecPointer(C%da, b, xb1d, xb, readonly=.True.)

      ! forward sweep through x
      do k=C%zs,C%ze-1
        if( atm%l1d(atmk(atm,k)) ) then
          do j=dy(1), dy(2), dy(3)
            do i=dx(1), dx(2), dx(3)

              do idst = 0, solver%dirtop%dof-1
                x0(idst, k+i1, i, j) = x0(idst, k, i, j) * atm%a33(atmk(atm,k),i,j)
              enddo
            enddo
          enddo
        else
          do j=dy(1), dy(2), dy(3)
            do i=dx(1), dx(2), dx(3)

              v(0:C%dof-1, 0:C%dof-1) => coeffs(:,k-C%zs+1,i-C%xs+1,j-C%ys+1)

              dst = 0
              do idst = 0, solver%dirtop%dof-1
                x0(dst, k+i1, i, j) = 0
                src = 0
                do isrc = 0, solver%dirtop%dof-1
                  x0(dst, k+i1, i, j) = x0(dst, k+i1, i, j) + x0(src, k, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                  x0(dst, k+i1, i, j) = x0(dst, k+i1, i, j) + x0(src, k, i+1-xinc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                  x0(dst, k+i1, i, j) = x0(dst, k+i1, i, j) + x0(src, k, i, j+1-yinc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%dirside%dof-1
                x0(dst, k, i+xinc, j) = 0
                src = 0
                do isrc = 0, solver%dirtop%dof-1
                  x0(dst, k, i+xinc, j) = x0(dst, k, i+xinc, j) + x0(src, k, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                  x0(dst, k, i+xinc, j) = x0(dst, k, i+xinc, j) + x0(src, k, i+1-xinc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                  x0(dst, k, i+xinc, j) = x0(dst, k, i+xinc, j) + x0(src, k, i, j+1-yinc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%dirside%dof-1
                x0(dst, k, i, j+yinc) = 0
                src = 0
                do isrc = 0, solver%dirtop%dof-1
                  x0(dst, k, i, j+yinc) = x0(dst, k, i, j+yinc) + x0(src, k, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                  x0(dst, k, i, j+yinc) = x0(dst, k, i, j+yinc) + x0(src, k, i+1-xinc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                  x0(dst, k, i, j+yinc) = x0(dst, k, i, j+yinc) + x0(src, k, i, j+1-yinc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

            enddo
          enddo
        endif
      enddo
      call restoreVecPointer(C%da, x, x01d, x0)
    end associate
  end subroutine

  !> @brief explicit loop to compute diffuse radiation
  subroutine explicit_ediff(solver, prefix, vb, solution, ierr)
  class(t_solver), intent(inout) :: solver
    character(len=*), intent(in) :: prefix
    type(tVec), intent(in) :: vb
    type(t_state_container), intent(inout) :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: x0=>null(), xg=>null()
    real(ireals),pointer,dimension(:) :: x01d=>null(), xg1d=>null()
    type(tVec) :: lvb, v0

    integer(iintegers) :: iter, isub, maxiter, sub_iter

    real(ireals), allocatable :: residual(:)
    real(ireals) :: residual_mmm(3), rel_residual, atol, rtol, omega

    real(ireals) :: ignore_max_it! Ignore max iter setting if time is less

    logical :: lskip_residual, lmonitor_residual, lconverged, lflg, lflg2
    logical :: laccept_incomplete_solve, lconverged_reason

    integer(mpiint) :: ierr

    ierr = 0

    call PetscLogEventBegin(solver%logs%compute_Ediff, ierr)
    associate( &
        & vediff => solution%ediff, &
        & atm => solver%atm, &
        & C   => solver%C_diff)


      maxiter=1000
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, prefix, &
        "-ksp_max_it", maxiter, lflg, ierr) ;call CHKERR(ierr)

      ignore_max_it=huge(ignore_max_it)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, prefix, "-ksp_ignore_max_it", &
        ignore_max_it, lflg , ierr) ;call CHKERR(ierr)
      if (solution%time(1).lt.ignore_max_it) then
         maxiter=1001
      endif
      allocate(residual(maxiter)) 
      sub_iter = 1
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, prefix, &
        "-sub_it", sub_iter, lflg, ierr) ;call CHKERR(ierr)

      lskip_residual = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
        "-ksp_skip_residual", lskip_residual, lflg , ierr) ;call CHKERR(ierr)

      lmonitor_residual = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
        "-ksp_monitor", lmonitor_residual, lflg , ierr) ;call CHKERR(ierr)
      if(.not.lmonitor_residual) then
        call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
          & "-ksp_monitor_true_residual", lmonitor_residual, lflg , ierr) ;call CHKERR(ierr)
      endif

      lconverged_reason = lmonitor_residual
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
        & "-ksp_converged_reason", lconverged_reason, lflg , ierr) ;call CHKERR(ierr)

      laccept_incomplete_solve = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
        "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)

      call determine_ksp_tolerances(C, atm%l1d, rtol, atol)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, prefix, "-ksp_atol", &
        atol, lflg , ierr) ;call CHKERR(ierr)
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, prefix, "-ksp_rtol", &
        rtol, lflg2, ierr) ;call CHKERR(ierr)
      if(lmonitor_residual) then
        if(lflg.or.lflg2) then
          if(solver%myid.eq.0) print *,'pprts_explicit_ediff setting rtol/atol', rtol, atol
        endif
      endif

      omega = 1
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, prefix, "-pc_sor_omega", &
        omega, lflg , ierr) ;call CHKERR(ierr)


      call DMGetLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(C%da, vediff, INSERT_VALUES, v0, ierr) ;call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C%da, vediff, INSERT_VALUES, v0, ierr) ;call CHKERR(ierr)

      call DMGetLocalVector(C%da, lvb, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(C%da, vb, INSERT_VALUES, lvb, ierr) ;call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C%da, vb, INSERT_VALUES, lvb, ierr) ;call CHKERR(ierr)

      do iter = 1, maxiter

        do isub=1, sub_iter
          if(approx(omega,1._ireals)) then
            call explicit_ediff_forward_sweep(&
              & solver, &
              & solver%diff2diff, &
              & dx = [C%xs, C%xe, i1], &
              & dy = [C%ys, C%ye, i1], &
              & dz = [C%zs, C%ze-1, i1], &
              & b=lvb, x=v0 )
            call explicit_ediff_forward_sweep(&
              & solver, &
              & solver%diff2diff, &
              & dx = [C%xe, C%xs, -i1], &
              & dy = [C%ye, C%ys, -i1], &
              & dz = [C%ze-1, C%zs, -i1], &
              & b=lvb, x=v0 )
          else
            call explicit_ediff_sor_sweep(&
              & solver, &
              & solver%diff2diff, &
              & dx = [C%xs, C%xe, i1], &
              & dy = [C%ys, C%ye, i1], &
              & dz = [C%zs, C%ze-1, i1], &
              & omega = omega, &
              & b=lvb, x=v0 )
            call explicit_ediff_sor_sweep(&
              & solver, &
              & solver%diff2diff, &
              & dx = [C%xe, C%xs, -i1], &
              & dy = [C%ye, C%ys, -i1], &
              & dz = [C%ze-1, C%zs, -i1], &
              & omega = omega, &
              & b=lvb, x=v0 )
          endif
        enddo

        call exchange_diffuse_boundary(solver, v0, ierr); call CHKERR(ierr)

        ! Residual computations
        if(.not.lskip_residual) then
          call getVecPointer(C%da, vediff, xg1d, xg)
          call getVecPointer(C%da, v0, x01d, x0, readonly=.True.)

          residual(iter) = norm2(xg - x0(:,:,C%xs:C%xe,C%ys:C%ye))
          xg = x0(:,:,C%xs:C%xe,C%ys:C%ye)

          call restoreVecPointer(C%da, vediff, xg1d, xg)
          call restoreVecPointer(C%da, v0, x01d, x0, readonly=.True.)

          if(residual(1).le.sqrt(tiny(residual))) then
            rel_residual = 0
          else
            rel_residual = residual(iter)/residual(1)
          endif
          if(solver%myid.eq.0.and.lmonitor_residual) then
            call imp_min_mean_max(solver%comm, residual(iter), residual_mmm)
            residual(iter) = residual_mmm(2)
            if(residual(1).le.sqrt(tiny(residual))) then
              rel_residual = 0
            else
              rel_residual = residual(iter)/residual(1)
            endif
            print *,trim(prefix), ' iter ', toStr(iter), ' residual (min/mean/max)', residual_mmm, &
              & 'rel res', rel_residual
          endif
          solution%diff_ksp_residual_history(min(size(solution%diff_ksp_residual_history, kind=iintegers), iter)) = residual(iter)

          lconverged = mpi_logical_and(solver%comm, residual(iter).lt.atol.or.rel_residual.lt.rtol)
          if(lconverged) then
            if(solver%myid.eq.0.and.lconverged_reason) &
              & print *,trim(prefix),' solve converged after', iter, 'iterations'
            exit
          endif
        endif

        if(iter.eq.maxiter) then
          if(.not.laccept_incomplete_solve) then
            call CHKERR(int(iter,mpiint), trim(prefix)//" did not converge")
          endif
        endif
      enddo ! iter

      call getVecPointer(C%da, vediff, xg1d, xg)
      call getVecPointer(C%da, v0, x01d, x0, readonly=.True.)

      ! update solution vec
      xg = x0(:,:,C%xs:C%xe,C%ys:C%ye)

      call restoreVecPointer(C%da, vediff, xg1d, xg)
      call restoreVecPointer(C%da, v0, x01d, x0, readonly=.True.)

      call DMRestoreLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMRestoreLocalVector(C%da, lvb, ierr); call CHKERR(ierr)

      call PetscObjectSetName(vediff,'debug_ediff',ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(vediff, PETSC_NULL_VEC, "-show_debug_ediff", ierr); call CHKERR(ierr)
    end associate

    solution%lchanged=.True.
    solution%lWm2_diff=.False.
    call PetscLogEventEnd(solver%logs%compute_Ediff, ierr)
  end subroutine

  subroutine exchange_diffuse_boundary(solver, v0, ierr)
  class(t_solver), intent(in) :: solver
    type(tVec), intent(inout) :: v0
    integer(mpiint), intent(out) :: ierr

    integer(mpiint), parameter :: tag_e=1, tag_w=2, tag_n=3, tag_s=4
    integer(mpiint) :: neigh_s, neigh_e, neigh_w, neigh_n
    integer(mpiint) :: requests(8), statuses(mpi_status_size,8)
    integer(iintegers) :: k,i,j,d1,d2,dof,idof

    real(ireals), pointer :: x0(:,:,:,:)=>null(), x01d(:)=>null()

    real(ireals), allocatable :: mpi_send_bfr_e(:,:,:), mpi_send_bfr_n(:,:,:)
    real(ireals), allocatable :: mpi_send_bfr_w(:,:,:), mpi_send_bfr_s(:,:,:)
    real(ireals), allocatable :: mpi_recv_bfr_e(:,:,:), mpi_recv_bfr_n(:,:,:)
    real(ireals), allocatable :: mpi_recv_bfr_w(:,:,:), mpi_recv_bfr_s(:,:,:)


    associate( &
        & C   => solver%C_diff)

      allocate(&
        & mpi_send_bfr_e(solver%diffside%dof/2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_send_bfr_w(solver%diffside%dof/2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_send_bfr_n(solver%diffside%dof/2, C%zs:C%ze, C%xs:C%xe), &
        & mpi_send_bfr_s(solver%diffside%dof/2, C%zs:C%ze, C%xs:C%xe), &
        & mpi_recv_bfr_e(solver%diffside%dof/2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_recv_bfr_w(solver%diffside%dof/2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_recv_bfr_n(solver%diffside%dof/2, C%zs:C%ze, C%xs:C%xe), &
        & mpi_recv_bfr_s(solver%diffside%dof/2, C%zs:C%ze, C%xs:C%xe) &
        & )


      call getVecPointer(C%da, v0, x01d, x0)

      neigh_s = int(C%neighbors( 4), mpiint)
      neigh_w = int(C%neighbors(10), mpiint)
      neigh_e = int(C%neighbors(16), mpiint)
      neigh_n = int(C%neighbors(22), mpiint)

      call MPI_Irecv(mpi_recv_bfr_w, size(mpi_recv_bfr_w, kind=mpiint), &
        & imp_ireals, neigh_w, tag_e, solver%comm, requests(1), ierr); call CHKERR(ierr)
      call MPI_Irecv(mpi_recv_bfr_e, size(mpi_recv_bfr_e, kind=mpiint), &
        & imp_ireals, neigh_e, tag_w, solver%comm, requests(2), ierr); call CHKERR(ierr)

      call MPI_Irecv(mpi_recv_bfr_s, size(mpi_recv_bfr_s, kind=mpiint), &
        & imp_ireals, neigh_s, tag_n, solver%comm, requests(3), ierr); call CHKERR(ierr)
      call MPI_Irecv(mpi_recv_bfr_n, size(mpi_recv_bfr_n, kind=mpiint), &
        & imp_ireals, neigh_n, tag_s, solver%comm, requests(4), ierr); call CHKERR(ierr)

      ! Boundary exchanges
      ! x direction scatters, east/west boundary
      do j=C%ys,C%ye
        do k=C%zs,C%ze
          d1=1; d2=1
          do idof=i0, solver%diffside%dof-1
            dof = solver%difftop%dof + idof
            if (solver%diffside%is_inward(i1+idof)) then ! to the right
              mpi_send_bfr_e(d1,k,j) = x0(dof, k, C%xe+1, j)
              d1 = d1+1
            else !leftward
              mpi_send_bfr_w(d2,k,j) = x0(dof, k, C%xs, j)
              d2 = d2+1
            endif
          enddo
        enddo
      enddo
      call MPI_Isend(mpi_send_bfr_w, size(mpi_send_bfr_w, kind=mpiint), &
        & imp_ireals, neigh_w, tag_w, solver%comm, requests(5), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_e, size(mpi_send_bfr_e, kind=mpiint), &
        & imp_ireals, neigh_e, tag_e, solver%comm, requests(6), ierr); call CHKERR(ierr)

      ! y direction scatters, south/north boundary
      do i=C%xs,C%xe
        do k=C%zs,C%ze
          d1=1; d2=1
          do idof=i0, solver%diffside%dof-1
            dof = solver%difftop%dof + solver%diffside%dof + idof
            if (solver%diffside%is_inward(i1+idof)) then ! forward
              mpi_send_bfr_n(d1,k,i) = x0(dof, k, i, C%ye+1)
              d1 = d1+1
            else ! backward
              mpi_send_bfr_s(d2,k,i) = x0(dof, k, i, C%ys)
              d2 = d2+1
            endif
          enddo
        enddo
      enddo
      call MPI_Isend(mpi_send_bfr_s, size(mpi_send_bfr_s, kind=mpiint), &
        & imp_ireals, neigh_s, tag_s, solver%comm, requests(7), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_n, size(mpi_send_bfr_n, kind=mpiint), &
        & imp_ireals, neigh_n, tag_n, solver%comm, requests(8), ierr); call CHKERR(ierr)

      call MPI_Waitall(8_mpiint, requests, statuses, ierr); call CHKERR(ierr)

      ! receive buffers
      do j=C%ys,C%ye
        do k=C%zs,C%ze
          d1=1; d2=1
          do idof=i0, solver%diffside%dof-1
            dof = solver%difftop%dof + idof
            if (solver%diffside%is_inward(i1+idof)) then ! to the right
              x0(dof, k, C%xs, j) = mpi_recv_bfr_w(d1,k,j)
              d1 = d1+1
            else ! leftward
              x0(dof, k, C%xe+1, j) = mpi_recv_bfr_e(d2,k,j)
              d2 = d2+1
            endif
          enddo
        enddo
      enddo

      do i=C%xs,C%xe
        do k=C%zs,C%ze
          d1=1; d2=1
          do idof=i0, solver%diffside%dof-1
            dof = solver%difftop%dof + solver%diffside%dof + idof
            if (solver%diffside%is_inward(i1+idof)) then
              x0(dof, k, i, C%ys) = mpi_recv_bfr_s(d1,k,i)
              d1 = d1+1
            else
              x0(dof, k, i, C%ye+1) = mpi_recv_bfr_n(d2,k,i)
              d2 = d2+1
            endif
          enddo
        enddo
      enddo

      call restoreVecPointer(C%da, v0, x01d, x0)
    end associate
    ierr = 0
  end subroutine

  subroutine explicit_ediff_forward_sweep(solver, coeffs, dx, dy, dz, b, x)
  class(t_solver), intent(inout) :: solver
    real(ireals), target, intent(in) :: coeffs(:,:,:,:)
    integer(iintegers), dimension(3), intent(in) :: dx, dy, dz ! start, end, increment for each dimension
    type(tVec), intent(in) :: b, x

    real(ireals),pointer,dimension(:,:,:,:) :: x0=>null(), xb=>null()
    real(ireals),pointer,dimension(:) :: x01d=>null(), xb1d=>null()
    integer(iintegers) :: k,i,j
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:,:) ! dim(src, dst)
    integer(iintegers) :: msrc, mdst

    associate( &
        & atm => solver%atm, &
        & C   => solver%C_diff)

      call getVecPointer(C%da, x, x01d, x0)
      call getVecPointer(C%da, b, xb1d, xb, readonly=.True.)


      if(dz(3).lt.0) then ! if going from bottom to top, we do it here at the beginning
        do j=dy(1), dy(2), dy(3)
          do i=dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof-1
              if (.not.solver%difftop%is_inward(i1+idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i,j)
              endif
            enddo
          enddo
        enddo
      endif

      ! forward sweep through v0
      do k=dz(1), dz(2), dz(3)
        if( atm%l1d(atmk(atm,k)) ) then
          do j=dy(1), dy(2), dy(3)
            do i=dx(1), dx(2), dx(3)
              do idst = 0, solver%difftop%dof-1
                if (solver%difftop%is_inward(i1+idst)) then ! edn
                  x0(idst, k+i1, i, j) = xb(idst, k+1, i, j) + &
                    & x0(        idst , k   , i, j) * atm%a11(atmk(atm,k),i,j) + &
                    & x0(inv_dof(idst), k+i1, i, j) * atm%a12(atmk(atm,k),i,j)
                else ! eup
                  x0(idst, k, i, j) = xb(idst, k, i, j) + &
                    & x0(        idst , k+i1, i, j) * atm%a11(atmk(atm,k),i,j) + &
                    & x0(inv_dof(idst), k   , i, j) * atm%a12(atmk(atm,k),i,j)
                endif
              enddo
            enddo
          enddo
        else
          do j=dy(1), dy(2), dy(3)
            do i=dx(1), dx(2), dx(3)

              v(0:C%dof-1, 0:C%dof-1) => coeffs(1:C%dof**2,k-C%zs+1,i-C%xs+1,j-C%ys+1)

              dst = 0
              do idst = 0, solver%difftop%dof-1
                mdst = merge(k+1, k, solver%difftop%is_inward(i1+idst))
                x0(dst, mdst, i, j) = xb(dst, mdst, i, j)
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  x0(dst, mdst, i, j) = x0(dst, mdst, i, j) + x0(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  x0(dst, mdst, i, j) = x0(dst, mdst, i, j) + x0(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  x0(dst, mdst, i, j) = x0(dst, mdst, i, j) + x0(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%diffside%dof-1
                mdst = merge(i+1, i, solver%diffside%is_inward(i1+idst))
                x0(dst, k, mdst, j) = xb(dst, k, mdst, j)
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  x0(dst, k, mdst, j) = x0(dst, k, mdst, j) + x0(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  x0(dst, k, mdst, j) = x0(dst, k, mdst, j) + x0(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  x0(dst, k, mdst, j) = x0(dst, k, mdst, j) + x0(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%diffside%dof-1
                mdst = merge(j+1, j, solver%diffside%is_inward(i1+idst))
                x0(dst, k, i, mdst) = xb(dst, k, i, mdst)
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  x0(dst, k, i, mdst) = x0(dst, k, i, mdst) + x0(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  x0(dst, k, i, mdst) = x0(dst, k, i, mdst) + x0(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  x0(dst, k, i, mdst) = x0(dst, k, i, mdst) + x0(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

            enddo
          enddo
        endif
      enddo

      if(dz(3).gt.0) then ! if going from top to bottom, we do it here
        do j=dy(1), dy(2), dy(3)
          do i=dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof-1
              if (.not.solver%difftop%is_inward(i1+idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i,j)
              endif
            enddo
          enddo
        enddo
      endif

      call restoreVecPointer(C%da, x, x01d, x0)
      call restoreVecPointer(C%da, b, xb1d, xb, readonly=.True.)
    end associate

  contains

    pure function inv_dof(dof) ! returns the dof that is the same stream but the opposite direction
      integer(iintegers), intent(in) :: dof
      integer(iintegers) :: inv_dof, inc
      if(solver%difftop%is_inward(1)) then ! starting with downward streams
        inc = 1
      else
        inc = -1
      endif
      if(solver%difftop%is_inward(i1+dof)) then ! downward stream
        inv_dof = dof + inc
      else
        inv_dof = dof - inc
      endif
    end function
  end subroutine

  subroutine explicit_ediff_sor_sweep(solver, coeffs, dx, dy, dz, omega, b, x)
  class(t_solver), intent(inout) :: solver
    real(ireals), target, intent(in) :: coeffs(:,:,:,:)
    integer(iintegers), dimension(3), intent(in) :: dx, dy, dz ! start, end, increment for each dimension
    real(ireals), intent(in) :: omega
    type(tVec), intent(in) :: b, x

    real(ireals),pointer,dimension(:,:,:,:) :: x0=>null(), xb=>null()
    real(ireals),pointer,dimension(:) :: x01d=>null(), xb1d=>null()
    integer(iintegers) :: k,i,j
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:,:) ! dim(src, dst)
    integer(iintegers) :: msrc, mdst

    real(ireals), parameter :: diag = 1
    real(ireals) :: sigma

    associate( &
        & atm => solver%atm, &
        & C   => solver%C_diff)

      call getVecPointer(C%da, x, x01d, x0)
      call getVecPointer(C%da, b, xb1d, xb, readonly=.True.)

      if(dz(3).lt.0) then ! if going from bottom to top, we do it here at the beginning
        do j=dy(1), dy(2), dy(3)
          do i=dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof-1
              if (.not.solver%difftop%is_inward(i1+idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i,j)
              endif
            enddo
          enddo
        enddo
      endif

      ! forward sweep through v0
      do j=dy(1), dy(2), dy(3)
        do i=dx(1), dx(2), dx(3)
          do k=dz(1), dz(2), dz(3)
            if( atm%l1d(atmk(atm,k)) ) then
              do idst = 0, solver%difftop%dof-1
                if (solver%difftop%is_inward(i1+idst)) then ! edn
                  x0(idst, k+i1, i, j) = xb(idst, k+1, i, j) + &
                    & x0(        idst , k   , i, j) * atm%a11(atmk(atm,k),i,j) + &
                    & x0(inv_dof(idst), k+i1, i, j) * atm%a12(atmk(atm,k),i,j)
                else ! eup
                  x0(idst, k, i, j) = xb(idst, k, i, j) + &
                    & x0(        idst , k+i1, i, j) * atm%a11(atmk(atm,k),i,j) + &
                    & x0(inv_dof(idst), k   , i, j) * atm%a12(atmk(atm,k),i,j)
                endif
              enddo
            else

              v(0:C%dof-1, 0:C%dof-1) => coeffs(1:C%dof**2,k-C%zs+1,i-C%xs+1,j-C%ys+1)

              dst = 0
              do idst = 0, solver%difftop%dof-1
                mdst = merge(k+1, k, solver%difftop%is_inward(i1+idst))
                sigma = 0
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  sigma = sigma + x0(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  sigma = sigma + x0(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  sigma = sigma + x0(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                x0(dst, mdst, i, j) = (one - omega) * x0(dst, mdst, i, j) + (omega / diag) * (xb(dst, mdst, i, j) + sigma)
                dst = dst+1
              enddo

              do idst = 0, solver%diffside%dof-1
                mdst = merge(i+1, i, solver%diffside%is_inward(i1+idst))
                sigma = 0
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  sigma = sigma + x0(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  sigma = sigma + x0(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  sigma = sigma + x0(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                x0(dst, k, mdst, j) = (one - omega) * x0(dst, k, mdst, j) + (omega / diag) * (xb(dst, k, mdst, j) + sigma)
                dst = dst+1
              enddo

              do idst = 0, solver%diffside%dof-1
                mdst = merge(j+1, j, solver%diffside%is_inward(i1+idst))
                sigma = 0
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  sigma = sigma + x0(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  sigma = sigma + x0(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  sigma = sigma + x0(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                x0(dst, k, i, mdst) = (one - omega) * x0(dst, k, i, mdst) + (omega / diag) * (xb(dst, k, i, mdst) + sigma)
                dst = dst+1
              enddo

            endif ! endif l1d
          enddo
        enddo
      enddo

      if(dz(3).gt.0) then ! if going from top to bottom, we do it here
        do j=dy(1), dy(2), dy(3)
          do i=dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof-1
              if (.not.solver%difftop%is_inward(i1+idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i,j)
              endif
            enddo
          enddo
        enddo
      endif

      call restoreVecPointer(C%da, x, x01d, x0)
      call restoreVecPointer(C%da, b, xb1d, xb, readonly=.True.)
    end associate

  contains

    pure function inv_dof(dof) ! returns the dof that is the same stream but the opposite direction
      integer(iintegers), intent(in) :: dof
      integer(iintegers) :: inv_dof, inc
      if(solver%difftop%is_inward(1)) then ! starting with downward streams
        inc = 1
      else
        inc = -1
      endif
      if(solver%difftop%is_inward(i1+dof)) then ! downward stream
        inv_dof = dof + inc
      else
        inv_dof = dof - inc
      endif
    end function
  end subroutine
end module

