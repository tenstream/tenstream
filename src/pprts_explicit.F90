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

  use m_data_parameters, only: &
    & i0, i1, &
    & iintegers, mpiint, &
    & ireals, imp_ireals, &
    & zero, one

  use m_helper_functions, only: &
    & approx, &
    & CHKERR, &
    & get_petsc_opt, &
    & imp_min_mean_max, &
    & mpi_logical_and, &
    & toStr

  use m_pprts_base, only: &
    & atmk, &
    & determine_ksp_tolerances, &
    & setup_incSolar, &
    & t_solver, &
    & t_state_container

  use m_petsc_helpers, only: &
    getVecPointer, restoreVecPointer

  use m_adaptive_spectral_integration, only: polyfit

  implicit none

  private
  public :: &
    & exchange_diffuse_boundary, &
    & exchange_direct_boundary, &
    & explicit_ediff, &
    & explicit_ediff_sor_sweep, &
    & explicit_edir, &
    & explicit_edir_forward_sweep

  logical, parameter :: ldebug = .false.

contains

  !> @brief explicit loop to compute direct radiation
  subroutine explicit_edir(solver, prefix, edirTOA, solution, ierr)
    class(t_solver), target, intent(in) :: solver
    character(len=*), intent(in) :: prefix
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), intent(inout) :: solution

    real(ireals), pointer, dimension(:, :, :, :) :: x0 => null(), xg => null()
    real(ireals), pointer, dimension(:) :: x01d => null(), xg1d => null()
    type(tVec) :: v0, b, lb

    integer(iintegers), dimension(3) :: dx, dy ! start, end, increment for each dimension
    integer(iintegers) :: iter, maxiter, maxit_ignore

    real(ireals), allocatable :: residual(:)
    real(ireals) :: residual_mmm(3), rel_residual, atol, rtol, atol_default, rtol_default

    integer(iintegers), parameter :: default_max_it = 1000
    real(ireals) :: ignore_max_it! Ignore max iter setting if time is less

    logical :: lksp_view, lcomplete_initial_run
    logical :: lsun_north, lsun_east, lpermute, lskip_residual, lmonitor_residual, lflg
    logical :: laccept_incomplete_solve, lconverged_atol, lconverged_rtol, lconverged_reason

    integer(mpiint) :: ierr

    ierr = 0

    call PetscLogEventBegin(solver%logs%compute_Edir, ierr)
    associate ( &
        & vedir => solution%edir, &
        & sun => solver%sun, &
        & atm => solver%atm, &
        & C => solver%C_dir)

      maxiter = default_max_it
      call get_petsc_opt(prefix, "-ksp_max_it", maxiter, lflg, ierr); call CHKERR(ierr)

      lskip_residual = .false.
      call get_petsc_opt(prefix, "-ksp_skip_residual", lskip_residual, lflg, ierr); call CHKERR(ierr)

      ignore_max_it = -huge(ignore_max_it)
      call get_petsc_opt(prefix, "-ksp_ignore_max_it", ignore_max_it, lflg, ierr); call CHKERR(ierr)
      if (solution%time(1) .lt. ignore_max_it) then
        maxiter = default_max_it + 1
        lskip_residual = .false.
      end if

      call determine_ksp_tolerances(C, atm%unconstrained_fraction, &
        & rtol_default, atol_default, maxit_ignore)
      rtol = rtol_default
      atol = atol_default
      call get_petsc_opt(prefix, "-ksp_rtol", rtol, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(prefix, "-ksp_atol", atol, lflg, ierr); call CHKERR(ierr)

      lcomplete_initial_run = .true.
      call get_petsc_opt(prefix, "-ksp_complete_initial_run", lcomplete_initial_run, lflg, ierr); call CHKERR(ierr)
      if (lcomplete_initial_run .and. solution%dir_ksp_residual_history(1) .lt. 0) then
        maxiter = default_max_it + 2
        lskip_residual = .false.
        rtol = min(rtol, rtol_default)
        atol = min(atol, atol_default)
      end if

      allocate (residual(maxiter))

      lmonitor_residual = .false.
      call get_petsc_opt(prefix, "-ksp_monitor", lmonitor_residual, lflg, ierr); call CHKERR(ierr)
      if (.not. lmonitor_residual) then
        call get_petsc_opt(prefix, "-ksp_monitor_true_residual", lmonitor_residual, lflg, ierr); call CHKERR(ierr)
      end if

      lconverged_reason = lmonitor_residual
      call get_petsc_opt(prefix, "-ksp_converged_reason", lconverged_reason, lflg, ierr); call CHKERR(ierr)

      laccept_incomplete_solve = .false.
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(prefix, "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)

      lsun_north = sun%yinc .eq. i0
      lsun_east = sun%xinc .eq. i0

      dx = [C%xs, C%xe, i1]
      dy = [C%ys, C%ye, i1]

      lpermute = .true.
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-explicit_edir_permute", lpermute, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(prefix, "-explicit_edir_permute", lpermute, lflg, ierr); call CHKERR(ierr)
      if (lpermute) then
        if (lsun_east) dx = [dx(2), dx(1), -dx(3)]
        if (lsun_north) dy = [dy(2), dy(1), -dy(3)]
      end if

      lksp_view = .false.
      call get_petsc_opt(prefix, "-ksp_view", lksp_view, lflg, ierr); call CHKERR(ierr)
      if (solver%myid .eq. 0 .and. lksp_view) then
        print *, '* Using pprts explicit solver for prefix <'//trim(prefix)//'>'
        print *, '  -'//trim(prefix)//'ksp_max_it '//toStr(maxiter)
        print *, '  -'//trim(prefix)//'ksp_atol ', atol
        print *, '  -'//trim(prefix)//'ksp_rtol ', rtol
        print *, '  -'//trim(prefix)//'skip_residual '//toStr(lskip_residual)
        print *, '  -'//trim(prefix)//'accept_incomplete_solve '//toStr(laccept_incomplete_solve)
        print *, '  -'//trim(prefix)//'explicit_edir_permute '//toStr(lpermute)//&
          & ' horizontal-xy-iterator: ['//toStr(dx(3))//', '//toStr(dy(3))//']'
      end if

      call DMGetGlobalVector(C%da, b, ierr); call CHKERR(ierr)
      call VecSet(b, zero, ierr); call CHKERR(ierr)
      call setup_incSolar(solver, edirTOA, b)

      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call DMGlobalToLocal(C%da, b, INSERT_VALUES, lb, ierr); call CHKERR(ierr)

      call DMGetLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMGlobalToLocal(C%da, vedir, INSERT_VALUES, v0, ierr); call CHKERR(ierr)

      do iter = 1, maxiter

        call explicit_edir_forward_sweep(solver, solver%dir2dir, dx, dy, lb, v0)

        call exchange_direct_boundary(solver, lsun_north, lsun_east, v0, ierr); call CHKERR(ierr)

        ! Residual computations
        if (.not. lskip_residual) then
          call getVecPointer(C%da, vedir, xg1d, xg)
          call getVecPointer(C%da, v0, x01d, x0, readonly=.true.)
          residual(iter) = max(tiny(one), norm2(xg - x0(:, :, C%xs:C%xe, C%ys:C%ye)))
          xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)
          call restoreVecPointer(C%da, vedir, xg1d, xg)
          call restoreVecPointer(C%da, v0, x01d, x0, readonly=.true.)
          if (residual(1) .le. sqrt(tiny(residual))) then
            rel_residual = 0
          else
            rel_residual = residual(iter) / residual(1)
          end if
          if (solver%myid .eq. 0 .and. lmonitor_residual) then
            call imp_min_mean_max(solver%comm, residual(iter), residual_mmm)
            residual(iter) = residual_mmm(2)
            if (residual(1) .le. sqrt(tiny(residual))) then
              rel_residual = 0
            else
              rel_residual = residual(iter) / residual(1)
            end if
            print *, trim(prefix)//" iter "//toStr(iter)//' residual (min/mean/max)', residual_mmm, &
              & 'rel res', rel_residual
          end if
          solution%dir_ksp_residual_history(min(size(solution%dir_ksp_residual_history, kind=iintegers), iter)) = residual(iter)

          lconverged_atol = mpi_logical_and(solver%comm, residual(iter) .lt. atol)
          if (.not. lconverged_atol) then ! only reduce converged_RTOL if necessary
            lconverged_rtol = mpi_logical_and(solver%comm, rel_residual .lt. rtol)
          else
            lconverged_rtol = .false.
          end if
          if (lconverged_atol .or. lconverged_rtol) then
            if (solver%myid .eq. 0 .and. lconverged_reason) then
              if (lconverged_atol) then
                print *, trim(prefix)//' solve converged due to CONVERGED_ATOL iterations', iter
              else
                print *, trim(prefix)//' solve converged due to CONVERGED_RTOL iterations', iter
              end if
            end if
            exit
          end if
        else
          solution%dir_ksp_residual_history(min(size(solution%dir_ksp_residual_history, kind=iintegers), iter)) = zero
        end if

        if (iter .eq. maxiter) then
          if (.not. laccept_incomplete_solve) then
            call CHKERR(int(iter, mpiint), trim(prefix)//" did not converge")
          end if
        end if
      end do ! iter

      ! update solution vec
      call getVecPointer(C%da, vedir, xg1d, xg)
      call getVecPointer(C%da, v0, x01d, x0, readonly=.true.)
      xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)
      call restoreVecPointer(C%da, vedir, xg1d, xg)
      call restoreVecPointer(C%da, v0, x01d, x0, readonly=.true.)

      call DMRestoreLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call DMRestoreGlobalVector(C%da, b, ierr); call CHKERR(ierr)

      call PetscObjectSetName(vedir, 'debug_edir', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(vedir, PETSC_NULL_VEC, "-show_debug_edir", ierr); call CHKERR(ierr)
    end associate

    solution%lchanged = .true.
    solution%lWm2_dir = .false.
    call PetscLogEventEnd(solver%logs%compute_Edir, ierr)

  end subroutine

  subroutine exchange_direct_boundary(solver, lsun_north, lsun_east, x, ierr)
    class(t_solver), intent(in) :: solver
    logical, intent(in) :: lsun_north, lsun_east
    type(tVec), intent(inout) :: x
    integer(mpiint), intent(out) :: ierr

    real(ireals), pointer :: x0(:, :, :, :) => null(), x01d(:) => null()

    integer(mpiint), parameter :: tag_x = 1, tag_y = 2
    integer(mpiint) :: neigh_s, neigh_r, requests(4), statuses(mpi_status_size, 4)
    integer(iintegers) :: dofstart, dofend

    real(ireals), allocatable :: mpi_send_bfr_x(:, :, :), mpi_send_bfr_y(:, :, :)
    real(ireals), allocatable :: mpi_recv_bfr_x(:, :, :), mpi_recv_bfr_y(:, :, :)

    associate ( &
        & C => solver%C_dir)

      allocate (mpi_send_bfr_x(solver%dirside%dof, C%zm, C%ys:C%ye), mpi_recv_bfr_x(solver%dirside%dof, C%zm, C%ys:C%ye))
      allocate (mpi_send_bfr_y(solver%dirside%dof, C%zm, C%xs:C%xe), mpi_recv_bfr_y(solver%dirside%dof, C%zm, C%xs:C%xe))

      call getVecPointer(C%da, x, x01d, x0)

      ! Boundary exchanges
      ! x direction scatters
      dofstart = solver%dirtop%dof
      dofend = -1 + solver%dirtop%dof + solver%dirside%dof
      if (lsun_east) then
        if (solver%lopen_bc .and. C%xs .eq. i0) then
          mpi_send_bfr_x = 0
        else
          mpi_send_bfr_x = x0(dofstart:dofend, :, C%xs, C%ys:C%ye)
        end if
        neigh_s = int(C%neighbors(10), mpiint) ! neigh west
        neigh_r = int(C%neighbors(16), mpiint) ! neigh east
      else
        if (solver%lopen_bc .and. C%xe + 1 .eq. C%glob_xm) then
          mpi_send_bfr_x = 0
        else
          mpi_send_bfr_x = x0(dofstart:dofend, :, C%xe + 1, C%ys:C%ye)
        end if
        neigh_s = int(C%neighbors(16), mpiint) ! neigh east
        neigh_r = int(C%neighbors(10), mpiint) ! neigh west
      end if
      call MPI_Irecv(mpi_recv_bfr_x, size(mpi_recv_bfr_x, kind=mpiint), &
        & imp_ireals, neigh_r, tag_x, solver%comm, requests(1), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_x, size(mpi_send_bfr_x, kind=mpiint), &
        & imp_ireals, neigh_s, tag_x, solver%comm, requests(2), ierr); call CHKERR(ierr)

      ! y direction scatters
      dofstart = solver%dirtop%dof + solver%dirside%dof
      dofend = -1 + solver%dirtop%dof + solver%dirside%dof * 2
      if (lsun_north) then
        if (solver%lopen_bc .and. C%ys .eq. i0) then
          mpi_send_bfr_y = 0
        else
          mpi_send_bfr_y = x0(dofstart:dofend, :, C%xs:C%xe, C%ys)
        end if
        neigh_s = int(C%neighbors(4), mpiint) ! neigh south
        neigh_r = int(C%neighbors(22), mpiint) ! neigh north
      else
        if (solver%lopen_bc .and. C%ye + 1 .eq. C%glob_ym) then
          mpi_send_bfr_y = 0
        else
          mpi_send_bfr_y = x0(dofstart:dofend, :, C%xs:C%xe, C%ye + 1)
        end if
        neigh_s = int(C%neighbors(22), mpiint) ! neigh north
        neigh_r = int(C%neighbors(4), mpiint) ! neigh south
      end if
      call MPI_Irecv(mpi_recv_bfr_y, size(mpi_recv_bfr_y, kind=mpiint), &
        & imp_ireals, neigh_r, tag_y, solver%comm, requests(3), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_y, size(mpi_send_bfr_y, kind=mpiint), &
        & imp_ireals, neigh_s, tag_y, solver%comm, requests(4), ierr); call CHKERR(ierr)

      call MPI_Waitall(4_mpiint, requests, statuses, ierr); call CHKERR(ierr)

      dofstart = solver%dirtop%dof
      dofend = -1 + solver%dirtop%dof + solver%dirside%dof
      if (lsun_east) then
        x0(dofstart:dofend, :, C%xe + 1, C%ys:C%ye) = mpi_recv_bfr_x
      else
        x0(dofstart:dofend, :, C%xs, C%ys:C%ye) = mpi_recv_bfr_x
      end if

      dofstart = solver%dirtop%dof + solver%dirside%dof
      dofend = -1 + solver%dirtop%dof + solver%dirside%dof * 2
      if (lsun_north) then
        x0(dofstart:dofend, :, C%xs:C%xe, C%ye + 1) = mpi_recv_bfr_y
      else
        x0(dofstart:dofend, :, C%xs:C%xe, C%ys) = mpi_recv_bfr_y
      end if
      call restoreVecPointer(C%da, x, x01d, x0)
    end associate
    ierr = 0
  end subroutine

  subroutine explicit_edir_forward_sweep(solver, coeffs, dx, dy, b, x)
    class(t_solver), intent(in) :: solver
    real(ireals), target, intent(in) :: coeffs(:, :, :, :)
    integer(iintegers), dimension(3), intent(in) :: dx, dy ! start, end, increment for each dimension
    type(tVec), intent(in) :: b, x

    real(ireals), pointer :: x0(:, :, :, :) => null(), x01d(:) => null()
    real(ireals), pointer :: xb(:, :, :, :) => null(), xb1d(:) => null()

    integer(iintegers) :: i, j, k
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:, :) ! dim(src, dst)
    logical :: lsun_north, lsun_east

    associate ( &
        & atm => solver%atm, &
        & C => solver%C_dir, &
        & xinc => solver%sun%xinc, &
        & yinc => solver%sun%yinc)

      call getVecPointer(C%da, x, x01d, x0)
      call getVecPointer(C%da, b, xb1d, xb, readonly=.true.)
      x0(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, C%ys:C%ye) = xb(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, C%ys:C%ye)

      if (solver%lopen_bc) then
        lsun_north = yinc .eq. i0
        lsun_east = xinc .eq. i0

        if (lsun_north) then
          x0(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ye + 1) = &
            & xb(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ye + 1)
        else
          x0(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ys) = &
            & xb(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ys)
        end if
        if (lsun_east) then
          x0(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xe + 1, C%ys:C%ye) = &
            & xb(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xe + 1, C%ys:C%ye)
        else
          x0(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xs, C%ys:C%ye) = &
            & xb(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xs, C%ys:C%ye)
        end if
      end if

      ! forward sweep through x
      do k = C%zs, C%ze - 1
        if (atm%l1d(atmk(atm, k))) then
          do j = dy(1), dy(2), dy(3)
            do i = dx(1), dx(2), dx(3)

              do idst = 0, solver%dirtop%dof - 1
                x0(idst, k + i1, i, j) = x0(idst, k, i, j) * atm%a33(atmk(atm, k), i, j)
              end do
            end do
          end do
        else
          do j = dy(1), dy(2), dy(3)
            do i = dx(1), dx(2), dx(3)

              v(0:C%dof - 1, 0:C%dof - 1) => coeffs(:, k - C%zs + 1, i - C%xs + 1, j - C%ys + 1)

              dst = 0
              do idst = 0, solver%dirtop%dof - 1
                x0(dst, k + i1, i, j) = 0
                src = 0
                do isrc = 0, solver%dirtop%dof - 1
                  x0(dst, k + i1, i, j) = x0(dst, k + i1, i, j) + x0(src, k, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  x0(dst, k + i1, i, j) = x0(dst, k + i1, i, j) + x0(src, k, i + 1 - xinc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  x0(dst, k + i1, i, j) = x0(dst, k + i1, i, j) + x0(src, k, i, j + 1 - yinc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%dirside%dof - 1
                x0(dst, k, i + xinc, j) = 0
                src = 0
                do isrc = 0, solver%dirtop%dof - 1
                  x0(dst, k, i + xinc, j) = x0(dst, k, i + xinc, j) + x0(src, k, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  x0(dst, k, i + xinc, j) = x0(dst, k, i + xinc, j) + x0(src, k, i + 1 - xinc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  x0(dst, k, i + xinc, j) = x0(dst, k, i + xinc, j) + x0(src, k, i, j + 1 - yinc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%dirside%dof - 1
                x0(dst, k, i, j + yinc) = 0
                src = 0
                do isrc = 0, solver%dirtop%dof - 1
                  x0(dst, k, i, j + yinc) = x0(dst, k, i, j + yinc) + x0(src, k, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  x0(dst, k, i, j + yinc) = x0(dst, k, i, j + yinc) + x0(src, k, i + 1 - xinc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  x0(dst, k, i, j + yinc) = x0(dst, k, i, j + yinc) + x0(src, k, i, j + 1 - yinc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

            end do
          end do
        end if
      end do
      call restoreVecPointer(C%da, b, xb1d, xb, readonly=.true.)
      call restoreVecPointer(C%da, x, x01d, x0)
    end associate
  end subroutine

  !> @brief explicit loop to compute diffuse radiation
  subroutine explicit_ediff(solver, prefix, vb, solution, ierr)
    class(t_solver), intent(inout) :: solver
    character(len=*), intent(in) :: prefix
    type(tVec), intent(in) :: vb
    type(t_state_container), intent(inout) :: solution

    real(ireals), pointer, dimension(:, :, :, :) :: x0 => null(), xg => null()
    real(ireals), pointer, dimension(:) :: x01d => null(), xg1d => null()
    type(tVec) :: lvb, v0

    integer(iintegers) :: iter, isub, maxiter, sub_iter, maxit_ignore

    real(ireals), allocatable :: residual(:)
    real(ireals) :: residual_mmm(3), rel_residual, atol, rtol, atol_default, rtol_default

    integer(iintegers), parameter :: default_max_it = 10000
    real(ireals) :: ignore_max_it! Ignore max iter setting if time is less

    logical :: lksp_view, lcomplete_initial_run
    logical :: lskip_residual, lmonitor_residual, lconverged_atol, lconverged_rtol, lflg
    logical :: laccept_incomplete_solve, lconverged_reason

    logical :: lomega_set, ladaptive_omega
    real(ireals) :: omega, omega_adaptive, omega_increment, omega_min, omega_max

    integer(iintegers), parameter :: omega_pfN = 10
    integer(iintegers) :: i
    real(ireals) :: omega_pfx(omega_pfN), omega_pfy(omega_pfN), pf(2)

    integer(mpiint) :: ierr

    ierr = 0

    call PetscLogEventBegin(solver%logs%compute_Ediff, ierr)
    associate ( &
        & vediff => solution%ediff, &
        & atm => solver%atm, &
        & C => solver%C_diff)

      maxiter = default_max_it
      call get_petsc_opt(prefix, "-ksp_max_it", maxiter, lflg, ierr); call CHKERR(ierr)

      lskip_residual = .false.
      call get_petsc_opt(prefix, "-ksp_skip_residual", lskip_residual, lflg, ierr); call CHKERR(ierr)

      ignore_max_it = -huge(ignore_max_it)
      call get_petsc_opt(prefix, "-ksp_ignore_max_it", ignore_max_it, lflg, ierr); call CHKERR(ierr)
      if (solution%time(1) .lt. ignore_max_it) then
        maxiter = default_max_it + 1
        lskip_residual = .false.
      end if

      call determine_ksp_tolerances(C, atm%unconstrained_fraction, &
        & rtol_default, atol_default, maxit_ignore)
      rtol = rtol_default
      atol = atol_default
      call get_petsc_opt(prefix, "-ksp_rtol", rtol, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(prefix, "-ksp_atol", atol, lflg, ierr); call CHKERR(ierr)

      laccept_incomplete_solve = .false.
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(prefix, "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)

      omega = 1
      call get_petsc_opt(prefix, "-pc_sor_omega", omega, lomega_set, ierr); call CHKERR(ierr)
      omega_adaptive = omega
      ladaptive_omega = .true.
      call get_petsc_opt(prefix, "-pc_sor_omega_adaptive", ladaptive_omega, lflg, ierr); call CHKERR(ierr)
      omega_increment = .1_ireals
      call get_petsc_opt(prefix, "-pc_sor_omega_increment", omega_increment, lflg, ierr); call CHKERR(ierr)
      omega_min = 1._ireals
      call get_petsc_opt(prefix, "-pc_sor_omega_min", omega_min, lflg, ierr); call CHKERR(ierr)
      omega_max = 1.25_ireals
      call get_petsc_opt(prefix, "-pc_sor_omega_max", omega_max, lflg, ierr); call CHKERR(ierr)

      lcomplete_initial_run = .true.
      call get_petsc_opt(prefix, "-ksp_complete_initial_run", lcomplete_initial_run, lflg, ierr); call CHKERR(ierr)
      if (lcomplete_initial_run .and. (solution%diff_ksp_residual_history(1) .lt. 0)) then
        maxiter = default_max_it + 2
        lskip_residual = .false.
        laccept_incomplete_solve = .false.
        rtol = min(rtol, rtol_default)
        atol = min(atol, atol_default)
      end if

      allocate (residual(maxiter))
      sub_iter = 1
      call get_petsc_opt(prefix, "-pc_sub_it", sub_iter, lflg, ierr); call CHKERR(ierr)

      lmonitor_residual = .false.
      call get_petsc_opt(prefix, "-ksp_monitor", lmonitor_residual, lflg, ierr); call CHKERR(ierr)
      if (.not. lmonitor_residual) then
        call get_petsc_opt(prefix, "-ksp_monitor_true_residual", lmonitor_residual, lflg, ierr); call CHKERR(ierr)
      end if

      lconverged_reason = lmonitor_residual
      call get_petsc_opt(prefix, "-ksp_converged_reason", lconverged_reason, lflg, ierr); call CHKERR(ierr)

      lksp_view = .false.
      call get_petsc_opt(prefix, "-ksp_view", lksp_view, lflg, ierr); call CHKERR(ierr)
      if (solver%myid .eq. 0 .and. lksp_view) then
        print *, '* Using pprts explicit solver for prefix <'//trim(prefix)//'>'
        print *, '  -'//trim(prefix)//'ksp_max_it '//toStr(maxiter)
        print *, '  -'//trim(prefix)//'pc_sub_it '//toStr(sub_iter)
        print *, '  -'//trim(prefix)//'ksp_atol ', atol
        print *, '  -'//trim(prefix)//'ksp_rtol ', rtol
        print *, '  -'//trim(prefix)//'pc_sor_omega '//toStr(omega)
        print *, '  -'//trim(prefix)//'skip_residual '//toStr(lskip_residual)
        print *, '  -'//trim(prefix)//'accept_incomplete_solve '//toStr(laccept_incomplete_solve)
      end if

      call DMGetLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(C%da, vediff, INSERT_VALUES, v0, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd(C%da, vediff, INSERT_VALUES, v0, ierr); call CHKERR(ierr)

      call DMGetLocalVector(C%da, lvb, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(C%da, vb, INSERT_VALUES, lvb, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd(C%da, vb, INSERT_VALUES, lvb, ierr); call CHKERR(ierr)

      do iter = 1, maxiter
        do isub = 1, sub_iter
          if (approx(omega_adaptive, 1._ireals)) then
            if (modulo(iter + isub, 2) .eq. 0) then
              call explicit_ediff_forward_sweep(&
                & solver, &
                & solver%diff2diff, &
                & dx=[C%xs, C%xe, i1], &
                & dy=[C%ys, C%ye, i1], &
                & dz=[C%zs, C%ze - 1, i1], &
                & b=lvb, x=v0)
            else
              call explicit_ediff_forward_sweep(&
                & solver, &
                & solver%diff2diff, &
                & dx=[C%xe, C%xs, -i1], &
                & dy=[C%ye, C%ys, -i1], &
                & dz=[C%ze - 1, C%zs, -i1], &
                & b=lvb, x=v0)
            end if
          else
            if (modulo(iter + isub, 2) .eq. 0) then
              call explicit_ediff_sor_sweep(&
                & solver, &
                & solver%diff2diff, &
                & dx=[C%xs, C%xe, i1], &
                & dy=[C%ys, C%ye, i1], &
                & dz=[C%zs, C%ze - 1, i1], &
                & omega=omega_adaptive, &
                & b=lvb, x=v0)
            else
              call explicit_ediff_sor_sweep(&
                & solver, &
                & solver%diff2diff, &
                & dx=[C%xe, C%xs, -i1], &
                & dy=[C%ye, C%ys, -i1], &
                & dz=[C%ze - 1, C%zs, -i1], &
                & omega=omega_adaptive, &
                & b=lvb, x=v0)
            end if
          end if
        end do

        call exchange_diffuse_boundary(solver, v0, ierr); call CHKERR(ierr)

        ! Residual computations
        if (.not. lskip_residual) then
          call getVecPointer(C%da, vediff, xg1d, xg)
          call getVecPointer(C%da, v0, x01d, x0, readonly=.true.)

          residual(iter) = norm2(xg - x0(:, :, C%xs:C%xe, C%ys:C%ye))
          xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)

          call restoreVecPointer(C%da, vediff, xg1d, xg)
          call restoreVecPointer(C%da, v0, x01d, x0, readonly=.true.)

          if (residual(1) .le. sqrt(tiny(residual))) then
            rel_residual = 0
          else
            rel_residual = residual(iter) / residual(1)
          end if
          if (solver%myid .eq. 0 .and. lmonitor_residual) then
            call imp_min_mean_max(solver%comm, residual(iter), residual_mmm)
            residual(iter) = residual_mmm(2)
            if (residual(1) .le. sqrt(tiny(residual))) then
              rel_residual = 0
            else
              rel_residual = residual(iter) / residual(1)
            end if
            print *, trim(prefix), ' iter ', toStr(iter), ' residual (min/mean/max)', residual_mmm, &
              & 'rel res', rel_residual, 'omega', omega_adaptive
          end if
          solution%diff_ksp_residual_history(min(size(solution%diff_ksp_residual_history, kind=iintegers), iter)) = residual(iter)

          lconverged_atol = mpi_logical_and(solver%comm, residual(iter) .lt. atol)
          if (.not. lconverged_atol) then ! only reduce converged_RTOL if necessary
            lconverged_rtol = mpi_logical_and(solver%comm, rel_residual .lt. rtol)
          else
            lconverged_rtol = .false.
          end if
          if (lconverged_atol .or. lconverged_rtol) then
            if (solver%myid .eq. 0 .and. lconverged_reason) then
              if (lconverged_atol) then
                print *, trim(prefix)//' solve converged due to CONVERGED_ATOL iterations', iter
              else
                print *, trim(prefix)//' solve converged due to CONVERGED_RTOL iterations', iter
              end if
            end if
            exit
          end if

          if (ladaptive_omega) then
            do i = 1, omega_pfN
              omega_pfx(i) = real(i)
              omega_pfy(i) = residual(max(i1, iter - omega_pfN + i)) / residual(1)
            end do
            pf(:) = polyfit(omega_pfx, omega_pfy, 1_iintegers, ierr)
            !print *,'polyfit', omega_pfx, omega_pfy, '=>', pf
            if (ierr .eq. 0) then
              if (pf(2) .lt. -0.1_ireals) then ! converging very well
                ! keep on as we have, we're doing fine
              else if (pf(2) .lt. -0.001_ireals) then ! converging ok
                omega_adaptive = omega_adaptive + omega_increment
              else if (pf(2) .lt. 0._ireals) then ! converging slowly
                omega_adaptive = 1._ireals
              else if (pf(2) .gt. 1._ireals) then ! diverging extremely
                omega_adaptive = 1._ireals
              else if (pf(2) .gt. 0._ireals) then ! diverging slowly
                omega_adaptive = omega_adaptive - omega_increment
              end if
            end if
            omega_adaptive = min(max(omega_adaptive, omega_min), omega_max)
            if (residual(iter) .lt. 10 * atol) omega_adaptive = 1._ireals ! nearly converged, use stable omega
          end if

        else
          solution%diff_ksp_residual_history(min(size(solution%diff_ksp_residual_history, kind=iintegers), iter)) = zero
        end if

        if (iter .eq. maxiter) then
          if (.not. laccept_incomplete_solve) then
            call CHKERR(int(iter, mpiint), trim(prefix)//" did not converge")
          end if
        end if
      end do ! iter

      call getVecPointer(C%da, vediff, xg1d, xg)
      call getVecPointer(C%da, v0, x01d, x0, readonly=.true.)

      ! update solution vec
      xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)

      call restoreVecPointer(C%da, vediff, xg1d, xg)
      call restoreVecPointer(C%da, v0, x01d, x0, readonly=.true.)

      call DMRestoreLocalVector(C%da, v0, ierr); call CHKERR(ierr)
      call DMRestoreLocalVector(C%da, lvb, ierr); call CHKERR(ierr)

      call PetscObjectSetName(vediff, 'debug_ediff', ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(vediff, PETSC_NULL_VEC, "-show_debug_ediff", ierr); call CHKERR(ierr)
    end associate

    solution%lchanged = .true.
    solution%lWm2_diff = .false.
    call PetscLogEventEnd(solver%logs%compute_Ediff, ierr)
  end subroutine

  subroutine exchange_diffuse_boundary(solver, v0, ierr)
    class(t_solver), intent(in) :: solver
    type(tVec), intent(inout) :: v0
    integer(mpiint), intent(out) :: ierr

    integer(mpiint), parameter :: tag_e = 1, tag_w = 2, tag_n = 3, tag_s = 4
    integer(mpiint) :: neigh_s, neigh_e, neigh_w, neigh_n
    integer(mpiint) :: requests(8), statuses(mpi_status_size, 8)
    integer(iintegers) :: k, i, j, d1, d2, dof, idof

    real(ireals), pointer :: x0(:, :, :, :) => null(), x01d(:) => null()

    real(ireals), allocatable :: mpi_send_bfr_e(:, :, :), mpi_send_bfr_n(:, :, :)
    real(ireals), allocatable :: mpi_send_bfr_w(:, :, :), mpi_send_bfr_s(:, :, :)
    real(ireals), allocatable :: mpi_recv_bfr_e(:, :, :), mpi_recv_bfr_n(:, :, :)
    real(ireals), allocatable :: mpi_recv_bfr_w(:, :, :), mpi_recv_bfr_s(:, :, :)

    associate ( &
        & C => solver%C_diff)

      allocate (&
        & mpi_send_bfr_e(solver%diffside%dof / 2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_send_bfr_w(solver%diffside%dof / 2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_send_bfr_n(solver%diffside%dof / 2, C%zs:C%ze, C%xs:C%xe), &
        & mpi_send_bfr_s(solver%diffside%dof / 2, C%zs:C%ze, C%xs:C%xe), &
        & mpi_recv_bfr_e(solver%diffside%dof / 2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_recv_bfr_w(solver%diffside%dof / 2, C%zs:C%ze, C%ys:C%ye), &
        & mpi_recv_bfr_n(solver%diffside%dof / 2, C%zs:C%ze, C%xs:C%xe), &
        & mpi_recv_bfr_s(solver%diffside%dof / 2, C%zs:C%ze, C%xs:C%xe) &
        & )

      call getVecPointer(C%da, v0, x01d, x0)

      neigh_s = int(C%neighbors(4), mpiint)
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
      do j = C%ys, C%ye
        do k = C%zs, C%ze
          d1 = 1; d2 = 1
          do idof = i0, solver%diffside%dof - 1
            dof = solver%difftop%dof + idof
            if (solver%diffside%is_inward(i1 + idof)) then ! to the right
              mpi_send_bfr_e(d1, k, j) = x0(dof, k, C%xe + 1, j)
              d1 = d1 + 1
            else !leftward
              mpi_send_bfr_w(d2, k, j) = x0(dof, k, C%xs, j)
              d2 = d2 + 1
            end if
          end do
        end do
      end do
      call MPI_Isend(mpi_send_bfr_w, size(mpi_send_bfr_w, kind=mpiint), &
        & imp_ireals, neigh_w, tag_w, solver%comm, requests(5), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_e, size(mpi_send_bfr_e, kind=mpiint), &
        & imp_ireals, neigh_e, tag_e, solver%comm, requests(6), ierr); call CHKERR(ierr)

      ! y direction scatters, south/north boundary
      do i = C%xs, C%xe
        do k = C%zs, C%ze
          d1 = 1; d2 = 1
          do idof = i0, solver%diffside%dof - 1
            dof = solver%difftop%dof + solver%diffside%dof + idof
            if (solver%diffside%is_inward(i1 + idof)) then ! forward
              mpi_send_bfr_n(d1, k, i) = x0(dof, k, i, C%ye + 1)
              d1 = d1 + 1
            else ! backward
              mpi_send_bfr_s(d2, k, i) = x0(dof, k, i, C%ys)
              d2 = d2 + 1
            end if
          end do
        end do
      end do
      call MPI_Isend(mpi_send_bfr_s, size(mpi_send_bfr_s, kind=mpiint), &
        & imp_ireals, neigh_s, tag_s, solver%comm, requests(7), ierr); call CHKERR(ierr)
      call MPI_Isend(mpi_send_bfr_n, size(mpi_send_bfr_n, kind=mpiint), &
        & imp_ireals, neigh_n, tag_n, solver%comm, requests(8), ierr); call CHKERR(ierr)

      call MPI_Waitall(8_mpiint, requests, statuses, ierr); call CHKERR(ierr)

      ! receive buffers
      do j = C%ys, C%ye
        do k = C%zs, C%ze
          d1 = 1; d2 = 1
          do idof = i0, solver%diffside%dof - 1
            dof = solver%difftop%dof + idof
            if (solver%diffside%is_inward(i1 + idof)) then ! to the right
              x0(dof, k, C%xs, j) = mpi_recv_bfr_w(d1, k, j)
              d1 = d1 + 1
            else ! leftward
              x0(dof, k, C%xe + 1, j) = mpi_recv_bfr_e(d2, k, j)
              d2 = d2 + 1
            end if
          end do
        end do
      end do

      do i = C%xs, C%xe
        do k = C%zs, C%ze
          d1 = 1; d2 = 1
          do idof = i0, solver%diffside%dof - 1
            dof = solver%difftop%dof + solver%diffside%dof + idof
            if (solver%diffside%is_inward(i1 + idof)) then
              x0(dof, k, i, C%ys) = mpi_recv_bfr_s(d1, k, i)
              d1 = d1 + 1
            else
              x0(dof, k, i, C%ye + 1) = mpi_recv_bfr_n(d2, k, i)
              d2 = d2 + 1
            end if
          end do
        end do
      end do

      call restoreVecPointer(C%da, v0, x01d, x0)
    end associate
    ierr = 0
  end subroutine

  subroutine explicit_ediff_forward_sweep(solver, coeffs, dx, dy, dz, b, x)
    class(t_solver), intent(inout) :: solver
    real(ireals), target, intent(in) :: coeffs(:, :, :, :)
    integer(iintegers), dimension(3), intent(in) :: dx, dy, dz ! start, end, increment for each dimension
    type(tVec), intent(in) :: b, x

    real(ireals), pointer, dimension(:, :, :, :) :: x0 => null(), xb => null()
    real(ireals), pointer, dimension(:) :: x01d => null(), xb1d => null()
    integer(iintegers) :: k, i, j
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:, :) ! dim(src, dst)
    integer(iintegers) :: msrc, mdst

    associate ( &
        & atm => solver%atm, &
        & C => solver%C_diff)

      call getVecPointer(C%da, x, x01d, x0)
      call getVecPointer(C%da, b, xb1d, xb, readonly=.true.)

      if (dz(3) .lt. 0) then ! if going from bottom to top, we do it here at the beginning
        do j = dy(1), dy(2), dy(3)
          do i = dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof - 1
              if (.not. solver%difftop%is_inward(i1 + idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i, j)
              end if
            end do
          end do
        end do
      end if

      ! forward sweep through v0
      do k = dz(1), dz(2), dz(3)
        if (atm%l1d(atmk(atm, k))) then
          do j = dy(1), dy(2), dy(3)
            do i = dx(1), dx(2), dx(3)
              do idst = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + idst)) then ! edn
                  x0(idst, k + i1, i, j) = xb(idst, k + 1, i, j) + &
                    & x0(idst, k, i, j) * atm%a11(atmk(atm, k), i, j) + &
                    & x0(inv_dof(idst), k + i1, i, j) * atm%a12(atmk(atm, k), i, j)
                else ! eup
                  x0(idst, k, i, j) = xb(idst, k, i, j) + &
                    & x0(idst, k + i1, i, j) * atm%a11(atmk(atm, k), i, j) + &
                    & x0(inv_dof(idst), k, i, j) * atm%a12(atmk(atm, k), i, j)
                end if
              end do
            end do
          end do
        else
          do j = dy(1), dy(2), dy(3)
            do i = dx(1), dx(2), dx(3)

              v(0:C%dof - 1, 0:C%dof - 1) => coeffs(1:C%dof**2, k - C%zs + 1, i - C%xs + 1, j - C%ys + 1)

              dst = 0
              do idst = 0, solver%difftop%dof - 1
                mdst = merge(k + 1, k, solver%difftop%is_inward(i1 + idst))
                x0(dst, mdst, i, j) = xb(dst, mdst, i, j)
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  x0(dst, mdst, i, j) = x0(dst, mdst, i, j) + x0(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  x0(dst, mdst, i, j) = x0(dst, mdst, i, j) + x0(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  x0(dst, mdst, i, j) = x0(dst, mdst, i, j) + x0(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%diffside%dof - 1
                mdst = merge(i + 1, i, solver%diffside%is_inward(i1 + idst))
                x0(dst, k, mdst, j) = xb(dst, k, mdst, j)
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  x0(dst, k, mdst, j) = x0(dst, k, mdst, j) + x0(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  x0(dst, k, mdst, j) = x0(dst, k, mdst, j) + x0(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  x0(dst, k, mdst, j) = x0(dst, k, mdst, j) + x0(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%diffside%dof - 1
                mdst = merge(j + 1, j, solver%diffside%is_inward(i1 + idst))
                x0(dst, k, i, mdst) = xb(dst, k, i, mdst)
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  x0(dst, k, i, mdst) = x0(dst, k, i, mdst) + x0(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  x0(dst, k, i, mdst) = x0(dst, k, i, mdst) + x0(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  x0(dst, k, i, mdst) = x0(dst, k, i, mdst) + x0(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

            end do
          end do
        end if
      end do

      if (dz(3) .gt. 0) then ! if going from top to bottom, we do it here
        do j = dy(1), dy(2), dy(3)
          do i = dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof - 1
              if (.not. solver%difftop%is_inward(i1 + idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i, j)
              end if
            end do
          end do
        end do
      end if

      call restoreVecPointer(C%da, x, x01d, x0)
      call restoreVecPointer(C%da, b, xb1d, xb, readonly=.true.)
    end associate

  contains

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
  end subroutine

  subroutine explicit_ediff_sor_sweep(solver, coeffs, dx, dy, dz, omega, b, x)
    class(t_solver), intent(inout) :: solver
    real(ireals), target, intent(in) :: coeffs(:, :, :, :)
    integer(iintegers), dimension(3), intent(in) :: dx, dy, dz ! start, end, increment for each dimension
    real(ireals), intent(in) :: omega
    type(tVec), intent(in) :: b, x

    real(ireals), pointer, dimension(:, :, :, :) :: x0 => null(), xb => null()
    real(ireals), pointer, dimension(:) :: x01d => null(), xb1d => null()
    integer(iintegers) :: k, i, j
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:, :) ! dim(src, dst)
    integer(iintegers) :: msrc, mdst

    real(ireals), parameter :: diag = 1
    real(ireals) :: sigma

    associate ( &
        & atm => solver%atm, &
        & C => solver%C_diff)

      call getVecPointer(C%da, x, x01d, x0)
      call getVecPointer(C%da, b, xb1d, xb, readonly=.true.)

      if (dz(3) .lt. 0) then ! if going from bottom to top, we do it here at the beginning
        do j = dy(1), dy(2), dy(3)
          do i = dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof - 1
              if (.not. solver%difftop%is_inward(i1 + idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i, j)
              end if
            end do
          end do
        end do
      end if

      ! forward sweep through v0
      do j = dy(1), dy(2), dy(3)
        do i = dx(1), dx(2), dx(3)
          do k = dz(1), dz(2), dz(3)
            if (atm%l1d(atmk(atm, k))) then
              do idst = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + idst)) then ! edn
                  x0(idst, k + i1, i, j) = xb(idst, k + 1, i, j) + &
                    & x0(idst, k, i, j) * atm%a11(atmk(atm, k), i, j) + &
                    & x0(inv_dof(idst), k + i1, i, j) * atm%a12(atmk(atm, k), i, j)
                else ! eup
                  x0(idst, k, i, j) = xb(idst, k, i, j) + &
                    & x0(idst, k + i1, i, j) * atm%a11(atmk(atm, k), i, j) + &
                    & x0(inv_dof(idst), k, i, j) * atm%a12(atmk(atm, k), i, j)
                end if
              end do
            else

              v(0:C%dof - 1, 0:C%dof - 1) => coeffs(1:C%dof**2, k - C%zs + 1, i - C%xs + 1, j - C%ys + 1)

              dst = 0
              do idst = 0, solver%difftop%dof - 1
                mdst = merge(k + 1, k, solver%difftop%is_inward(i1 + idst))
                sigma = 0
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                x0(dst, mdst, i, j) = (one - omega) * x0(dst, mdst, i, j) + (omega / diag) * (xb(dst, mdst, i, j) + sigma)
                dst = dst + 1
              end do

              do idst = 0, solver%diffside%dof - 1
                mdst = merge(i + 1, i, solver%diffside%is_inward(i1 + idst))
                sigma = 0
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                x0(dst, k, mdst, j) = (one - omega) * x0(dst, k, mdst, j) + (omega / diag) * (xb(dst, k, mdst, j) + sigma)
                dst = dst + 1
              end do

              do idst = 0, solver%diffside%dof - 1
                mdst = merge(j + 1, j, solver%diffside%is_inward(i1 + idst))
                sigma = 0
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  sigma = sigma + x0(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                x0(dst, k, i, mdst) = (one - omega) * x0(dst, k, i, mdst) + (omega / diag) * (xb(dst, k, i, mdst) + sigma)
                dst = dst + 1
              end do

            end if ! endif l1d
          end do
        end do
      end do

      if (dz(3) .gt. 0) then ! if going from top to bottom, we do it here
        do j = dy(1), dy(2), dy(3)
          do i = dx(1), dx(2), dx(3)
            do idst = 0, solver%difftop%dof - 1
              if (.not. solver%difftop%is_inward(i1 + idst)) then ! Eup
                x0(idst, C%ze, i, j) = xb(idst, C%ze, i, j) + &
                  & x0(inv_dof(idst), C%ze, i, j) * atm%albedo(i, j)
              end if
            end do
          end do
        end do
      end if

      call restoreVecPointer(C%da, x, x01d, x0)
      call restoreVecPointer(C%da, b, xb1d, xb, readonly=.true.)
    end associate

  contains

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
  end subroutine
end module
