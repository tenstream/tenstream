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

module m_pprts_shell

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & i0, i1, &
    & ireals, iintegers, mpiint, &
    & zero, one

  use m_helper_functions, only: &
    & approx, &
    & CHKERR

  use m_pprts_base, only: &
    & atmk, &
    & t_coord, &
    & t_pprts_shell_ctx, &
    & t_solver

  use m_pprts_explicit, only: &
    & exchange_diffuse_boundary, &
    & exchange_direct_boundary, &
    & explicit_ediff_sor_sweep, &
    & explicit_edir_forward_sweep

  use m_petsc_helpers, only: &
    getVecPointer, restoreVecPointer

  implicit none

  private
  public :: &
    & op_mat_getdiagonal, &
    & op_mat_mult_ediff, &
    & op_mat_mult_edir, &
    & op_mat_sor_ediff, &
    & op_mat_sor_edir, &
    & setup_matshell

  abstract interface
    subroutine mat_mult_sub(A, x, b, ierr)
      import tMat, tVec, mpiint
      type(tMat), intent(in) :: A
      type(tVec), intent(in) :: x
      type(tVec), intent(inout) :: b
      integer(mpiint), intent(out) :: ierr
    end subroutine
  end interface
  abstract interface
    subroutine mat_sor_sub(A, b, omega, sortype, shift, its, lits, x, ierr)
      use petsc
      import ireals, iintegers, mpiint
      type(tMat), intent(in) :: A
      type(tVec), intent(in) :: b
      real(ireals), intent(in) :: omega
      MatSORType, intent(in) :: sortype
      real(ireals), intent(in) :: shift
      integer(iintegers), intent(in) :: its
      integer(iintegers), intent(in) :: lits
      type(tVec), intent(inout) :: x
      integer(mpiint), intent(out) :: ierr
    end subroutine
  end interface
  abstract interface
    subroutine mat_getdiagonal_sub(A, d, ierr)
      import tMat, tVec, mpiint
      type(tMat), intent(in) :: A
      type(tVec), intent(inout) :: d
      integer(mpiint), intent(out) :: ierr
    end subroutine
  end interface

contains

  !> @brief call the matrix assembly and petsc solve routines for pprts solvers
  subroutine setup_matshell(solver, C, A, mat_mult_subroutine, mat_sor_subroutine, mat_getdiagonal_subroutine)
    class(t_solver), target, intent(inout) :: solver
    type(t_coord), intent(in) :: C
    type(tMat), allocatable, intent(inout) :: A
    procedure(mat_mult_sub) :: mat_mult_subroutine
    procedure(mat_sor_sub) :: mat_sor_subroutine
    procedure(mat_getdiagonal_sub) :: mat_getdiagonal_subroutine

    integer(iintegers) :: Nlocal, Nglobal
    integer(mpiint) :: ierr

    solver%shell_ctx%solver => solver

    if (.not. allocated(A)) then
      allocate (A)
      Nlocal = C%zm * C%xm * C%ym * C%dof
      Nglobal = C%glob_zm * C%glob_xm * C%glob_ym * C%dof
      call MatCreateShell(solver%comm, Nlocal, Nlocal, Nglobal, Nglobal, solver%shell_ctx, A, ierr); call CHKERR(ierr)
      call MatShellSetContext(A, solver%shell_ctx, ierr); call CHKERR(ierr)
      call MatShellSetOperation(A, MATOP_MULT, mat_mult_subroutine, ierr); call CHKERR(ierr)
      call MatShellSetOperation(A, MATOP_SOR, mat_sor_subroutine, ierr); call CHKERR(ierr)
      call MatShellSetOperation(A, MATOP_GET_DIAGONAL, mat_getdiagonal_subroutine, ierr); call CHKERR(ierr)
    end if
  end subroutine

  !> @brief define the matmult function for dir2dir computations (used for MatShell)
  subroutine op_mat_mult_edir(A, x, b, ierr)
    type(tMat), intent(in) :: A
    type(tVec), intent(in) :: x
    type(tVec), intent(inout) :: b
    integer(mpiint), intent(out) :: ierr

    type(t_pprts_shell_ctx), pointer :: ctx_ptr
    class(t_solver), pointer :: solver
    real(ireals), pointer, dimension(:, :, :, :) :: xx => null(), xb => null()
    real(ireals), pointer, dimension(:) :: xx1d => null(), xb1d => null()
    real(ireals), pointer :: v(:, :) ! dim(src, dst)
    integer(iintegers) :: k, i, j
    integer(iintegers) :: idst, isrc, src, dst, xinc, yinc
    type(tVec) :: lb, lx

    nullify (ctx_ptr)
    call MatShellGetContext(A, ctx_ptr, ierr); call CHKERR(ierr)

    nullify (solver)
    solver => ctx_ptr%solver
    if (.not. associated(solver)) call CHKERR(1_mpiint, 'mat shell context has not been set!')

    call PetscObjectViewFromOptions(PetscObjectCast(x), PETSC_NULL_OBJECT, "-show_shell_x", ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(PetscObjectCast(b), PETSC_NULL_OBJECT, "-show_shell_b", ierr); call CHKERR(ierr)

    associate (sun => solver%sun, &
               atm => solver%atm, &
               C => solver%C_dir)

      call DMGetLocalVector(C%da, lx, ierr); call CHKERR(ierr)
      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call VecSet(lb, zero, ierr); call CHKERR(ierr)

      call DMGlobalToLocalBegin(C%da, x, INSERT_VALUES, lx, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd(C%da, x, INSERT_VALUES, lx, ierr); call CHKERR(ierr)

      call getVecPointer(C%da, lx, xx1d, xx, readonly=.true.)
      call getVecPointer(C%da, lb, xb1d, xb)

      do j = C%ys, C%ye
        do i = C%xs, C%xe
          do k = C%zs, C%ze - 1

            if (atm%l1d(atmk(atm, k))) then
              do idst = 0, solver%dirtop%dof - 1
                xb(idst, k + i1, i, j) = xb(idst, k + i1, i, j) - xx(idst, k, i, j) * atm%a33(atmk(atm, k), i, j)
              end do
            else

              xinc = sun%xinc
              yinc = sun%yinc

              v(0:C%dof - 1, 0:C%dof - 1) => solver%dir2dir(:, k, i, j)

              dst = 0
              do idst = 0, solver%dirtop%dof - 1
                src = 0
                do isrc = 0, solver%dirtop%dof - 1
                  xb(dst, k + i1, i, j) = xb(dst, k + i1, i, j) - xx(src, k, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  xb(dst, k + i1, i, j) = xb(dst, k + i1, i, j) - xx(src, k, i + 1 - xinc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  xb(dst, k + i1, i, j) = xb(dst, k + i1, i, j) - xx(src, k, i, j + 1 - yinc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%dirside%dof - 1
                src = 0
                do isrc = 0, solver%dirtop%dof - 1
                  xb(dst, k, i + xinc, j) = xb(dst, k, i + xinc, j) - xx(src, k, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  xb(dst, k, i + xinc, j) = xb(dst, k, i + xinc, j) - xx(src, k, i + 1 - xinc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  xb(dst, k, i + xinc, j) = xb(dst, k, i + xinc, j) - xx(src, k, i, j + 1 - yinc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%dirside%dof - 1
                src = 0
                do isrc = 0, solver%dirtop%dof - 1
                  xb(dst, k, i, j + yinc) = xb(dst, k, i, j + yinc) - xx(src, k, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  xb(dst, k, i, j + yinc) = xb(dst, k, i, j + yinc) - xx(src, k, i + 1 - xinc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%dirside%dof - 1
                  xb(dst, k, i, j + yinc) = xb(dst, k, i, j + yinc) - xx(src, k, i, j + 1 - yinc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do
            end if

          end do
        end do
      end do

      call restoreVecPointer(C%da, lb, xb1d, xb)
      call restoreVecPointer(C%da, lx, xx1d, xx, readonly=.true.)
      call DMRestoreLocalVector(C%da, lx, ierr); call CHKERR(ierr)

      call VecSet(b, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(C%da, lb, ADD_VALUES, b, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd(C%da, lb, ADD_VALUES, b, ierr); call CHKERR(ierr)

      ! vec copy simulates diagonal entries
      call VecAXPY(b, one, x, ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
    end associate
    call PetscObjectViewFromOptions(PetscObjectCast(b), PETSC_NULL_OBJECT, "-show_shell_rb", ierr); call CHKERR(ierr)

    ierr = 0
  end subroutine

  !> @brief define the matsor function for dir2dir computations (used for MatShell)
  subroutine op_mat_sor_edir(A, b, omega, sortype, shift, its, lits, x, ierr)
    type(tMat), intent(in) :: A
    type(tVec), intent(in) :: b
    real(ireals), intent(in) :: omega
    MatSORType, intent(in) :: sortype
    real(ireals), intent(in) :: shift
    integer(iintegers), intent(in) :: its
    integer(iintegers), intent(in) :: lits
    type(tVec), intent(inout) :: x
    integer(mpiint), intent(out) :: ierr

    type(t_pprts_shell_ctx), pointer :: ctx_ptr
    class(t_solver), pointer :: solver
    type(tVec) :: lx, lb

    real(ireals), pointer :: xx(:, :, :, :) => null(), xx1d(:) => null()
    real(ireals), pointer :: xg(:, :, :, :) => null(), xg1d(:) => null()
    real(ireals), pointer :: xb(:, :, :, :) => null(), xb1d(:) => null()

    integer(iintegers), dimension(3) :: dx, dy ! start, end, increment for each dimension
    integer(iintegers) :: i_its, i_lits
    logical :: lsweep_forward, lsweep_backward, lsweep_symmetric, lzero_initial_guess
    logical :: lsun_north, lsun_east

    call CHKERR(1_mpiint, 'Not yet implemented') ! well this seems implemented but its wrong.... explicit version should be better anyway

    nullify (ctx_ptr)
    call MatShellGetContext(A, ctx_ptr, ierr); call CHKERR(ierr)

    nullify (solver)
    solver => ctx_ptr%solver
    if (.not. associated(solver)) call CHKERR(1_mpiint, 'mat shell context has not been set!')

    if (.not. approx(shift, zero)) call CHKERR(1_mpiint, 'Shift has to be zero')
    if (.not. approx(omega, one)) call CHKERR(1_mpiint, 'Omega has to be one')

    lsweep_symmetric = iand(sortype%v, SOR_LOCAL_SYMMETRIC_SWEEP%v) .eq. SOR_LOCAL_SYMMETRIC_SWEEP%v .or. &
      & iand(sortype%v, SOR_SYMMETRIC_SWEEP%v) .eq. SOR_SYMMETRIC_SWEEP%v
    lsweep_backward = iand(sortype%v, SOR_LOCAL_BACKWARD_SWEEP%v) .eq. SOR_LOCAL_BACKWARD_SWEEP%v .or. &
      & iand(sortype%v, SOR_BACKWARD_SWEEP%v) .eq. SOR_BACKWARD_SWEEP%v
    lsweep_forward = iand(sortype%v, SOR_LOCAL_FORWARD_SWEEP%v) .eq. SOR_LOCAL_FORWARD_SWEEP%v .or. &
      & iand(sortype%v, SOR_FORWARD_SWEEP%v) .eq. SOR_FORWARD_SWEEP%v
    lzero_initial_guess = iand(sortype%v, SOR_ZERO_INITIAL_GUESS%v) .eq. SOR_ZERO_INITIAL_GUESS%v

    associate ( &
        & sun => solver%sun, &
        & C => solver%C_dir)

      lsun_north = sun%yinc .eq. i0
      lsun_east = sun%xinc .eq. i0

      dx = [C%xs, C%xe, i1]
      dy = [C%ys, C%ye, i1]
      if (lsun_east) dx = [dx(2), dx(1), -dx(3)]
      if (lsun_north) dy = [dy(2), dy(1), -dy(3)]

      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call DMGlobalToLocal(C%da, b, INSERT_VALUES, lb, ierr); call CHKERR(ierr)

      call DMGetLocalVector(C%da, lx, ierr); call CHKERR(ierr)
      if (lzero_initial_guess) then
        call VecSet(lx, zero, ierr); call CHKERR(ierr)
      else
        call DMGlobalToLocal(C%da, x, INSERT_VALUES, lx, ierr); call CHKERR(ierr)
      end if

      call getVecPointer(C%da, lx, xx1d, xx)
      call getVecPointer(C%da, lb, xb1d, xb, readonly=.true.)
      xx(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, C%ys:C%ye) = xb(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, C%ys:C%ye)
      call restoreVecPointer(C%da, lx, xx1d, xx)
      call restoreVecPointer(C%da, lb, xb1d, xb, readonly=.true.)

      do i_its = 1, its

        do i_lits = 1, lits
          call explicit_edir_forward_sweep(solver, solver%dir2dir, dx, dy, lb, lx)
        end do

        call exchange_direct_boundary(solver, lsun_north, lsun_east, lx, ierr); call CHKERR(ierr)

        ! update solution vec
        call getVecPointer(C%da, x, xg1d, xg)
        call getVecPointer(C%da, lx, xx1d, xx, readonly=.true.)
        xg = xx(:, :, C%xs:C%xe, C%ys:C%ye)
        call restoreVecPointer(C%da, x, xg1d, xg)
        call restoreVecPointer(C%da, lx, xx1d, xx, readonly=.true.)
      end do

      call DMRestoreLocalVector(C%da, lx, ierr); call CHKERR(ierr)
      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
    end associate
    ierr = 0
  end subroutine

  !> @brief define the mat_getdiagonal function for dir2dir computations (used for MatShell)
  subroutine op_mat_getdiagonal(A, d, ierr)
    type(tMat), intent(in) :: A
    type(tVec), intent(inout) :: d
    integer(mpiint), intent(out) :: ierr
    call VecSet(d, one, ierr); call CHKERR(ierr)
    if (.false.) print *, 'ignoring Mat argument', A
  end subroutine

  subroutine op_mat_mult_ediff(A, x, b, ierr)
    type(tMat), intent(in) :: A
    type(tVec), intent(in) :: x
    type(tVec), intent(inout) :: b
    integer(mpiint), intent(out) :: ierr

    type(t_pprts_shell_ctx), pointer :: ctx_ptr
    class(t_solver), pointer :: solver
    real(ireals), pointer, dimension(:, :, :, :) :: xx => null(), xb => null()
    real(ireals), pointer, dimension(:) :: xx1d => null(), xb1d => null()
    real(ireals), pointer :: v(:, :) ! dim(src,dst)
    integer(iintegers) :: k, i, j
    integer(iintegers) :: idst, isrc, dst, src
    integer(iintegers) :: mdst, msrc
    integer(iintegers), allocatable :: row(:, :), col(:, :)
    type(tVec) :: lb, lx

    nullify (ctx_ptr)
    call MatShellGetContext(A, ctx_ptr, ierr); call CHKERR(ierr)

    nullify (solver)
    solver => ctx_ptr%solver
    if (.not. associated(solver)) call CHKERR(1_mpiint, 'mat shell context has not been set!')

    call PetscObjectViewFromOptions(PetscObjectCast(x), PETSC_NULL_OBJECT, "-show_shell_x", ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(PetscObjectCast(b), PETSC_NULL_OBJECT, "-show_shell_b", ierr); call CHKERR(ierr)

    associate (atm => solver%atm, &
               C => solver%C_diff)

      allocate (row(4, 0:C%dof - 1), col(4, 0:C%dof - 1))

      call DMGetLocalVector(C%da, lx, ierr); call CHKERR(ierr)
      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call VecSet(lb, zero, ierr); call CHKERR(ierr)

      call DMGlobalToLocalBegin(C%da, x, INSERT_VALUES, lx, ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd(C%da, x, INSERT_VALUES, lx, ierr); call CHKERR(ierr)

      call getVecPointer(C%da, lx, xx1d, xx, readonly=.true.)
      call getVecPointer(C%da, lb, xb1d, xb)

      do j = C%ys, C%ye
        do i = C%xs, C%xe
          do k = C%zs, C%ze - 1

            if (atm%l1d(atmk(atm, k))) then

              do idst = 0, solver%difftop%dof - 1
                if (solver%difftop%is_inward(i1 + idst)) then ! edn
                  xb(idst, k + i1, i, j) = xb(idst, k + i1, i, j) - xx(idst, k, i, j) * atm%a11(atmk(atm, k), i, j)
                  xb(idst, k + i1, i, j) = xb(idst, k + i1, i, j) - xx(inv_dof(idst), k + i1, i, j) * atm%a12(atmk(atm, k), i, j)
                else ! eup
                  xb(idst, k, i, j) = xb(idst, k, i, j) - xx(idst, k + i1, i, j) * atm%a11(atmk(atm, k), i, j)
                  xb(idst, k, i, j) = xb(idst, k, i, j) - xx(inv_dof(idst), k, i, j) * atm%a12(atmk(atm, k), i, j)
                end if
              end do

            else

              v(0:C%dof - 1, 0:C%dof - 1) => solver%diff2diff(1:C%dof**2, k, i, j)

              dst = 0
              do idst = 0, solver%difftop%dof - 1
                mdst = merge(k + 1, k, solver%difftop%is_inward(i1 + idst))
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  xb(dst, mdst, i, j) = xb(dst, mdst, i, j) - xx(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  xb(dst, mdst, i, j) = xb(dst, mdst, i, j) - xx(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  xb(dst, mdst, i, j) = xb(dst, mdst, i, j) - xx(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%diffside%dof - 1
                mdst = merge(i + 1, i, solver%diffside%is_inward(i1 + idst))
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  xb(dst, k, mdst, j) = xb(dst, k, mdst, j) - xx(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  xb(dst, k, mdst, j) = xb(dst, k, mdst, j) - xx(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  xb(dst, k, mdst, j) = xb(dst, k, mdst, j) - xx(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

              do idst = 0, solver%diffside%dof - 1
                mdst = merge(j + 1, j, solver%diffside%is_inward(i1 + idst))
                src = 0
                do isrc = 0, solver%difftop%dof - 1
                  msrc = merge(k, k + 1, solver%difftop%is_inward(i1 + isrc))
                  xb(dst, k, i, mdst) = xb(dst, k, i, mdst) - xx(src, msrc, i, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(i, i + 1, solver%diffside%is_inward(i1 + isrc))
                  xb(dst, k, i, mdst) = xb(dst, k, i, mdst) - xx(src, k, msrc, j) * v(src, dst)
                  src = src + 1
                end do
                do isrc = 0, solver%diffside%dof - 1
                  msrc = merge(j, j + 1, solver%diffside%is_inward(i1 + isrc))
                  xb(dst, k, i, mdst) = xb(dst, k, i, mdst) - xx(src, k, i, msrc) * v(src, dst)
                  src = src + 1
                end do
                dst = dst + 1
              end do

            end if
          end do

          ! Albedo:
          do idst = 0, solver%difftop%dof - 1
            if (.not. solver%difftop%is_inward(i1 + idst)) then ! Eup
              xb(idst, C%ze, i, j) = xb(idst, C%ze, i, j) - xx(inv_dof(idst), C%ze, i, j) * solver%atm%albedo(i, j)
            end if
          end do

        end do
      end do

      call restoreVecPointer(C%da, lb, xb1d, xb)
      call restoreVecPointer(C%da, lx, xx1d, xx, readonly=.true.)
      call DMRestoreLocalVector(C%da, lx, ierr); call CHKERR(ierr)

      call VecSet(b, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(C%da, lb, ADD_VALUES, b, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd(C%da, lb, ADD_VALUES, b, ierr); call CHKERR(ierr)

      ! vec copy simulates diagonal entries
      call VecAXPY(b, one, x, ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
    end associate
    call PetscObjectViewFromOptions(PetscObjectCast(b), PETSC_NULL_OBJECT, "-show_shell_rb", ierr); call CHKERR(ierr)

    ierr = 0
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

  !> @brief define the matsor function for diff2diff computations (used for MatShell)
  subroutine op_mat_sor_ediff(A, b, omega, sortype, shift, its, lits, x, ierr)
    type(tMat), intent(in) :: A
    type(tVec), intent(in) :: b
    real(ireals), intent(in) :: omega
    MatSORType, intent(in) :: sortype
    real(ireals), intent(in) :: shift
    integer(iintegers), intent(in) :: its
    integer(iintegers), intent(in) :: lits
    type(tVec), intent(inout) :: x
    integer(mpiint), intent(out) :: ierr

    type(t_pprts_shell_ctx), pointer :: ctx_ptr
    class(t_solver), pointer :: solver
    type(tVec) :: lb, lx

    real(ireals), pointer, dimension(:, :, :, :) :: xx => null(), xg => null()
    real(ireals), pointer, dimension(:) :: xx1d => null(), xg1d => null()

    integer(iintegers) :: i_its, i_lits
    logical :: lsweep_forward, lsweep_backward, lsweep_symmetric, lzero_initial_guess

    nullify (ctx_ptr)
    call MatShellGetContext(A, ctx_ptr, ierr); call CHKERR(ierr)

    nullify (solver)
    solver => ctx_ptr%solver
    if (.not. associated(solver)) call CHKERR(1_mpiint, 'mat shell context has not been set!')

    if (.not. approx(shift, zero)) call CHKERR(1_mpiint, 'Shift has to be zero')
    if (.not. approx(omega, one)) call CHKERR(1_mpiint, 'Omega has to be one')

    lsweep_symmetric = iand(sortype%v, SOR_LOCAL_SYMMETRIC_SWEEP%v) .eq. SOR_LOCAL_SYMMETRIC_SWEEP%v .or. &
      & iand(sortype%v, SOR_SYMMETRIC_SWEEP%v) .eq. SOR_SYMMETRIC_SWEEP%v
    lsweep_backward = iand(sortype%v, SOR_LOCAL_BACKWARD_SWEEP%v) .eq. SOR_LOCAL_BACKWARD_SWEEP%v .or. &
      & iand(sortype%v, SOR_BACKWARD_SWEEP%v) .eq. SOR_BACKWARD_SWEEP%v
    lsweep_forward = iand(sortype%v, SOR_LOCAL_FORWARD_SWEEP%v) .eq. SOR_LOCAL_FORWARD_SWEEP%v .or. &
      & iand(sortype%v, SOR_FORWARD_SWEEP%v) .eq. SOR_FORWARD_SWEEP%v
    lzero_initial_guess = iand(sortype%v, SOR_ZERO_INITIAL_GUESS%v) .eq. SOR_ZERO_INITIAL_GUESS%v

    associate ( &
        & C => solver%C_diff)

      call DMGetLocalVector(C%da, lx, ierr); call CHKERR(ierr)
      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call DMGlobalToLocal(C%da, b, INSERT_VALUES, lb, ierr); call CHKERR(ierr)

      if (lzero_initial_guess) then
        call VecSet(lx, zero, ierr); call CHKERR(ierr)
      else
        call DMGlobalToLocal(C%da, x, INSERT_VALUES, lx, ierr); call CHKERR(ierr)
      end if

      do i_its = 1, its

        do i_lits = 1, lits
          if (lsweep_forward .or. lsweep_symmetric) then
            call explicit_ediff_sor_sweep(&
              & solver, &
              & solver%diff2diff, &
              & dx=[C%xs, C%xe, i1], &
              & dy=[C%ys, C%ye, i1], &
              & dz=[C%zs, C%ze - 1, i1], &
              & omega=omega, &
              & b=lb, x=lx)
          end if
          if (lsweep_backward .or. lsweep_symmetric) then
            call explicit_ediff_sor_sweep(&
              & solver, &
              & solver%diff2diff, &
              & dx=[C%xe, C%xs, -i1], &
              & dy=[C%ye, C%ys, -i1], &
              & dz=[C%ze - 1, C%zs, -i1], &
              & omega=omega, &
              & b=lb, x=lx)
          end if
        end do

        call exchange_diffuse_boundary(solver, lx, ierr); call CHKERR(ierr)

        call getVecPointer(C%da, x, xg1d, xg)
        call getVecPointer(C%da, lx, xx1d, xx, readonly=.true.)

        ! update solution vec
        xg = xx(:, :, C%xs:C%xe, C%ys:C%ye)

        call restoreVecPointer(C%da, x, xg1d, xg)
        call restoreVecPointer(C%da, lx, xx1d, xx, readonly=.true.)
      end do

      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call DMRestoreLocalVector(C%da, lx, ierr); call CHKERR(ierr)
    end associate
    ierr = 0
  end subroutine

end module
