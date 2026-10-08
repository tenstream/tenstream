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

  use mpi

  use m_data_parameters, only: &
    & i0, i1, &
    & iintegers, mpiint, &
    & ireals, imp_ireals, &
    & zero, one

  use m_helper_functions, only: &
    & approx, &
    & CHKERR, &
    & get_petsc_opt, &
    & imp_allreduce_mean, &
    & mpi_logical_and, &
    & toStr

  use m_pprts_base, only: &
    & atmk, &
    & determine_ksp_tolerances, &
    & t_coord, &
    & t_solver, &
    & t_state_container

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

  subroutine explicit_edir(solver, prefix, edirTOA, vedir, lb, v0, solution, ierr)
    class(t_solver), target, intent(in) :: solver
    character(len=*), intent(in) :: prefix
    real(ireals), intent(in) :: edirTOA
    real(ireals), target, contiguous, intent(inout) :: vedir(:, :, :, :) ! (0:dof-1, zs:ze, xs:xe, ys:ye)
    real(ireals), target, contiguous, intent(in) :: lb(:, :, :, :) ! incSolar, ghost-extended
    real(ireals), target, contiguous, intent(inout) :: v0(:, :, :, :) ! working copy, ghost-extended
    type(t_state_container), target, intent(inout) :: solution
    integer(mpiint), intent(out) :: ierr

    real(ireals), pointer, dimension(:, :, :, :) :: x0, xg

    integer(iintegers), dimension(3) :: dx, dy ! start, end, increment for each dimension
    integer(iintegers) :: iter, maxiter, maxit_ignore

    real(ireals), allocatable :: residual(:)
    real(ireals) :: residual_mean, rel_residual, atol, rtol, atol_default, rtol_default

    integer(iintegers), parameter :: default_max_it = 1000
    real(ireals) :: ignore_max_it ! Ignore max iter setting if time is less

    logical :: lksp_view, lcomplete_initial_run
    logical :: lsun_north, lsun_east, lpermute, lskip_residual, lmonitor_residual, lflg
    logical :: laccept_incomplete_solve, lconverged_atol, lconverged_rtol, lconverged_reason

    x0 => null()
    xg => null()
    ierr = 0

    associate ( &
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
      call get_petsc_opt('', "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(prefix, "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)

      lsun_north = sun%yinc .eq. i0
      lsun_east = sun%xinc .eq. i0

      dx = [C%xs, C%xe, i1]
      dy = [C%ys, C%ye, i1]

      lpermute = .true.
      call get_petsc_opt('', "-explicit_edir_permute", lpermute, lflg, ierr); call CHKERR(ierr)
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
        print *, '  -'//trim(prefix)//'explicit_edir_permute '//toStr(lpermute)// &
          & ' horizontal-xy-iterator: ['//toStr(dx(3))//', '//toStr(dy(3))//']'
      end if

      ! warm-start ghost from previous solution before the iteration loop
      call exchange_direct_boundary(solver, lsun_north, lsun_east, v0, ierr); call CHKERR(ierr)

      ! 2D open boundaries iterate on the inflow, start with the inflow of the single column solves
      if (solver%lopen_bc .and. solver%lopen_bc_2d) call set_open_bc_inflow(solver, lb, v0)

      do iter = 1, maxiter

        call explicit_edir_forward_sweep(solver, solver%dir2dir, dx, dy, lb, v0)

        call exchange_direct_boundary(solver, lsun_north, lsun_east, v0, ierr); call CHKERR(ierr)

        ! Residual computations
        if (.not. lskip_residual) then
          xg(0:C%dof - 1, C%zs:C%ze, C%xs:C%xe, C%ys:C%ye) => vedir
          x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => v0
          residual(iter) = max(tiny(one), norm2(xg - x0(:, :, C%xs:C%xe, C%ys:C%ye)))
          xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)
          nullify (xg, x0)

          call imp_allreduce_mean(solver%comm, residual(iter), residual_mean)
          residual(iter) = residual_mean

          if (residual(1) .le. sqrt(tiny(residual))) then
            rel_residual = 0
          else
            rel_residual = residual(iter) / residual(1)
          end if

          if (solver%myid .eq. 0 .and. lmonitor_residual) then
            print *, trim(prefix)//" iter "//toStr(iter)//' residual', residual_mean, 'rel res', rel_residual
          end if
          solution%dir_ksp_residual_history(min(size(solution%dir_ksp_residual_history, kind=iintegers), iter)) = residual(iter)

          lconverged_atol = residual(iter) .lt. atol
          lconverged_rtol = rel_residual .lt. rtol

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
      xg(0:C%dof - 1, C%zs:C%ze, C%xs:C%xe, C%ys:C%ye) => vedir
      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => v0
      xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)
      nullify (xg, x0)
    end associate

    solution%lchanged = .true.
    solution%lWm2_dir = .false.

  end subroutine

  subroutine exchange_direct_boundary(solver, lsun_north, lsun_east, x, ierr)
    class(t_solver), intent(in) :: solver
    logical, intent(in) :: lsun_north, lsun_east
    real(ireals), target, contiguous, intent(inout) :: x(:, :, :, :)
    integer(mpiint), intent(out) :: ierr

    real(ireals), pointer :: x0(:, :, :, :)

    integer(mpiint), parameter :: tag_x = 1, tag_y = 2
    integer(mpiint) :: neigh_s, neigh_r, requests(4), statuses(mpi_status_size, 4)
    integer(iintegers) :: dofstart, dofend

    real(ireals), allocatable :: mpi_send_bfr_x(:, :, :), mpi_send_bfr_y(:, :, :)
    real(ireals), allocatable :: mpi_recv_bfr_x(:, :, :), mpi_recv_bfr_y(:, :, :)

    x0 => null()

    associate ( &
        & C => solver%C_dir)

      ! including the ghost rows/columns, with open boundaries they hold the ghost cells outside of the sunward edges
      allocate (mpi_send_bfr_x(solver%dirside%dof, C%zm, C%gys:C%gye), mpi_recv_bfr_x(solver%dirside%dof, C%zm, C%gys:C%gye))
      allocate (mpi_send_bfr_y(solver%dirside%dof, C%zm, C%gxs:C%gxe), mpi_recv_bfr_y(solver%dirside%dof, C%zm, C%gxs:C%gxe))

      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => x

      ! Boundary exchanges
      ! x direction scatters
      dofstart = solver%dirtop%dof
      dofend = -1 + solver%dirtop%dof + solver%dirside%dof
      if (lsun_east) then
        if (solver%lopen_bc_x .and. C%xs .eq. i0) then
          mpi_send_bfr_x = 0
        else
          mpi_send_bfr_x = x0(dofstart:dofend, :, C%xs, C%gys:C%gye)
        end if
        neigh_s = int(C%neighbors(10), mpiint) ! neigh west
        neigh_r = int(C%neighbors(16), mpiint) ! neigh east
      else
        if (solver%lopen_bc_x .and. C%xe + 1 .eq. C%glob_xm) then
          mpi_send_bfr_x = 0
        else
          mpi_send_bfr_x = x0(dofstart:dofend, :, C%xe + 1, C%gys:C%gye)
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
        if (solver%lopen_bc_y .and. C%ys .eq. i0) then
          mpi_send_bfr_y = 0
        else
          mpi_send_bfr_y = x0(dofstart:dofend, :, C%gxs:C%gxe, C%ys)
        end if
        neigh_s = int(C%neighbors(4), mpiint) ! neigh south
        neigh_r = int(C%neighbors(22), mpiint) ! neigh north
      else
        if (solver%lopen_bc_y .and. C%ye + 1 .eq. C%glob_ym) then
          mpi_send_bfr_y = 0
        else
          mpi_send_bfr_y = x0(dofstart:dofend, :, C%gxs:C%gxe, C%ye + 1)
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
      ! with open boundaries, the inflow face at the domain edge holds the boundary condition, dont overwrite it.
      ! With -pprts_open_bc_2d, the ghost cells outside of the domain edge set it in the sweep
      if (lsun_east) then
        if (.not. (solver%lopen_bc_x .and. C%xe + 1 .eq. C%glob_xm)) then
          x0(dofstart:dofend, :, C%xe + 1, C%gys:C%gye) = mpi_recv_bfr_x
        end if
      else
        if (.not. (solver%lopen_bc_x .and. C%xs .eq. i0)) then
          x0(dofstart:dofend, :, C%xs, C%gys:C%gye) = mpi_recv_bfr_x
        end if
      end if

      dofstart = solver%dirtop%dof + solver%dirside%dof
      dofend = -1 + solver%dirtop%dof + solver%dirside%dof * 2
      if (lsun_north) then
        if (.not. (solver%lopen_bc_y .and. C%ye + 1 .eq. C%glob_ym)) then
          x0(dofstart:dofend, :, C%gxs:C%gxe, C%ye + 1) = mpi_recv_bfr_y
        end if
      else
        if (.not. (solver%lopen_bc_y .and. C%ys .eq. i0)) then
          x0(dofstart:dofend, :, C%gxs:C%gxe, C%ys) = mpi_recv_bfr_y
        end if
      end if
      nullify (x0)
    end associate
    ierr = 0
  end subroutine

  !> @brief copy the inflow at the sunward open domain edges from the source term b into x
  subroutine set_open_bc_inflow(solver, b, x)
    class(t_solver), intent(in) :: solver
    real(ireals), target, contiguous, intent(in) :: b(:, :, :, :)
    real(ireals), target, contiguous, intent(inout) :: x(:, :, :, :)

    real(ireals), pointer :: x0(:, :, :, :), xb(:, :, :, :)
    logical :: lsun_north, lsun_east

    associate ( &
        & C => solver%C_dir, &
        & xinc => solver%sun%xinc, &
        & yinc => solver%sun%yinc)

      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => x
      xb(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => b

      lsun_north = yinc .eq. i0
      lsun_east = xinc .eq. i0

      if (solver%lopen_bc_y .and. lsun_north .and. C%ye + 1 .eq. C%glob_ym) then
        x0(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ye + 1) = &
          & xb(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ye + 1)
      end if
      if (solver%lopen_bc_y .and. .not. lsun_north .and. C%ys .eq. i0) then
        x0(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ys) = &
          & xb(solver%dirtop%dof + solver%dirside%dof:C%dof - 1, :, C%xs:C%xe, C%ys)
      end if
      if (solver%lopen_bc_x .and. lsun_east .and. C%xe + 1 .eq. C%glob_xm) then
        x0(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xe + 1, C%ys:C%ye) = &
          & xb(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xe + 1, C%ys:C%ye)
      end if
      if (solver%lopen_bc_x .and. .not. lsun_east .and. C%xs .eq. i0) then
        x0(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xs, C%ys:C%ye) = &
          & xb(solver%dirtop%dof:solver%dirtop%dof + solver%dirside%dof - 1, :, C%xs, C%ys:C%ye)
      end if
      nullify (x0, xb)
    end associate
  end subroutine

  subroutine explicit_edir_forward_sweep(solver, coeffs, dx, dy, b, x)
    class(t_solver), intent(in) :: solver
    real(ireals), target, intent(in) :: coeffs(:, :, :, :)
    integer(iintegers), dimension(3), intent(in) :: dx, dy ! start, end, increment for each dimension
    real(ireals), target, contiguous, intent(in) :: b(:, :, :, :)
    real(ireals), target, contiguous, intent(inout) :: x(:, :, :, :)

    real(ireals), pointer :: x0(:, :, :, :), xb(:, :, :, :)

    integer(iintegers) :: i, j, k
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:, :) ! dim(src, dst)
    logical :: lghost_x, lghost_y, lghost_corner
    integer(iintegers) :: ig, jg, ie, je

    x0 => null()
    xb => null()

    associate ( &
        & atm => solver%atm, &
        & C => solver%C_dir, &
        & xinc => solver%sun%xinc, &
        & yinc => solver%sun%yinc)

      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => x
      xb(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => b

      x0(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, C%ys:C%ye) = xb(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, C%ys:C%ye)

      if (solver%lopen_bc .and. .not. solver%lopen_bc_2d) call set_open_bc_inflow(solver, b, x)
      nullify (xb)

      ! ghost cells outside of the sunward open domain edges
      lghost_x = solver%lopen_bc_2d .and. allocated(solver%dir2dir_ghost_x)
      lghost_y = solver%lopen_bc_2d .and. allocated(solver%dir2dir_ghost_y)
      ig = merge(C%xs - 1, C%xe + 1, xinc .eq. i1) ! ghost column
      jg = merge(C%ys - 1, C%ye + 1, yinc .eq. i1) ! ghost row
      ie = merge(C%xs, C%xe, xinc .eq. i1) ! edge column next to it
      je = merge(C%ys, C%ye, yinc .eq. i1)
      lghost_corner = solver%lopen_bc_2d .and. allocated(solver%dir2dir_ghost_corner)
      if (lghost_x) x0(0:solver%dirtop%dof - 1, C%zs, ig, C%ys:C%ye) = x0(0:solver%dirtop%dof - 1, C%zs, ie, C%ys:C%ye)
      if (lghost_y) x0(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, jg) = x0(0:solver%dirtop%dof - 1, C%zs, C%xs:C%xe, je)
      if (lghost_corner) x0(0:solver%dirtop%dof - 1, C%zs, ig, jg) = x0(0:solver%dirtop%dof - 1, C%zs, ie, je)

      ! forward sweep through x
      do k = C%zs, C%ze - 1
        if (atm%l1d(atmk(atm, k))) then
          if (lghost_x) then
            do j = C%ys, C%ye
              x0(0:solver%dirtop%dof - 1, k + i1, ig, j) = x0(0:solver%dirtop%dof - 1, k, ig, j) * atm%a33(atmk(atm, k), ie, j)
            end do
          end if
          if (lghost_y) then
            do i = C%xs, C%xe
              x0(0:solver%dirtop%dof - 1, k + i1, i, jg) = x0(0:solver%dirtop%dof - 1, k, i, jg) * atm%a33(atmk(atm, k), i, je)
            end do
          end if
          if (lghost_corner) then
            x0(0:solver%dirtop%dof - 1, k + i1, ig, jg) = x0(0:solver%dirtop%dof - 1, k, ig, jg) * atm%a33(atmk(atm, k), ie, je)
          end if
          do j = dy(1), dy(2), dy(3)
            do i = dx(1), dx(2), dx(3)

              do idst = 0, solver%dirtop%dof - 1
                x0(idst, k + i1, i, j) = x0(idst, k, i, j) * atm%a33(atmk(atm, k), i, j)
              end do
            end do
          end do
        else
          ! the ghost cells: zero gradient across the domain edge, i.e. what enters a ghost cell from the outside
          ! is what leaves it towards the domain. Along the edge, they see their neighbours,
          ! in the sunward corner of the domain this is the corner ghost cell
          if (lghost_corner) then
            call ghost_cell(solver, solver%dir2dir_ghost_corner(:, k), x0, k, ig, jg, ig + xinc, jg + yinc)
          end if
          if (lghost_x) then
            do j = dy(1), dy(2), dy(3)
              call ghost_cell(solver, solver%dir2dir_ghost_x(:, k, j), x0, k, ig, j, ig + xinc, j + 1 - yinc)
            end do
          end if
          if (lghost_y) then
            do i = dx(1), dx(2), dx(3)
              call ghost_cell(solver, solver%dir2dir_ghost_y(:, k, i), x0, k, i, jg, i + 1 - xinc, jg + yinc)
            end do
          end if
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
      nullify (x0)
    end associate

  end subroutine

  !> @brief transport through a ghost cell i,j of the explicit direct sweep
  !> @details side inflow is read from faces ix (x) and iy (y). If a face is also the outflow face of the ghost cell
  !> (zero gradient across a domain edge), inflow and outflow are equal and we solve for them directly
  subroutine ghost_cell(solver, coeff, x0, k, i, j, ix, iy)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(in) :: coeff(:)
    real(ireals), intent(inout) :: x0(0:, solver%C_dir%zs:, solver%C_dir%gxs:, solver%C_dir%gys:)
    integer(iintegers), intent(in) :: k, i, j, ix, iy

    real(ireals) :: c(solver%C_dir%dof, solver%C_dir%dof) ! dim(src, dst)
    real(ireals) :: xin(solver%C_dir%dof), xout(solver%C_dir%dof)
    real(ireals), allocatable :: A(:, :), rhs(:)
    logical :: lself(solver%C_dir%dof)
    integer(iintegers), allocatable :: idx(:)
    integer(iintegers) :: ntop, nside, dof, xinc, yinc, n, m, p, q

    ntop = solver%dirtop%dof
    nside = solver%dirside%dof
    dof = solver%C_dir%dof
    xinc = solver%sun%xinc
    yinc = solver%sun%yinc
    c = reshape(coeff, [dof, dof])

    lself = .false.
    lself(ntop + 1:ntop + nside) = ix .eq. i + xinc
    lself(ntop + nside + 1:dof) = iy .eq. j + yinc

    xin(1:ntop) = x0(0:ntop - 1, k, i, j)
    xin(ntop + 1:ntop + nside) = x0(ntop:ntop + nside - 1, k, ix, j)
    xin(ntop + nside + 1:dof) = x0(ntop + nside:dof - 1, k, i, iy)

    if (any(lself)) then
      ! outflow of the self coupled streams: (I - C_ss^T) o_s = C_ns^T x_n
      idx = pack([(m, m=1, dof)], lself)
      n = size(idx)
      allocate (A(n, n), rhs(n))
      do p = 1, n
        rhs(p) = zero
        do m = 1, dof
          if (.not. lself(m)) rhs(p) = rhs(p) + xin(m) * c(m, idx(p))
        end do
        do q = 1, n
          A(p, q) = -c(idx(q), idx(p))
        end do
        A(p, p) = A(p, p) + one
      end do
      call solve_small(A, rhs)
      xin(idx) = rhs
    end if

    xout = matmul(xin, c)
    x0(0:ntop - 1, k + i1, i, j) = xout(1:ntop)
    x0(ntop:ntop + nside - 1, k, i + xinc, j) = xout(ntop + 1:ntop + nside)
    x0(ntop + nside:dof - 1, k, i, j + yinc) = xout(ntop + nside + 1:dof)
  contains
    !> gaussian elimination, the system is diagonally dominant
    subroutine solve_small(A, b)
      real(ireals), intent(inout) :: A(:, :), b(:)
      integer(iintegers) :: r, s
      real(ireals) :: f
      do r = 1, size(b)
        do s = r + 1, size(b)
          f = A(s, r) / A(r, r)
          A(s, r:) = A(s, r:) - f * A(r, r:)
          b(s) = b(s) - f * b(r)
        end do
      end do
      do r = size(b), 1, -1
        b(r) = (b(r) - dot_product(A(r, r + 1:), b(r + 1:))) / A(r, r)
      end do
    end subroutine
  end subroutine

  subroutine explicit_ediff(solver, prefix, vb, vediff, solution, ierr)
    class(t_solver), intent(inout) :: solver
    character(len=*), intent(in) :: prefix
    real(ireals), target, intent(in) :: vb(:, :, :, :)
    real(ireals), target, contiguous, intent(inout) :: vediff(:, :, :, :) ! (0:dof-1, zs:ze, xs:xe, ys:ye)
    type(t_state_container), target, intent(inout) :: solution
    integer(mpiint), intent(out) :: ierr

    real(ireals), pointer, dimension(:, :, :, :) :: x0, xg
    real(ireals), allocatable, target :: lvb(:, :, :, :), v0(:, :, :, :)

    integer(iintegers) :: iter, isub, maxiter, sub_iter, maxit_ignore

    real(ireals), allocatable :: residual(:)
    real(ireals) :: residual_mean, rel_residual, atol, rtol, atol_default, rtol_default

    integer(iintegers), parameter :: default_max_it = 10000
    real(ireals) :: ignore_max_it ! Ignore max iter setting if time is less

    logical :: lksp_view, lcomplete_initial_run
    logical :: lskip_residual, lmonitor_residual, lconverged_atol, lconverged_rtol, lflg
    logical :: laccept_incomplete_solve, lconverged_reason

    logical :: lomega_set, ladaptive_omega, lomega_frozen
    real(ireals) :: omega, omega_adaptive, omega_increment, omega_min, omega_max, omega_stagnation
    real(ireals) :: omega_dir, omega_step, omega_save, log_rate, log_rate_prev
    real(ireals) :: best_residual_solve
    integer(iintegers) :: iter_at_best
    integer(iintegers), parameter :: stagnation_window = 50

    x0 => null()
    xg => null()
    ierr = 0

    associate ( &
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
      call get_petsc_opt('', "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)
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
      omega_stagnation = 0.5_ireals
      call get_petsc_opt(prefix, "-pc_sor_omega_stagnation", omega_stagnation, lflg, ierr); call CHKERR(ierr)
      omega_dir = 1._ireals
      omega_step = omega_increment * 0.5_ireals
      log_rate_prev = 0._ireals
      if (ladaptive_omega .and. solution%diff_sor_omega .gt. zero) then
        omega_adaptive = solution%diff_sor_omega ! warm-start from previous solve for this g-point
      end if
      omega_save = omega_adaptive
      best_residual_solve = huge(best_residual_solve)
      iter_at_best = 1
      lomega_frozen = .false.

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

      allocate (v0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye)); v0 = zero
      v0(:, :, C%xs:C%xe, C%ys:C%ye) = vediff
      ! warm-start ghost from previous solution
      call fill_ghost(solver%comm, C, v0, ierr); call CHKERR(ierr)
      call set_open_bc_diffuse(solver, v0)

      allocate (lvb(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye)); lvb = zero
      lvb(:, :, C%xs:C%xe, C%ys:C%ye) = vb
      ! RHS ghost cells are used by the sweep for incoming diffuse at boundaries
      call fill_ghost(solver%comm, C, lvb, ierr); call CHKERR(ierr)

      sorloop: do iter = 1, maxiter
        do isub = 1, sub_iter
          ! ghost cells outside of the open domain edges, they feed the domain
          if (solver%lopen_bc .and. solver%lopen_bc_2d) call diffuse_ghost_sweep(solver, lvb, v0)
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
        end do

        call exchange_diffuse_boundary(solver, v0, ierr); call CHKERR(ierr)

        ! Residual computations
        if (.not. lskip_residual) then
          xg(0:C%dof - 1, C%zs:C%ze, C%xs:C%xe, C%ys:C%ye) => vediff
          x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => v0

          residual(iter) = norm2(xg - x0(:, :, C%xs:C%xe, C%ys:C%ye))
          xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)

          nullify (xg, x0)

          call imp_allreduce_mean(solver%comm, residual(iter), residual_mean)
          residual(iter) = residual_mean

          if (residual(1) .le. sqrt(tiny(residual))) then
            rel_residual = 0
          else
            rel_residual = residual(iter) / residual(1)
          end if

          if (solver%myid .eq. 0 .and. lmonitor_residual) then
            print *, trim(prefix), ' iter ', toStr(iter), &
              & ' residual', residual_mean, &
              & 'rel res', rel_residual, &
              & 'omega', omega_adaptive
          end if
          solution%diff_ksp_residual_history(min(size(solution%diff_ksp_residual_history, kind=iintegers), iter)) = residual(iter)

          lconverged_atol = residual(iter) .lt. atol
          lconverged_rtol = rel_residual .lt. rtol

          if (lconverged_atol .or. lconverged_rtol) then
            if (solver%myid .eq. 0 .and. lconverged_reason) then
              if (lconverged_atol) then
                print *, trim(prefix)//' solve converged due to CONVERGED_ATOL iterations', iter
              else
                print *, trim(prefix)//' solve converged due to CONVERGED_RTOL iterations', iter
              end if
            end if
            exit sorloop
          end if

          if (residual(iter) .lt. best_residual_solve) then
            best_residual_solve = residual(iter)
            iter_at_best = iter
          end if

          if (ladaptive_omega .and. iter .ge. 3) then
            ! Stagnation guard: if no new best for stagnation_window iters while omega > 1,
            ! under-relax (omega < 1) for the rest of the solve to escape a limit cycle
            if (.not. lomega_frozen &
              & .and. omega_adaptive .gt. omega_min &
              & .and. (iter - iter_at_best) .gt. stagnation_window) then
              omega_adaptive = omega_stagnation
              omega_save = omega_adaptive
              lomega_frozen = .true.
            end if
            if (.not. lomega_frozen) then
              if (residual(iter) .gt. zero .and. residual(iter - 2) .gt. zero) then
                log_rate = 0.5_ireals * log(residual(iter) / residual(iter - 2))
                if (log_rate .lt. log_rate_prev) then
                  omega_step = min(omega_step * 1.3_ireals, omega_max - omega_min)
                else
                  omega_dir = -omega_dir
                  omega_step = max(omega_step * 0.5_ireals, 0.01_ireals)
                end if
                log_rate_prev = log_rate
                omega_adaptive = min(max(omega_adaptive + omega_dir * omega_step, omega_min), omega_max)
                omega_save = omega_adaptive
              end if
            end if
          end if

        else
          solution%diff_ksp_residual_history(min(size(solution%diff_ksp_residual_history, kind=iintegers), iter)) = zero
        end if

        if (iter .eq. maxiter) then
          if (.not. laccept_incomplete_solve) then
            call CHKERR(int(iter, mpiint), trim(prefix)//" did not converge")
          end if
        end if
      end do sorloop

      if (ladaptive_omega) solution%diff_sor_omega = omega_save

      xg(0:C%dof - 1, C%zs:C%ze, C%xs:C%xe, C%ys:C%ye) => vediff
      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => v0

      ! update solution vec
      xg = x0(:, :, C%xs:C%xe, C%ys:C%ye)

      ! keep the fluxes through the east/north open boundary, they live on a ghost face and are not part of ediff.
      ! A periodic halo exchange would replace them with the fluxes through the opposite edge
      if (allocated(solution%ediff_open_bc_x)) deallocate (solution%ediff_open_bc_x)
      if (allocated(solution%ediff_open_bc_y)) deallocate (solution%ediff_open_bc_y)
      if (solver%lopen_bc_x .and. C%xe + 1 .eq. C%glob_xm) then
        allocate (solution%ediff_open_bc_x(0:C%dof - 1, C%zs:C%ze, C%ys:C%ye), source=x0(:, :, C%xe + 1, C%ys:C%ye))
      end if
      if (solver%lopen_bc_y .and. C%ye + 1 .eq. C%glob_ym) then
        allocate (solution%ediff_open_bc_y(0:C%dof - 1, C%zs:C%ze, C%xs:C%xe), source=x0(:, :, C%xs:C%xe, C%ye + 1))
      end if

      nullify (xg, x0)
      deallocate (v0, lvb)
    end associate

    solution%lchanged = .true.
    solution%lWm2_diff = .false.
  end subroutine

  subroutine exchange_diffuse_boundary(solver, v0, ierr)
    class(t_solver), intent(in) :: solver
    real(ireals), target, contiguous, intent(inout) :: v0(:, :, :, :)
    integer(mpiint), intent(out) :: ierr

    integer(mpiint), parameter :: tag_e = 1, tag_w = 2, tag_n = 3, tag_s = 4
    integer(mpiint) :: neigh_s, neigh_e, neigh_w, neigh_n
    integer(mpiint) :: requests(8), statuses(mpi_status_size, 8)
    integer(iintegers) :: k, i, j, d1, d2, dof, idof

    real(ireals), pointer :: x0(:, :, :, :)

    real(ireals), allocatable :: mpi_send_bfr_e(:, :, :), mpi_send_bfr_n(:, :, :)
    real(ireals), allocatable :: mpi_send_bfr_w(:, :, :), mpi_send_bfr_s(:, :, :)
    real(ireals), allocatable :: mpi_recv_bfr_e(:, :, :), mpi_recv_bfr_n(:, :, :)
    real(ireals), allocatable :: mpi_recv_bfr_w(:, :, :), mpi_recv_bfr_s(:, :, :)

    logical :: lopen_w, lopen_e, lopen_s, lopen_n

    x0 => null()

    associate ( &
        & C => solver%C_diff)

      lopen_w = solver%lopen_bc_x .and. C%xs .eq. i0
      lopen_e = solver%lopen_bc_x .and. C%xe + 1 .eq. C%glob_xm
      lopen_s = solver%lopen_bc_y .and. C%ys .eq. i0
      lopen_n = solver%lopen_bc_y .and. C%ye + 1 .eq. C%glob_ym

      allocate ( &
        & mpi_send_bfr_e(solver%diffside%dof / 2, C%zs:C%ze, C%gys:C%gye), &
        & mpi_send_bfr_w(solver%diffside%dof / 2, C%zs:C%ze, C%gys:C%gye), &
        & mpi_send_bfr_n(solver%diffside%dof / 2, C%zs:C%ze, C%gxs:C%gxe), &
        & mpi_send_bfr_s(solver%diffside%dof / 2, C%zs:C%ze, C%gxs:C%gxe), &
        & mpi_recv_bfr_e(solver%diffside%dof / 2, C%zs:C%ze, C%gys:C%gye), &
        & mpi_recv_bfr_w(solver%diffside%dof / 2, C%zs:C%ze, C%gys:C%gye), &
        & mpi_recv_bfr_n(solver%diffside%dof / 2, C%zs:C%ze, C%gxs:C%gxe), &
        & mpi_recv_bfr_s(solver%diffside%dof / 2, C%zs:C%ze, C%gxs:C%gxe) &
        & )
      mpi_recv_bfr_e = zero; mpi_recv_bfr_w = zero
      mpi_recv_bfr_n = zero; mpi_recv_bfr_s = zero

      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => v0

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
      do j = C%gys, C%gye
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
      do i = C%gxs, C%gxe
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
      do j = C%gys, C%gye
        do k = C%zs, C%ze
          d1 = 1; d2 = 1
          do idof = i0, solver%diffside%dof - 1
            dof = solver%difftop%dof + idof
            ! the inflow at open domain edges is set by the open boundary condition
            if (solver%diffside%is_inward(i1 + idof)) then ! to the right
              if (.not. lopen_w) x0(dof, k, C%xs, j) = mpi_recv_bfr_w(d1, k, j)
              d1 = d1 + 1
            else ! leftward
              if (.not. lopen_e) x0(dof, k, C%xe + 1, j) = mpi_recv_bfr_e(d2, k, j)
              d2 = d2 + 1
            end if
          end do
        end do
      end do

      do i = C%gxs, C%gxe
        do k = C%zs, C%ze
          d1 = 1; d2 = 1
          do idof = i0, solver%diffside%dof - 1
            dof = solver%difftop%dof + solver%diffside%dof + idof
            if (solver%diffside%is_inward(i1 + idof)) then
              if (.not. lopen_s) x0(dof, k, i, C%ys) = mpi_recv_bfr_s(d1, k, i)
              d1 = d1 + 1
            else
              if (.not. lopen_n) x0(dof, k, i, C%ye + 1) = mpi_recv_bfr_n(d2, k, i)
              d2 = d2 + 1
            end if
          end do
        end do
      end do

      nullify (x0)
    end associate

    call set_open_bc_diffuse(solver, v0)
    ierr = 0
  end subroutine

  !> @brief open boundaries for diffuse radiation without ghost cells (-pprts_open_bc_2d no): zero gradient across the domain edges
  !> @details what enters an edge cell through the domain edge is what leaves this cell in the same direction,
  !> i.e. the edge columns continue outwards but they do see their neighbours along the edge
  subroutine set_open_bc_diffuse(solver, v0)
    class(t_solver), intent(in) :: solver
    real(ireals), target, contiguous, intent(inout) :: v0(:, :, :, :)

    real(ireals), pointer :: x0(:, :, :, :)
    integer(iintegers) :: idof, dof

    if (.not. solver%lopen_bc) return
    if (solver%lopen_bc_2d) return ! the ghost cells set the inflow, see diffuse_ghost_sweep

    associate (C => solver%C_diff)
      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => v0

      if (solver%lopen_bc_x) then
        do idof = i0, solver%diffside%dof - 1
          dof = solver%difftop%dof + idof
          if (solver%diffside%is_inward(i1 + idof)) then ! to the right, enters at the west edge
            if (C%xs .eq. i0) x0(dof, :, C%xs, C%ys:C%ye) = x0(dof, :, C%xs + 1, C%ys:C%ye)
          else ! leftward, enters at the east edge
            if (C%xe + 1 .eq. C%glob_xm) x0(dof, :, C%xe + 1, C%ys:C%ye) = x0(dof, :, C%xe, C%ys:C%ye)
          end if
        end do
      end if

      if (solver%lopen_bc_y) then
        do idof = i0, solver%diffside%dof - 1
          dof = solver%difftop%dof + solver%diffside%dof + idof
          if (solver%diffside%is_inward(i1 + idof)) then ! forward, enters at the south edge
            if (C%ys .eq. i0) x0(dof, :, C%xs:C%xe, C%ys) = x0(dof, :, C%xs:C%xe, C%ys + 1)
          else ! backward, enters at the north edge
            if (C%ye + 1 .eq. C%glob_ym) x0(dof, :, C%xs:C%xe, C%ye + 1) = x0(dof, :, C%xs:C%xe, C%ye)
          end if
        end do
      end if
      nullify (x0)
    end associate
  end subroutine

  !> @brief one pass through the diffuse ghost cells outside of all open domain edges and corners
  !> @details the ghost cells continue the edge columns outwards (see alloc_coeff_diff2diff_ghost) and use the sources
  !> of the edge cells. Across the domain edge, zero gradient: what enters a ghost cell from the outside is what it emits
  !> in the same direction towards the domain. Along the edge, ghost cells see their neighbours,
  !> at the corners the corner ghost cells.
  subroutine diffuse_ghost_sweep(solver, b, x)
    class(t_solver), intent(in) :: solver
    real(ireals), target, contiguous, intent(in) :: b(:, :, :, :)
    real(ireals), target, contiguous, intent(inout) :: x(:, :, :, :)

    real(ireals), pointer :: x0(:, :, :, :), xb(:, :, :, :)
    integer(iintegers) :: i, j

    associate (C => solver%C_diff)
      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => x
      xb(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => b

      if (allocated(solver%diff2diff_ghost_w)) then
        do j = C%ys, C%ye
          call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_w(:, :, j), C%xs - 1, j, i1, i0, 1_iintegers, &
                                    0_iintegers)
        end do
      end if
      if (allocated(solver%diff2diff_ghost_e)) then
        do j = C%ys, C%ye
          call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_e(:, :, j), C%xe + 1, j, -i1, i0, 2_iintegers, &
                                    0_iintegers)
        end do
      end if
      if (allocated(solver%diff2diff_ghost_s)) then
        do i = C%xs, C%xe
          call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_s(:, :, i), i, C%ys - 1, i0, i1, 0_iintegers, &
                                    1_iintegers)
        end do
      end if
      if (allocated(solver%diff2diff_ghost_n)) then
        do i = C%xs, C%xe
          call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_n(:, :, i), i, C%ye + 1, i0, -i1, 0_iintegers, &
                                    2_iintegers)
        end do
      end if
      if (allocated(solver%diff2diff_ghost_sw)) &
        & call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_sw, C%xs - 1, C%ys - 1, i1, i1, 1_iintegers, 1_iintegers)
      if (allocated(solver%diff2diff_ghost_se)) &
       & call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_se, C%xe + 1, C%ys - 1, -i1, i1, 2_iintegers, 1_iintegers)
      if (allocated(solver%diff2diff_ghost_nw)) &
       & call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_nw, C%xs - 1, C%ye + 1, i1, -i1, 1_iintegers, 2_iintegers)
      if (allocated(solver%diff2diff_ghost_ne)) &
      & call diffuse_ghost_column(solver, x0, xb, solver%diff2diff_ghost_ne, C%xe + 1, C%ye + 1, -i1, -i1, 2_iintegers, 2_iintegers)
      nullify (x0, xb)
    end associate

  end subroutine

  !> @brief diffuse ghost column i,j, (di,dj) points to the edge column,
  !> outer = 1: the outer face is the low face, 2: the high face, 0: none
  subroutine diffuse_ghost_column(solver, x0, xb, coeffs, i, j, di, dj, outer_x, outer_y)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(inout) :: x0(0:, solver%C_diff%zs:, solver%C_diff%gxs:, solver%C_diff%gys:)
    real(ireals), intent(in) :: xb(0:, solver%C_diff%zs:, solver%C_diff%gxs:, solver%C_diff%gys:)
    real(ireals), intent(in) :: coeffs(:, solver%C_diff%zs:)
    integer(iintegers), intent(in) :: i, j, di, dj, outer_x, outer_y
    integer(iintegers) :: k, idof

    associate (C => solver%C_diff, atm => solver%atm)
      do idof = 0, solver%difftop%dof - 1
        if (solver%difftop%is_inward(i1 + idof)) x0(idof, C%zs, i, j) = x0(idof, C%zs, i + di, j + dj)
      end do
      do k = C%zs, C%ze - 1
        call diffuse_ghost_cell(solver, x0, xb, coeffs(:, k), k, i, j, di, dj, outer_x, outer_y)
      end do
      do idof = 0, solver%difftop%dof - 1
        if (.not. solver%difftop%is_inward(i1 + idof)) then
          x0(idof, C%ze, i, j) = xb(idof, C%ze, i + di, j + dj) &
                             & + x0(diff_inv_dof(solver, idof), C%ze, i, j) * atm%albedo(i + di, j + dj)
        end if
      end do
      do k = C%ze - 1, C%zs, -1
        call diffuse_ghost_cell(solver, x0, xb, coeffs(:, k), k, i, j, di, dj, outer_x, outer_y)
      end do
    end associate
  end subroutine

  !> @brief diffuse ghost cell, see diffuse_ghost_column
  subroutine diffuse_ghost_cell(solver, x0, xb, coeff, k, i, j, di, dj, outer_x, outer_y)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(inout) :: x0(0:, solver%C_diff%zs:, solver%C_diff%gxs:, solver%C_diff%gys:)
    real(ireals), intent(in) :: xb(0:, solver%C_diff%zs:, solver%C_diff%gxs:, solver%C_diff%gys:)
    real(ireals), intent(in) :: coeff(:)
    integer(iintegers), intent(in) :: k, i, j, di, dj, outer_x, outer_y

    integer(iintegers) :: d, s, n, p, q, r, ntop, nside, ndof, ak
    integer(iintegers), allocatable, dimension(:) :: kin, iin, jin, kout, iout, jout, idx
    logical, allocatable :: lself(:)
    real(ireals), allocatable :: c(:, :), xin(:), xout(:), A(:, :), rhs(:)
    real(ireals) :: f

    ntop = solver%difftop%dof
    nside = solver%diffside%dof
    ndof = solver%C_diff%dof
    ak = atmk(solver%atm, k)

    if (solver%atm%l1d(ak)) then
      do d = 0, ntop - 1
        if (solver%difftop%is_inward(i1 + d)) then
          x0(d, k + 1, i, j) = xb(d, k + 1, i + di, j + dj) &
            & + x0(d, k, i, j) * solver%atm%a11(ak, i + di, j + dj) &
            & + x0(diff_inv_dof(solver, d), k + 1, i, j) * solver%atm%a12(ak, i + di, j + dj)
        else
          x0(d, k, i, j) = xb(d, k, i + di, j + dj) &
            & + x0(d, k + 1, i, j) * solver%atm%a11(ak, i + di, j + dj) &
            & + x0(diff_inv_dof(solver, d), k, i, j) * solver%atm%a12(ak, i + di, j + dj)
        end if
      end do
      return
    end if

    ! where each stream enters and leaves the ghost cell
    allocate (kin(0:ndof - 1), source=k)
    allocate (kout(0:ndof - 1), source=k)
    allocate (iin(0:ndof - 1), source=i)
    allocate (iout(0:ndof - 1), source=i)
    allocate (jin(0:ndof - 1), source=j)
    allocate (jout(0:ndof - 1), source=j)
    allocate (lself(0:ndof - 1), source=.false.)
    do d = 0, ntop - 1
      if (solver%difftop%is_inward(i1 + d)) then
        kout(d) = k + 1
      else
        kin(d) = k + 1
      end if
    end do
    do s = 0, nside - 1
      d = ntop + s
      if (solver%diffside%is_inward(i1 + s)) then ! to the right
        iout(d) = i + 1
        lself(d) = outer_x .eq. 1
      else
        iin(d) = i + 1
        lself(d) = outer_x .eq. 2
      end if
      d = ntop + nside + s
      if (solver%diffside%is_inward(i1 + s)) then ! forward
        jout(d) = j + 1
        lself(d) = outer_y .eq. 1
      else
        jin(d) = j + 1
        lself(d) = outer_y .eq. 2
      end if
    end do

    allocate (c(0:ndof - 1, 0:ndof - 1), xin(0:ndof - 1), xout(0:ndof - 1))
    c = reshape(coeff, [ndof, ndof]) ! dim(src, dst)
    xin = zero
    do d = 0, ndof - 1
      if (.not. lself(d)) xin(d) = x0(d, kin(d), iin(d), jin(d))
    end do

    ! streams that enter through the outer faces are equal to the ones that leave towards the domain
    n = count(lself)
    if (n .gt. 0) then
      allocate (idx(n), A(n, n), rhs(n))
      p = 0
      do d = 0, ndof - 1
        if (lself(d)) then
          p = p + 1
          idx(p) = d
        end if
      end do
      do p = 1, n
        rhs(p) = xb(idx(p), kout(idx(p)), iout(idx(p)) + di, jout(idx(p)) + dj)
        do d = 0, ndof - 1
          if (.not. lself(d)) rhs(p) = rhs(p) + xin(d) * c(d, idx(p))
        end do
        do q = 1, n
          A(p, q) = -c(idx(q), idx(p))
        end do
        A(p, p) = A(p, p) + one
      end do
      do r = 1, n ! gaussian elimination, the system is diagonally dominant
        do q = r + 1, n
          f = A(q, r) / A(r, r)
          A(q, r:n) = A(q, r:n) - f * A(r, r:n)
          rhs(q) = rhs(q) - f * rhs(r)
        end do
      end do
      do r = n, 1, -1
        do q = r + 1, n
          rhs(r) = rhs(r) - A(r, q) * rhs(q)
        end do
        rhs(r) = rhs(r) / A(r, r)
      end do
      do p = 1, n
        xin(idx(p)) = rhs(p)
      end do
    end if

    xout = matmul(xin, c)
    do d = 0, ndof - 1
      ! what leaves through the outer face on the high side has no place in the array and is not needed
      if (iout(d) .gt. solver%C_diff%gxe .or. jout(d) .gt. solver%C_diff%gye) cycle
      x0(d, kout(d), iout(d), jout(d)) = xb(d, kout(d), iout(d) + di, jout(d) + dj) + xout(d)
    end do
  end subroutine

  !> @brief returns the diffuse dof that is the same stream but the opposite direction
  pure function diff_inv_dof(solver, dof) result(inv_dof)
    class(t_solver), intent(in) :: solver
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

  subroutine explicit_ediff_sor_sweep(solver, coeffs, dx, dy, dz, omega, b, x)
    class(t_solver), intent(inout) :: solver
    real(ireals), target, intent(in) :: coeffs(:, :, :, :)
    integer(iintegers), dimension(3), intent(in) :: dx, dy, dz ! start, end, increment for each dimension
    real(ireals), intent(in) :: omega
    real(ireals), target, contiguous, intent(in) :: b(:, :, :, :)
    real(ireals), target, contiguous, intent(inout) :: x(:, :, :, :)

    real(ireals), pointer, dimension(:, :, :, :) :: x0, xb
    integer(iintegers) :: k, i, j
    integer(iintegers) :: idst, isrc, src, dst
    real(ireals), pointer :: v(:, :) ! dim(src, dst)
    integer(iintegers) :: msrc, mdst

    real(ireals) :: sigma

    x0 => null()
    xb => null()

    associate ( &
        & atm => solver%atm, &
        & C => solver%C_diff)

      x0(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => x
      xb(0:C%dof - 1, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye) => b

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
                x0(dst, mdst, i, j) = (one - omega) * x0(dst, mdst, i, j) + omega * (xb(dst, mdst, i, j) + sigma)
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
                x0(dst, k, mdst, j) = (one - omega) * x0(dst, k, mdst, j) + omega * (xb(dst, k, mdst, j) + sigma)
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
                x0(dst, k, i, mdst) = (one - omega) * x0(dst, k, i, mdst) + omega * (xb(dst, k, i, mdst) + sigma)
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

      nullify (x0, xb)
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

  !> Fill ghost cells of v from neighboring MPI ranks.
  !> v must have ghost extent 1: shape (dof, zm, gxm=xm+2, gym=ym+2).
  subroutine fill_ghost(comm, C, v, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_coord), intent(in) :: C
    real(ireals), intent(inout) :: v(:, :, :, :)
    integer(mpiint), intent(out) :: ierr

    ! v is (1:dof, 1:zm, 1:gxm, 1:gym) in 1-based assumed-shape indexing:
    !   interior x: 2..gxm-1, ghost west: 1, ghost east: gxm
    !   interior y: 2..gym-1, ghost south: 1, ghost north: gym

    real(ireals), allocatable :: se(:, :, :), sw(:, :, :), sn(:, :, :), ss(:, :, :)
    real(ireals), allocatable :: re(:, :, :), rw(:, :, :), rn(:, :, :), rs(:, :, :)
    integer(mpiint) :: rqs(8), sts(MPI_STATUS_SIZE, 8)
    integer(mpiint) :: nw, ne, ns, nn
    integer(iintegers) :: gxm, gym

    nw = int(C%neighbors(10), mpiint)
    ne = int(C%neighbors(16), mpiint)
    ns = int(C%neighbors(4), mpiint)
    nn = int(C%neighbors(22), mpiint)

    gxm = size(v, 3, kind=iintegers)
    gym = size(v, 4, kind=iintegers)

    allocate (se(size(v, 1), size(v, 2), size(v, 4) - 2)); se = zero
    allocate (sw(size(v, 1), size(v, 2), size(v, 4) - 2)); sw = zero
    allocate (re(size(v, 1), size(v, 2), size(v, 4) - 2)); re = zero
    allocate (rw(size(v, 1), size(v, 2), size(v, 4) - 2)); rw = zero
    allocate (sn(size(v, 1), size(v, 2), size(v, 3) - 2)); sn = zero
    allocate (ss(size(v, 1), size(v, 2), size(v, 3) - 2)); ss = zero
    allocate (rn(size(v, 1), size(v, 2), size(v, 3) - 2)); rn = zero
    allocate (rs(size(v, 1), size(v, 2), size(v, 3) - 2)); rs = zero

    se = v(:, :, gxm - 1, 2:gym - 1) ! east interior edge, interior y only
    sw = v(:, :, 2, 2:gym - 1)        ! west interior edge, interior y only
    sn = v(:, :, 2:gxm - 1, gym - 1) ! north interior edge
    ss = v(:, :, 2:gxm - 1, 2)       ! south interior edge

    ! Tag convention: sender uses its own neighbor index as tag so receiver
    ! can post Irecv with the same tag.
    call MPI_Irecv(rw, size(rw, kind=mpiint), imp_ireals, nw, 16, comm, rqs(1), ierr)
    call MPI_Irecv(re, size(re, kind=mpiint), imp_ireals, ne, 10, comm, rqs(2), ierr)
    call MPI_Irecv(rs, size(rs, kind=mpiint), imp_ireals, ns, 22, comm, rqs(3), ierr)
    call MPI_Irecv(rn, size(rn, kind=mpiint), imp_ireals, nn, 4, comm, rqs(4), ierr)
    call MPI_Isend(se, size(se, kind=mpiint), imp_ireals, ne, 16, comm, rqs(5), ierr)
    call MPI_Isend(sw, size(sw, kind=mpiint), imp_ireals, nw, 10, comm, rqs(6), ierr)
    call MPI_Isend(sn, size(sn, kind=mpiint), imp_ireals, nn, 22, comm, rqs(7), ierr)
    call MPI_Isend(ss, size(ss, kind=mpiint), imp_ireals, ns, 4, comm, rqs(8), ierr)
    call MPI_Waitall(8_mpiint, rqs, sts, ierr)

    v(:, :, 1, 2:gym - 1) = rw
    v(:, :, gxm, 2:gym - 1) = re
    v(:, :, 2:gxm - 1, 1) = rs
    v(:, :, 2:gxm - 1, gym) = rn
  end subroutine

end module
