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

module m_adaptive_spectral_integration
  use mpi, only: mpi_comm_rank
  use m_pprts_base, only: t_state_container
  use m_data_parameters, only: iintegers, ireals, default_str_len, mpiint, zero, one, nil
  use m_tenstream_options, only: options_max_solution_err, options_max_solution_time
  use m_helper_functions, only: approx, CHKERR, CHKWARN, toStr

  implicit none

  private
  public need_new_solution, polyfit

  logical, parameter :: ldebug = .false.

  real(ireals), parameter :: time_debug_solutions = zero  ! if enabled, we calculate new solutions but do not return the update solutions to the host model.(set to zero to disable)

contains

  function need_new_solution(comm, solution, time, lenable_solutions_err_estimates)
    integer(mpiint), intent(in) :: comm
    type(t_state_container) :: solution
    real(ireals), intent(in), optional :: time
    logical, intent(in) :: lenable_solutions_err_estimates
    logical :: need_new_solution

    integer, parameter :: Nfit = 3   ! Number of used residuals
    integer, parameter :: Nporder = 2 ! max order of polynomial
    real(ireals) :: t(Nfit), tm(Nfit), dt(Nfit - 1), err(2, 2 * (Nfit - 1)), error_estimate
    real(ireals) :: polyc(Nporder + 1), estimate(Nporder)

    character(default_str_len) :: reason
    integer, parameter :: out_unit = 20

    integer(iintegers) :: k, ipoly
    integer(mpiint) :: ierr, myid

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    ! Make time an optional argument here for
    ! convenience of the interface --
    ! otherwise the user needs to check if he
    ! opt_time is present and so on...
    if (.not. present(time)) then
      need_new_solution = .true.
      write (reason, *) 'time not provided'
      if (ldebug .and. myid .eq. 0) print *, 'new calc', need_new_solution, ' bc ', trim(reason), ' t', time, solution%uid
      return
    end if

    if ((.not. lenable_solutions_err_estimates) .or. &
        options_max_solution_time .le. zero .or. &
        options_max_solution_err .le. zero) then
      need_new_solution = .true.
      write (reason, *) 'solutions_err_estimates disabled'
      if (ldebug .and. myid .eq. 0) print *, 'new calc', need_new_solution, ' bc ', trim(reason), ' t', time, solution%uid
      return
    end if

    if (.not. solution%lset) then !if we did not store a solution, return immediately
      need_new_solution = .true.
      write (reason, *) 'no solution yet -- not initialized'
      if (ldebug .and. myid .eq. 0) print *, 'new calc', need_new_solution, ' bc ', trim(reason), ' t', time, solution%uid
      return
    end if

    do k = 1, Nfit
      t(k) = solution%time(Nfit - k + 1)
    end do

    ! t is pts where the solution got updated
    ! tm is in between those time
    ! dt is the weight for the integral
    do k = 1, Nfit - 1
      tm(k) = (t(k) + t(k + 1))*.5_ireals
    end do
    tm(Nfit) = (t(Nfit) + time)*.5_ireals

    do k = 1, Nfit - 1
      dt(k) = tm(k + 1) - tm(k)
    end do

    ! setup error timeseries on which to fit -- first index is time, second is error
    do k = 1, Nfit - 1
      err(1, 2 * k - 1) = tm(k)
      err(1, 2 * k) = t(k + 1)
      err(2, 2 * k - 1) = zero
      err(2, 2 * k) = solution%maxnorm(Nfit - k)
    end do

    ! try several polynomials and find max error:
    do ipoly = 1, Nporder
      polyc(1:ipoly + 1) = polyfit(err(1, :), err(2, :), ipoly, ierr) ! e.g. second order polynomial has 3 coefficients
      if (ierr .ne. 0) then
        need_new_solution = .true.
        write (reason, *) 'problem fitting error curve', ierr
        if (ldebug .and. myid .eq. 0) print *, 'new calc', need_new_solution, ' bc ', trim(reason), ' t', time, solution%uid
        return
      end if
      estimate(ipoly) = zero
      do k = 1, ipoly + 1
        estimate(ipoly) = estimate(ipoly) + polyc(k) * time**(k - 1)
      end do
    end do

    error_estimate = maxval(abs(estimate))

    if (myid .eq. 0 .and. ldebug) then ! .and. error_estimate.le.zero) then
      print *, 'DEBUG t', t
      print *, 'DEBUG tm', tm
      print *, 'DEBUG dt', dt
      print *, 'DEBUG err', err(1, :), '::', err(2, :)
      print *, 'DEBUG err_est', error_estimate, '::', estimate
      print *, 'DEBUG polyc', polyc
    end if

    if (error_estimate .le. options_max_solution_err) then
      need_new_solution = .false.
      write (reason, *) 'ERR_TOL_IN_BOUND'
    else
      need_new_solution = .true.
      write (reason, *) 'ERR_TOL_EXCEEDED'
    end if

    if (any(t .lt. zero)) then
      need_new_solution = .true.
      write (reason, *) 'FEW_SOLUTIONS'
    end if

    if (time - solution%time(1) .gt. options_max_solution_time) then
      need_new_solution = .true.
      write (reason, *) 'MIN_TIME_EXCEEDED'
    end if

    if (time_debug_solutions .gt. zero) then
      if (.not. need_new_solution) then
        need_new_solution = .true. ! overwrite it and calculate anyway
        write (reason, *) 'MANUAL OVERRIDE'
        ! Hack to monitor error growth...
        ! We tell the user that he has to calculate radiation again.
        ! We will calculate and update the solution vectors...
        ! But once he tries to load em, we give him the old, un-updated solution.
        ! This happens in restore_solution i.e. we overwrite solution with old one.
      end if
    end if

    !            if(ldebug .and. myid.eq.0 .and. .not. need_new_solution ) &
    if (ldebug .and. myid .eq. 0) then
      print *, ''
      print *, ''
      print *, 'new calc', need_new_solution, ' bc ', trim(reason), ' t', time, solution%uid, &
        '    ::     est.', error_estimate, '[W]', error_estimate * 86.1, '[K/d]'
      print *, ' dir residuals _solver ::', solution%dir_ksp_residual_history(1:4)
      print *, ' diff residuals _solver ::', solution%diff_ksp_residual_history(1:4)
      if (solution%uid .eq. 501) then
        open (unit=out_unit, file="residuals.log", action="write", status="replace")
      else
        open (unit=out_unit, file="residuals.log", action="readwrite", status="unknown", position="append")
      end if

      write (out_unit, *) solution%uid, solution%maxnorm
      write (out_unit, *) solution%uid, time
      close (out_unit)
    end if
  end function need_new_solution

  function polyfit(vx, vy, d, ierr) !Rosetta Code http://rosettacode.org/wiki/Polynomial_regression#Fortran
    implicit none
    integer(iintegers), intent(in) :: d
    real(ireals), dimension(d + 1) :: polyfit
    real(ireals), dimension(:), intent(in) :: vx, vy
    integer(mpiint), intent(out) :: ierr

    real(ireals), dimension(size(vx), d + 1) :: X
    real(ireals), dimension(d + 1, size(vx)) :: XT
    real(ireals), dimension(d + 1, d + 1) :: XTX

    integer(iintegers) :: i, j

    integer :: n, lda, lwork
    integer :: info
    integer, dimension(d + 1) :: ipiv
    real(ireals), dimension(d + 1) :: work

    ierr = 0

    n = int(d + 1, kind(n))
    lda = n
    lwork = n

    do i = 1, size(vx)
      if (any(approx(vx(i), vx(i + 1:size(vx))))) then ! polyfit cannot cope with same x values --> matrix gets singular
        polyfit = 0
        polyfit(1) = nil
        ierr = 1
        if(ldebug) &
          & call CHKWARN(ierr, 'polyfit cannot cope with same x values --> matrix gets singular'//toStr(vx)//' :: '//toStr(vy))
        return
      end if
    end do

    ! prepare the matrix
    do i = 0, d
      do j = 1, size(vx)
        X(j, i + 1) = vx(j)**i
      end do
    end do

    XT = transpose(X)
    XTX = matmul(XT, X)

    ! calls to LAPACK subs DGETRF and DGETRI
    if (sizeof(one) .eq. 4) then !single precision
      call SGETRF(n, n, XTX, lda, ipiv, info)
    else if (sizeof(one) .eq. 8) then !double precision
      call DGETRF(n, n, XTX, lda, ipiv, info)
    else
      call CHKERR(1_mpiint, 'dont know which lapack routine to call for reals with sizeof=='//toStr(sizeof(one)))
    end if
    if (info /= 0) then
      ierr = 2
      if (ldebug) call CHKWARN(ierr, "problem with lapack lsqr :: 1 :: info"//toStr(info))
      return
    end if
    if (sizeof(one) .eq. 4) then !single precision
      call SGETRI(n, XTX, lda, ipiv, work, lwork, info)
    else if (sizeof(one) .eq. 8) then !double precision
      call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    else
      call CHKERR(1_mpiint, 'dont know which lapack routine to call for reals with sizeof=='//toStr(sizeof(one)))
    end if
    if (info /= 0) then
      ierr = 3
      call CHKWARN(ierr, "problem with lapack lsqr :: 2 :: info"//toStr(info))
      return
    end if

    polyfit = matmul(matmul(XTX, XT), vy)
  end function
end module
