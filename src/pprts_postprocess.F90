!-------------------------------------------------------------------------
! This file is part of the TenStream solver.
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

module m_pprts_postprocess
  use mpi, only: mpi_comm_rank
  use m_data_parameters, only: iintegers, ireals, zero, one, i0, i1, i2, i3, i4, i5, i6, mpiint
  use m_pprts_base, only: t_solver, convolve_ediff_srfc, halo_fill_5pt
  use m_helper_functions, only: &
    & CHKERR, &
    & cross_3d, &
    & get_petsc_opt, &
    & imp_allreduce_max, &
    & imp_allreduce_mean, &
    & imp_allreduce_min, &
    & toStr
  use m_tenstream_log, only: t_ts_log_event, ts_log_begin, ts_log_end, ts_log_event_register

  implicit none

  private
  public :: smooth_surface_fluxes, slope_correction_fluxes, ts_nca, nca_wrapper

  logical, parameter :: ldebug = .false.
  type(t_ts_log_event), save :: log_smooth = t_ts_log_event(ts_id=-1)

contains

  subroutine smooth_surface_fluxes(solver, edn, eup)
    class(t_solver), intent(inout) :: solver
    real(ireals), allocatable, dimension(:, :, :), intent(inout) :: edn, eup
    integer(iintegers) :: i, kernel_width, Niter
    real(ireals) :: radius, mflx_up, mflx_dn
    logical :: lflg
    integer(mpiint) :: myid, ierr

    if (log_smooth%ts_id < 0) then
      call ts_log_event_register('pprts_smooth_surface_fluxes', log_smooth, ierr)
      call CHKERR(ierr)
    end if
    call ts_log_begin(log_smooth, ierr); call CHKERR(ierr)
    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    radius = 0
    call get_petsc_opt(solver%prefix, "-pprts_smooth_srfc_flx", radius, lflg, ierr); call CHKERR(ierr)

    if (lflg) then

      if (radius .lt. zero) then
        call imp_allreduce_mean(solver%comm, edn(ubound(edn, 1), :, :), mflx_dn); edn(ubound(edn, 1), :, :) = mflx_dn
        call imp_allreduce_mean(solver%comm, eup(ubound(eup, 1), :, :), mflx_up); eup(ubound(eup, 1), :, :) = mflx_up
        if (ldebug .and. myid .eq. 0) &
          print *, 'Smoothing diffuse srfc fluxes over the entire domain '// &
          'mean downward flx', mflx_dn, 'up', mflx_up

      else

        call find_iter_and_kernelwidth(solver, radius, Niter, kernel_width)
        do i = 1, Niter
          call convolve_ediff_srfc(solver%comm, solver%C_diff, kernel_width, edn(ubound(edn, 1):ubound(edn, 1), :, :))
          call convolve_ediff_srfc(solver%comm, solver%C_diff, kernel_width, eup(ubound(eup, 1):ubound(eup, 1), :, :))
        end do
      end if
    end if
    call ts_log_end(log_smooth, ierr); call CHKERR(ierr)
  end subroutine

  subroutine find_iter_and_kernelwidth(solver, radius, Niter, kernel_width)
    class(t_solver), intent(in) :: solver
    real(ireals), intent(in) :: radius
    integer(iintegers), intent(out) :: kernel_width, Niter
    integer(iintegers), parameter :: test_Ni = 500
    integer(iintegers), parameter :: Ni_limit = 250
    integer(iintegers) :: i, min_iter
    real(ireals) :: test_k(test_Ni), residuals(test_Ni)
    real(ireals) :: radius_in_pixel
    integer(mpiint) :: myid, ierr

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    radius_in_pixel = radius / ((solver%atm%dx + solver%atm%dy) / 2)
    ! Try a couple of number of iterations and determine optimal kernel_width
    min_iter = 1
    do i = 1, test_Ni
      test_k(i) = (sqrt((12 * radius_in_pixel**2 + real(i, ireals)) / real(i, ireals) + 1) - 1) / 2
      if (nint(test_k(i)) .ge. min(solver%C_diff%xm, solver%C_diff%ym)) min_iter = i + 1
    end do
    if (min_iter .gt. Ni_limit) &
      call CHKERR(int(min_iter, mpiint), 'the smoothing iteration count would be larger than '// &
                  toStr(Ni_limit)//'... this could become expensive, you can set the value higher but I'// &
                  'suspect that you are trying something weird')
    ! We want it
    ! as close as possible to an integer
    ! and as big as possible
    ! but not bigger than the local domain size
    ! or a certain max value
    residuals = abs(test_k - nint(test_k))

    Niter = min_iter - 1 + minloc(residuals(min_iter:test_Ni), dim=1)
    kernel_width = nint(test_k(Niter))

    call imp_allreduce_max(solver%comm, Niter, i); Niter = i
    call imp_allreduce_min(solver%comm, kernel_width, i); kernel_width = i

    if (kernel_width .eq. 0) Niter = 0

    if (ldebug .and. myid .eq. 0) then
      do i = 1, test_Ni
        print *, 'iter', i, 'test_k', test_k(i), residuals(i)
      end do
      print *, 'Smoothing diffuse srfc fluxes with radius', solver%atm%dx, radius, radius_in_pixel, &
        'Niter', Niter, 'kwidth', kernel_width
    end if

  end subroutine

  subroutine slope_correction_fluxes(solver, edir)
    class(t_solver), intent(inout) :: solver
    real(ireals), allocatable, dimension(:, :, :), intent(inout) :: edir
    logical :: lslope_correction, latm_correction, lflg
    integer(mpiint) :: myid, ierr

    real(ireals), pointer :: grad(:, :, :, :)
    real(ireals) :: fac, n(3)
    integer(iintegers) :: i, j, k

    grad => null()

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    lslope_correction = .false.
    call get_petsc_opt(solver%prefix, "-pprts_slope_correction", lslope_correction, lflg, ierr); call CHKERR(ierr)

    latm_correction = .false.
    call get_petsc_opt(solver%prefix, "-pprts_atm_correction", latm_correction, lflg, ierr); call CHKERR(ierr)

    if (.not. lslope_correction .and. .not. latm_correction) return

    associate ( &
      atm => solver%atm, &
      sun => solver%sun, &
      C_dir => solver%C_dir, &
      C_two1 => solver%C_two1)

      grad => atm%hgrad

      if (latm_correction) then
        do j = C_two1%ys, C_two1%ye
          do i = C_two1%xs, C_two1%xe
            do k = C_two1%zs, C_two1%ze

              n = cross_3d([one, zero, grad(i0, k, i, j)], [zero, one, grad(i1, k, i, j)])
              n = n / norm2(n)

              fac = dot_product(solver%sun%sundir, [zero, zero, one]) / dot_product(solver%sun%sundir, n)
              edir(k - C_two1%zs + 1, i - C_two1%xs + 1, j - C_two1%ys + 1) = &
                & edir(k - C_two1%zs + 1, i - C_two1%xs + 1, j - C_two1%ys + 1) * fac
            end do
          end do
        end do
      end if

      if (lslope_correction) then
        k = C_two1%ze
        do j = C_two1%ys, C_two1%ye
          do i = C_two1%xs, C_two1%xe

            n = cross_3d([one, zero, grad(i0, k, i, j)], [zero, one, grad(i1, k, i, j)])
            n = n / norm2(n)

            fac = dot_product(solver%sun%sundir, n) / dot_product(solver%sun%sundir, [zero, zero, one])
            fac = max(0._ireals, fac)
            edir(k - C_two1%zs + 1, i - C_two1%xs + 1, j - C_two1%ys + 1) = &
              & edir(k - C_two1%zs + 1, i - C_two1%xs + 1, j - C_two1%ys + 1) * fac
          end do
        end do
      end if

      nullify (grad)
    end associate
  end subroutine

  !> @brief Pack fluxes + atmosphere into a ghost-extended buffer and call ts_nca.
  !> @details ediff and abso are local arrays with no ghost cells (interior only).
  !> The wrapper allocates a 1-cell ghost halo, fills it with halo_fill_5pt,
  !> calls ts_nca, and writes heating-rate / dz back into abso.
  subroutine nca_wrapper(solver, ediff, abso)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: ediff(0:, :, :, :)    ! (0:dof-1, zs:ze, xs:xe, ys:ye) interior only
    real(ireals), intent(inout) :: abso(0:, :, :, :)  ! (0:0, zs:ze, xs:xe, ys:ye) interior only

    integer(iintegers), parameter :: E_up = i0, E_dn = i1, idz = i2, iplanck = i3, ikabs = i4, ihr = i5
    integer(iintegers) :: k
    integer(mpiint) :: ierr
    real(ireals), allocatable :: lnca(:, :, :, :)

    associate ( &
      atm => solver%atm, &
      C_diff => solver%C_diff, &
      C_one => solver%C_one)

      if (C_diff%dof < i6) &
        call CHKERR(1_mpiint, 'For NCA, need a solver with at least 6 diffuse streams to copy over some data')

      ! ghost-extended buffer: (0:5, zm, gxm=xm+2, gym=ym+2)
      allocate (lnca(0:5, C_diff%zs:C_diff%ze, &
                     C_diff%xs - 1:C_diff%xe + 1, &
                     C_diff%ys - 1:C_diff%ye + 1), source=zero)

      lnca(idz, C_diff%zs:C_diff%ze - 1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = atm%dz
      lnca(iplanck, C_diff%zs:C_diff%ze, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = atm%planck
      lnca(ikabs, C_diff%zs:C_diff%ze - 1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = atm%kabs
      lnca(E_up, :, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = ediff(E_up, :, :, :)
      lnca(E_dn, :, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = ediff(E_dn, :, :, :)

      call halo_fill_5pt(solver%comm, C_diff, lnca, ierr); call CHKERR(ierr)

      call ts_nca(atm%dx, atm%dy, &
                  lnca(idz, :, :, :), &
                  lnca(iplanck, :, :, :), &
                  lnca(ikabs, :, :, :), &
                  lnca(E_dn, :, :, :), &
                  lnca(E_up, :, :, :), &
                  lnca(ihr, :, :, :))

      do k = C_one%zs, C_one%ze
        abso(i0, k, :, :) = lnca(ihr, k, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) / &
                             lnca(idz, k, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye)
      end do
    end associate
  end subroutine nca_wrapper

  ! Neighbouring Column Approximation
  ! RTE Solver for the thermal spectral range, calculation of heating rates.
  ! Klinger and Mayer (2015), doi:10.1016/j.jqsrt.2015.04.002
  ! Contact Carolin Klinger for an implementation.
  subroutine ts_nca(dx, dy, dz, B, kabs_3d, Edn, Eup, hr)
    real(ireals), intent(in) :: dx, dy
    real(ireals), intent(in), dimension(:, :, :) :: dz, kabs_3d  ! (Nlevel,Nx,Ny) including ghosts
    real(ireals), intent(in), dimension(:, :, :) :: B, Edn, Eup  ! (Nlevel,Nx,Ny) including ghosts
    real(ireals), intent(out), dimension(:, :, :) :: hr           ! (Nlevel,Nx,Ny) including ghosts

    hr = dz + kabs_3d + B + Edn + Eup * (dx + dy) ! suppress unused-variable warnings
    hr = -9999999._ireals
    stop 'NCA not freely available.... please consult Carolin Klinger for an implementation.'
  end subroutine ts_nca

end module m_pprts_postprocess
