module test_pprts_distorted_mesh
  ! Do distorted meshes give reasonable results?
  ! A periodic hill, defined by a surface pressure deficit, with terrain following layers.
  ! Compare pprts against the MonteCarlo solver RayLi

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & zero, one, pi

  use pfunit_mod

  use m_helper_functions, only: &
    & CHKERR, &
    & deg2rad, &
    & insert_petsc_opt, &
    & spherical_2_cartesian, &
    & toStr

  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, pprts_get_result_toZero
  use m_pprts_external_solvers, only: destroy_rayli_info

  implicit none

  integer(iintegers), parameter :: Nx = 16, Ny = 16, Nz = 8
  real(ireals), parameter :: dx = 100, dz = 100
  real(ireals), parameter :: incSolar = 1000
  real(ireals), parameter :: tau_clearsky = .2_ireals
  real(ireals), parameter :: p_srfc = 1000, scale_height = 8000 ! [hPa], [m]
  real(ireals), parameter :: dP_hill = 30 ! [hPa] pressure deficit at the top of the hill, about 245m
  real(ireals), parameter :: max_slope = 30 ! [deg]

contains
  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    call destroy_rayli_info()
    call finalize_mpi(&
      & this%getMpiCommunicator(), &
      & lfinalize_mpi=.false., &
      & lfinalize_petsc=.true.)
  end subroutine teardown

  ! surface pressure deficit [hPa] of a periodic hill in the middle of the domain, i,j are global indices starting at 1
  pure function pressure_deficit(i, j) result(dP)
    integer(iintegers), intent(in) :: i, j
    real(ireals) :: dP, cx, cy
    cx = (one - cos(2 * pi * (real(i, ireals) - .5_ireals) / real(Nx, ireals))) / 2
    cy = (one - cos(2 * pi * (real(j, ireals) - .5_ireals) / real(Ny, ireals))) / 2
    dP = dP_hill * cx * cy
  end function

  ! layer thicknesses [m] of a column with terrain following pressure levels in a isothermal atmosphere
  pure function column_dz(i, j, lflat) result(dz1d)
    integer(iintegers), intent(in) :: i, j
    logical, intent(in) :: lflat
    real(ireals) :: dz1d(Nz)
    real(ireals) :: ptop, ps, plev(Nz + 1)
    integer(iintegers) :: k
    ptop = p_srfc * exp(-dz * real(Nz, ireals) / scale_height)
    ps = p_srfc
    if (.not. lflat) ps = ps - pressure_deficit(i, j)
    do k = 1, Nz + 1
      plev(k) = ptop + (ps - ptop) * real(k - 1, ireals) / real(Nz, ireals)
    end do
    do k = 1, Nz
      dz1d(k) = scale_height * log(plev(k + 1) / plev(k))
    end do
  end function

  pure function surface_height(i, j) result(h)
    integer(iintegers), intent(in) :: i, j
    real(ireals) :: h
    h = dz * real(Nz, ireals) - sum(column_dz(modulo(i - 1, Nx) + 1, modulo(j - 1, Ny) + 1, .false.))
  end function

  ! results have global shape and are only allocated on rank 0
  subroutine solve(comm, solvername, lgeometric, lflat, phi0, theta0, w0, albedo, edir, edn, eup, abso)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: solvername
    logical, intent(in) :: lgeometric, lflat
    real(ireals), intent(in) :: phi0, theta0, w0, albedo
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: edir, edn, eup, abso

    class(t_solver), allocatable :: solver
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g, dz3d
    real(ireals) :: dz1d(Nz)
    integer(iintegers) :: i, j, xs, xm, ys, ym
    integer(mpiint) :: ierr

    call insert_petsc_opt('-pprts_geometric_coeffs '//merge('yes', 'no ', lgeometric), ierr); call CHKERR(ierr)

    ! dz3d has to have the shape of the local domain, ask a solver on a regular mesh for the domain decomposition
    dz1d = dz
    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)
    call init_pprts(comm, Nz, Nx, Ny, dx, dx, spherical_2_cartesian(phi0, theta0), solver, dz1d=dz1d)
    xs = solver%C_one%xs; xm = solver%C_one%xm
    ys = solver%C_one%ys; ym = solver%C_one%ym
    call destroy_pprts(solver, lfinalizepetsc=.false.)
    deallocate (solver)

    allocate (dz3d(Nz, xm, ym))
    do j = 1, ym
      do i = 1, xm
        dz3d(:, i, j) = column_dz(xs + i, ys + j, lflat)
      end do
    end do

    call allocate_pprts_solver_from_commandline(solver, solvername, ierr); call CHKERR(ierr)
    call init_pprts(comm, Nz, Nx, Ny, dx, dx, spherical_2_cartesian(phi0, theta0), solver, dz3d=dz3d)

    allocate (kabs(Nz, xm, ym), source=(one - w0) * tau_clearsky / (dz * real(Nz, ireals)))
    allocate (ksca(Nz, xm, ym), source=w0 * tau_clearsky / (dz * real(Nz, ireals)))
    allocate (g(Nz, xm, ym), source=zero)

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call solve_pprts(solver, lthermal=.false., lsolar=.true., edirTOA=incSolar)
    call pprts_get_result_toZero(solver, edn, eup, abso, edir)
    call destroy_pprts(solver, lfinalizepetsc=.false.)
    call insert_petsc_opt('-pprts_geometric_coeffs no', ierr); call CHKERR(ierr)
  end subroutine

  subroutine stats(name, a, ref)
    character(len=*), intent(in) :: name
    real(ireals), intent(in) :: a(:, :), ref(:, :)
    print '(a,a30,6(a,f9.2))', '   ', name, ' mean', sum(a) / size(a), ' ref', sum(ref) / size(ref), &
      & ' bias', sum(a - ref) / size(a), ' rmse', sqrt(sum((a - ref)**2) / size(a)), &
      & ' max|d|', maxval(abs(a - ref)), ' ref range', maxval(ref) - minval(ref)
  end subroutine

  @test(npes=[1])
  subroutine test_distorted_mesh_vs_rayli(this)
    class(MpiTestMethod), intent(inout) :: this

    real(ireals), parameter :: phis(3) = [real(ireals) :: 270, 0, 225]
    real(ireals), parameter :: thetas(2) = [real(ireals) :: 40, 60]
    real(ireals), parameter :: w0s(2) = [real(ireals) :: 0, .9]
    integer(iintegers), parameter :: Ncases = size(phis) * size(thetas) * size(w0s)

    type t_r
      real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso
    end type
    type(t_r) :: ray(Ncases), geo(Ncases), lut(Ncases), flat(Ncases)
    real(ireals) :: albedo, steepest, h(0:Nx + 1, 0:Ny + 1), aratio(Nx, Ny)
    real(ireals), allocatable :: ca(:, :), cb(:, :)
    integer(iintegers) :: icase, iphi, itheta, iw, i, j, k
    integer(mpiint) :: comm, myid

    comm = this%getMpiCommunicator()
    myid = this%getProcessRank()

    icase = 0
    do iw = 1, size(w0s)
      albedo = merge(zero, .2_ireals, iw .eq. 1)
      do itheta = 1, size(thetas)
        do iphi = 1, size(phis)
          icase = icase + 1
          associate (p => phis(iphi), t => thetas(itheta), w => w0s(iw))
            call solve(comm, 'rayli', .false., .false., p, t, w, albedo, ray(icase)%edir, ray(icase)%edn, ray(icase)%eup, &
                       ray(icase)%abso)
            call solve(comm, '3_10', .true., .false., p, t, w, albedo, geo(icase)%edir, geo(icase)%edn, geo(icase)%eup, &
                       geo(icase)%abso)
            call solve(comm, '3_10', .false., .false., p, t, w, albedo, lut(icase)%edir, lut(icase)%edn, lut(icase)%eup, &
                       lut(icase)%abso)
            call solve(comm, '3_10', .false., .true., p, t, w, albedo, flat(icase)%edir, flat(icase)%edn, flat(icase)%eup, &
                       flat(icase)%abso)
          end associate
        end do
      end do
    end do

    if (myid .ne. 0) return

    do j = 0, Ny + 1
      do i = 0, Nx + 1
        h(i, j) = surface_height(i, j)
      end do
    end do
    print *, 'surface height [m]'
    do j = Ny, 1, -1
      print '(*(i4))', (nint(h(i, j)), i=1, Nx)
    end do
    steepest = max(maxval(abs(h(1:Nx + 1, :) - h(0:Nx, :))), maxval(abs(h(:, 1:Ny + 1) - h(:, 0:Ny)))) / dx
    print *, 'steepest slope [deg]', atan(steepest) * 180 / pi
    ! ratio of horizontal to actual surface area, from centered height differences
    do j = 1, Ny
      do i = 1, Nx
        aratio(i, j) = one / sqrt(one + ((h(i + 1, j) - h(i - 1, j)) / (2 * dx))**2 + ((h(i, j + 1) - h(i, j - 1)) / (2 * dx))**2)
      end do
    end do
    print *, 'mean ratio of horizontal to surface area', sum(aratio) / size(aratio)
    @assertTrue(steepest .le. tan(deg2rad(max_slope)), 'slopes must not exceed 30 degrees')

    icase = 0
    do iw = 1, size(w0s)
      do itheta = 1, size(thetas)
        do iphi = 1, size(phis)
          icase = icase + 1
          print *, ''
          print *, '=== phi0', phis(iphi), 'theta0', thetas(itheta), 'w0', w0s(iw)
          k = Nz + 1
          print *, 'surface edir rayli | geometric | lut'
          do j = Ny, 1, -1
            print '(*(i4))', (nint(ray(icase)%edir(k, i, j)), i=1, Nx), -1, (nint(geo(icase)%edir(k, i, j)), i=1, Nx), -1, &
              & (nint(lut(icase)%edir(k, i, j)), i=1, Nx)
          end do
          print *, 'energy budget [W/m2]: incoming', incSolar * cos(deg2rad(thetas(itheta))), &
    & 'rayli edir+edn-eup at srfc', sum(ray(icase)%edir(k, :, :) + ray(icase)%edn(k, :, :) - ray(icase)%eup(k, :, :)) / (Nx * Ny), &
            & 'rayli toa eup', sum(ray(icase)%eup(1, :, :)) / (Nx * Ny), &
    & 'pprts edir+edn-eup at srfc', sum(geo(icase)%edir(k, :, :) + geo(icase)%edn(k, :, :) - geo(icase)%eup(k, :, :)) / (Nx * Ny), &
            & 'pprts toa eup', sum(geo(icase)%eup(1, :, :)) / (Nx * Ny)
          call stats('srfc edir geometric', geo(icase)%edir(k, :, :), ray(icase)%edir(k, :, :))
          call stats('srfc edir geo vs rayli/nz', geo(icase)%edir(k, :, :), ray(icase)%edir(k, :, :) / aratio)
          call stats('srfc edir lut vs rayli/nz', lut(icase)%edir(k, :, :), ray(icase)%edir(k, :, :) / aratio)
          call stats('srfc edir flat vs rayli/nz', flat(icase)%edir(k, :, :), ray(icase)%edir(k, :, :) / aratio)
          call stats('mid level edir geo', geo(icase)%edir(Nz / 2 + 1, :, :), ray(icase)%edir(Nz / 2 + 1, :, :))
          call stats('toa edir geo', geo(icase)%edir(1, :, :), ray(icase)%edir(1, :, :))
          call stats('srfc edir lut', lut(icase)%edir(k, :, :), ray(icase)%edir(k, :, :))
          call stats('srfc edir flat mesh', flat(icase)%edir(k, :, :), ray(icase)%edir(k, :, :))
          call stats('srfc edn geometric', geo(icase)%edn(k, :, :), ray(icase)%edn(k, :, :))
          call stats('srfc edn lut', lut(icase)%edn(k, :, :), ray(icase)%edn(k, :, :))
          call stats('toa eup geometric', geo(icase)%eup(1, :, :), ray(icase)%eup(1, :, :))
          call stats('toa eup lut', lut(icase)%eup(1, :, :), ray(icase)%eup(1, :, :))
          ! column absorption [W/m2] needs the layer thickness, use flux divergence instead: abso*dz summed
          allocate (ca(Nx, Ny), cb(Nx, Ny))
          do j = 1, Ny
            do i = 1, Nx
              ca(i, j) = sum(geo(icase)%abso(:, i, j) * column_dz(i, j, .false.))
              cb(i, j) = sum(ray(icase)%abso(:, i, j) * column_dz(i, j, .false.))
            end do
          end do
          call stats('column abso geometric', ca, cb)
          do j = 1, Ny
            do i = 1, Nx
              ca(i, j) = sum(lut(icase)%abso(:, i, j) * column_dz(i, j, .false.))
            end do
          end do
          call stats('column abso lut', ca, cb)
          deallocate (ca, cb)
        end do
      end do
    end do
  end subroutine
end module
