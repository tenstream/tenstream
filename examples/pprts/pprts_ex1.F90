module m_examples_pprts_ex1
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, zero, pi
  use m_helper_functions, only: CHKERR, spherical_2_cartesian, get_arg, cstr, toStr
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, pprts_get_result, set_angles
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts
  use m_tenstream_options, only: read_commandline_options

  implicit none

contains
  subroutine pprts_ex1( &
    & comm, &
    & lthermal, &
    & lsolar, &
    & nxp, nyp, nv, &
    & dx, dy, &
    & phi0, theta0, &
    & albedo, dz, &
    & incSolar, &
    & Bplck, &
    & Bplck_srfc, &
    & dtau_clearsky, w0_clearsky, g_clearsky, &
    & cld_layer_idx, &
    & dtau_cloud, w0_cloud, g_cloud, &
    & fdir, fdn, fup, fdiv, &
    & lverbose)

    integer(mpiint), intent(in) :: comm

    logical, intent(in) :: lthermal, lsolar
    integer(iintegers), intent(in) :: nxp, nyp, nv
    real(ireals), intent(in) :: dx, dy
    real(ireals), intent(in) :: phi0, theta0
    real(ireals), intent(in) :: albedo, dz
    real(ireals), intent(in) :: incSolar, Bplck, Bplck_srfc
    real(ireals), intent(in) :: dtau_clearsky, w0_clearsky, g_clearsky
    integer(iintegers), intent(in) :: cld_layer_idx(2)
    real(ireals), intent(in) :: dtau_cloud, w0_cloud, g_cloud
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: fdir, fdn, fup, fdiv
    logical, intent(in), optional :: lverbose

    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, g, planck
    real(ireals), allocatable, dimension(:, :) :: planck_srfc
    real(ireals) :: dz1d(nv)
    real(ireals) :: sundir(3)

    class(t_solver), allocatable :: solver

    integer(mpiint) :: myid, numnodes, ierr
    integer(iintegers) :: i, k, Ncld

    call init_mpi_data_parameters(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    dz1d = dz

    sundir = spherical_2_cartesian(phi0, theta0)

    call init_pprts(comm, nv, nxp, nyp, dx, dy, sundir, solver, dz1d)

    associate ( &
        & Catm => solver%C_one_atm, &
        & Catm1 => solver%C_one_atm1)

      allocate (kabs(Catm%zm, Catm%xm, Catm%ym))
      allocate (ksca(Catm%zm, Catm%xm, Catm%ym))
      allocate (g(Catm%zm, Catm%xm, Catm%ym))

      if (lthermal) then
        allocate (planck(Catm1%zm, Catm1%xm, Catm1%ym), source=Bplck)
        allocate (planck_srfc(Catm1%xm, Catm1%ym), source=Bplck_srfc)
      end if
    end associate

    kabs = dtau_clearsky / dz / real(nv, ireals) * (1._ireals - w0_clearsky)
    ksca = dtau_clearsky / dz / real(nv, ireals) * w0_clearsky
    g = g_clearsky

    Ncld = 1 + cld_layer_idx(2) - cld_layer_idx(1)
    if (myid .eq. 0 .and. lverbose) print *, 'Have '//tostr(Ncld)//' cloud layer(s) between '//toStr(cld_layer_idx)
    do k = cld_layer_idx(1), cld_layer_idx(2)
      kabs(k, :, :) = kabs(k, :, :) + dtau_cloud / real(Ncld, ireals) / dz * (1._ireals - w0_cloud)
      ksca(k, :, :) = ksca(k, :, :) + dtau_cloud / real(Ncld, ireals) / dz * w0_cloud
      g(k, :, :) = &
        & (dtau_clearsky * w0_clearsky * g_clearsky &
        & + dtau_cloud * w0_cloud * g_cloud) &
        & / (dtau_clearsky * w0_clearsky + dtau_cloud * w0_cloud)
    end do

    if (lthermal) then
      call set_optical_properties(solver, albedo, kabs, ksca, g, &
        & planck=planck, planck_srfc=planck_srfc)
    else
      call set_optical_properties(solver, albedo, kabs, ksca, g)
    end if
    call set_angles(solver, sundir)

    call solve_pprts(solver, &
      & lthermal=lthermal, &
      & lsolar=lsolar, &
      & edirTOA=incSolar)

    call pprts_get_result(solver, fdn, fup, fdiv, fdir)

    if (get_arg(.false., lverbose)) then
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
      do i = 0, numnodes - 1
        if (myid .eq. i) then
          print *, ''
          print *, cstr('rank '//toStr(myid), 'red')
          do k = 1, nv + 1
            print *, k, cstr('fdir x=1', 'red'), fdir(k, 1, :), 'y=1', fdir(k, :, 1)
          end do
          print *, ''
          do k = 1, nv + 1
            print *, k, cstr('fdn  x=1', 'blue'), fdn(k, 1, :), 'y=1', fdn(k, :, 1)
          end do
          print *, ''
          do k = 1, nv + 1
            print *, k, cstr('fup  x=1', 'red'), fup(k, 1, :), 'y=1', fup(k, :, 1)
          end do
          print *, ''
          do k = 1, nv
            print *, k, cstr('fdiv x=1', 'blue'), fdiv(k, 1, :), 'y=1', fdiv(k, :, 1)
          end do
        end if
        call mpi_barrier(comm, ierr); call CHKERR(ierr)
      end do
    end if

    call destroy_pprts(solver, .false.)
  end subroutine
end module
