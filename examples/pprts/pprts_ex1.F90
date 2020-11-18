module m_examples_pprts_ex1
    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, pi
    use m_helper_functions, only : CHKERR, spherical_2_cartesian, get_arg, cstr, toStr
    use m_pprts, only : init_pprts, set_optical_properties, solve_pprts, pprts_get_result, set_angles
    use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline, destroy_pprts
    use m_tenstream_options, only: read_commandline_options

  implicit none

  contains
    subroutine pprts_ex1( &
      & comm, &
      & nxp, nyp, nv, &
      & dx, dy, &
      & phi0, theta0, &
      & albedo, dz, &
      & incSolar, &
      & dtau_clearsky, w0_clearsky, g_clearsky, &
      & cld_layer_idx, &
      & dtau_cloud, w0_cloud, g_cloud, &
      & fdir,fdn,fup,fdiv, &
      & lverbose )

      integer(mpiint), intent(in) :: comm

      integer(iintegers),intent(in) :: nxp,nyp,nv
      real(ireals), intent(in) :: dx,dy
      real(ireals), intent(in) :: phi0, theta0
      real(ireals), intent(in) :: albedo, dz
      real(ireals), intent(in) :: incSolar
      real(ireals), intent(in) :: dtau_clearsky, w0_clearsky, g_clearsky
      integer(iintegers), intent(in) :: cld_layer_idx(2)
      real(ireals), intent(in) :: dtau_cloud, w0_cloud, g_cloud
      real(ireals), allocatable,dimension(:,:,:), intent(out) :: fdir,fdn,fup,fdiv
      logical, intent(in), optional :: lverbose

      real(ireals), allocatable,dimension(:,:,:) :: kabs,ksca,g
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

      call init_pprts(comm, nv, nxp, nyp, dx,dy, sundir, solver, dz1d)

      allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
      allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
      allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

      kabs = dtau_clearsky/(dz*nv) * (1._ireals - w0_clearsky)
      ksca = dtau_clearsky/(dz*nv) * w0_clearsky
      g    = g_clearsky

      Ncld = 1 + cld_layer_idx(2) - cld_layer_idx(1)
      if(myid.eq.0 .and. lverbose) print *,'Have '//tostr(Ncld)//' cloud layer(s) between '//toStr(cld_layer_idx)
      do k = cld_layer_idx(1), cld_layer_idx(2)
        kabs(k, :, :) = kabs(k,:,:) + dtau_cloud / Ncld / dz * (1._ireals - w0_cloud)
        ksca(k, :, :) = ksca(k,:,:) + dtau_cloud / Ncld / dz * w0_cloud
        g   (k, :, :) = &
          & ( dtau_clearsky * w0_clearsky * g_clearsky &
          & + dtau_cloud * w0_cloud * g_cloud ) &
          & / (dtau_clearsky * w0_clearsky + dtau_cloud * w0_cloud)
      enddo

      call set_optical_properties(solver, albedo, kabs, ksca, g)
      call set_angles(solver, sundir)

      call solve_pprts(solver, &
        & lthermal=.False., &
        & lsolar=.True., &
        & edirTOA=incSolar )

      call pprts_get_result(solver, fdn,fup,fdiv,fdir)

      if(get_arg(.False., lverbose)) then
        call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
        do i = 0, numnodes-1
          if(myid.eq.i) then
            print *,''
            print *,cstr('rank '//toStr(myid), 'red')
            do k = 1, nv+1
              print *,k, cstr('fdir x=1', 'red'), fdir(k,1,:), 'y=1', fdir(k,:,1)
            enddo
            print *,''
            do k = 1, nv+1
              print *,k, cstr('fdn  x=1', 'blue'), fdn (k,1,:), 'y=1', fdn (k,:,1)
            enddo
            print *,''
            do k = 1, nv+1
              print *,k, cstr('fup  x=1', 'red'), fup (k,1,:), 'y=1', fup (k,:,1)
            enddo
            print *,''
            do k = 1, nv
              print *,k, cstr('fdiv x=1', 'blue'), fdiv(k,1,:), 'y=1', fdiv(k,:,1)
            enddo
          endif
          call mpi_barrier(comm, ierr); call CHKERR(ierr)
        enddo
      endif

      call destroy_pprts(solver, .False.)
    end subroutine
  end module
