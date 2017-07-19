@test(npes =[2,1]) 
subroutine error_growth_tracking(this)

    use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, one, zero, myid

    use m_tenstream, only : init_tenstream, set_angles, need_new_solution, &
        set_optical_properties, solve_tenstream, destroy_tenstream, &
        tenstream_get_result, getvecpointer, restorevecpointer, &
        t_coord,C_diff,C_one

    use m_tenstream_options, only: read_commandline_options

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: iter, k, numnodes, comm
    !integer(iintegers) :: myid

    integer(iintegers),parameter :: nxp=3,nyp=3,nv=3
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: albedo=0, dz=dx
    real(ireals),parameter :: incSolar = 1000
    real(ireals) :: dz1d(nv)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    real(ireals) :: time
    logical :: lneed

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()

    call init_mpi_data_parameters(comm)
    call init_tenstream(comm, nv, nxp, nyp, dx, dy, phi0, theta0, dz1d=dz1d)

    allocate(kabs(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(ksca(C_one%zm , C_one%xm,  C_one%ym ))
    allocate(g   (C_one%zm , C_one%xm,  C_one%ym ))

    ! First solve for solar radiation

    kabs = 1._ireals/nv
    ksca = 1e-8
    g    = zero

    allocate(fdn  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fup  (C_diff%zm, C_diff%xm, C_diff%ym))
    allocate(fdiv (C_one%zm, C_one%xm, C_one%ym))

    do iter=1,5
      do k=1,2
        time = iter*one

        lneed = need_new_solution(k, time)
        print *, myid, 'Need_new_solution?', k, time, ' :: ', lneed

        if(iter.le.3) then
            call assertTrue(lneed)
        else
            call assertFalse(lneed)
        endif

        call set_optical_properties(albedo, kabs, ksca, g)
        call solve_tenstream(incSolar, opt_solution_uid=k, opt_solution_time=iter*one)

        allocate(fdir (C_diff%zm, C_diff%xm, C_diff%ym))
        call tenstream_get_result(fdir, fdn, fup, fdiv, opt_solution_uid=k)
        deallocate(fdir)

      enddo
    enddo

    call destroy_tenstream(.True.)
end subroutine
