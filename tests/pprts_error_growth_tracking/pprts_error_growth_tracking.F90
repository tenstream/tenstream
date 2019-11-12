module pprts_error_growth_tracking
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, one, zero

  use m_adaptive_spectral_integration, only: need_new_solution
  use m_pprts_base, only : t_coord, t_solver_3_10
  use m_pprts, only : init_pprts, set_angles, &
    set_optical_properties, solve_pprts, destroy_pprts, &
    pprts_get_result

  use m_tenstream_options, only: read_commandline_options

  use pfunit_mod

  implicit none

  type(t_solver_3_10) :: solver

contains
  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    ! Tidy up
    call destroy_pprts(solver, lfinalizepetsc=.True.)
  end subroutine teardown

  @test(npes =[2,1])
  subroutine error_growth_tracking(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: iter, k
    integer(mpiint) :: comm, myid, numnodes

    integer(iintegers),parameter :: nxp=3,nyp=3,nlyr=3
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals),parameter :: phi0=0, theta0=60
    real(ireals),parameter :: albedo=0, dz=dx
    real(ireals),parameter :: incSolar = 1000
    real(ireals) :: dz1d(nlyr)

    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    real(ireals) :: time
    logical :: lneed

    dz1d = dz

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)
    call init_pprts(comm, nlyr, nxp, nyp, dx, dy, phi0, theta0, solver, dz1d=dz1d)

    allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
    allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

    ! First solve for solar radiation

    kabs = 1._ireals/nlyr
    ksca = 1e-8
    g    = zero

    !allocate(fdn  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    !allocate(fup  (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    !allocate(fdiv (solver%C_one%zm,  solver%C_one%xm,  solver%C_one%ym))

    do iter=1,5
      do k=1,2
        time = real(iter, ireals)

        lneed = need_new_solution(comm, solver%solutions(k), time, solver%lenable_solutions_err_estimates)
        print *, myid, 'Need_new_solution?', k, time, ' :: ', lneed

        if(iter.le.3) then
          call assertTrue(lneed)
        else
          call assertFalse(lneed)
        endif

        call set_optical_properties(solver, albedo, kabs, ksca, g)
        call solve_pprts(solver, incSolar, opt_solution_uid=k, opt_solution_time=real(iter,ireals))

        !allocate(fdir (solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
        call pprts_get_result(solver, fdn, fup, fdiv, fdir, opt_solution_uid=k)
        deallocate(fdir, fdn, fup, fdiv)

      enddo
    enddo
  end subroutine

end module
