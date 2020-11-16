module pprts_error_growth_tracking
  use m_data_parameters, only : &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & one, zero

  use m_helper_functions, only : spherical_2_cartesian

  use m_adaptive_spectral_integration, only: need_new_solution

  use m_pprts_base, only : &
    & t_coord, &
    & t_solver, &
    & t_solver_3_10, &
    & t_solver_8_16, &
    & destroy_pprts

  use m_pprts, only : &
    & init_pprts, &
    & set_angles, &
    & set_optical_properties, &
    & solve_pprts, &
    & pprts_get_result

  use m_tenstream_options, only: read_commandline_options

  use pfunit_mod

  implicit none

contains
  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
    call read_commandline_options(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    call finalize_mpi(&
      & comm=this%getMpiCommunicator(), &
      & lfinalize_mpi=.False., &
      & lfinalize_petsc=.True.)
  end subroutine teardown

  @test(npes =[1,2])
  subroutine error_growth_tracking(this)
    class (MpiTestMethod), intent(inout) :: this

    type(t_solver_3_10) :: solver_3_10
    type(t_solver_8_16) :: solver_8_16

    call this_test(solver_3_10)
    call this_test(solver_8_16)

  contains
    subroutine this_test(solver)
      class(t_solver), intent(inout) :: solver
      integer(iintegers),parameter :: nxp=3,nyp=3,nlyr=3
      real(ireals),parameter :: dx=100,dy=dx
      real(ireals),parameter :: phi0=0, theta0=60
      real(ireals),parameter :: albedo=0.1, dz=dx
      real(ireals),parameter :: incSolar = 1000

      integer(iintegers) :: iter, k
      integer(mpiint) :: comm, myid

      real(ireals) :: dz1d(nlyr)

      real(ireals) :: sundir(3)
      real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
      real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

      real(ireals) :: time
      logical :: lneed

      dz1d = dz

      comm     = this%getMpiCommunicator()
      myid     = this%getProcessRank()

      sundir = spherical_2_cartesian(phi0, theta0)

      call init_pprts(comm, nlyr, nxp, nyp, dx, dy, sundir, solver, dz1d=dz1d)

      allocate(kabs(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
      allocate(ksca(solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))
      allocate(g   (solver%C_one%zm , solver%C_one%xm,  solver%C_one%ym ))

      ! First solve for solar radiation

      kabs = 1._ireals/nlyr
      ksca = 1e-8
      g    = zero

      do iter=1,5
        do k=1,2
          time = real(iter, ireals)

          lneed = need_new_solution(comm, solver%solutions(k), time, solver%lenable_solutions_err_estimates)
          print *, myid, 'Need_new_solution?', k, time, ' :: ', lneed

          if(iter.le.3) then
            @mpiassertTrue(lneed, 'in the first couple of iterations we expect to need a new solution')
          else
            @mpiassertFalse(lneed, 'after a couple of iterations and no error growth, we should not need a new one')
          endif

          call set_optical_properties(solver, albedo, kabs, ksca, g)
          call solve_pprts(solver, &
            & lthermal=.False., &
            & lsolar=.True., &
            & edirTOA=incSolar, &
            & opt_solution_uid=k, &
            & opt_solution_time=real(iter,ireals))

          call pprts_get_result(solver, fdn, fup, fdiv, fdir, opt_solution_uid=k)
          deallocate(fdir, fdn, fup, fdiv)
        enddo
      enddo
      call destroy_pprts(solver, lfinalizepetsc=.False.)
    end subroutine
  end subroutine
end module
