module test_ranlux
  use m_data_parameters, only: ireals, mpiint, iintegers, init_mpi_data_parameters
  use m_ranlux, only: ranlux, rluxgo

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: comm
    comm = this%getMpiCommunicator()
    call init_mpi_data_parameters(comm)
  end subroutine setup

  @after
  subroutine teardown(this)
    class(MpiTestMethod), intent(inout) :: this
    logical :: lpetsc_is_initialized
    integer(mpiint) :: ierr
    call PetscInitialized(lpetsc_is_initialized, ierr)
    if (lpetsc_is_initialized) call PetscFinalize(ierr)
  end subroutine teardown

  @test(npes=[1])
  subroutine test_ranlux_call(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(iintegers), parameter :: N = 100000, iter = 3
    integer :: ilvl, i
    real, allocatable :: R(:)
    double precision :: s, e

    allocate(R(N))

    do ilvl = 1, 4
      call RLUXGO(ilvl, 1, 0, 0)
      call cpu_time(s)
      do i = 1, iter
        call RANLUX(R, int(N))
      end do
      call cpu_time(e)
      @assertGreaterThan(1, R)
      @assertLessThan(0, R)
      print *, 'Ranlux level', ilvl, 'time', e - s
    end do
  end subroutine
end module
