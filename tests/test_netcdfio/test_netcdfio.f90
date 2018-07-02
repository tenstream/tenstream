module test_netcdfio
    use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters

    use m_helper_functions, only: CHKERR, itoa

    use m_netcdfIO, only: ncwrite, ncload, acquire_file_lock, release_file_lock

    use pfunit_mod

  implicit none

  character(len=*), parameter :: fname='pfunit_test.nc'

  contains

@test(npes =[2])
subroutine test_file_locks(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, flockunit

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    if(myid.eq.0) then
      call acquire_file_lock(fname, flockunit, ierr); call CHKERR(ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call release_file_lock(flockunit, ierr); call CHKERR(ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
    else
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call acquire_file_lock(fname, flockunit, ierr, blocking=.False.)
      @assertFalse(ierr.eq.0)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
      call acquire_file_lock(fname, flockunit, ierr, blocking=.False.)
      @assertTrue(ierr.eq.0)
      call release_file_lock(flockunit, ierr); call CHKERR(ierr)
    endif
    call release_file_lock(flockunit, ierr)
    @assertEqual(2,ierr)
end subroutine

@test(npes =[2])
subroutine test_netcdf_load_write(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, rank
    real(ireals),allocatable :: a1d(:)
    real(ireals),allocatable :: a2d(:,:)

    integer(mpiint),parameter :: N=10

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    allocate(a1d(N), source=real(myid, ireals))
    call ncwrite([fname, 'rank'//itoa(myid), 'a1d'], a1d, ierr); call CHKERR(ierr, 'Could not write 1d array to nc file')
    call mpi_barrier(comm, ierr); call CHKERR(ierr)
    do rank= 0, numnodes-1
      deallocate(a1d)
      call ncload([fname, 'rank'//itoa(rank), 'a1d'], a1d, ierr); call CHKERR(ierr, 'Could not read 1d array from nc file')
      @assertEqual(real(rank, ireals), a1d(1))
      @assertEqual(N, size(a1d))
    enddo

    allocate(a2d(N,N), source=real(myid, ireals))
    call ncwrite([fname, 'rank'//itoa(myid), 'a2d'], a2d, ierr); call CHKERR(ierr, 'Could not write 2d array to nc file')
    call mpi_barrier(comm, ierr); call CHKERR(ierr)
    do rank= 0, numnodes-1
      deallocate(a2d)
      call ncload([fname, 'rank'//itoa(rank), 'a2d'], a2d, ierr); call CHKERR(ierr, 'Could not read 2d array from nc file')
      @assertEqual(real(rank, ireals), a2d(1,1))
      @assertEqual(N, size(a2d, dim=1))
      @assertEqual(N, size(a2d, dim=2))
    enddo
end subroutine

end module
