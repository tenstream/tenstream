module test_netcdfio
    use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters, default_str_len
    use iso_fortran_env, only: REAL32, REAL64

    use m_helper_functions, only: CHKERR, itoa

    use m_netcdfIO, only: ncwrite, ncload, acquire_file_lock, release_file_lock
    use m_c_syscall_wrappers, only: acquire_flock_lock, release_flock_lock

    use pfunit_mod

  implicit none

  character(len=*), parameter :: fname='pfunit_test.nc'

  contains

!@test(npes =[2]) ! not passing with a deadlock at the moment
subroutine test_c_lockf(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    if(myid.eq.0) then
      call acquire_flock_lock(fname, ierr)
      @assertEqual(0,ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call release_flock_lock(fname, ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
    else
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call acquire_flock_lock(fname, ierr)
      @assertEqual(0,ierr)
      call release_flock_lock(fname, ierr)
      @assertEqual(0,ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
    endif
    call release_flock_lock(fname, ierr)
    @assertEqual(0,ierr)
end subroutine

@test(npes =[2])
subroutine test_file_locks(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, flockunit

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    if(myid.eq.0) then
      call acquire_file_lock(fname, flockunit, ierr, waittime=1); call CHKERR(ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call release_file_lock(flockunit, ierr); call CHKERR(ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
    else
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call acquire_file_lock(fname, flockunit, ierr, waittime=1, blocking=.False.)
      @assertFalse(ierr.eq.0)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
      call acquire_file_lock(fname, flockunit, ierr, waittime=1, blocking=.False.)
      @assertTrue(ierr.eq.0)
      call release_file_lock(flockunit, ierr); call CHKERR(ierr)
    endif
    call release_file_lock(flockunit, ierr)
    @assertEqual(2,ierr)
end subroutine

@test(npes =[2])
subroutine test_netcdf_load_write_r32(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, rank
    real(REAL32),allocatable :: a1d(:)
    real(REAL32),allocatable :: a2d(:,:)

    integer(mpiint),parameter :: N=10
    character(len=default_str_len) :: groups(3)

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)//'_kind_'//itoa(kind(a1d))
    groups(3) = 'a1d'

    allocate(a1d(N), source=real(myid, kind(a1d)))
    do rank = 0, numnodes-1
      if(rank.eq.myid) then
        call ncwrite(groups, a1d, ierr); call CHKERR(ierr, 'Could not write 1d array to nc file')
      endif
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    enddo
    ! Now everyone has written his stuff into the netcdf file
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    ! Try reading it without barriers
    do rank = 0, numnodes-1
      deallocate(a1d)
      groups(2) = 'rank'//itoa(rank)
      call ncload(groups, a1d, ierr); call CHKERR(ierr, 'Could not read 1d array from nc file')
      @mpiassertEqual(real(rank, kind(a1d)), a1d(1))
      @mpiassertEqual(N, size(a1d))
    enddo

    allocate(a2d(N,N), source=real(myid, kind(a2d)))
    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)
    groups(3) = 'a2d'

    do rank = 0, numnodes-1
      if(rank.eq.myid) then
        call ncwrite(groups, a2d, ierr); call CHKERR(ierr, 'Could not write 2d array to nc file')
      endif
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    enddo
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    do rank= 0, numnodes-1
      deallocate(a2d)
      groups(2) = 'rank'//itoa(rank)
      call ncload(groups, a2d, ierr); call CHKERR(ierr, 'Could not read 2d array from nc file')
      @mpiassertEqual(real(rank, kind(a2d)), a2d(1,1))
      @mpiassertEqual(N, size(a2d, dim=1))
      @mpiassertEqual(N, size(a2d, dim=2))
    enddo
end subroutine


@test(npes =[2])
subroutine test_netcdf_load_write_r64(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, rank
    real(REAL64),allocatable :: a1d(:)
    real(REAL64),allocatable :: a2d(:,:)

    integer(mpiint),parameter :: N=10
    character(len=default_str_len) :: groups(3)

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)//'_kind_'//itoa(kind(a1d))
    groups(3) = 'a1d'

    allocate(a1d(N), source=real(myid, kind(a1d)))
    do rank = 0, numnodes-1
      if(rank.eq.myid) then
        call ncwrite(groups, a1d, ierr); call CHKERR(ierr, 'Could not write 1d array to nc file')
      endif
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    enddo
    ! Now everyone has written his stuff into the netcdf file
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    ! Try reading it without barriers
    do rank = 0, numnodes-1
      deallocate(a1d)
      groups(2) = 'rank'//itoa(rank)
      call ncload(groups, a1d, ierr); call CHKERR(ierr, 'Could not read 1d array from nc file')
      @mpiassertEqual(real(rank, kind(a1d)), a1d(1))
      @mpiassertEqual(N, size(a1d))
    enddo

    allocate(a2d(N,N), source=real(myid, kind(a2d)))
    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)
    groups(3) = 'a2d'

    do rank = 0, numnodes-1
      if(rank.eq.myid) then
        call ncwrite(groups, a2d, ierr); call CHKERR(ierr, 'Could not write 2d array to nc file')
      endif
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    enddo
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    do rank= 0, numnodes-1
      deallocate(a2d)
      groups(2) = 'rank'//itoa(rank)
      call ncload(groups, a2d, ierr); call CHKERR(ierr, 'Could not read 2d array from nc file')
      @mpiassertEqual(real(rank, kind(a2d)), a2d(1,1))
      @mpiassertEqual(N, size(a2d, dim=1))
      @mpiassertEqual(N, size(a2d, dim=2))
    enddo
end subroutine

end module
