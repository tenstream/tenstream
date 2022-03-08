module test_netcdfio
  use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters, default_str_len
  use iso_fortran_env, only: int32, int64, real32, real64

  use m_helper_functions, only: CHKERR, itoa, char_arr_to_str

  use m_netcdfIO, only: ncwrite, ncload, acquire_file_lock, release_file_lock, &
                        get_global_attribute, set_global_attribute
  use m_c_syscall_wrappers, only: acquire_flock_lock, release_flock_lock

  use pfunit_mod

  implicit none

  character(len=*), parameter :: fname = 'pfunit_test.nc'

contains

!@test(npes =[2]) ! not passing with a deadlock at the moment
  subroutine test_c_lockf(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    if (myid .eq. 0) then
      call acquire_flock_lock(fname, ierr)
      @assertEqual(0, ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call release_flock_lock(fname, ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
    else
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call acquire_flock_lock(fname, ierr)
      @assertEqual(0, ierr)
      call release_flock_lock(fname, ierr)
      @assertEqual(0, ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
    end if
    call release_flock_lock(fname, ierr)
    @assertEqual(0, ierr)
  end subroutine

  @test(npes=[2])
  subroutine test_file_locks(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, flockunit

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    if (myid .eq. 0) then
      call acquire_file_lock(fname, flockunit, ierr, waittime=1); call CHKERR(ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call release_file_lock(flockunit, ierr); call CHKERR(ierr)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
    else
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 1
      call acquire_file_lock(fname, flockunit, ierr, waittime=1, blocking=.false.)
      @assertFalse(ierr .eq. 0)
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 2
      call mpi_barrier(comm, ierr); call CHKERR(ierr) ! 3
      call acquire_file_lock(fname, flockunit, ierr, waittime=1, blocking=.false.)
      @assertTrue(ierr .eq. 0)
      call release_file_lock(flockunit, ierr); call CHKERR(ierr)
    end if
    call release_file_lock(flockunit, ierr)
    @assertEqual(2, ierr)
  end subroutine

  @test(npes=[2])
  subroutine test_netcdf_load_write_r32(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, rank
    real(real32), allocatable :: a1d(:)
    real(real32), allocatable :: a2d(:, :)

    integer(mpiint), parameter :: N = 10
    character(len=default_str_len) :: groups(3)

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)//'_kind_'//itoa(kind(a1d))
    groups(3) = 'a1d'

    allocate (a1d(N), source=real(myid, kind(a1d)))
    do rank = 0, numnodes - 1
      if (rank .eq. myid) then
        call ncwrite(groups, a1d, ierr); call CHKERR(ierr, 'Could not write 1d array to nc file')
      end if
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    end do
    ! Now everyone has written his stuff into the netcdf file
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    ! Try reading it without barriers
    do rank = 0, numnodes - 1
      deallocate (a1d)
      groups(2) = 'rank'//itoa(rank)//'_kind_'//itoa(kind(a1d))
      call ncload(groups, a1d, ierr); call CHKERR(ierr, 'Could not read 1d array from nc file: '//char_arr_to_str(groups, ' / '))
      @mpiassertEqual(real(rank, kind(a1d)), a1d(1))
      @mpiassertEqual(N, size(a1d))
    end do

    allocate (a2d(N, N), source=real(myid, kind(a2d)))
    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)//'_kind_'//itoa(kind(a2d))
    groups(3) = 'a2d'

    do rank = 0, numnodes - 1
      if (rank .eq. myid) then
        call ncwrite(groups, a2d, ierr); call CHKERR(ierr, 'Could not write 2d array to nc file')
      end if
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    end do
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    do rank = 0, numnodes - 1
      deallocate (a2d)
      groups(2) = 'rank'//itoa(rank)//'_kind_'//itoa(kind(a2d))
      call ncload(groups, a2d, ierr); call CHKERR(ierr, 'Could not read 2d array from nc file: '//char_arr_to_str(groups, ' / '))
      @mpiassertEqual(real(rank, kind(a2d)), a2d(1, 1))
      @mpiassertEqual(N, size(a2d, dim=1))
      @mpiassertEqual(N, size(a2d, dim=2))
    end do
  end subroutine

  @test(npes=[2])
  subroutine test_netcdf_load_write_r64(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr, rank
    real(real64), allocatable :: a1d(:)
    real(real64), allocatable :: a2d(:, :)

    integer(mpiint), parameter :: N = 10
    character(len=default_str_len) :: groups(3)

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)//'_kind_'//itoa(kind(a1d))
    groups(3) = 'a1d'

    allocate (a1d(N), source=real(myid, kind(a1d)))
    do rank = 0, numnodes - 1
      if (rank .eq. myid) then
        call ncwrite(groups, a1d, ierr); call CHKERR(ierr, 'Could not write 1d array to nc file')
      end if
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    end do
    ! Now everyone has written his stuff into the netcdf file
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    ! Try reading it without barriers
    do rank = 0, numnodes - 1
      deallocate (a1d)
      groups(2) = 'rank'//itoa(rank)//'_kind_'//itoa(kind(a1d))
      call ncload(groups, a1d, ierr); call CHKERR(ierr, 'Could not read 1d array from nc file: '//char_arr_to_str(groups, ' / '))
      @mpiassertEqual(real(rank, kind(a1d)), a1d(1))
      @mpiassertEqual(N, size(a1d))
    end do

    allocate (a2d(N, N), source=real(myid, kind(a2d)))
    groups(1) = trim(fname)
    groups(2) = 'rank'//itoa(myid)//'_kind_'//itoa(kind(a2d))
    groups(3) = 'a2d'

    do rank = 0, numnodes - 1
      if (rank .eq. myid) then
        call ncwrite(groups, a2d, ierr); call CHKERR(ierr, 'Could not write 2d array to nc file')
      end if
      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    end do
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    do rank = 0, numnodes - 1
      deallocate (a2d)
      groups(2) = 'rank'//itoa(rank)//'_kind_'//itoa(kind(a2d))
      call ncload(groups, a2d, ierr); call CHKERR(ierr, 'Could not read 2d array from nc file: '//char_arr_to_str(groups, ' / '))
      @mpiassertEqual(real(rank, kind(a2d)), a2d(1, 1))
      @mpiassertEqual(N, size(a2d, dim=1))
      @mpiassertEqual(N, size(a2d, dim=2))
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_netcdf_write_hyperslab_1d(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr
    integer(mpiint) :: i
    real(ireals), allocatable :: a1d(:)

    integer(mpiint), parameter :: N = 10
    character(len=default_str_len) :: groups(2)

    groups(1) = 'pfunit_hyperslab_test.nc'
    groups(2) = 'hyperslab1d'

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    allocate (a1d(N))
    do i = 1, N
      a1d(i) = i
    end do

    call ncwrite(groups, a1d(1:N / 2), ierr, &
                 arr_shape=shape(a1d), startp=[1], countp=[N / 2])

    call ncwrite(groups, a1d(N / 2 + 1:N), ierr, &
                 arr_shape=shape(a1d), startp=[N / 2 + 1], countp=[N / 2])

    deallocate (a1d)

    call ncload(groups, a1d, ierr)

    do i = 1, N
      @assertEqual(i, a1d(i))
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_netcdf_write_hyperslab_2d(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid, ierr
    integer(mpiint) :: i, j
    real(ireals), allocatable :: a2d(:, :)

    integer(mpiint), parameter :: N = 10
    character(len=default_str_len) :: groups(2)

    groups(1) = 'pfunit_hyperslab_test.nc'
    groups(2) = 'hyperslab2d'

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    allocate (a2d(N, N))
    do i = 1, N
      a2d(:, i) = i
    end do

    do i = 1, N
      call ncwrite(groups, a2d(i, 1:N), ierr, &
                   arr_shape=shape(a2d), startp=[1, i], countp=[N])
    end do

    deallocate (a2d)

    call ncload(groups, a2d, ierr)

    do j = 1, N
      do i = 1, N
        @assertEqual(i, a2d(i, j))
      end do
    end do
  end subroutine

  @test(npes=[1])
  subroutine test_get_set_global_attributes(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid
    real(real32) :: attr_r32
    real(real64) :: attr_r64
    integer(int32) :: attr_i32
    integer(int64) :: attr_i64
    character(len=default_str_len) :: attr_str
    integer(mpiint) :: ierr

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    attr_r32 = 7._real32
    call set_global_attribute(fname, 'r32_test_attr', attr_r32, ierr)
    call get_global_attribute(fname, 'r32_test_attr', attr_r32, ierr)
    @assertEqual(7._real32, attr_r32)

    attr_r64 = 7._real64
    call set_global_attribute(fname, 'r64_test_attr', attr_r64, ierr)
    call get_global_attribute(fname, 'r64_test_attr', attr_r64, ierr)
    @assertEqual(7._real64, attr_r64)

    attr_i32 = 8_int32
    call set_global_attribute(fname, 'i32_test_attr', attr_i32, ierr)
    call get_global_attribute(fname, 'i32_test_attr', attr_i32, ierr)
    @assertEqual(8_int32, attr_i32)

    attr_i64 = 9_int64
    call set_global_attribute(fname, 'i64_test_attr', attr_i64, ierr)
    call get_global_attribute(fname, 'i64_test_attr', attr_i64, ierr)
    @assertEqual(9_int64, attr_i64)

    attr_str = 'this is a test string'
    call set_global_attribute(fname, 'str_test_attr', attr_str, ierr)
    call get_global_attribute(fname, 'str_test_attr', attr_str, ierr)
    @assertEqual('this is a test string', trim(attr_str))
  end subroutine
end module
