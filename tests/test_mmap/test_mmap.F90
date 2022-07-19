module test_mmap
  use m_data_parameters, only: irealLUT, iintegers, mpiint, init_mpi_data_parameters
  use m_mmap, only: arr_to_mmap, munmap_mmap_ptr, arr_to_binary_datafile, binary_file_to_mmap, load_mmap_array

  use m_helper_functions, only: CHKERR

  use pfunit_mod

  implicit none

contains

  @test(npes=[2, 1])
  subroutine test_mpi_write_and_read_mmaps(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(iintegers), parameter :: N = 3
    character(len=*), parameter :: fname = 'pfunit_test.mmap'

    integer(iintegers) :: i, j
    integer(mpiint) :: numnodes, comm, myid, ierr

    real(irealLUT), pointer :: arr(:, :)
    real(irealLUT), pointer :: mmap_ptr(:, :)

    arr => null()
    mmap_ptr => null()

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    print *, 'Testing mmap with ', numnodes, 'ranks ...'

    call init_mpi_data_parameters(comm)

    if (myid .eq. 0) then
      allocate (arr(N, N))
      do i = 1, N
        do j = 1, N
          arr(i, j) = real((i - 1) * N + j, irealLUT)
        end do
      end do
    end if

    call arr_to_mmap(comm, fname, mmap_ptr, ierr, arr)
    call CHKERR(ierr, 'Err creating mmap')
    if (associated(arr)) deallocate (arr)

    arr => mmap_ptr
    do i = 1, N
      do j = 1, N
        @assertEqual(real((i - 1) * N + j, irealLUT), arr(i, j))
      end do
    end do

    call munmap_mmap_ptr(arr, ierr); call CHKERR(ierr, 'Err unmapping mmap')
    print *, 'Testing mmap with ', numnodes, 'ranks ... done'
  end subroutine

  @test(npes=[2, 1])
  subroutine test_write_read_mmap(this)
    class(MpiTestMethod), intent(inout) :: this

    integer(iintegers), parameter :: N = 5
    character(len=*), parameter :: fname = 'pfunit_test.mmap'

    integer(iintegers) :: i, j
    integer(mpiint) :: numnodes, comm, myid, ierr

    real(irealLUT), pointer :: arr(:, :)
    real(irealLUT), pointer :: mmap_ptr(:, :)

    arr => null()
    mmap_ptr => null()

    comm = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid = this%getProcessRank()

    print *, 'Testing mmap with ', numnodes, 'ranks ...'

    call init_mpi_data_parameters(comm)

    if (myid .eq. 0) then
      allocate (arr(N, N))
      do i = 1, N
        do j = 1, N
          arr(i, j) = real((i - 1) * N + j, irealLUT)
        end do
      end do
      call arr_to_binary_datafile(arr, fname, ierr); call CHKERR(ierr)
    end if

    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    call load_mmap_array('<a non existing file>', mmap_ptr, ierr)
    @assertTrue(ierr .ne. 0)

    call load_mmap_array(fname, mmap_ptr, ierr)

    do i = 1, N
      do j = 1, N
        @assertEqual(real((i - 1) * N + j, irealLUT), mmap_ptr(i, j))
      end do
    end do

    call munmap_mmap_ptr(mmap_ptr, ierr); call CHKERR(ierr, 'Err unmapping mmap')
  end subroutine

end module
