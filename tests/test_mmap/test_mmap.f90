module test_mmap
    use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters
    use m_mmap, only: arr_to_mmap, munmap_mmap_ptr

    use m_helper_functions, only: CHKERR

    use pfunit_mod

  implicit none

  contains

@test(npes =[2,1])
subroutine test_mpi_functions(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers), parameter :: N=3
    character(len=*), parameter :: fname='pfunit_test.mmap'

    integer(iintegers) :: i, j
    integer(mpiint) :: numnodes, comm, myid, ierr

    real(ireals), pointer :: arr(:,:)
    real(ireals), pointer :: mmap_ptr(:,:)

    arr=>NULL()
    mmap_ptr=>NULL()

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    print *,'Testing mmap with ', numnodes,'ranks ...'

    call init_mpi_data_parameters(comm)

    if(myid.eq.0) then
      allocate(arr(N,N))
      do i = 1, N
        do j = 1, N
          arr(i,j) = (i-1)*N+j
        enddo
      enddo
    endif

    call arr_to_mmap(comm, fname, mmap_ptr, ierr, arr)
    call CHKERR(ierr, 'Err creating mmap')
    if(associated(arr)) deallocate(arr)

    arr => mmap_ptr
    do i = 1, N
      do j = 1, N
        @assertEqual((i-1)*N+j, arr(i,j))
      enddo
    enddo

    call munmap_mmap_ptr(arr, ierr)
    call CHKERR(ierr, 'Err unmapping mmap')
    print *,'Testing mmap with ', numnodes,'ranks ... done'
end subroutine

end module
