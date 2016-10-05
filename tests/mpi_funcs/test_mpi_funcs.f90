@test(npes =[1,2,4])
subroutine test_dtypes(this)

    use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters
    use m_helper_functions, only : imp_allgather_int_inplace

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid, i

    integer(iintegers),allocatable :: arr(:)

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    allocate(arr(numnodes), source=-1)

    arr(myid+1) = myid

    call imp_allgather_int_inplace(comm, arr)

    do i=1,numnodes
      @assertEqual(i-1, arr(i))
    enddo

end subroutine

