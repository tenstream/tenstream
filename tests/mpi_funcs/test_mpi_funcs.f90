@test(npes =[1,2,4])
subroutine test_dtypes(this)

    use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters
    use m_helper_functions, only : imp_bcast, imp_allgather_int_inplace, mpi_logical_and, mpi_logical_or

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid, i

    integer(iintegers),allocatable :: arr(:)
    integer(iintegers),allocatable :: bcast_arr(:)
    integer(iintegers) :: bcast_scalar

    real(ireals),allocatable :: bcast_1d_arr(:)
    real(ireals),allocatable :: bcast_2d_arr(:,:)
    real(ireals),allocatable :: bcast_3d_arr(:,:,:)
    real(ireals),allocatable :: bcast_5d_arr(:,:,:,:,:)
    real(ireals) :: bcast_real_scalar

    logical :: l_all_true, l_all_false, l_even_true

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

    ! Check if scalar and array bcasts work
    if(myid.eq.0) then
      bcast_scalar = 1234
      bcast_real_scalar = 1234._ireals
      allocate(bcast_arr(10), source=1234)
      allocate(bcast_1d_arr(2), source=1234._ireals)
      allocate(bcast_2d_arr(2,1), source=1234._ireals)
      allocate(bcast_3d_arr(2,1,1), source=1234._ireals)
      allocate(bcast_5d_arr(2,1,1,1,1), source=1234._ireals)
    else
      bcast_scalar = -1
      bcast_real_scalar = -1
    endif

    call imp_bcast(comm, bcast_arr, 0_mpiint)
    @assertTrue(allocated(bcast_arr), 'int array broadcast wrong, allocation failed')
    do i=1,size(bcast_arr)
      @assertEqual(1234, bcast_arr(i), 'int array broadcast wrong')
    enddo

    call imp_bcast(comm, bcast_1d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_1d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_1d_arr(1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_1d_arr(2), 'real array broadcast wrong')

    call imp_bcast(comm, bcast_2d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_2d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_2d_arr(1,1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_2d_arr(2,1), 'real array broadcast wrong')

    call imp_bcast(comm, bcast_3d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_3d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_3d_arr(1,1,1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_3d_arr(2,1,1), 'real array broadcast wrong')

    call imp_bcast(comm, bcast_5d_arr, 0_mpiint)
    @assertTrue(allocated(bcast_5d_arr), 'real array broadcast wrong, allocation failed')
    @assertEqual(1234, bcast_5d_arr(1,1,1,1,1), 'real array broadcast wrong')
    @assertEqual(1234, bcast_5d_arr(2,1,1,1,1), 'real array broadcast wrong')


    ! Scalar Bcasts:
    call imp_bcast(comm, bcast_scalar, 0_mpiint)
    @assertEqual(1234, bcast_scalar, 'int scalar broadcast wrong')

    call imp_bcast(comm, bcast_real_scalar, 0_mpiint)
    @assertEqual(1234, bcast_real_scalar, 'real scalar broadcast wrong')

    ! Logical Bcasts:
    if (myid.eq.0) then
      l_all_true = .True.
    else
      l_all_true = .False.
    endif
    call imp_bcast(comm, l_all_true, 0_mpiint)
    @assertEqual(.True., l_all_true, 'logical bcast wrong')


    ! Check for logical reductions
    l_all_true  = .True.
    l_all_false = .False.

    if (modulo(myid, 2_mpiint).eq.0) then
      l_even_true = .True.
    else
      l_even_true = .False.
    endif

    @assertEqual(.True., mpi_logical_and(comm, l_all_true), 'mpi_logical_and is wrong')
    @assertEqual(.False., mpi_logical_and(comm, l_all_false), 'mpi_logical_and is wrong')

    if (numnodes.gt.1) then
      @assertEqual(.False., mpi_logical_and(comm, l_even_true), 'mpi_logical_and is wrong')
    else
      @assertEqual(.True., mpi_logical_and(comm, l_even_true), 'mpi_logical_and is wrong')
    endif

    @assertEqual(.True., mpi_logical_or(comm, l_all_true), 'mpi_logical_or is wrong')
    @assertEqual(.False., mpi_logical_or(comm, l_all_false), 'mpi_logical_or is wrong')

    @assertEqual(.True., mpi_logical_or(comm, l_even_true), 'mpi_logical_or is wrong')


end subroutine

