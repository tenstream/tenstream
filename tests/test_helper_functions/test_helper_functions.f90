@test(npes =[1,2,4])
subroutine test_mpi_functions(this)

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

@test(npes =[1])
subroutine test_triangle_functions()

    use m_data_parameters, only: ireals
    use m_helper_functions, only : compute_normal_3d, hit_plane, pnt_in_triangle, norm, distance_to_edge

    use pfunit_mod

    implicit none
    real(ireals),parameter :: zero=0, one=1, dx = 100
    real(ireals),parameter :: A(2) = [zero, zero]
    real(ireals),parameter :: B(2) = [dx, zero]
    real(ireals),parameter :: C(2) = [dx/2.,sqrt(dx**2 - (dx/2)**2)]
    real(ireals) :: P(2), distance

    real(ireals) :: normal(3), new_loc(3)

    ! Tests determining the distance of a point to a 2D line/edge
    @assertEqual(zero, distance_to_edge(A,B,A), 'from point on line, the distance to same line should be zero distance_to_edge1')
    @assertEqual(zero, distance_to_edge(A,B,[dx/2,zero]), 'from point on line, the distance to same line should be zero distance_to_edge2')
    @assertEqual(one, distance_to_edge(A,B,[dx/2,one]), 'here point line <-> distance should be one distance_to_edge3')
    @assertEqual(sqrt(epsilon(dx)), distance_to_edge(A,B,[dx/2,sqrt(epsilon(dx))]), 'here point line <-> distance should be different distance_to_edge4')
    @assertEqual(epsilon(dx), distance_to_edge(A,B,[dx/2,epsilon(dx)]), 'here point line <-> distance should be different distance_to_edge5')
    @assertEqual(epsilon(dx), distance_to_edge(A,B,[dx/2,-epsilon(dx)]), 'here point line <-> distance should be different distance_to_edge6')
    @assertEqual(one, distance_to_edge(A,B,B+[zero,one]), 'here point line <-> distance should be one test distance_to_edge7')


    ! Checks if points lie in a triangle
    new_loc = [0.38475394248962402_ireals, zero, zero]
    @assertTrue(pnt_in_triangle(A,B,C, [new_loc(1), new_loc(2)]), 'custom edge case point should be in triangle!')

    normal = compute_normal_3d([A(1),A(2),zero], [B(1),B(2),zero], [C(1),C(2),zero])
    @assertEqual([zero,zero,one], normal, '3D normal not as expected')

    normal = compute_normal_3d([A(1),A(2),zero], [C(1),C(2),zero], [B(1),B(2),zero])
    @assertEqual([zero,zero,one], -normal, '3D normal not as expected')

    @assertEqual(one, norm(normal), 'returned normal is not normed to one')

    ! Check if we can determine if a point is in a triangle
    @assertTrue(pnt_in_triangle(A,B,C, A), 'pnt_in_triangle wrong for edge case in A')
    @assertTrue(pnt_in_triangle(A,B,C, B), 'pnt_in_triangle wrong for edge case in B')
    @assertTrue(pnt_in_triangle(A,B,C, C), 'pnt_in_triangle wrong for edge case in C')
    @assertTrue(pnt_in_triangle(A,B,C, [one/2, one/2]), 'pnt_in_triangle wrong for center of triangle')

    @assertTrue(pnt_in_triangle(A,B,C, A+(C-A)/2), 'pnt_in_triangle wrong for edge case on line between A and C')
    @assertTrue(pnt_in_triangle(A,B,C, A+(B-A)/2), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertTrue(pnt_in_triangle(A,B,C, C+(B-C)/2), 'pnt_in_triangle wrong for edge case on line between B and C')

    @assertFalse(pnt_in_triangle(A,B,C, A-[one,one ]), 'pnt_in_triangle wrong for outside case 1')
    @assertFalse(pnt_in_triangle(A,B,C, B+[one,zero]), 'pnt_in_triangle wrong for outside case 2')
    @assertFalse(pnt_in_triangle(A,B,C, C+[one,one] ), 'pnt_in_triangle wrong for outside case 3')


    ! vector from C to pnt halfway between (AB):
    P = A+(B-A)/2 - C
    @assertTrue(pnt_in_triangle(A,B,C, C+P), 'pnt_in_triangle wrong for edge case on line between A and B')
    @assertFalse(pnt_in_triangle(A,B,C, C+(one+sqrt(epsilon(one)))*P), 'pnt_in_triangle wrong for edge case epsilon after line between A and B')

    ! Check if distance caluclations are OK
    distance = hit_plane([A(1),A(2),one], [zero,zero,-one], [A(1),A(2),zero], normal)
    @assertEqual(one, distance, 'distance calculation not correct 1')

    distance = hit_plane([A(1),A(2),one], [zero,zero,+one], [A(1),A(2),zero], normal)
    @assertEqual(-one, distance, 'distance calculation not correct 2')


    distance = hit_plane([A(1),A(2),one], [zero,zero,-one], [C(1),C(2),zero], normal)
    @assertEqual(one, distance, 'distance calculation not correct 3')

    distance = hit_plane([A(1),A(2),one], [zero,zero,+one], [C(1),C(2),zero], normal)
    @assertEqual(-one, distance, 'distance calculation not correct 4')
end subroutine
