module test_helper_functions
  use iso_c_binding
  use m_data_parameters, only: ireals, iintegers, mpiint, init_mpi_data_parameters
  use m_helper_functions, only : imp_bcast, imp_allgather_int_inplace, mpi_logical_and, mpi_logical_or, &
    compute_normal_3d, hit_plane, pnt_in_triangle, norm, distance_to_edge, determine_normal_direction, &
    cumprod, reverse, rotation_matrix_around_axis_vec, deg2rad, char_arr_to_str, cstr, &
    search_sorted_bisection

  use pfunit_mod

  implicit none

contains

@test(npes =[2,1])
subroutine test_mpi_functions(this)
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

    real(ireals),pointer :: bcast_2d_arr_ptr(:,:)=>NULL()

    logical :: l_all_true, l_all_false, l_even_true

    integer(iintegers),parameter :: repetitions=10000
    integer(iintegers) :: rep

    integer(c_size_t) :: large_size_t_int

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    call init_mpi_data_parameters(comm)

    do rep=1,repetitions
      if(.not.allocated(arr)) allocate(arr(numnodes), source=-1_iintegers)
      arr(myid+1) = myid

      call imp_allgather_int_inplace(comm, arr)

      do i=1,numnodes
        @assertEqual(i-1, arr(i))
      enddo

      ! Check if c_size_t broadcasts work
      if(myid.eq.0) then
        large_size_t_int = 8589934592_c_size_t ! 2**33
      endif
      call imp_bcast(comm, large_size_t_int, 0_mpiint)
      @assertEqual(8589934592_c_size_t, large_size_t_int, 'c_size_t broadcast wrong')

      ! Check if scalar and array bcasts work
      if(myid.eq.0) then
        bcast_scalar = 1234
        bcast_real_scalar = 1234._ireals
        if(.not.allocated(bcast_arr))    allocate(bcast_arr(10), source=1234_iintegers)
        if(.not.allocated(bcast_1d_arr)) allocate(bcast_1d_arr(2), source=1234._ireals)
        if(.not.allocated(bcast_2d_arr)) allocate(bcast_2d_arr(2,1), source=1234._ireals)
        if(.not.allocated(bcast_3d_arr)) allocate(bcast_3d_arr(2,1,1), source=1234._ireals)
        if(.not.allocated(bcast_5d_arr)) allocate(bcast_5d_arr(2,1,1,1,1), source=1234._ireals)

        if(.not.associated(bcast_2d_arr_ptr)) allocate(bcast_2d_arr_ptr(2,1), source=1234._ireals)
      else
        bcast_scalar = -1
        bcast_real_scalar = -1
        if(allocated(bcast_arr))    deallocate(bcast_arr)
        if(allocated(bcast_1d_arr)) deallocate(bcast_1d_arr)
        if(allocated(bcast_2d_arr)) deallocate(bcast_2d_arr)
        if(allocated(bcast_3d_arr)) deallocate(bcast_3d_arr)
        if(allocated(bcast_5d_arr)) deallocate(bcast_5d_arr)

        if(associated(bcast_2d_arr_ptr)) deallocate(bcast_2d_arr_ptr)
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

      ! Check pointer Bcasts:
      call imp_bcast(comm, bcast_2d_arr_ptr, 0_mpiint)
      @assertTrue(associated(bcast_2d_arr_ptr), 'real pointer array broadcast wrong, allocation failed')
      @assertEqual(1234, bcast_2d_arr_ptr(1,1), 'real pointer array broadcast wrong')
      @assertEqual(1234, bcast_2d_arr_ptr(2,1), 'real pointer array broadcast wrong')


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
    enddo ! repetitions
end subroutine

@test(npes =[1])
subroutine test_triangle_functions(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals),parameter :: zero=0, one=1, dx = 100
    real(ireals),parameter :: A(2) = [zero, zero]
    real(ireals),parameter :: B(2) = [dx, zero]
    real(ireals),parameter :: C(2) = [dx/2.,sqrt(dx**2 - (dx/2)**2)]
    real(ireals) :: P(2), distance

    real(ireals), dimension(3) :: normal, new_loc, center_face, center_cell
    integer(iintegers) :: normal_direction

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

    ! test routines that determine direction of normal
    normal = [zero, zero, one]
    center_face = [zero, zero, zero]
    center_cell = [zero, zero, one]
    normal_direction = determine_normal_direction(normal, center_face, center_cell)
    @assertEqual(1_iintegers, normal_direction, 'direction of normal not towards cell center (case 1)')

    center_cell = [zero, zero, -one]
    normal_direction = determine_normal_direction(normal, center_face, center_cell)
    @assertEqual(-1_iintegers, normal_direction, 'direction of normal not towards cell center (case 2)')
end subroutine

@test(npes=[1])
subroutine test_cumprod(this)
  class (MpiTestMethod), intent(inout) :: this
  integer(iintegers),parameter :: iarr(3) = [1,2,3]
  real(ireals),parameter :: rarr(3) = [1,2,3]
  @assertEqual([1,2,6], cumprod(iarr))
  @assertEqual(real([1,2,6], ireals), cumprod(rarr))
end subroutine

@test(npes=[1])
subroutine test_reverse(this)
  class (MpiTestMethod), intent(inout) :: this
  real(ireals),parameter :: arr(3) = [1,2,3]
  real(ireals) :: arr2d(2,3)
  @assertEqual([3,2,1], reverse(arr))

  arr2d(1,:) = arr
  arr2d(2,:) = arr+1

  arr2d = reverse(arr2d)

  @assertEqual(arr+1, arr2d(1,:))
  @assertEqual(arr, arr2d(2,:))

  arr2d = reverse(arr2d, dim=2_iintegers)

  @assertEqual(reverse(arr+1), arr2d(1,:))
  @assertEqual(reverse(arr), arr2d(2,:))
end subroutine

@test(npes=[1])
subroutine test_rotation_matrix_around_axis_vec(this)
  class (MpiTestMethod), intent(inout) :: this
  real(ireals), dimension(3) :: ex, ey, ez, x1
  real(ireals) :: rot_angle, Mrot(3,3)
  real(ireals), parameter :: eps=sqrt(epsilon(eps))
  integer(iintegers) :: i

  ex = [1,0,0]
  ey = [0,1,0]
  ez = [0,0,1]
  x1 = [1,0,0]

  Mrot = rotation_matrix_around_axis_vec(deg2rad(180._ireals), ez)
  @assertEqual(-x1, matmul(Mrot, x1), eps)

  Mrot = rotation_matrix_around_axis_vec(deg2rad(360._ireals), ez)
  @assertEqual(x1, matmul(Mrot, x1), eps)

  do i = 0, 360
    rot_angle = real(i, ireals)
    Mrot = rotation_matrix_around_axis_vec(deg2rad(rot_angle), ex)
    @assertEqual(x1, matmul(Mrot, x1), eps)
  enddo

  Mrot = rotation_matrix_around_axis_vec(deg2rad(90._ireals), ey)
  @assertEqual(-ez, matmul(Mrot, x1), eps)

  Mrot = rotation_matrix_around_axis_vec(deg2rad(270._ireals), ey)
  @assertEqual(ez, matmul(Mrot, x1), eps)
end subroutine

@test(npes=[1])
subroutine test_char_arr_to_str(this)
  class (MpiTestMethod), intent(inout) :: this
  character(len=4)  :: a(3)
  a(1) = '1'
  a(2) = '23'
  a(3) = '456'

  @assertEqual('1, 23, 456', char_arr_to_str(a))
  @assertEqual('1 :: 23 :: 456', char_arr_to_str(a, ' :: '))
end subroutine

@test(npes=[1])
subroutine test_search_sorted_bisection(this)
  class (MpiTestMethod), intent(inout) :: this
  real(ireals), parameter :: A(3) = [-10, 0, 2]
  !real(ireals) :: x

  @assertEqual(1.0, search_sorted_bisection(A, -20._ireals))
  @assertEqual(1.0, search_sorted_bisection(A, -10._ireals))
  @assertEqual(1.5, search_sorted_bisection(A,  -5._ireals))
  @assertEqual(2.0, search_sorted_bisection(A,  -0._ireals))
  @assertEqual(2.0, search_sorted_bisection(A,   0._ireals))
  @assertEqual(2.5, search_sorted_bisection(A,   1._ireals))
  @assertEqual(3.0, search_sorted_bisection(A,   2._ireals))
  @assertEqual(3.0, search_sorted_bisection(A,   3._ireals))
end subroutine
end module
