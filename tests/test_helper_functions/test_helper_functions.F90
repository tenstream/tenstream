module test_helper_functions
  use iso_fortran_env, only: REAL32, REAL64
  use iso_c_binding
  use m_data_parameters, only: ireals, ireal128, iintegers, mpiint, init_mpi_data_parameters
  use m_helper_functions, only : &
    & approx, &
    & char_arr_to_str, &
    & compute_normal_3d, &
    & cstr, &
    & cumprod, &
    & deg2rad, &
    & determine_normal_direction, &
    & distance_to_edge, &
    & expm1, &
    & imp_allgather_int_inplace, &
    & imp_allreduce_mean, &
    & imp_allreduce_sum, &
    & imp_bcast, &
    & imp_reduce_mean, &
    & imp_reduce_sum, &
    & ind_1d_to_nd, &
    & ind_nd_to_1d, &
    & is_between, &
    & is_inrange, &
    & mpi_logical_all_same, &
    & mpi_logical_and, &
    & mpi_logical_or, &
    & ndarray_offsets, &
    & normalize_vec, &
    & read_ascii_file_2d, &
    & resize_arr, &
    & reverse, &
    & rotation_matrix_around_axis_vec, &
    & rotation_matrix_local_basis_to_world, &
    & rotation_matrix_world_to_local_basis, &
    & solve_quadratic, &
    & cartesian_2_spherical, &
    & spherical_2_cartesian, &
    & toStr

  use pfunit_mod

  implicit none

contains

@before
subroutine setup(this)
class (MpiTestMethod), intent(inout) :: this
  integer(mpiint) :: comm
  comm     = this%getMpiCommunicator()
  call init_mpi_data_parameters(comm)
end subroutine setup

@after
subroutine teardown(this)
class (MpiTestMethod), intent(inout) :: this
  logical :: lpetsc_is_initialized
  integer(mpiint) :: ierr
  call PetscInitialized(lpetsc_is_initialized, ierr)
  if(lpetsc_is_initialized) call PetscFinalize(ierr)
end subroutine teardown
@test(npes =[1,2,3])
subroutine test_mpi_char_bcast(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(mpiint) :: numnodes, comm, myid

    character(len=10) :: s, s_target
    character(len=10), allocatable :: s1d(:)

    integer(iintegers), parameter :: N=10
    integer(iintegers), parameter :: repetitions=100
    integer(iintegers) :: rep, i

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    if(myid.eq.0) allocate(s1d(N))

    do rep = 1, repetitions
      ! single char
      s_target = "Rep"//toStr(rep)
      if (myid.eq.0) then
        s = s_target
      endif
      call imp_bcast(comm, s, 0_mpiint)
      @mpiAssertEqual(s_target, s)

      ! array of strings
      if(myid.eq.0) then
        do i = 1, size(s1d)
          s1d(i) = "Rep"//toStr(rep)//':'//toStr(i)
        enddo
      endif
      call imp_bcast(comm, s1d, 0_mpiint)
      do i = 1, size(s1d)
        s_target = "Rep"//toStr(rep)//':'//toStr(i)
        @mpiAssertEqual(s_target, s1d(i))
      enddo

    enddo
end subroutine

@test(npes =[1,2,3])
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

      ! check for `all the same` operations
      @assertTrue(mpi_logical_all_same(comm, l_all_true), 'mpi_logical_all_same(l_all_true) is wrong')
      @assertTrue(mpi_logical_all_same(comm, l_all_false), 'mpi_logical_all_same(l_all_false) is wrong')
      if(numnodes.gt.1) then
        @assertFalse(mpi_logical_all_same(comm, l_even_true), 'mpi_logical_all_same(l_even_true) is wrong')
      else
        @assertTrue(mpi_logical_all_same(comm, l_even_true), 'mpi_logical_all_same(l_all_true) is wrong')
      endif
    enddo ! repetitions
end subroutine

@test(npes =[1,2,3])
subroutine test_mpi_reductions(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: numnodes, comm, myid
    real(ireals)   :: v_ireals, sireals
    real(ireal128) :: v_128, s128
    real(ireals)  , allocatable :: m_ireals(:), mean
    real(ireal128), allocatable :: m_128   (:)
    integer(iintegers) :: i, s

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    v_ireals = real(myid, ireals)
    v_128    = real(myid, ireal128)

    sireals = (numnodes-1)*numnodes/2
    s128    = (numnodes-1)*numnodes/2
    print *,myid, numnodes, 'v128', v_128

    call imp_reduce_sum(comm, v_ireals)
    call imp_reduce_sum(comm, v_128   )
    if(myid.eq.0) then
      @assertEqual(sireals, v_ireals, epsilon(v_ireals), 'ireals reduce_sum is not correct')
      sireals  = real(s128 , ireals)
      v_ireals = real(v_128, ireals)
      @assertEqual(sireals, v_ireals, epsilon(v_ireals), '128 bit reduce_sum is not correct Nranks:'//toStr(numnodes))
    endif

    call imp_allreduce_sum(comm, real(myid, ireals  ), v_ireals)
    @assertEqual(sireals, v_ireals, epsilon(v_ireals), 'all reduce sum is not correct')


    ! allreduce_mean scalar
    call imp_allreduce_mean(comm, real(myid, ireals), sireals)
    @assertEqual(real(numnodes-1, ireals)/2, sireals, epsilon(v_ireals), 'ireals reduce_mean scalar is not correct')

    allocate(m_ireals(myid+1)); m_ireals = real(myid, kind(m_ireals))
    allocate(m_128   (myid+1)); m_128    = real(myid, kind(m_128   ))

    mean = 0
    s = 0
    do i=1, numnodes
      s    = s + i
      mean = mean + real(i*(i-1), ireals) ! i entries with the value i-1 in them
    enddo
    mean = mean / real(s, ireals)

    call imp_allreduce_mean(comm, m_ireals, sireals)
    @assertEqual(mean, sireals, epsilon(v_ireals), 'ireals reduce_mean is not correct Nranks:'//toStr(numnodes))
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
subroutine test_rotation_matrix_world_to_local(this)
  class (MpiTestMethod), intent(inout) :: this
  real(ireals), dimension(3) :: ex, ey, ez, x1
  real(ireals) :: Mrot(3,3)
  real(ireals), parameter :: eps=sqrt(epsilon(eps))

  ex = [1,0,0]
  ey = [0,1,0]
  ez = [0,0,1]
  x1 = [1,2,3]

  ! because ex, ey, ez are just world coords, this transformation does not do anything
  Mrot = rotation_matrix_world_to_local_basis(ex, ey, ez)
  @assertEqual(x1, matmul(Mrot, x1), eps)

  ! flip the coord system
  Mrot = rotation_matrix_world_to_local_basis(-ex, -ey, -ez)
  @assertEqual(x1, -matmul(Mrot, x1), eps)


  ex = [0,0,1]
  ey = [0,1,0]
  ez = [1,0,0]
  x1 = [1,2,3]

  Mrot = rotation_matrix_world_to_local_basis(ex, ey, ez)
  @assertEqual([3,2,1], matmul(Mrot, x1), eps)


  ! have a local basis with 45 deg rotation in the xy-plane
  ex = [ 1,1,0]; ex = ex/norm2(ex)
  ey = [-1,1,0]; ey = ey/norm2(ey)
  ez = [ 0,0,1]; ez = ez/norm2(ez)
  x1 = [ 1, 1, 1]

  Mrot = rotation_matrix_world_to_local_basis(ex, ey, ez)
  @assertEqual([ sqrt(2._ireals),0._ireals,1._ireals], matmul(Mrot, x1), eps)
end subroutine

@test(npes=[1])
subroutine test_rotation_matrix_local_basis_to_world(this)
  class (MpiTestMethod), intent(inout) :: this
  real(ireals), dimension(3) :: ex, ey, ez, x1
  real(ireals) :: Mrot(3,3)
  real(ireals), parameter :: eps=sqrt(epsilon(eps))

  ex = [1,0,0]
  ey = [0,1,0]
  ez = [0,0,1]
  x1 = [1,2,3]

  ! because ex, ey, ez are the same as world coords, this transformation does not do anything
  Mrot = rotation_matrix_local_basis_to_world(ex, ey, ez)
  @assertEqual(x1, matmul(Mrot, x1), eps)

  ! flip the coord system
  Mrot = rotation_matrix_local_basis_to_world(-ex, -ey, -ez)
  @assertEqual(x1, -matmul(Mrot, x1), eps)


  ! have a local basis with 45 deg rotation in the xy-plane
  ex = [ 1,1,0]; ex = ex/norm2(ex)
  ey = [-1,1,0]; ey = ey/norm2(ey)
  ez = [ 0,0,1]; ez = ez/norm2(ez)
  x1 = [ sqrt(2._ireals),0._ireals,1._ireals]

  Mrot = rotation_matrix_local_basis_to_world(ex, ey, ez)
  @assertEqual([1,1,1], matmul(Mrot, x1), eps)
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
subroutine test_solve_quadratic(this)
  class (MpiTestMethod), intent(inout) :: this
  integer(mpiint) :: ierr
  real(ireals) :: x(2)
  real(ireals), parameter :: eps=sqrt(epsilon(eps))

  call solve_quadratic(1._ireals,-2._ireals,-3._ireals, x, ierr)
  @assertEqual(0, ierr)
  @assertEqual([-1._ireals, 3._ireals], x, eps)

  call solve_quadratic(2._ireals, 4._ireals, -4._ireals, x, ierr)
  @assertEqual(0, ierr)
  @assertEqual([-1._ireals - sqrt(3._ireals), -1._ireals + sqrt(3._ireals)], x, eps)

  call solve_quadratic(1._ireals, 2._ireals, 3._ireals, x, ierr)
  @assertFalse(ierr.eq.0)

  call solve_quadratic(1._ireals, 2._ireals, 1._ireals, x, ierr)
  @assertEqual(0, ierr)
  @assertEqual([-1._ireals, -1._ireals], x, eps)
end subroutine

@test(npes=[1])
subroutine test_is_between(this)
  class (MpiTestMethod), intent(inout) :: this

  @assertTrue(is_between(0.5_ireals, 0.0_ireals, 1.0_ireals))
  @assertTrue(is_between(0.5_ireals, 1.0_ireals, 0.0_ireals))
  @assertTrue(is_between(-0.5_ireals, 0.0_ireals, -1.0_ireals))
  @assertTrue(is_between(-0.5_ireals, -1.0_ireals, 0.0_ireals))

  @assertFalse(is_between(-0.5_ireals, 0.0_ireals, 1.0_ireals))
  @assertFalse(is_between(-0.5_ireals, 1.0_ireals, 0.0_ireals))
  @assertFalse(is_between(1.5_ireals, 0.0_ireals, 1.0_ireals))
  @assertFalse(is_between(1.5_ireals, 1.0_ireals, 0.0_ireals))

  @assertTrue(is_between(0.5_ireals, 1.0_ireals, 0.0_ireals))

  @assertTrue(is_between(780,0,2875391))

  @assertFalse(is_between(0,0,2))
  @assertTrue (is_between(1,0,2))
  @assertFalse(is_between(2,0,2))

  @assertFalse(is_between(0,2,0))
  @assertTrue (is_between(1,2,0))
  @assertFalse(is_between(2,2,0))
end subroutine

subroutine test_is_inrange(this)
  class (MpiTestMethod), intent(inout) :: this

  @assertTrue(is_inrange(0,0,2))
  @assertTrue(is_inrange(1,0,2))
  @assertTrue(is_inrange(2,0,2))

  @assertTrue(is_inrange(0,2,0))
  @assertTrue(is_inrange(1,2,0))
  @assertTrue(is_inrange(2,2,0))

  @assertFalse(is_inrange(-1,0,2))
  @assertFalse(is_inrange( 4,0,2))

  @assertFalse(is_inrange(-1,2,0))
  @assertFalse(is_inrange( 4,2,0))

  @assertTrue(is_inrange(780,0,2875391))
end subroutine

@test(npes=[1])
subroutine test_normalize_vec(this)
  class (MpiTestMethod), intent(inout) :: this
  real(REAL32) :: v1(2), v2(2)
  real(REAL64) :: p1(2), p2(2)
  integer(mpiint) :: ierr
  v1 = [2,0]
  call normalize_vec(v1, v2, ierr)
  @assertEqual(0_mpiint, ierr)
  @assertEqual([1,0], v2, epsilon(v2))

  call normalize_vec(v1, ierr)
  @assertEqual(0_mpiint, ierr)
  @assertEqual([1,0], v1, epsilon(v1))

  v1 = [0,0]
  call normalize_vec(v1, ierr)
  @assertTrue(0_mpiint.ne.ierr)

  p1 = [2,0]
  call normalize_vec(p1, p2, ierr)
  @assertEqual(0_mpiint, ierr)
  @assertEqual([1,0], p2, epsilon(p2))

  call normalize_vec(p1, ierr)
  @assertEqual(0_mpiint, ierr)
  @assertEqual([1,0], p1, epsilon(p1))

  p1 = [0,0]
  call normalize_vec(p1, ierr)
  @assertTrue(0_mpiint.ne.ierr)
end subroutine

@test(npes=[1])
subroutine test_resize_arr(this)
  class (MpiTestMethod), intent(inout) :: this
  integer(iintegers), allocatable :: i1d(:)
  real(ireals), allocatable :: i2d(:,:), i3d(:,:,:)

  allocate(i1d(3), source=int([0,1,2], kind(i1d)))
  call resize_arr(2_iintegers, i1d)
  @assertEqual([2], shape(i1d), 'wrong dimension after shrinking')

  call resize_arr(1_iintegers, i1d)
  @assertEqual([1], shape(i1d), 'wrong dimension after shrinking')

  call resize_arr(2_iintegers, i1d, fillVal=1_iintegers)
  @assertEqual([2], shape(i1d), 'wrong dimension after shrinking')
  @assertEqual([0,1], i1d, 'wrong values after repeat')

  call resize_arr(4_iintegers, i1d, lrepeat=.True.)
  @assertEqual([4], shape(i1d), 'wrong dimension after shrinking')
  @assertEqual([0,1,0,1], i1d, 'wrong values after repeat')

  call resize_arr(5_iintegers, i1d, lrepeat=.False., fillVal=5_iintegers)
  @assertEqual([0,1,0,1,5], i1d, 'wrong values after repeat')

  allocate(i2d(2,3))
  i2d(:,1) = real([1,2], ireals)
  i2d(:,2) = real([3,4], ireals)
  i2d(:,3) = real([5,6], ireals)

  call resize_arr(1_iintegers, i2d, dim=1)
  @assertEqual([1,3], shape(i2d), 'wrong dimension after shrinking')

  call resize_arr(2_iintegers, i2d, dim=1, lrepeat=.False., fillVal=8._ireals)
  @assertEqual([2,3], shape(i2d), 'wrong dimension after shrinking')
  @assertEqual(8, i2d(2,:), 'wrong values after fill in dim 1')

  call resize_arr(2_iintegers, i2d, dim=2)
  @assertEqual([2,2], shape(i2d), 'wrong dimension after shrinking')

  call resize_arr(3_iintegers, i2d, dim=2, lrepeat=.False., fillVal=9._ireals)
  @assertEqual([2,3], shape(i2d), 'wrong dimension after shrinking')
  @assertEqual(9, i2d(:,3), 'wrong values after fill in dim 2')


  allocate(i3d(2,1,3))
  i3d(:,1,1) = [1,2]
  i3d(:,1,2) = [3,4]
  i3d(:,1,3) = [5,6]

  call resize_arr(1_iintegers, i3d, dim=1)
  @assertEqual([1,1,3], shape(i3d), 'wrong dimension after shrinking')

  call resize_arr(1_iintegers, i3d, dim=3)
  @assertEqual([1,1,1], shape(i3d), 'wrong dimension after shrinking')

  call resize_arr(2_iintegers, i3d, dim=2, lrepeat=.False., fillVal=5._ireals)
  @assertEqual([1,2,1], shape(i3d), 'wrong dimension after shrinking')
  @assertEqual(5, i3d(:,2,:), 'wrong values after fill in dim 2')

  call resize_arr(4_iintegers, i3d, dim=2, lrepeat=.True.)
  @assertEqual([1,4,1], shape(i3d), 'wrong dimension after shrinking')
  @assertEqual(5, i3d(:,2,:), 'wrong values after fill in dim 2')
  @assertEqual(5, i3d(:,4,:), 'wrong values after fill in dim 2')
end subroutine

@test(npes=[1])
subroutine test_approx(this)
  class (MpiTestMethod), intent(inout) :: this
  @assertTrue(     approx(0._ireals, 0.0_ireals, epsilon(0._ireals)))
  @assertTrue(.not.approx(0._ireals, 1.0_ireals, epsilon(0._ireals)))

  @assertTrue(     approx(0._ireals, 0.9_ireals, 1._ireals))
  @assertTrue(.not.approx(0._ireals, 1.1_ireals, 1._ireals))

  @assertTrue(     approx(0._ireals,-0.9_ireals, 1._ireals))
  @assertTrue(.not.approx(0._ireals,-1.1_ireals, 1._ireals))

  @assertTrue(     approx(0._ireals, 1.0_ireals-epsilon(1.0_ireals), 1._ireals))
  @assertTrue(.not.approx(0._ireals, 1.0_ireals+epsilon(1.0_ireals), 1._ireals))
  @assertTrue(     approx(0._ireals,-1.0_ireals+epsilon(1.0_ireals), 1._ireals))
  @assertTrue(.not.approx(0._ireals,-1.0_ireals-epsilon(1.0_ireals), 1._ireals))
end subroutine

@test(npes=[1])
subroutine test_read_ascii_file(this)
  class (MpiTestMethod), intent(inout) :: this
  character(len=*), parameter :: test_fname="tenstream_test_ascii_file.txt"
  real(ireals), allocatable :: arr(:,:)
  integer(mpiint) :: ierr
  integer :: funit

  open (newunit=funit, file=test_fname, action="write", status="replace", iostat=ierr)
  @assertEqual(0, ierr, 'error when trying to write a test ascii file')
  write (funit,*) 1, 2, 3
  write (funit,*)"# skip this line"
  write (funit,*) "     # also skip this one"
  write (funit,*) 4, 5, 6
  write (funit,*) 7, 8, 9
  close (funit)

  call read_ascii_file_2d(test_fname, arr, ierr, verbose=.True.)
  @assertEqual(0_mpiint, ierr)
  @assertEqual(3, size(arr,dim=1), 'line count not correct')
  @assertEqual(3, size(arr,dim=2), 'column count not correct')
  deallocate(arr)

  call read_ascii_file_2d(test_fname, arr, ierr, &
    skiplines=2_iintegers, verbose=.True.)
  @assertEqual(0_mpiint, ierr)
  @assertEqual(2, size(arr,dim=1), 'line count not correct')
  @assertEqual(3, size(arr,dim=2), 'column count not correct')
  deallocate(arr)

  open (newunit=funit, file=test_fname, status="old", position="append", iostat=ierr)
  @assertEqual(0, ierr, 'error when trying to append to a test ascii file')
  write (funit,*) 10, 11, 12, 13
  close (funit)
  call read_ascii_file_2d(test_fname, arr, ierr, verbose=.True.)
  @assertFalse(0_mpiint.eq.ierr)
end subroutine

@test(npes=[1])
subroutine test_ndarray_idx(this)
  class (MpiTestMethod), intent(inout) :: this
  integer(iintegers) :: arr_shape(3)
  integer(iintegers) :: nd_offsets(3)
  integer(iintegers) :: i1d, ind(3)

  arr_shape = [2,3,2]
  call ndarray_offsets(arr_shape, nd_offsets)

  @assertEqual(1_iintegers, nd_offsets(1))
  @assertEqual(2_iintegers, nd_offsets(2))
  @assertEqual(6_iintegers, nd_offsets(3))

  i1d = ind_nd_to_1d(nd_offsets, [integer(iintegers) :: 1, 1, 1])
  @assertEqual(1_iintegers, i1d)

  call ind_1d_to_nd(nd_offsets, i1d, ind)
  @assertEqual(1_iintegers, ind(1))
  @assertEqual(1_iintegers, ind(2))
  @assertEqual(1_iintegers, ind(3))

  i1d = ind_nd_to_1d(nd_offsets, [integer(iintegers) :: 2, 3, 2])
  @assertEqual(12_iintegers, i1d)

  call ind_1d_to_nd(nd_offsets, i1d, ind)
  @assertEqual(2_iintegers, ind(1))
  @assertEqual(3_iintegers, ind(2))
  @assertEqual(2_iintegers, ind(3))

  do i1d = 1, 12
    call ind_1d_to_nd(nd_offsets, i1d, ind)
    @assertEqual(i1d, ind_nd_to_1d(nd_offsets, ind))
  enddo
end subroutine

@test(npes=[1])
subroutine test_ndarray_idx_c_style(this)
  class (MpiTestMethod), intent(inout) :: this
  integer(iintegers) :: arr_shape(3)
  integer(iintegers) :: nd_offsets(3)
  integer(iintegers) :: i1d, ind(3)

  arr_shape = [2,3,2]
  call ndarray_offsets(arr_shape, nd_offsets)

  @assertEqual(1_iintegers, nd_offsets(1))
  @assertEqual(2_iintegers, nd_offsets(2))
  @assertEqual(6_iintegers, nd_offsets(3))

  i1d = ind_nd_to_1d(nd_offsets, [integer(iintegers) :: 0, 0, 0], cstyle=.True.)
  @assertEqual(0_iintegers, i1d)

  call ind_1d_to_nd(nd_offsets, i1d, ind, cstyle=.True.)
  @assertEqual(0_iintegers, ind(1))
  @assertEqual(0_iintegers, ind(2))
  @assertEqual(0_iintegers, ind(3))

  i1d = ind_nd_to_1d(nd_offsets, [integer(iintegers) :: 1, 2, 1], cstyle=.True.)
  @assertEqual(11_iintegers, i1d)

  call ind_1d_to_nd(nd_offsets, i1d, ind, cstyle=.True.)
  @assertEqual(1_iintegers, ind(1))
  @assertEqual(2_iintegers, ind(2))
  @assertEqual(1_iintegers, ind(3))

  do i1d = 0, 11
    call ind_1d_to_nd(nd_offsets, i1d, ind, cstyle=.True.)
    @assertEqual(i1d, ind_nd_to_1d(nd_offsets, ind, cstyle=.True.))
  enddo
end subroutine

@test(npes=[1])
subroutine test_spherical2cartesian(this)
  class (MpiTestMethod), intent(inout) :: this
  real(ireals), parameter :: eps=max(sqrt(epsilon(eps)), 1e-5_ireals)
  real(ireals) :: sundir(3), sundir2(3), phi, theta
  integer(iintegers) :: iter
  integer(mpiint) :: ierr

  sundir = [real(ireals) :: 0, -tiny(sundir), -1]
  call cartesian_2_spherical(sundir, phi, theta, ierr)
  @assertEqual(0._ireals, phi, eps)
  @assertEqual(0._ireals, theta, eps)

  sundir = [real(ireals) :: 0, +tiny(sundir), -1]
  call cartesian_2_spherical(sundir, phi, theta, ierr)
  @assertEqual(180._ireals, abs(phi), eps)
  @assertEqual(0._ireals, theta, eps)

  sundir = [real(ireals) :: -tiny(sundir), 0, -1]
  call cartesian_2_spherical(sundir, phi, theta, ierr)
  @assertEqual(90._ireals, phi, eps)
  @assertEqual(0._ireals, theta, eps)

  sundir = [real(ireals) :: +tiny(sundir), 0, -1]
  call cartesian_2_spherical(sundir, phi, theta, ierr)
  @assertEqual(-90._ireals, phi, eps)
  @assertEqual(0._ireals, theta, eps)

  sundir = [real(ireals) :: 0, -1, -1]
  call cartesian_2_spherical(sundir, phi, theta, ierr)
  @assertEqual(0._ireals, phi, eps)
  @assertEqual(45._ireals, theta, eps)

  sundir = [real(ireals) :: -1, -1, -1]
  call cartesian_2_spherical(sundir, phi, theta, ierr)
  @assertEqual(45._ireals, phi, eps)
  @assertEqual(54.73561_ireals, theta, eps)

  sundir = [real(ireals) :: 1, 1, -1]
  call cartesian_2_spherical(sundir, phi, theta, ierr)
  @assertEqual(-135._ireals, phi, eps)
  @assertEqual(54.73561_ireals, theta, eps)

  ! roundtrips
  do iter = 1, 1000
    call random_number(sundir)
    call normalize_vec(sundir, ierr)
    @assertEqual(0_mpiint, ierr)
    call cartesian_2_spherical(sundir, phi, theta, ierr)
    @assertEqual(0_mpiint, ierr)
    sundir2 = spherical_2_cartesian(phi, theta)
    if(.not.all(approx(sundir, sundir2, eps))) print *,'iter', iter, 'sundir', sundir, 'phi/theta', phi, theta, 'sundir_rev', sundir2
    @assertEqual(sundir, sundir2, eps)
  enddo
end subroutine


@test(npes=[1])
subroutine test_expm1(this)
  class (MpiTestMethod), intent(inout) :: this
  real(ireals), parameter :: eps=sqrt(epsilon(eps))
  real(ireals) :: x, y

  x = 1e-20_ireals
  !y = (1.-exp(-x)) / x ! the trivial equation is not numerically stable around 0
  y = - expm1(-x) / x
  @assertEqual(1._ireals, y, eps)
end subroutine

end module
