module test_interp
    use m_data_parameters, only: iintegers, ireals, default_str_len, one, zero
    use m_helper_functions, only: ndarray_offsets, mean
    use m_tenstream_interpolation, only: interp_1d, interp_vec_1d, interp_2d, &
      interp_vec_simplex_nd

    use pfunit_mod

    implicit none

    real(ireals),parameter :: tol=sqrt(epsilon(one))*10
contains
  @test(npes =[1])
  subroutine test_interp_vec_simplex_1d_2d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals) :: db_2d(1,2,2), db(1,4)
    real(ireals) :: res(1)
    integer(iintegers) :: db_shape(3)
    integer(iintegers) :: db_offsets(2)

    integer :: i, j

    character(default_str_len) :: msg
    write(msg,*) "2D linear interpolation not as expected: "

    ! 2D LUT in this example @ coordinates: (1,1), (2,1), (1,2), (2,2)
    !
    ! 2 --- 4
    ! |     |
    ! |     |
    ! 0 --- 1

    db_2d(1,1,1) = 0
    db_2d(1,2,1) = 1
    db_2d(1,1,2) = 2
    db_2d(1,2,2) = 4

    db(1,:) = pack(db_2d, .True.)

    db_shape = shape(db_2d, iintegers)
    db_offsets = ndarray_offsets(shape(db_2d(1,:,:), iintegers))

    ! 0-dim
    do i=1,2
      do j=1,2
        call interp_vec_simplex_nd([one*i, one*j], db, db_offsets, res)
        @assertEqual(db_2d(1,i,j) , res(1),  tol, msg)
      enddo
    enddo

    ! 1-dim
    call interp_vec_simplex_nd([1.5_ireals, 1._ireals], db, db_offsets, res)
    @assertEqual(.5_ireals , res(1),  tol, msg)
    call interp_vec_simplex_nd([1.75_ireals, 1._ireals], db, db_offsets, res)
    @assertEqual(.75_ireals , res(1),  tol, msg)

    call interp_vec_simplex_nd([1.5_ireals, 2._ireals], db, db_offsets, res)
    @assertEqual(3._ireals , res(1),  tol, msg)
    call interp_vec_simplex_nd([1.75_ireals, 2._ireals], db, db_offsets, res)
    @assertEqual(3.5_ireals , res(1),  tol, msg)

    call interp_vec_simplex_nd([1._ireals, 1.5_ireals], db, db_offsets, res)
    @assertEqual(1._ireals , res(1),  tol, msg)
    call interp_vec_simplex_nd([1._ireals, 1.75_ireals], db, db_offsets, res)
    @assertEqual(1.5_ireals , res(1),  tol, msg)

    call interp_vec_simplex_nd([2._ireals, 1.5_ireals], db, db_offsets, res)
    @assertEqual(2.5_ireals , res(1),  tol, msg)
    call interp_vec_simplex_nd([2._ireals, 1.75_ireals], db, db_offsets, res)
    @assertEqual(3.25_ireals , res(1),  tol, msg)

    ! 2-dim Checks:

    ! 2 --- 4
    ! |     |
    ! |     |
    ! 0 --- 2

    db_2d(1,1,1) = 0
    db_2d(1,2,1) = 2
    db_2d(1,1,2) = 2
    db_2d(1,2,2) = 4
    db(1,:) = pack(db_2d, .True.)

    call interp_vec_simplex_nd([1.5_ireals, 1.5_ireals], db, db_offsets, res)
    @assertEqual(2._ireals , res(1),  tol, msg)

    ! 2D LUT in this example @ coordinates: (1,1), (2,1), (1,2), (2,2)
    !
    ! 2 --- 2
    ! |     |
    ! |     |
    ! 0 --- 0

    db_2d(1,1,1) = 0
    db_2d(1,2,1) = 0
    db_2d(1,1,2) = 2
    db_2d(1,2,2) = 2
    db(1,:) = pack(db_2d, .True.)

    do j = 0 , 9
      do i = 0 , 9
        call interp_vec_simplex_nd([one+i*.1, one+j*.1], db, db_offsets, res)
        @assertEqual(j*.1*2 , res(1),  tol, msg)
      enddo
    enddo
  end subroutine

  @test(npes =[1])
  subroutine test_interp_vec_simplex_3d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals) :: db_3d(1,2,2,2), db(1,8)
    real(ireals) :: res(1)
    integer(iintegers) :: db_shape(4)
    integer(iintegers) :: db_offsets(3)

    integer :: i, j, k

    character(default_str_len) :: msg
    write(msg,*) "3D linear interpolation not as expected: "

    ! 2D LUT in this example @ coordinates: (1,1), (2,1), (1,2), (2,2)
    !
    ! 2 --- 4
    ! |     |
    ! |     |
    ! 0 --- 2

    db_3d(1,1,1,1) = 0
    db_3d(1,2,1,1) = 2
    db_3d(1,1,2,1) = 2
    db_3d(1,2,2,1) = 4

    db_3d(1,1,1,2) = 10
    db_3d(1,2,1,2) = 12
    db_3d(1,1,2,2) = 12
    db_3d(1,2,2,2) = 14

    db(1,:) = pack(db_3d, .True.)

    db_shape = shape(db_3d, iintegers)
    db_offsets = ndarray_offsets(shape(db_3d(1,:,:,:), iintegers))

    ! 0-dim, just check corners again, this should actually not even end up in 3dim interp but just go to 0D
    do i=1,2
      do j=1,2
        do k=1,2
          call interp_vec_simplex_nd([one*i, one*j, one*k], db, db_offsets, res)
          @assertEqual(db_3d(1,i,j,k) , res(1),  tol, msg)
        enddo
      enddo
    enddo

    ! Center point
    call interp_vec_simplex_nd([one, one, one]+.5_ireals, db, db_offsets, res)
    @assertEqual(mean(db), res(1),  tol, msg)

  end subroutine

  @test(npes =[1])
  subroutine test_interp_vec_simplex_4d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals) :: db_4d(1,2,2,2,2), db(1,16)
    real(ireals) :: res(1)
    integer(iintegers) :: db_shape(5)
    integer(iintegers) :: db_offsets(4)

    integer :: i, j, k, l

    character(default_str_len) :: msg
    write(msg,*) "3D linear interpolation not as expected: "

    ! 2D LUT in this example @ coordinates: (1,1), (2,1), (1,2), (2,2)
    !
    ! 2 --- 4
    ! |     |
    ! |     |
    ! 0 --- 2

    db_4d(1,1,1,1,1) = 0
    db_4d(1,2,1,1,1) = 2
    db_4d(1,1,2,1,1) = 2
    db_4d(1,2,2,1,1) = 4

    db_4d(1,1,1,2,1) = 10
    db_4d(1,2,1,2,1) = 12
    db_4d(1,1,2,2,1) = 12
    db_4d(1,2,2,2,1) = 14

    db_4d(1,1,1,1,2) = 0
    db_4d(1,2,1,1,2) = 2
    db_4d(1,1,2,1,2) = 2
    db_4d(1,2,2,1,2) = 4

    db_4d(1,1,1,2,2) = 10
    db_4d(1,2,1,2,2) = 12
    db_4d(1,1,2,2,2) = 12
    db_4d(1,2,2,2,2) = 14

    db(1,:) = pack(db_4d, .True.)

    db_shape = shape(db_4d, iintegers)
    db_offsets = ndarray_offsets(shape(db_4d(1,:,:,:,:), iintegers))

    ! 0-dim, just check corners again, this should actually not even end up in 3dim interp but just go to 0D
    do i=1,2
      do j=1,2
        do k=1,2
          do l=1,2
            call interp_vec_simplex_nd([one*i, one*j, one*k, one*l], db, db_offsets, res)
            @assertEqual(db_4d(1,i,j,k,l) , res(1),  tol, msg)
          enddo
        enddo
      enddo
    enddo

    ! Center point
    call interp_vec_simplex_nd([one, one, one, one]+.5_ireals, db, db_offsets, res)
    @assertEqual(mean(db), res(1),  tol, msg)

  end subroutine

  @test(npes =[1])
  subroutine test_interp2d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals) :: db_2d(1,2,2)
    real(ireals) :: res(1)

    character(default_str_len) :: msg
    write(msg,*) "2D linear interpolation not as expected: "

    ! 2D LUT in this example @ coordinates: (1,1), (2,1), (1,2), (2,2)
    !
    ! 2 --- 4
    ! |     |
    ! |     |
    ! 0 --- 1
    !
    db_2d = reshape([ 0, 1, 2, 4 ], [1, 2, 2])
    !do i=1,2
    !  print *,'2D LUT Data:',i , db_2d(1,i,:)
    !enddo

    call interp_2d([1._ireals, 1.5_ireals], db_2d, res)
    @assertEqual(1_ireals , res(1),  tol, msg)

    call interp_2d([2._ireals, 1.5_ireals], db_2d, res)
    @assertEqual(2.5_ireals , res(1),  tol, msg)

    call interp_2d([1.5_ireals, 1._ireals], db_2d, res)
    @assertEqual(0.5_ireals , res(1),  tol, msg)

    call interp_2d([1.5_ireals, 2._ireals], db_2d, res)
    @assertEqual(3_ireals , res(1),  tol, msg)

    call interp_2d([1.5_ireals, 1.5_ireals], db_2d, res)
    @assertEqual(1.75_ireals , res(1),  tol, msg)
  end subroutine

  @test(npes =[1])
  subroutine test_interp1d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals) :: db_1d(1,5)
    real(ireals) :: res

    integer :: i

    character(default_str_len) :: msg
    write(msg,*) "1D linear interpolation not as expected: "

    ! 1D LUT in this example @ coordinates:
    !
    !    ^
    !    |
    ! 10 |              x
    !    |
    !  1 |      x   x
    !    |
    ! -1 |  x               x
    !    --------------------->
    !       1   2   3   4   5

    db_1d(1,:) = [ -1, 1, 1, 10, -1 ]

    res = interp_1d(1.5_ireals, db_1d(1,:))
    @assertEqual(zero , res,  tol, msg)

    res = interp_1d(2.4_ireals, db_1d(1,:))
    @assertEqual(one , res,  tol, msg)
    res = interp_1d(2.5_ireals, db_1d(1,:))
    @assertEqual(one , res,  tol, msg)
    res = interp_1d(2.6_ireals, db_1d(1,:))
    @assertEqual(one , res,  tol, msg)

    res = interp_1d(3.5_ireals, db_1d(1,:))
    @assertEqual(5.5_ireals, res,  tol, msg)

    res = interp_1d(3.7_ireals, db_1d(1,:))
    @assertEqual(7.3_ireals, res,  tol, msg)

    res = interp_1d(4.1_ireals, db_1d(1,:))
    @assertEqual(8.9_ireals, res,  tol, msg)

    res = interp_1d(4.9_ireals, db_1d(1,:))
    @assertEqual(0.1_ireals, res,  tol, msg)

    do i=1,size(db_1d, dim=2)
      res = interp_1d(one*i, db_1d(1,:))
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
    do i=2,size(db_1d, dim=2)
      res = interp_1d(one*i-epsilon(one), db_1d(1,:))
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
    do i=1,size(db_1d, dim=2)-1
      res = interp_1d(one*i+epsilon(one), db_1d(1,:))
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
  end subroutine
  @test(npes =[1])
  subroutine test_interp1d_vec(this)
    class (MpiTestMethod), intent(inout) :: this

    real(ireals) :: db_1d(2,5)
    real(ireals) :: res(2)

    integer :: i

    character(default_str_len) :: msg
    write(msg,*) "1D linear interpolation not as expected: "

    ! 1D LUT in this example @ coordinates:
    !
    !    ^
    !    |
    ! 10 |              x
    !    |
    !  1 |      x   x
    !    |
    ! -1 |  x               x
    !    --------------------->
    !       1   2   3   4   5

    db_1d(1,:) = [ -1, 1, 1, 10, -1 ]
    db_1d(2,:) = db_1d(1,:)

    res = interp_vec_1d(1.5_ireals, db_1d)
    @assertEqual(zero , res,  tol, msg)

    res = interp_vec_1d(2.4_ireals, db_1d)
    @assertEqual(one , res,  tol, msg)
    res = interp_vec_1d(2.5_ireals, db_1d)
    @assertEqual(one , res,  tol, msg)
    res = interp_vec_1d(2.6_ireals, db_1d)
    @assertEqual(one , res,  tol, msg)

    res = interp_vec_1d(3.5_ireals, db_1d)
    @assertEqual(5.5_ireals, res,  tol, msg)

    res = interp_vec_1d(3.7_ireals, db_1d)
    @assertEqual(7.3_ireals, res,  tol, msg)

    res = interp_vec_1d(4.1_ireals, db_1d)
    @assertEqual(8.9_ireals, res,  tol, msg)

    res = interp_vec_1d(4.9_ireals, db_1d)
    @assertEqual(0.1_ireals, res,  tol, msg)

    do i=1,size(db_1d, dim=2)
      res = interp_vec_1d(one*i, db_1d)
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
    do i=2,size(db_1d, dim=2)
      res = interp_vec_1d(one*i-epsilon(one), db_1d)
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
    do i=1,size(db_1d, dim=2)-1
      res = interp_vec_1d(one*i+epsilon(one), db_1d)
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
  end subroutine
end module
