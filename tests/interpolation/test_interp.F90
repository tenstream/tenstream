module test_interp
    use m_data_parameters, only: iintegers, irealLUT, default_str_len
    use m_helper_functions, only: ndarray_offsets, meanval
    use m_tenstream_interpolation, only: interp_1d, interp_2d, &
      interp_vec_simplex_nd

    use pfunit_mod

    implicit none

    real(irealLUT),parameter :: one=1, zero=0, tol=sqrt(epsilon(one))*10

contains
  @test(npes =[1])
  subroutine test_interp_vec_simplex_1d_2d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(irealLUT) :: db_2d(1,2,2), db(1,4)
    real(irealLUT) :: res(1)
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
    call ndarray_offsets(shape(db_2d(1,:,:), iintegers), db_offsets)

    ! 0-dim
    do i=1,2
      do j=1,2
        call interp_vec_simplex_nd([one*i, one*j], db, db_offsets, res)
        @assertEqual(db_2d(1,i,j) , res(1),  tol, msg)
      enddo
    enddo

    ! 1-dim
    call interp_vec_simplex_nd([1.5_irealLUT, 1._irealLUT], db, db_offsets, res)
    @assertEqual(.5_irealLUT , res(1),  tol, msg)
    call interp_vec_simplex_nd([1.75_irealLUT, 1._irealLUT], db, db_offsets, res)
    @assertEqual(.75_irealLUT , res(1),  tol, msg)

    call interp_vec_simplex_nd([1.5_irealLUT, 2._irealLUT], db, db_offsets, res)
    @assertEqual(3._irealLUT , res(1),  tol, msg)
    call interp_vec_simplex_nd([1.75_irealLUT, 2._irealLUT], db, db_offsets, res)
    @assertEqual(3.5_irealLUT , res(1),  tol, msg)

    call interp_vec_simplex_nd([1._irealLUT, 1.5_irealLUT], db, db_offsets, res)
    @assertEqual(1._irealLUT , res(1),  tol, msg)
    call interp_vec_simplex_nd([1._irealLUT, 1.75_irealLUT], db, db_offsets, res)
    @assertEqual(1.5_irealLUT , res(1),  tol, msg)

    call interp_vec_simplex_nd([2._irealLUT, 1.5_irealLUT], db, db_offsets, res)
    @assertEqual(2.5_irealLUT , res(1),  tol, msg)
    call interp_vec_simplex_nd([2._irealLUT, 1.75_irealLUT], db, db_offsets, res)
    @assertEqual(3.25_irealLUT , res(1),  tol, msg)

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

    call interp_vec_simplex_nd([1.5_irealLUT, 1.5_irealLUT], db, db_offsets, res)
    @assertEqual(2._irealLUT , res(1),  tol, msg)

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

    real(irealLUT) :: db_3d(1,2,2,2), db(1,8)
    real(irealLUT) :: res(1)
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
    call ndarray_offsets(shape(db_3d(1,:,:,:), iintegers), db_offsets)

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
    call interp_vec_simplex_nd([one, one, one]+.5_irealLUT, db, db_offsets, res)
    @assertEqual(meanval(db), res(1),  tol, msg)

  end subroutine

  @test(npes =[1])
  subroutine test_interp_vec_simplex_4d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(irealLUT) :: db_4d(1,2,2,2,2), db(1,16)
    real(irealLUT) :: res(1)
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
    call ndarray_offsets(shape(db_4d(1,:,:,:,:), iintegers), db_offsets)

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
    call interp_vec_simplex_nd([one, one, one, one]+.5_irealLUT, db, db_offsets, res)
    @assertEqual(meanval(db), res(1),  tol, msg)

  end subroutine

  @test(npes =[1])
  subroutine test_interp2d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(irealLUT) :: db_2d(1,2,2)
    real(irealLUT) :: res(1)

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

    call interp_2d([1._irealLUT, 1.5_irealLUT], db_2d, res)
    @assertEqual(1_irealLUT , res(1),  tol, msg)

    call interp_2d([2._irealLUT, 1.5_irealLUT], db_2d, res)
    @assertEqual(2.5_irealLUT , res(1),  tol, msg)

    call interp_2d([1.5_irealLUT, 1._irealLUT], db_2d, res)
    @assertEqual(0.5_irealLUT , res(1),  tol, msg)

    call interp_2d([1.5_irealLUT, 2._irealLUT], db_2d, res)
    @assertEqual(3_irealLUT , res(1),  tol, msg)

    call interp_2d([1.5_irealLUT, 1.5_irealLUT], db_2d, res)
    @assertEqual(1.75_irealLUT , res(1),  tol, msg)
  end subroutine

  @test(npes =[1])
  subroutine test_interp1d(this)
    class (MpiTestMethod), intent(inout) :: this

    real(irealLUT) :: db_1d(1,5)
    real(irealLUT) :: res

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

    res = interp_1d(1.5_irealLUT, db_1d(1,:))
    @assertEqual(zero , res,  tol, msg)

    res = interp_1d(2.4_irealLUT, db_1d(1,:))
    @assertEqual(one , res,  tol, msg)
    res = interp_1d(2.5_irealLUT, db_1d(1,:))
    @assertEqual(one , res,  tol, msg)
    res = interp_1d(2.6_irealLUT, db_1d(1,:))
    @assertEqual(one , res,  tol, msg)

    res = interp_1d(3.5_irealLUT, db_1d(1,:))
    @assertEqual(5.5_irealLUT, res,  tol, msg)

    res = interp_1d(3.7_irealLUT, db_1d(1,:))
    @assertEqual(7.3_irealLUT, res,  tol, msg)

    res = interp_1d(4.1_irealLUT, db_1d(1,:))
    @assertEqual(8.9_irealLUT, res,  tol, msg)

    res = interp_1d(4.9_irealLUT, db_1d(1,:))
    @assertEqual(0.1_irealLUT, res,  tol, msg)

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

    real(irealLUT) :: db_1d(2,5)
    real(irealLUT) :: res(2)
    integer(iintegers) :: offsets(1)

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

    call ndarray_offsets(shape(db_1d(1,:), kind=iintegers), offsets)

    call interp_vec_simplex_nd([1.5_irealLUT], db_1d, offsets, res)
    @assertEqual(zero , res,  tol, msg)

    call interp_vec_simplex_nd([2.4_irealLUT], db_1d, offsets, res)
    @assertEqual(one , res,  tol, msg)

    call interp_vec_simplex_nd([2.5_irealLUT], db_1d, offsets, res)
    @assertEqual(one , res,  tol, msg)

    call interp_vec_simplex_nd([2.6_irealLUT], db_1d, offsets, res)
    @assertEqual(one , res,  tol, msg)

    call interp_vec_simplex_nd([3.5_irealLUT], db_1d, offsets, res)
    @assertEqual(5.5_irealLUT, res,  tol, msg)

    call interp_vec_simplex_nd([3.7_irealLUT], db_1d, offsets, res)
    @assertEqual(7.3_irealLUT, res,  tol, msg)

    call interp_vec_simplex_nd([4.1_irealLUT], db_1d, offsets, res)
    @assertEqual(8.9_irealLUT, res,  tol, msg)

    call interp_vec_simplex_nd([4.9_irealLUT], db_1d, offsets, res)
    @assertEqual(0.1_irealLUT, res,  tol, msg)

    do i=1,size(db_1d, dim=2)
    call interp_vec_simplex_nd([one*i], db_1d, offsets, res)
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
    do i=2,size(db_1d, dim=2)
    call interp_vec_simplex_nd([one*i-epsilon(one)], db_1d, offsets, res)
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
    do i=1,size(db_1d, dim=2)-1
    call interp_vec_simplex_nd([one*i+epsilon(one)], db_1d, offsets, res)
      @assertEqual(db_1d(1,i) , res,  tol, msg)
    enddo
  end subroutine


  @test(npes =[1])
  subroutine test_interp_vec_6d(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers), parameter :: Nv=4
    real(irealLUT) :: db_6d(2,Nv,Nv,Nv,Nv,Nv,Nv), db(2,Nv**6)
    real(irealLUT) :: res(2)
    integer(iintegers) :: db_offsets(6)

    integer :: i, j, k, l, m, n

    character(default_str_len) :: msg
    write(msg,*) "6D linear interpolation not as expected: "

    do n=1,Nv
    do m=1,Nv
    do l=1,Nv
    do k=1,Nv
    do j=1,Nv
    do i=1,Nv
      db_6d(1, i,j,k,l,m,n) = x(1.*i,1.*j,1.*k,1.*l,1.*m,1.*n)
      db_6d(2, i,j,k,l,m,n) = x(2.*i,2.*j,2.*k,2.*l,2.*m,2.*n)
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo

    db(1,:) = pack(db_6d(1,:,:,:,:,:,:), .True.)
    db(2,:) = pack(db_6d(2,:,:,:,:,:,:), .True.)

    call ndarray_offsets(shape(db_6d(1,:,:,:,:,:,:), iintegers), db_offsets)

    do n=1,Nv-1
    do m=1,Nv-1
    do l=1,Nv-1
    do k=1,Nv-1
    do j=1,Nv-1
    do i=1,Nv-1
      ! 0-dim, just check corners, this should actually not even end up in high dimensional interp but just go to 0D
      call interp_vec_simplex_nd(real([i, j, k, l, m, n], irealLUT), db, db_offsets, res)
      !print *,i, j, k, l, m, n, 'res', res, ':', x(1.*i,1.*j,1.*k,1.*l,1.*m,1.*n), x(2.*i,2.*j,2.*k,2.*l,2.*m,2.*n)
      @assertEqual(x(1.*i,1.*j,1.*k,1.*l,1.*m,1.*n), res(1),  msg)
      @assertEqual(x(2.*i,2.*j,2.*k,2.*l,2.*m,2.*n), res(2),  msg)

      ! Then check in between
      call interp_vec_simplex_nd(real([i, j, k, l, m, n], irealLUT)+.5, db, db_offsets, res)
      @assertEqual(x(i+.5,j+.5,k+.5,l+.5,m+.5,n+.5) , res(1),  msg)
      @assertEqual(x(2*(i+.5),2*(j+.5),2*(k+.5),2*(l+.5),2*(m+.5),2*(n+.5)) , res(2),  msg)
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo

    contains
      function x(i,j,k,l,m,n)
        real(irealLUT) :: i,j,k,l,m,n
        real(irealLUT) :: x
        x = &
        2**0 * i + &
        2**1 * j + &
        2**3 * k + &
        2**4 * l + &
        2**5 * m + &
        2**6 * n
      end function
  end subroutine

end module
