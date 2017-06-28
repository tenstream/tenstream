@test(npes =[1])
subroutine test_interp(this)

    use m_data_parameters, only: ireals, default_str_len
    use m_tenstream_interpolation, only: interp_2d

    use pfunit_mod

    implicit none

    class (MpiTestMethod), intent(inout) :: this

    real(ireals),parameter :: tol=1e-6
    real(ireals) :: db_2d(1,2,2)
    real(ireals) :: res(1)

    integer :: i

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
    do i=1,2
      print *,'2D LUT Data:',i , db_2d(1,i,:)
    enddo

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
