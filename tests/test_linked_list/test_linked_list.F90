module test_linked_list
  use iso_c_binding
  use m_data_parameters, only: iintegers, mpiint, ireal_dp
  use m_linked_list_iintegers, only: t_list_iintegers, t_node

  use pfunit_mod

  implicit none

contains

  @test(npes=[1])
  subroutine test_list_iintegers(this)
    class(MpiTestMethod), intent(inout) :: this

    type(t_list_iintegers) :: list
    integer(iintegers) :: i
    integer(mpiint) :: ierr

    call list%add(1_iintegers)
    @assertEqual(1_iintegers, list%get_first())
    @assertEqual(1_iintegers, list%get_last())

    call list%add(2_iintegers)
    @assertEqual(1_iintegers, list%get_first())
    @assertEqual(2_iintegers, list%get_last())

    call list%for_each(square)
    @assertEqual(1_iintegers, list%get_first())
    @assertEqual(4_iintegers, list%get_last())

    @assertEqual(2_iintegers, list%len())
    call list%add(3_iintegers)
    @assertEqual(3_iintegers, list%len())

    call list%get_nth(2_iintegers, i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(4_iintegers, i)

    call list%get_nth(-2_iintegers, i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(4_iintegers, i)

    call list%get_nth(3_iintegers, i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(3_iintegers, i)

    call list%get_nth(-1_iintegers, i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(3_iintegers, i)

    call list%pop(i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(3_iintegers, i)
    call list%pop(i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(4_iintegers, i)
    call list%pop(i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(1_iintegers, i)

    call list%finalize()
    @assertEqual(0_iintegers, list%len())

    call list%add(1_iintegers)
    call list%add(2_iintegers)
    call list%add(3_iintegers)

    call list%pop_nth(5_iintegers, i, ierr)
    @assertLessThan(0_iintegers, ierr)

    call list%view()

    call list%pop_nth(2_iintegers, i, ierr)
    call list%view()
    @assertEqual(0_iintegers, ierr)
    @assertEqual(2_iintegers, list%len())

    @assertEqual(1_iintegers, list%get_first())
    @assertEqual(3_iintegers, list%get_last())

    call list%pop_nth(2_iintegers, i, ierr)
    call list%view()
    @assertEqual(0_iintegers, ierr)
    @assertEqual(1_iintegers, list%len())

    @assertEqual(1_iintegers, list%get_first())
    @assertEqual(1_iintegers, list%get_last())

    call list%pop_nth(1_iintegers, i, ierr)
    @assertEqual(0_iintegers, ierr)
    @assertEqual(0_iintegers, list%len())

    call list%pop_nth(1_iintegers, i, ierr)
    @assertLessThan(0_iintegers, ierr)

    call list%add(1_iintegers)
    call list%add(2_iintegers)
    call list%add(3_iintegers)
    call list%add(4_iintegers)
    call list%add(5_iintegers)

    call list%view()

    call list%pop_nth(1_iintegers, i, ierr)
    call list%pop_nth(4_iintegers, i, ierr)

    call list%view()
    call list%pop_nth(1_iintegers, i, ierr)
    call list%pop_nth(1_iintegers, i, ierr)
    call list%pop_nth(1_iintegers, i, ierr)
    call list%view()

    call list%add(1_iintegers)
    call list%add(2_iintegers)
    call list%add(3_iintegers)
    call list%add(4_iintegers)
    call list%add(5_iintegers)
    call list%view()
    call list%for_each(delete_uneven)
    call list%view()

  contains
    subroutine square(idx, node, val)
      integer(iintegers), intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(iintegers), intent(inout) :: val
      print *, 'squaring idx', idx, val
      val = val**2
    end subroutine

    recursive subroutine delete_uneven(idx, node, val)
      integer(iintegers), intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(iintegers), intent(inout) :: val
      if (modulo(val, 2_iintegers) .ne. 0) then
        print *, 'deleting', idx, val
        call list%del_node(node, ierr)
      end if
    end subroutine
  end subroutine

  @test(npes=[1])
  subroutine test_list_performance(this)
    class(MpiTestMethod), intent(inout) :: this

    type(t_list_iintegers) :: list
    integer(iintegers), parameter :: N = 100000
    integer(iintegers) :: i, s
    integer(mpiint) :: ierr

    real(ireal_dp) :: tstart, tend

    call cpu_time(tstart)
    do i = 1, N
      call list%add(i)
    end do
    call cpu_time(tend)
    print *, 'Time for add:', tend - tstart, 's'

    call cpu_time(tstart)
    s = 0
    do i = 1, N
      s = s + list%get_first()
    end do
    call cpu_time(tend)
    print *, 'Time for get_first():', tend - tstart, 's'

    call cpu_time(tstart)
    s = 0
    do i = 1, N
      s = s + list%get_last()
    end do
    call cpu_time(tend)
    print *, 'Time for get_last():', tend - tstart, 's'

    call cpu_time(tstart)
    call list%get_nth(N - 1_iintegers, s, ierr)
    call cpu_time(tend)
    print *, 'Time for get_nth(N-1):', tend - tstart, 's'

    call cpu_time(tstart)
    call list%get_nth(-1_iintegers, s, ierr)
    call cpu_time(tend)
    print *, 'Time for get_nth(-1):', tend - tstart, 's'

    call cpu_time(tstart)
    call list%for_each(delete_uneven)
    call cpu_time(tend)
    print *, 'Time for delete_uneven():', tend - tstart, 's'

  contains

    recursive subroutine delete_uneven(idx, node, val)
      integer(iintegers), intent(in) :: idx
      type(t_node), pointer, intent(inout) :: node
      integer(iintegers), intent(inout) :: val
      if (modulo(val, 2) .ne. 0) then
        call list%del_node(node, ierr)
      end if
    end subroutine
  end subroutine

end module
