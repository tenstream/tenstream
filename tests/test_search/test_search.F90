module test_search
  use iso_c_binding
  use m_data_parameters, only: ireals, ireal_dp, iintegers, mpiint, init_mpi_data_parameters
  use m_search, only: search_sorted_bisection, find_real_location_petsc, find_real_location_linear, &
                      find_real_location

  use pfunit_mod

  implicit none

contains

  @test(npes=[1])
  subroutine test_search_sorted_bisection(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: A(3) = [-10, 0, 2]

    @assertEqual(1.0, search_sorted_bisection(A, -20._ireals))
    @assertEqual(1.0, search_sorted_bisection(A, -10._ireals))
    @assertEqual(1.25, search_sorted_bisection(A, -7.5_ireals))
    @assertEqual(1.5, search_sorted_bisection(A, -5._ireals))
    @assertEqual(2.0, search_sorted_bisection(A, -0._ireals))
    @assertEqual(2.0, search_sorted_bisection(A, 0._ireals))
    @assertEqual(2.5, search_sorted_bisection(A, 1._ireals))
    @assertEqual(3.0, search_sorted_bisection(A, 2._ireals))
    @assertEqual(3.0, search_sorted_bisection(A, 3._ireals))
    @assertEqual(1.0, search_sorted_bisection([1._ireals], -1._ireals))
    @assertEqual(1.0, search_sorted_bisection([1._ireals], 1._ireals))
    @assertEqual(1.0, search_sorted_bisection([1._ireals], 2._ireals))
    @assertEqual(1.0, search_sorted_bisection([0._ireals, 1._ireals], -1._ireals))
    @assertEqual(2.0, search_sorted_bisection([0._ireals, 1._ireals], 1._ireals))
    @assertEqual(2.0, search_sorted_bisection([0._ireals, 1._ireals], 2._ireals))
  end subroutine

  @test(npes=[1])
  subroutine test_search_sorted_bisection_reversed(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: A(3) = [10, 0, -2]

    @assertEqual(1.0, search_sorted_bisection(A, 20._ireals))
    @assertEqual(1.0, search_sorted_bisection(A, 10._ireals))
    @assertEqual(1.25, search_sorted_bisection(A, 7.5_ireals))
    @assertEqual(1.5, search_sorted_bisection(A, 5._ireals))
    @assertEqual(2.0, search_sorted_bisection(A, 0._ireals))
    @assertEqual(2.0, search_sorted_bisection(A, -0._ireals))
    @assertEqual(2.5, search_sorted_bisection(A, -1._ireals))
    @assertEqual(3.0, search_sorted_bisection(A, -2._ireals))
    @assertEqual(3.0, search_sorted_bisection(A, -3._ireals))
  end subroutine

  @test(npes=[1])
  subroutine test_find_real_location_petsc(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: A(3) = [-10, 0, 2]

    @assertEqual(1.0, find_real_location_petsc(A, -20._ireals))
    @assertEqual(1.0, find_real_location_petsc(A, -10._ireals))
    @assertEqual(1.25, find_real_location_petsc(A, -7.5_ireals))
    @assertEqual(1.5, find_real_location_petsc(A, -5._ireals))
    @assertEqual(2.0, find_real_location_petsc(A, -0._ireals))
    @assertEqual(2.0, find_real_location_petsc(A, 0._ireals))
    @assertEqual(2.5, find_real_location_petsc(A, 1._ireals))
    @assertEqual(3.0, find_real_location_petsc(A, 2._ireals))
    @assertEqual(3.0, find_real_location_petsc(A, 3._ireals))
  end subroutine

  @test(npes=[1])
  subroutine test_find_real_location(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: A(3) = [-10, 0, 2]

    @assertEqual(1.0, find_real_location(A, -20._ireals))
    @assertEqual(1.0, find_real_location(A, -10._ireals))
    @assertEqual(1.25, find_real_location(A, -7.5_ireals))
    @assertEqual(1.5, find_real_location(A, -5._ireals))
    @assertEqual(2.0, find_real_location(A, -0._ireals))
    @assertEqual(2.0, find_real_location(A, 0._ireals))
    @assertEqual(2.5, find_real_location(A, 1._ireals))
    @assertEqual(3.0, find_real_location(A, 2._ireals))
    @assertEqual(3.0, find_real_location(A, 3._ireals))
  end subroutine

  @test(npes=[1])
  subroutine test_find_real_location_reversed(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: A(3) = [10, 0, -2]

    @assertEqual(1.0, find_real_location(A, 20._ireals))
    @assertEqual(1.0, find_real_location(A, 10._ireals))
    @assertEqual(1.25, find_real_location(A, 7.5_ireals))
    @assertEqual(1.5, find_real_location(A, 5._ireals))
    @assertEqual(2.0, find_real_location(A, 0._ireals))
    @assertEqual(2.0, find_real_location(A, -0._ireals))
    @assertEqual(2.5, find_real_location(A, -1._ireals))
    @assertEqual(3.0, find_real_location(A, -2._ireals))
    @assertEqual(3.0, find_real_location(A, -3._ireals))
  end subroutine

  @test(npes=[1])
  subroutine test_find_real_location_linear(this)
    class(MpiTestMethod), intent(inout) :: this
    real(ireals), parameter :: A(3) = [-10, 0, 2]

    @assertEqual(1.0, find_real_location_linear(A, -20._ireals))
    @assertEqual(1.0, find_real_location_linear(A, -10._ireals))
    @assertEqual(1.25, find_real_location_linear(A, -7.5_ireals))
    @assertEqual(1.5, find_real_location_linear(A, -5._ireals))
    @assertEqual(2.0, find_real_location_linear(A, -0._ireals))
    @assertEqual(2.0, find_real_location_linear(A, 0._ireals))
    @assertEqual(2.5, find_real_location_linear(A, 1._ireals))
    @assertEqual(3.0, find_real_location_linear(A, 2._ireals))
    @assertEqual(3.0, find_real_location_linear(A, 3._ireals))
  end subroutine

  @test(npes=[1])
  subroutine test_search_runtime(this)
    class(MpiTestMethod), intent(inout) :: this
    integer(iintegers), parameter :: Nsize = 200, Niter = 100000 * Nsize
    integer(iintegers) :: i, s, itest
    real(ireals), allocatable :: A(:), r(:)
    real(ireal_dp) :: tstart, tend
    real(ireals) :: sum_res(4), time(4)
    print *, 'Running performance test for search routines'
    print *, '                    find_real_location, '// &
      & 'find_real_location_linear, '// &
      & 'search_sorted_bisection, '// &
      & 'find_real_location_petsc'
    allocate (A(Niter), r(Niter))
    do i = 1, size(A)
      A(i) = real(i - 1, ireals)
    end do
    do i = 1, Niter
      call random_number(r(i))
      r(i) = r(i) * size(A)
    end do

    do s = 5, Nsize, 5
      associate (arr => A(1:s))

        sum_res = 0

        itest = 1
        call cpu_time(tstart)
        do i = 1, Niter / s
          sum_res(itest) = sum_res(itest) + find_real_location(arr, r(i))
        end do
        call cpu_time(tend)
        time(itest) = real(tend - tstart, ireals)

        itest = 2
        call cpu_time(tstart)
        do i = 1, Niter / s
          sum_res(itest) = sum_res(itest) + find_real_location_linear(arr, r(i))
        end do
        call cpu_time(tend)
        time(itest) = real(tend - tstart, ireals)

        itest = 3
        call cpu_time(tstart)
        do i = 1, Niter / s
          sum_res(itest) = sum_res(itest) + search_sorted_bisection(arr, r(i))
        end do
        call cpu_time(tend)
        time(itest) = real(tend - tstart, ireals)

        itest = 4
        call cpu_time(tstart)
        do i = 1, Niter / s
          sum_res(itest) = sum_res(itest) + find_real_location_petsc(arr, r(i))
        end do
        call cpu_time(tend)
        time(itest) = real(tend - tstart, ireals)

        do itest = 2, size(sum_res)
          @assertEqual(sum_res(1), sum_res(itest))
        end do
        print *, 'arr size', s, 'time', time, ':', time / minval(time), &
          & ': auto_select algorithm', time(1), '( ratio=', time(1) / minval(time), ')'
        ! make sure that the generic version of search algorithm is selected well,
        ! i.e. that we are not off by a factor of 3
        @assertTrue(time(1) .lt. minval(time) * 2, 'the auto selected search algorithm was the wrong one')
      end associate
    end do
  end subroutine

end module
