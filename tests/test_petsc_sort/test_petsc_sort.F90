module test_petsc_sort

#include "petsc/finclude/petsc.h"
  use petsc
  use m_data_parameters, only : iintegers, ireals, mpiint, zero, one

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine teardown


  @Test(npes=[1])
  subroutine petsc_sort_int_arrays(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers),parameter :: N=3, a=-1, b=2, c=5
    integer(iintegers) :: x(N), x2(N)
    integer(iintegers), parameter :: t(N)=[a,b,c]
    integer(mpiint) :: ierr

    x = [b,c,a]
    call PetscSortInt(N, x, ierr)
    @assertEqual(t,x)

    x = [c,b,a]
    x2= [a,b,c]
    call PetscSortIntWithArray(N, x, x2, ierr)
    @assertEqual(t,x)
    @assertEqual([c,b,a],x2)
  end subroutine

  @Test(npes=[1])
  subroutine petsc_sort_int_array_with_type(this)
    class (MpiTestMethod), intent(inout) :: this
    integer(iintegers),parameter :: N=3, a=-1, b=2, c=5
    integer(mpiint) :: ierr
    integer(iintegers), parameter :: t(N)=[a,b,c]
    type t_test
      integer(iintegers) :: i,j
    end type
    type(t_test) :: x(N)

    x(:)%i = [c,b,a]
    x(:)%j = [a,b,c]
    call PetscSortIntWithArray(N, x(:)%i, x(:)%j, ierr)
    @assertEqual(t,x(:)%i)
    @assertEqual([c,b,a],x(:)%j)
  end subroutine

end module
