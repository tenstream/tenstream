module m_search
  use iso_fortran_env, only: REAL32, REAL64
  use m_data_parameters,only : iintegers, mpiint, ireals
  use m_helper_functions, only: CHKERR

  implicit none

  interface find_real_location
    module procedure find_real_location_r32, find_real_location_r64
  end interface
  interface find_real_location_linear
    module procedure find_real_location_linear_r32, find_real_location_linear_r64
  end interface
  interface search_sorted_bisection
    module procedure search_sorted_bisection_r32, search_sorted_bisection_r64
  end interface

  integer(iintegers), parameter :: LINEAR_SEARCH_LIMIT=40

  contains

    ! return index+residula i where val is between arr(i) and arr(i+1)
    ! chooses different implementations depending on input array size
    pure function find_real_location_r32(arr, val) result(loc)
      real(REAL32),intent(in) :: arr(:)
      real(REAL32),intent(in) :: val
      real(REAL32)            :: loc
      loc = search_sorted_bisection(arr,val)
      ! Seems that the Petsc Func is always the slowest implementation but unnoticeable after N>100
      !select case (size(arr))
      !case(:LINEAR_SEARCH_LIMIT)
      !  loc = find_real_location_linear(arr, val)
      !case default
      !  loc = search_sorted_bisection(arr,val)
      !end select
    end function
    pure function find_real_location_r64(arr, val) result(loc)
      real(REAL64),intent(in) :: arr(:)
      real(REAL64),intent(in) :: val
      real(REAL64)            :: loc
      loc = search_sorted_bisection(arr,val)
      !select case (size(arr))
      !case(:LINEAR_SEARCH_LIMIT)
      !  loc = find_real_location_linear(arr, val)
      !case default
      !  loc = search_sorted_bisection(arr,val)
      !end select
    end function

    ! return index+residula i where val is between arr(i) and arr(i+1)
    pure function find_real_location_linear_r32(arr, val) result(loc)
      real(REAL32),intent(in) :: arr(:)
      real(REAL32),intent(in) :: val
      real(REAL32)            :: loc
      real(REAL32)            :: residual
      integer(iintegers)      :: i, N

      N = size(arr, kind=iintegers)

      if(arr(1).ge.val) then
        loc = 1
        return
      endif
      do i=2,N
        if(arr(i).ge.val) then
          residual = (val - arr(i-1)) / (arr(i) - arr(i-1))
          loc = real(i-1, real32) + residual
          return
        endif
      enddo
      loc = real(N, real32)
    end function
    pure function find_real_location_linear_r64(arr, val) result(loc)
      real(REAL64),intent(in) :: arr(:)
      real(REAL64),intent(in) :: val
      real(REAL64)            :: loc
      real(REAL64)            :: residual
      integer(iintegers)      :: i, N

      N = size(arr, kind=iintegers)

      if(arr(1).ge.val) then
        loc = 1
        return
      endif
      do i=2,N
        if(arr(i).ge.val) then
          residual = (val - arr(i-1)) / (arr(i) - arr(i-1))
          loc = real(i-1, real64) + residual
          return
        endif
      enddo
      loc = real(N, real64)
    end function

    ! return index+residula i where val is between arr(i) and arr(i+1)
    function find_real_location_petsc(arr, val) result(loc)
      real(ireals),intent(in) :: arr(:)
      real(ireals),intent(in) :: val
      real(ireals)            :: loc
      real(ireals)            :: residual
      integer(iintegers)      :: lb, N
      real(ireals), parameter :: eps=epsilon(eps)
      integer(mpiint)         :: ierr

      N = size(arr, kind=iintegers)
      call PetscFindReal(val, N, arr, eps, lb, ierr); call CHKERR(ierr)
      if(lb.eq.-1_iintegers) then
        loc = 1
        return
      else if(lb.eq.-N-1) then
        loc = real(N, ireals)
        return
      else if(lb.lt.-1._iintegers) then ! we have not hit the real value exactly, lets add the offset
        lb = -lb - 1
        residual = (val - arr(lb)) / (arr(lb+1) - arr(lb))
        loc = real(lb, ireals) + residual
        return
      else
        loc = real(lb+1) ! PetscFind returns the c-notation index
      endif
    end function

    ! return index+residula i where val is between arr(i) and arr(i+1)
    pure function search_sorted_bisection_r32(arr,val) result(res)
      real(REAL32),intent(in) :: arr(:)
      real(REAL32),intent(in) :: val
      real(REAL32) :: res
      real(REAL32) :: loc_increment
      integer(iintegers) :: i,j,k, N

      N = size(arr, kind=iintegers)
      if(N.le.LINEAR_SEARCH_LIMIT) then ! for small arrays it is quicker to do a linear search
        if(arr(1).ge.val) then
          res = 1
          return
        endif
        do i=2,N
          if(arr(i).ge.val) then
            loc_increment = (val - arr(i-1)) / (arr(i) - arr(i-1))
            res = real(i-1, kind=kind(res)) + loc_increment
            return
          endif
        enddo
        res = real(N, kind=kind(res))
        return
      endif

      i=lbound(arr,1)
      j=ubound(arr,1)

      if(arr(i).le.arr(j)) then ! ascending order
        do
          k=(i+j)/2
          if (val < arr(k)) then
            j=k
          else
            i=k
          endif
          if (i+1 >= j) then ! only single or tuple left
            ! i is left bound and j is right bound index
            if(i.eq.j) then
              loc_increment = 0
            else
              loc_increment = (val - arr(i)) / ( arr(j) - arr(i) )
            endif
            res= min(max(real(lbound(arr,1), REAL32), real(i, REAL32) + loc_increment), real(ubound(arr,1), REAL32)) ! return `real-numbered` location of val
            return
          endif
        end do
      else !descending order
        do
          k=(i+j)/2
          if (val > arr(k)) then
            j=k
          else
            i=k
          endif
          if (i+1 >= j) then ! only single or tuple left
            ! i is left bound and j is right bound index
            if(i.eq.j) then
              loc_increment = 0
            else
              loc_increment = (val - arr(j)) / ( arr(i) - arr(j) )
            endif
            res = min(max(real(lbound(arr,1), REAL32), real(j, REAL32) - loc_increment), real(ubound(arr,1), REAL32)) ! return `real-numbered` location of val
            return
          endif
        end do
      endif
    end function
    pure function search_sorted_bisection_r64(arr,val) result(res)
      real(REAL64),intent(in) :: arr(:)
      real(REAL64),intent(in) :: val
      real(REAL64) :: res
      real(REAL64) :: loc_increment
      integer(iintegers) :: i,j,k,N

      N = size(arr, kind=iintegers)
      if(N.le.LINEAR_SEARCH_LIMIT) then ! for small arrays it is quicker to do a linear search
        if(arr(1).ge.val) then
          res = 1
          return
        endif
        do i=2,N
          if(arr(i).ge.val) then
            loc_increment = (val - arr(i-1)) / (arr(i) - arr(i-1))
            res = real(i-1, kind=kind(res)) + loc_increment
            return
          endif
        enddo
        res = real(N, kind=kind(res))
        return
      endif

      i=lbound(arr,1)
      j=ubound(arr,1)

      if(arr(i).le.arr(j)) then ! ascending order
        do
          k=(i+j)/2
          if (val < arr(k)) then
            j=k
          else
            i=k
          endif
          if (i+1 >= j) then ! only single or tuple left
            ! i is left bound and j is right bound index
            if(i.eq.j) then
              loc_increment = 0
            else
              loc_increment = (val - arr(i)) / ( arr(j) - arr(i) )
            endif
            res= min(max(1._REAL64*lbound(arr,1), i + loc_increment), 1._REAL64*ubound(arr,1)) ! return `real-numbered` location of val
            return
          endif
        end do
      else !descending order
        do
          k=(i+j)/2
          if (val > arr(k)) then
            j=k
          else
            i=k
          endif
          if (i+1 >= j) then ! only single or tuple left
            ! i is left bound and j is right bound index
            if(i.eq.j) then
              loc_increment = 0
            else
              loc_increment = (val - arr(j)) / ( arr(i) - arr(j) )
            endif
            res = min(max(1._REAL64*lbound(arr,1), j - loc_increment), 1._REAL64*ubound(arr,1)) ! return `real-numbered` location of val
            return
          endif
        end do
      endif
    end function
end module
