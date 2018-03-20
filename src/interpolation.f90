!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

module m_tenstream_interpolation
  use m_data_parameters, only: iintegers, ireals, mpiint, zero,one, i1
  use m_helper_functions, only: approx, CHKERR, itoa, ftoa, &
    triangle_area_by_vertices, ind_nd_to_1d, ind_1d_to_nd, pnt_in_triangle
  implicit none

  private
  public :: interp_4d, interp_4d_recursive, interp_2d, &
    interp_1d, interp_vec_1d, &
    interp_vec_simplex_nd

  ! a has the bounds on axes
  ! t has the distance weights

  integer(iintegers) :: permu2d(2,2**2)
  DATA permu2d(:, 1) / 0,  0 /
  DATA permu2d(:, 2) / 1,  0 /
  DATA permu2d(:, 3) / 0,  1 /
  DATA permu2d(:, 4) / 1,  1 /

  integer(iintegers) :: permu4d(4,2**4)
  DATA permu4d(:, 1) / 0,  0,  0,  0 /
  DATA permu4d(:, 2) / 1,  0,  0,  0 /
  DATA permu4d(:, 3) / 0,  1,  0,  0 /
  DATA permu4d(:, 4) / 1,  1,  0,  0 /
  DATA permu4d(:, 5) / 0,  0,  1,  0 /
  DATA permu4d(:, 6) / 1,  0,  1,  0 /
  DATA permu4d(:, 7) / 0,  1,  1,  0 /
  DATA permu4d(:, 8) / 1,  1,  1,  0 /
  DATA permu4d(:, 9) / 0,  0,  0,  1 /
  DATA permu4d(:,10) / 1,  0,  0,  1 /
  DATA permu4d(:,11) / 0,  1,  0,  1 /
  DATA permu4d(:,12) / 1,  1,  0,  1 /
  DATA permu4d(:,13) / 0,  0,  1,  1 /
  DATA permu4d(:,14) / 1,  0,  1,  1 /
  DATA permu4d(:,15) / 0,  1,  1,  1 /
  DATA permu4d(:,16) / 1,  1,  1,  1 /

  integer(iintegers) :: permu6d(6,2**6)
  DATA permu6d(:, 1 )  / 0,  0,  0,  0,  0,  0 /
  DATA permu6d(:, 2 )  / 1,  0,  0,  0,  0,  0 /
  DATA permu6d(:, 3 )  / 0,  1,  0,  0,  0,  0 /
  DATA permu6d(:, 4 )  / 1,  1,  0,  0,  0,  0 /
  DATA permu6d(:, 5 )  / 0,  0,  1,  0,  0,  0 /
  DATA permu6d(:, 6 )  / 1,  0,  1,  0,  0,  0 /
  DATA permu6d(:, 7 )  / 0,  1,  1,  0,  0,  0 /
  DATA permu6d(:, 8 )  / 1,  1,  1,  0,  0,  0 /
  DATA permu6d(:, 9 )  / 0,  0,  0,  1,  0,  0 /
  DATA permu6d(:, 10 ) / 1,  0,  0,  1,  0,  0 /
  DATA permu6d(:, 11 ) / 0,  1,  0,  1,  0,  0 /
  DATA permu6d(:, 12 ) / 1,  1,  0,  1,  0,  0 /
  DATA permu6d(:, 13 ) / 0,  0,  1,  1,  0,  0 /
  DATA permu6d(:, 14 ) / 1,  0,  1,  1,  0,  0 /
  DATA permu6d(:, 15 ) / 0,  1,  1,  1,  0,  0 /
  DATA permu6d(:, 16 ) / 1,  1,  1,  1,  0,  0 /
  DATA permu6d(:, 17 ) / 0,  0,  0,  0,  1,  0 /
  DATA permu6d(:, 18 ) / 1,  0,  0,  0,  1,  0 /
  DATA permu6d(:, 19 ) / 0,  1,  0,  0,  1,  0 /
  DATA permu6d(:, 20 ) / 1,  1,  0,  0,  1,  0 /
  DATA permu6d(:, 21 ) / 0,  0,  1,  0,  1,  0 /
  DATA permu6d(:, 22 ) / 1,  0,  1,  0,  1,  0 /
  DATA permu6d(:, 23 ) / 0,  1,  1,  0,  1,  0 /
  DATA permu6d(:, 24 ) / 1,  1,  1,  0,  1,  0 /
  DATA permu6d(:, 25 ) / 0,  0,  0,  1,  1,  0 /
  DATA permu6d(:, 26 ) / 1,  0,  0,  1,  1,  0 /
  DATA permu6d(:, 27 ) / 0,  1,  0,  1,  1,  0 /
  DATA permu6d(:, 28 ) / 1,  1,  0,  1,  1,  0 /
  DATA permu6d(:, 29 ) / 0,  0,  1,  1,  1,  0 /
  DATA permu6d(:, 30 ) / 1,  0,  1,  1,  1,  0 /
  DATA permu6d(:, 31 ) / 0,  1,  1,  1,  1,  0 /
  DATA permu6d(:, 32 ) / 1,  1,  1,  1,  1,  0 /
  DATA permu6d(:, 33 ) / 0,  0,  0,  0,  0,  1 /
  DATA permu6d(:, 34 ) / 1,  0,  0,  0,  0,  1 /
  DATA permu6d(:, 35 ) / 0,  1,  0,  0,  0,  1 /
  DATA permu6d(:, 36 ) / 1,  1,  0,  0,  0,  1 /
  DATA permu6d(:, 37 ) / 0,  0,  1,  0,  0,  1 /
  DATA permu6d(:, 38 ) / 1,  0,  1,  0,  0,  1 /
  DATA permu6d(:, 39 ) / 0,  1,  1,  0,  0,  1 /
  DATA permu6d(:, 40 ) / 1,  1,  1,  0,  0,  1 /
  DATA permu6d(:, 41 ) / 0,  0,  0,  1,  0,  1 /
  DATA permu6d(:, 42 ) / 1,  0,  0,  1,  0,  1 /
  DATA permu6d(:, 43 ) / 0,  1,  0,  1,  0,  1 /
  DATA permu6d(:, 44 ) / 1,  1,  0,  1,  0,  1 /
  DATA permu6d(:, 45 ) / 0,  0,  1,  1,  0,  1 /
  DATA permu6d(:, 46 ) / 1,  0,  1,  1,  0,  1 /
  DATA permu6d(:, 47 ) / 0,  1,  1,  1,  0,  1 /
  DATA permu6d(:, 48 ) / 1,  1,  1,  1,  0,  1 /
  DATA permu6d(:, 49 ) / 0,  0,  0,  0,  1,  1 /
  DATA permu6d(:, 50 ) / 1,  0,  0,  0,  1,  1 /
  DATA permu6d(:, 51 ) / 0,  1,  0,  0,  1,  1 /
  DATA permu6d(:, 52 ) / 1,  1,  0,  0,  1,  1 /
  DATA permu6d(:, 53 ) / 0,  0,  1,  0,  1,  1 /
  DATA permu6d(:, 54 ) / 1,  0,  1,  0,  1,  1 /
  DATA permu6d(:, 55 ) / 0,  1,  1,  0,  1,  1 /
  DATA permu6d(:, 56 ) / 1,  1,  1,  0,  1,  1 /
  DATA permu6d(:, 57 ) / 0,  0,  0,  1,  1,  1 /
  DATA permu6d(:, 58 ) / 1,  0,  0,  1,  1,  1 /
  DATA permu6d(:, 59 ) / 0,  1,  0,  1,  1,  1 /
  DATA permu6d(:, 60 ) / 1,  1,  0,  1,  1,  1 /
  DATA permu6d(:, 61 ) / 0,  0,  1,  1,  1,  1 /
  DATA permu6d(:, 62 ) / 1,  0,  1,  1,  1,  1 /
  DATA permu6d(:, 63 ) / 0,  1,  1,  1,  1,  1 /
  DATA permu6d(:, 64 ) / 1,  1,  1,  1,  1,  1 /

  logical, parameter :: ldebug=.True.

contains

  recursive subroutine interpn(n,a,t,i,res)
    real(ireals),intent(in) :: a(:,:),t(:)
    integer(iintegers),intent(in) :: n,i
    real(ireals),intent(out) :: res(:)

    real(ireals) :: a0(size(res)),a1(size(res))

    if(n.eq.1) then
      res = spline( t(1), a(:,i), a(:,i+1) )
    else
      call interpn(n-1,a,t,i,a0)
      call interpn(n-1,a,t,i+2**(n-1),a1)
      res = spline( t(n), a0, a1 )
    endif
  end subroutine

  pure elemental function spline(t,a0,a1)
    real(ireals),intent(in) :: t,a0,a1 ! t is weighting distance from a0
    real(ireals) :: spline
    real(ireals) :: s
    logical,parameter :: lspline = .False.
    !        logical,parameter :: lspline = .True.

    if(lspline) then
      if( a0 .gt. a1 ) then
        s = t**2
      else if (a0.le.a1) then
        s = sqrt(t)
      endif
      !          s = t**2*(3_ireals - 2_ireals*t)
    else
      s = t
    endif

    !        spline = s*a1 + (one-s)*a0
    spline = a0 + s * ( a1-a0 )
    !      print *,'interpolating t,a0,a1',t,a0,a1,'==>',f
  end function

  function interp_1d(t,a0)
    real(ireals),intent(in) :: t, a0(:) ! t is weighting distance from [1, size(a0)]
    real(ireals) :: interp_1d
    integer(iintegers) :: i
    real(ireals) :: offset

    if(t.lt.one .or. t.gt.size(a0)) call CHKERR(1_mpiint, 'Cannot use interp_1d with weights outside of [0,1]')
    i = floor(t)
    offset = modulo(t,one)
    if(approx(offset,zero)) then
      interp_1d = a0(i)
    else
      interp_1d = (one-offset) * a0(i) + offset * a0(min(i+1, size(a0)))
    endif
  end function

  function interp_vec_1d(t,a0)
    real(ireals),intent(in) :: t       ! t is weighting distance from [1, size(a0)]
    real(ireals),intent(in) :: a0(:,:) ! first dimension is vector which is interpolated for
    real(ireals) :: interp_vec_1d(size(a0,dim=1))
    integer(iintegers) :: i
    real(ireals) :: offset

    if(t.lt.one .or. t.gt.size(a0, dim=2)) call CHKERR(1_mpiint, &
      'Cannot use interp_1d with weights outside of the array bounds '//ftoa(t) &
      //' : '//itoa(size(a0,dim=1,kind=iintegers))//','//itoa(size(a0,dim=2,kind=iintegers)))
    i = floor(t)
    offset = modulo(t,one)
    if(approx(offset,zero)) then
      interp_vec_1d = a0(:,i)
    else
      !interp_vec_1d = (one-offset) * a0(:,i) + offset * a0(:,min(i+1, size(a0, dim=2)))
      interp_vec_1d = spline(offset, a0(:,i), a0(:,i+1))
    endif
  end function

  subroutine interp_6d_recursive(pti, db, C)
    integer(iintegers),parameter :: Ndim=6
    real(ireals),intent(in) :: pti(Ndim), db(:,:,:,:,:,:,:)
    real(ireals),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    real(ireals) :: weights(Ndim)
    real(ireals) :: bound_vals(size(C),2**Ndim)
    integer(iintegers) :: i,d

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    !        print *,'interp6d',pti,'weights',weights
    do i=1,2**Ndim
      indices(:,i) = permu6d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db,d+i1), indices(d,:) ) )
    enddo

    ! Then get the corner values of hypercube
    do i=1,2**Ndim
      bound_vals(:,i) = db(:, indices(1,i),indices(2,i),indices(3,i),indices(4,i), indices(5,i), indices(6,i) )
      !          print *,'interp6d',i,'ind',indices(:,i),'bounds',bound_vals(:,i)
    enddo

    ! And plug bound_vals and weights into recursive interpolation...
    call interpn(Ndim,bound_vals,weights,i1, C)
  end subroutine
  pure subroutine interp_6d(pti, db, C)
    integer(iintegers),parameter :: Ndim=6
    real(ireals),intent(in) :: pti(Ndim), db(:,:,:,:,:,:,:)
    real(ireals),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    integer(iintegers) :: i,d,ind(6)

    real(ireals) :: weights(Ndim)
    real(ireals) :: db6(size(C),2**Ndim)
    real(ireals) :: db5(size(C),2**(Ndim-1))
    real(ireals) :: db4(size(C),2**(Ndim-2))
    real(ireals) :: db3(size(C),2**(Ndim-3))
    real(ireals) :: db2(size(C),2**(Ndim-4))
    real(ireals) :: db1(size(C),2**(Ndim-5))

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    !        print *,'interp6d',pti,'weights',weights
    do i=1,2**Ndim
      indices(:,i) = permu6d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db,d+i1), indices(d,:) ) )
    enddo
    ! Then get the corner values of hypercube
    do i=1,2**Ndim
      ind = indices(:,i)
      db6(:,i) = db(:, ind(1),ind(2),ind(3),ind(4), ind(5), ind(6) )
    enddo

    ! Permutations for 1st axis
    do i=1,2**(Ndim-1)
      db5(:,i)  = spline( weights(1), db6(:,2*i-1), db6(:,2*i) )
    enddo
    do i=1,2**(Ndim-2)
      db4(:,i)  = spline( weights(2), db5(:,2*i-1), db5(:,2*i) )
    enddo
    do i=1,2**(Ndim-3)
      db3(:,i)  = spline( weights(3), db5(:,4*i-1), db5(:,4*i) )
    enddo
    do i=1,2**(Ndim-4)
      db2(:,i)  = spline( weights(4), db5(:,3*i-1), db3(:,2*i) )
    enddo
    do i=1,2**(Ndim-5)
      db1(:,i)  = spline( weights(5), db2(:,2*i-1), db2(:,2*i) )
    enddo

    C(:)  = spline( weights(6), db1(:,2), db1(:,1) )
  end subroutine
  pure subroutine interp_2d(pti, db, C)
    integer(iintegers),parameter :: Ndim=2
    real(ireals),intent(in) :: pti(Ndim), db(:,:,:)
    real(ireals),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    integer(iintegers) :: i,d
    real(ireals) :: weights(Ndim)
    real(ireals) :: db2(size(C),2**(Ndim ))
    real(ireals) :: db1(size(C),2**(Ndim-1))

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    do i=1,2**Ndim
      indices(:,i) = permu2d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db,d+i1), indices(d,:) ) )
    enddo

    ! Then get the corner values of hypercube
    do i=1,2**Ndim
      db2(:,i) = db(:, indices(1,i), indices(2,i))
    enddo

    ! Permutations for 1st axis
    do i=1,2**(Ndim-1)
      db1(:,i)  =  db2(:,2*i-1) + ( db2(:,2*i) - db2(:,2*i-1) ) * (weights(1))
    enddo

    C(:)        =  db1(:,2  -1) + ( db1(:,2  ) - db1(:,2  -1) ) * (weights(2))
  end subroutine
  pure subroutine interp_4d(pti, db, C)
    integer(iintegers),parameter :: Ndim=4
    real(ireals),intent(in) :: pti(Ndim), db(:,:,:,:,:)
    real(ireals),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    integer(iintegers) :: i,d
    real(ireals) :: weights(Ndim)
    real(ireals) :: db4(size(C),2**(Ndim ))
    real(ireals) :: db3(size(C),2**(Ndim-1))
    real(ireals) :: db2(size(C),2**(Ndim-2))
    real(ireals) :: db1(size(C),2**(Ndim-3))

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    do i=1,2**Ndim
      indices(:,i) = permu4d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db,d+i1), indices(d,:) ) )
    enddo

    ! Then get the corner values of hypercube
    do i=1,2**Ndim
      db4(:,i) = db(:, indices(1,i),indices(2,i),indices(3,i),indices(4,i) )
    enddo

    ! Permutations for 1st axis
    do i=1,2**(Ndim-1)
      db3(:,i)  =  db4(:,2*i-1) + ( db4(:,2*i) - db4(:,2*i-1) ) * (weights(1))
    enddo
    do i=1,2**(Ndim-2)
      db2(:,i)  =  db3(:,2*i-1) + ( db3(:,2*i) - db3(:,2*i-1) ) * (weights(2))
    enddo
    do i=1,2**(Ndim-3)
      db1(:,i)  =  db2(:,2*i-1) + ( db2(:,2*i) - db2(:,2*i-1) ) * (weights(3))
    enddo

    C(:)        =  db1(:,2  -1) + ( db1(:,2  ) - db1(:,2  -1) ) * (weights(4))
  end subroutine
  subroutine interp_4d_recursive(pti, db, C)
    integer(iintegers),parameter :: Ndim=4
    real(ireals),intent(in) :: pti(Ndim), db(:,:,:,:,:)
    real(ireals),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    real(ireals) :: weights(Ndim)
    real(ireals) :: bound_vals(size(C),2**Ndim)
    integer(iintegers) :: i,d

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    do i=1,2**Ndim
      indices(:,i) = permu4d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db,d+i1), indices(d,:) ) )
    enddo

    ! Then get the corner values of hypercube
    do i=1,2**Ndim
      bound_vals(:,i) = db(:, indices(1,i),indices(2,i),indices(3,i),indices(4,i) )
    enddo

    ! And plug bound_vals and weights into recursive interpolation...
    call interpn(Ndim,bound_vals,weights,i1, C)
  end subroutine

  ! http://www.hpl.hp.com/techreports/2002/HPL-2002-320.pdf
  subroutine interp_vec_simplex_nd(pti, db, db_offsets, Cres)
    real(ireals),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    real(ireals),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(ireals),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))

    integer(iintegers), allocatable :: interp_dims(:)
    integer(iintegers) :: nd_indices(size(pti))

    if(ldebug) then
      nd_indices = ind_1d_to_nd(db_offsets, size(db, dim=2, kind=iintegers))
      if(any(pti.lt.one).or.any(pti.gt.nd_indices)) then
        print *,'db dimensions', nd_indices
        print *,'pti', pti
        call CHKERR(1_mpiint, 'called with pti that does not fit database dimensions')
      endif
    endif

    interp_dims = get_dims_that_need_interpolation(pti)

    select case (size(interp_dims))
    case(0)
      Cres = db(:, ind_nd_to_1d(db_offsets, nint(pti, iintegers)))
    case(1)
      call interp_vec_simplex_1d(pti, interp_dims(1), db, db_offsets, Cres)
    case(2)
      call interp_vec_simplex_2d(pti, interp_dims, db, db_offsets, Cres)
    case(3:10)
      !call interp_vec_bilinear_recursive(size(interp_dims), pti, interp_dims, db, db_offsets, Cres)
      call interp_vec_simplex_recursive(size(interp_dims, kind=iintegers), pti, interp_dims, db, db_offsets, Cres)
    case default
      call CHKERR(1_mpiint, 'interp_vec_simplex not implemented for '//itoa(size(interp_dims, kind=iintegers))//' dimensions')
    end select
  end subroutine

  subroutine interp_vec_simplex_1d(pti, interp_dim, db, db_offsets, Cres)
    real(ireals),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db dim(Ndimensions)
    integer(iintegers),intent(in) :: interp_dim ! dimension in which the 1D interpolation should happen
    real(ireals),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(ireals),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))

    real(ireals) :: A, B, wgt
    integer(iintegers) :: nd_indices(size(db_offsets)), indA, indB

    nd_indices = nint(pti, iintegers)

    A = floor(pti(interp_dim))
    B = ceiling(pti(interp_dim))
    wgt = modulo(pti(interp_dim), one)

    nd_indices(interp_dim) = nint(A)
    indA = ind_nd_to_1d(db_offsets, nd_indices)

    nd_indices(interp_dim) = nint(B)
    indB = ind_nd_to_1d(db_offsets, nd_indices)

    Cres = spline(wgt, db(:,indA), db(:,indB)) ! (one-wgt) * db(:,indA) + wgt * db(:,indB)
  end subroutine

  subroutine interp_vec_simplex_2d(pti, interp_dims, db, db_offsets, Cres)
    use m_helper_functions, only : distance, triangle_area_by_vertices
    real(ireals),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    integer(iintegers),intent(in) :: interp_dims(:) ! dimensions in which the interpolation should happen
    real(ireals),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(ireals),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))
    real(ireals), dimension(2,5) :: points ! i.e. A, B, C, D, P
    real(ireals) :: areas(3) ! areas between vertices ABP, BCP, ACP

    integer(iintegers) :: i, ipnt
    integer(iintegers) :: nd_indices(size(pti))
    real(ireals) :: wgt_1d, interp_values(size(db,dim=1),3), interp_areas(4)

    associate( A => points(:,1), &
        B => points(:,2), &
        C => points(:,3), &
        D => points(:,4), &
        P => points(:,5) )

      P = [pti(interp_dims(1)), pti(interp_dims(2))]

      A = [floor  (P(1)), floor  (P(2))]
      B = [ceiling(P(1)), floor  (P(2))]
      C = [floor  (P(1)), ceiling(P(2))]
      D = [ceiling(P(1)), ceiling(P(2))]

      nd_indices = nint(pti, iintegers)

      if(pnt_in_triangle(B,C,D,P)) then
        A = D
      endif

      do ipnt = 1, 3
        nd_indices(interp_dims) = nint(points(:,ipnt))
        interp_values(:,ipnt) = db(:, ind_nd_to_1d(db_offsets, nd_indices))
        interp_areas(ipnt) = triangle_area_by_vertices(points(:,modulo(ipnt,3)+1), points(:,modulo(ipnt+1,3)+1), P)
      enddo
      do ipnt = 1, 3
        ! if the point lies on the newly made up line, quick exit with 1D interpolation
        if(approx(interp_areas(ipnt), zero, 10*sqrt(epsilon(zero)))) then
          associate( p1 => modulo(ipnt,3)+1, p2 => modulo(ipnt+1,3)+1 )
            wgt_1d = interp_areas(p2) / (interp_areas(p1) + interp_areas(p2))
            Cres = spline(wgt_1d, interp_values(:,p1), interp_values(:,p2))
          end associate
          return
        endif
      enddo

      interp_areas(4) = triangle_area_by_vertices(A,B,C)

      if(sum(interp_areas(1:3))-10*sqrt(epsilon(one)).ge.interp_areas(4)) then
        print *,'interp values', int(interp_values), 'areas', interp_areas, &
          '->', sum(interp_areas(1:3))-sqrt(epsilon(one)*10)
        call CHKERR(1_mpiint, 'whoops, inner triangle areas are bigger than total. This means we`ve done something wrong')
      endif

      do i = 1, size(db, dim=1)
        Cres(i) = dot_product(interp_values(i,:), interp_areas(1:3)) / interp_areas(4)
      enddo
    end associate
  end subroutine

  recursive subroutine interp_vec_bilinear_recursive(Ndim, pti, interp_dims, db, db_offsets, Cres)
    use m_helper_functions, only : distance, triangle_area_by_vertices
    integer(iintegers),intent(in) :: Ndim
    real(ireals),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    integer(iintegers),intent(in) :: interp_dims(:) ! dimensions in which the interpolation should happen
    real(ireals),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(ireals),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))

    integer(iintegers) :: Ninterpdim
    real(ireals) :: pti_intermediate(size(pti)), db_intermediate(size(db,dim=1), 2), wgt_1d

    Ninterpdim=size(interp_dims)
    ! Interpolate first two dimensions with simplex and then one 1d interpolation

    ! Simplex:
    pti_intermediate = pti
    pti_intermediate(interp_dims(Ninterpdim)) = floor(pti(interp_dims(Ninterpdim)))
    if(Ndim.eq.2) then
      call interp_vec_simplex_1d(pti_intermediate, interp_dims(1), db, db_offsets, db_intermediate(:,1))
    else
      call interp_vec_bilinear_recursive(Ndim-1, pti_intermediate, interp_dims(1:Ninterpdim-1), db, db_offsets, db_intermediate(:,1))
    endif

    pti_intermediate(interp_dims(Ninterpdim)) = ceiling(pti(interp_dims(Ninterpdim)))
    if(Ndim.eq.2) then
      call interp_vec_simplex_1d(pti_intermediate, interp_dims(1), db, db_offsets, db_intermediate(:,2))
    else
      call interp_vec_bilinear_recursive(Ndim-1, pti_intermediate, interp_dims(1:Ninterpdim-1), db, db_offsets, db_intermediate(:,2))
    endif

    wgt_1d = one + modulo(pti(interp_dims(Ninterpdim)), one)
    Cres = interp_vec_1d(wgt_1d, db_intermediate)

  end subroutine

  recursive subroutine interp_vec_simplex_recursive(Ndim, pti, interp_dims, db, db_offsets, Cres)
    use m_helper_functions, only : distance, triangle_area_by_vertices
    integer(iintegers),intent(in) :: Ndim
    real(ireals),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    integer(iintegers),intent(in) :: interp_dims(:) ! dimensions in which the interpolation should happen
    real(ireals),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(ireals),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))

    integer(iintegers) :: Ninterpdim
    real(ireals) :: pti_intermediate(size(pti)), db_intermediate(size(db,dim=1), 2), wgt_1d

    Ninterpdim=size(interp_dims)
    ! Interpolate first two dimensions with simplex and then one 1d interpolation

    ! Simplex:
    pti_intermediate = pti
    pti_intermediate(interp_dims(Ninterpdim)) = floor(pti(interp_dims(Ninterpdim)))
    if(Ndim.eq.3) then
      call interp_vec_simplex_2d(pti_intermediate, interp_dims(1:Ninterpdim-1), db, db_offsets, db_intermediate(:,1))
    else
      call interp_vec_simplex_recursive(Ndim-1, pti_intermediate, interp_dims(1:Ninterpdim-1), db, db_offsets, db_intermediate(:,1))
    endif

    pti_intermediate(interp_dims(Ninterpdim)) = ceiling(pti(interp_dims(Ninterpdim)))
    if(Ndim.eq.3) then
      call interp_vec_simplex_2d(pti_intermediate, interp_dims(1:Ninterpdim-1), db, db_offsets, db_intermediate(:,2))
    else
      call interp_vec_simplex_recursive(Ndim-1, pti_intermediate, interp_dims(1:Ninterpdim-1), db, db_offsets, db_intermediate(:,2))
    endif

    wgt_1d = one + modulo(pti(interp_dims(Ninterpdim)), one)
    Cres = interp_vec_1d(wgt_1d, db_intermediate)

  end subroutine

  pure function get_dims_that_need_interpolation(pti)
    real(ireals),intent(in) :: pti(:)
    integer(iintegers), allocatable :: get_dims_that_need_interpolation(:)
    integer(iintegers) :: i
    integer(iintegers) :: rank(size(pti))
    logical :: linterp(size(pti))
    linterp = dim_needs_interpolation(pti)
    rank = [ (i, i=1,size(pti)) ]

    allocate(get_dims_that_need_interpolation(count(linterp)))
    get_dims_that_need_interpolation = pack(rank, mask=linterp)
  end function

  pure elemental function dim_needs_interpolation(pti)
    real(ireals),intent(in) :: pti
    logical :: dim_needs_interpolation
    dim_needs_interpolation = .not.approx(modulo(pti, one), zero, sqrt(epsilon(zero)))
  end function


end module

!      program main
!      use interpolation
!
!      integer(iintegers),parameter :: Ndim=3
!      real(ireals) :: a(2**Ndim),t(Ndim)
!
!      a = [3,3,2,2,1,1,0,0]
!      t = [.5,.5,.5]
!      print *,interpn(Ndim,a,t,1),'should be 1.5'
!
!      a=[3,3,9,9,2,2,11,21] ; t=[.7,.5,.5]
!      print *,interpn(Ndim,a,t,1),'should be 8'
!
!      a=[3,3,9,9,2,2,11,21] ; t=[1,1,1]
!      print *,interpn(Ndim,a,t,1),'should be ',a(2**Ndim)
!
!      a=[ 0.29738766,  0.69095583,  0.1003632 ,  0.68999427,  0.0239015 , 0.88141948,  0.85422108,  0.78470368]
!      t = [.5,.5,.5]
!      print *,interpn(Ndim,a,t,1),'should be 0.54'
!
!      end program
