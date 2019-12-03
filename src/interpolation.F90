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
  use iso_fortran_env, only: REAL32, REAL64
  use m_data_parameters, only: iintegers, irealLUT, mpiint, i1
  use m_helper_functions, only: approx, CHKERR, itoa, ftoa, &
    triangle_area_by_vertices, ind_nd_to_1d, ind_1d_to_nd, pnt_in_triangle
  implicit none

  private
  public :: interp_4d, interp_4d_recursive, interp_2d, &
    interp_1d, interp_vec_simplex_nd

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

  logical, parameter :: ldebug=.False.

  real(irealLUT), parameter :: interpolation_lattice_snapping=max(1e-3_irealLUT, epsilon(interpolation_lattice_snapping))
  real(irealLUT), parameter :: zero=0, one=1

  interface interp_1d
    module procedure interp_1d_r32, interp_1d_r64
  end interface

contains

  recursive subroutine interpn(n,a,t,i,res)
    real(irealLUT),intent(in) :: a(:,:),t(:)
    integer(iintegers),intent(in) :: n,i
    real(irealLUT),intent(out) :: res(:)

    real(irealLUT) :: a0(size(res)),a1(size(res))

    if(n.eq.1) then
      res = spline( t(1), a(:,i), a(:,i+1) )
    else
      call interpn(n-1,a,t,i,a0)
      call interpn(n-1,a,t,i+2**(n-1),a1)
      res = spline( t(n), a0, a1 )
    endif
  end subroutine

  pure elemental function spline(t,a0,a1)
    real(irealLUT),intent(in) :: t,a0,a1 ! t is weighting distance from a0
    real(irealLUT) :: spline
    real(irealLUT) :: s
    logical,parameter :: lspline = .False.
    ! logical,parameter :: lspline = .True.

    if(lspline) then
      if( a0 .gt. a1 ) then
        s = t**2
      else if (a0.le.a1) then
        s = sqrt(t)
      endif
      ! s = t**2*(3_irealLUT - 2_irealLUT*t)
    else
      s = t
    endif

    ! spline = s*a1 + (one-s)*a0
    spline = a0 + s * ( a1-a0 )
  end function

  function interp_1d_r32(t,a0) result(res)
    real(REAL32),intent(in) :: t, a0(:) ! t is weighting distance from [1, size(a0)]
    real(REAL32) :: res
    integer(iintegers) :: i
    real(REAL32) :: offset

    if(t.lt.1._REAL32 .or. t.gt.size(a0)) call CHKERR(1_mpiint, 'Cannot use interp_1d with weights outside of [0,1]')
    i = floor(t)
    offset = modulo(t,1._REAL32)
    if(approx(offset,0._REAL32)) then
      res = a0(i)
    else
      res = (1._REAL32-offset) * a0(i) + offset * a0(min(i+1, size(a0,kind=iintegers)))
    endif
  end function
  function interp_1d_r64(t,a0) result(res)
    real(REAL64),intent(in) :: t, a0(:) ! t is weighting distance from [1, size(a0)]
    real(REAL64) :: res
    integer(iintegers) :: i
    real(REAL64) :: offset

    if(t.lt.1._REAL64 .or. t.gt.size(a0)) call CHKERR(1_mpiint, 'Cannot use interp_1d with weights outside of [0,1]')
    i = floor(t)
    offset = modulo(t,1._REAL64)
    if(approx(offset,0._REAL64)) then
      res = a0(i)
    else
      res = (1._REAL64-offset) * a0(i) + offset * a0(min(i+1, size(a0,kind=iintegers)))
    endif
  end function

  pure subroutine interp_2d(pti, db, C)
    integer(iintegers),parameter :: Ndim=2
    real(irealLUT),intent(in) :: pti(Ndim), db(:,:,:)
    real(irealLUT),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    integer(iintegers) :: i,d
    real(irealLUT) :: weights(Ndim)
    real(irealLUT) :: db2(size(C),2**(Ndim ))
    real(irealLUT) :: db1(size(C),2**(Ndim-1))

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    do i=1,2**Ndim
      indices(:,i) = permu2d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db,d+i1,iintegers), indices(d,:) ) )
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
    real(irealLUT),intent(in) :: pti(Ndim), db(:,:,:,:,:)
    real(irealLUT),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    integer(iintegers) :: i,d
    real(irealLUT) :: weights(Ndim)
    real(irealLUT) :: db4(size(C),2**(Ndim ))
    real(irealLUT) :: db3(size(C),2**(Ndim-1))
    real(irealLUT) :: db2(size(C),2**(Ndim-2))
    real(irealLUT) :: db1(size(C),2**(Ndim-3))

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    do i=1,2**Ndim
      indices(:,i) = permu4d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db, d+i1, iintegers), indices(d,:) ) )
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
    real(irealLUT),intent(in) :: pti(Ndim), db(:,:,:,:,:)
    real(irealLUT),intent(out) :: C(:)

    integer(iintegers) :: indices(Ndim,2**Ndim),fpti(Ndim)
    real(irealLUT) :: weights(Ndim)
    real(irealLUT) :: bound_vals(size(C),2**Ndim)
    integer(iintegers) :: i,d

    ! First determine the array indices, where to look.
    fpti = floor(pti)
    weights = modulo(pti, one)

    do i=1,2**Ndim
      indices(:,i) = permu4d(:,i) + fpti
    enddo
    ! Make sure we dont recall a value outside of array dimensions
    do d=1,Ndim
      indices(d,:) = max( i1, min( ubound(db,d+i1, iintegers), indices(d,:) ) )
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
    real(irealLUT),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    real(irealLUT),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(irealLUT),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))

    integer(iintegers), allocatable :: nd_indices(:)

    if(ldebug) then
      allocate(nd_indices(size(pti)))
      call ind_1d_to_nd(db_offsets, size(db, dim=2, kind=iintegers), nd_indices)
      if(any(pti.lt.one).or.any(pti.gt.real(nd_indices, irealLUT))) then
        print *,'db dimensions', nd_indices
        print *,'pti', pti
        call CHKERR(1_mpiint, 'called with pti that does not fit database dimensions')
      endif
    endif

    call interp_vec_bilinear_recursive(pti, db, db_offsets, Cres)
  end subroutine

  pure subroutine interp_vec_simplex_1d(pti, interp_dim, db, db_offsets, Cres)
    real(irealLUT),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db dim(Ndimensions)
    integer(iintegers),intent(in) :: interp_dim ! dimension in which the 1D interpolation should happen
    real(irealLUT),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(irealLUT),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))

    real(irealLUT) :: A, B, wgt
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
    real(irealLUT),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    integer(iintegers),intent(in) :: interp_dims(:) ! dimensions in which the interpolation should happen
    real(irealLUT),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(irealLUT),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))
    real(irealLUT), dimension(2,5) :: points ! i.e. A, B, C, D, P

    integer(iintegers) :: i, ipnt
    integer(iintegers) :: nd_indices(size(pti))
    real(irealLUT) :: wgt_1d, interp_values(size(db,dim=1),3), interp_areas(4)

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
        interp_areas(ipnt) = triangle_area_by_vertices(  &
                points(:, modulo(ipnt,3_iintegers)+1  ), &
                points(:, modulo(ipnt+1,3_iintegers)+1), P)
      enddo
      do ipnt = 1, 3
        ! if the point lies on the newly made up line, quick exit with 1D interpolation
        if(approx(interp_areas(ipnt), 0._irealLUT, 10*sqrt(epsilon(0._irealLUT)))) then
          associate( p1 => modulo(ipnt,3_iintegers)+1, p2 => modulo(ipnt+1,3_iintegers)+1 )
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

  recursive pure subroutine interp_vec_bilinear_recursive(pti, db, db_offsets, Cres)
    real(irealLUT),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    real(irealLUT),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(irealLUT),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))
    Cres = 0
    call interp_vec_bilinear_recursive_(size(pti, kind=iintegers), pti, &
      db, db_offsets, 1_iintegers, 1._irealLUT, Cres)
  end subroutine

  recursive pure subroutine interp_vec_bilinear_recursive_(Ndim, pti, db, db_offsets, ofs, weight, Cres)
    integer(iintegers),intent(in) :: Ndim ! Number of dimensions that still need interpolation
    real(irealLUT),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    real(irealLUT),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    integer(iintegers),intent(in) :: ofs ! current database offset, in the end, will give raveled index
    real(irealLUT), intent(in) :: weight
    real(irealLUT),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))
    integer(iintegers) :: ind

    real(irealLUT) :: wgt_1d

    if(Ndim.eq.0) then
      do ind=1,size(db, dim=1)
        Cres(ind) = Cres(ind) + weight * db(ind,ofs)
      enddo
      return
    endif

    if(dim_needs_interpolation(pti(Ndim))) then
      ind = int(pti(Ndim))
      wgt_1d = pti(Ndim) - real(ind, irealLUT) ! === modulo(pti(Ndim), one)
      call interp_vec_bilinear_recursive_(Ndim-1, pti, &
        db, db_offsets, ofs + db_offsets(Ndim) * (ind-1), weight * (one-wgt_1d), Cres)
      call interp_vec_bilinear_recursive_(Ndim-1, pti, &
        db, db_offsets, ofs + db_offsets(Ndim) * ind, weight * wgt_1d, Cres)
    else ! snap to nearest val
      call interp_vec_bilinear_recursive_(Ndim-1, pti, &
        db, db_offsets, ofs + db_offsets(Ndim) * (nint(pti(Ndim))-1), weight, Cres)
    endif
  end subroutine

  recursive subroutine interp_vec_simplex_recursive(Ndim, pti, interp_dims, db, db_offsets, Cres)
    use m_helper_functions, only : distance, triangle_area_by_vertices
    integer(iintegers),intent(in) :: Ndim
    real(irealLUT),intent(in) :: pti(:) ! weigths/indices in the respective unraveled db, dim(Ndimensions)
    integer(iintegers),intent(in) :: interp_dims(:) ! dimensions in which the interpolation should happen
    real(irealLUT),intent(in) :: db(:,:) ! first dimension is the vector dimension, ie if just one scalar should be interpolated, call it with shape [1, ravel(db)]
    integer(iintegers),intent(in) :: db_offsets(:) ! offsets of the db dim(Ndimensions)
    real(irealLUT),intent(out) :: Cres(:) ! output, has the dimension(size(db,dim=1))

    integer(iintegers) :: Ninterpdim
    real(irealLUT) :: pti_intermediate(size(pti)), db_intermediate(size(db,dim=1), 2), wgt_1d

    Ninterpdim=size(interp_dims)
    ! Interpolate first two dimensions with simplex and then one 1d interpolation

    ! Simplex:
    pti_intermediate = pti
    pti_intermediate(interp_dims(Ninterpdim)) = floor(pti(interp_dims(Ninterpdim)))
    select case (Ndim)
    case(2)
      call interp_vec_simplex_1d(pti_intermediate, interp_dims(1), &
              db, db_offsets, db_intermediate(:,1))
    case(3)
      call interp_vec_simplex_2d(pti_intermediate, interp_dims(1:Ninterpdim-1), &
              db, db_offsets, db_intermediate(:,1))
    case(4,5,6,7)
      call interp_vec_simplex_recursive(Ndim-1, pti_intermediate, interp_dims(1:Ninterpdim-1), &
              db, db_offsets, db_intermediate(:,1))
    case default
      call CHKERR(1_mpiint, 'interp_vec_simplex_recursive not implemented for '// &
              itoa(size(interp_dims, kind=iintegers))//' dimensions')
    end select

    pti_intermediate(interp_dims(Ninterpdim)) = ceiling(pti(interp_dims(Ninterpdim)))
    select case (Ndim)
    case(2)
      call interp_vec_simplex_1d(pti_intermediate, interp_dims(1), &
              db, db_offsets, db_intermediate(:,2))
    case(3)
      call interp_vec_simplex_2d(pti_intermediate, interp_dims(1:Ninterpdim-1), &
              db, db_offsets, db_intermediate(:,2))
    case(4,5,6,7)
      call interp_vec_simplex_recursive(Ndim-1, pti_intermediate, interp_dims(1:Ninterpdim-1), &
              db, db_offsets, db_intermediate(:,2))
    case default
      call CHKERR(1_mpiint, 'interp_vec_simplex_recursive not implemented for '// &
              itoa(size(interp_dims, kind=iintegers))//' dimensions')
    end select

    wgt_1d = modulo(pti(interp_dims(Ninterpdim)), one)
    Cres = spline(wgt_1d, db_intermediate(:,1), db_intermediate(:,2))
  end subroutine

  elemental function dim_needs_interpolation(pti)
    real(irealLUT),intent(in) :: pti
    logical :: dim_needs_interpolation
    if( pti - int(pti) .lt. interpolation_lattice_snapping ) then
      dim_needs_interpolation = .False.
    elseif(pti - int(pti) .gt. 1._irealLUT - interpolation_lattice_snapping) then
      dim_needs_interpolation = .False.
    else
      dim_needs_interpolation = .True.
    endif
  end function


end module

!      program main
!      use interpolation
!
!      integer(iintegers),parameter :: Ndim=3
!      real(irealLUT) :: a(2**Ndim),t(Ndim)
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
