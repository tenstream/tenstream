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
      use m_data_parameters, only: iintegers, ireals,zero,one
      implicit none

      private
      public :: interp_4d, interp_4d_recursive, interp_2d, interp_1d

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

      integer(iintegers),parameter :: i1=1

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

pure function spline(t,a0,a1)
        real(ireals),intent(in) :: t,a0(:),a1(:) ! t is weighting distance from a0
        real(ireals) :: spline(size(a1))
        real(ireals) :: s(size(a1))
        logical,parameter :: lspline = .False.
!        logical,parameter :: lspline = .True.

        if(lspline) then
          where( a0 .gt. a1 )
            s = t**2
          elsewhere(a0.le.a1)
            s = sqrt(t)
          endwhere
!          s = t**2*(3_ireals - 2_ireals*t)
        else
          s = t
        endif

!        spline = s*a1 + (one-s)*a0
        spline = a0 + s * ( a1-a0 )
        !      print *,'interpolating t,a0,a1',t,a0,a1,'==>',f
      end function

pure function interp_1d(t,a0)
        real(ireals),intent(in) :: t,a0(:) ! t is weighting distance from a0
        real(ireals) :: interp_1d
        integer(iintegers) :: i
        real(ireals) :: offset

        i = floor(t)
        offset = modulo(t,one)

        interp_1d = (one-offset) * a0(i) + offset * a0(min(i+1, size(a0)))
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
!          db5(:,i)  =  db6(:,2*i-1) + ( db6(:,2*i) - db6(:,2*i-1) ) * spline(weights(1))
          db5(:,i)  = spline( weights(1), db6(:,2*i-1), db6(:,2*i) )
        enddo
        do i=1,2**(Ndim-2)
!          db4(:,i)  =  db5(:,2*i-1) + ( db5(:,2*i) - db5(:,2*i-1) ) * spline(weights(2))
          db4(:,i)  = spline( weights(2), db5(:,2*i-1), db5(:,2*i) )
        enddo
        do i=1,2**(Ndim-3)
!          db3(:,i)  =  db4(:,2*i-1) + ( db4(:,2*i) - db4(:,2*i-1) ) * spline(weights(3))
          db3(:,i)  = spline( weights(3), db5(:,4*i-1), db5(:,4*i) )
        enddo
        do i=1,2**(Ndim-4)
!          db2(:,i)  =  db3(:,2*i-1) + ( db3(:,2*i) - db3(:,2*i-1) ) * spline(weights(4))
          db2(:,i)  = spline( weights(4), db5(:,3*i-1), db3(:,2*i) )
        enddo
        do i=1,2**(Ndim-5)
!          db1(:,i)  =  db2(:,2*i-1) + ( db2(:,2*i) - db2(:,2*i-1) ) * spline(weights(5))
          db1(:,i)  = spline( weights(5), db2(:,2*i-1), db2(:,2*i) )
        enddo

!        C(:)        =  db1(:,2  -1) + ( db1(:,2  ) - db1(:,2  -1) ) * spline(weights(6))
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
