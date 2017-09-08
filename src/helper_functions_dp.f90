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

module m_helper_functions_dp
      use m_data_parameters,only : iintegers,ireal_dp,imp_real_dp,imp_int,imp_logical,mpiint, pi_dp
      use m_helper_functions, only: CHKERR
      use mpi

      implicit none

      private
      public imp_bcast,norm,deg2rad,rad2deg,rmse,mean,approx,rel_approx,delta_scale_optprop,delta_scale,cumsum,inc, &
          mpi_logical_and,mpi_logical_or,imp_allreduce_min,imp_allreduce_max,imp_reduce_sum,                &
          pnt_in_triangle, compute_normal_3d, hit_plane, spherical_2_cartesian, distance_to_edge, distance_to_triangle_edges, rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec

      interface imp_bcast
        module procedure imp_bcast_real_1d,imp_bcast_real_2d,imp_bcast_real_3d,imp_bcast_real_5d,imp_bcast_int_1d,imp_bcast_int_2d,imp_bcast_int,imp_bcast_real,imp_bcast_logical
      end interface

      integer(mpiint) :: mpierr
      real(ireal_dp),parameter :: zero=0, one=1

    contains
      pure elemental subroutine inc(x,i)
          real(ireal_dp),intent(inout) :: x
          real(ireal_dp),intent(in) :: i
          x=x+i
      end subroutine

      pure function norm(v)
          real(ireal_dp) :: norm
          real(ireal_dp),intent(in) :: v(:)
          norm = sqrt(dot_product(v,v))
      end function

      elemental function deg2rad(deg)
          real(ireal_dp) :: deg2rad
          real(ireal_dp),intent(in) :: deg
          deg2rad = deg * pi_dp / 180
      end function
      elemental function rad2deg(rad)
        real(ireal_dp) :: rad2deg
        real(ireal_dp),intent(in) :: rad
        rad2deg = rad / pi_dp * 180
      end function

      pure function rmse(a,b)
          real(ireal_dp) :: rmse(2)
          real(ireal_dp),intent(in) :: a(:),b(:)
          rmse(1) = sqrt( mean( (a-b)**2 ) )
          rmse(2) = rmse(1)/max( mean(b), epsilon(rmse) )
      end function

      pure function mean(arr)
          real(ireal_dp) :: mean
          real(ireal_dp),intent(in) :: arr(:)
          mean = sum(arr)/size(arr)
      end function

      elemental logical function approx(a,b,precision)
          real(ireal_dp),intent(in) :: a,b
          real(ireal_dp),intent(in),optional :: precision
          real(ireal_dp) :: factor
          if(present(precision) ) then
            factor = precision
          else
            factor = 10 * epsilon(b)
          endif
          if( a.le.b+factor .and. a.ge.b-factor ) then
            approx = .True.
          else
            approx = .False.
          endif
      end function
      elemental logical function rel_approx(a,b,precision)
          real(ireal_dp),intent(in) :: a,b
          real(ireal_dp),intent(in),optional :: precision
          real(ireal_dp) :: factor,rel_error
          if(present(precision) ) then
            factor = precision
          else
            factor = 10*epsilon(b)
          endif
          rel_error = abs( (a-b)/ max(epsilon(a), ( (a+b)*.5_ireal_dp ) ) )

          if( rel_error .lt. precision ) then
            rel_approx = .True.
          else
            rel_approx = .False.
          endif
      end function


      function mpi_logical_and(comm, lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_and
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_and, 1_mpiint, imp_logical, MPI_LAND, comm, mpierr); call CHKERR(mpierr)
      end function
      function mpi_logical_or(comm, lval)
          integer(mpiint),intent(in) :: comm
          logical :: mpi_logical_or
          logical,intent(in) :: lval
          call mpi_allreduce(lval, mpi_logical_or, 1_mpiint, imp_logical, MPI_LOR, comm, mpierr); call CHKERR(mpierr)
      end function

      subroutine imp_allreduce_min(comm, v,r)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(in) :: v
          real(ireal_dp),intent(out) :: r
          call mpi_allreduce(v,r,1,imp_real_dp, MPI_MIN,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_allreduce_max(comm, v,r)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(in) :: v
          real(ireal_dp),intent(out) :: r
          call mpi_allreduce(v,r,1,imp_real_dp, MPI_MAX,comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine imp_reduce_sum(comm, v, myid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(inout) :: v
          integer(mpiint),intent(in) :: myid
          integer(mpiint) :: commsize

          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return 

          if(myid.eq.0) then
            call mpi_reduce(MPI_IN_PLACE, v, 1, imp_real_dp, MPI_SUM, 0, comm, mpierr); call CHKERR(mpierr)
          else
            call mpi_reduce(v, MPI_IN_PLACE, 1, imp_real_dp, MPI_SUM, 0, comm, mpierr); call CHKERR(mpierr)
          endif
      end subroutine

      subroutine  imp_bcast_logical(comm, val, sendid)
          integer(mpiint),intent(in) :: comm
          logical,intent(inout) :: val
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val, 1_mpiint, imp_logical, sendid, comm, mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int(comm, val,sendid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val,1_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_1d(comm, arr,sendid,myid)
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot

          if(sendid.eq.myid) Ntot = size(arr)
          call mpi_bcast(Ntot,1_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_int_2d(comm, arr,sendid,myid)!
          integer(mpiint),intent(in) :: comm
          integer(iintegers),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot(2)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_int,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real(comm, val,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),intent(inout) :: val
          integer(mpiint),intent(in) :: sendid

          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return

          call mpi_bcast(val,1_mpiint,imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_1d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = size(arr)
          call imp_bcast_int(comm, Ntot, sendid)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_2d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(2)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,2_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2)) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_3d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(3)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,3_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3) ) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine
      subroutine  imp_bcast_real_5d(comm, arr,sendid)
          integer(mpiint),intent(in) :: comm
          real(ireal_dp),allocatable,intent(inout) :: arr(:,:,:,:,:)
          integer(mpiint),intent(in) :: sendid
          integer(mpiint) :: myid

          integer(iintegers) :: Ntot(5)
          integer(mpiint) :: commsize
          call MPI_Comm_size( comm, commsize, mpierr); call CHKERR(mpierr)
          if(commsize.le.1) return
          call MPI_Comm_rank( comm, myid, mpierr); call CHKERR(mpierr)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,5_mpiint,imp_int,sendid,comm,mpierr); call CHKERR(mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3), Ntot(4), Ntot(5) ) )
          call mpi_bcast(arr,size(arr),imp_real_dp,sendid,comm,mpierr); call CHKERR(mpierr)
      end subroutine

      elemental subroutine delta_scale( kabs,ksca,g,factor )
          real(ireal_dp),intent(inout) :: kabs,ksca,g ! kabs, ksca, g
          real(ireal_dp),intent(in),optional :: factor
          real(ireal_dp) :: dtau, w0
          dtau = max( kabs+ksca, epsilon(dtau) )
          w0   = ksca/dtau
          g    = g

          if(present(factor)) then
            call delta_scale_optprop( dtau, w0, g, factor)
          else
            call delta_scale_optprop( dtau, w0, g)
          endif

          kabs= dtau * (one-w0)
          ksca= dtau * w0
      end subroutine
      elemental subroutine delta_scale_optprop( dtau, w0, g,factor)
          real(ireal_dp),intent(inout) :: dtau,w0,g
          real(ireal_dp),intent(in),optional :: factor
          real(ireal_dp) :: f

          g = min( g, one-epsilon(g)*10)
          if(present(factor)) then
            f = factor
          else
            f = g**2
          endif
          dtau = dtau * ( one - w0 * f )
          g    = ( g - f ) / ( one - f )
          w0   = w0 * ( one - f ) / ( one - f * w0 )
      end subroutine

      function cumsum(arr)
          real(ireal_dp),intent(in) :: arr(:)
          real(ireal_dp) :: cumsum(size(arr))
          integer :: i
          cumsum(1) = arr(1)
          do i=2,size(arr)
            cumsum(i) = cumsum(i-1) + arr(i)
          enddo
      end function

      pure function compute_normal_3d(p1,p2,p3)
        ! for a triangle p1, p2, p3, if the vector U = p2 - p1 and the vector V = p3 - p1 then the normal
        ! N = U X V and can be calculated by:
        real(ireal_dp), intent(in) :: p1(3), p2(3), p3(3)
        real(ireal_dp) :: compute_normal_3d(3)
        real(ireal_dp) :: U(3), V(3)

        U = p2-p1
        V = p3-p1

        compute_normal_3d(1) = U(2)*V(3) - U(3)*V(2)
        compute_normal_3d(2) = U(3)*V(1) - U(1)*V(3)
        compute_normal_3d(3) = U(1)*V(2) - U(2)*V(1)

        compute_normal_3d = compute_normal_3d / norm(compute_normal_3d)
      end function

    !> @brief For local azimuth and zenith angles, return the local cartesian vectors phi azimuth, theta zenith angles, angles are input in degrees.
    !> @details theta == 0 :: z = -1, i.e. downward
    !> @details azimuth == 0 :: vector going toward minus y, i.e. sun shines from the north
    !> @details azimuth == 90 :: vector going toward minus x, i.e. sun shines from the east
    pure function spherical_2_cartesian(phi, theta, r)
      real(ireal_dp), intent(in) :: phi, theta
      real(ireal_dp), intent(in), optional :: r

      real(ireal_dp) :: spherical_2_cartesian(3)

      spherical_2_cartesian(1) = -sin(deg2rad(theta)) * sin(deg2rad(phi))
      spherical_2_cartesian(2) = -sin(deg2rad(theta)) * cos(deg2rad(phi))
      spherical_2_cartesian(3) = -cos(deg2rad(theta))

      if(present(r)) spherical_2_cartesian = spherical_2_cartesian*r
    end function

    pure function angle_between_two_vec(p1,p2)
      real(ireal_dp),intent(in) :: p1(:), p2(:)
      real(ireal_dp) :: angle_between_two_vec
      real(ireal_dp) :: n1, n2
      n1 = norm(p1)
      n2 = norm(p2)
      angle_between_two_vec = acos(dot_product(p1/norm(p1), p2/norm(p2)))
    end function

      !> @brief determine distance where a photon p intersects with a plane
      !> @details inputs are the location and direction of a photon aswell as the origin and surface normal of the plane
      pure function hit_plane(p_loc, p_dir, po, pn)
        real(ireal_dp) :: hit_plane
        real(ireal_dp),intent(in) :: p_loc(3), p_dir(3)
        real(ireal_dp),intent(in) :: po(3), pn(3)
        real(ireal_dp) :: discr
        discr = dot_product(p_dir,pn)
        if( ( discr.le. epsilon(discr) ) .and. ( discr.gt.-epsilon(discr)  ) ) then
          hit_plane = huge(hit_plane)
        else
          hit_plane = dot_product(po-p_loc, pn) / discr
        endif
      end function

      !> @brief determine if point is inside a triangle p1,p2,p3
      pure function pnt_in_triangle(p1,p2,p3, p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2,p3, p
        logical :: pnt_in_triangle
        real(ireal_dp),parameter :: eps = epsilon(eps), eps2 = 100*eps
        real(ireal_dp) :: a, b, c, edge_dist

        ! First check on rectangular bounding box
        if ( p(1).lt.minval([p1(1),p2(1),p3(1)])-eps2 .or. p(1).gt.maxval([p1(1),p2(1),p3(1)])+eps2 ) then ! outside of xrange
            pnt_in_triangle=.False.
            !print *,'pnt_in_triangle, bounding box check failed:', p
            return
        endif
        if ( p(2).lt.minval([p1(2),p2(2),p3(2)])-eps2 .or. p(2).gt.maxval([p1(2),p2(2),p3(2)])+eps2 ) then ! outside of yrange
            pnt_in_triangle=.False.
            !print *,'pnt_in_triangle, bounding box check failed:', p
            return
        endif

        ! Then check for sides
        a = ((p2(2)- p3(2))*(p(1) - p3(1)) + (p3(1) - p2(1))*(p(2) - p3(2))) / ((p2(2) - p3(2))*(p1(1) - p3(1)) + (p3(1) - p2(1))*(p1(2) - p3(2)))
        b = ((p3(2) - p1(2))*(p(1) - p3(1)) + (p1(1) - p3(1))*(p(2) - p3(2))) / ((p2(2) - p3(2))*(p1(1) - p3(1)) + (p3(1) - p2(1))*(p1(2) - p3(2)))
        c = one - (a + b)

        pnt_in_triangle = all([a,b,c].ge.zero)

        if(.not.pnt_in_triangle) then ! Compute distances to each edge and allow the check to be positive if the distance is small
          edge_dist = distance_to_triangle_edges(p1,p2,p3,p)
          if(edge_dist.le.sqrt(eps)) pnt_in_triangle=.True.
        endif
        !print *,'pnt_in_triangle final:', pnt_in_triangle,'::',a,b,c,':',p,'edgedist',distance_to_triangle_edges(p1,p2,p3,p),distance_to_triangle_edges(p1,p2,p3,p).le.eps
      end function

      pure function distance_to_triangle_edges(p1,p2,p3,p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2,p3, p
        real(ireal_dp) :: distance_to_triangle_edges
        distance_to_triangle_edges = distance_to_edge(p1,p2,p)
        distance_to_triangle_edges = min(distance_to_triangle_edges, distance_to_edge(p2,p3,p))
        distance_to_triangle_edges = min(distance_to_triangle_edges, distance_to_edge(p1,p3,p))
      end function

      pure function distance_to_edge(p1,p2,p)
        real(ireal_dp), intent(in), dimension(2) :: p1,p2, p
        real(ireal_dp) :: distance_to_edge

        distance_to_edge = abs( (p2(2)-p1(2))*p(1) - (p2(1)-p1(1))*p(2) + p2(1)*p1(2) - p2(2)*p1(1) ) / norm(p2-p1)
      end function

      pure function rotate_angle_x(v,angle)
        real(ireal_dp) :: rotate_angle_x(3)
        real(ireal_dp),intent(in) :: v(3), angle
        real(ireal_dp) :: M(3,3),s,c
        s=sin(deg2rad(angle))
        c=cos(deg2rad(angle))

        M(1,:)=[one ,zero ,zero]
        M(2,:)=[zero, c   , s  ]
        M(3,:)=[zero,-s   , c  ]

        rotate_angle_x = matmul(M,v)
      end function
      pure function rotate_angle_y(v,angle)
        real(ireal_dp) :: rotate_angle_y(3)
        real(ireal_dp),intent(in) :: v(3), angle
        real(ireal_dp) :: M(3,3),s,c
        s=sin(deg2rad(angle))
        c=cos(deg2rad(angle))

        M(1,:)=[ c  ,zero , -s ]
        M(2,:)=[zero, one ,zero]
        M(3,:)=[ s  , zero, c  ]

        rotate_angle_y = matmul(M,v)
      end function
      pure function rotate_angle_z(v,angle)
        real(ireal_dp) :: rotate_angle_z(3)
        real(ireal_dp),intent(in) :: v(3), angle
        real(ireal_dp) :: M(3,3),s,c
        s=sin(deg2rad(angle))
        c=cos(deg2rad(angle))

        M(1,:)=[ c  , s   ,zero]
        M(2,:)=[-s  , c   ,zero]
        M(3,:)=[zero, zero, one]

        rotate_angle_z = matmul(M,v)
      end function
    end module
