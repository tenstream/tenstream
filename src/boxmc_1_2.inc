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


subroutine update_dir_stream_1_2(bmc,vertices,p,T)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    type(t_photon),intent(in) :: p
    real(irealbmc),intent(inout) :: T(:)

    select case (p%side)
    case(1:2)
      T(1) = T(1)+p%weight
      return
    case default
      call CHKERR(1_mpiint, 'Dont know what to do with this p%side')
      print *,'remove unused compiler warning', vertices
    end select
    select type(bmc)
    end select
end subroutine
subroutine update_diff_stream_1_2(bmc,vertices,p,S)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    type(t_photon),intent(in) :: p
    real(irealbmc),intent(inout) :: S(:)

    !         _______1_______
    !        |           10  |
    !      5 |            8  | 6
    !        |               |
    !        |   9           |
    !        |   7           |
    !      3 |_______________| 4
    !                2
    if(p%side.eq.1) then
      S(1) = S(1)+p%weight
      return
    else if(p%side.eq.2) then
      S(2) = S(2)+p%weight
      return
    else
      print *,'Couldnt find a stream on to which I can put the photon weight on?!'
      call print_photon(p)
      call exit
      print *,'remove unused compiler warning', vertices
    endif
    select type(bmc)
    end select

end subroutine

subroutine init_dir_photon_1_2(bmc, p, src, ldirect, initial_dir, vertices, ierr)
    class(t_boxmc_1_2) :: bmc
    type(t_photon),intent(inout) :: p
    real(ireal_dp),intent(in) :: initial_dir(:), vertices(:)
    integer(iintegers),intent(in) :: src
    logical,intent(in) :: ldirect
    integer(mpiint),intent(out) :: ierr
    real(ireal_dp) :: dx, dy, dz

    call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

    if(src.eq.1) then
      p%loc = (/L(dx), L(dy), dz/)
    else
      print *,'Dont know what to do with source spec:',src
      call exit
    endif

    p%weight=one
    p%alive = .True.
    p%direct= ldirect
    p%side = int(nil)
    p%src  = src
    p%scattercnt  = 0
    p%dir = initial_dir

    select type(bmc)
    end select

    ierr = 0
end subroutine

subroutine init_diff_photon_1_2(bmc, p, src, vertices, ierr)
    class(t_boxmc_1_2) :: bmc
    type(t_photon),intent(inout) :: p
    real(ireal_dp),intent(in) :: vertices(:)
    integer(iintegers),intent(in) :: src
    integer(mpiint),intent(out) :: ierr

    real(ireal_dp) :: phi,mu,mu1,mu2
    real(ireal_dp) :: dx, dy, dz

    call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

    mu2 = R()
    mu1 = sqrt(R())
    select case(src)
    case (-1) ! thermal emission
      p%loc = (/L(dx), L(dy), L(dz) /)
      mu = 2._ireal_dp*mu2 - one
      phi= deg2rad( R()*360._ireal_dp )
    case (2)
      p%loc = (/L(dx), L(dy),    dz  /)
      mu = -mu1 !mu = interv_R( deg2mu(90._ireal_dp),deg2mu(180._ireal_dp))
      phi= deg2rad( R()*360._ireal_dp )
    case (1)
      p%loc = (/L(dx), L(dy),    zero/)
      mu = mu1 !mu = interv_R( deg2mu(0._ireal_dp),deg2mu(90._ireal_dp))
      phi= deg2rad( R()*360._ireal_dp )
    case default
      print *,'Dont know what to do with source spec:',src
      call exit
    end select

    p%dir = (/cos(phi)*sin(acos(mu)) , sin(phi)*sin(acos(mu)) , mu  /)
    p%weight=one
    p%alive = .True.
    p%direct= .False.
    p%side = int(nil)
    p%src  = src
    p%scattercnt  = 0
    if(.not. approx( norm2(p%dir) ,one ) ) stop 'init_diff_photon_1_2 :: direction was not normed?!'

    select type(bmc)
    end select

    ierr = 0
end subroutine

subroutine intersect_distance_1_2(bmc,vertices,p,max_dist)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    type(t_photon),intent(inout) :: p
    real(ireal_dp),intent(out) :: max_dist

    real(ireal_dp) :: x,y
    real(ireal_dp) :: dx, dy, dz

    dx = 0; dy = 0

    call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

    !crossing with bottom and top plane:
    if(p%dir(3).ge.zero) then
      max_dist = hit_plane(p%loc, p%dir,[zero,zero,dz ],[zero,zero,one])
      p%side=1
      x = p%loc(1)+p%dir(1)*max_dist
      y = p%loc(2)+p%dir(2)*max_dist
      return
    endif
    if(p%dir(3).le.zero) then
      max_dist = hit_plane(p%loc, p%dir,[zero,zero,zero ],[zero,zero,one])
      p%side=2
      x = p%loc(1)+p%dir(1)*max_dist
      y = p%loc(2)+p%dir(2)*max_dist
      return
    endif

    print *,'should actually not be here at the end of crossings in intersect distance!'// &
	    ' - however, please check if distance makes sense?:', &
      max_dist, ':', vertices
    call print_photon(p)

    select type(bmc)
    end select
end subroutine

subroutine box_halfspaces_1_2(bmc, vertices, origins, normals)
    class(t_boxmc_1_2) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    real(ireal_dp),allocatable,intent(out) :: origins(:,:), normals(:,:)
    call box_halfspaces(vertices, origins, normals)
    select type(bmc)
    end select
end subroutine
