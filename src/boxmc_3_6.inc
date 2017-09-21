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

! Distribution Code for direct streams/targets:
!
!   C ------- D
!   |         |
!  2|   1     |
!   |         |
!   |         |
!   A ------- B
!        3
! Labeling for destinations is symmetric horizontally and vertically
!
! Face numbering is
! 1 -> top/bot face
! 2 -> side face between (A,C)/(A,D)
! 3 -> side face between (A,B)/(C,D)

! Distribution Code for diffuse streams/targets:
!        6
!   C ------- D
!   |         |
!  3|   1     |4
!   |    2    |
!   |         |
!   A ------- B
!        5
!
! Face numbering is
! 1 -> top
! 2 -> bottom face
! 3 -> side face between (A,C)
! 5 -> side face between (A,B)

subroutine init_dir_photon_3_6(bmc,p,src,direct,initial_dir,dx,dy,dz)
  class(t_boxmc_3_6) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: dx,dy,dz,initial_dir(3)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: direct

  select case(src)
  case (1)
    p%loc = [R()*dx,R()*dy, dz]
  case (2)
    p%loc = [zero, R()*dy, R()*dz]
  case (3)
    p%loc = [R()*dx, zero, R()*dz]
  case default
    print *,'Dont know what to do with source spec:',src
    call exit
  end select

  p%weight=one
  p%dx   = dx
  p%dy   = dy
  p%dz   = dz
  p%alive = .True.
  p%direct= direct
  p%side = int(nil)
  p%src  = src
  p%scattercnt = 0
  p%dir = initial_dir
  select type(bmc)
  end select
end subroutine

subroutine update_dir_stream_3_6(bmc,p,T)
  class(t_boxmc_3_6) :: bmc
  type(photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: T(:)
  select case(p%side)
   case(1,2)
     T(1)=T(1)+p%weight
   case(3,4)
     T(2)=T(2)+p%weight
   case(5,6)
     T(3)=T(3)+p%weight
   case default
     stop 'called update dirstream with incorrect source'
  end select
  select type(bmc)
  end select
end subroutine

subroutine update_diff_stream_3_6(bmc,p,S)
  class(t_boxmc_3_6) :: bmc
  type(photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: S(:)

  S(p%side) = S(p%side) + p%weight
  select type(bmc)
  end select
end subroutine

subroutine init_diff_photon_3_6(bmc,p,src,dx,dy,dz)
  use m_helper_functions_dp, only : pnt_in_triangle,rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, rad2deg
  class(t_boxmc_3_6) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: dx,dy,dz
  integer(iintegers),intent(in) :: src

  real(ireal_dp) :: phi,mu

  !print *,'---------- new diffuse photon --------------'

  select case(src)
  case (1)
    p%loc = [R()*dx,R()*dy, dz]
  case (2)
    p%loc = [R()*dx, R()*dy, zero]
  case (3)
    p%loc = [zero, R()*dy, R()*dz]
  case (4)
    p%loc = [dx, R()*dy, R()*dz]
  case (5)
    p%loc = [R()*dx, zero, R()*dz]
  case (6)
    p%loc = [R()*dx, dy, R()*dz]
  case default
    print *,'Dont know what to do with source spec:',src
    call exit
  end select

  phi = deg2rad( R()*360 )
  mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
  !mu = R() ! evenly distribution between 0 and 90 degree

  p%dir = [cos(phi)*sin(acos(mu)) , sin(phi)*sin(acos(mu)) , mu]

 select case(src)
  case(1)
    p%dir = rotate_angle_y(p%dir, 180*one)

  case(2)
    continue

  case(3)
    p%dir = rotate_angle_y(p%dir, 270*one)

  case(4)
    p%dir = rotate_angle_y(p%dir, 90*one)

  case(5)
    p%dir = rotate_angle_x(p%dir, 90*one)

  case(6)
    p%dir = rotate_angle_x(p%dir, 270*one)

  case default
    print*,'Sont knwo what to do with source spec:',src
    call exit
  end select

  p%weight=one
  p%dx   = dx
  p%dy   = dy
  p%dz   = dz
  p%alive = .True.
  p%direct= .False.
  p%side = int(nil)
  p%src  = src
  p%scattercnt = 0
  if(.not. approx( norm(p%dir) ,one ) ) stop 'init_diff_photon_3_10 :: direction was not normed?!'

  select type(bmc)
  end select
end subroutine

subroutine intersect_distance_3_6(bmc,p,max_dist)
  use m_helper_functions_dp, only: pnt_in_rectangle, distance_to_edge, distance_to_triangle_edges
  class(t_boxmc_3_6) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(out) :: max_dist

  integer(iintegers) :: i

  real(ireal_dp) :: dist(6) ! distances to all the faces
  real(ireal_dp) :: new_loc(3)

  logical :: l_in_rectangle(2)

  dist = huge(dist)
  !crossing with bottom and top plane:
  if(p%dir(3).ge.zero) then
    dist(1) = hit_plane(p%loc, p%dir,[zero, zero ,p%dz ],[zero,zero,one])
  else if(p%dir(3).le.zero) then
    dist(2) = hit_plane(p%loc, p%dir,[zero, zero, zero ],[zero,zero,one])
  endif

  !crossing with side planes
  if(p%dir(1).le.zero) then
    dist(3) = hit_plane(p%loc, p%dir, [zero, zero, zero], [one, zero, zero])
  else if(p%dir(1).ge.zero) then
    dist(4) = hit_plane(p%loc, p%dir, [p%dx, zero, zero], [one, zero, zero])
  endif

  !crossing with other side planes
  if(p%dir(2).le.zero) then
    dist(5) = hit_plane(p%loc, p%dir, [zero, zero, zero], [zero, one, zero])
  else if(p%dir(2).ge.zero) then
    dist(6) = hit_plane(p%loc, p%dir, [zero, p%dy, zero], [zero, one, zero])
  endif

  max_dist = huge(dist)
  p%side = -i1
  do i=1,size(dist)
    if(dist(i).lt.zero) cycle ! dont travel backwards

    if(dist(i).le.max_dist) then
      max_dist = dist(i)
      p%side = i
    endif
  enddo

  if(p%side.eq.-i1) then
    print *,'should actually not be here at the end of crossings in intersect distance! - however, please check if distance makes sense?:',max_dist
    call print_photon(p)
    stop 'debug'
  endif

  new_loc = p%loc + max_dist * p%dir
  l_in_rectangle(1) = pnt_in_rectangle([zero,zero], [p%dx,zero], [zero,p%dy], [new_loc(1), new_loc(2)])
  l_in_rectangle(2) = pnt_in_rectangle([zero,zero], [zero,p%dz], [p%dx,p%dz], [new_loc(1), new_loc(3)])
  if(any(.not.l_in_rectangle)) then
    print *,'photon not in rectangle', l_in_rectangle
    print *,'dists to faces', dist
    print *,'new point not in rectangle', new_loc, 'side', p%side, 'dist', max_dist
    call print_photon(p)
    stop 'debug'
  endif

  select type(bmc)
  end select
end subroutine
