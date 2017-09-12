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

! Distribution Code for Wedges:
!
!        C
!       / \
!     3/   \4
!     /  1  \
!    /    5  \
!   A ------- B
!        2
!
! Face numbering is
! 1 -> top
! 2 -> side face between (A,B)
! 3 -> side face between (A,C)
! 5 -> bottom face
!
! We always assume the triangle to have dx edge length along y=0 (A,B) and dy is the edge length between (A,C)
! Distribute Photons on triangles: https://doi.org/10.1145/571647.571648

subroutine define_wedge_triangle_pts(dx,dy,A,B,C)
  real(ireal_dp),intent(in) :: dx,dy
  real(ireal_dp),dimension(2),intent(out) :: A, B, C ! points on triangle [A,B,C]

  A = [zero, zero]
  B = [dx, zero]
  C = [dx/2, sqrt(dy**2 - (dx/2)**2)]
end subroutine

subroutine init_dir_photon_wedge_5_5(bmc,p,src,direct,initial_dir,dx,dy,dz)
  use m_helper_functions_dp, only: pnt_in_triangle, determine_normal_direction
  class(t_boxmc_wedge_5_5) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: dx,dy,dz,initial_dir(3)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: direct
  real(ireal_dp),dimension(2) :: A, B, C, PNT, CENTER ! point PNT on triangle [A,B,C]
  real(ireal_dp) :: r1, r2

  call define_wedge_triangle_pts(dx,dy,A,B,C)
  CENTER = (A+B+C)/3

  select case(src)
  case (1)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), dz]
  case (2)
    PNT = A + (R() * (B-A))
    p%loc = [PNT(1), PNT(2), L(dz)]
  case (3)
    PNT = A + (R() * (C-A))
    p%loc = [PNT(1), PNT(2), L(dz)]
  case (4)
    PNT = C + (R() * (B-C))
    p%loc = [PNT(1), PNT(2), L(dz)]

  case (5)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), zero]

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

  if(.not. pnt_in_triangle(A,B,C, PNT)) then
    print *,'UUPS, the initial location of the photon seems to be outside of the wedge triangle :: init_dir_photon_wedge_5_5'
    call print_photon(p)
    stop 'initial location outside of the wedge triangle'
  endif

end subroutine

subroutine update_dir_stream_wedge_5_5(bmc,p,T)
  class(t_boxmc_wedge_5_5) :: bmc
  type(photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: T(:)

  T(p%side) = T(p%side) + p%weight
  select type(bmc)
  end select
end subroutine

subroutine update_diff_stream_wedge_5_5(bmc,p,S)
  class(t_boxmc_wedge_5_5) :: bmc
  type(photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: S(:)

  S(p%side) = S(p%side) + p%weight
  select type(bmc)
  end select

  !print *,'updating S',S,'(',p%side,')','+',p%weight
  !call print_photon(p)
end subroutine

subroutine init_diff_photon_wedge_5_5(bmc,p,src,dx,dy,dz)
  use m_helper_functions_dp, only : pnt_in_triangle,rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_5) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: dx,dy,dz
  integer(iintegers),intent(in) :: src

  real(ireal_dp) :: phi,mu
  real(ireal_dp),dimension(2) :: A, B, C, PNT ! point PNT on triangle [A,B,C]
  real(ireal_dp) :: r1, r2

  !print *,'---------- new diffuse photon --------------'
  call define_wedge_triangle_pts(dx,dy,A,B,C)

  select case (src)
  case (1)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), dz]
  case (2)
    PNT = A + (B - A) * R()
    p%loc = [PNT(1), PNT(2), R()*dz]
  case (3)
    PNT = A + (C - A) * R()
    p%loc = [PNT(1), PNT(2), R()*dz]
  case (4)
    PNT = C + (B - C) * R()
    p%loc = [PNT(1), PNT(2), R()*dz]
  case (5)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), zero]
  case default
    stop 'init_diff_photon_wedge_5_5 called with unknown src'
  end select

  phi = deg2rad( R()*360 )
  ! mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
  mu = R() ! evenly distribution between 0 and 90 degree

  p%dir = (/cos(phi)*sin(acos(mu)) , sin(phi)*sin(acos(mu)) , mu  /)

  select case(src)
  case(1)
    p%dir = rotate_angle_y(p%dir,180*one) ! flip from upward to downward direction

  case(2)
    p%dir = rotate_angle_x(p%dir, 90*one)
    ! rotate towards horizon, pointing in +y and tilt the angle a bit more in accordance to plane defined from (AB).
    ! the trick with 90 deg back helps to retain direction because
    ! arccos does not distinguish between postitive or negative angles
    r1 = rad2deg(angle_between_two_vec([zero, one], B-A)) - 90
    p%dir = rotate_angle_z(p%dir, r1)

  case(3)
    p%dir = rotate_angle_x(p%dir, 90*one)

    r1 = rad2deg(angle_between_two_vec([zero, one], C-A)) + 90
    p%dir = rotate_angle_z(p%dir, r1)

  case(4)
    p%dir = rotate_angle_x(p%dir, 90*one)

    r1 = rad2deg(angle_between_two_vec([zero, one], B-C)) + 90
    p%dir = rotate_angle_z(p%dir, r1)

  case (5)
    continue ! already in correct direction
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

  !call print_photon(p)

  if(.not. pnt_in_triangle(A,B,C, PNT)) then
    print *,'UUPS, the initial location of the photon seems to be outside of the wedge triangle :: init_diff_photon_wedge_5_5'
    call print_photon(p)
    stop 'initial location outside of the wedge triangle'
  endif
  !print *,'--------- end new diffuse photon ----------'
  !stop 'not yet implemented :: init_diff_photon_wedge_5_5'
end subroutine

subroutine intersect_distance_wedge_5_5(bmc,p,max_dist)
  use m_helper_functions_dp, only: pnt_in_triangle, distance_to_edge, distance_to_triangle_edges, determine_normal_direction, angle_between_two_vec
  class(t_boxmc_wedge_5_5) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(out) :: max_dist

  integer(iintegers) :: i

  real(ireal_dp) :: dist(5) ! distances to all the faces
  real(ireal_dp) :: normals(3,5), new_loc(3)
  real(ireal_dp),dimension(2) :: A, B, C, CENTER ! points on triangle [A,B,C]

  logical :: l_in_triangle

  call define_wedge_triangle_pts(p%dx,p%dy,A,B,C)

  CENTER = (A+B+C)/3

  normals(:,1) = [zero,zero,-one]
  normals(:,2) = compute_normal_3d([A(1),A(2),zero], [A(1),A(2),p%dz], [B(1),B(2),zero])
  normals(:,3) = compute_normal_3d([A(1),A(2),zero], [A(1),A(2),p%dz], [C(1),C(2),zero])
  normals(:,4) = compute_normal_3d([C(1),C(2),zero], [C(1),C(2),p%dz], [B(1),B(2),zero])
  normals(:,5) = [zero,zero, one]

  i=2; normals(:,i) = normals(:,i) * determine_normal_direction(normals(1:2,i), A+(B-A)/2, CENTER)
  i=3; normals(:,i) = normals(:,i) * determine_normal_direction(normals(1:2,i), A+(C-A)/2, CENTER)
  i=4; normals(:,i) = normals(:,i) * determine_normal_direction(normals(1:2,i), C+(B-C)/2, CENTER)

  !crossing with bottom and top plane:
  dist(1) = hit_plane(p%loc, p%dir,[A(1),A(2),p%dz ],normals(:,1))
  dist(5) = hit_plane(p%loc, p%dir,[A(1),A(2),zero ],normals(:,5))

  !crossing with side planes:
  ! plane 2, along y=0
  dist(2) = hit_plane(p%loc, p%dir, [A(1), A(2), zero], normals(:,2))
  dist(3) = hit_plane(p%loc, p%dir, [C(1), C(2), zero], normals(:,3))
  dist(4) = hit_plane(p%loc, p%dir, [C(1), C(2), zero], normals(:,4))

  max_dist = huge(dist)
  p%side = -i1
  do i=1,size(dist)
    if(p%scattercnt.eq.0 .and. i.eq.p%src) cycle ! we just started this, dont allow that it hits the start face

    if(dist(i).ge.-sqrt(epsilon(zero)) .and. dist(i).le.zero) then ! have photon that is very close to the edge, this could be a numerics issue... lets check in which direction the photon travels
      if( angle_between_two_vec(p%dir, -normals(:,i)) .le. pi_dp/2 ) then ! photon is flying towards this edge but intersection code said that the crossing would be behind us, this can only mean that we are actually already outside of the triangle...
        dist(i) = zero
      endif
    endif

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
  l_in_triangle = pnt_in_triangle(A,B,C, [new_loc(1), new_loc(2)])
  if(.not.l_in_triangle) then
    print *,'distance:', distance_to_triangle_edges(A,B,C, [new_loc(1), new_loc(2)]), ':: dists to faces',dist
    print *,'called pnt_in_triangle(', A, B, C, [new_loc(1), new_loc(2)], ')'
    print *,'new point not in triangle', new_loc, 'side', p%side, 'dist', max_dist
    call print_photon(p)
    stop 'debug'
  endif


  select type(bmc)
  end select

end subroutine
! ******************* BOXMC wedge 3 5 ************************

