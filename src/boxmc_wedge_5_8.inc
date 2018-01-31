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

!defined in wedge_5_5 subroutine define_wedge_triangle_pts(dx,dy,A,B,C, nAB, nBC, nCA)

subroutine init_dir_photon_wedge_5_8(bmc,p,src,direct,initial_dir,dx,dy,dz)
  use m_helper_functions_dp, only: pnt_in_triangle, determine_normal_direction, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_8) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: dx,dy,dz,initial_dir(3)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: direct
  real(ireal_dp),dimension(2) :: A, B, C, PNT, nAB, nBC, nCA ! point PNT on triangle [A,B,C]
  real(ireal_dp) :: r1, r2, normal(3), triangle_center(3), dir_p_to_triangle_center(3)

  call define_wedge_triangle_pts(dx,dy,A,B,C, nAB, NBC, nCA)

  select case(src)
  case (1)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), dz]
    normal = [zero,zero,-one]
  case (2)
    PNT = A + (R() * (B-A))
    p%loc = [PNT(1), PNT(2), L(dz)]
    normal = [nAB(1), nAB(2), zero]
  case (3)
    PNT = A + (R() * (C-A))
    p%loc = [PNT(1), PNT(2), L(dz)]
    normal = [nCA(1), nCA(2), zero]
  case (4)
    PNT = C + (R() * (B-C))
    p%loc = [PNT(1), PNT(2), L(dz)]
    normal = [nBC(1), nBC(2), zero]
  case (5)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), zero]
    normal = [zero,zero,one]

  case default
    print *,'Dont know what to do with source spec:',src
    call exit
  end select

  p%src_side = src

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

  if(.not. pnt_in_triangle(A,B,C, p%loc)) then
    print *,'UUPS, the initial location of the photon seems to be outside of the wedge triangle :: init_dir_photon_wedge_5_8'
    call print_photon(p)
    stop 'initial location outside of the wedge triangle'
  endif

  if(angle_between_two_vec(p%dir, normal).ge.pi_dp/2-sqrt(epsilon(pi_dp))) then
    if(angle_between_two_vec(p%dir, normal).le.pi_dp/2+sqrt(epsilon(pi_dp))) then
      triangle_center(1:2) = (A+B+C)/3
      triangle_center(3) = dz/2
      dir_p_to_triangle_center = triangle_center - p%loc
      dir_p_to_triangle_center = dir_p_to_triangle_center / norm(dir_p_to_triangle_center)
      !print *,'the difference is really small so I will just move the photon a bit towards the center of the triangle and continue'
      !print *,'triangle_center', triangle_center, 'ploc', p%loc(1:2)
      !print *,'dir',  dir_p_to_triangle_center,'::',normal, sqrt(epsilon(pi_dp)) * dir_p_to_triangle_center

      p%loc = p%loc + sqrt(epsilon(pi_dp)) * dir_p_to_triangle_center
    else
      print *,'UUPS, are we starting in a wrong direction? :: init_dir_photon_wedge_5_8', rad2deg(angle_between_two_vec(p%dir, normal))
      print *,'normal',normal,'direction',p%dir
      call print_photon(p)
      stop 'initial direction of the wedge triangle'
    endif
  endif
end subroutine

subroutine update_dir_stream_wedge_5_8(bmc,p,T)
  class(t_boxmc_wedge_5_8) :: bmc
  type(photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: T(:)

  T(p%side) = T(p%side) + p%weight
  select type(bmc)
  end select
end subroutine

subroutine update_diff_stream_wedge_5_8(bmc,p,S)
  class(t_boxmc_wedge_5_8) :: bmc
  type(photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: S(:)

  select case (p%side)
  case (1)
    S(1) = S(1)+p%weight
  case (2)
    if(p%dir(3).le.zero) then
      S(2) = S(2)+p%weight
    else
      S(3) = S(3)+p%weight
    endif
  case (3)
    if(p%dir(3).le.zero) then
      S(4) = S(4)+p%weight
    else
      S(5) = S(5)+p%weight
    endif
  case (4)
    if(p%dir(3).le.zero) then
      S(6) = S(6)+p%weight
    else
      S(7) = S(7)+p%weight
    endif
  case (5)
    S(8) = S(8)+p%weight
  case default
    call print_photon(p)
    stop 'update_diff_photon_wedge_5_8 called with unknown p%side'
  end select

  select type(bmc)
  end select

  !print *,'updating S',S,'(',p%side,')','+',p%weight
  !call print_photon(p)
  !stop 'debug update_diff_stream_wedge_5_8'
end subroutine

subroutine init_diff_photon_wedge_5_8(bmc,p,src,dx,dy,dz)
  use m_helper_functions_dp, only : pnt_in_triangle,rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_8) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: dx,dy,dz
  integer(iintegers),intent(in) :: src

  real(ireal_dp) :: phi,mu
  real(ireal_dp),dimension(2) :: A, B, C, PNT, nAB, NBC, nCA ! point PNT on triangle [A,B,C]
  real(ireal_dp) :: r1, r2, intersec_dist

  !print *,'---------- new diffuse photon --------------'
  call define_wedge_triangle_pts(dx,dy,A,B,C, nAB, NBC, nCA)

  select case (src)
  case (1)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), dz]
    p%src_side = 1
  case (2:3)
    PNT = A + (B - A) * R()
    p%loc = [PNT(1), PNT(2), R()*dz]
    p%src_side = 2
  case (4:5)
    PNT = A + (C - A) * R()
    p%loc = [PNT(1), PNT(2), R()*dz]
    p%src_side = 3
  case (6:7)
    PNT = C + (B - C) * R()
    p%loc = [PNT(1), PNT(2), R()*dz]
    p%src_side = 4
  case (8)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), zero]
    p%src_side = 5
  case default
    stop 'init_diff_photon_wedge_5_8 called with unknown src'
  end select

  select case(src)
  case(1,8)
    phi= deg2rad( R()*360 )
  case(2:7)
    phi= deg2rad( interv_R(0._ireal_dp, 180._ireal_dp) ) ! from +x to -x through +y
  end select

  ! mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
  mu = R() ! evenly distribution between 0 and 90 degree

  p%dir = (/cos(phi)*sin(acos(mu)) , sin(phi)*sin(acos(mu)) , mu  /)

  select case(src)
  case(1)
    p%dir = rotate_angle_y(p%dir,180*one) ! flip from upward to downward direction

  case(2:3)
    if(src.eq.3) p%dir = rotate_angle_z(p%dir, 180*one) ! rotate by 180 so that phi looks towards the other direction
    p%dir = rotate_angle_x(p%dir, 90*one)
    ! rotate towards horizon, pointing in +y and tilt the angle a bit more in accordance to plane defined from (AB).
    ! the trick with 90 deg back helps to retain direction because
    ! arccos does not distinguish between postitive or negative angles
    r1 = rad2deg(angle_between_two_vec([zero, one], B-A)) - 90
    p%dir = rotate_angle_z(p%dir, r1)

  case(4:5)
    if(src.eq.5) p%dir = rotate_angle_z(p%dir, 180*one) ! rotate by 180 so that phi looks towards the other direction
    p%dir = rotate_angle_x(p%dir, 90*one)
    r1 = rad2deg(angle_between_two_vec([zero, one], C-A)) + 90
    p%dir = rotate_angle_z(p%dir, r1)

  case(6:7)
    if(src.eq.7) p%dir = rotate_angle_z(p%dir, 180*one) ! rotate by 180 so that phi looks towards the other direction
    p%dir = rotate_angle_x(p%dir, 90*one)
    r1 = rad2deg(angle_between_two_vec([zero, one], B-C)) + 90
    p%dir = rotate_angle_z(p%dir, r1)

  case (8)
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

  if(.not. pnt_in_triangle(A,B,C, PNT)) then
    print *,'UUPS, the initial location of the photon seems to be outside of the wedge triangle :: init_diff_photon_wedge_5_8'
    call print_photon(p)
    stop 'initial location outside of the wedge triangle'
  endif


  !call print_photon(p)
  !print *,'--------- end new diffuse photon ----------'
  !stop 'not yet implemented :: init_diff_photon_wedge_5_8'
end subroutine

subroutine intersect_distance_wedge_5_8(bmc,p,max_dist)
  use m_helper_functions_dp, only: pnt_in_triangle, distance_to_edge, &
    determine_normal_direction, angle_between_two_vec, distances_to_triangle_edges
  class(t_boxmc_wedge_5_8) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(out) :: max_dist

  real(ireal_dp),dimension(2) :: nAB, nBC, nCA
  real(ireal_dp),dimension(2) :: A, B, C ! points on triangle [A,B,C]
  real(ireals),dimension(3) :: Ab, Bb, Cb ! points on triangle [A,B,C] at bottom
  real(ireals),dimension(3) :: At, Bt, Ct ! points on triangle [A,B,C] at top

  logical :: l_in_triangle
  logical :: lhit
  real(ireals) :: hit(4)

  call define_wedge_triangle_pts(p%dx, p%dy, A, B, C, nAB, NBC, nCA)

  Ab = [A(1), A(2), zero]
  Bb = [B(1), B(2), zero]
  Cb = [C(1), C(2), zero]

  At = [A(1), A(2), p%dz]
  Bt = [B(1), B(2), p%dz]
  Ct = [C(1), C(2), p%dz]

  !crossing with bottom and top plane:
  if(p%dir(3).gt.zero) then
    if(p%scattercnt.gt.0 .or. p%src_side.ne.1) then
      call triangle_intersection(p%loc, p%dir, At, Bt, Ct, lhit, hit)
      if(lhit) p%side = 1
    endif
  else if (p%dir(3).lt.zero) then
    if(p%scattercnt.gt.0 .or. p%src_side.ne.5) then
      call triangle_intersection(p%loc, p%dir, Ab, Bb, Cb, lhit, hit)
      if(lhit) p%side = 5
    endif
  endif

  !crossing with side planes:
  ! plane 2, along y=0
  if(.not.lhit) then
    if(p%scattercnt.gt.0 .or. p%src_side.ne.2) then
      call square_intersection(p%loc, p%dir, Ab, Bb, Bt, At, lhit, hit)
      if(lhit) p%side = 2
    endif
  endif
  if(.not.lhit) then
    if(p%scattercnt.gt.0 .or. p%src_side.ne.3) then
      call square_intersection(p%loc, p%dir, Ab, Cb, Ct, At, lhit, hit)
      if(lhit) p%side = 3
    endif
  endif
  if(.not.lhit) then
    if(p%scattercnt.gt.0 .or. p%src_side.ne.4) then
      call square_intersection(p%loc, p%dir, Bb, Cb, Ct, Bt, lhit, hit)
      if(lhit) p%side = 4
    endif
  endif

  if(.not.lhit) then
    print *,'should actually not be here at the end of crossings in intersect distance!'
    call print_photon(p)
    stop 'ERROR in Raytracer, didnt hit anything!'
  endif

  l_in_triangle = pnt_in_triangle(A,B,C, hit(1:2))
  if(.not.l_in_triangle) then
    print *,'I hit something...', p%side, '::', hit
    print *,'distance to edges:', distances_to_triangle_edges(A,B,C, hit(1:2)), ':: dists to faces',hit(4)
    print *,'called pnt_in_triangle(', A, B, C, hit(1:2), ')'
    print *,'target point not in triangle', hit, 'side', p%side, 'dist', hit(4)
    call print_photon(p)
    stop 'debug'
  endif

  select type(bmc)
  end select

  max_dist = hit(4)
end subroutine
! ******************* BOXMC wedge 5 5 ************************

