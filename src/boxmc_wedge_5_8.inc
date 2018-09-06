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
function rand_pnt_on_triangle(A,B,C) result(pnt)
  real(ireal_dp), dimension(:), intent(in) :: A, B, C
  real(ireal_dp) :: pnt(size(A))
  real(ireal_dp) :: r1, r2
  r1 = R()
  r2 = R()
  pnt = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
end function

subroutine init_dir_photon_wedge_5_8(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  use m_helper_functions_dp, only: pnt_in_triangle, determine_normal_direction, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_8) :: bmc
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: ldirect
  real(ireal_dp) :: normal(3), triangle_center(3), dir_p_to_triangle_center(3)
  real(ireal_dp),parameter :: edge_start_range=.999_ireal_dp
  integer(mpiint),intent(out) :: ierr

  associate( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18) )

  select case(src)
  case (1)
    p%loc = rand_pnt_on_triangle(pD,pE,pF)
    normal = compute_normal_3d(pD,pF,pE)
  case (2)
    if(R().le..5_ireals) then
      p%loc = rand_pnt_on_triangle(pA,pB,pD)
      normal = compute_normal_3d(pA,pD,pB)
    else
      p%loc = rand_pnt_on_triangle(pB,pD,pE)
      normal = compute_normal_3d(pD,pE,pB)
    endif
  case (3)
    if(R().le..5_ireals) then
      p%loc = rand_pnt_on_triangle(pA,pC,pD)
      normal = compute_normal_3d(pA,pC,pD)
    else
      p%loc = rand_pnt_on_triangle(pC,pD,pF)
      normal = compute_normal_3d(pC,pF,pD)
    endif
  case (4)
    if(R().le..5_ireals) then
      p%loc = rand_pnt_on_triangle(pB,pC,pE)
      normal = compute_normal_3d(pB,pE,pC)
    else
      p%loc = rand_pnt_on_triangle(pC,pE,pF)
      normal = compute_normal_3d(pC,pE,pF)
    endif
  case (5)
    p%loc = rand_pnt_on_triangle(pA,pB,pC)
    normal = compute_normal_3d(pA,pC,pB)
  case default
    call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
    call exit
  end select

  triangle_center = mean([pA,pB,pC,pD,pE,pF])
  end associate

  p%src_side = src

  p%weight=one
  p%alive = .True.
  p%direct= ldirect
  p%side = int(nil)
  p%src  = src
  p%scattercnt = 0
  p%dir = initial_dir

  ! always nudge the photon a little bit inside the volume
  dir_p_to_triangle_center = triangle_center - p%loc
  dir_p_to_triangle_center = dir_p_to_triangle_center / norm(dir_p_to_triangle_center)

  p%loc = p%loc + sqrt(epsilon(pi_dp)) * dir_p_to_triangle_center

  select type(bmc)
  end select

  !print *,'triangle center', triangle_center
  !print *,'normal', normal
  !call print_photon(p)

  if(dot_product(p%dir, normal).lt.zero) then
    ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
    if(ldebug) then
      print *,'UUPS, are we starting in a wrong direction? :: init_dir_photon_wedge_5_8', rad2deg(angle_between_two_vec(p%dir, normal))
      print *,'normal',normal,'direction',p%dir
      call print_photon(p)
    endif
    !stop 'initial direction of the wedge triangle'
    return
  endif
  ierr = 0

end subroutine

subroutine update_dir_stream_wedge_5_8(bmc,vertices,p,T)
  class(t_boxmc_wedge_5_8) :: bmc
  real(ireal_dp),intent(in) :: vertices(:)
  type(t_photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: T(:)

  T(p%side) = T(p%side) + p%weight
  select type(bmc)
  end select
  if(.False.) print *, 'to prevent unused warnings', vertices
end subroutine

subroutine update_diff_stream_wedge_5_8(bmc,p,S)
  class(t_boxmc_wedge_5_8) :: bmc
  type(t_photon),intent(in) :: p
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
    call CHKERR(1_mpiint, 'update_diff_photon_wedge_5_8 called with unknown p%side')
  end select

  select type(bmc)
  end select

  !print *,'updating S',S,'(',p%side,')','+',p%weight
  !call print_photon(p)
  !stop 'debug update_diff_stream_wedge_5_8'
end subroutine

subroutine init_diff_photon_wedge_5_8(bmc, p, src, vertices, ierr)
  use m_helper_functions_dp, only : pnt_in_triangle,rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_8) :: bmc
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: vertices(:)
  integer(iintegers),intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi,mu
  real(ireal_dp),dimension(2) :: A, B, C, nAB, nBC, nCA ! point PNT on triangle [A,B,C]
  real(ireal_dp) :: triangle_center(3), dir_p_to_triangle_center(3)
  real(ireal_dp) :: r1, dx, dy, dz

  associate( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18) )

  select case(src)
  case (1)
    p%loc = rand_pnt_on_triangle(pD,pE,pF)
    p%src_side = 1
  case (2:3)
    if(R().le..5_ireals) then
      p%loc = rand_pnt_on_triangle(pA,pB,pD)
    else
      p%loc = rand_pnt_on_triangle(pB,pD,pE)
    endif
    p%src_side = 2
  case (4:5)
    if(R().le..5_ireals) then
      p%loc = rand_pnt_on_triangle(pA,pC,pD)
    else
      p%loc = rand_pnt_on_triangle(pC,pD,pF)
    endif
    p%src_side = 3
  case (6:7)
    if(R().le..5_ireals) then
      p%loc = rand_pnt_on_triangle(pB,pC,pE)
    else
      p%loc = rand_pnt_on_triangle(pC,pE,pF)
    endif
    p%src_side = 4
  case (8)
    p%loc = rand_pnt_on_triangle(pA,pB,pC)
    p%src_side = 5
  case default
    call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
    call exit
  end select

  triangle_center = mean([pA,pB,pC,pD,pE,pF])
  end associate

  !print *,'---------- new diffuse photon --------------'
  call setup_wedge_coords_from_vertices(vertices, A, B, C, nAB, NBC, nCA, dx, dy, dz)

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
  p%alive = .True.
  p%direct= .False.
  p%side = int(nil)
  p%src  = src
  p%scattercnt = 0

  ! always nudge the photon a little bit inside the volume
  dir_p_to_triangle_center = triangle_center - p%loc
  dir_p_to_triangle_center = dir_p_to_triangle_center / norm(dir_p_to_triangle_center)

  p%loc = p%loc + sqrt(epsilon(pi_dp)) * dir_p_to_triangle_center

  if(.not. approx( norm(p%dir) ,one ) ) call CHKERR(1_mpiint, 'init_diff_photon_3_10 :: direction was not normed?!')

  ierr = 0
  !call print_photon(p)
  !print *,'--------- end new diffuse photon ----------'
  !stop 'not yet implemented :: init_diff_photon_wedge_5_8'
  select type(bmc)
  end select
end subroutine

subroutine intersect_distance_wedge_5_8(bmc, vertices, p, max_dist)
  use m_boxmc_geometry, only : intersect_wedge
  class(t_boxmc_wedge_5_8) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist

  call intersect_wedge(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    p%side, p%weight, max_dist)

  select type(bmc)
  end select
end subroutine
! ******************* BOXMC wedge 5 8 ************************

