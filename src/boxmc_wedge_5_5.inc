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

subroutine init_dir_photon_wedge_5_5(bmc,p,src,ldirect,initial_dir,vertices,ierr)
  use m_helper_functions_dp, only: pnt_in_triangle, determine_normal_direction, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_5) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: ldirect
  integer(mpiint),intent(out) :: ierr

  real(ireal_dp),dimension(2) :: A, B, C, PNT, nAB, nBC, nCA ! point PNT on triangle [A,B,C]
  real(ireal_dp) :: r1, r2, normal(3), triangle_center(3), dir_p_to_triangle_center(3)

  call setup_wedge_coords_from_vertices(vertices, A, B, C, nAB, NBC, nCA, p%dx, p%dy, p%dz)

  select case(src)
  case (1)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), p%dz]
    normal = [zero,zero,-one]
  case (2)
    PNT = A + (R() * (B-A))
    p%loc = [PNT(1), PNT(2), L(p%dz)]
    normal = [nAB(1), nAB(2), zero]
  case (3)
    PNT = A + (R() * (C-A))
    p%loc = [PNT(1), PNT(2), L(p%dz)]
    normal = [nCA(1), nCA(2), zero]
  case (4)
    PNT = C + (R() * (B-C))
    p%loc = [PNT(1), PNT(2), L(p%dz)]
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

  p%weight=one
  p%alive = .True.
  p%direct= ldirect
  p%side = int(nil)
  p%src  = src
  p%scattercnt = 0
  p%dir = initial_dir
  select type(bmc)
  end select

  if(.not. pnt_in_triangle(A,B,C, p%loc)) then
    print *,'UUPS, the initial location of the photon seems to be outside of the wedge triangle :: init_dir_photon_wedge_5_5'
    call print_photon(p)
    stop 'initial location outside of the wedge triangle'
  endif

  if(angle_between_two_vec(p%dir, normal).gt.pi_dp/2-sqrt(epsilon(pi_dp))) then
    call print_photon(p)
    if(angle_between_two_vec(p%dir, normal).le.pi_dp/2+10*epsilon(pi_dp)) then
      !print *,'the difference is really small so I will just move the photon a bit towards the center of the triangle and continue'
      triangle_center(1:2) = (A+B+C)/3
      triangle_center(3) = p%dz/2
      dir_p_to_triangle_center = triangle_center - p%loc
      dir_p_to_triangle_center = dir_p_to_triangle_center / norm(dir_p_to_triangle_center)
      p%loc = p%loc + sqrt(epsilon(pi_dp)) * dir_p_to_triangle_center
    else
      ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
      if(ldebug) print *,'UUPS, are we starting in a wrong direction?:: init_dir_photon_wedge_5_5', rad2deg(angle_between_two_vec(p%dir, normal))
      if(ldebug) print *,'normal',normal,'direction',p%dir
      return
      !stop 'initial direction of the wedge triangle'
    endif
  endif
  ierr = 0
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

subroutine init_diff_photon_wedge_5_5(bmc, p, src, vertices, ierr)
  use m_helper_functions_dp, only : pnt_in_triangle,rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_5) :: bmc
  type(photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: vertices(:)
  integer(iintegers),intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi,mu
  real(ireal_dp),dimension(2) :: A, B, C, PNT, nAB, NBC, nCA ! point PNT on triangle [A,B,C]
  real(ireal_dp) :: r1, r2

  !print *,'---------- new diffuse photon --------------'
  call setup_wedge_coords_from_vertices(vertices, A, B, C, nAB, NBC, nCA, p%dx, p%dy, p%dz)

  select case (src)
  case (1)
    r1 = R()
    r2 = R()
    PNT = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    p%loc = [PNT(1), PNT(2), p%dz]
  case (2)
    PNT = A + (B - A) * R()
    p%loc = [PNT(1), PNT(2), R()*p%dz]
  case (3)
    PNT = A + (C - A) * R()
    p%loc = [PNT(1), PNT(2), R()*p%dz]
  case (4)
    PNT = C + (B - C) * R()
    p%loc = [PNT(1), PNT(2), R()*p%dz]
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
  ierr = 0
end subroutine

subroutine intersect_distance_wedge_5_5(bmc,vertices,p,max_dist)
  use m_helper_functions_dp, only: pnt_in_triangle, distance_to_edge, distances_to_triangle_edges, determine_normal_direction, angle_between_two_vec
  use m_boxmc_geometry, only : intersect_wedge
  class(t_boxmc_wedge_5_5) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist


  call intersect_wedge(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    p%side, p%weight, max_dist)

  select type(bmc)
  end select

end subroutine
! ******************* BOXMC wedge 5 5 ************************

