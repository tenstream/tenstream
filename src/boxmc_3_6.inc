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
!        5
!   C ------- D
!   |         |
!  4|   1     |3
!   |    2    |
!   |         |
!   A ------- B
!        6
!
! Face numbering is
! 2 -> top
! 1 -> bottom face
! 4 -> side face between (A,C)
! 6 -> side face between (A,B)

subroutine init_dir_photon_3_6(bmc,p,src,ldirect,initial_dir,vertices,ierr)
  class(t_boxmc_3_6) :: bmc
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: ldirect
  integer(mpiint),intent(out) :: ierr

  real(ireal_dp), parameter  :: e=1e-6_ireal_dp
  real(ireal_dp) :: dx, dy, dz

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  p%alive = .False.

  select case(src)
  case (1)
    p%loc = [R()*dx, R()*dy, dz*(one-e)]
  case (2)
    p%loc = [e, R()*dy, R()*dz]
  case (3)
    p%loc = [R()*dx, e, R()*dz]
  case default
    print *,'Dont know what to do with source spec:',src
    call exit
  end select

  p%weight     = one
  p%alive      = .True.
  p%direct     = ldirect
  p%side       = int(nil)
  p%src        = src
  p%scattercnt = 0
  p%dir        = initial_dir
  select type(bmc)
  end select

  ierr = 0
end subroutine

subroutine update_dir_stream_3_6(bmc,vertices,p,T)
  class(t_boxmc_3_6) :: bmc
  real(ireal_dp),intent(in) :: vertices(:)
  type(t_photon),intent(in) :: p
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
     print *,'remove unused compiler warning', vertices
  end select
  select type(bmc)
  end select
end subroutine

subroutine update_diff_stream_3_6(bmc,vertices,p,S)
  class(t_boxmc_3_6) :: bmc
  real(ireal_dp),intent(in) :: vertices(:)
  type(t_photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: S(:)

  select case(p%side)
    case(1)
      S(1) = S(1) + p%weight
    case(2)
      S(2) = S(2) + p%weight
    case(3)
      S(3) = S(3) + p%weight
    case(4)
      S(4) = S(4) + p%weight
    case(5)
      S(5) = S(5) + p%weight
    case(6)
      S(6) = S(6) + p%weight
    case default
      print *,'Couldnt find a stream on to which I can put the photon weight on?!'
      call print_photon(p)
      call exit
      print *,'remove unused compiler warning', vertices
  end select
  select type(bmc)
  end select
end subroutine

subroutine init_diff_photon_3_6(bmc,p,src,vertices,ierr)
  use m_helper_functions_dp, only : pnt_in_triangle,rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, rad2deg
  class(t_boxmc_3_6) :: bmc
  type(t_photon),intent(inout) :: p
  integer(iintegers),intent(in) :: src
  real(ireal_dp),intent(in) :: vertices(:)
  integer(mpiint),intent(out) :: ierr

  real(ireal_dp) :: phi,mu
  real(ireal_dp) :: dx, dy, dz

  !print *,'---------- new diffuse photon --------------'

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  select case(src)
  case (2)
    p%loc = [R()*dx,R()*dy, dz]
  case (1)
    p%loc = [R()*dx, R()*dy, zero]
  case (4)
    p%loc = [zero, R()*dy, R()*dz]
  case (3)
    p%loc = [dx, R()*dy, R()*dz]
  case (6)
    p%loc = [R()*dx, zero, R()*dz]
  case (5)
    p%loc = [R()*dx, dy, R()*dz]
  case default
    print *,'Dont know what to do with source spec:',src
    call exit
  end select

  phi = deg2rad( R()*360 )
  mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
  !mu = R() ! evenly distribution between 0 and 90 degree

  p%dir = [sin(phi)*sin(acos(mu)) , cos(phi)*sin(acos(mu)) , mu]

 select case(src)
  case(1)
    continue

  case(2)
    p%dir = rotate_angle_y(p%dir, 180*one)

  case(3)
    p%dir = rotate_angle_y(p%dir, 90*one)

  case(4)
    p%dir = rotate_angle_y(p%dir, 270*one)

  case(5)
    p%dir = rotate_angle_x(p%dir, 270*one)

  case(6)
    p%dir = rotate_angle_x(p%dir, 90*one)

  case default
    print*,'Dont know what to do with source spec:',src
    call exit
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

  ierr = 0
end subroutine

subroutine intersect_distance_3_6(bmc,vertices,p,max_dist)
  use m_helper_functions, only: cstr
  use m_helper_functions_dp, only: pnt_in_rectangle, angle_between_two_vec, pnt_in_cube
  class(t_boxmc_3_6) :: bmc
  real(ireal_dp),intent(in) :: vertices(:)
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(out) :: max_dist

  if(.not.pnt_in_cube(vertices, p%loc)) then
    print *,'Before intersection:', pnt_in_cube(vertices, p%loc), ':', vertices
    call print_photon(p)
    call CHKERR(1_mpiint, cstr('photon is not inside the local cube','red'))
  endif
  call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    p%side, max_dist)

  !print *,cstr('After intersection:','red'), pnt_in_cube(vertices, p%loc), ':', vertices
  !call print_photon(p)
  select type(bmc)
  end select
end subroutine
