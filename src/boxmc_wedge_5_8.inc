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
subroutine rand_pnt_on_plane(A,B,C,D, pnt, normal, U,V)
  use m_helper_functions_dp, only: cross_3d
  real(ireal_dp), dimension(:), intent(in) :: A, B, C, D
  real(ireal_dp), dimension(size(A)), intent(out) :: pnt, normal, U, V
  real(ireal_dp) :: area(2)

  area(1) = triangle_area_by_vertices(A,B,C)
  area(2) = triangle_area_by_vertices(A,C,D)

  if(R().lt.area(1)/sum(area)) then
    pnt = rand_pnt_on_triangle(A,B,C)
    normal = compute_normal_3d(A,B,C)
    U = (C-B)/norm(C-B)
    V = -cross_3d(U,normal)
  else
    pnt = rand_pnt_on_triangle(A,C,D)
    normal = compute_normal_3d(A,C,D)
    U = (D-A)/norm(D-A)
    V = -cross_3d(U,normal)
  endif
end subroutine

subroutine init_dir_photon_wedge_5_8(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  use m_helper_functions_dp, only: pnt_in_triangle, determine_normal_direction, angle_between_two_vec, rad2deg
  class(t_boxmc_wedge_5_8) :: bmc
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: ldirect
  real(ireal_dp), dimension(3) :: normal, triangle_center, dir_p_to_triangle_center, U, V
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
    call rand_pnt_on_plane(pA,pD,pE,pB, p%loc, normal, U, V)
  case (3)
    call rand_pnt_on_plane(pA,pC,pF,pD, p%loc, normal, U, V)
  case (4)
    call rand_pnt_on_plane(pB,pE,pF,pC, p%loc, normal, U, V)
  case (5)
    p%loc = rand_pnt_on_triangle(pA,pB,pC)
    normal = compute_normal_3d(pA,pB,pC)
  case default
    call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
    call exit
  end select

  triangle_center = mean([pA,pB,pC,pD,pE,pF])

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

  end associate

  select type(bmc)
  end select

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

subroutine update_diff_stream_wedge_5_8(bmc,vertices,p,S)
  class(t_boxmc_wedge_5_8) :: bmc
  real(ireal_dp),intent(in) :: vertices(:)
  type(t_photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: S(:)
  real(ireal_dp) :: normal(3)

  ! in case of spherically distorted elements,
  ! we should probably use the tilting angle of the side planes
  ! to separate streams instead of z-component of direction

  associate( &
      pA => vertices(1:3), &
      pB => vertices(4:6), &
      pC => vertices(7:9), &
      pD => vertices(10:12), &
      pE => vertices(13:15), &
      pF => vertices(16:18) )

    select case (p%side)
    case (1)
      S(1) = S(1)+p%weight
    case (2)
      normal = determine_normal(pA,pB,pE,pD)
      if(p%dir(3).le.normal(3)) then
        S(2) = S(2)+p%weight
      else
        S(3) = S(3)+p%weight
      endif
    case (3)
      normal = -determine_normal(pA,pC,pF,pD)
      if(p%dir(3).le.normal(3)) then
        S(4) = S(4)+p%weight
      else
        S(5) = S(5)+p%weight
      endif
    case (4)
      normal = determine_normal(pB,pC,pF,pE)
      if(p%dir(3).le.normal(3)) then
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
  end associate

  select type(bmc)
  end select

  !print *,'updating S',S,'(',p%side,')','+',p%weight
  !call print_photon(p)
  !stop 'debug update_diff_stream_wedge_5_8'
  contains
    function determine_normal(A,B,C,D)
      real(ireal_dp), intent(in) :: A(:), B(:), C(:), D(:)
      real(ireal_dp) :: determine_normal(size(A))

      if(p%subface.eq.1) then
        determine_normal = compute_normal_3d(A,B,C)
      elseif(p%subface.eq.2) then
        determine_normal = compute_normal_3d(A,C,D)
      else
        call CHKERR(1_mpiint, 'subface not determined?'//itoa(p%subface))
      endif
    end function
end subroutine

subroutine init_diff_photon_wedge_5_8(bmc, p, src, vertices, ierr)
  use m_helper_functions_dp, only : pnt_in_triangle,rotate_angle_x, rotate_angle_y, rotate_angle_z, angle_between_two_vec, rad2deg
  use m_helper_functions, only : rotation_matrix_local_basis_to_world
  class(t_boxmc_wedge_5_8) :: bmc
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: vertices(:)
  integer(iintegers),intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi,mu
  real(ireal_dp), dimension(3) :: normal, triangle_center, dir_p_to_triangle_center, rotvec, U, V
  real(ireal_dp) :: r1
  real(ireals) :: Mrot(3,3)

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
  case (2:3)
    call rand_pnt_on_plane(pA,pD,pE,pB, p%loc, normal, U, V)
    p%src_side = 2
  case (4:5)
    call rand_pnt_on_plane(pA,pD,pF,pC, p%loc, normal, U, V)
    p%src_side = 3
  case (6:7)
    call rand_pnt_on_plane(pB,pE,pF,pC, p%loc, normal, U, V)
    p%src_side = 4
  case (8)
    p%loc = rand_pnt_on_triangle(pA,pB,pC)
    normal = compute_normal_3d(pA,pB,pC)
    p%src_side = 5
  case default
    call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
    call exit
  end select

  triangle_center = mean([pA,pB,pC,pD,pE,pF])

  select case(src)
  case(1,8)
    phi= deg2rad( R()*360 )
  case(2:7)
    phi= deg2rad( interv_R(0._ireal_dp, 180._ireal_dp) ) ! from +x to -x through +y
  end select

  mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu --> Lambertian
  mu = max(cos(deg2rad(89.9_ireals)), mu)

  p%dir = (/cos(phi)*sin(acos(mu)) , sin(phi)*sin(acos(mu)) , mu  /)

  select case(src)
  case(1)
    p%dir = rotate_angle_y(p%dir,180*one) ! flip from upward to downward direction

  case(2)
    Mrot = rotation_matrix_local_basis_to_world( U, V, normal )
    p%dir = matmul(real(Mrot,ireal_dp), p%dir)

  case(3)
    !print *,'normal dir', rotate_angle_x(rotate_angle_z(p%dir, 180*one), 90*one)
    Mrot = rotation_matrix_local_basis_to_world( -U, -V, normal )
    p%dir = matmul(real(Mrot,ireal_dp), p%dir)

  case(4)
    !print *,'simple dir', rotate_angle_z(rotate_angle_x(p%dir, 90*one), 120*one)
    Mrot = rotation_matrix_local_basis_to_world( -U, V, -normal )
    p%dir = matmul(real(Mrot,ireal_dp), p%dir)

  case(5)
    !print *,'simple dir', rotate_angle_z(rotate_angle_x(rotate_angle_z(p%dir, 180*one), 90*one), 120*one)
    Mrot = rotation_matrix_local_basis_to_world( U, -V, -normal )
    p%dir = matmul(real(Mrot,ireal_dp), p%dir)

  case(6)
    !print *,'simple dir', rotate_angle_z(rotate_angle_x(rotate_angle_z(p%dir, 0*one), 90*one), 240*one)
    !print *, 'U          ', U
    !print *, 'V          ', V
    !print *, 'face normal', normal
    Mrot = rotation_matrix_local_basis_to_world( U, V, normal )
    p%dir = matmul(real(Mrot,ireal_dp), p%dir)
    !print *,'rot dir', p%dir

  case(7)
    Mrot = rotation_matrix_local_basis_to_world( -U, -V, normal )
    p%dir = matmul(real(Mrot,ireal_dp), p%dir)

  case (8)
    continue ! already in correct direction

  end select
  end associate

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
    p%side, p%weight, max_dist, p%subface)

  select type(bmc)
  end select
end subroutine
! ******************* BOXMC wedge 5 8 ************************

