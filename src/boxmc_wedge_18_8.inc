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
! Horizontal View on Top Face
!
!        F
!       / \
!      / 3 \
!     /--x--\
!    / 1 | 2 \
!   D ------- E
!
!
! View on the vertical side of the base face, left face and right face
!
!   D---------E  D---------F  F---------E
!   | 4  | 5  |  | 8  | 9  |  | 12 | 13 |
!   |____|____|  |____|____|  |____|____|
!   |    |    |  |    |    |  |    |    |
!   | 6  | 7  |  | 10 | 11 |  | 14 | 15 |
!   A---------B  A---------C  C---------B
!
! Horizontal View on Bot Face
!
!        C
!       / \
!      /18 \
!     /--x--\
!    / 16|17 \
!   A ------- B
!
!
! Face numbering is
! 1 -> top
! 2 -> side face between (A,B)
! 3 -> side face between (A,C)
! 5 -> bottom face
!
! We always assume the triangle to have dx edge length along y=0 (A,B) and dy is the edge length between (A,C)
! Distribute Photons on triangles: https://doi.org/10.1145/571647.571648

subroutine init_dir_photon_wedge_18_8(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  class(t_boxmc_wedge_18_8) :: bmc
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers),intent(in) :: src
  logical,intent(in) :: ldirect
  real(ireal_dp), dimension(3) :: normal, triangle_center, dir_p_to_triangle_center, U, V
  real(ireal_dp), dimension(3) :: top_center, bot_center, side_center_base, side_center_left, side_center_right
  real(ireal_dp), dimension(3) :: pmAB, pmAC, pmBC, pmDE, pmDF, pmEF
  real(ireal_dp), dimension(3) :: pmAD, pmBE, pmCF
  real(ireal_dp),parameter :: edge_start_range=.999_ireal_dp
  integer(mpiint),intent(out) :: ierr

  associate( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18) )

    bot_center = (pA+pB+pC)/3
    top_center = (pD+pE+pF)/3

    side_center_base = (pA+pB+pD+pE)/4
    side_center_left = (pA+pC+pD+pF)/4
    side_center_right= (pC+pB+pF+pE)/4

    pmAB = (pA+pB)/2
    pmAC = (pA+pC)/2
    pmBC = (pB+pC)/2
    pmDE = (pD+pE)/2
    pmDF = (pD+pF)/2
    pmEF = (pE+pF)/2

    pmAD = (pA+pD)/2
    pmBE = (pB+pE)/2
    pmCF = (pC+pF)/2

  select case(src)
  case (1)
    call rand_pnt_on_plane(pD,pmDF,top_center,pmDE, p%loc, normal, U, V)
  case (2)
    call rand_pnt_on_plane(pE,pmDE,top_center,pmEF, p%loc, normal, U, V)
  case (3)
    call rand_pnt_on_plane(pF,pmEF,top_center,pmDF, p%loc, normal, U, V)

  case (4)
    call rand_pnt_on_plane(pD,pmDE,side_center_base,pmAD, p%loc, normal, U, V)
  case (5)
    call rand_pnt_on_plane(pE,pmBE,side_center_base,pmDE, p%loc, normal, U, V)
  case (6)
    call rand_pnt_on_plane(pA,pmAD,side_center_base,pmAB, p%loc, normal, U, V)
  case (7)
    call rand_pnt_on_plane(pB,pmAB,side_center_base,pmBE, p%loc, normal, U, V)


  case (8)
    call rand_pnt_on_plane(pD,pmAD,side_center_left,pmDF, p%loc, normal, U, V)
  case (9)
    call rand_pnt_on_plane(pF,pmDF,side_center_left,pmCF, p%loc, normal, U, V)
  case (10)
    call rand_pnt_on_plane(pA,pmAC,side_center_left,pmAD, p%loc, normal, U, V)
  case (11)
    call rand_pnt_on_plane(pC,pmCF,side_center_left,pmAC, p%loc, normal, U, V)

  case (12)
    call rand_pnt_on_plane(pF,pmCF,side_center_right,pmEF, p%loc, normal, U, V)
  case (13)
    call rand_pnt_on_plane(pE,pmEF,side_center_right,pmBE, p%loc, normal, U, V)
  case (14)
    call rand_pnt_on_plane(pC,pmBC,side_center_right,pmCF, p%loc, normal, U, V)
  case (15)
    call rand_pnt_on_plane(pB,pmBE,side_center_right,pmBC, p%loc, normal, U, V)

  case (16)
    call rand_pnt_on_plane(pA,pmAB,bot_center,pmAC, p%loc, normal, U, V)
  case (17)
    call rand_pnt_on_plane(pB,pmBC,bot_center,pmAB, p%loc, normal, U, V)
  case (18)
    call rand_pnt_on_plane(pC,pmAC,bot_center,pmBC, p%loc, normal, U, V)

  case default
    call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
    call exit
  end select

  triangle_center = (pA+pB+pC+pD+pE+pF)/6._ireal_dp

  select case(src)
  case (1:3)
    p%src_side = i1
  case (4:7)
    p%src_side = i2
  case (8:11)
    p%src_side = i3
  case (12:15)
    p%src_side = i4
  case (16:18)
    p%src_side = i5
  end select

  p%weight=one
  p%alive = .True.
  p%direct= ldirect
  p%side = int(nil)
  p%src  = src
  p%scattercnt = 0
  p%dir = initial_dir

  ! always nudge the photon a little bit inside the volume
  dir_p_to_triangle_center = triangle_center - p%loc

  p%loc = p%loc + sqrt(epsilon(pi64)) * dir_p_to_triangle_center

  end associate

  select type(bmc)
  end select

  if(dot_product(p%dir, normal).lt.zero) then
    ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
    if(ldebug) then
      print *,'UUPS, are we starting in a wrong direction? :: init_dir_photon_wedge_18_8', &
        rad2deg(angle_between_two_vec(p%dir, normal))
      print *,'src', src, 'initial_dir', initial_dir
      print *,'normal',normal,'direction',p%dir
      call print_photon(p)
    endif
    !stop 'initial direction of the wedge triangle'
    return
  endif
  ierr = 0
end subroutine

subroutine update_dir_stream_wedge_18_8(bmc,vertices,p,T)
  class(t_boxmc_wedge_18_8) :: bmc
  real(ireal_dp),intent(in) :: vertices(:)
  type(t_photon),intent(in) :: p
  real(ireal_dp),intent(inout) :: T(:)

  real(ireal_dp), dimension(3) :: top_center, bot_center
  real(ireal_dp), dimension(3) :: side_center_base, side_center_left, side_center_right
  real(ireal_dp), dimension(3) :: pmAB, pmAC, pmBC, pmDE, pmDF, pmEF
  real(ireal_dp), dimension(3) :: pmAD, pmBE, pmCF

  real(ireal_dp) :: origin(3), hit(4)
  logical :: lhit
  integer(iintegers) :: isubface

  origin = p%loc - p%dir
  select type(bmc)
  end select

  associate( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18) )

    bot_center = (pA+pB+pC)/3
    top_center = (pD+pE+pF)/3

    side_center_base = (pA+pB+pD+pE)/4
    side_center_left = (pA+pC+pD+pF)/4
    side_center_right= (pC+pB+pF+pE)/4

    pmAB = (pA+pB)/2
    pmAC = (pA+pC)/2
    pmBC = (pB+pC)/2
    pmDE = (pD+pE)/2
    pmDF = (pD+pF)/2
    pmEF = (pE+pF)/2

    pmAD = (pA+pD)/2
    pmBE = (pB+pE)/2
    pmCF = (pC+pF)/2

    select case(p%side)
    case(1)
      call square_intersection(origin, p%dir, pD, pmDF, top_center, pmDE, lhit, hit, isubface)
      if(lhit) then
        T(1) = T(1) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pE, pmDE, top_center, pmEF, lhit, hit, isubface)
      if(lhit) then
        T(2) = T(2) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pF, pmEF, top_center, pmDF, lhit, hit, isubface)
      if(lhit) then
        T(3) = T(3) + p%weight
        return
      endif

    case(2)
      call square_intersection(origin, p%dir, pD, pmDE, side_center_base, pmAD, lhit, hit, isubface)
      if(lhit) then
        T(4) = T(4) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pE, pmBE, side_center_base, pmDE, lhit, hit, isubface)
      if(lhit) then
        T(5) = T(5) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pA, pmAD, side_center_base, pmAB, lhit, hit, isubface)
      if(lhit) then
        T(6) = T(6) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pB, pmAB, side_center_base, pmBE, lhit, hit, isubface)
      if(lhit) then
        T(7) = T(7) + p%weight
        return
      endif

    case(3)
      call square_intersection(origin, p%dir, pD, pmAD, side_center_left, pmDF, lhit, hit, isubface)
      if(lhit) then
        T(8) = T(8) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pF, pmDF, side_center_left, pmCF, lhit, hit, isubface)
      if(lhit) then
        T(9) = T(9) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pA, pmAC, side_center_left, pmAD, lhit, hit, isubface)
      if(lhit) then
        T(10) = T(10) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pC, pmCF, side_center_left, pmAC, lhit, hit, isubface)
      if(lhit) then
        T(11) = T(11) + p%weight
        return
      endif

    case(4)
      call square_intersection(origin, p%dir, pF, pmCF, side_center_right, pmEF, lhit, hit, isubface)
      if(lhit) then
        T(12) = T(12) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pE, pmEF, side_center_right, pmBE, lhit, hit, isubface)
      if(lhit) then
        T(13) = T(13) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pC, pmBC, side_center_right, pmCF, lhit, hit, isubface)
      if(lhit) then
        T(14) = T(14) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pB, pmBE, side_center_right, pmBC, lhit, hit, isubface)
      if(lhit) then
        T(15) = T(15) + p%weight
        return
      endif

    case(5)
      call square_intersection(origin, p%dir, pA, pmAB, bot_center, pmAC, lhit, hit, isubface)
      if(lhit) then
        T(16) = T(16) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pB, pmBC, bot_center, pmAB, lhit, hit, isubface)
      if(lhit) then
        T(17) = T(17) + p%weight
        return
      endif
      call square_intersection(origin, p%dir, pC, pmAC, bot_center, pmBC, lhit, hit, isubface)
      if(lhit) then
        T(18) = T(18) + p%weight
        return
      endif

    case default
      call CHKERR(1_mpiint, 'Dont know what to do with side spec: '//itoa(p%side))
      call exit
    end select

  end associate
  call print_photon(p)
  print *,'origin', origin
  call CHKERR(1_mpiint, 'Dont know where to put the energy! dst_side: '//itoa(p%side))
end subroutine

subroutine update_diff_stream_wedge_18_8(bmc,vertices,p,S)
  class(t_boxmc_wedge_18_8) :: bmc
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
      call CHKERR(1_mpiint, 'update_diff_photon_wedge_18_8 called with unknown p%side')
    end select
  end associate

  select type(bmc)
  end select

  !print *,'updating S',S,'(',p%side,')','+',p%weight
  !call print_photon(p)
  !stop 'debug update_diff_stream_wedge_18_8'
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

subroutine init_diff_photon_wedge_18_8(bmc, p, src, vertices, ierr)
  class(t_boxmc_wedge_18_8) :: bmc
  type(t_photon),intent(inout) :: p
  real(ireal_dp),intent(in) :: vertices(:)
  integer(iintegers),intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi,mu
  real(ireal_dp), dimension(3) :: normal, triangle_center, dir_p_to_triangle_center, U, V
  real(ireal_dp) :: Mrot(3,3)

  associate( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18) )

  select case(src)
  case (1)
    call rand_pnt_on_triangle(pD,pE,pF,p%loc)
    normal = compute_normal_3d(pD,pF,pE)
    p%src_side = 1
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
    call rand_pnt_on_triangle(pA,pB,pC,p%loc)
    normal = compute_normal_3d(pA,pB,pC)
    p%src_side = 5
  case default
    call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//itoa(src))
    call exit
  end select

  triangle_center = (pA+pB+pC+pD+pE+pF)/6._ireal_dp

  select case(src)
  case(1,8)
    phi= deg2rad( R()*360 )
  case(2:7)
    phi= deg2rad( interv_R(0._ireal_dp, 180._ireal_dp) ) ! from +x to -x through +y
  end select

  mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu --> Lambertian
  mu = max(cos(deg2rad(89.9_ireal_dp)), mu)

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

  p%loc = p%loc + sqrt(epsilon(pi64)) * dir_p_to_triangle_center

  if(.not. approx( norm(p%dir) ,one ) ) call CHKERR(1_mpiint, 'init_diff_photon_3_10 :: direction was not normed?!')

  ierr = 0
  !call print_photon(p)
  !print *,'--------- end new diffuse photon ----------'
  !stop 'not yet implemented :: init_diff_photon_wedge_18_8'
  select type(bmc)
  end select
end subroutine

subroutine intersect_distance_wedge_18_8(bmc, vertices, p, max_dist)
  class(t_boxmc_wedge_18_8) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist
  integer(mpiint) :: ierr

  call intersect_wedge(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    p%side, p%weight, max_dist, p%subface, ierr)

  if(ierr.ne.0_mpiint) then
    call print_photon(p)
  endif
  call CHKERR(ierr, 'Error in intersect_wedge')

  select type(bmc)
  end select
end subroutine
subroutine wedge_halfspaces_18_8(bmc, vertices, origins, normals)
    class(t_boxmc_wedge_18_8) :: bmc
    real(ireal_dp),intent(in) :: vertices(:)
    real(ireal_dp),allocatable,intent(out) :: origins(:,:), normals(:,:)
    if(.not.allocated(origins)) allocate(origins(3,5))
    if(.not.allocated(normals)) allocate(normals(3,5))
    call wedge_halfspaces(vertices, origins, normals)
    select type(bmc)
    end select
end subroutine
! ******************* BOXMC wedge 18 8 ************************

