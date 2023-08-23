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
! 2 -> side face between (A,C)/(B,D)
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
! Stream dst numbering is
! 2 -> top
! 1 -> bottom face
! 4 -> side face between (A,C)
! 6 -> side face between (A,B)

subroutine init_dir_photon_3_6(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  class(t_boxmc_3_6) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers), intent(in) :: src
  logical, intent(in) :: ldirect
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp), dimension(3) :: normal, U, V, cube_center, dir_p_to_cube_center

  p%dir = initial_dir
  p%alive = .false.

  call CHKERR(size(vertices) - 24, 'wrong size of vertices, expected 24 but found ', size(vertices))
  associate ( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18), &
    pG => vertices(19:21), &
    pH => vertices(22:24) &
    )
    select case (src)
    case (1)
      if (p%dir(3) .le. 0) then
        call rand_pnt_on_plane(pE, pG, pH, pF, p%loc, normal, U, V)
        p%src_side = 1
      else
        call rand_pnt_on_plane(pA, pB, pD, pC, p%loc, normal, U, V)
        p%src_side = 2
      end if
    case (2)
      if (p%dir(1) .ge. 0) then
        call rand_pnt_on_plane(pA, pC, pG, pE, p%loc, normal, U, V)
        p%src_side = 3
      else
        call rand_pnt_on_plane(pB, pF, pH, pD, p%loc, normal, U, V)
        p%src_side = 4
      end if
    case (3)
      if (p%dir(2) .ge. 0) then
        call rand_pnt_on_plane(pE, pF, pB, pA, p%loc, normal, U, V)
        p%src_side = 5
      else
        call rand_pnt_on_plane(pG, pC, pD, pH, p%loc, normal, U, V)
        p%src_side = 6
      end if
    case default
      print *, 'Dont know what to do with source spec:', src
      call exit
    end select

    ! always nudge the photon a little bit inside the volume
    cube_center = (pA + pB + pC + pD + pE + pF + pG + pH) / 8._ireal_dp
    dir_p_to_cube_center = cube_center - p%loc

    p%loc = p%loc + sqrt(epsilon(p%loc)) * dir_p_to_cube_center
  end associate

  p%weight = one
  p%alive = .true.
  p%direct = ldirect
  p%side = int(nil)
  p%src = src
  p%scattercnt = 0

  select type (bmc)
  end select

  if (dot_product(p%dir, normal) .lt. zero) then
    ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
    call CHKWARN(1_mpiint, 'initial direction of the photon is pointing outside of the cube')
    print *, 'UUPS, are we starting in a wrong direction?'// &
      & ' :: init_dir_photon_boxmc_3_6', rad2deg(angle_between_two_vec(p%dir, normal))
    print *, 'src', src, 'initial_dir', initial_dir
    print *, 'normal', normal, 'direction', p%dir
    call print_photon(p)
    call CHKERR(ierr)
    return
  end if

  ierr = 0
end subroutine

subroutine update_dir_stream_3_6(bmc, vertices, p, T)
  class(t_boxmc_3_6) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(in) :: p
  real(irealbmc), intent(inout) :: T(:)
  select case (p%side)
  case (1, 2)
    T(1) = T(1) + p%weight
  case (3, 4)
    T(2) = T(2) + p%weight
  case (5, 6)
    T(3) = T(3) + p%weight
  case default
    stop 'called update dirstream with incorrect source'
    print *, 'remove unused compiler warning', vertices
  end select
  select type (bmc)
  end select
end subroutine

subroutine update_diff_stream_3_6(bmc, vertices, p, S)
  class(t_boxmc_3_6) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(in) :: p
  real(irealbmc), intent(inout) :: S(:)

  select case (p%side)
  case (1)
    S(1) = S(1) + p%weight
  case (2)
    S(2) = S(2) + p%weight
  case (3)
    S(3) = S(3) + p%weight
  case (4)
    S(4) = S(4) + p%weight
  case (5)
    S(5) = S(5) + p%weight
  case (6)
    S(6) = S(6) + p%weight
  case default
    print *, 'Couldnt find a stream on to which I can put the photon weight on?!'
    call print_photon(p)
    call exit
    print *, 'remove unused compiler warning', vertices
  end select
  select type (bmc)
  end select
end subroutine

subroutine init_diff_photon_3_6(bmc, p, src, vertices, ierr)
  class(t_boxmc_3_6) :: bmc
  type(t_photon), intent(inout) :: p
  integer(iintegers), intent(in) :: src
  real(ireal_dp), intent(in) :: vertices(:)
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi, mu
  real(ireal_dp), dimension(3) :: normal, U, V, cube_center, dir_p_to_cube_center
  real(ireal_dp) :: Mrot(3, 3)

  call CHKERR(size(vertices) - 24, 'wrong size of vertices, expected 24 but found ', size(vertices))
  associate ( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18), &
    pG => vertices(19:21), &
    pH => vertices(22:24) &
    )
    select case (src)
    case (-1) ! thermal emission
      call CHKERR(-1_mpiint, 'thermal emission on distorted cubes currently not implemented')
    case (1) ! eup
      call rand_pnt_on_plane(pA, pB, pD, pC, p%loc, normal, U, V)
      p%src_side = PPRTS_BOT_FACE
    case (2) ! edn
      call rand_pnt_on_plane(pE, pG, pH, pF, p%loc, normal, U, V)
      p%src_side = PPRTS_TOP_FACE
    case (3)
      call rand_pnt_on_plane(pB, pF, pH, pD, p%loc, normal, U, V)
      p%src_side = PPRTS_RIGHT_FACE
    case (4)
      call rand_pnt_on_plane(pA, pC, pG, pE, p%loc, normal, U, V)
      p%src_side = PPRTS_LEFT_FACE
    case (5)
      call rand_pnt_on_plane(pG, pC, pD, pH, p%loc, normal, U, V)
      p%src_side = PPRTS_FRONT_FACE
    case (6)
      call rand_pnt_on_plane(pE, pF, pB, pA, p%loc, normal, U, V)
      p%src_side = PPRTS_REAR_FACE
    case default
      call CHKERR(1_mpiint, 'Dont know what to do with source spec: ', src)
    end select
    ! always nudge the photon a little bit inside the volume
    cube_center = (pA + pB + pC + pD + pE + pF + pG + pH) / 8._ireal_dp
    dir_p_to_cube_center = cube_center - p%loc

    p%loc = p%loc + sqrt(epsilon(p%loc)) * dir_p_to_cube_center
  end associate

  phi = deg2rad(R() * 360)
  mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
  !mu = R() ! evenly distribution between 0 and 90 degree

  p%dir = [sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu]

  select case (src)
  case (1)
    Mrot = rotation_matrix_local_basis_to_world(U, V, normal)
    p%dir = matmul(real(Mrot, ireal_dp), p%dir)
  case (2)
    Mrot = rotation_matrix_local_basis_to_world(U, V, normal)
    p%dir = matmul(real(Mrot, ireal_dp), p%dir)
  case (3)
    p%dir = rotate_angle_z(p%dir, -90._ireal_dp)
    p%dir = rotate_angle_y(p%dir, -90._ireal_dp)
    p%dir(1) = -p%dir(1)
  case (4)
    p%dir = rotate_angle_z(p%dir, -90._ireal_dp)
    p%dir = rotate_angle_y(p%dir, -90._ireal_dp)
  case (5)
    p%dir = rotate_angle_x(p%dir, -90._ireal_dp)
  case (6)
    p%dir = rotate_angle_x(p%dir, -90._ireal_dp)
    p%dir(2) = -p%dir(2)
  end select

  p%weight = one
  p%alive = .true.
  p%direct = .false.
  p%side = int(nil)
  p%src = src
  p%scattercnt = 0

  if (.not. approx(norm2(p%dir), one)) stop 'init_diff_photon_3_6 :: direction was not normed?!'

  select type (bmc)
  end select

  if (dot_product(p%dir, normal) .lt. zero) then
    ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
    call CHKWARN(1_mpiint, 'initial direction of the photon is pointing outside of the cube')
    print *, 'UUPS, are we starting in a wrong direction?'// &
      & ' :: init_diff_photon_bmx_3_6', rad2deg(angle_between_two_vec(p%dir, normal))
    print *, 'src', src
    print *, 'normal', normal, 'direction', p%dir
    call print_photon(p)
    call CHKERR(ierr)
    return
  end if

  ierr = 0
end subroutine

subroutine intersect_distance_3_6(bmc, vertices, p, max_dist)
  class(t_boxmc_3_6) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist
  integer(mpiint) :: i, ierr

  if (.not. pnt_in_cube(vertices, p%loc)) then
    print *, 'Before intersection: pnt is in cube? ', pnt_in_cube(vertices, p%loc)
    do i = 1, 8
      print *, 'vertex ', i, ':', vertices(1 + (i - 1) * 3:i * 3)
    end do
    call print_photon(p)
    call CHKERR(1_mpiint, cstr('photon is not inside the local cube', 'red'))
  end if

  call intersect_cube_general(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    & p%side, p%weight, max_dist, p%subface, ierr)

  !print *,cstr('After intersection:','red'), pnt_in_cube(vertices, p%loc), ':', vertices
  !call print_photon(p)
  select type (bmc)
  end select
end subroutine
