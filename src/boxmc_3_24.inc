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
!
! 5 streams per face
! one center stream which has angles close to the surface normal
! then 4 streams for each quadrant
!
! Src numbering
!         _2_4_6_8_______
!        |               |
!      9 |     18..24..2 |10
!     11 |               |12
!     13 |               |14
!     15 |  17..23..2    |16
!        |               |
!        |_______________|
!          1 3 5 7
!
! i.e. src numbering is
!  1..8..2 -> bot face, i.e. upwelling
!  2..8..2 -> top face, i.e. down welling
!  9..16..2 -> side faces in x-plane left to right
! 10..16..2 -> side faces in x-plane right to left
! 17..24..2 -> side faces in y-plane backward to forward
! 18..24..2 -> side faces in y-plane forward to back

! Dst numbering
!         _1_3_5_7_______
!        |               |
!     10 |     17..24..2 | 9
!     12 |               |11
!     14 |               |13
!     16 |  18..23..2    |15
!        |               |
!        |_______________|
!          2 5 6 8
!
! i.e. numbering is
!  1..8..2 -> bot face, i.e. upwelling
!  2..8..2 -> top face, i.e. down welling
!  9..16..2 -> side faces in x-plane left to right
! 10..16..2 -> side faces in x-plane right to left
! 17..24..2 -> side faces in y-plane backward to forward
! 18..24..2 -> side faces in y-plane forward to back
!
! Angular Distribution on side 1 (top) and on side 2 (bot)
!                                              !
!  y                                           !  y
!  ^  ________________                         !  ^  ________________
!  | |        |        |                       !  | |        |        |
!  | |        |        |                       !  | |        |        |
!    |    1   |   3    |                       !    |    2   |   4    |
!    |        |        |                       !    |        |        |
!    |-----------------|                       !    |-----------------|
!    |        |        |                       !    |        |        |
!    |    5   |   7    |                       !    |    6   |   8    |
!    |        |        |                       !    |        |        |
!    | _______|_______ |                       !    | _______|_______ |
!                       -->  x                 !                       -->  x
!                                              !

! Angular Distribution on side 3 (x-axis left) and on side 4 (x-axis right)
!                                              !
!  z                                           !  z
!  ^  ________________                         !  ^  ________________
!  | |        |        |                       !  | |        |        |
!  | |        |        |                       !  | |        |        |
!    |   10   |   12   |                       !    |    9   |   11   |
!    |        |        |                       !    |        |        |
!    |-----------------|                       !    |-----------------|
!    |        |        |                       !    |        |        |
!    |   14   |   16   |                       !    |   13   |   15   |
!    |        |        |                       !    |        |        |
!    | _______|_______ |                       !    | _______|_______ |
!                       -->  y                 !                       -->  y
!                                                                  !

! Angular Distribution on side 5 (y-axis backward) and on side 6 (y-axis forward)
!                                              !
!  z                                           !  z
!  ^  ________________                         !  ^  ________________
!  | |        |        |                       !  | |        |        |
!  | |        |        |                       !  | |        |        |
!    |   18   |   20   |                       !    |   17   |   19   |
!    |        |        |                       !    |        |        |
!    |-----------------|                       !    |-----------------|
!    |        |        |                       !    |        |        |
!    |   22   |   24   |                       !    |   21   |   23   |
!    |        |        |                       !    |        |        |
!    | _______|_______ |                       !    | _______|_______ |
!                       -->  x                 !                       -->  x
!                                              !

subroutine init_dir_photon_3_24(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  class(t_boxmc_3_24) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers), intent(in) :: src
  logical, intent(in) :: ldirect
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp), dimension(3) :: normal, U, V, cube_center, dir_p_to_cube_center

  p%dir = initial_dir
  p%alive = .false.

  call CHKERR(size(vertices) - 24, 'wrong size of vertices, expected 24 but found '//toStr(size(vertices)))
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
      & ' :: init_dir_photon_wedge_5_8', rad2deg(angle_between_two_vec(p%dir, normal))
    print *, 'src', src, 'initial_dir', initial_dir
    print *, 'normal', normal, 'direction', p%dir
    call print_photon(p)
    call CHKERR(ierr)
    return
  end if

  ierr = 0
end subroutine

subroutine update_dir_stream_3_24(bmc, vertices, p, T)
  class(t_boxmc_3_24) :: bmc
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
    call CHKERR(1_mpiint, 'called update dirstream with incorrect source'//toStr(p%side))
    print *, 'remove unused compiler warning', vertices
  end select
  select type (bmc)
  end select
end subroutine

subroutine init_diff_photon_3_24(bmc, p, src, vertices, ierr)
  use m_helper_functions, only: rotate_angle_x, rotate_angle_y, rotate_angle_z
  class(t_boxmc_3_24) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: vertices(:)
  integer(iintegers), intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi, mu
  real(ireal_dp), dimension(3) :: normal, U, V, cube_center, dir_p_to_cube_center

  call CHKERR(size(vertices) - 24, 'wrong size of vertices, expected 24 but found '//toStr(size(vertices)))
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
    case (1, 3, 5, 7)
      call rand_pnt_on_plane(pA, pB, pD, pC, p%loc, normal, U, V)
      p%src_side = 2
    case (2, 4, 6, 8)
      call rand_pnt_on_plane(pE, pG, pH, pF, p%loc, normal, U, V)
      p%src_side = 1
    case (9, 11, 13, 15)
      call rand_pnt_on_plane(pA, pC, pG, pE, p%loc, normal, U, V)
      p%src_side = 3
    case (10, 12, 14, 16)
      call rand_pnt_on_plane(pB, pF, pH, pD, p%loc, normal, U, V)
      p%src_side = 4
    case (17, 19, 21, 23)
      call rand_pnt_on_plane(pE, pF, pB, pA, p%loc, normal, U, V)
      p%src_side = 5
    case (18, 20, 22, 24)
      call rand_pnt_on_plane(pG, pC, pD, pH, p%loc, normal, U, V)
      p%src_side = 6
    case default
      call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
    end select

    select case (src)
    case (1, 2, 9, 10, 17, 18)
      phi = deg2rad(interv_R(-90._ireal_dp, 0._ireal_dp))
    case (3, 4, 11, 12, 19, 20)
      phi = deg2rad(interv_R(0._ireal_dp, 90._ireal_dp))
    case (5, 6, 13, 14, 21, 22)
      phi = deg2rad(interv_R(180._ireal_dp, 270._ireal_dp))
    case (7, 8, 15, 16, 23, 24)
      phi = deg2rad(interv_R(90._ireal_dp, 180._ireal_dp))
    case default
      phi = -999
      call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
    end select

    select case (src)
    case (1:24)
      mu = sqrt(interv_R(zero, one))
    case default
      mu = -1
      call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
    end select

    p%dir = (/sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu/)
    ! phi = 0 is y-axes

    select case (src)
    case (1, 3, 5, 7)
      continue
    case (2, 4, 6, 8)
      p%dir(3) = -p%dir(3)
    case (10, 12, 14, 16)
      p%dir = rotate_angle_z(p%dir, -90._ireal_dp)
      p%dir = rotate_angle_y(p%dir, -90._ireal_dp)
      p%dir(1) = -p%dir(1)
    case (9, 11, 13, 15)
      p%dir = rotate_angle_z(p%dir, -90._ireal_dp)
      p%dir = rotate_angle_y(p%dir, -90._ireal_dp)
    case (18, 20, 22, 24)
      p%dir = rotate_angle_x(p%dir, -90._ireal_dp)
    case (17, 19, 21, 23)
      p%dir = rotate_angle_x(p%dir, -90._ireal_dp)
      p%dir(2) = -p%dir(2)
    case default
      call CHKERR(1_mpiint, 'Dont know what to do with source spec: '//toStr(src))
    end select

    ! always nudge the photon a little bit inside the volume
    cube_center = (pA + pB + pC + pD + pE + pF + pG + pH) / 8._ireal_dp
    dir_p_to_cube_center = cube_center - p%loc

    p%loc = p%loc + sqrt(epsilon(p%loc)) * dir_p_to_cube_center
  end associate

  p%weight = one
  p%alive = .true.
  p%direct = .false.
  p%side = int(nil)
  p%src = src
  p%scattercnt = 0

  if (.not. approx(norm2(p%dir), one)) then
    print *, 'p%dir', p%dir, 'norm', norm2(p%dir)
    stop 'init_diff_photon_3_24 :: direction was not normed?!'
  end if

  select type (bmc)
  end select

  if (dot_product(p%dir, normal) .lt. zero) then
    ierr = 1000_mpiint + int(angle_between_two_vec(p%dir, normal), kind=mpiint)
    call CHKWARN(1_mpiint, 'initial direction of the photon is pointing outside of the cube')
    print *, 'UUPS, are we starting in a wrong direction?'// &
      & ' :: init_diff_photon_bmx_3_24', rad2deg(angle_between_two_vec(p%dir, normal))
    print *, 'src', src
    print *, 'normal', normal, 'direction', p%dir
    call print_photon(p)
    call CHKERR(ierr)
    return
  end if

  ierr = 0
  !call print_photon(p)
  !call CHKERR(1_mpiint, 'DEBUG')
end subroutine

subroutine update_diff_stream_3_24(bmc, vertices, p, S)
  class(t_boxmc_3_24) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(in) :: p
  real(irealbmc), intent(inout) :: S(:)

  if (p%side .eq. 1) then
    if (p%dir(1) .gt. 0) then
      if (p%dir(2) .gt. 0) then
        S(3) = S(3) + p%weight
      else
        S(7) = S(7) + p%weight
      end if
    else
      if (p%dir(2) .gt. 0) then
        S(1) = S(1) + p%weight
      else
        S(5) = S(5) + p%weight
      end if
    end if
    return

  else if (p%side .eq. 2) then
    if (p%dir(1) .gt. 0) then
      if (p%dir(2) .gt. 0) then
        S(4) = S(4) + p%weight
      else
        S(8) = S(8) + p%weight
      end if
    else
      if (p%dir(2) .gt. 0) then
        S(2) = S(2) + p%weight
      else
        S(6) = S(6) + p%weight
      end if
    end if
    return

  else if (p%side .eq. 3) then
    if (p%dir(2) .gt. 0) then
      if (p%dir(3) .gt. 0) then
        S(12) = S(12) + p%weight
      else
        S(16) = S(16) + p%weight
      end if
    else
      if (p%dir(3) .gt. 0) then
        S(10) = S(10) + p%weight
      else
        S(14) = S(14) + p%weight
      end if
    end if
    return

  else if (p%side .eq. 4) then
    if (p%dir(2) .gt. 0) then
      if (p%dir(3) .gt. 0) then
        S(11) = S(11) + p%weight
      else
        S(15) = S(15) + p%weight
      end if
    else
      if (p%dir(3) .gt. 0) then
        S(9) = S(9) + p%weight
      else
        S(13) = S(13) + p%weight
      end if
    end if
    return

  else if (p%side .eq. 5) then
    if (p%dir(1) .gt. 0) then
      if (p%dir(3) .gt. 0) then
        S(20) = S(20) + p%weight
      else
        S(24) = S(24) + p%weight
      end if
    else
      if (p%dir(3) .gt. 0) then
        S(18) = S(18) + p%weight
      else
        S(22) = S(22) + p%weight
      end if
    end if
    return

  else if (p%side .eq. 6) then
    if (p%dir(1) .gt. 0) then
      if (p%dir(3) .gt. 0) then
        S(19) = S(19) + p%weight
      else
        S(23) = S(23) + p%weight
      end if
    else
      if (p%dir(3) .gt. 0) then
        S(17) = S(17) + p%weight
      else
        S(21) = S(21) + p%weight
      end if
    end if
    return

  else
    call print_photon(p)
    call CHKERR(1_mpiint, 'Couldnt find a stream on to which I can put the photon weight on?!')
    print *, 'remove unused compiler warning', vertices
  end if
  select type (bmc)
  end select
end subroutine

subroutine intersect_distance_3_24(bmc, vertices, p, max_dist)
  class(t_boxmc_3_24) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist
  integer(mpiint) :: ierr
  call intersect_cube_general(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    & p%side, p%weight, max_dist, p%subface, ierr)

  !call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
  !  & p%side, max_dist)

  select type (bmc)
  end select
end subroutine
subroutine box_halfspaces_3_24(bmc, vertices, origins, normals)
  class(t_boxmc_3_24) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  real(ireal_dp), allocatable, intent(out) :: origins(:, :), normals(:, :)
  call box_halfspaces(vertices, origins, normals)
  select type (bmc)
  end select
end subroutine
! ******************* BOXMC 3 30 ************************

