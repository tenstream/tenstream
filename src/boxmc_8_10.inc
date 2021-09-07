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

subroutine update_dir_stream_8_10(bmc, vertices, p, T)
  class(t_boxmc_8_10) :: bmc
  type(t_photon), intent(in) :: p
  real(ireal_dp), intent(in) :: vertices(:)
  real(irealbmc), intent(inout) :: T(:)
  real(ireal_dp), parameter :: h = .5_ireal_dp
  real(ireal_dp) :: dx, dy, dz

!          /     3       4 /
!         / 1       2     /
!        /_______________/
!        |               |
!      6 |           8   |
!        |               |
!        |               |
!        |    7          |
!      5 |_______________|

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  select case (p%side)
  case (1:2)
    if (p%loc(1) .le. dx * h .and. p%loc(2) .le. dy * h) then
      T(1) = T(1) + p%weight
      return
    else if (p%loc(1) .gt. dx * h .and. p%loc(2) .le. dy * h) then
      T(2) = T(2) + p%weight
      return
    else if (p%loc(1) .le. dx * h .and. p%loc(2) .gt. dy * h) then
      T(3) = T(3) + p%weight
      return
    else if (p%loc(1) .gt. dx * h .and. p%loc(2) .gt. dy * h) then
      T(4) = T(4) + p%weight
      return
    else
      print *, 'Couldnt find a stream on to which I can put the photon weight on?!'
      call print_photon(p)
    end if

  case (3:4)

    if (p%loc(3) .le. dz * h) then
      T(5) = T(5) + p%weight
      return
    else if (p%loc(3) .gt. dz * h) then
      T(6) = T(6) + p%weight
      return
    else
      print *, 'Couldnt find a stream on to which I can put the photon weight on?!'
      call print_photon(p)
    end if
  case (5:6)

    if (p%loc(3) .le. dz * h) then
      T(7) = T(7) + p%weight
      return
    else if (p%loc(3) .gt. dz * h) then
      T(8) = T(8) + p%weight
      return
    else
      print *, 'Couldnt find a stream on to which I can put the photon weight on?!'
      call print_photon(p)
    end if
  case default
    print *, 'Dont know what to do with this p%side'
    call print_photon(p)
  end select
  select type (bmc)
  end select

end subroutine
subroutine update_diff_stream_8_10(bmc, vertices, p, S)
  class(t_boxmc_8_10) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(in) :: p
  real(irealbmc), intent(inout) :: S(:)

!         _______2_______
!        |            9  |
!      6 |            7  | 5
!        |               |
!        |   10          |
!        |   8           |
!      4 |_______________| 3
!                1

  if (p%side .eq. 1) then
    S(1) = S(1) + p%weight
    return

  else if (p%side .eq. 2) then
    S(2) = S(2) + p%weight
    return

  else if (p%side .eq. 3 .and. p%dir(3) .le. zero) then
    S(3) = S(3) + p%weight
    return
  else if (p%side .eq. 3 .and. p%dir(3) .gt. zero) then
    S(5) = S(5) + p%weight
    return

  else if (p%side .eq. 4 .and. p%dir(3) .le. zero) then
    S(4) = S(4) + p%weight
    return
  else if (p%side .eq. 4 .and. p%dir(3) .gt. zero) then
    S(6) = S(6) + p%weight
    return

  else if (p%side .eq. 5 .and. p%dir(3) .le. zero) then
    S(7) = S(7) + p%weight
    return
  else if (p%side .eq. 5 .and. p%dir(3) .gt. zero) then
    S(9) = S(9) + p%weight
    return

  else if (p%side .eq. 6 .and. p%dir(3) .le. zero) then
    S(8) = S(8) + p%weight
    return
  else if (p%side .eq. 6 .and. p%dir(3) .gt. zero) then
    S(10) = S(10) + p%weight
    return
  else
    print *, 'Couldnt find a stream on to which I can put the photon weight on?!'
    call print_photon(p)
    call exit
    print *, 'remove unused compiler warning', vertices
  end if
  select type (bmc)
  end select

end subroutine

subroutine init_dir_photon_8_10(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  class(t_boxmc_8_10) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers), intent(in) :: src
  logical, intent(in) :: ldirect
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp), dimension(3) :: normal, U, V, cube_center, dir_p_to_cube_center

  p%dir = initial_dir
  p%alive = .false.

!          /     3       4 /
!         / 1       2     /
!        /_______________/
!        |               |
!      6 |           8   |
!        |               |
!        |               |
!        |    7          |
!      5 |_______________|

  call CHKERR(size(vertices) - 24, 'wrong size of vertices, expected 24 but found '//toStr(size(vertices)))
  associate ( &
    pA => vertices(1:3), &
    pB => vertices(4:6), &
    pC => vertices(7:9), &
    pD => vertices(10:12), &
    pE => vertices(13:15), &
    pF => vertices(16:18), &
    pG => vertices(19:21), &
    pH => vertices(22:24), &
    mid_bot => (vertices(1:3) + vertices(4:6) + vertices(7:9) + vertices(10:12)) / 4, &
    mid_top => (vertices(13:15) + vertices(16:18) + vertices(19:21) + vertices(22:24)) / 4, &
    AB => (vertices(1:3) + vertices(4:6)) / 2, &
    CD => (vertices(7:9) + vertices(10:12)) / 2, &
    AC => (vertices(1:3) + vertices(7:9)) / 2, &
    BD => (vertices(4:6) + vertices(10:12)) / 2, &
    EF => (vertices(13:15) + vertices(16:18)) / 2, &
    GH => (vertices(19:21) + vertices(22:24)) / 2, &
    EG => (vertices(13:15) + vertices(19:21)) / 2, &
    FH => (vertices(16:18) + vertices(22:24)) / 2, &
    AE => (vertices(1:3) + vertices(13:15)) / 2, &
    BF => (vertices(4:6) + vertices(16:18)) / 2, &
    CG => (vertices(7:9) + vertices(19:21)) / 2, &
    DH => (vertices(10:12) + vertices(22:24)) / 2 &
    )
    select case (src)
    case (1)
      !p%loc = [ L(dx)*h     , L(dy)*h     , dz*(one-e)   ]
      if (p%dir(3) .le. 0) then
        call rand_pnt_on_plane(pE, EG, mid_top, EF, p%loc, normal, U, V)
        p%src_side = 1
      else
        call rand_pnt_on_plane(pA, AC, mid_bot, AB, p%loc, normal, U, V)
        p%src_side = 2
      end if
    case (2)
      !p%loc = [ dx*h+L(dx)*h, L(dy)*h     , dz*(one-e)   ]
      if (p%dir(3) .le. 0) then
        call rand_pnt_on_plane(EF, mid_top, FH, pF, p%loc, normal, U, V)
        p%src_side = 1
      else
        call rand_pnt_on_plane(AB, mid_bot, BD, pB, p%loc, normal, U, V)
        p%src_side = 2
      end if
    case (3)
      !p%loc = [ L(dx)*h     , dy*h+L(dy)*h, dz*(one-e)   ]
      if (p%dir(3) .le. 0) then
        call rand_pnt_on_plane(EG, pG, GH, mid_top, p%loc, normal, U, V)
        p%src_side = 1
      else
        call rand_pnt_on_plane(AC, pC, CD, mid_bot, p%loc, normal, U, V)
        p%src_side = 2
      end if
    case (4)
      !p%loc = [ dx*h+L(dx)*h, dy*h+L(dy)*h, dz*(one-e)   ]
      if (p%dir(3) .le. 0) then
        call rand_pnt_on_plane(GH, pH, FH, mid_top, p%loc, normal, U, V)
        p%src_side = 1
      else
        call rand_pnt_on_plane(CD, pD, BD, mid_bot, p%loc, normal, U, V)
        p%src_side = 2
      end if
    case (5)
      if (p%dir(1) .ge. 0) then
        p%src_side = 3
        !p%loc = [ zero+dx*e   , L(dy)       , L(dz)*h      ]
        call rand_pnt_on_plane(pA, pC, CG, AE, p%loc, normal, U, V)
      else
        p%src_side = 4
        !p%loc = [ dx  -dx*e   , L(dy)       , L(dz)*h      ]
        call rand_pnt_on_plane(pB, BF, DH, pD, p%loc, normal, U, V)
      end if
    case (6)
      if (p%dir(1) .ge. 0) then
        p%src_side = 3
        !p%loc = [ zero+dx*e   , L(dy)       , dz*h+L(dz)*h ]
        call rand_pnt_on_plane(AE, CG, pG, pE, p%loc, normal, U, V)
      else
        p%src_side = 4
        !p%loc = [ dx  -dx*e   , L(dy)       , dz*h+L(dz)*h ]
        call rand_pnt_on_plane(BF, pF, pH, DH, p%loc, normal, U, V)
      end if
    case (7)
      if (p%dir(2) .ge. 0) then
        p%src_side = 5
        !p%loc = [ L(dx)       , zero+dy*e   , L(dz)*h      ]
        call rand_pnt_on_plane(pA, AE, BF, pB, p%loc, normal, U, V)
      else
        p%src_side = 6
        !p%loc = [ L(dx)       , dy  -dy*e   , L(dz)*h      ]
        call rand_pnt_on_plane(pC, pD, DH, CG, p%loc, normal, U, V)
      end if
    case (8)
      if (p%dir(2) .ge. 0) then
        p%src_side = 5
        !p%loc = [ L(dx)       , zero+dy*e   , dz*h+L(dz)*h ]
        call rand_pnt_on_plane(AE, pE, pF, BF, p%loc, normal, U, V)
      else
        p%src_side = 6
        !p%loc = [ L(dx)       , dy  -dy*e   , dz*h+L(dz)*h ]
        call rand_pnt_on_plane(pG, CG, DH, pH, p%loc, normal, U, V)
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

  ierr = 0
end subroutine

subroutine init_diff_photon_8_10(bmc, p, src, vertices, ierr)
  class(t_boxmc_8_10) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: vertices(:)
  integer(iintegers), intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi, mu, mu1, mu2
  real(ireal_dp) :: dx, dy, dz

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  mu1 = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian
  mu2 = R() ! evenly distribution between 0 and 90 degree

  mu = mu1
  select case (src)
  case (-1) ! thermal emission
    p%loc = (/L(dx), L(dy), L(dz)/)
    mu = 2 * mu2 - one
    phi = deg2rad(R() * 360)
  case (1)
    p%loc = (/L(dx), L(dy), zero/)
    phi = deg2rad(R() * 360)
  case (2)
    p%loc = (/L(dx), L(dy), dz/)
    phi = deg2rad(R() * 360)
  case (3)
    p%loc = (/dx, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (4)
    p%loc = (/zero, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (5)
    p%loc = (/dx, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (6)
    p%loc = (/zero, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (7)
    p%loc = (/L(dx), dy, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (8)
    p%loc = (/L(dx), zero, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (9)
    p%loc = (/L(dx), dy, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (10)
    p%loc = (/L(dx), zero, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case default
    print *, 'Dont know what to do with source spec:', src
    call exit
  end select
  select case (src)
  case (1:2)
    phi = deg2rad(R() * 360)
  case (3:10)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  end select

  p%dir = (/sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu/)
  ! phi = 0 is y-axes

  select case (src)
  case (2)
    p%dir = rotate_angle_y(p%dir, 180 * one)
  case (3)
    p%dir = rotate_angle_y(p%dir, 90 * one)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (4)
    p%dir = rotate_angle_y(p%dir, -90 * one)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (5)
    p%dir = rotate_angle_y(p%dir, 90 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (6)
    p%dir = rotate_angle_y(p%dir, -90 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (7)
    p%dir = rotate_angle_z(p%dir, 180 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (8)
    p%dir = rotate_angle_x(p%dir, 90 * one)

  case (9)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (10)
    p%dir = rotate_angle_z(p%dir, 180 * one)
    p%dir = rotate_angle_x(p%dir, 90 * one)

  end select
  p%weight = one
  p%alive = .true.
  p%direct = .false.
  p%side = int(nil)
  p%src = src
  p%scattercnt = 0
  if (.not. approx(norm2(p%dir), one)) then
    print *, 'p%dir', p%dir, 'norm', norm2(p%dir)
    stop 'init_diff_photon_8_10 :: direction was not normed?!'
  end if
  select type (bmc)
  end select

  ierr = 0
end subroutine

subroutine intersect_distance_8_10(bmc, vertices, p, max_dist)
  class(t_boxmc_8_10) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist
  integer(mpiint) :: ierr
  call intersect_cube_general(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
    & p%side, p%weight, max_dist, p%subface, ierr)

  !call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
  !  p%side, max_dist)

  select type (bmc)
  end select

end subroutine
subroutine box_halfspaces_8_10(bmc, vertices, origins, normals)
  class(t_boxmc_8_10) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  real(ireal_dp), allocatable, intent(out) :: origins(:, :), normals(:, :)
  call box_halfspaces(vertices, origins, normals)
  select type (bmc)
  end select
end subroutine
! ******************* BOXMC 8 10 ************************

