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

subroutine update_dir_stream_8_12(bmc, vertices, p, T)
  class(t_boxmc_8_12) :: bmc
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
subroutine update_diff_stream_8_12(bmc, vertices, p, S)
  class(t_boxmc_8_12) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(in) :: p
  real(irealbmc), intent(inout) :: S(:)
  real(ireal_dp) :: mu

!         _____1_3_______
!        |           12  |
!      7 |           10  | 8
!        |               |
!        |   11          |
!        |   9           |
!      5 |_______________| 6
!              2 4

  if (p%side .eq. 1) then
    mu = dot_product(p%dir, [zero, zero, one])
    if (mu .lt. 0) then
      call print_photon(p)
      call CHKERR(1_mpiint, 'wrong p%side?')
    end if
    if (mu .gt. 0.5_ireal_dp) then ! angles less than 60 deg go into stream 1
      S(1) = S(1) + p%weight
    else ! larger angles into stream 3
      S(3) = S(3) + p%weight
    end if
    return

  else if (p%side .eq. 2) then
    mu = dot_product(p%dir, [zero, zero, -one])
    if (mu .lt. 0) then
      call print_photon(p)
      call CHKERR(1_mpiint, 'wrong p%side?')
    end if
    if (mu .gt. 0.5_ireal_dp) then ! angles less than 60 deg go into stream 1
      S(2) = S(2) + p%weight
    else ! larger angles into stream 3
      S(4) = S(4) + p%weight
    end if
    return

  else if (p%side .eq. 3 .and. p%dir(3) .le. zero) then
    S(5) = S(5) + p%weight
    return
  else if (p%side .eq. 3 .and. p%dir(3) .gt. zero) then
    S(7) = S(7) + p%weight
    return

  else if (p%side .eq. 4 .and. p%dir(3) .le. zero) then
    S(6) = S(6) + p%weight
    return
  else if (p%side .eq. 4 .and. p%dir(3) .gt. zero) then
    S(8) = S(8) + p%weight
    return

  else if (p%side .eq. 5 .and. p%dir(3) .le. zero) then
    S(9) = S(9) + p%weight
    return
  else if (p%side .eq. 5 .and. p%dir(3) .gt. zero) then
    S(11) = S(11) + p%weight
    return

  else if (p%side .eq. 6 .and. p%dir(3) .le. zero) then
    S(10) = S(10) + p%weight
    return
  else if (p%side .eq. 6 .and. p%dir(3) .gt. zero) then
    S(12) = S(12) + p%weight
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

subroutine init_dir_photon_8_12(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  class(t_boxmc_8_12) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers), intent(in) :: src
  logical, intent(in) :: ldirect
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp), parameter :: h = .5_ireal_dp
  real(ireal_dp), parameter :: e = 1e-6_ireal_dp
  real(ireal_dp) :: dx, dy, dz

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

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
  select case (src)
  case (1)
    p%loc = [L(dx) * h, L(dy) * h, dz * (one - e)]
  case (2)
    p%loc = [dx * h + L(dx) * h, L(dy) * h, dz * (one - e)]
  case (3)
    p%loc = [L(dx) * h, dy * h + L(dy) * h, dz * (one - e)]
  case (4)
    p%loc = [dx * h + L(dx) * h, dy * h + L(dy) * h, dz * (one - e)]
  case (5)
    p%loc = [zero + dx * e, L(dy), L(dz) * h]
  case (6)
    p%loc = [zero + dx * e, L(dy), dz * h + L(dz) * h]
  case (7)
    p%loc = [L(dx), zero + dy * e, L(dz) * h]
  case (8)
    p%loc = [L(dx), zero + dy * e, dz * h + L(dz) * h]
  case default
    print *, 'Dont know what to do with source spec:', src
    call exit
  end select

  p%weight = one
  p%alive = .true.
  p%direct = ldirect
  p%side = int(nil)
  p%src = src
  p%scattercnt = 0
  p%dir = initial_dir
  select type (bmc)
  end select

  ierr = 0
end subroutine

subroutine init_diff_photon_8_12(bmc, p, src, vertices, ierr)
  class(t_boxmc_8_12) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: vertices(:)
  integer(iintegers), intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi, mu
  real(ireal_dp) :: dx, dy, dz

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian

  select case (src)
  case (1)
    p%loc = (/L(dx), L(dy), zero/)
    phi = deg2rad(R() * 360)
    mu = sqrt(interv_R(.25_ireal_dp, one))
  case (2)
    p%loc = (/L(dx), L(dy), dz/)
    phi = deg2rad(R() * 360)
    mu = sqrt(interv_R(.25_ireal_dp, one))
  case (3)
    p%loc = (/L(dx), L(dy), zero/)
    phi = deg2rad(R() * 360)
    mu = sqrt(interv_R(zero, .25_ireal_dp))
  case (4)
    p%loc = (/L(dx), L(dy), dz/)
    phi = deg2rad(R() * 360)
    mu = sqrt(interv_R(zero, .25_ireal_dp))
  case (5)
    p%loc = (/dx, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (6)
    p%loc = (/zero, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (7)
    p%loc = (/dx, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (8)
    p%loc = (/zero, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (9)
    p%loc = (/L(dx), dy, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (10)
    p%loc = (/L(dx), zero, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (11)
    p%loc = (/L(dx), dy, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (12)
    p%loc = (/L(dx), zero, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case default
    print *, 'Dont know what to do with source spec:', src
    call exit
  end select
  select case (src)
  case (1:4)
    phi = deg2rad(R() * 360)
  case (5:12)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  end select

  p%dir = (/sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu/)

  select case (src)
  case (1, 3)
    continue
  case (2, 4)
    p%dir = rotate_angle_y(p%dir, 180 * one)
  case (5)
    p%dir = rotate_angle_y(p%dir, 90 * one)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (6)
    p%dir = rotate_angle_y(p%dir, -90 * one)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (7)
    p%dir = rotate_angle_y(p%dir, 90 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (8)
    p%dir = rotate_angle_y(p%dir, -90 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (9)
    p%dir = rotate_angle_z(p%dir, 180 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (10)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (11)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (12)
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
    stop 'init_diff_photon_8_12 :: direction was not normed?!'
  end if
  select type (bmc)
  end select

  ierr = 0
end subroutine

subroutine intersect_distance_8_12(bmc, vertices, p, max_dist)
  class(t_boxmc_8_12) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist

  call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
                      p%side, max_dist)

  select type (bmc)
  end select

end subroutine
subroutine box_halfspaces_8_12(bmc, vertices, origins, normals)
  class(t_boxmc_8_12) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  real(ireal_dp), allocatable, intent(out) :: origins(:, :), normals(:, :)
  call box_halfspaces(vertices, origins, normals)
  select type (bmc)
  end select
end subroutine
! ******************* BOXMC 8 12 ************************

