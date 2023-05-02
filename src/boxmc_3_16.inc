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

subroutine update_dir_stream_3_16(bmc, vertices, p, T)
  class(t_boxmc_3_16) :: bmc
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

subroutine update_diff_stream_3_16(bmc, vertices, p, S)
  class(t_boxmc_3_16) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(in) :: p
  real(irealbmc), intent(inout) :: S(:)
  integer(iintegers), parameter :: soff = 8 ! offset for sideward streams

  if (p%side .eq. 1) then
    if (abs(p%dir(2)) .gt. abs(p%dir(1))) then
      if (p%dir(2) .gt. zero) then
        S(1) = S(1) + p%weight
      else
        S(5) = S(5) + p%weight
      end if
    else
      if (p%dir(1) .gt. zero) then
        S(7) = S(7) + p%weight
      else
        S(3) = S(3) + p%weight
      end if
    end if

  else if (p%side .eq. 2) then
    if (abs(p%dir(2)) .gt. abs(p%dir(1))) then
      if (p%dir(2) .gt. zero) then
        S(2) = S(2) + p%weight
      else
        S(6) = S(6) + p%weight
      end if
    else
      if (p%dir(1) .gt. zero) then
        S(8) = S(8) + p%weight
      else
        S(4) = S(4) + p%weight
      end if
    end if

  else if (p%side .eq. 3 .and. p%dir(3) .le. zero) then
    S(soff + 1) = S(soff + 1) + p%weight

  else if (p%side .eq. 3 .and. p%dir(3) .gt. zero) then
    S(soff + 3) = S(soff + 3) + p%weight

  else if (p%side .eq. 4 .and. p%dir(3) .le. zero) then
    S(soff + 2) = S(soff + 2) + p%weight
  else if (p%side .eq. 4 .and. p%dir(3) .gt. zero) then
    S(soff + 4) = S(soff + 4) + p%weight

  else if (p%side .eq. 5 .and. p%dir(3) .le. zero) then
    S(soff + 5) = S(soff + 5) + p%weight

  else if (p%side .eq. 5 .and. p%dir(3) .gt. zero) then
    S(soff + 7) = S(soff + 7) + p%weight

  else if (p%side .eq. 6 .and. p%dir(3) .le. zero) then
    S(soff + 6) = S(soff + 6) + p%weight

  else if (p%side .eq. 6 .and. p%dir(3) .gt. zero) then
    S(soff + 8) = S(soff + 8) + p%weight

  else
    print *, 'Couldnt find a stream on to which I can put the photon weight on?!'
    call print_photon(p)
    call exit
    print *, 'remove unused compiler warning', vertices
  end if

  select type (bmc)
  end select

end subroutine

subroutine init_dir_photon_3_16(bmc, p, src, ldirect, initial_dir, vertices, ierr)
  class(t_boxmc_3_16) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: initial_dir(:), vertices(:)
  integer(iintegers), intent(in) :: src
  logical, intent(in) :: ldirect
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp), parameter :: e = 1e-6_ireal_dp
  real(ireal_dp) :: dx, dy, dz

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  p%dir = initial_dir
  p%alive = .false.

  select case (src)
  case (1)
    p%loc = [R() * dx, R() * dy, dz * (one - e)]
  case (2)
    if (p%dir(1) .ge. 0) then
      p%loc = [e, R() * dy, R() * dz]
    else
      p%loc = [dx - e, R() * dy, R() * dz]
    end if
  case (3)
    if (p%dir(2) .ge. 0) then
      p%loc = [R() * dx, e, R() * dz]
    else
      p%loc = [R() * dx, dy - e, R() * dz]
    end if
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
  select type (bmc)
  end select

  ierr = 0
end subroutine

subroutine init_diff_photon_3_16(bmc, p, src, vertices, ierr)
  class(t_boxmc_3_16) :: bmc
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(in) :: vertices(:)
  integer(iintegers), intent(in) :: src
  integer(mpiint), intent(out) :: ierr

  real(ireal_dp) :: phi, mu
  real(ireal_dp) :: dx, dy, dz
  integer(iintegers), parameter :: soff = 8

  call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

  mu = sqrt(R()) ! distribution between 0 and 90 degree with more towards 90 , i.e. mu=1 --> Lambertian

  select case (src)
  case (1)
    p%loc = (/L(dx), L(dy), zero/)
    phi = deg2rad(interv_R(-45._ireal_dp, 45._ireal_dp))
  case (2)
    p%loc = (/L(dx), L(dy), dz/)
    phi = deg2rad(interv_R(-45._ireal_dp, 45._ireal_dp))
  case (3)
    p%loc = (/L(dx), L(dy), zero/)
    phi = deg2rad(interv_R(225._ireal_dp, 315._ireal_dp))
  case (4)
    p%loc = (/L(dx), L(dy), dz/)
    phi = deg2rad(interv_R(225._ireal_dp, 315._ireal_dp))
  case (5)
    p%loc = (/L(dx), L(dy), zero/)
    phi = deg2rad(interv_R(135._ireal_dp, 225._ireal_dp))
  case (6)
    p%loc = (/L(dx), L(dy), dz/)
    phi = deg2rad(interv_R(135._ireal_dp, 225._ireal_dp))
  case (7)
    p%loc = (/L(dx), L(dy), zero/)
    phi = deg2rad(interv_R(45._ireal_dp, 135._ireal_dp))
  case (8)
    p%loc = (/L(dx), L(dy), dz/)
    phi = deg2rad(interv_R(45._ireal_dp, 135._ireal_dp))
  case (soff + 1)
    p%loc = (/dx, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (soff + 2)
    p%loc = (/zero, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (soff + 3)
    p%loc = (/dx, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (soff + 4)
    p%loc = (/zero, L(dy), L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (soff + 5)
    p%loc = (/L(dx), dy, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (soff + 6)
    p%loc = (/L(dx), zero, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (soff + 7)
    p%loc = (/L(dx), dy, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case (soff + 8)
    p%loc = (/L(dx), zero, L(dz)/)
    phi = deg2rad(interv_R(-90._ireal_dp, 90._ireal_dp))
  case default
    print *, 'Dont know what to do with source spec:', src
    call exit
  end select

  ! phi 0(north) 90(east) 180(south) 270(west)
  p%dir = (/sin(phi) * sin(acos(mu)), cos(phi) * sin(acos(mu)), mu/)

  select case (src)
  case (1, 3, 5, 7)
    continue
  case (2, 4, 6, 8)
    p%dir(3) = -p%dir(3)
  case (soff + 1)
    p%dir = rotate_angle_y(p%dir, 90 * one)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (soff + 2)
    p%dir = rotate_angle_y(p%dir, -90 * one)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (soff + 3)
    p%dir = rotate_angle_y(p%dir, 90 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (soff + 4)
    p%dir = rotate_angle_y(p%dir, -90 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (soff + 5)
    p%dir = rotate_angle_z(p%dir, 180 * one)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (soff + 6)
    p%dir = rotate_angle_x(p%dir, 90 * one)
  case (soff + 7)
    p%dir = rotate_angle_x(p%dir, -90 * one)
  case (soff + 8)
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
    stop 'init_diff_photon_3_16 :: direction was not normed?!'
  end if
  select type (bmc)
  end select

  ierr = 0
end subroutine

subroutine intersect_distance_3_16(bmc, vertices, p, max_dist)
  class(t_boxmc_3_16) :: bmc
  real(ireal_dp), intent(in) :: vertices(:)
  type(t_photon), intent(inout) :: p
  real(ireal_dp), intent(out) :: max_dist

  call intersect_cube(vertices, p%loc, p%dir, p%scattercnt, p%src_side, &
                      p%side, max_dist)

  select type (bmc)
  end select

end subroutine
! ******************* BOXMC 3 16 ************************
