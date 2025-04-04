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

module m_boxmc_geometry
  use iso_fortran_env, only: real32, real64
  use m_data_parameters, only: mpiint, iintegers, ireals, ireal_dp, one, zero
  use m_helper_functions, only: &
    & angle_between_two_vec, &
    & approx, &
    & CHKERR, &
    & compute_normal_3d, &
    & cross_3d, &
    & determine_normal_direction, &
    & distances_to_triangle_edges, &
    & distance_to_edge, &
    & meanval, &
    & toStr, &
    & triangle_area_by_vertices, &
    & triangle_inner_circle_center

  use m_intersection, only: &
    & hit_plane, &
    & pnt_in_triangle, &
    & square_intersection, &
    & triangle_intersection

  implicit none

  integer(iintegers), parameter :: &
    & PPRTS_TOP_FACE = 1, & ! z+0
    & PPRTS_BOT_FACE = 2, & ! z+1
    & PPRTS_LEFT_FACE = 3, & ! x+0
    & PPRTS_RIGHT_FACE = 4, & ! x+1
    & PPRTS_REAR_FACE = 5, & ! y+0
    & PPRTS_FRONT_FACE = 6    ! y+1

  interface setup_default_unit_cube_geometry
    module procedure setup_default_unit_cube_geometry_r32, setup_default_unit_cube_geometry_r64
  end interface
  interface setup_default_wedge_geometry
    module procedure setup_default_wedge_geometry_r32, setup_default_wedge_geometry_r64
  end interface

contains
  ! Defining related geometric variables given the following vertex coordinates
  !
  !          G________________H
  !          /|              /|
  !         / |             / |
  !       E/__|___________F/  |
  !        |  |            |  |
  !        |  |            |  |
  !        |  |____________|__|
  !        |  /C           |  /D
  !        | /             | /
  !        |/______________|/
  !        A               B

  subroutine setup_cube_coords_from_vertices(vertices, dx, dy, dz)
    real(ireal_dp), intent(in) :: vertices(:)
    real(ireal_dp), intent(out) :: dx, dy, dz
    real(ireal_dp), dimension(3) :: A, B, C, D, E, F, G, H
    logical :: ladvanced = .false.

    if (ladvanced) then
      if (size(vertices) .eq. 4 * 2) then ! given are vertices on the top of a cube(x,y)
        A(1:2) = vertices(1:2); A(3) = 0
        B(1:2) = vertices(3:4); B(3) = 0
        C(1:2) = vertices(5:6); C(3) = 0
        D(1:2) = vertices(7:8); D(3) = 0

        dx = meanval([norm2(B - A), norm2(D - C)])
        dy = meanval([norm2(C - A), norm2(D - B)])
        dz = one
      elseif (size(vertices) .eq. 2 * 4 * 3) then ! 3D coords
        A = vertices(1:3)
        B = vertices(4:6)
        C = vertices(7:9)
        D = vertices(10:12)
        E = vertices(13:15)
        F = vertices(16:18)
        G = vertices(19:21)
        H = vertices(22:24)

        dx = meanval([norm2(B - A), norm2(D - C), norm2(F - E), norm2(H - G)])
        dy = meanval([norm2(C - A), norm2(D - B), norm2(G - E), norm2(H - F)])
        dz = meanval([norm2(E - A), norm2(F - B), norm2(H - D), norm2(G - C)])
      else
        call CHKERR(1_mpiint, 'dont know how to handle coords with '//&
          & toStr(size(vertices, kind=iintegers))//' vertex entries')
      end if
    else
      if (size(vertices) .eq. 4 * 2) then ! given are vertices on the top of a cube(x,y)
        dx = vertices(3) - vertices(1)
        dy = vertices(6) - vertices(2)
        dz = one
      elseif (size(vertices) .eq. 2 * 4 * 3) then ! 3D coords
        dx = vertices(4) - vertices(1)
        dy = vertices(8) - vertices(2)
        dz = vertices(15) - vertices(3)
      else
        call CHKERR(1_mpiint, 'dont know how to handle coords with '//&
          & toStr(size(vertices, kind=iintegers))//' vertex entries')
      end if
    end if
  end subroutine

  pure subroutine setup_default_cube_geometry(A, B, C, D, dz, vertices)
    real(ireals), intent(in) :: A(2), B(2), C(2), D(2), dz
    real(ireals), intent(out) :: vertices(24)

    vertices(1:2) = A
    vertices(4:5) = B
    vertices(7:8) = C
    vertices(10:11) = D
    vertices([3, 6, 9, 12]) = zero

    vertices(13:14) = A
    vertices(16:17) = B
    vertices(19:20) = C
    vertices(22:23) = D
    vertices([15, 18, 21, 24]) = dz
  end subroutine
  pure subroutine setup_default_unit_cube_geometry_r32(dx, dy, dz, vertices)
    real(real32), intent(in) :: dx, dy, dz
    real(real32), intent(out) :: vertices(24)
    real(kind=kind(dx)), parameter :: zero = 0

    vertices(1:2) = [zero, zero]
    vertices(4:5) = [dx, zero]
    vertices(7:8) = [zero, dy]
    vertices(10:11) = [dx, dy]
    vertices([3, 6, 9, 12]) = zero

    vertices(13:14) = [zero, zero]
    vertices(16:17) = [dx, zero]
    vertices(19:20) = [zero, dy]
    vertices(22:23) = [dx, dy]
    vertices([15, 18, 21, 24]) = dz
  end subroutine
  pure subroutine setup_default_unit_cube_geometry_r64(dx, dy, dz, vertices)
    real(real64), intent(in) :: dx, dy, dz
    real(real64), intent(out) :: vertices(24)
    real(kind=kind(dx)), parameter :: zero = 0

    vertices(1:2) = [zero, zero]
    vertices(4:5) = [dx, zero]
    vertices(7:8) = [zero, dy]
    vertices(10:11) = [dx, dy]
    vertices([3, 6, 9, 12]) = zero

    vertices(13:14) = [zero, zero]
    vertices(16:17) = [dx, zero]
    vertices(19:20) = [zero, dy]
    vertices(22:23) = [dx, dy]
    vertices([15, 18, 21, 24]) = dz
  end subroutine

  ! Distribution Code for Wedges:
  !
  !        F
  !       / \
  !    dy/   \
  !     /     \
  !    /   dx  \
  !   D________ E
  !   |         |
  !   |    |    |
  !   |    |    |
  !   |    |    |
  !   |         |
  !   |    C    |
  !   |   / \   |
  !   |  /   \  |
  !   | /     \ |
  !   |/       \|
  !   A _______ B

  ! We always assume the triangle to have dx edge length along y=0 (A,B) and dy is the edge length between (A,C)
  ! Distribute Photons on triangles: https://doi.org/10.1145/571647.571648

  subroutine setup_wedge_coords_from_vertices(vertices, A, B, C, nAB, nBC, nCA, dx, dy, dz)
    real(ireal_dp), intent(in) :: vertices(:) ! should be the vertex coordinates for A, B, C, D, E, F in 3D
    real(ireal_dp), dimension(:), intent(out) :: A, B, C ! points on triangle [A,B,C]
    real(ireal_dp), dimension(:), intent(out) :: nAB, nBC, nCA ! normals on triangle [A,B,C], pointing towards center
    real(ireal_dp), intent(out) :: dx, dy, dz

    real(ireal_dp), dimension(2) :: D, E, F ! points on triangle above [A,B,C]
    integer :: i

    if (size(A) .ne. 2 .or. size(B) .ne. 2 .or. size(C) .ne. 2) call CHKERR(1_mpiint, 'Coordinates have to be 2D coordinates')
    if (size(nAB) .ne. 2 .or. size(nBC) .ne. 2 .or. size(nCA) .ne. 2) call CHKERR(1_mpiint, 'Normals have to be 2D coordinates')
    if (size(vertices) .ne. 2 * 3 * 3) call CHKERR(1_mpiint, 'Size of vertices have to be 3D coordinates for 6 points')

    A = vertices(1:2)
    B = vertices(4:5)
    C = vertices(7:8)

    D = vertices(10:11)
    E = vertices(13:14)
    F = vertices(16:17)

    nAB = (A - B); nAB = nAB([2, 1]); nAB = nAB*[one, -one] / norm2(nAB)
    nBC = (B - C); nBC = nBC([2, 1]); nBC = nBC*[one, -one] / norm2(nBC)
    nCA = (C - A); nCA = nCA([2, 1]); nCA = nCA*[one, -one] / norm2(nCA)

    dx = norm2(B - A)
    dy = norm2(C - A)
    dz = vertices(12) - vertices(3)
    if (any(approx([dx, dy, dz], 0._ireal_dp))) then
      do i = 1, 6
        print *, 'vertices', i, vertices((i - 1) * 3 + 1:i * 3)
      end do
      print *, 'A', A
      print *, 'B', B
      print *, 'C', C
      print *, 'D', D
      print *, 'E', E
      print *, 'F', F
      print *, 'nA', nAB
      print *, 'nB', nBC
      print *, 'nC', nCA
      print *, 'dx/dy/dz', dx, dy, dz
      call CHKERR(1_mpiint, 'bad differential box lengths')
    end if
  end subroutine

  subroutine setup_default_wedge_geometry_r32(A, B, C, dz, vertices, sphere_radius)
    real(real32), intent(in) :: A(2), B(2), C(2), dz
    real(real32), intent(out) :: vertices(18)
    real(real32), intent(in), optional :: sphere_radius
    real(real32) :: s, center(3)

    vertices(1:2) = A
    vertices(4:5) = B
    vertices(7:8) = C
    vertices([3, 6, 9]) = zero

    vertices(10:11) = A
    vertices(13:14) = B
    vertices(16:17) = C
    vertices([12, 15, 18]) = dz

    !print *,'A', vertices(1:3)
    !print *,'B', vertices(4:6)
    !print *,'C', vertices(7:9)
    !print *,'D', vertices(10:12)
    !print *,'E', vertices(13:15)
    !print *,'F', vertices(16:18)

    if (present(sphere_radius)) then
      if (sphere_radius .gt. zero) then
        s = (sphere_radius + dz) / sphere_radius
        associate (D => vertices(10:12), E => vertices(13:15), F => vertices(16:18))
          call triangle_inner_circle_center(D, E, F, center)

          D = center + s * (D - center)
          E = center + s * (E - center)
          F = center + s * (F - center)
        end associate
        !print *,'center', center, 's', s
        !print *,'sA', vertices(1:3)
        !print *,'sB', vertices(4:6)
        !print *,'sC', vertices(7:9)
        !print *,'sD', vertices(10:12)
        !print *,'sE', vertices(13:15)
        !print *,'sF', vertices(16:18)
      end if
    end if
  end subroutine
  subroutine setup_default_wedge_geometry_r64(A, B, C, dz, vertices, sphere_radius)
    real(real64), intent(in) :: A(2), B(2), C(2), dz
    real(real64), intent(out) :: vertices(18)
    real(real64), intent(in), optional :: sphere_radius
    real(real64) :: s, center(3)

    vertices(1:2) = A
    vertices(4:5) = B
    vertices(7:8) = C
    vertices([3, 6, 9]) = zero

    vertices(10:11) = A
    vertices(13:14) = B
    vertices(16:17) = C
    vertices([12, 15, 18]) = dz

    !print *,'A', vertices(1:3)
    !print *,'B', vertices(4:6)
    !print *,'C', vertices(7:9)
    !print *,'D', vertices(10:12)
    !print *,'E', vertices(13:15)
    !print *,'F', vertices(16:18)

    if (present(sphere_radius)) then
      if (sphere_radius .gt. zero) then
        s = (sphere_radius + dz) / sphere_radius
        associate (D => vertices(10:12), E => vertices(13:15), F => vertices(16:18))
          call triangle_inner_circle_center(D, E, F, center)
          D = center + s * (D - center)
          E = center + s * (E - center)
          F = center + s * (F - center)
        end associate
        !print *,'center', center, 's', s
        !print *,'sA', vertices(1:3)
        !print *,'sB', vertices(4:6)
        !print *,'sC', vertices(7:9)
        !print *,'sD', vertices(10:12)
        !print *,'sE', vertices(13:15)
        !print *,'sF', vertices(16:18)
      end if
    end if
  end subroutine
  subroutine setup_default_unit_wedge_geometry(dx, dy, dz, vertices)
    real(ireal_dp), intent(in) :: dx, dy, dz
    real(ireal_dp), intent(out) :: vertices(18)
    real(ireal_dp), parameter :: zero = 0

    vertices(1:2) = [zero, zero]
    vertices(4:5) = [dx, zero]
    vertices(7:8) = [dx / 2, sqrt(dy**2 - (dx / 2)**2)]
    vertices([3, 6, 9]) = zero

    vertices(10:11) = [zero, zero]
    vertices(13:14) = [dx, zero]
    vertices(16:17) = [dx / 2, sqrt(dy**2 - (dx / 2)**2)]
    vertices([12, 15, 18]) = dz
  end subroutine

  subroutine intersect_cube(vertices, ploc, pdir, pscattercnt, psrc_side, &
                            pside, max_dist)
    real(ireal_dp), intent(in) :: vertices(:)
    real(ireal_dp), intent(in) :: pdir(:), ploc(:)
    integer(iintegers), intent(in) :: pscattercnt, psrc_side
    integer(iintegers), intent(inout) :: pside
    real(ireal_dp), intent(out) :: max_dist

    real(ireal_dp) :: x, y, z, dx, dy, dz
    integer(iintegers) :: i, sides(3)

    real(ireal_dp) :: dist(3)
    real(ireal_dp), parameter :: zero = 0, one = 1

    real(ireal_dp), parameter :: nx(3) = [one, zero, zero]
    real(ireal_dp), parameter :: ny(3) = [zero, one, zero]
    real(ireal_dp), parameter :: nz(3) = [zero, zero, one]
    real(ireal_dp), parameter :: o0(3) = [zero, zero, zero]

    call setup_cube_coords_from_vertices(vertices, dx, dy, dz)

    !crossing with bottom and top plane:
    if (pdir(3) .ge. zero) then
      max_dist = hit_plane(ploc, pdir, [zero, zero, dz], nz)
      pside = 1
      x = ploc(1) + pdir(1) * max_dist
      y = ploc(2) + pdir(2) * max_dist
      if ((x .gt. zero .and. x .lt. dx) .and. (y .gt. zero .and. y .lt. dy)) return
      dist(1) = max_dist; sides(1) = 1
    end if
    if (pdir(3) .le. zero) then
      max_dist = hit_plane(ploc, pdir, o0, nz)
      pside = 2
      x = ploc(1) + pdir(1) * max_dist
      y = ploc(2) + pdir(2) * max_dist
      if ((x .gt. zero .and. x .lt. dx) .and. (y .gt. zero .and. y .lt. dy)) return
      dist(1) = max_dist; sides(1) = 2
    end if

    !crossing with left and right plane:
    if (pdir(1) .le. zero) then
      max_dist = hit_plane(ploc, pdir, o0, nx)
      pside = 3
      y = ploc(2) + pdir(2) * max_dist
      z = ploc(3) + pdir(3) * max_dist
      if ((y .gt. zero .and. y .lt. dy) .and. (z .gt. zero .and. z .lt. dz)) return
      dist(2) = max_dist; sides(2) = 3
    end if
    if (pdir(1) .ge. zero) then
      max_dist = hit_plane(ploc, pdir, [dx, zero, zero], nx)
      pside = 4
      y = ploc(2) + pdir(2) * max_dist
      z = ploc(3) + pdir(3) * max_dist
      if ((y .gt. zero .and. y .lt. dy) .and. (z .gt. zero .and. z .lt. dz)) return
      dist(2) = max_dist; sides(2) = 4
    end if

    !crossing with back and forward plane:
    if (pdir(2) .le. zero) then
      max_dist = hit_plane(ploc, pdir, o0, ny)
      pside = 5
      x = ploc(1) + pdir(1) * max_dist
      z = ploc(3) + pdir(3) * max_dist
      if ((x .gt. zero .and. x .lt. dx) .and. (z .gt. zero .and. z .lt. dz)) return
      dist(3) = max_dist; sides(3) = 5
    end if
    if (pdir(2) .ge. zero) then
      max_dist = hit_plane(ploc, pdir, [zero, dy, zero], ny)
      pside = 6
      x = ploc(1) + pdir(1) * max_dist
      z = ploc(3) + pdir(3) * max_dist
      if ((x .gt. zero .and. x .lt. dx) .and. (z .gt. zero .and. z .lt. dz)) return
      dist(3) = max_dist; sides(3) = 6
    end if

    !Ohhh there was a problem.. maybe with numerics, seems that it may happen that we dont find a solution if norm of pdir is not equal to one....
    max_dist = huge(max_dist)
    do i = 1, 3
      if (.not. approx(pdir(i), zero)) then
        if (pscattercnt .eq. 0 .and. pside .eq. psrc_side) cycle
        if (dist(i) .le. max_dist) then
          pside = sides(i)
          max_dist = dist(i)
        end if
      end if
    end do

    if (max_dist .gt. norm2([dx, dy, dz])) then
      print *, 'should actually not be here at the end of crossings in intersect_cube! '// &
        '- however, please check if distance makes sense?:', &
        max_dist, norm2([dx, dy, dz]), '::', dist, ':', vertices, &
        'pdir', pdir, 'ploc side', ploc, psrc_side, 'target_side', pside
      call CHKERR(1_mpiint, 'DEBUG')
    end if

  end subroutine

  !> @brief: computes intersection of a distorted cube,
  !> i.e. the generalized version of intersect_cube
  subroutine intersect_cube_general(vertices, ploc, pdir, pscattercnt, psrc_side, &
                                    pside, pweight, max_dist, psubface, ierr)
    real(ireal_dp), intent(in) :: vertices(:)
    real(ireal_dp), intent(in) :: pdir(:), ploc(:)
    integer(iintegers), intent(in) :: pscattercnt, psrc_side
    integer(iintegers), intent(inout) :: pside
    real(ireal_dp), intent(inout) :: pweight
    real(ireal_dp), intent(out) :: max_dist
    integer(iintegers), intent(out) :: psubface
    integer(mpiint), intent(out) :: ierr

    logical :: lhit(6)
    real(ireal_dp) :: hit(6, 4)
    integer(iintegers) :: i, iface(6)

    ierr = 0

    associate ( &
      A => vertices(1:3), &
      B => vertices(4:6), &
      C => vertices(7:9), &
      D => vertices(10:12), &
      E => vertices(13:15), &
      F => vertices(16:18), &
      G => vertices(19:21), &
      H => vertices(22:24))

      lhit = .false.
      hit = huge(hit)
      iface(:) = -1
      !crossing with bottom and top plane:
      call square_intersection(ploc, pdir, E, F, H, G, lhit(PPRTS_TOP_FACE), hit(PPRTS_TOP_FACE, :), iface(PPRTS_TOP_FACE))
      call square_intersection(ploc, pdir, A, C, D, B, lhit(PPRTS_BOT_FACE), hit(PPRTS_BOT_FACE, :), iface(PPRTS_BOT_FACE))

      !crossing with side planes:
      call square_intersection(ploc, pdir, A, E, G, C, lhit(PPRTS_LEFT_FACE), hit(PPRTS_LEFT_FACE, :), iface(PPRTS_LEFT_FACE))
      call square_intersection(ploc, pdir, B, D, H, F, lhit(PPRTS_RIGHT_FACE), hit(PPRTS_RIGHT_FACE, :), iface(PPRTS_RIGHT_FACE))
      call square_intersection(ploc, pdir, A, B, F, E, lhit(PPRTS_REAR_FACE), hit(PPRTS_REAR_FACE, :), iface(PPRTS_REAR_FACE))
      call square_intersection(ploc, pdir, C, G, H, D, lhit(PPRTS_FRONT_FACE), hit(PPRTS_FRONT_FACE, :), iface(PPRTS_FRONT_FACE))

      pside = 0
      max_dist = huge(max_dist)
      do i = 1, 6
        if (hit(i, 4) .lt. zero) cycle
        if (pscattercnt .eq. 0 .and. i .eq. psrc_side) cycle
        if (hit(i, 4) .lt. max_dist) then
          max_dist = hit(i, 4)
          pside = i
          psubface = iface(i)
        end if
      end do

      ! If we did not hit anything else, I assume that we point towards the src side.
      ! We collect it there but set energy to 0
      ! if(pscattercnt.eq.0 .and. pside.eq.psrc_side .and. lhit(psrc_side) ) then
      !   max_dist = hit(psrc_side,4)
      !   pside = psrc_side
      !   pweight = zero ! we dont allow energy to hit the src face, at least not right after it started!
      !   psubface = iface(psrc_side)
      ! endif

      !print *,pside,'hit',hit(pside,1:3)

      if (count(lhit) .eq. 0 .or. pside .eq. 0) then
        print *, 'should actually not be here at the end of crossings in intersect distance!'
        print *, 'max dist, pside', max_dist, pside, 'src_side', psrc_side
        print *, 'ploc', ploc
        print *, 'pdir', pdir
        print *, 'pwgt', pweight
        print *, 'A', A
        print *, 'B', B
        print *, 'C', C
        print *, 'D', D
        print *, 'E', E
        print *, 'F', F
        print *, 'G', G
        print *, 'H', H
        print *, 'lhit', lhit
        print *, 'hit1', hit(1, :)
        print *, 'hit2', hit(2, :)
        print *, 'hit3', hit(3, :)
        print *, 'hit4', hit(4, :)
        print *, 'hit5', hit(5, :)
        print *, 'hit6', hit(6, :)
        ierr = 1_mpiint
        call CHKERR(1_mpiint, 'ERROR in Raytracer, didnt hit anything!')
      end if
    end associate
  end subroutine

  subroutine intersect_wedge(vertices, ploc, pdir, pscattercnt, psrc_side, &
                             pside, pweight, max_dist, psubface, ierr)
    real(ireal_dp), intent(in) :: vertices(:)
    real(ireal_dp), intent(in) :: pdir(:), ploc(:)
    integer(iintegers), intent(in) :: pscattercnt, psrc_side
    integer(iintegers), intent(inout) :: pside
    real(ireal_dp), intent(inout) :: pweight
    real(ireal_dp), intent(out) :: max_dist
    integer(iintegers), intent(out) :: psubface
    integer(mpiint), intent(out) :: ierr

    real(ireal_dp), parameter :: rng(2) = [0._ireal_dp, huge(rng)]
    logical :: l_in_triangle
    logical :: lhit(5)
    real(ireal_dp) :: hit(4, 5)
    integer(iintegers) :: i, iface(5)
    ierr = 0

    associate ( &
      Ab => vertices(1:3), &
      Bb => vertices(4:6), &
      Cb => vertices(7:9), &
      At => vertices(10:12), &
      Bt => vertices(13:15), &
      Ct => vertices(16:18))

      lhit = .false.
      hit = huge(hit)
      iface = [1, -1, -1, -1, 1]
      !crossing with bottom and top plane:
      if (pdir(3) .ge. zero) then
        call triangle_intersection(ploc, pdir, At, Bt, Ct, rng, lhit(1), hit(:, 1))
        lhit(5) = .false.
      end if
      if (pdir(3) .le. zero) then
        call triangle_intersection(ploc, pdir, Ab, Bb, Cb, rng, lhit(5), hit(:, 5))
        lhit(1) = .false.
      end if

      !crossing with side planes:
      ! plane 2, along y=0
      call square_intersection(ploc, pdir, Ab, Bb, Bt, At, lhit(2), hit(:, 2), iface(2))
      call square_intersection(ploc, pdir, Ab, Cb, Ct, At, lhit(3), hit(:, 3), iface(3))
      call square_intersection(ploc, pdir, Bb, Cb, Ct, Bt, lhit(4), hit(:, 4), iface(4))

      pside = 0
      max_dist = huge(max_dist)
      do i = 1, 5
        if (hit(4, i) .lt. zero) cycle
        if (pscattercnt .eq. 0 .and. i .eq. psrc_side) cycle
        if (hit(4, i) .lt. max_dist) then
          max_dist = hit(4, i)
          pside = i
          psubface = iface(i)
        end if
      end do

      ! If we did not hit anything else, I assume that we point towards the src side.
      ! We collect it there but set energy to 0
      if (pscattercnt .eq. 0 .and. pside .eq. psrc_side .and. lhit(psrc_side)) then
        max_dist = hit(4, psrc_side)
        pside = psrc_side
        pweight = zero ! we dont allow energy to hit the src face, at least not right after it started!
        psubface = iface(psrc_side)
      end if

      !print *,pside,'hit',hit(1:3, pside)

      if (count(lhit) .eq. 0) then
        print *, 'should actually not be here at the end of crossings in intersect distance!'
        print *, 'max dist, pside', max_dist, pside, 'src_side', psrc_side
        print *, 'ploc', ploc
        print *, 'pdir', pdir
        print *, 'At', At
        print *, 'Bt', Bt
        print *, 'Ct', Ct
        print *, 'Ab', Ab
        print *, 'Bb', Bb
        print *, 'Cb', Cb
        print *, 'lhit', lhit
        print *, 'hit1', hit(:, 1)
        print *, 'hit2', hit(:, 2)
        print *, 'hit3', hit(:, 3)
        print *, 'hit4', hit(:, 4)
        print *, 'hit5', hit(:, 5)
        ierr = 1_mpiint
        !call CHKERR(1_mpiint, 'ERROR in Raytracer, didnt hit anything!')
      end if

      select case (pside)
      case (1)
        l_in_triangle = pnt_in_triangle(At, Bt, Ct, hit(1:2, pside))
      case (5)
        l_in_triangle = pnt_in_triangle(Ab, Bb, Cb, hit(1:2, pside))
      case default
        l_in_triangle = .true.
      end select
      if (.not. l_in_triangle) then
        print *, 'max dist, pside', max_dist, pside, 'src_side', psrc_side
        print *, 'scattercnt', pscattercnt
        print *, 'ploc', ploc
        print *, 'pdir', pdir
        print *, 'At', At
        print *, 'Bt', Bt
        print *, 'Ct', Ct
        print *, 'Ab', Ab
        print *, 'Bb', Bb
        print *, 'Cb', Cb
        print *, 'lhit', lhit
        print *, 'hit1', hit(:, 1)
        print *, 'hit2', hit(:, 2)
        print *, 'hit3', hit(:, 3)
        print *, 'hit4', hit(:, 4)
        print *, 'hit5', hit(:, 5)
        print *, 'target point not in triangle', hit, 'side', pside, 'dist', hit(4, pside)
        ierr = 2_mpiint
        !call CHKERR(1_mpiint, 'Photon not inside the triangle')
      end if
    end associate
  end subroutine

  subroutine box_halfspaces(vertices, origins, normals)
    real(ireal_dp), intent(in) :: vertices(:)
    real(ireal_dp), intent(out) :: origins(3, 6), normals(3, 6) ! size 3,6

    if (size(vertices) .ne. 2 * 4 * 3) call CHKERR(1_mpiint, 'did not expect that')

    associate ( &
      A => vertices(1:3), &
      B => vertices(4:6), &
      C => vertices(7:9), &
      D => vertices(10:12), &
      E => vertices(13:15), &
      F => vertices(16:18), &
      G => vertices(19:21), &
      H => vertices(22:24))

      origins(:, 1) = H
      origins(:, 2) = A
      origins(:, 3) = A
      origins(:, 4) = H
      origins(:, 5) = A
      origins(:, 6) = H
      normals(:, 1) = compute_normal_3d(H, F, G)
      normals(:, 2) = compute_normal_3d(A, B, C)
      normals(:, 3) = compute_normal_3d(A, E, B)
      normals(:, 4) = compute_normal_3d(H, G, D)
      normals(:, 5) = compute_normal_3d(A, C, E)
      normals(:, 6) = compute_normal_3d(H, D, F)
    end associate
  end subroutine

  subroutine wedge_halfspaces(vertices, origins, normals)
    real(ireal_dp), intent(in) :: vertices(:)
    real(ireal_dp), intent(out) :: origins(:, :), normals(:, :) ! size 3,5

    associate ( &
      A => vertices(1:3), &
      B => vertices(4:6), &
      C => vertices(7:9), &
      D => vertices(10:12), &
      E => vertices(13:15), &
      F => vertices(16:18))

      origins(:, 1) = E
      origins(:, 2) = E
      origins(:, 3) = A
      origins(:, 4) = E
      origins(:, 5) = A
      normals(:, 1) = compute_normal_3d(E, D, F)
      normals(:, 2) = compute_normal_3d(E, B, D)
      normals(:, 3) = compute_normal_3d(A, C, D)
      normals(:, 4) = compute_normal_3d(E, F, B)
      normals(:, 5) = compute_normal_3d(A, B, C)
    end associate
  end subroutine

  ! Distribute Photons on triangles: https://doi.org/10.1145/571647.571648
  subroutine rand_pnt_on_triangle(A, B, C, pnt, eps)
    real(ireal_dp), dimension(3), intent(in) :: A, B, C
    real(ireal_dp), dimension(3), intent(out) :: pnt
    real(ireal_dp) :: r1, r2
    real(ireal_dp), dimension(3) :: center, cA, cB, cC
    real(ireal_dp), intent(in), optional :: eps ! move all points a wee bit towards the center

    call random_number(r1)
    call random_number(r2)

    if (present(eps)) then
      call triangle_inner_circle_center(A, B, C, center)
      cA = A + (center - A) * eps
      cB = B + (center - B) * eps
      cC = C + (center - C) * eps
      pnt = (one - sqrt(r1)) * cA + sqrt(r1) * (one - r2) * cB + sqrt(r1) * r2 * cC
    else
      pnt = (one - sqrt(r1)) * A + sqrt(r1) * (one - r2) * B + sqrt(r1) * r2 * C
    end if
  end subroutine

  subroutine rand_pnt_on_plane(A, B, C, D, pnt, normal, U, V)
    real(ireal_dp), dimension(3), intent(in) :: A, B, C, D
    real(ireal_dp), dimension(3), intent(out) :: pnt, normal, U, V
    real(ireal_dp) :: r, area(2)

    area(1) = triangle_area_by_vertices(A, B, C)
    area(2) = triangle_area_by_vertices(A, C, D)

    call random_number(r)
    if (r .lt. area(1) / sum(area)) then
      call rand_pnt_on_triangle(A, B, C, pnt)
      normal = compute_normal_3d(A, B, C)
      U = (C - B)
      U = U / norm2(U)
      V = -cross_3d(U, normal)
    else
      call rand_pnt_on_triangle(A, C, D, pnt)
      normal = compute_normal_3d(A, C, D)
      U = (D - A)
      U = U / norm2(U)
      V = -cross_3d(U, normal)
    end if
  end subroutine

end module
