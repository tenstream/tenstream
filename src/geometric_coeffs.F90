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

module m_geometric_coeffs

use m_data_parameters, only : ireals, mpiint, iintegers, zero, one, Pi
use m_helper_functions, only : pentagon_area_by_vertices, quadrangle_area_by_vertices, &
  triangle_area_by_vertices, compute_normal_3d, cross_3d, cstr, CHKERR, toStr, expm1
use m_intersection, only: hit_plane, line_intersection_3d

implicit none
private
public :: dir2dir3_geometric_coeffs

logical, parameter :: ldebug= .False.
contains

  subroutine dir2dir3_geometric_coeffs(vertices, sundir, extinction_coeff, coeffs, num_intervals)
    real(ireals), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: vertices(:), sundir(:), extinction_coeff
    integer(iintegers), intent(in), optional :: num_intervals
    real(ireals), dimension(3) :: d_p, b_p, a_p, h_p, f_p, e_p, g_p
    real(ireals) :: sun_up_down
    real(ireals), parameter :: small=tiny(extinction_coeff)

    if (extinction_coeff .lt. small) &
      & call CHKERR(1_mpiint, 'Extinction coeff too small: '//toStr(extinction_coeff)//' < min = '//toStr(small))

    associate ( &
      a => vertices( 1: 3), &
      b => vertices( 4: 6), &
      c => vertices( 7: 9), &
      d => vertices(10:12), &
      e => vertices(13:15), &
      f => vertices(16:18), &
      g => vertices(19:21), &
      h => vertices(22:24)  &
      )

    if (ldebug) then
      print *, cstr('box vertices', 'green')
      print *, 'a', a
      print *, 'b', b
      print *, 'c', c
      print *, 'd', d
      print *, 'e', e
      print *, 'f', f
      print *, 'g', g
      print *, 'h', h
      print *, '_________________________________________________________________'
      print *, cstr('sundir', 'yellow'), sundir
    endif

    if (ldebug) then
      print *, cstr('sun_up_down dot_product', 'red')
      print *, dot_product( - sundir, compute_normal_3d(e, f, h))
    endif

    sun_up_down = min(sign(one, dot_product(- sundir, compute_normal_3d(e, f, h))), zero)

    if(ldebug) then
      print *, '_________________________________________________________________'
      print *, cstr('src z', 'blue')
    endif
    call create_proj_copies(h, g, e, f, h_p, g_p, e_p, f_p)
    call project_points(sundir, d, compute_normal_3d(d, b, a), h_p, g_p, e_p, f_p)
    call rearange_projections(d, c, a, b, h_p, g_p, e_p, f_p)
    call gomtrc_coeffs( &
      d   , c   , a   , b,    g, & ! fixed
      h_p , g_p , e_p , f_p,     & ! projected
      extinction_coeff, &
      [integer(iintegers) :: 1, 7, 4], &
      [integer(iintegers) :: 2, 3, 1], &
      coeffs, & ! slice of relevant coefficients , and coefficient array
      sun_up_down, &
      num_intervals=num_intervals &
      )

    if (ldebug) then
      print *, '_________________________________________________________________'
      print *, cstr('src x', 'blue')
    endif
    call create_proj_copies(h, d, b, f, h_p, d_p, b_p, f_p)
    call project_points(sundir, a, compute_normal_3d(c, a, e), h_p, d_p, b_p, f_p)
    call rearange_projections(g, c, a, e, h_p, d_p, b_p, f_p)
    call gomtrc_coeffs( &
      g   , c   , a   , e,    d, & ! fixed
      h_p , d_p , b_p , f_p,     & ! projected
      extinction_coeff, &
      [integer(iintegers) :: 5, 8, 2], &
      [integer(iintegers) :: 2, 1, 3], &
      coeffs, & ! slice of relevant coefficients , and coefficient array
      sun_up_down, &
      num_intervals=num_intervals &
      )

    if (ldebug) then
      print *, '_________________________________________________________________'
      print *, cstr('src y', 'blue')
    endif
    call create_proj_copies(f, b, a, e, f_p, b_p, a_p, e_p)
    call project_points(sundir, d, compute_normal_3d(c, d, h), f_p, b_p, a_p, e_p)
    call rearange_projections(h, d, c, g, f_p, b_p, a_p, e_p)
    call gomtrc_coeffs( &
      h   , d   , c   , g,    b, & ! fixed
      f_p , b_p , a_p , e_p,     & ! projected
      extinction_coeff, &
      [integer(iintegers) :: 9, 6, 3], &
      [integer(iintegers) ::1, 2, 3], &
      coeffs,       & ! slice of relevant coefficients , and coefficient array
      sun_up_down, &
      num_intervals=num_intervals &
      )
  end associate

  if (ldebug) print *, '_________________________________________________________________'

  contains

    subroutine gomtrc_coeffs( &
        & f1, f2, f3, f4, f5, &
        & v1, v2, v3, v4, &
        & c_ext, &
        & slice, &
        & other_slice, &
        & coeffs, &
        & sun_up_down, &
        & num_intervals &
        )
      real(ireals), intent(in) :: sun_up_down
      real(ireals), intent(in) :: c_ext
      real(ireals), intent(in), dimension(3) :: f1, f2, f3, f4, f5, v1, v2, v3, v4
      integer(iintegers), intent(in), optional :: num_intervals
      integer(iintegers), dimension(3), intent(in) :: slice, other_slice
      real(ireals), intent(inout) :: coeffs(9)
      real(ireals) :: area1, area2, area3, sin_theta, cos_src_trgt, aq, at, s
      real(ireals), dimension(3) :: pl, pb, pt, pr, normal
      real(ireals), parameter :: small = sqrt(epsilon(one))

      if (ldebug) then
        print *, cstr('fixed', 'yellow')
        print *, 'f1', f1
        print *, 'f2', f2
        print *, 'f3', f3
        print *, 'f4', f4
        print *, '_________________________________________________________________'
      endif

      s =  norm2(hit_plane(f1, sundir, f5, compute_normal_3d(f1, f2, f3)) * sundir) ! |f1 - f1rp|
      normal = compute_normal_3d(f1, f2, f5)
      sin_theta = max(sin(abs(atan(sundir(other_slice(1)) / &
        sqrt(sundir(other_slice(2))**2 + sundir(other_slice(3))**2)))), tiny(sin_theta))
      cos_src_trgt = cos(acos(dot_product(f1 - f2, f1 - f4) / (norm2(f1 - f2) * norm2(f1 - f4))) - Pi / 2)

      if (norm2(v1-f1) .gt. small) then
        if (ldebug) print *, cstr('v1', 'red')
        call proj_var_to_edges(f1, f2, f3, f4, v1, pl, pb, pt, pr)
        area1 = quadrangle_area_by_vertices(v1, pb, f3, pl) * exp( - c_ext * s)

        s = norm2(hit_plane(v1, sundir, f1, normal) * sundir) ! prp - v1
        aq = quadrangle_area_by_vertices(v1, pr, f2, pb) * ext(s * c_ext)
        at = num_dst(s, norm2(f1 - pr) * cos_src_trgt, norm2(pr - v1), c_ext, num_intervals)
        area2 = aq + at

        normal = compute_normal_3d(f3, f2, f5)
        s = norm2(hit_plane(v1, sundir, f1, normal) * sundir) ! prp - v1
        aq = quadrangle_area_by_vertices(v1, pl, f4, pt) * ext(s * c_ext)
        at = num_dst(s, norm2(f1 - pt) * cos_src_trgt, norm2(pt - v1), c_ext, num_intervals)
      else if (norm2(v2-f2) .gt. small) then
        if (ldebug) print *, cstr('v2', 'red')
        call proj_var_to_edges(f2, f1, f4, f3, v2, pl, pt, pb, pr)
        area1 = quadrangle_area_by_vertices(v2, pl, f4, pt) * exp( - c_ext * s)

        s = norm2(hit_plane(v2, sundir, f2, normal) * sundir) ! prp - v2
        aq = quadrangle_area_by_vertices(v2, pr, f1, pt) * ext(s * c_ext)
        at = num_dst(s, norm2(f2 - pr) * cos_src_trgt, norm2(pr - v2), c_ext, num_intervals)
        area2 = aq + at

        normal = compute_normal_3d(f3, f2, f5)
        s = norm2(hit_plane(v2, sundir, f2, normal) * sundir) ! prp - v2
        aq = quadrangle_area_by_vertices(v2, pl, f3, pb) * ext(s * c_ext)
        at = num_dst(s, norm2(f2 - pb) * cos_src_trgt, norm2(pb - v2), c_ext, num_intervals)
      else if (norm2(v3-f3) .gt. small) then
        call proj_var_to_edges(f4, f3, f2, f1, v3, pr, pb, pt, pl)
        if (ldebug) print *, cstr('v3', 'red')
        area1 = quadrangle_area_by_vertices(v3, pt, f1, pr) * exp( - c_ext * s)

        s = norm2(hit_plane(v3, sundir, f3, normal) * sundir) ! prp - v3
        aq = quadrangle_area_by_vertices(v3, pl, f4, pt) * ext(s * c_ext)
        at = num_dst(s, norm2(f3 - pl) * cos_src_trgt, norm2(pl - v3), c_ext, num_intervals)
        area2 = aq + at

        normal = compute_normal_3d(f3, f2, f5)
        s = norm2(hit_plane(v3, sundir, f3, normal) * sundir) ! prp - v3
        aq = quadrangle_area_by_vertices(v3, pr, f2, pb) * ext(s * c_ext)
        at = num_dst(s, norm2(f3 - pb) * cos_src_trgt, norm2(pb - v3), c_ext, num_intervals)
      else if (norm2(v4-f4) .gt. small) then
        call proj_var_to_edges(f4, f3, f2, f1, v4, pr, pb, pt, pl)
        if (ldebug) print *, cstr('v4', 'red')
        area1 = quadrangle_area_by_vertices(v4, pr, f2, pb) * exp( - c_ext * s)

        s = norm2(hit_plane(v4, sundir, f4, normal) * sundir) ! prp - v4
        aq = quadrangle_area_by_vertices(v4, pl, f3, pb) * ext(s * c_ext)
        at = num_dst(s, norm2(f4 - pl) * cos_src_trgt, norm2(pl - v4), c_ext, num_intervals)
        area2 = aq + at

        normal = compute_normal_3d(f3, f2, f5)
        s = norm2(hit_plane(v4, sundir, f4, normal) * sundir) ! prp - v4
        aq = quadrangle_area_by_vertices(v4, pr, f1, pt) * ext(s * c_ext)
        at = num_dst(s, norm2(f4 - pt) * cos_src_trgt, norm2(pt - v4), c_ext, num_intervals)
      else
        ! initializing non initialized in case of no other fullfilled case
        area1 = quadrangle_area_by_vertices(f1, f2, f3, f4) * exp( - c_ext * s)
        area2 = zero
        aq = zero
        at = zero
      endif

      area3 = aq + at

      if (ldebug) then
        print *, 'a3q', aq
        print *, 'a3t', at
        print *, 'atot', quadrangle_area_by_vertices(f1, f2, f3, f4)
        print *, 'a1', area1
        print *, 'a2', area2
        print *, 'a3', area3
      endif


      if (ldebug) print *, 'sun_up_down', sun_up_down
      area1 = area1 - sun_up_down * area3
      area3 = area3 + sun_up_down * area3

      coeffs(slice) = [max(area1, zero), max(area2, zero), max(area3, zero)]
      coeffs(slice) = coeffs(slice) / max(quadrangle_area_by_vertices(f1, f2, f3, f4), sum(coeffs(slice)))

    end subroutine

    real(ireals) function num_dst(s0, l0, h0, c_ext, num_intervals)
      real(ireals), intent(in) :: s0, l0, h0, c_ext
      integer(iintegers), intent(in), optional :: num_intervals
      integer(iintegers) :: n
      integer(iintegers) :: i
      real(ireals) :: s, dh, ds, dl, h, j

      if (present(num_intervals)) then
        n = num_intervals
      else
        n = 10
      endif

      dl = l0 / real(n, ireals)
      dh = h0 / real(n, ireals)
      ds = s0 / real(n, ireals)

      num_dst = zero
      do i=1, n
        j = real(i, ireals) - 0.5_ireals
        s = j * ds
        h = j * dh
        num_dst = num_dst + dl * h * ext(s * c_ext)
      enddo
    end function

    real(ireals) function ext(beta_s)
      real(ireals), intent(in) :: beta_s
      real(ireals) :: x

      x = max(beta_s, tiny(x))
      ext = expm1(-x) / (-x)
    end function
  end subroutine

  subroutine create_proj_copies(p1, p2, p3, p4, c1, c2, c3, c4)
    real(ireals), intent(in), dimension(3) :: p1, p2, p3, p4
    real(ireals), intent(out), dimension(3) :: c1, c2, c3, c4

    c1 = p1
    c2 = p2
    c3 = p3
    c4 = p4
  end subroutine

  subroutine project_points(sundir, origin, normal, v1, v2, v3, v4)
    real(ireals), intent(in), dimension(3) :: sundir, normal, origin
    real(ireals), intent(inout), dimension(3) :: v1, v2, v3, v4
    real(ireals) :: sundir_proj(3)
    real(ireals), parameter :: big = 1e5_ireals !  THIS IS PROBLEMATIC IN SP

    if (ldebug) then
      print *, cstr('unprojected', 'yellow')
      print *, 'v1', v1
      print *, 'v2', v2
      print *, 'v3', v3
      print *, 'v4', v4
      print *, '_________________________________________________________________'
    endif

    if (ldebug) then
      print *, 'dot_product', dot_product(sundir, normal)
      print *, 'normal', normal
      print *, 'epsilon', epsilon(sundir)
      print *, 'sqrt', sqrt(epsilon(sundir))
      print *, 'big', big
    endif

    if (abs(dot_product(sundir, normal)) < sqrt(epsilon(sundir))) then
      sundir_proj = -(abs(normal) - [one, one, one]) * sundir
      v1 = v1 + hit_plane(v1, normal, origin, normal) * normal - big * sundir_proj
      v2 = v2 + hit_plane(v2, normal, origin, normal) * normal - big * sundir_proj
      v3 = v3 + hit_plane(v3, normal, origin, normal) * normal - big * sundir_proj
      v4 = v4 + hit_plane(v4, normal, origin, normal) * normal - big * sundir_proj
    else
      v1 = v1 + hit_plane(v1, sundir, origin, normal) * sundir
      v2 = v2 + hit_plane(v2, sundir, origin, normal) * sundir
      v3 = v3 + hit_plane(v3, sundir, origin, normal) * sundir
      v4 = v4 + hit_plane(v4, sundir, origin, normal) * sundir
    endif

    if (ldebug) then
      print *, cstr('projections', 'yellow')
      print *, 'v1', v1
      print *, 'v2', v2
      print *, 'v3', v3
      print *, 'v4', v4
      print *, '_________________________________________________________________'
    endif
  end subroutine

  subroutine rearange_point(origin, direction, coefficient, point)
    real(ireals), intent(inout) :: point(3)
    real(ireals), intent(in) :: origin(3), direction(3), coefficient

    point = origin + coefficient * direction
  end subroutine

  subroutine proj_var_to_edges( &
      c1, c2, c3, c4, & ! side's corners
      v, &              ! point to be projected
      p1, p2, p3, p4 &  ! projections
      ) ! c: corner
    real(ireals), intent(in), dimension(3) :: c1, c2, c3, c4, v
    real(ireals), intent(out), dimension(3) :: p1, p2, p3, p4
    real(ireals) :: c, t
    integer(mpiint) :: ierr

    call line_intersection_3d(v, c4-c1, c3, c4-c3, c, t, ierr)
    call rearange_point(v, c4-c1, c, p1)
    call line_intersection_3d(v, c2-c1, c2, c3-c2, c, t, ierr)
    call rearange_point(v, c2-c1, c, p2)
    call line_intersection_3d(v, c1-c2, c4, c1-c4, c, t, ierr)
    call rearange_point(v, c1-c2, c, p3)
    call line_intersection_3d(v, c4-c1, c2, c1-c2, c, t, ierr)
    call rearange_point(v, c3-c2, c, p4)
  end subroutine

  subroutine rearange_projections(f1, f2, f3, f4, v1, v2, v3, v4)
    real(ireals), intent(in), dimension(3) :: f1, f2, f3, f4
    real(ireals), intent(inout), dimension(3) :: v1, v2, v3, v4
    real(ireals), parameter :: eps=epsilon(f1)

    if (norm2(f1-v1) .gt. eps) call rearange_projection(f1-v1, f3, f4-f3, f2, f3-f2, v1)
    if (norm2(f2-v2) .gt. eps) call rearange_projection(f2-v2, f3, f4-f3, f4, f1-f4, v2)
    if (norm2(f3-v3) .gt. eps) call rearange_projection(f3-v3, f2, f1-f2, f4, f1-f4, v3)
    if (norm2(f4-v4) .gt. eps) call rearange_projection(f4-v4, f2, f1-f2, f2, f3-f2, v4)

    if (ldebug) then
      print *, '_________________________________________________________________'
      print *, cstr('rearangements', 'yellow')
      print *, 'v1', v1
      print *, 'v2', v2
      print *, 'v3', v3
      print *, 'v4', v4
      print *, '_________________________________________________________________'
    endif
  end subroutine

  subroutine rearange_projection(direction1, origin2, direction2, origin3, direction3, origin1)
    real(ireals), intent(in), dimension(3) :: direction1, origin2, direction2, origin3, direction3
    real(ireals), intent(inout), dimension(3) :: origin1
    real(ireals) :: coeff21, coeff22, coeff31, coeff32
    integer(mpiint) :: ierr1, ierr2

    call line_intersection_3d(origin1, direction1, origin2, direction2, coeff21, coeff22, ierr1)
    call line_intersection_3d(origin1, direction1, origin3, direction3, coeff31, coeff32, ierr2)
    call rearange_point(origin1, direction1, min(max(coeff21, coeff31, zero), one), origin1)
  end subroutine
end module
