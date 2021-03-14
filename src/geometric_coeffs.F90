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

use m_data_parameters, only : irealLUT, ireals, mpiint, iintegers, zero, one
use m_helper_functions, only : pentagon_area_by_vertices, quadrangle_area_by_vertices, &
  triangle_area_by_vertices, compute_normal_3d
use m_intersection, only: hit_plane, line_intersection_3d

implicit none
private
public :: dir2dir3_geometric_coeff_corr

logical, parameter :: lDEBUG_geometric_coeff_correction = .true.
contains

  subroutine dir2dir3_geometric_coeff_corr(verts, sundir, optical_props, coeffs)
    real(irealLUT), intent(inout) :: coeffs(9)
    real(ireals), intent(in) :: verts(24), sundir(3), optical_props(3)
    real(ireals), dimension(3) :: d_p, b_p, a_p, h_p, f_p, e_p, g_p

    associate ( &
      a => verts( 1: 3), &
      b => verts( 4: 6), &
      c => verts( 7: 9), &
      d => verts(10:12), &
      e => verts(13:15), &
      f => verts(16:18), &
      g => verts(19:21), &
      h => verts(22:24)  &
      )

    print *, 'a', a
    print *, 'b', b
    print *, 'c', c
    print *, 'd', d
    print *, 'e', e
    print *, 'f', f
    print *, 'g', g
    print *, 'h', h

    if (lDEBUG_geometric_coeff_correction) print *, 'sundir', sundir

    coeffs=0._irealLUT

    if(lDEBUG_geometric_coeff_correction) print *, 'src z'
    call create_proj_copies(h, g, e, f, h_p, g_p, e_p, f_p)
    call project_points(sundir, d, compute_normal_3d(d, b, a), h_p, g_p, e_p, f_p)
    call rearange_projections(d, c, a, b, h_p, g_p, e_p, f_p)
    call correct_coeffs( &
      d   , c   , a   , b,    g, e, & ! fixed
      h_p , g_p , e_p , f_p,     & ! projected
      optical_props(1), [1      , 7       , 4], [2, 3, 1]      , coeffs       & ! slice of relevant coefficients , and coefficient array
      )

    if (lDEBUG_geometric_coeff_correction)  print *, 'src x'
    call create_proj_copies(h, d, b, f, h_p, d_p, b_p, f_p)
    call project_points(sundir, a, compute_normal_3d(c, a, e), h_p, d_p, b_p, f_p)
    call rearange_projections(g, c, a, e, h_p, d_p, b_p, f_p)
    call correct_coeffs( &
      g   , c   , a   , e,    d, b, & ! fixed
      h_p , d_p , b_p , f_p,     & ! projected
      optical_props(1), [5      , 8       , 2], [2, 1, 3]      , coeffs       & ! slice of relevant coefficients , and coefficient array
      )

    if (lDEBUG_geometric_coeff_correction) print *, 'src y'
    call create_proj_copies(f, b, a, e, f_p, b_p, a_p, e_p)
    call project_points(sundir, d, compute_normal_3d(c, d, h), f_p, b_p, a_p, e_p)
    call rearange_projections(h, d, c, g, f_p, b_p, a_p, e_p)
    call correct_coeffs( &
      h   , d   , c   , g,    b, a, & ! fixed
      f_p , b_p , a_p , e_p,     & ! projected
      optical_props(1), [9      , 6       , 3], [1, 2, 3]      , coeffs       & ! slice of relevant coefficients , and coefficient array
      )
  end associate

  contains

    subroutine correct_coeffs( &
        f1, f2, f3, f4, f5, f6, &
        v1, v2, v3, v4, &
        extinction_coeff, slice, other_slice, coeffs &
        )
      real(ireals), intent(in) :: extinction_coeff
      real(ireals), intent(in), dimension(3) ::  &
        f1, f2, f3, f4, f5, f6, v1, v2, v3, v4
      integer(iintegers), intent(in) :: slice(3), other_slice(3)
      real(irealLUT), intent(inout) :: coeffs(9)
      real(ireals) :: &
        area1, area2, area3, area_total_src, areas(3), &
        s1, s21, s22, s23, s24, s31, s32, s33, s34, sin_theta
      real(ireals), dimension(3) :: &
        p1l, p1b, p1t, p1r, p2l, p2t, p2b, p2r, p3r, p3t, p3b, p3l, p4r, p4b, p4t, p4l, &
        v1_p, n
      integer(iintegers) :: coord_is(3)

      call proj_vars_to_edges( &
        f1, f2, f3, f4, &
        v1, v2, v3, v4, &
        p1l, p1b, p1t, p1r, &
        p2l, p2t, p2b, p2r, &
        p3r, p3t, p3b, p3l, &
        p4r, p4b, p4t, p4l  &
        )

      area_total_src = quadrangle_area_by_vertices(f1, f2, f3, f4)

      area2 = compute_quadrangle_areas( &
        v1, p1b, f2, f1, &
        v2, p2t, f1, f2, &
        v3, p3t, f4, f3, &
        v4, p4b, f4, f3  &
        )
      area3 = compute_pentagon_areas( &
        v1, p1l, f4, p1t, f1, &
        v2, f2, p2b, f3, p2l, &
        v3, p3r, f2, p3b, f3, &
        v4, f4, p4t, f1, p4r  &
        )
      area1 = area_total_src - area2 - area3

      n = compute_normal_3d(f5, f2, f3)
      v1_p = v1 + hit_plane(v1, sundir, f5, compute_normal_3d(f5, f2, f3)) * sundir


      s1 =  norm2(f1 - (f1 + hit_plane(f1, sundir, f5, compute_normal_3d(f1, f2, f3)) * sundir))

      area1 = area1 * exp( - extinction_coeff * s1)
      area1 = 0._irealLUT
      print *, 'area', area1

      area2 = compute_quadrangle_areas( &
        v1, p1b, f2, f1, &
        v2, p2t, f1, f2, &
        v3, p3t, f4, f3, &
        v4, p4b, f4, f3  &
        )

      print *, 'v1', v1
      print *, 'v2', v2
      print *, 'v3', v3
      print *, 'v4', v4

      coord_is = other_slice

      !sin_theta = sin(abs(atan(sundir(2) / sqrt(sundir(1)**2 + sundir(3)**2)))) !x, z
      !sin_theta = sin(abs(atan(sundir(1) / sqrt(sundir(3)**2 + sundir(2)**2))))  !y
      sin_theta = max(sin(abs(atan(sundir(coord_is(1)) / sqrt(sundir(coord_is(2))**2 + sundir(coord_is(3))**2)))), tiny(sin_theta))

      s21 = norm2(v1 - (v1 + hit_plane(v1, sundir, f1, compute_normal_3d(f1, f2, f5)) * sundir))
      s22 = norm2(v2 - (v2 + hit_plane(v2, sundir, f2, compute_normal_3d(f1, f2, f5)) * sundir))
      s23 = norm2(v3 - (v3 + hit_plane(v3, sundir, f3, compute_normal_3d(f4, f3, f6)) * sundir))
      s24 = norm2(v4 - (v4 + hit_plane(v4, sundir, f4, compute_normal_3d(f4, f3, f6)) * sundir))

      area2 = max( &
        quadrangle_area_by_vertices(v1, p1r, f2, p1b) * &
        (one - exp( - extinction_coeff * s21)) / max(tiny(area2), (extinction_coeff * s21)) + &
        num(f1(coord_is(3)) - v1(coord_is(3)), f1(coord_is(1)) - v1(coord_is(1)), extinction_coeff, sin_theta) &
        , &
        quadrangle_area_by_vertices(v2, p2r, f1, p2t) * &
        (one - exp( - extinction_coeff * s22)) / max(tiny(area2), (extinction_coeff * s22)) + &
        num(v2(coord_is(3)), f2(coord_is(1)) - v2(coord_is(1)), extinction_coeff, sin_theta) &
        , &
        quadrangle_area_by_vertices(v3, p3l, f4, p3t) * &
        (one - exp( - extinction_coeff * s23)) / max(tiny(area2), (extinction_coeff * s23)) + &
        num(v3(coord_is(3)), v3(coord_is(1)), extinction_coeff, sin_theta) &
        , &
        quadrangle_area_by_vertices(v4, p4l, f3, p4b) * &
        (one - exp( - extinction_coeff * s24)) / max(tiny(area2), (extinction_coeff * s24)) + &
        num(f4(coord_is(3)) - v4(coord_is(3)), v4(coord_is(1)), extinction_coeff, sin_theta) &
        )

      sin_theta = max(sin(abs(atan(sundir(coord_is(3)) / sqrt(sundir(coord_is(1))**2 + sundir(coord_is(2))**2)))), tiny(sin_theta))
      ! 90 - dotproduct(sundir, facenormal), (dot product = cos)

      s31 = norm2(v1 - (v1 + hit_plane(v1, sundir, f1, compute_normal_3d(f3, f2, f5)) * sundir))
      s32 = norm2(v2 - (v2 + hit_plane(v2, sundir, f2, compute_normal_3d(f3, f2, f5)) * sundir))
      s33 = norm2(v3 - (v3 + hit_plane(v3, sundir, f3, compute_normal_3d(f3, f2, f5)) * sundir))
      s34 = norm2(v4 - (v4 + hit_plane(v4, sundir, f4, compute_normal_3d(f3, f2, f5)) * sundir))

      area3 = max( &
        quadrangle_area_by_vertices(v1, p1l, f4, p1t) * &
        (one - exp( - extinction_coeff * s31)) / max(tiny(area3), (extinction_coeff * s31)) + &
        num(f1(coord_is(1)) - v1(coord_is(1)), f1(coord_is(3)) - v1(coord_is(3)), extinction_coeff, sin_theta) & ! src y
        , &
        quadrangle_area_by_vertices(v2, p2l, f3, p2b) * &
        (one - exp(- extinction_coeff * s32)) / max(tiny(area3), (extinction_coeff * s32)) + &
        num(f2(coord_is(1)) - v2(coord_is(1)), v2(coord_is(3)), extinction_coeff, sin_theta) & ! z
        , &
        quadrangle_area_by_vertices(v3, p3r, f2, p3b) * &
        (one - exp(- extinction_coeff * s33)) / max(tiny(area3), (extinction_coeff * s33)) + &
        num(v3(coord_is(1)), v3(coord_is(3)), extinction_coeff, sin_theta) &
        , &
        quadrangle_area_by_vertices(v4, p4r, f1, p4t) * &
        (one - exp( - extinction_coeff * s34)) / max(tiny(area3), (extinction_coeff * s34)) + &
        num(v4(coord_is(1)), f4(coord_is(3)) - v4(coord_is(3)), extinction_coeff, sin_theta) &
        )
      area3 = 0._irealLUT

      areas = max([area1, area2, area3], zero)

      coeffs(slice) = real(areas / area_total_src, irealLUT)

    end subroutine

    function num(l0, h0, extinction_coeff, sin_theta)
      real(ireals), intent(in) :: l0, h0, extinction_coeff, sin_theta
      integer(iintegers), parameter :: n = 20
      integer(iintegers) :: i
      real(ireals) :: dl, num, l, h

      dl = max(tiny(dl), l0 / n)
      num = zero
      do i=1,n
        l = (i - 0.5_ireals) * dl
        h = l * max(h0, zero) / max(tiny(l0), l0)
        num = num + dl * f(h, extinction_coeff, sin_theta)
      enddo
    end function

    real(ireals) function f(z, extinction_coeff, sin_theta)
      real(ireals), intent(in) :: z, extinction_coeff, sin_theta

      f = z * (one - exp(-extinction_coeff * z / sin_theta)) / max(tiny(f), (extinction_coeff * z / sin_theta))
    end function

    function compute_quadrangle_areas( &
        a1, a2, a3, a4, &
        b1, b2, b3, b4, &
        c1, c2, c3, c4, &
        d1, d2, d3, d4  &
        )
      real(ireals), intent(in), dimension(3) :: &
        a1, a2, a3, a4, &
        b1, b2, b3, b4, &
        c1, c2, c3, c4, &
        d1, d2, d3, d4
      real(ireals) :: compute_quadrangle_areas

      compute_quadrangle_areas = max( &
        quadrangle_area_by_vertices(a1, a2, a3, a4), &
        quadrangle_area_by_vertices(b1, b2, b3, b4), &
        quadrangle_area_by_vertices(c1, c2, c3, c4), &
        quadrangle_area_by_vertices(d1, d2, d3, d4)  &
        )
    end function

    function compute_pentagon_areas( &
        a1, a2, a3, a4, a5, &
        b1, b2, b3, b4, b5, &
        c1, c2, c3, c4, c5, &
        d1, d2, d3, d4, d5  &
        )
      real(ireals), intent(in), dimension(3) :: &
        a1, a2, a3, a4, a5, &
        b1, b2, b3, b4, b5, &
        c1, c2, c3, c4, c5, &
        d1, d2, d3, d4, d5
      real(ireals) :: compute_pentagon_areas

      compute_pentagon_areas = max( &
        pentagon_area_by_vertices(a1, a2, a3, a4, a5), &
        pentagon_area_by_vertices(b1, b2, b3, b4, b5), &
        pentagon_area_by_vertices(c1, c2, c3, c4, c5), &
        pentagon_area_by_vertices(d1, d2, d3, d4, d5)  &
        )
    end function
  end subroutine dir2dir3_geometric_coeff_corr

  subroutine create_proj_copies(p1, p2, p3, p4, c1, c2, c3, c4)
    real(ireals), intent(in), dimension(3) :: p1, p2, p3, p4
    real(ireals), intent(out), dimension(3) :: c1, c2, c3, c4

    c1 = p1
    c2 = p2
    c3 = p3
    c4 = p4
  end subroutine

  subroutine project_points(sundir, origin, normal, v1, v2, v3, v4)
    real(ireals), intent(in) :: sundir(3), normal(3), origin(3)
    real(ireals), intent(inout) :: v1(3), v2(3), v3(3), v4(3)
    real(ireals), parameter :: eps = 1._ireals / sqrt(epsilon(eps)) ! try to delete this

    v1 = v1 + min(hit_plane(v1, sundir, origin, normal), eps) * sundir
    v2 = v2 + min(hit_plane(v2, sundir, origin, normal), eps) * sundir
    v3 = v3 + min(hit_plane(v3, sundir, origin, normal), eps) * sundir
    v4 = v4 + min(hit_plane(v4, sundir, origin, normal), eps) * sundir

    v1 = v1 + hit_plane(v1, normal, origin, normal) * normal
    v2 = v2 + hit_plane(v2, normal, origin, normal) * normal
    v3 = v3 + hit_plane(v3, normal, origin, normal) * normal
    v4 = v4 + hit_plane(v4, normal, origin, normal) * normal
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

  subroutine proj_vars_to_edges( &
      f1, f2, f3, f4, &
      v1, v2, v3, v4, &
      p1l, p1b, p1t, p1r, &
      p2l, p2t, p2b, p2r, &
      p3r, p3t, p3b, p3l, &
      p4r, p4b, p4t, p4l  &
      )
    real(ireals), intent(in), dimension(3) :: f1, f2, f3, f4, v1, v2, v3, v4
    real(ireals), intent(out), dimension(3) :: p1l, p1b, p1t, p1r, p2l, p2t, p2b, p2r, p3r, p3t, p3b, p3l, p4r, p4b, p4t, p4l

    call proj_var_to_edges(f1, f2, f3, f4, v1, p1l, p1b, p1t, p1r)
    call proj_var_to_edges(f2, f1, f4, f3, v2, p2l, p2t, p2b, p2r)
    call proj_var_to_edges(f4, f3, f2, f1, v3, p3r, p3b, p3t, p3l)
    call proj_var_to_edges(f4, f3, f2, f1, v4, p4r, p4b, p4t, p4l)
  end subroutine
  subroutine rearange_projections(f1, f2, f3, f4, v1, v2, v3, v4)
    real(ireals), intent(in) :: f1(3), f2(3), f3(3), f4(3)
    real(ireals), intent(inout) :: v1(3), v2(3), v3(3), v4(3)

    call rearange_projection(f1-v1, f3, f4-f3, f2, f3-f2, v1)
    call rearange_projection(f2-v2, f3, f4-f3, f4, f1-f4, v2)
    call rearange_projection(f3-v3, f2, f1-f2, f4, f1-f4, v3)
    call rearange_projection(f4-v4, f2, f1-f2, f2, f3-f2, v4)
  end subroutine

  subroutine rearange_projection(direction1, origin2, direction2, origin3, direction3, origin1)
    real(ireals), intent(in), dimension(3) :: direction1, origin2, direction2, origin3, direction3
    real(ireals), intent(inout), dimension(3) :: origin1
    real(ireals) :: coeff21, coeff22, coeff31, coeff32
    integer(mpiint) :: ierr

    call line_intersection_3d(origin1, direction1, origin2, direction2, coeff21, coeff22, ierr)
    call line_intersection_3d(origin1, direction1, origin3, direction3, coeff31, coeff32, ierr)
    call rearange_point(origin1, direction1, min(max(coeff21, coeff31, zero), one), origin1)
  end subroutine
end module
