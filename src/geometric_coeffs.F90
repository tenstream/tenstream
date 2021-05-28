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
  triangle_area_by_vertices, compute_normal_3d, cross_3d, cstr, CHKERR, toStr, stripSpaces
use m_intersection, only: hit_plane, line_intersection_3d

implicit none
private
public :: dir2dir3_geometric_coeffs

logical, parameter :: ldebug= .True.
contains

  subroutine dir2dir3_geometric_coeffs(verts, sundir, bg, coeffs)
    real(ireals), intent(inout) :: coeffs(:)
    real(ireals), intent(in) :: verts(:), sundir(:), bg(:)
    real(ireals), dimension(3) :: d_p, b_p, a_p, h_p, f_p, e_p, g_p, a_p_r, b_p_r, d_p_r, e_p_r, f_p_r, g_p_r, h_p_r
    real(ireals) :: sun_up_down, extinction_coeff
    real(ireals), parameter :: small=sqrt(epsilon(extinction_coeff))

    extinction_coeff = bg(1) + bg(2)
    if (extinction_coeff .lt. small) &
      & call CHKERR(1_mpiint, 'Extinction coeff too small. bg1='// &
        & toStr(bg(1))//'; bg2='//toStr(bg(2)))

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
    call rearange_projections(d, c, a, b, h_p, g_p, e_p, f_p, h_p_r, g_p_r, e_p_r, f_p_r)
    call gomtrc_coeffs( &
      d   , c   , a   , b,    g, & ! fixed
      h_p , g_p , e_p , f_p,     & ! projections
      h_p_r , g_p_r , e_p_r , f_p_r,     & ! projections rearanged
      sun_up_down, &
      sundir, &
      extinction_coeff, &
      [integer(iintegers) :: 1, 7, 4], &
      [integer(iintegers) :: 2, 3, 1], &
      coeffs & ! slice of relevant coefficients , and coefficient array
      )


    if (ldebug) then
      print *, '_________________________________________________________________'
      print *, cstr('src x', 'blue')
    endif
    print *, 'normal comp'
    print *, 'c', c
    print *, 'a', a
    print *, 'e', e
    print *, compute_normal_3d(c,a,e)
    print *, 'done'
    call create_proj_copies(h, d, b, f, h_p, d_p, b_p, f_p)
    call project_points(sundir, a, compute_normal_3d(c, a, e), h_p, d_p, b_p, f_p)
    call rearange_projections(g, c, a, e, h_p, d_p, b_p, f_p, h_p_r, d_p_r, b_p_r, f_p_r)
    call gomtrc_coeffs( &
      g   , c   , a   , e,    d, & ! fixed
      h_p , d_p , b_p , f_p,     & ! projections
      h_p_r, d_p_r, b_p_r, f_p_r,     & ! projections rearanged
      sun_up_down, &
      sundir, &
      extinction_coeff, &
      [integer(iintegers) :: 5, 8, 2], &
      [integer(iintegers) :: 2, 1, 3], &
      coeffs & ! slice of relevant coefficients , and coefficient array
      )

    if (ldebug) then
      print *, '_________________________________________________________________'
      print *, cstr('src y', 'blue')
    endif
    call create_proj_copies(f, b, a, e, f_p, b_p, a_p, e_p)
    call project_points(sundir, d, compute_normal_3d(c, d, h), f_p, b_p, a_p, e_p)
    call rearange_projections(h, d, c, g, f_p, b_p, a_p, e_p, f_p_r, b_p_r, a_p_r, e_p_r)
    call gomtrc_coeffs( &
      h   , d   , c   , g,    b, & ! fixed
      f_p , b_p , a_p , e_p,     & ! projections
      f_p_r , b_p_r , a_p_r , e_p_r,     & ! projections rearanged
      sun_up_down, &
      sundir, &
      extinction_coeff, &
      [integer(iintegers) :: 9, 6, 3], &
      [integer(iintegers) :: 1, 2, 3], &
      coeffs & ! slice of relevant coefficients , and coefficient array
      )
  end associate

  if (ldebug) print *, '_________________________________________________________________'

  end subroutine

  subroutine gomtrc_coeffs( &
      f1, f2, f3, f4, f5, &
      v1, v2, v3, v4, &
      v1r,  v2r, v3r, v4r, &
      sun_up_down, &
      sundir, &
      extinction_coeff, &
      slice, &
      other_slice, &
      coeffs &
      )
    real(ireals), intent(in) :: sun_up_down
    real(ireals), intent(in) :: extinction_coeff
    real(ireals), intent(in), dimension(3) ::  f1, f2, f3, f4, f5, v1, v2, v3, v4, v1r, v2r, v3r, v4r, sundir
    integer(iintegers), intent(in) :: slice(3), other_slice(3)
    real(ireals), intent(inout) :: coeffs(9)
    real(ireals) :: area1, area2, area3, area_total_src, areas(3), sin_theta, cos_src_trgt_t, cos_src_trgt_b, s, &
      a2q, a2t, a3q, a3t, c, t, s_max, s_min, h_max, h_min
    real(ireals), dimension(3) :: pl, pr, pt, pb, normal2, normal3, prp, vrp
    real(ireals), parameter :: small = sqrt(epsilon(one))
    integer(mpiint) :: ierr

    area_total_src = quadrangle_area_by_vertices(f1, f2, f3, f4)
    normal2 = compute_normal_3d(f1, f2, f5)
    normal3 = compute_normal_3d(f3, f2, f5)
    sin_theta = max(sin(abs(atan(sundir(other_slice(1)) / &
      sqrt(sundir(other_slice(2))**2 + sundir(other_slice(3))**2)))), tiny(sin_theta))
    cos_src_trgt_b = cos(acos(dot_product(f1 - f2, f1 - f4) / (norm2(f1 - f2) * norm2(f1 - f4))) - Pi / 2)
    cos_src_trgt_t = cos(acos(dot_product(f1 - f2, f2 - f3) / (norm2(f1 - f2) * norm2(f2 - f3))) - Pi / 2)

    if (norm2(v1r-f1) .gt. small) then
      print *, cstr('v1', 'red')

      call line_intersection_3d(v1r, v1-v4, f3, f4-f3, c, t, ierr)
      call rearange_point(v1r, v1-v4, c, pl)
      call line_intersection_3d(pl, f4-f3, f3, f2-f3, c, t, ierr)
      call rearange_point(pl, f4-f3, max(c, zero), pl)

      call line_intersection_3d(v1r, f3-f2, f2, f1-f2, c, t, ierr)
      call rearange_point(v1r, f3-f2, c, pr)
      call line_intersection_3d(pr, f1-f2, f3, f2-f3, c, t, ierr)
      call rearange_point(pr, f1-f2, max(c, zero), pr)

      call line_intersection_3d(v1r, f4-f3, f1, f4-f1, c, t, ierr)
      call rearange_point(v1r, f4-f3, c, pt)

      call line_intersection_3d(v1r, f4-f3, f2, f3-f2, c, t, ierr)
      call rearange_point(v1r, f4-f3, c, pb)

      print *, 'fixed'
      print *, 'f1', f1
      print *, 'f2', f2
      print *, 'f3', f3
      print *, 'f4', f4
      print *, 'v1r', v1r

      print *, 'pr', pr
      print *, 'pt', pt
      print *, 'pl', pl

      a2q = quadrangle_area_by_vertices(v1r, pr, f2, pb)
      a2t = triangle_area_by_vertices(v1r, f1, pr)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v1r, pl, f4, pt)
      a3t = triangle_area_by_vertices(v1r, pt, f1)
      area3 = a3q + a3t
      area1 = area_total_src - area2 - area3
      if (ldebug) then
        print *, 'area1 pre extinction='//toStr(area1)
        print *, 'area2 pre extinction='//toStr(area2)//' with a2q='//toStr(a2q)//' and a2t='//toStr(a2t)
        print *, 'area3 pre extinction='//toStr(area3)//' with a3q='//toStr(a3q)//' and a3t='//toStr(a3t)
      endif
      vrp = v1r + hit_plane(v1r, sundir, f1, normal2) * sundir
      s = norm2(vrp - v1r)

      a2q = a2q * f_dst(s, extinction_coeff)
      a2t = num_dst(s, norm2(f1 - pr) * cos_src_trgt_t, norm2(pr - v1r), extinction_coeff)

      prp = pl + hit_plane(pl, sundir, f1, normal3) * sundir
      vrp = v1 + hit_plane(v1r, sundir, f1, normal3) * sundir
      s_max = norm2(prp - v1r)
      s_min = norm2(vrp - v1r)
      h_max = norm2(f4 - pl)
      h_min = norm2(pt - v1r)

      a3q = num_dst_q(s_max, s_min, norm2(f4-pt) * cos_src_trgt_b, h_max, h_min, extinction_coeff)
      a3t = num_dst(s_max, norm2(f1 - pt) * cos_src_trgt_t, norm2(pt - v1r), extinction_coeff)

    else if (norm2(v2r - f2) .gt. small) then
      print *, cstr('v2', 'red')

      call line_intersection_3d(v2r, v3-v2, f3, f4-f3, c, t, ierr)
      call rearange_point(v2r, v3-v2, c, pl)
      call line_intersection_3d(pl, f4-f3, f3, f2-f3, c, t, ierr)
      call rearange_point(pl, f4-f3, max(c, zero), pl)
      call line_intersection_3d(pl, f3-f4, f4, f4-f1, c, t, ierr)
      call rearange_point(pl, f3-f4, max(c, zero), pl)

      call line_intersection_3d(v2r, f1-f4, f2, f1-f2, c, t, ierr)
      call rearange_point(v2r, f1-f4, c, pr)
      call line_intersection_3d(pr, f1-f2, f3, f2-f3, c, t, ierr)
      call rearange_point(pr, f1-f2, max(c, zero), pr)
      call line_intersection_3d(pr, f2-f1, f4, f1-f4, c, t, ierr)
      call rearange_point(pr, f2-f1, max(c, zero), pr)

      call line_intersection_3d(v2r, v1-v2, f1, f4-f1, c, t, ierr)
      call rearange_point(v2r, v1-v2, c, pt)

      call line_intersection_3d(v2r, f4-f3, f2, f3-f2, c, t, ierr)
      call rearange_point(v2r, f4-f3, c, pb)

      a2q = quadrangle_area_by_vertices(v2r, pr, f1, pt)
      a2t = triangle_area_by_vertices(v2r, f2, pr)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v2r, pl, f3, pb)
      a3t = triangle_area_by_vertices(v2r, pb, f2)
      area3 = a3q + a3t
      area1 = area_total_src - area2 - area3
      if (ldebug) then
        print *, 'area1 pre extinction='//toStr(area1)
        print *, 'area2 pre extinction='//toStr(area2)//' with a2q='//toStr(a2q)//' and a2t='//toStr(a2t)
        print *, 'area3 pre extinction='//toStr(area3)//' with a3q='//toStr(a3q)//' and a3t='//toStr(a3t)
      endif
      prp = v2 + hit_plane(v2r, sundir, f2, normal2) * sundir
      s = norm2(prp - v2r)
      a2q = a2q * f_dst(s, extinction_coeff)
      a2t = num_dst(s, norm2(f2 - pr) * cos_src_trgt_b, norm2(pr - v2r), extinction_coeff)

      prp = pl + hit_plane(pl, sundir, f2, normal3) * sundir
      vrp = v2r + hit_plane(v2r, sundir, f2, normal3) * sundir
      s_max = norm2(prp - v2r)
      s_min = norm2(vrp - v2r)
      h_max = norm2(f3 - pl)
      h_min = norm2(pb - v2r)
      a3q = num_dst_q(s_max, s_min, norm2(f3 - pb) * cos_src_trgt_t, h_max, h_min, extinction_coeff)
      a3t = num_dst(s, norm2(f2 - pb) * cos_src_trgt_b, norm2(pb - v2r), extinction_coeff)

    else if (norm2(v3r - f3) .gt. small) then
      print *, cstr('v3', 'red')

      call line_intersection_3d(v3r, f4-f1, f3, f4-f3, c, t, ierr)
      call rearange_point(v3r, f4-f1, c, pl)
      call line_intersection_3d(pl, f4-f3, f3, f2-f3, c, t, ierr)
      call rearange_point(pl, f4-f3, max(c, zero), pl)
      call line_intersection_3d(pl, f3-f4, f4, f4-f1, c, t, ierr)
      call rearange_point(pl, f3-f4, max(c, zero), pl)

      call line_intersection_3d(v3r, v2-f3, f2, f1-f2, c, t, ierr)
      call rearange_point(v3r, v2-v3, c, pr)
      call line_intersection_3d(pr, f1-f2, f3, f2-f3, c, t, ierr)
      call rearange_point(pr, f1-f2, max(c, zero), pr)
      call line_intersection_3d(pr, f2-f1, f4, f1-f4, c, t, ierr)
      call rearange_point(pr, f2-f1, max(c, zero), pr)

      call line_intersection_3d(v3r, v4-v3, f1, f4-f1, c, t, ierr)
      call rearange_point(v3r, v4-v3, c, pt)

      call line_intersection_3d(v3r, f2-f1, f2, f3-f2, c, t, ierr)
      call rearange_point(v3r, f2-f1, c, pb)

      a2q = quadrangle_area_by_vertices(v3r, pl, f4, pt)
      a2t = triangle_area_by_vertices(v3r, f3, pl)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v3r, pr, f2, pb)
      a3t = triangle_area_by_vertices(v3r, pb, f4)
      area3 = a3q + a3t
      area1 = area_total_src - area2 - area3
      if (ldebug) then
        print *, 'area1 pre extinction='//toStr(area1)
        print *, 'area2 pre extinction='//toStr(area2)//' with a2q='//toStr(a2q)//' and a2t='//toStr(a2t)
        print *, 'area3 pre extinction='//toStr(area3)//' with a3q='//toStr(a3q)//' and a3t='//toStr(a3t)
      endif
      prp = v3 + hit_plane(v3r, sundir, f3, normal2) * sundir
      s = norm2(prp - v3r)
      a2q = a2q * f_dst(s, extinction_coeff)
      a2t = num_dst(s, norm2(f3 - pl) * cos_src_trgt_b, norm2(pl - v3r), extinction_coeff)

      prp = pr + hit_plane(pr, sundir, f3, normal3) * sundir
      vrp = v3r + hit_plane(v3r, sundir, f3, normal3) * sundir
      s_max = norm2(prp - v3r)
      s_min = norm2(vrp - v3r)
      h_max = norm2(f2 - pr)
      h_min = norm2(pb - v3r)

      a3q = num_dst_q(s_max, s_min, norm2(f2 - pb) * cos_src_trgt_t, h_max, h_min, extinction_coeff)
      a3t = num_dst(s, norm2(f3 - pb) * cos_src_trgt_b, norm2(pb - v3r), extinction_coeff)

    else if (norm2(v4r - f4) .gt. small) then
      print *, cstr('v4', 'red')

      call line_intersection_3d(v4r, f3-f2, f3, f4-f3, c, t, ierr)
      call rearange_point(v3r, f3-f2, c, pl)
      call line_intersection_3d(pl, f4-f3, f3, f2-f3, c, t, ierr)
      call rearange_point(pl, f4-f3, max(c, zero), pl)
      call line_intersection_3d(pl, f3-f4, f4, f4-f1, c, t, ierr)
      call rearange_point(pl, f3-f4, max(c, zero), pl)

      call line_intersection_3d(v4r, v1-f4, f2, f1-f2, c, t, ierr)
      call rearange_point(v4r, v1-v4, c, pr)
      call line_intersection_3d(pr, f1-f2, f3, f2-f3, c, t, ierr)
      call rearange_point(pr, f1-f2, max(c, zero), pr)
      call line_intersection_3d(pr, f2-f1, f4, f1-f4, c, t, ierr)
      call rearange_point(pr, f2-f1, max(c, zero), pr)

      call line_intersection_3d(v4r, f1-f2, f1, f4-f1, c, t, ierr)
      call rearange_point(v3r, f1-f2, c, pt)

      call line_intersection_3d(v4r, v3-v4, f2, f3-f2, c, t, ierr)
      call rearange_point(v4r, v3-v4, c, pb)

      a2q = quadrangle_area_by_vertices(v4r, pl, f3, pb)
      a2t = triangle_area_by_vertices(v4r, f4, pl)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v4r, pr, f1, pt)
      a3t = triangle_area_by_vertices(v4r, pt, f4)
      area3 = a3q + a3t
      area1 = area_total_src - area2 - area3
      if (ldebug) then
        print *, 'area1 pre extinction='//toStr(area1)
        print *, 'area2 pre extinction='//toStr(area2)//' with a2q='//toStr(a2q)//' and a2t='//toStr(a2t)
        print *, 'area3 pre extinction='//toStr(area3)//' with a3q='//toStr(a3q)//' and a3t='//toStr(a3t)
      endif
      prp = v4 + hit_plane(v4r, sundir, f4, normal2) * sundir
      s = norm2(prp - v4r)
      a2q = a2q * f_dst(s, extinction_coeff)
      a2t = num_dst(s, norm2(f4 - pr) * cos_src_trgt_t, norm2(pr - v4r), extinction_coeff)

      prp = pr + hit_plane(pr, sundir, f4, normal3) * sundir
      vrp = v4r + hit_plane(v3r, sundir, f4, normal3) * sundir
      s_max = norm2(prp - v4r)
      s_min = norm2(vrp - v4r)
      h_max = norm2(f1 - pr)
      h_min = norm2(pt - v4r)

      a3q = num_dst_q(s_max, s_min, norm2(f1 - pt) * cos_src_trgt_t, h_max, h_min, extinction_coeff)
      a3t = num_dst(s, norm2(f4 - pt) * cos_src_trgt_t, norm2(pt - v4r), extinction_coeff)
    endif

    area1 = area1 * &
      exp( - extinction_coeff * norm2(f1 - (f1 + hit_plane(f1, sundir, f5, compute_normal_3d(f1,f2,f3)) * sundir)))
    area2 = a2q + a2t
    area3 = a3q + a3t
    if (ldebug) then
      print *, 'area1 post extinction='//toStr(area1)
      print *, 'area2 post extinction='//toStr(area2)//' with a2q='//toStr(a2q)//' and a2t='//toStr(a2t)
      print *, 'area3 post extinction='//toStr(area3)//' with a3q='//toStr(a3q)//' and a3t='//toStr(a3t)
    endif

    area1 = area1 - sun_up_down * area3
    area3 = area3 + sun_up_down * area3
    if (ldebug) then
      print *, cstr('areas extinction included', 'red')
      print *, area1, area2, area3
    endif

    areas = max([area1, area2, area3], zero)
    coeffs(slice) = areas / area_total_src
    if (ldebug) print *, 'areas', coeffs(slice)

  end subroutine

  real(ireals) function num_dst(s0, l0, h0, extinction_coeff)
    real(ireals), intent(in) :: s0, l0, h0, extinction_coeff
    integer(iintegers), parameter :: n = 20
    integer(iintegers) :: i
    real(ireals) :: s, dh, ds, dl, h, j

    dl = l0 / n
    dh = h0 / n
    ds = s0 / n

    num_dst = zero
    do i=1, n
      j = i - 0.5_ireals
      s = j * ds
      h = j * dh
      num_dst = num_dst + dl * h * f_dst(s, extinction_coeff)
    enddo
  end function

  real(ireals) function num_dst_q(s_max, s_min, l, h_max, h_min, extinction_coeff)
    real(ireals), intent(in) :: s_max, s_min, l, h_max, h_min, extinction_coeff
    integer(iintegers), parameter :: n = 20
    integer(iintegers) :: i
    real(ireals) :: s, dh, ds, dl, h, j

    dl = l / n
    dh = (h_max - h_min) / n
    ds = (s_max - s_min) / n
    print *, 'dh', dh
    print *, 'ds', ds

    ! if dh .le. epsilon then use only one iteration
    num_dst_q = zero
    do i=1, n
      j = i - 0.5_ireals
      s = s_min + j * ds
      h = h_min + j * dh
      num_dst_q = num_dst_q + dl * h * f_dst(s, extinction_coeff)
    enddo
  end function

  real(ireals) function f_dst(s, extinction_coeff)
    real(ireals), intent(in) :: s, extinction_coeff
    real(ireals), parameter :: small = sqrt(tiny(f_dst))

    f_dst = (one - exp( - extinction_coeff * s)) / max(small, (extinction_coeff * s))
  end function

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
    real(ireals), parameter :: big = 1e5_ireals ! might be problematic

    if (ldebug) then
      print *, cstr('unprojected', 'yellow')
      print *, 'v1', v1
      print *, 'v2', v2
      print *, 'v3', v3
      print *, 'v4', v4
      print *, '_________________________________________________________________'
    endif

    if (abs(dot_product(sundir, normal)) < sqrt(epsilon(sundir))) then
      print *, cstr('sundir_proj being used', 'red')
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

  subroutine rearange_projections(f1, f2, f3, f4, v1, v2, v3, v4, v1r, v2r, v3r, v4r)
    real(ireals), intent(in), dimension(3) :: f1, f2, f3, f4, v1, v2, v3, v4
    real(ireals), intent(out), dimension(3) :: v1r, v2r, v3r, v4r

    call create_proj_copies(v1, v2, v3, v4, v1r, v2r, v3r, v4r)

    if (ldebug) then
      print *, cstr('rearangement coeffs', 'green')
      print *, 'v1'
    endif

    call rearange_projection(f1-v1, f3, f4-f3, f2, f3-f2, v1r)
    if (ldebug) print *, 'v2'
    call rearange_projection(f2-v2, f3, f4-f3, f4, f1-f4, v2r)
    if (ldebug) print *, 'v3'
    call rearange_projection(f3-v3, f2, f1-f2, f4, f1-f4, v3r)
    if (ldebug) print *, 'v4'
    call rearange_projection(f4-v4, f2, f1-f2, f2, f3-f2, v4r)

    if (ldebug) then
      print *, '_________________________________________________________________'
      print *, cstr('rearangements', 'yellow')
      print *, 'v1r', v1r
      print *, 'v2r', v2r
      print *, 'v3r', v3r
      print *, 'v4r', v4r
      print *, '_________________________________________________________________'
    endif
  end subroutine

  subroutine rearange_projection(direction1, origin2, direction2, origin3, direction3, origin1)
    real(ireals), intent(in), dimension(3) :: direction1, origin2, direction2, origin3, direction3
    real(ireals), intent(inout), dimension(3) :: origin1
    real(ireals) :: coeff21, coeff22, coeff31, coeff32
    integer(mpiint) :: ierr1, ierr2

    if (all(abs(direction1) .le. sqrt(epsilon(direction1)))) then
      origin1 = origin1
    else
      call line_intersection_3d(origin1, direction1, origin2, direction2, coeff21, coeff22, ierr1)
      call line_intersection_3d(origin1, direction1, origin3, direction3, coeff31, coeff32, ierr2)
      if (ldebug) then
        print *, coeff21, coeff31
        print *, ierr1, ierr2
        print *, '_________________________________________________________________'
      endif

      call rearange_point(origin1, direction1, min(max(coeff21, coeff31, zero), one), origin1)
    endif

  end subroutine
end module
