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
    real(ireals) :: area1, area2, area3, area_total_src, areas(3), sin_theta, cos_src_trgt, s1, s2, s3, s4, &
      a21q, a22q, a23q, a24q, a21t, a22t, a23t, a24t, a31q, a32q, a33q, a34q, a31t, a32t, a33t, a34t
    real(ireals), dimension(3) :: &
      p1l, p1b, p1t, p1r, p2l, p2t, p2b, p2r, p3r, p3t, p3b, p3l, p4r, p4b, p4t, p4l, normal, &
      p1rp, p2rp, p3rp, p4rp
    real(ireals), parameter :: small = sqrt(espilon(one))

    area_total_src = quadrangle_area_by_vertices(f1, f2, f3, f4)
    normal2 = compute_normal_3d(f1,f2,f5)
    normal3 = compute_normal_3d(f3,f2,f5)
    sin_theta = max(sin(abs(atan(sundir(other_slice(1)) / &
      sqrt(sundir(other_slice(2))**2 + sundir(other_slice(3))**2)))), tiny(sin_theta))
    cos_src_trgt = cos(acos(dot_product(f1 - f2, f1 - f4) / (norm2(f1 - f2) * norm2(f1 - f4))) - Pi / 2)

    if (norm2(p1r-f1) .gt. small) then
      call proj_var_to_edges(f1, f2, f3, f4, v1, v2, v3, v4, v1r, p1l, p1b, p1t, p1r)
      a2q = quadrangle_area_by_vertices(v1r, p1r, f2, p1b)
      a2t = triangle_area_by_vertices(v1r, f1, p1r)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v1r, p1l, f4, p1t)
      a3t = triangle_area_by_vertices(v1r, p1t, f1)
      area3 = a3q + a3t
      area1 = area_total_src - area2 - area3
      if (ldebug) then
        print *, 'area1 pre extinction='//toStr(area1)
        print *, 'area2 pre extinction='//toStr(area2)//' with a2q='//toStr(a2q)//' and a2t='//toStr(a2t)
        print *, 'area3 pre extinction='//toStr(area3)//' with a3q='//toStr(a3q)//' and a3t='//toStr(a3t)
      endif
      prp = v1 + hit_plane(v1r, sundir, f1, normal2) * sundir
      s = norm2(prp - v1r)

      a2q = a2q * f_dst(s, extinction_coeff)
      a2t = num_dst(s, norm2(f1 - p1r) * cos_src_trgt, norm2(p1r - v1r), extinction_coeff)

      prp = v1 + hit_plane(v1r, sundir, f1, normal3) * sundir
      s = norm2(prp - v1r)
      a3q = a3q * f_dst(s, extinction_coeff)
      a3t = num_dst(s1, norm2(f1 - p1t) * cos_src_trgt, norm2(p1t - v1r), extinction_coeff)

    else if (norm2(v2r - f2) .gt. small) then
      call proj_var_to_edges(f2, f1, f4, f3, v1, v2, v3, v4, v2r, p2l, p2t, p2b, p2r)
      a2q = quadrangle_area_by_vertices(v2r, p2r, f1, p2t)
      a2t = triangle_area_by_vertices(v2r, f2, p2r)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v2r, p2l, f3, p2b)
      a3t = triangle_area_by_vertices(v2r, p2b, f2)
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
      a2t = num_dst(s, norm2(f2 - p2r) * cos_src_trgt, norm2(p2r - v2r), extinction_coeff)
      prp = v2 + hit_plane(v2r, sundir, f1, normal3) * sundir
      s = norm2(prp - v2r)
      a3q = a3q * f_dst(s, extinction_coeff)
      a3t = num_dst(s, norm2(f2 - p2t) * cos_src_trgt, norm2(p2t - v2r), extinction_coeff)

    else if (norm2(v3r - f3) .gt. small) then
      call proj_var_to_edges(f4, f3, f2, f1, v1, v2, v3, v4, v3r, p3r, p3b, p3t, p3l)
      a2q = quadrangle_area_by_vertices(v3r, p3l, f4, p3t)
      a2t = triangle_area_by_vertices(v3r, f3, p3l)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v3r, p3r, f2, p3b)
      a3t = triangle_area_by_vertices(v3r, p3b, f4)
      area3 = a3q + a3t
      area1 = area_total_src - area2 - area3
      if (ldebug) then
        print *, 'area1 pre extinction='//toStr(area1)
        print *, 'area2 pre extinction='//toStr(area2)//' with a2q='//toStr(a2q)//' and a2t='//toStr(a2t)
        print *, 'area3 pre extinction='//toStr(area3)//' with a3q='//toStr(a3q)//' and a3t='//toStr(a3t)
      endif
      prp = v3 + hit_plane(v3r, sundir, f3, normal2) * sundir
      s = norm2(prp - v2r)
      a2q = a2q * f_dst(s, extinction_coeff)
      a2t = num_dst(s, norm2(f3 - p3r) * cos_src_trgt, norm2(p3r - v3r), extinction_coeff)
      prp = v3 + hit_plane(v3r, sundir, f1, normal3) * sundir
      s = norm2(prp - v3r)
      a3q = a3q * f_dst(s, extinction_coeff)
      a3t = num_dst(s, norm2(f3 - p3t) * cos_src_trgt, norm2(p3t - v3r), extinction_coeff)

    else if (norm2(v4r - f4) .gt. small) then
      call proj_var_to_edges(f4, f3, f2, f1, v1, v2, v3, v4, v4r, p4r, p4b, p4t, p4l)
      a2q = quadrangle_area_by_vertices(v4r, p4l, f3, p4b)
      a2t = triangle_area_by_vertices(v4r, f4, p4l)
      area2 = a2q + a2t
      a3q = quadrangle_area_by_vertices(v4r, p4r, f1, p4t)
      a3t = triangle_area_by_vertices(v4r, p4t, f4)
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
      a2t = num_dst(s, norm2(f4 - p4r) * cos_src_trgt, norm2(p4r - v4r), extinction_coeff)
      prp = v4 + hit_plane(v4r, sundir, f1, normal3) * sundir
      s = norm2(prp - v4r)
      a3q = a3q * f_dst(s, extinction_coeff)
      a3t = num_dst(s, norm2(f4 - p4t) * cos_src_trgt, norm2(p4t - v4r), extinction_coeff)
    endif

    area1 = area1 * &
      exp( - extinction_coeff * norm2(f1 - (f1 + hit_plane(f1, sundir, f5, compute_normal_3d(f1,f2,f3)) * sundir)))
    area2 = a2q + a2t
    area3 = a2q + a2t
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
    !END

    area1 = area1 - sun_up_down * area3
    area3 = area3 + sun_up_down * area3
    if (ldebug) then
      print *, cstr('areas extinction included', 'red')
      print *, area1, area2, area3
    endif
    areas = max([area1, area2, area3], zero)
    coeffs(slice) = areas / area_total_src

    ! end new

    a23q = quadrangle_area_by_vertices(v3r, p3l, f4, p3t)
    a24q = quadrangle_area_by_vertices(v4r, p4l, f3, p4b)

    a23t = triangle_area_by_vertices(v3r, f3, p3l)
    a24t = triangle_area_by_vertices(v4r, f4, p4l)


    a33q = quadrangle_area_by_vertices(v3r, p3r, f2, p3b)
    a34q = quadrangle_area_by_vertices(v4r, p4r, f1, p4t)

    a33t = triangle_area_by_vertices(v3r, p3b, f3)
    a34t = triangle_area_by_vertices(v4r, p4t, f4)

    area3 = max(a31q+a31t, a32q+a32t, a33q+a33t, a34q+a34t)
    area1 = area_total_src - area2 - area3

    s1 =  norm2(f1 - (f1 + hit_plane(f1, sundir, f5, compute_normal_3d(f1, f2, f3)) * sundir))

    area1 = area1 * exp( - extinction_coeff * s1)

    normal = compute_normal_3d(f1, f2, f5)

    p1rp = v1 + hit_plane(v1r, sundir, f1, normal) * sundir
    p2rp = v2 + hit_plane(v2r, sundir, f2, normal) * sundir
    p3rp = v3 + hit_plane(v3r, sundir, f3, normal) * sundir
    p4rp = v4 + hit_plane(v4r, sundir, f4, normal) * sundir

    s1 = norm2(p1rp - v1r)
    s2 = norm2(p2rp - v2r)
    s3 = norm2(p3rp - v3r)
    s4 = norm2(p4rp - v4r)

    sin_theta = max(sin(abs(atan(sundir(other_slice(1)) / &
      sqrt(sundir(other_slice(2))**2 + sundir(other_slice(3))**2)))), tiny(sin_theta))
    cos_src_trgt = cos(acos(dot_product(f1 - f2, f1 - f4) / (norm2(f1 - f2) * norm2(f1 - f4))) - Pi / 2)

    a21q = a21q * f_dst(s1, extinction_coeff)
    a22q = a22q * f_dst(s2, extinction_coeff)
    a23q = a23q * f_dst(s3, extinction_coeff)
    a24q = a24q * f_dst(s4, extinction_coeff)

    a21t = num_dst(s1, norm2(f1 - p1r) * cos_src_trgt, norm2(p1r - v1r), extinction_coeff)
    a22t = num_dst(s2, norm2(f2 - p2r) * cos_src_trgt, norm2(p2r - v2r), extinction_coeff)
    a23t = num_dst(s3, norm2(f3 - p3l) * cos_src_trgt, norm2(p3l - v3r), extinction_coeff)
    a24t = num_dst(s4, norm2(f4 - p4l) * cos_src_trgt, norm2(p4l - v4r), extinction_coeff)

    area2 = max(a21q + a21t, a22q + a22t, a23q + a23t, a24q + a24t)

    if (ldebug) then
      print *, cstr('area2 computation', 'yellow')
      print *, cstr('quadrangle areas, triangle areas', 'yellow')
      print *, stripSpaces('v1, s1='//toStr(s1)//', a21q='//toStr(a21q)//', a21t='//toStr(a21t))
      print *, stripSpaces('v2, s2='//toStr(s2)//', a22q='//toStr(a22q)//', a22t='//toStr(a21t))
      print *, stripSpaces('v3, s3='//toStr(s3)//', a23q='//toStr(a23q)//', a23t='//toStr(a23t))
      print *, stripSpaces('v4, s4='//toStr(s4)//', a24q='//toStr(a24q)//', a24t='//toStr(a24t))
      print *, '_________________________________________________________________'
    endif

    normal = compute_normal_3d(f3, f2, f5)

    p1rp = v1 + hit_plane(v1r, sundir, f1, normal) * sundir
    p2rp = v2 + hit_plane(v2r, sundir, f2, normal) * sundir
    p3rp = v3 + hit_plane(v3r, sundir, f3, normal) * sundir
    p4rp = v4 + hit_plane(v4r, sundir, f4, normal) * sundir

    s1 = norm2(p1rp - v1r)
    s2 = norm2(p2rp - v2r)
    s3 = norm2(p3rp - v3r)
    s4 = norm2(p4rp - v4r)

    a31q = a31q * f_dst(s1, extinction_coeff)
    a32q = a32q * f_dst(s2, extinction_coeff)
    a33q = a33q * f_dst(s3, extinction_coeff)
    a34q = a34q * f_dst(s4, extinction_coeff)

    a31t = num_dst(s1, norm2(f1 - p1t) * cos_src_trgt, norm2(p1t - v1r), extinction_coeff)
    a32t = num_dst(s2, norm2(f2 - p2b) * cos_src_trgt, norm2(p2b - v2r), extinction_coeff)
    a33t = num_dst(s3, norm2(f3 - p3b) * cos_src_trgt, norm2(p3b - v3r), extinction_coeff)
    a34t = num_dst(s4, norm2(f4 - p4t) * cos_src_trgt, norm2(p4t - v4r), extinction_coeff)

    area3 = max(a31q + a31t, a32q + a32t, a33q + a33t, a34q + a34t)

    if (ldebug) then
      print *, cstr('area3 computation', 'yellow')
      print *, cstr('quadrangle areas, triangle areas', 'yellow')
      print *, 'v1, s1='//stripSpaces(toStr(s1))//&
        &', a31q='//stripSpaces(toStr(a31q))//', f_dst='//stripSpaces(toStr(f_dst(s1, extinction_coeff)))//&
        &', a31t='//stripSpaces(toStr(a31t))
      print *, 'v2, s2='//stripSpaces(toStr(s2))//&
        &', a32q='//stripSpaces(toStr(a32q))//', f_dst='//stripSpaces(toStr(f_dst(s2, extinction_coeff)))//&
        &', a32t='//stripSpaces(toStr(a31t))
      print *, 'v3, s3='//stripSpaces(toStr(s3))//&
        &', a33q='//stripSpaces(toStr(a33q))//', f_dst='//stripSpaces(toStr(f_dst(s3, extinction_coeff)))//&
        &', a33t='//stripSpaces(toStr(a33t))
      print *, 'v4, s4='//stripSpaces(toStr(s4))//&
        &', a34q='//stripSpaces(toStr(a34q))//', f_dst='//stripSpaces(toStr(f_dst(s4, extinction_coeff)))//&
        &', a34t='//stripSpaces(toStr(a34t))
      print *, '_________________________________________________________________'
    endif

    if (ldebug) then
      print *, cstr('sun_up_down', 'green')
      print *, sun_up_down
      print *, 'area1', area1, ' -> ', area1 - sun_up_down * area3
      print *, 'area3', area3, ' -> ', area3 + sun_up_down * area3
      print *, 'sundir', sundir
    endif

    area1 = area1 - sun_up_down * area3
    area3 = area3 + sun_up_down * area3

    ! do not recompute areas, but just multiply partwise with extincten / replace triangle

    if (ldebug) then
      print *, cstr('areas extinction included', 'red')
      print *, area1, area2, area3
    endif

    areas = max([area1, area2, area3], zero)

    coeffs(slice) = areas / area_total_src

  end subroutine

  subroutine gomtrc_coeffs_comp( &
      f1, f2, f3, f4, &
      v1, v2, v3, v4, &
      vr, pr, pl, pb, pt &
      )
    real(ireals), dimension(3), intent(in) :: f1, f2, f3, f4, v1, v2, v3, v4, vr, pr, pl, pb, pt
    real(ireals), intent(out) :: areas(3)

    if (ldebug) then
      print *, cstr('fixed', 'yellow')
      print *, 'f1', f1
      print *, 'f2', f2
      print *, 'f3', f3
      print *, 'f4', f4
      print *, cstr('variable', 'green')
      print *, cstr('v1r', 'green'), v1r
      print *, 'p1l', p1l
      print *, 'p1r', p1r
      print *, 'p1b', p1b
      print *, 'p1t', p1t
    endif

    area1 = area1 * &
      exp( - extinction_coeff * norm2(f1 - (f1 + hit_plane(f1, sundir, f5, compute_normal_3d(f1,f2,f3)) * sundir)))

    normal = compute_normal_3d(f1,f2,f5)
    prp = vr + hit_plane(vr, sundir, f1, normal) * sundir
    s = norm2(prp - vr)
    sin_theta = max(sin(abs(atan(sundir(other_slice(1)) / &
      sqrt(sundir(other_slice(2))**2 + sundir(other_slice(3))**2)))), tiny(sin_theta))
    cos_src_trgt = cos(acos(dot_product(f1 - f2, f1 - f4) / (norm2(f1 - f2) * norm2(f1 - f4))) - Pi / 2)

    a2q = a2q * f_dst(s, extinction_coeff)

    a21t = num_dst(s1, norm2(f1 - p1r) * cos_src_trgt, norm2(p1r - v1r), extinction_coeff)
    a22t = num_dst(s2, norm2(f2 - p2r) * cos_src_trgt, norm2(p2r - v2r), extinction_coeff)
    a23t = num_dst(s3, norm2(f3 - p3l) * cos_src_trgt, norm2(p3l - v3r), extinction_coeff)
    a24t = num_dst(s4, norm2(f4 - p4l) * cos_src_trgt, norm2(p4l - v4r), extinction_coeff)


    normal = compute_normal_3d(f1,f2,f5)
    p1rp = v1 + hit_plane(v1r, sundir, f1, normal) * sundir
    s1 = norm2(p1rp - v1r)
    sin_theta = max(sin(abs(atan(sundir(other_slice(1)) / &
      sqrt(sundir(other_slice(2))**2 + sundir(other_slice(3))**2)))), tiny(sin_theta))
    cos_src_trgt = cos(acos(dot_product(f1 - f2, f1 - f4) / (norm2(f1 - f2) * norm2(f1 - f4))) - Pi / 2)
    a21q = a21q * f_dst(s1, extinction_coeff)
    a21t = num_dst(s1, norm2(f1 - p1r) * cos_src_trgt, norm2(p1r - v1r), extinction_coeff)
    area2 = a21q + a21t
    if (ldebug) print *, 'area2 pre extinction='//toStr(area)//' with a21q='//toStr(a21q)//' and a21t='//toStr(a21t)
    normal = compute_normal_3d(f3, f2, f5)
    p1rp = v1 + hit_plane(v1r, sundir, f1, normal) * sundir
    s1 = norm2(p1rp - v1r)
    a31q = a31q * f_dst(s1, extinction_coeff)
    a31t = num_dst(s1, norm2(f1 - p1t) * cos_src_trgt, norm2(p1t - v1r), extinction_coeff)
    area3 = a31q + a31t
    if (ldebug) print *, 'area3 pre extinction='//toStr(area)//' with a31q='//toStr(a31q)//' and a31t='//toStr(a31t)

    area1 = area1 - sun_up_down * area3
    area3 = area3 + sun_up_down * area3
    if (ldebug) then
      print *, cstr('areas extinction included', 'red')
      print *, area1, area2, area3
    endif
    areas = max([area1, area2, area3], zero)
    coeffs(slice) = areas / area_total_src
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
      f1, f2, f3, f4, & ! side's corners
      v1, v2, v3, v4, & ! var side's corners
      vr, &             ! point to be projected
      p1, p2, p3, p4 &  ! projections
      )
    real(ireals), intent(in), dimension(3) :: f1, f2, f3, f4, vr, v1, v2, v3, v4
    real(ireals), intent(out), dimension(3) :: p1, p2, p3, p4
    real(ireals) :: c, t
    integer(mpiint) :: ierr

    call line_intersection_3d(vr, v4-v1, f3, f4-f3, c, t, ierr)
    call rearange_point(vr, v4-v1, c, p1)
    call line_intersection_3d(vr, v2-v1, f2, f3-f2, c, t, ierr)
    call rearange_point(vr, v2-v1, c, p2)
    call line_intersection_3d(vr, v1-v2, f4, f1-f4, c, t, ierr)
    call rearange_point(vr, v1-v2, c, p3)
    call line_intersection_3d(vr, v4-v1, f2, f1-f2, c, t, ierr)
    call rearange_point(vr, v3-v2, c, p4)

    if (ldebug) then
      print *, cstr('fixed', 'yellow')
      print *, 'f1', f1
      print *, 'f2', f2
      print *, 'f3', f3
      print *, 'f4', f4
      print *, cstr('projections', 'green')
      print *, 'v1', v1
      print *, 'v2', v2
      print *, 'v3', v3
      print *, 'v4', v4
      print *, cstr('responsible', 'red')
      print *, cstr('v2r', 'green'), v2r
      print *, 'p1', p1
      print *, 'p2', p2
      print *, 'p3', p3
      print *, 'p4', p4
    endif
  end subroutine

  subroutine proj_vars_to_edges( &
      f1, f2, f3, f4, &
      v1, v2, v3, v4, &
      v1r, v2r, v3r, v4r, &
      p1l, p1b, p1t, p1r, &
      p2l, p2t, p2b, p2r, &
      p3r, p3t, p3b, p3l, &
      p4r, p4b, p4t, p4l  &
      )
    real(ireals), intent(in), dimension(3) :: f1, f2, f3, f4, v1, v2, v3, v4, v1r, v2r, v3r, v4r
    real(ireals), intent(out), dimension(3) :: p1l, p1b, p1t, p1r, p2l, p2t, p2b, p2r, p3r, p3t, p3b, p3l, p4r, p4b, p4t, p4l

    call proj_var_to_edges(f1, f2, f3, f4, v1, v2, v3, v4, v1r, p1l, p1b, p1t, p1r)
    call proj_var_to_edges(f2, f1, f4, f3, v1, v2, v3, v4, v2r, p2l, p2t, p2b, p2r)
    call proj_var_to_edges(f4, f3, f2, f1, v1, v2, v3, v4, v3r, p3r, p3b, p3t, p3l)
    call proj_var_to_edges(f4, f3, f2, f1, v1, v2, v3, v4, v4r, p4r, p4b, p4t, p4l)
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

    if (ldebug) then
      print *, 'v_i', origin1
      print *, 'dir1', direction1
      print *, 'o2', origin2
      print*, 'dir2', direction2
      print *, 'o3', origin3
      print*, 'dir3', direction3
    endif

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
