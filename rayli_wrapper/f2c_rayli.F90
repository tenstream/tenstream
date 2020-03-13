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

module m_f2c_rayli
  use iso_c_binding
  use iso_fortran_env, only: REAL64, INT32
  use m_data_parameters, only: mpiint
  use m_helper_functions, only: CHKERR

  implicit none

#ifdef HAVE_RAYLI
  interface
    integer(c_int) function rfft_wedgeF90(&
        Nphotons, Nwedges, Nfaces, Nverts, cyclic, &
        verts_of_face, faces_of_wedges, vert_coords, &
        kabs, ksca, g, albedo_on_faces, sundir, diffuse_point_origin, &
        flx_through_faces_edir, flx_through_faces_ediff, abso_in_cells) &
        bind(c, name='rfft_wedge')
      use iso_c_binding
      integer(c_size_t), value :: Nphotons
      integer(c_size_t), value :: Nwedges
      integer(c_size_t), value :: Nfaces
      integer(c_size_t), value :: Nverts
      integer(c_int)   , value :: cyclic
      integer(c_size_t) :: verts_of_face(1:4,1:Nfaces)
      integer(c_size_t) :: faces_of_wedges(1:5,1:Nfaces)
      real(c_double) :: vert_coords(1:3,1:Nverts)
      real(c_float ) :: kabs(1:Nwedges), ksca(1:Nwedges), g(1:Nwedges)
      real(c_float) :: albedo_on_faces(1:Nfaces)
      real(c_float) :: sundir(1:3)
      real(c_float) :: diffuse_point_origin(1:3)
      real(c_float) :: flx_through_faces_edir(1:Nfaces)
      real(c_float) :: flx_through_faces_ediff(1:Nfaces)
      real(c_float) :: abso_in_cells(1:Nwedges)
    end function
  end interface
  interface
    integer(c_int) function rpt_img_wedgeF90(&
        img_Nx, img_Ny, &
        Nphotons, Nwedges, Nfaces, Nverts, &
        verts_of_face, faces_of_wedges, vert_coords, &
        kabs, ksca, g, albedo_on_faces, &
        sundir, &
        cam_location, cam_viewing_dir, cam_up_vec, &
        fov_width, fov_height, &
        img) bind(c, name='rpt_img_wedge')
      use iso_c_binding
      integer(c_size_t), value :: img_Nx, img_Ny
      integer(c_size_t), value :: Nphotons
      integer(c_size_t), value :: Nwedges
      integer(c_size_t), value :: Nfaces
      integer(c_size_t), value :: Nverts
      integer(c_size_t) :: verts_of_face(1:4,1:Nfaces)
      integer(c_size_t) :: faces_of_wedges(1:5,1:Nfaces)
      real(c_double) :: vert_coords(1:3,1:Nverts)
      real(c_float ) :: kabs(1:Nwedges), ksca(1:Nwedges), g(1:Nwedges)
      real(c_float) :: albedo_on_faces(1:Nfaces)
      real(c_float) :: sundir(1:3)
      real(c_float), dimension(1:3) :: cam_location, cam_viewing_dir, cam_up_vec
      real(c_float), value :: fov_width, fov_height
      real(c_float) :: img(1:img_Nx, 1:img_Ny)
    end function
  end interface

#else

contains
  integer(c_int) function rfft_wedgeF90(&
        Nphotons, Nwedges, Nfaces, Nverts, cyclic, &
        verts_of_face, faces_of_wedges, vert_coords, &
        kabs, ksca, g, albedo_on_faces, sundir, diffuse_point_origin, &
        flx_through_faces_edir, flx_through_faces_ediff, abso_in_cells)
      use iso_c_binding
      integer(c_size_t), value :: Nphotons
      integer(c_size_t), value :: Nwedges
      integer(c_size_t), value :: Nfaces
      integer(c_size_t), value :: Nverts
      integer(c_int),    value :: cyclic
      integer(c_size_t), intent(in) :: verts_of_face(:,:)
      integer(c_size_t), intent(in) :: faces_of_wedges(:,:)
      real(c_double), intent(in) :: vert_coords(:,:)
      real(c_float ), intent(in) :: kabs(:), ksca(:), g(:)
      real(c_float ), intent(in) :: albedo_on_faces(1:Nfaces)
      real(c_float ), intent(in) :: sundir(:), diffuse_point_origin(:)
      real(c_float ), intent(out) :: flx_through_faces_edir(:)
      real(c_float ), intent(out) :: flx_through_faces_ediff(:)
      real(c_float ), intent(out) :: abso_in_cells(:)

      rfft_wedgeF90 = 1
      call CHKERR(1_mpiint, "You tried calling The RayLi Monte Carlo solver "// &
        "but the Tenstream package was not compiled to use it.."// &
        " try to export RAYLI_DIR=<rayli-root>/build/package")

      if(.False.) then ! unused var warnings
        flx_through_faces_edir(1) = real(Nphotons+Nwedges+Nfaces+Nverts+verts_of_face(1,1)+faces_of_wedges(1,1), c_float) + &
          real(cyclic, c_float)
        flx_through_faces_ediff(1) = real(vert_coords(1,1), c_float) + kabs(1) + ksca(1) + g(1) + sundir(1) + &
          albedo_on_faces(1) + diffuse_point_origin(1)
        abso_in_cells(1) = 0
      endif
    end function

    integer(c_int) function rpt_img_wedgeF90(&
        img_Nx, img_Ny, &
        Nphotons, Nwedges, Nfaces, Nverts, &
        verts_of_face, faces_of_wedges, vert_coords, &
        kabs, ksca, g, sundir, &
        albedo_on_faces, &
        cam_location, cam_viewing_dir, cam_up_vec, &
        fov_width, fov_height, img) bind(c, name='rpt_img_wedge')
      use iso_c_binding
      integer(c_size_t), value :: img_Nx, img_Ny
      integer(c_size_t), value :: Nphotons
      integer(c_size_t), value :: Nwedges
      integer(c_size_t), value :: Nfaces
      integer(c_size_t), value :: Nverts
      integer(c_size_t) :: verts_of_face(1:4,1:Nfaces)
      integer(c_size_t) :: faces_of_wedges(1:5,1:Nfaces)
      real(c_double) :: vert_coords(1:3,1:Nverts)
      real(c_float ) :: kabs(1:Nwedges), ksca(1:Nwedges), g(1:Nwedges)
      real(c_float) :: albedo_on_faces(1:Nfaces)
      real(c_float) :: sundir(1:3)
      real(c_float), dimension(1:3) :: cam_location, cam_viewing_dir, cam_up_vec
      real(c_float) :: fov_width, fov_height
      real(c_float) :: img(1:img_Nx, 1:img_Ny)

      rpt_img_wedgeF90 = 1
      call CHKERR(1_mpiint, "You tried calling The RayLi Monte Carlo solver "// &
        "but the Tenstream package was not compiled to use it.."// &
        " try to export RAYLI_DIR=<rayli-root>/build/package")

      if(.False.) then ! unused var warnings
        img(1,1) = real(Nphotons+Nwedges+Nfaces+Nverts+verts_of_face(1,1)+faces_of_wedges(1,1), c_float)
        img(2,1) = real(vert_coords(1,1), c_float) + kabs(1) + ksca(1) + g(1) + sundir(1) + albedo_on_faces(1)
        img(3,1) = cam_location(1) + cam_viewing_dir(1) + cam_up_vec(1) + fov_width + fov_height
        img(4,1) = real(img_Nx+img_Ny, c_float)
      endif
    end function
#endif
end module

