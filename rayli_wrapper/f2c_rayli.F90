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

#ifdef __HAVE_RAYLI__
  interface
    integer(c_int) function rfft_wedgeF90(&
        Nphotons, Nwedges, Nfaces, Nverts, &
        verts_of_face, wedges_of_face, vert_coords, &
        kabs, ksca, g, sundir, &
        flx_through_faces_edir, flx_through_faces_ediff) bind(c, name='rfft_wedge')
      use iso_c_binding
      integer(c_size_t), value :: Nphotons
      integer(c_size_t), value :: Nwedges
      integer(c_size_t), value :: Nfaces
      integer(c_size_t), value :: Nverts
      integer(c_size_t) :: verts_of_face(1:3,1:Nfaces)
      integer(c_size_t) :: wedges_of_face(1:2,1:Nfaces)
      real(c_double) :: vert_coords(1:3,1:Nverts)
      real(c_double) :: kabs(1:Nwedges), ksca(1:Nwedges), g(1:Nwedges)
      real(c_double) :: sundir(1:3)
      real(c_double) :: flx_through_faces_edir(1:Nfaces)
      real(c_double) :: flx_through_faces_ediff(1:Nfaces)
    end function
  end interface

#else

contains
  integer(c_int) function rfft_wedgeF90(&
        Nphotons, Nwedges, Nfaces, Nverts, &
        verts_of_face, wedges_of_face, vert_coords, &
        kabs, ksca, g, sundir, &
        flx_through_faces)
      use iso_c_binding
      integer(c_size_t), value :: Nphotons
      integer(c_size_t), value :: Nwedges
      integer(c_size_t), value :: Nfaces
      integer(c_size_t), value :: Nverts
      integer(c_size_t), intent(in) :: verts_of_face(:,:)
      integer(c_size_t), intent(in) :: wedges_of_face(:,:)
      real(c_double), intent(in) :: vert_coords(:,:)
      real(c_double), intent(in) :: kabs(:), ksca(:), g(:)
      real(c_double), intent(in) :: sundir(:)
      real(c_double), intent(out) :: flx_through_faces_edir(:)
      real(c_double), intent(out) :: flx_through_faces_ediff(:)

      rfft_wedgeF90 = 1
      call CHKERR(1_mpiint, "You tried calling The RayLi Monte Carlo solver "// &
        "but the Tenstream package was not compiled to use it.."// &
        " try to export RAYLI_DIR=<rayli-root>/build/package")

      if(.False.) then ! unused var warnings
        flx_through_faces_edir(1) = real(Nphotons+Nwedges+Nfaces+Nverts+verts_of_face(1,1)+wedges_of_face(1,1), c_double)
        flx_through_faces_ediff(1) = vert_coords(1,1) + kabs(1) + ksca(1) + g(1) + sundir(1)
      endif
    end function

#endif
end module

