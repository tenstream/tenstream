!----------------------------------------------------------------------------
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
! Copyright (C) 2010-2019  Carolin Klinger, <carolin.klinger@physik.lmu.de>
!----------------------------------------------------------------------------
! Neighbouring Column Approximation
! RTE Solver for the thermal spectral range, calculation of heating rates
! Klinger and Mayer, 2019, **
! carolin.klinger@physik.lmu.de
!----------------------------------------------------------------------------

module m_plexrt_nca
  use m_helper_functions, only: CHKERR
  use m_data_parameters, only : ireals, iintegers, mpiint
  implicit none

  private
  public :: plexrt_nca_init, plexrt_nca

  real(ireals), dimension(:,:), allocatable :: eps_tab_side, eps_tab_top
  real(ireals), dimension(:,:), allocatable :: corr_tab_side, corr_tab_top


contains
  subroutine plexrt_nca_init()
    integer(iintegers) :: i,j
    integer :: funit

    integer(iintegers), parameter :: ntau=16, n1=9, n2=36
    character(len=*), parameter :: eps_tab_side_fname = "Lookup/lookup_side_triangle.dat"
    character(len=*), parameter :: eps_tab_top_fname  = "Lookup/lookup_top_triangle.dat"
    character(len=*), parameter :: cor_tab_side_fname = "Lookup/lookup_correct_triangle_side.dat"
    character(len=*), parameter :: cor_tab_top_fname  = "Lookup/lookup_correct_triangle_top.dat"

    ! lookup table for emissivity
    call load_table(eps_tab_side_fname, ntau, ntau, eps_tab_side)
    call load_table(eps_tab_top_fname, ntau, ntau, eps_tab_top)

    ! lookup table for correction
    call load_table(cor_tab_side_fname, n1, n2, corr_tab_side)
    call load_table(cor_tab_top_fname, n1, n2, corr_tab_top)

  contains
    subroutine load_table(fname, n1, n2, arr)
      character(len=*), intent(in) :: fname
      integer(iintegers), intent(in) :: n1, n2
      real(ireals), allocatable, intent(inout) :: arr(:,:)
      logical :: lexists
      real(ireals) :: tmp

      if(.not.allocated(arr)) then
        allocate(arr(n1,n2))

        inquire(file=trim(fname), exist=lexists)
        if(.not.lexists) call CHKERR(1_mpiint, 'Could not find NCA LUT : '//trim(fname))

        open (newunit=funit, file=trim(fname), status="old", action="read")
        do i = 1, n1
          do j = 1, n2
            read(funit,*) tmp, tmp, arr(i,j)
          end do
        end do
        close(funit)
      endif
    end subroutine
  end subroutine


    ! NCA geometry for Wedges with vertices (ABC, DEF):
    !
    !
    !         F
    !        / \
    !    dx2/   \dx3
    !      /     \
    !     /  atop \
    !    /   dx1   \
    !   D __________E
    !   |           |
    !   |     |     |
    !   | a2  |  a3 |
    !   |     |     |
    !   |           |
    !   |     C     |
    !   |    / \    |
    !   |   /   \   |
    !   |  /     \  |
    !   | /       \ |
    !   |/  abot   \|
    !   A _________ B

    ! NCA flux/optprop info for Wedges with vertices (ABC, DEF), box on top has vertices (DEF, GHI):
    !      Base Wedge(1)      Neighbor Wedge
    !
    !         I                         F
    !        / \                       / \
    !       /   \                     /   \
    !      /     \                   /     \
    !     /       \                 /       \
    !    /         \               /         \
    !   G __________H             D __________E
    !
    !   |  kabs_top |
    !   |           |
    !   |     F     |                   F
    !   |    / \    |                  / \
    !   |   /   \   |                 /   \
    !   |  /     \  |                /     \
    !     /       \                 /       \
    !    /   Etop  \               /   Etop  \
    !   D __________E             D __________E
    !   |           |             |           |
    !   |   kabs    |             |   kabs    |
    !   |     |     |             |     |     |
    !   |     |     |             |     |     |
    !   |           |             |           |
    !   |     C     |             |     C     |
    !   |    / \    |             |    / \    |
    !   |   /   \   |             |   /   \   |
    !   |  /     \  |             |  /     \  |
    !   | /       \ |             | /       \ |
    !   |/   Ebot  \|             |/   Ebot  \|
    !   A _________ B             A _________ B


    subroutine plexrt_nca (dx1, dx2, dx3, dz, atop, abot, a1, a2, a3, v, &
        base_info, side_info, hr)
      real(ireals), intent(in) :: dx1, dx2, dx3          ! edge lengths of triangle: dx1, dx2, dx3
      real(ireals), intent(in) :: a1, a2, a3, atop, abot ! area of side faces and top and bot faces
      real(ireals), intent(in) :: dz, v                  ! height of grid box and volume

      ! info on this voxel dim(9) ( Edn_top, Eup_top, Btop, Edn_bot, Eup_bot, Bbot, kabs, kabs_top, kabs_bot)
      real(ireals), intent(in), dimension(:) :: base_info

      ! info for each of the side voxels dim(3 * 5) ( Edn_top, Eup_top, Edn_bot, Eup_bot, kabs)
      real(ireals), intent(in), dimension(:) :: side_info
      real(ireals), intent(out) :: hr        ! new 3D heating rate in the voxel

      ! ############## Definition of variables ##################

      if(.not. all( [ &
        allocated(eps_tab_side), &
        allocated(eps_tab_top), &
        allocated(corr_tab_side), &
        allocated(corr_tab_top) ] ) ) &
          call CHKERR(1_mpiint, 'called NCA but seems the NCA LUT`s havent been loaded ...'// &
          ' Try to call plexrt_nca_init() first')

      ! Dummy computations to get rid of unused compiler warnings
      hr = dx1 + dx2 + dx3 + a1 + a2 + a3 + atop + abot + dz + v + base_info(1) + side_info(1)
    end subroutine
end module m_plexrt_nca
