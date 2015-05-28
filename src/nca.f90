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
! Copyright (C) 2010-2015  Carolin Klinger, <carolin.klinger@physik.lmu.de>
!----------------------------------------------------------------------------
! Neighbouring Column Approximation
! RTE Solver for the thermal spectral range, calculation of heating rates
! Klinger and Mayer, 2015, **
! carolin.klinger@physik.lmu.de
!----------------------------------------------------------------------------

module m_ts_nca
  use m_data_parameters, only : ireals,iintegers
  implicit none

  private
  public :: ts_nca

contains
  subroutine ts_nca (dx, dy, dz, B, kabs_3d, Edn, Eup, hr)

    real(ireals), intent(in) :: dx, dy
    real(ireals), intent(in), dimension(:,:,:) :: dz, kabs_3d   ! dimensions including ghost values (Nlevel,Nx,Ny)
    real(ireals), intent(in), dimension(:,:,:) :: B,Edn,Eup     ! dimensions including ghost values (Nlevel,Nx,Ny)
    real(ireals), intent(out),dimension(:,:,:) :: hr            ! dimensions including ghost values (Nlevel,Nx,Ny)

    ! ############## Definition of variables ##################

    integer(iintegers)        :: ilyr = 0         ! counter z levels
    integer(iintegers)        :: ixx = 0          ! counter x grid boxes
    integer(iintegers)        :: iyy = 0          ! counter y grid boxes
    integer(iintegers)        :: iface = 0        ! counter side faces/directions of grid box

    real(ireals)           :: Absup = 0        ! Upwelling, Absorption, lower/upper face
    real(ireals)           :: Absdn = 0        ! Downwelling Absorption, lower/upper face
    real(ireals)           :: Absup1 = 0       ! Upwelling, Absorption, face 1
    real(ireals)           :: Absup2 = 0       ! Upwelling, Absorption, face 2
    real(ireals)           :: Absup3 = 0       ! Upwelling, Absorption, face 3
    real(ireals)           :: Absup4 = 0       ! Upwelling, Absorption, face 4
    real(ireals)           :: Absdn1 = 0       ! Downwelling Absorption, face 1  
    real(ireals)           :: Absdn2 = 0       ! Downwelling Absorption, face 2
    real(ireals)           :: Absdn3 = 0       ! Downwelling Absorption, face 3
    real(ireals)           :: Absdn4 = 0       ! Downwelling Absorption, face 4

    real(ireals)           :: Emup = 0         ! Upwelling, Emission, lower/upper face
    real(ireals)           :: Emdn = 0         ! Downwelling Emission, lower/upper face
    real(ireals)           :: Emup1 = 0        ! Upwelling, Emission, face 1
    real(ireals)           :: Emup2 = 0        ! Upwelling, Emission, face 2
    real(ireals)           :: Emup3 = 0        ! Upwelling, Emission, face 3
    real(ireals)           :: Emup4 = 0        ! Upwelling, Emission, face 4
    real(ireals)           :: Emdn1 = 0        ! Downwelling Emission, face 1  
    real(ireals)           :: Emdn2 = 0        ! Downwelling Emission, face 2  
    real(ireals)           :: Emdn3 = 0        ! Downwelling Emission, face 3  
    real(ireals)           :: Emdn4 = 0        ! Downwelling Emission, face 4  
   
    real(ireals)           :: HR_up = 0        ! Downwelling Emission, face 1  
    real(ireals)           :: HR_dn = 0        ! Downwelling Emission, face 2  
    real(ireals)           :: HR_up_s = 0      ! Downwelling Emission, face 3  
    real(ireals)           :: HR_dn_s = 0      ! Downwelling Emission, face 4  

    real(ireals)           :: mu = 0           ! zenith angle

    real(ireals)           :: areaweight = 0   ! area weight for non-cubic grid boxes

    real(ireals)           :: ax = 0           ! integration boarder for side contributions
    real(ireals)           :: bx = 0           ! integration boarder for side contributions 
    real(ireals)           :: cx = 0           ! integration boarder for side contributions

    real(ireals)           :: az = 0           ! integration boarder for side contributions
    real(ireals)           :: bz = 0           ! integration boarder for side contributions 
    real(ireals)           :: cz = 0           ! integration boarder for side contributions 

    real(ireals)           :: pi=3.141592653589793
 
    real(ireals)           :: factor = 0
    real(ireals)           :: l = 0.0
    real(ireals)           :: Trans = 0.0

    real(ireals)           :: L_up_3d(lbound(Eup,1):ubound(Eup,1),lbound(Eup,2):ubound(Eup,2),lbound(Eup,3):ubound(Eup,3))
    real(ireals)           :: L_dn_3d(lbound(Edn,1):ubound(Edn,1),lbound(Edn,2):ubound(Edn,2),lbound(Edn,3):ubound(Edn,3))
    real(ireals)           :: B1 = 0
    real(ireals)           :: B2 = 0
    real(ireals)           :: B1_1 = 0
    real(ireals)           :: B2_1 = 0

    integer(iintegers) :: is,ie,js,je,nlyr, Nx,Ny

    nlyr = ubound(hr,1)-1
    is = lbound(hr,2) +1
    ie = ubound(hr,2) -1
    js = lbound(hr,3) +1
    je = ubound(hr,3) -1

    Nx = ie-1
    Ny = je-1

    hr(1:nlyr,is:ie,js:je) = Edn(1:nlyr,is:ie,js:je) + Eup(2:nlyr+1,is:ie,js:je) - Eup(1:nlyr,is:ie,js:je) - Edn(2:nlyr+1,is:ie,js:je)
!    stop 'NCA is not freely available.... please consult Carolin Klinger for an implementation.'

    return
  end subroutine

end module m_ts_nca
