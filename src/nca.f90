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
! Klinger and Mayer, 2015, The Neighbouring Column Approximation - A fast approach for the calculation of 3D thermal heating rates, JQSRT
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

    integer(iintegers)     :: ilyr          ! counter z levels
    integer(iintegers)     :: ixx           ! counter x grid boxes
    integer(iintegers)     :: iyy           ! counter y grid boxes

    real(ireals)           :: Absup         ! Upwelling, Absorption, bottom
    real(ireals)           :: Absdn         ! Downwelling Absorption, top
    real(ireals)           :: Absup_s       ! Upwelling, Absorption, side
    real(ireals)           :: Absdn_s       ! Downwelling Absorption, side  
    real(ireals)           :: Emup          ! Upwelling, Emission, bottom
    real(ireals)           :: Emdn          ! Downwelling Emission, top
    real(ireals)           :: Emup_s        ! Upwelling, Emission, side
    real(ireals)           :: Emdn_s        ! Downwelling Emission, side  

    real(ireals)           :: HR_T          ! Heating Rate Top Contribution  
    real(ireals)           :: HR_s          ! Heating Rate Side Contribution 

    real(ireals)           :: mu            ! zenith angle

    real(ireals)           :: ax,ay         ! integration boarder for side contributions
    real(ireals)           :: bx,by         ! integration boarder for side contributions 
    real(ireals)           :: cx,cy         ! integration boarder for side contributions
    real(ireals)           :: az1,az        ! integration boarder for side contributions
    real(ireals)           :: bz1,bz        ! integration boarder for side contributions 
    real(ireals)           :: cz1,cz        ! integration boarder for side contributions 


    real(ireals)           :: l      ! Averaged Transmitted Flux, Top Contribution 
    real(ireals)           :: Trans  ! Transmitted Flux, Top Contribution

    real(ireals)           :: L_up_3d(lbound(Eup,1):ubound(Eup,1),lbound(Eup,2):ubound(Eup,2),lbound(Eup,3):ubound(Eup,3))
    real(ireals)           :: L_dn_3d(lbound(Edn,1):ubound(Edn,1),lbound(Edn,2):ubound(Edn,2),lbound(Edn,3):ubound(Edn,3))
    real(ireals)           :: B2   ! Planck General
    real(ireals)           :: B2_1 ! Planck Top Contribution

    integer(iintegers) :: is,ie,js,je,nlyr

    ! ## Fitting Parameters
    real(ireals)           :: afit1 
    real(ireals)           :: bfit1 
    real(ireals)           :: bfit2 
    real(ireals)           :: afit2 
    real(ireals)           :: cfit2 
    real(ireals)           :: factor2 
    real(ireals)           :: factor1 

    real(ireals),parameter :: eps = 0.0001
    real(ireals),parameter :: pi=3.141592653589793

    nlyr = ubound(hr,1)-1
    is = lbound(hr,2) +1
    ie = ubound(hr,2) -1
    js = lbound(hr,3) +1
    je = ubound(hr,3) -1
    

 
    hr = -9999999
    stop 'NCA not freely available.... please consult Carolin Klinger for an implementation.'

  end subroutine ts_nca


end module m_ts_nca
