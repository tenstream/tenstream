!-------------------------------------------------------------------------
! This file is part of the TenStream solver.
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

!> \page Routines to call tenstream with optical properties from RRTM
!! The routines in this file take care of concatenating the model input
!! with a background profile which is needed to do radiative transfer up
!! until TopOfAtmosphere
!!

module m_dyn_atm_to_rrtmg

  use iso_fortran_env, only: REAL32, REAL64
  use m_tenstr_parkind_sw, only: im => kind_im, rb => kind_rb

  use m_data_parameters, only : iintegers, mpiint, ireals, default_str_len, &
    zero

  use m_helper_functions, only: CHKERR, search_sorted_bisection

  use m_tenstream_interpolation, only : interp_1d

  implicit none

  logical,parameter :: ldebug=.True.
  !logical,parameter :: ldebug=.False.

  interface
    real function PLKINT(WVLLO, WVLHI, T)
      real :: WVLLO, WVLHI, T
    end function
  end interface

  interface hydrostat_dz
    module procedure hydrostat_dz_real32, hydrostat_dz_real64
  end interface

  interface rev_last_dim
    module procedure rev1d, rev2d
  end interface

  type t_atm
    real(ireals),allocatable :: plev   (:) ! dim(nlay+1)
    real(ireals),allocatable :: tlev   (:) !
    real(ireals),allocatable :: zt     (:) !
    real(ireals),allocatable :: h2o_lev(:) !
    real(ireals),allocatable :: o3_lev (:) !
    real(ireals),allocatable :: co2_lev(:) !
    real(ireals),allocatable :: ch4_lev(:) !
    real(ireals),allocatable :: n2o_lev(:) !
    real(ireals),allocatable :: o2_lev (:) !

    real(ireals),allocatable :: play   (:) ! dim(nlay)
    real(ireals),allocatable :: zm     (:) !
    real(ireals),allocatable :: dz     (:) !
    real(ireals),allocatable :: tlay   (:) !
    real(ireals),allocatable :: h2o_lay(:) !
    real(ireals),allocatable :: o3_lay (:) !
    real(ireals),allocatable :: co2_lay(:) !
    real(ireals),allocatable :: ch4_lay(:) !
    real(ireals),allocatable :: n2o_lay(:) !
    real(ireals),allocatable :: o2_lay (:) !
  end type

  contains

  pure elemental function hydrostat_dz_real32(dp, p, T)
    real(REAL32), intent(in) :: dp, p, T
    real(REAL32) :: hydrostat_dz_real32, rho
    rho = p / 287.058_REAL32 / T
    hydrostat_dz_real32 = dp / rho / 9.8065_REAL32
  end function
  pure elemental function hydrostat_dz_real64(dp, p, T)
    real(REAL64), intent(in) :: dp, p, T
    real(REAL64) :: hydrostat_dz_real64, rho
    rho = p / 287.058_REAL64 / T
    hydrostat_dz_real64 = dp / rho / 9.8065_REAL64
  end function

  subroutine hydrostat_lev(plev,tlay, hsrfc, hhl, dz)
    ! Integrate vertical height profile hydrostatically. arrays start at bottom(surface)
    real(ireals),intent(in) :: plev(:),tlay(:)
    real(ireals),intent(in) :: hsrfc
    real(ireals),intent(out) :: hhl(size(plev))
    real(ireals),intent(out) :: dz(size(tlay))
    integer(im) :: k
    hhl(1) = hsrfc
    do k=1,size(tlay)
      dz(k) = hydrostat_dz(abs(plev(k+1)-plev(k)), (plev(k+1)+plev(k))/2, tlay(k))
      hhl(k+1) = hhl(k) + dz(k)
    enddo
    if(any(dz.le.zero)) then
      print *,'plev',plev
      print *,'tlay',tlay
      print *,'dz',dz
      print *,'hhl',hhl
      stop 'error in dz'
    endif
  end subroutine


  subroutine sanitize_input(plev, tlev, tlay)
    real(ireals),intent(in),dimension(:) :: plev, tlev
    real(ireals),intent(in),dimension(:),optional :: tlay

    integer(mpiint) :: errcnt
    integer(iintegers) :: k
    logical :: lerr

    errcnt = 0
    lerr = maxval(plev) .gt. 1050
    if(lerr) then
      print *,'Pressure above 1050 hPa -- are you sure this is earth?', maxval(plev)
      errcnt = errcnt+1
    endif

    lerr = minval(plev) .lt. zero
    if(lerr) then
      print *,'Pressure negative -- are you sure this is physically correct?', minval(plev)
      errcnt = errcnt+1
    endif

    lerr = minval(tlev) .lt. 180
    if(lerr) then
      print *,'Temperature is very low -- are you sure RRTMG can handle that?', minval(tlev)
      errcnt = errcnt+1
    endif

    lerr = maxval(tlev) .gt. 400
    if(lerr) then
      print *,'Temperature is very high -- are you sure RRTMG can handle that?', maxval(tlev)
      errcnt = errcnt+1
    endif

    if(present(tlay) .and. ldebug) then
      do k=lbound(tlay,1), ubound(tlay,1)
        lerr = (tlay(k)-tlev(k).ge.zero) .eqv. (tlay(k)-tlev(k+1).gt.zero) ! different sign says its in between
        if(lerr) then
          print *,'Layer Temperature not between level temps?', k, tlev(k), '|', tlay(k), '|', tlev(k+1)
          errcnt = errcnt+1
        endif
      enddo
    endif


    if(errcnt.gt.0) then
      print *,'Found wonky input to pprts_rrtm_lw -- please check! -- will abort now.'
      call CHKERR(errcnt)
    endif

  end subroutine

  function rev2d(inp) ! reverse second dimension
    real(ireals),intent(in) :: inp(:,:)
    real(rb) :: rev2d(size(inp,1),size(inp,2))
    rev2d = inp(:,ubound(inp,2):lbound(inp,2):-1)
  end function

  function rev1d(inp) ! reverse array
    real(ireals),intent(in) :: inp(:)
    real(ireals) :: rev1d(size(inp,1))
    rev1d = inp(ubound(inp,1):lbound(inp,1):-1)
  end function

  ! merge the dynamics grid and the background profile together at lvl atm_ke
  ! NOTE! Only use with variables on layer
  subroutine merge_grid_var(a_hhl, d_hhl, atm_ke, a_lay, a_lev, col_var, d_var)
    integer(iintegers),intent(in) :: atm_ke
    real(ireals),intent(in) :: a_hhl(:), d_hhl(:), a_lay(:), a_lev(:) ! a_arr is from atm%, d_arr corresponds to dynamics grids
    real(rb),intent(out) :: col_var(:)
    real(ireals),intent(in),optional :: d_var(:)
    integer(iintegers) :: k, kt ! kt is reverse index
    real(ireals) :: h

    ! Top of atmosphere layers are always given by background profile
    do k=1,atm_ke
      kt = size(col_var)-k+1
      col_var(kt) = a_lay(k)
    enddo

    if(present(d_var)) then ! dynamics grid variable is provided, use that
      do k=1,size(d_var)
        col_var(k) = d_var(k)
      enddo
    else ! we may still use atmospheric grid file instead...
      do k=1,size(col_var)-atm_ke
        h = (d_hhl(k+1) + d_hhl(k)) / 2
        col_var(k) = interp_1d(search_sorted_bisection(a_hhl, h), a_lev)
      enddo
    endif

  end subroutine
end module

