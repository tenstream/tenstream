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

module m_optprop

#ifdef _XLF
      use ieee_arithmetic
#define isnan ieee_is_nan
#endif

use m_optprop_parameters, only : ldebug_optprop, coeff_mode
use m_helper_functions, only : rmse, CHKERR, itoa, ftoa, approx
use m_data_parameters, only: ireals,iintegers,one,zero,i0,i1,inil,mpiint
use m_optprop_LUT, only : t_optprop_LUT, t_optprop_LUT_1_2,t_optprop_LUT_8_10, t_optprop_LUT_3_6, t_optprop_LUT_3_10, &
  t_optprop_LUT_wedge_5_8
use m_optprop_ANN, only : ANN_init, ANN_get_dir2dir, ANN_get_dir2diff, ANN_get_diff2diff
use m_boxmc_geometry, only : setup_default_unit_cube_geometry

use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast

implicit none

private
public :: t_optprop, t_optprop_1_2,t_optprop_8_10, t_optprop_3_6, t_optprop_3_10, t_optprop_wedge_5_8

type,abstract :: t_optprop
  logical :: optprop_debug=ldebug_optprop
  class(t_optprop_LUT), allocatable :: OPP_LUT
  contains
    procedure :: init
    procedure :: get_coeff
    procedure :: get_coeff_bmc
    procedure :: dir2dir_coeff_symmetry
    procedure :: dir2diff_coeff_symmetry
    procedure :: destroy
end type

type,extends(t_optprop) :: t_optprop_1_2
end type

type,extends(t_optprop) :: t_optprop_8_10
end type

type,extends(t_optprop) :: t_optprop_3_6
end type

type,extends(t_optprop) :: t_optprop_3_10
end type

type,extends(t_optprop) :: t_optprop_wedge_5_8
end type

contains

  subroutine init(OPP, comm)
      class(t_optprop), intent(inout) :: OPP
      integer(mpiint) ,intent(in) :: comm
      integer(mpiint) :: ierr

      select case (coeff_mode)
          case(i0) ! LookUpTable Mode
            select type(OPP)
              class is (t_optprop_1_2)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_1_2::OPP%OPP_LUT)

              class is (t_optprop_8_10)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_10::OPP%OPP_LUT)

              class is (t_optprop_3_6)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_6::OPP%OPP_LUT)

              class is (t_optprop_3_10)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_3_10::OPP%OPP_LUT)

              class is (t_optprop_wedge_5_8)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_wedge_5_8::OPP%OPP_LUT)
              class default
                stop ' init optprop : unexpected type for optprop object!'
            end select
            call OPP%OPP_LUT%init(comm)

          case(i1) ! ANN
            call ANN_init(comm, ierr)
  !          stop 'ANN not yet implemented'
          case default
            stop 'coeff mode optprop initialization not defined '
        end select
  end subroutine
  subroutine destroy(OPP)
      class(t_optprop) :: OPP
      if(allocated(OPP%OPP_LUT)) then
          call OPP%OPP_LUT%destroy()
          deallocate(OPP%OPP_LUT)
      endif
  end subroutine

  subroutine get_coeff_bmc(OPP, vertices, tauz, w0, g, dir, C, angles)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: tauz, w0, g, vertices(:)
      logical,intent(in) :: dir
      real(ireals),intent(out):: C(:)
      real(ireals),intent(in),optional :: angles(2)

      real(ireals) :: S_diff(OPP%OPP_LUT%diff_streams),T_dir(OPP%OPP_LUT%dir_streams)
      real(ireals) :: S_tol (OPP%OPP_LUT%diff_streams),T_tol(OPP%OPP_LUT%dir_streams)
      integer(iintegers) :: isrc


      if(present(angles)) then
        if(dir) then !dir2dir
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, vertices, tauz, w0, g, .True., angles(1), angles(2), -1_mpiint, S_diff, T_dir, S_tol, T_tol)
            C((isrc-1)*OPP%OPP_LUT%dir_streams+1:isrc*OPP%OPP_LUT%dir_streams) = T_dir
          enddo
        else ! dir2diff
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, vertices, tauz, w0, g, .True., angles(1), angles(2), -1_mpiint, S_diff, T_dir, S_tol, T_tol)
            C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams) = S_diff
          enddo
        endif
      else
        ! diff2diff
        do isrc=1,OPP%OPP_LUT%diff_streams
          call OPP%OPP_LUT%bmc_wrapper(isrc, vertices, tauz, w0, g, .False., zero, zero, -1_mpiint, S_diff, T_dir, S_tol, T_tol)
          C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams) = S_diff
        enddo
      endif ! angles_present

  end subroutine

  subroutine get_coeff(OPP, tauz, w0, g, aspect_zx, dir, C, inp_angles, lswitch_east, lswitch_north, wedge_coords)
        class(t_optprop)                  :: OPP
        logical,intent(in)                :: dir
        real(ireals),intent(in)           :: tauz, w0, g, aspect_zx
        real(ireals),intent(in),optional  :: inp_angles(:)         ! phi and azimuth in degree
        logical,intent(in),optional       :: lswitch_east, lswitch_north
        real(ireals),intent(in),optional  :: wedge_coords(:)       ! 6 coordinates of wedge triangle, only used for wedge OPP types
        real(ireals),intent(out)          :: C(:)

        select type (OPP)
        class is (t_optprop_1_2)
          call boxmc_lut_call(OPP, tauz, w0, g, aspect_zx, dir, C, inp_angles, lswitch_east, lswitch_north)

        class is (t_optprop_3_6)
          call boxmc_lut_call(OPP, tauz, w0, g, aspect_zx, dir, C, inp_angles, lswitch_east, lswitch_north)

        class is (t_optprop_3_10)
          call boxmc_lut_call(OPP, tauz, w0, g, aspect_zx, dir, C, inp_angles, lswitch_east, lswitch_north)

        class is (t_optprop_8_10)
          call boxmc_lut_call(OPP, tauz, w0, g, aspect_zx, dir, C, inp_angles, lswitch_east, lswitch_north)

        class is (t_optprop_wedge_5_8)
          call wedge_lut_call(OPP, tauz, w0, g, aspect_zx, dir, C, inp_angles, wedge_coords)

        class default
          call CHKERR(1_mpiint, 'initialize LUT: unexpected type for optprop object!')
      end select

      contains
        subroutine wedge_lut_call(OPP, tauz, w0, g, aspect_zx, dir, C, inp_angles, wedge_coords)
          class(t_optprop)                  :: OPP
          logical,intent(in)                :: dir
          real(ireals),intent(in)           :: tauz, w0, g, aspect_zx
          real(ireals),intent(in),optional  :: inp_angles(:)
          real(ireals),intent(in),optional  :: wedge_coords(:)       ! 6 coordinates of wedge triangle, only used for wedge OPP types, have to be in local coords already (i.e. A=[0,0], B=[0,1], C=[...]
          real(ireals),intent(out)          :: C(:)

          real(ireals) :: angles(2), C_pnt(2)

          if(ldebug_optprop) then
            if(.not.present(wedge_coords)) call CHKERR(1_mpiint, 'If you use an wedge OPP object, I highly recommend that you provide the wedge coordinates')
          endif

          if(ldebug_optprop) then
            if(.not.all(approx(wedge_coords([1,2,4]), zero)) &
              .or. .not.approx(wedge_coords(3), one)) then
              print *,'wedge_coords:', wedge_coords
              call CHKERR(1_mpiint, 'provided wedge coords have to in local bmc coordinate system!')
            endif
            if(wedge_coords(5).lt.0.35_ireals .or. wedge_coords(5).gt.0.65_ireals) &
              call CHKERR(1_mpiint, 'wedge_coords(5) is outside the range we expected... take care of it somehow, do some magic! '//ftoa(wedge_coords(5)))
            call check_inp(tauz, w0, g, aspect_zx, dir, C)
          endif

          C_pnt = wedge_coords(5:6)

          select case (coeff_mode)

          case(i0) ! LookUpTable Mode

            if(present(inp_angles)) then ! obviously we want the direct coefficients
              if(dir) then ! dir2dir
                call OPP%OPP_LUT%LUT_get_dir2dir([tauz, w0, g, aspect_zx, C_pnt(1), C_pnt(2), angles(1), angles(2)], C)
              else         ! dir2diff
                call OPP%OPP_LUT%LUT_get_dir2diff([tauz, w0, g, aspect_zx, C_pnt(1), C_pnt(2), angles(1), angles(2)], C)
              endif
            else
              ! diff2diff
              call OPP%OPP_LUT%LUT_get_diff2diff([tauz, w0, g, aspect_zx, C_pnt(1), C_pnt(2)], C)
            endif

          case default
            call CHKERR(1_mpiint, 'particular value of coeff mode in optprop_parameters is not defined: '//itoa(coeff_mode))
          end select

        end subroutine

        subroutine boxmc_lut_call(OPP, tauz, w0, g, aspect_zx, dir, C, angles, lswitch_east, lswitch_north)
          class(t_optprop)                  :: OPP
          logical,intent(in)                :: dir
          real(ireals),intent(in)           :: tauz, w0, g, aspect_zx
          real(ireals),intent(in),optional  :: angles(:)
          logical,intent(in),optional       :: lswitch_east, lswitch_north
          real(ireals),intent(out)          :: C(:)

          logical,parameter :: compute_coeff_online=.False.
          real(ireals), allocatable :: vertices(:)

          if(compute_coeff_online) then
            call setup_default_unit_cube_geometry(one, one, aspect_zx, vertices)
            call get_coeff_bmc(OPP, vertices, tauz, w0, g, dir, C, inp_angles)
            return
          endif

          if(ldebug_optprop) call check_inp(tauz, w0, g, aspect_zx, dir, C)

          !if(present(inp_angles)) then
          !  angles = inp_angles
          !  if(inp_angles(2).le.1e-3_ireals) angles(1) = zero ! if sza is close to 0, azimuth is symmetric -> dont need to distinguish
          !endif

          select case (coeff_mode)

          case(i0) ! LookUpTable Mode

            if(present(angles)) then ! obviously we want the direct coefficients
              if(dir) then ! dir2dir
                call OPP%OPP_LUT%LUT_get_dir2dir([tauz, w0, g, aspect_zx, angles(1), angles(2)], C)
                call OPP%dir2dir_coeff_symmetry(C, lswitch_east, lswitch_north)
              else         ! dir2diff
                call OPP%OPP_LUT%LUT_get_dir2diff([tauz, w0, g, aspect_zx, angles(1), angles(2)], C)
                call OPP%dir2diff_coeff_symmetry(C, lswitch_east, lswitch_north)
              endif
            else
              ! diff2diff
              call OPP%OPP_LUT%LUT_get_diff2diff([tauz, w0, g, aspect_zx], C)
            endif


          case(i1) ! ANN

            if(present(angles)) then ! obviously we want the direct coefficients
              if(dir) then ! specifically the dir2dir
                call ANN_get_dir2dir(tauz, w0, g, aspect_zx, angles(1), angles(2), C)
              else ! dir2diff
                call ANN_get_dir2diff(tauz, w0, g, aspect_zx, angles(1), angles(2), C)
              endif
            else
              ! diff2diff
              call ANN_get_diff2diff(tauz, w0, g, aspect_zx, C)
            endif

          case default
            call CHKERR(1_mpiint, 'particular value of coeff mode in optprop_parameters is not defined: '//itoa(coeff_mode))
          end select

        end subroutine

        subroutine check_inp(tauz, w0, g, aspect_zx, dir, C)
            real(ireals),intent(in) :: tauz, w0, g, aspect_zx
            logical,intent(in) :: dir
            real(ireals),intent(in):: C(:)
            if(OPP%optprop_debug) then
              if( (any([aspect_zx, tauz, w0, g].lt.zero)) .or. (any(isnan([aspect_zx, tauz, w0, g]))) ) then
                print *,'optprop_lookup_coeff :: corrupt optical properties: bg:: ',[aspect_zx, tauz, w0, g]
                call exit
              endif
            endif
            if(present(inp_angles)) then
              if(dir .and. size(C).ne. OPP%OPP_LUT%dir_streams**2) then
                print *,'direct called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%dir_streams**2
              endif
              if(.not.dir .and. size(C).ne. OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams) then
                print *,'dir2diffuse called get_coeff with wrong shaped output array:',size(C),'should be',OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams
              endif
            else
              if(dir .and. size(C).ne. OPP%OPP_LUT%diff_streams) then
                print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams
              endif
              if(.not.dir .and. size(C).ne. OPP%OPP_LUT%diff_streams**2) then
                print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams**2
              endif
            endif
        end subroutine

  end subroutine

  subroutine dir2diff_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)

    class(t_optprop) :: OPP
    logical, intent(in),optional :: lswitch_east, lswitch_north
    real(ireals),intent(inout) :: coeff(:)

    integer(iintegers)              :: dof
    real(ireals)                    :: newcoeff(size(coeff))

    select type(OPP)
      class is (t_optprop_1_2)
        continue

      class is (t_optprop_3_6)
        !for solver_3_6 only the offset is changing in those sides which should be switched
        dof = 3
        if(present(lswitch_east)) then
          if(lswitch_east) then
            newcoeff = coeff
            coeff(7:9)   = newcoeff([1, 2, 3] + dof*3)
            coeff(10:12) = newcoeff([1, 2, 3] + dof*2)
          endif
        endif
        if(present(lswitch_north)) then
          if(lswitch_north) then
            newcoeff = coeff
            coeff(13:15) = newcoeff([1, 2, 3] + dof*5)
            coeff(16:18) = newcoeff([1, 2, 3] + dof*4)
          endif
        endif

      class is (t_optprop_8_10)
        !for solver_8_10 the offset is chaning and the destination order
        dof = 8
        if(present(lswitch_east)) then
          if(lswitch_east) then
            newcoeff = coeff
            coeff(1:8)   = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]        )
            coeff(9:16)  = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*1 )
            coeff(17:24) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*3 )
            coeff(25:32) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*2 )
            coeff(33:40) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*5 )
            coeff(41:48) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*4 )
            coeff(49:56) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*6 )
            coeff(57:64) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*7 )
            coeff(65:72) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*8 )
            coeff(73:80) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8] +dof*9 )
          endif
        endif
        if(present(lswitch_north)) then
          if (lswitch_north) then
            newcoeff =coeff
            coeff(1:8)   = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]        )
            coeff(9:16)  = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*1 )
            coeff(17:24) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*3 )
            coeff(25:32) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*2 )
            coeff(33:40) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*5 )
            coeff(41:48) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*4 )
            coeff(49:56) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*6 )
            coeff(57:64) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*7 )
            coeff(65:72) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*8 )
            coeff(73:80) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8] +dof*9 )
          endif
        endif

      class is (t_optprop_3_10)
        !for solver_3_10 the offset is chaning and the destination order
        dof = 3
        if(present(lswitch_east)) then
          if(lswitch_east) then
            newcoeff = coeff
            coeff(1:3)   = newcoeff([1, 2, 3]        )
            coeff( 4: 6) = newcoeff([1, 2, 3] + dof*1)
            coeff( 7: 9) = newcoeff([1, 2, 3] + dof*3)
            coeff(10:12) = newcoeff([1, 2, 3] + dof*2)
            coeff(13:15) = newcoeff([1, 2, 3] + dof*5)
            coeff(16:18) = newcoeff([1, 2, 3] + dof*4)
            coeff(19:21) = newcoeff([1, 2, 3] + dof*6)
            coeff(22:24) = newcoeff([1, 2, 3] + dof*7)
            coeff(25:27) = newcoeff([1, 2, 3] + dof*8)
            coeff(28:30) = newcoeff([1, 2, 3] + dof*9)
          endif
        endif
        if(present(lswitch_north)) then
          if(lswitch_north) then
            newcoeff = coeff
            coeff(1:3)   = newcoeff([1, 2, 3]        )
            coeff( 4: 6) = newcoeff([1, 2, 3] + dof*1)
            coeff( 7: 9) = newcoeff([1, 2, 3] + dof*3)
            coeff(10:12) = newcoeff([1, 2, 3] + dof*2)
            coeff(13:15) = newcoeff([1, 2, 3] + dof*5)
            coeff(16:18) = newcoeff([1, 2, 3] + dof*4)
            coeff(19:21) = newcoeff([1, 2, 3] + dof*6)
            coeff(22:24) = newcoeff([1, 2, 3] + dof*7)
            coeff(25:27) = newcoeff([1, 2, 3] + dof*8)
            coeff(28:30) = newcoeff([1, 2, 3] + dof*9)
          endif
        endif
    end select

  end subroutine

  subroutine dir2dir_coeff_symmetry(OPP, coeff, lswitch_east, lswitch_north)
    class(t_optprop)            :: OPP
    logical, intent(in),optional:: lswitch_east, lswitch_north
    real(ireals),intent(inout)  :: coeff(:)

    integer(iintegers) :: dof
    real(ireals)       :: newcoeff(size(coeff))

    select type(OPP)
      class is (t_optprop_1_2)
        continue

      class is (t_optprop_3_6)
        !nothing to do because of symmetrie
        continue

      class is (t_optprop_3_10)
        !nothing to do because of symmetrie
        continue

      class is (t_optprop_8_10)
        dof = 8
        if(present(lswitch_east)) then
          if(lswitch_east) then
            newcoeff = coeff
            coeff(1:8)   = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof  )
            coeff(9:16)  = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]      )
            coeff(17:24) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*3)
            coeff(25:32) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*2)
            coeff(33:40) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*4)
            coeff(41:48) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*5)
            coeff(49:56) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*6)
            coeff(57:64) = newcoeff([2, 1, 4, 3, 5, 6, 7, 8]+dof*7)
          endif
        endif
        if(present(lswitch_north)) then
          if (lswitch_north) then
            newcoeff = coeff
            coeff(1:8)   = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*2)
            coeff(9:16)  = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*3)
            coeff(17:24) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]      )
            coeff(25:32) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof  )
            coeff(33:40) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*4)
            coeff(41:48) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*5)
            coeff(49:56) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*6)
            coeff(57:64) = newcoeff([3, 4, 1, 2, 5, 6, 7, 8]+dof*7)
          endif
        endif
    end select
  end subroutine
end module
