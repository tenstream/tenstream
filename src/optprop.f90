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

use m_optprop_parameters, only : ldebug_optprop, Ndir_8_10,Ndiff_8_10, Ndir_1_2,Ndiff_1_2, coeff_mode
use m_helper_functions, only : rmse
use m_data_parameters, only: ireals,iintegers,one,zero,i0,i1,inil,mpiint
use m_optprop_LUT, only : t_optprop_LUT, t_optprop_LUT_1_2,t_optprop_LUT_8_10
use m_optprop_ANN, only : ANN_init, ANN_get_dir2dir, ANN_get_dir2diff, ANN_get_diff2diff

use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast

implicit none

private
public :: t_optprop_1_2,t_optprop_8_10

type,abstract :: t_optprop
  logical :: optprop_debug=ldebug_optprop
  integer(iintegers) :: dir_streams=inil,diff_streams=inil
  class(t_optprop_LUT),allocatable :: OPP_LUT
  contains
    procedure :: init
    procedure :: get_coeff
    procedure :: get_coeff_bmc
    procedure :: coeff_symmetry
    procedure :: destroy
end type

type,extends(t_optprop) :: t_optprop_1_2
end type

type,extends(t_optprop) :: t_optprop_8_10
end type

contains

  subroutine init(OPP, azis, szas, comm)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: szas(:),azis(:)
      integer(mpiint) ,intent(in) :: comm
      integer(mpiint) :: ierr

      select type(OPP)
        class is (t_optprop_1_2)
          OPP%dir_streams  =  Ndir_1_2
          OPP%diff_streams =  Ndiff_1_2

        class is (t_optprop_8_10)
          OPP%dir_streams  =  Ndir_8_10
          OPP%diff_streams =  Ndiff_8_10

        class default
        stop ' init optprop : unexpected type for optprop object!'
      end select

      select case (coeff_mode)
          case(i0) ! LookUpTable Mode
            select type(OPP)
              class is (t_optprop_1_2)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_1_2::OPP%OPP_LUT)

              class is (t_optprop_8_10)
               if(.not.allocated(OPP%OPP_LUT) ) allocate(t_optprop_LUT_8_10::OPP%OPP_LUT)

              class default
                stop ' init optprop : unexpected type for optprop object!'
            end select
            call OPP%OPP_LUT%init(azis, szas, comm)

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

  subroutine get_coeff_bmc(OPP, aspect, tauz, w0, g, dir, C, angles)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: aspect, tauz, w0, g
      logical,intent(in) :: dir
      real(ireals),intent(out):: C(:)
      real(ireals),intent(in),optional :: angles(2)

      real(ireals) :: S_diff(OPP%OPP_LUT%diff_streams),T_dir(OPP%OPP_LUT%dir_streams)
      real(ireals) :: S_tol (OPP%OPP_LUT%diff_streams),T_tol(OPP%OPP_LUT%dir_streams)
      integer(iintegers) :: isrc

      if(present(angles)) then
        if(dir) then !dir2dir
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, aspect, tauz, w0, g, .True., angles(1), angles(2), -1_mpiint, S_diff, T_dir, S_tol, T_tol)
            C((isrc-1)*OPP%OPP_LUT%dir_streams+1:isrc*OPP%OPP_LUT%dir_streams) = T_dir
          enddo
        else ! dir2diff
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc, aspect, tauz, w0, g, .True., angles(1), angles(2), -1_mpiint, S_diff, T_dir, S_tol, T_tol)
            C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams) = S_diff
          enddo
        endif
      else
        ! diff2diff
        do isrc=1,OPP%OPP_LUT%diff_streams
          call OPP%OPP_LUT%bmc_wrapper(isrc, aspect, tauz, w0, g, .False., zero, zero, -1_mpiint, S_diff, T_dir, S_tol, T_tol)
          C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams) = S_diff
        enddo
      endif ! angles_present

  end subroutine

  subroutine get_coeff(OPP, aspect, tauz, w0, g, dir, C, inp_angles)
        class(t_optprop) :: OPP
        logical,intent(in) :: dir
        real(ireals),intent(in) :: aspect, tauz, w0, g
        real(ireals),intent(in),optional :: inp_angles(2)
        real(ireals),intent(out):: C(:)

        real(ireals) :: angles(2)

        logical,parameter :: compute_coeff_online=.False.

        if(compute_coeff_online) then
          call get_coeff_bmc(OPP, aspect, tauz, w0, g, dir, C, inp_angles)
          return
        endif

        if(ldebug_optprop) call check_inp(aspect, tauz, w0, g, dir, C)

        if(present(inp_angles)) then
          angles = inp_angles
          if(inp_angles(2).le.1e-3_ireals) angles(1) = zero ! if sza is close to 0, azimuth is symmetric -> dont need to distinguish
        endif

        select case (coeff_mode)

          case(i0) ! LookUpTable Mode

            if(present(inp_angles)) then ! obviously we want the direct coefficients
              if(dir) then ! dir2dir
                call OPP%OPP_LUT%LUT_get_dir2dir(aspect, tauz, w0, g, angles(1), angles(2), C)
              else         ! dir2diff
                call OPP%OPP_LUT%LUT_get_dir2diff(aspect, tauz, w0, g, angles(1), angles(2), C)
              endif
            else
              ! diff2diff
              call OPP%OPP_LUT%LUT_get_diff2diff(aspect, tauz, w0, g, C)
            endif


          case(i1) ! ANN

            if(present(inp_angles)) then ! obviously we want the direct coefficients
              if(dir) then ! specifically the dir2dir
                call ANN_get_dir2dir(aspect, tauz, w0, g, angles(1), angles(2), C)
              else ! dir2diff
                call ANN_get_dir2diff(aspect, tauz, w0, g, angles(1), angles(2), C)
              endif
            else
              ! diff2diff
              call ANN_get_diff2diff(aspect, tauz, w0, g, C)
            endif

          case default
            stop 'coeff mode optprop initialization not defined '
        end select

      contains

        subroutine check_inp(aspect, tauz, w0, g, dir, C)
            real(ireals),intent(in) :: aspect, tauz, w0, g
            logical,intent(in) :: dir
            real(ireals),intent(in):: C(:)
            if(OPP%optprop_debug) then
              if( (any([aspect, tauz, w0, g].lt.zero)) .or. (any(isnan([aspect, tauz, w0, g]))) ) then
                print *,'optprop_lookup_coeff :: corrupt optical properties: bg:: ',[aspect, tauz, w0, g]
                call exit
              endif
            endif
            if(present(inp_angles)) then
              if(dir .and. size(C).ne. OPP%dir_streams**2) then
                print *,'direct called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%dir_streams**2
              endif
              if(.not.dir .and. size(C).ne. OPP%diff_streams*OPP%dir_streams) then
                print *,'dir2diffuse called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%diff_streams*OPP%dir_streams
              endif
            else
              if(dir .and. size(C).ne. OPP%diff_streams) then
                print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%diff_streams
              endif
              if(.not.dir .and. size(C).ne. OPP%diff_streams**2) then
                print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%diff_streams**2
              endif
            endif
        end subroutine

end subroutine
        function coeff_symmetry(OPP, isrc,coeff)
            class(t_optprop) :: OPP
            real(ireals) :: coeff_symmetry(OPP%diff_streams)
            real(ireals),intent(in) :: coeff(:)
            integer(iintegers),intent(in) :: isrc
            integer(iintegers),parameter :: l=1
            integer(iintegers),parameter :: k=5
            !TODO: this was just a simple test if we can enhance diffusion
            !      artificially by fiddling with the transport coefficients. -- only marginally improvement for Enet(surface) but increased
            !      rmse in atmosphere...
            real(ireals),parameter :: artificial_diffusion = zero

            ! integer,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9
            select type(OPP)

              class is (t_optprop_1_2)
                select case (isrc)
                  case(1)
                          coeff_symmetry = coeff([l+0, l+1])
                  case(2)
                          coeff_symmetry = coeff([l+1, l+0])
                  case default
                          stop 'cant call coeff_symmetry with isrc -- error happened for type optprop_1_2'
                end select

              class is (t_optprop_8_10)
                select case (isrc)
                  case(1)
                    coeff_symmetry = coeff([l+0, l+1, l+2, l+2, l+3, l+3, l+2, l+2, l+3, l+3])
                  case(2)
                    coeff_symmetry = coeff([l+1, l+0, l+3, l+3, l+2, l+2, l+3, l+3, l+2, l+2])
                  case(3)
                    coeff_symmetry = coeff([k+0, k+1, k+2, k+3, k+4, k+5, k+6, k+6, k+7, k+7])
                  case(4)
                    coeff_symmetry = coeff([k+0, k+1, k+3, k+2, k+5, k+4, k+6, k+6, k+7, k+7])
                  case(5)
                    coeff_symmetry = coeff([k+1, k+0, k+4, k+5, k+2, k+3, k+7, k+7, k+6, k+6])
                  case(6)
                    coeff_symmetry = coeff([k+1, k+0, k+5, k+4, k+3, k+2, k+7, k+7, k+6, k+6])
                  case(7)
                    coeff_symmetry = coeff([k+0, k+1, k+6, k+6, k+7, k+7, k+2, k+3, k+4, k+5])
                  case(8)
                    coeff_symmetry = coeff([k+0, k+1, k+6, k+6, k+7, k+7, k+3, k+2, k+5, k+4])
                  case(9)
                    coeff_symmetry = coeff([k+1, k+0, k+7, k+7, k+6, k+6, k+4, k+5, k+2, k+3])
                  case(10)
                    coeff_symmetry = coeff([k+1, k+0, k+7, k+7, k+6, k+6, k+5, k+4, k+3, k+2])
                  case default
                    stop 'cant call coeff_symmetry with isrc -- error happened for optprop_8_10'
                end select

                if(artificial_diffusion.gt.zero) then
                  if(isrc.eq.1 .or. isrc.eq.2) then
                    coeff_symmetry(3:4) = coeff_symmetry(3:4) + coeff_symmetry(2)*artificial_diffusion *.25_ireals
                    coeff_symmetry(7:8) = coeff_symmetry(7:8) + coeff_symmetry(2)*artificial_diffusion *.25_ireals
                    coeff_symmetry(2)   = coeff_symmetry(2)   - coeff_symmetry(2)*artificial_diffusion

                    coeff_symmetry(5: 6) = coeff_symmetry(5: 6) + coeff_symmetry(1)*artificial_diffusion *.25_ireals
                    coeff_symmetry(9:10) = coeff_symmetry(9:10) + coeff_symmetry(1)*artificial_diffusion *.25_ireals
                    coeff_symmetry(1)    = coeff_symmetry(1)    - coeff_symmetry(1)*artificial_diffusion
                  endif
                endif

              class default
                stop 'coeff_symmetry : unexpected type for OPP !'
            end select


            if(ldebug_optprop) then
              if(real(sum(coeff_symmetry)).gt.one+epsilon(one)*10._ireals) then
                print *,'sum of diffuse coeff_symmetrys bigger one!',sum(coeff_symmetry),'for src=',isrc,'coeff_symmetry:',coeff_symmetry
                call exit()
              endif
            endif
        end function

end module
