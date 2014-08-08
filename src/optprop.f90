module m_optprop
use m_optprop_parameters, only : ldebug_optprop
use m_helper_functions, only : rmse
use m_data_parameters, only: ireals,iintegers,one,zero,i0,i1,mpiint
use m_optprop_LUT, only : t_optprop_LUT, t_optprop_LUT_1_2,t_optprop_LUT_8_10

use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast

implicit none

private
public :: t_optprop_1_2,t_optprop_8_10

integer(mpiint) :: ierr

type,abstract :: t_optprop
  logical :: optprop_debug=.True.
  real(ireals) :: dx,dy
  integer(iintegers) :: coeff_mode=i0 ! 0 is LUT, 1 is Neural Network
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

  subroutine init(OPP, dx_inp,dy_inp, azis, szas, comm)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: szas(:),azis(:),dx_inp,dy_inp
      integer(mpiint) ,intent(in) :: comm

      OPP%dx=dx_inp
      OPP%dy=dy_inp

      select case (OPP%coeff_mode)
        case(i0) ! LookUpTable Mode
          select type(OPP)
            class is (t_optprop_1_2)
              allocate(t_optprop_LUT_1_2::OPP%OPP_LUT)

            class is (t_optprop_8_10)
              allocate(t_optprop_LUT_8_10::OPP%OPP_LUT)

            class default
              stop ' init optprop : unexpected type for optprop object!'
          end select
          call OPP%OPP_LUT%init(OPP%dx,OPP%dy, azis, szas, comm)

        case(i1) ! ANN
          stop 'ANN not yet implemented'
        case default
          stop 'coeff mode optprop initialization not defined ' 
      end select

  end subroutine
  subroutine destroy(OPP)
      class(t_optprop) :: OPP
      if(allocated(OPP%OPP_LUT)) deallocate(OPP%OPP_LUT)
  end subroutine

  subroutine get_coeff_bmc(OPP, dz,kabs,ksca,g,dir,C,angles)
      class(t_optprop) :: OPP
      real(ireals),intent(in) :: dz,g,kabs,ksca
      logical,intent(in) :: dir
      real(ireals),intent(out):: C(:)
      real(ireals),intent(in),optional :: angles(2)

      real(ireals) :: S_diff(OPP%OPP_LUT%diff_streams),T_dir(OPP%OPP_LUT%dir_streams)
      integer(iintegers) :: isrc

      if(present(angles)) then 
        if(dir) then !dir2dir
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper( isrc,OPP%dx,OPP%dy,dz,kabs,ksca,g,.True.,angles(1),angles(2),-1_mpiint,S_diff,T_dir)
            C((isrc-1)*OPP%OPP_LUT%dir_streams+1:isrc*OPP%OPP_LUT%dir_streams) = T_dir
          enddo
        else ! dir2diff
          do isrc=1,OPP%OPP_LUT%dir_streams
            call OPP%OPP_LUT%bmc_wrapper(isrc,OPP%dx,OPP%dy,dz,kabs,ksca,g,.True.,angles(1),angles(2),-1_mpiint,S_diff,T_dir)
            C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams) = S_diff
          enddo
        endif
      else
        ! diff2diff
        do isrc=1,OPP%OPP_LUT%diff_streams
          call OPP%OPP_LUT%bmc_wrapper(isrc,OPP%dx,OPP%dy,dz,kabs,ksca,g,.False.,zero,zero,-1_mpiint,S_diff,T_dir)
          C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams) = S_diff
        enddo
      endif ! angles_present

  end subroutine

  subroutine get_coeff(OPP, dz,kabs,ksca,g,dir,C,inp_angles)
        class(t_optprop) :: OPP
        logical,intent(in) :: dir
        real(ireals),intent(in) :: dz,g,kabs,ksca
        real(ireals),intent(in),optional :: inp_angles(2)
        real(ireals),intent(out):: C(:)

        real(ireals),allocatable :: C_diff(:)
        real(ireals) :: angles(2)
        integer(iintegers) :: isrc

! This enables on-line calculations of coefficients with bmc code. This takes FOREVER! - use this only to check if LUT is working correctly!
        logical,parameter :: determine_coeff_error=.False.
        logical,parameter :: compute_coeff_online=.False.
        real(ireals) :: S_diff(OPP%OPP_LUT%diff_streams),T_dir(OPP%OPP_LUT%dir_streams)
        real(ireals) :: frmse(2)
        real(ireals),parameter :: checking_limit=1e-1

        real(ireals) :: dx,dy,diff_streams,dir_streams

        if(compute_coeff_online) then
          call get_coeff_bmc(OPP, dz,kabs,ksca,g,dir,C,inp_angles)
          return
        endif

        dx = OPP%dx
        dy = OPP%dy
        diff_streams= OPP%OPP_LUT%diff_streams
        dir_streams = OPP%OPP_LUT%dir_streams

        if(ldebug_optprop) then
          if(OPP%optprop_debug) then
            if( (any([dz,kabs,ksca,g].lt.zero)) .or. (any(isnan([dz,kabs,ksca,g]))) ) then
              print *,'optprop_lookup_coeff :: corrupt optical properties: bg:: ',[dz,kabs,ksca,g]
              call exit
            endif
          endif
          if(present(inp_angles)) then
            if(dir .and. size(C).ne. OPP%OPP_LUT%dir_streams**2) then
              print *,'direct called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%dir_streams
            endif
            if(.not.dir .and. size(C).ne. OPP%OPP_LUT%diff_streams*OPP%OPP_LUT%dir_streams) then
              print *,'dir2diffuse called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams
            endif
          else
            if(size(C).ne. OPP%OPP_LUT%diff_streams**2) then
              print *,'diff2diff called get_coeff with wrong shaped output array:',size(C),'should be ',OPP%OPP_LUT%diff_streams
            endif
          endif
        endif

        if(present(inp_angles)) then
          angles = inp_angles
          if(inp_angles(2).le.1e-3_ireals) angles(1) = zero ! if sza is close to 0, azimuth is symmetric -> dont need to distinguish
        endif

          if(present(inp_angles)) then ! obviously we want the direct coefficients
            if(dir) then ! specifically the dir2dir
              call OPP%OPP_LUT%LUT_get_dir2dir(dz,kabs,ksca,g,angles(1),angles(2),C)
            else ! dir2diff
              call OPP%OPP_LUT%LUT_get_dir2diff(dz,kabs,ksca,g,angles(1),angles(2),C)
            endif
          else
            ! diff2diff
            call OPP%OPP_LUT%LUT_get_diff2diff(dz,kabs,ksca,g,C_diff)
          endif

        if(.not.present(inp_angles)) then
          do isrc=1,OPP%OPP_LUT%diff_streams
            C( (isrc-1)*OPP%OPP_LUT%diff_streams+i1 : isrc*OPP%OPP_LUT%diff_streams ) = OPP%coeff_symmetry(isrc, C_diff )
          enddo
          deallocate(C_diff)
        endif

        if(determine_coeff_error) then 
                call random_number(T_dir(1)) 
                if(T_dir(1).le.checking_limit) then
                        if(present(inp_angles)) then 
                                if(dir) then !dir2dir
                                        do isrc=1,OPP%OPP_LUT%dir_streams
                                                call OPP%OPP_LUT%bmc_wrapper( isrc,dx,dy,dz,kabs,ksca,g,.True.,angles(1),angles(2),-1_mpiint,S_diff,T_dir)
                                                frmse = RMSE(C((isrc-1)*OPP%OPP_LUT%dir_streams+1:isrc*OPP%OPP_LUT%dir_streams), T_dir)
                                                print "('check ',i1,' dir2dir ',6e10.2,' :: RMSE ',2e13.4,' coeff ',8f7.4,' bmc ',8f7.4)",&
                                                        isrc,dz,kabs,ksca,g,angles(1),angles(2),frmse, C((isrc-1)*OPP%OPP_LUT%dir_streams+1:isrc*OPP%OPP_LUT%dir_streams),T_dir
                                                if(all(frmse.gt..1_ireals) ) stop 'Something terrible happened... I checked the coefficients and they differ more than they should!'
                                        enddo
                                else ! dir2diff
                                        do isrc=1,OPP%OPP_LUT%dir_streams
                                                call OPP%OPP_LUT%bmc_wrapper(isrc,dx,dy,dz,kabs,ksca,g,.True.,angles(1),angles(2),-1_mpiint,S_diff,T_dir)
                                                frmse = RMSE(C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams), S_diff)
                                                print "('check ',i1,' dir2diff ',6e10.2,' :: RMSE ',2e13.4,' coeff ',10e10.2)",&
                                                        isrc,dz,kabs,ksca,g,angles(1),angles(2),frmse ,abs(C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams)-S_diff)
                                                print "(a10,e13.4,10e13.4)",'C_interp',sum(C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams)),C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams)
                                                print "(a10,e13.4,10e13.4)",'C_bmc   ',sum(S_diff),                  S_diff
                                                print *,''
                                                if(all(frmse.gt..1_ireals) ) stop 'Something terrible happened... I checked the coefficients and they differ more than they should!'
                                        enddo
                                endif
                        else
                                ! diff2diff
                                do isrc=1,OPP%OPP_LUT%diff_streams
                                        call OPP%OPP_LUT%bmc_wrapper(isrc,dx,dy,dz,kabs,ksca,g,.False.,zero,zero,-1_mpiint,S_diff,T_dir)
                                        frmse = RMSE(C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams), S_diff)
                                        print "('check ',i1,' diff2diff ',4e10.2,' :: RMSE ',2e13.4,' coeff err',10e10.2)",&
                                                isrc,dz,kabs,ksca,g,frmse,abs(C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams)-S_diff)
                                        print "(a10,e13.4,10e13.4)",'C_interp',sum(C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams)),C((isrc-1)*OPP%OPP_LUT%diff_streams+1:isrc*OPP%OPP_LUT%diff_streams)
                                        print "(a10,e13.4,10e13.4)",'C_bmc   ',sum(S_diff),S_diff
                                        print *,''
                                        if(all(frmse.gt..1_ireals) ) stop 'Something terrible happened... I checked the coefficients and they differ more than they should!'
                                enddo
                        endif ! angles_present
                endif ! is in checking_limit
        endif ! want to check
end subroutine
        function coeff_symmetry(OPP, isrc,coeff)
            class(t_optprop) :: OPP
            real(ireals) :: coeff_symmetry(OPP%OPP_LUT%diff_streams)
            real(ireals),intent(in) :: coeff(:)
            integer(iintegers),intent(in) :: isrc
            integer(iintegers),parameter :: l=1
            integer(iintegers),parameter :: k=5
            !               integer,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9
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

              class default
                stop 'coeff_symmetry : unexpected type for OPP !'
            end select

            if(ldebug_optprop) then
              if(real(sum(coeff_symmetry)).gt.one) then
                print *,'sum of diffuse coeff_symmetrys bigger one!',sum(coeff_symmetry),'for src=',isrc,'coeff_symmetry:',coeff_symmetry
                call exit()
              endif
            endif
        end function

end module
