module tenstream_optprop_8_10
use data_parameters, only: ireals,iintegers,one,zero,i0,i1,mpiint
use boxmc_parameters_8_10, only : delta_scale,dir_streams,diff_streams

use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast
use tenstream_optprop_LUT_8_10
use tenstream_optprop_ANN_8_10

implicit none

private
public :: init_optprop, optprop_lookup_coeff,optprop_debug
integer(mpiint) :: ierr

logical :: optprop_debug=.False.

real(ireals) :: dx,dy
integer(iintegers),parameter :: coeff_mode=i0 ! 0 is LUT, 1 is Neural Network

contains

  subroutine init_optprop(dx_inp,dy_inp, azis, szas, comm)
      real(ireals),intent(in) :: szas(:),azis(:),dx_inp,dy_inp
      integer(mpiint) ,intent(in) :: comm

      dx=dx_inp
      dy=dy_inp

      select case (coeff_mode)
        case(i0)
          call init_LUT(dx,dy, azis, szas, comm)
        case(i1)
          call ANN_init(dx,dy)
        case default
          print *,'coeff mode optprop initialization not defined for ',coeff_mode
          call exit(-1)
      end select
      print *,'Optical properties initialized!'

  end subroutine
  pure elemental logical function approx(a,b)
    real(ireals),intent(in) :: a,b
    real(ireals),parameter :: eps=1e-5
    if( a.le.b+eps .and. a.ge.b-eps ) then
      approx = .True.
    else
      approx = .False.
    endif
  end function

  subroutine optprop_lookup_coeff(dz,kabs,ksca,g,dir,C,inp_angles)
        logical,intent(in) :: dir
        real(ireals),intent(in) :: dz,g,kabs,ksca
        real(ireals),intent(in),optional :: inp_angles(2)
        real(ireals),intent(out):: C(:)

        real(ireals),allocatable :: C_diff(:)

        real(ireals) :: S_diff(diff_streams),T_dir(dir_streams)
        logical,parameter :: determine_coeff_error=.True.
        real(ireals),parameter :: checking_limit=1e-5
        real(ireals) :: angles(2)
        integer(iintegers) :: isrc

        if(optprop_debug) then
          if( (any([dz,kabs,ksca,g].lt.zero)) .or. (any(isnan([dz,kabs,ksca,g]))) ) then
            print *,'optprop_lookup_coeff :: corrupt optical properties: bg:: ',[dz,kabs,ksca,g]
            call exit
          endif
        endif

        if(present(inp_angles)) then
          angles = inp_angles
          if(inp_angles(2).le.1e-3_ireals) angles(1) = zero ! if sza is close to 0, azimuth is symmetric -> dont need to distinguish
        endif

        if(coeff_mode.eq.i0) then
          if(present(inp_angles)) then ! obviously we want the direct coefficients
            if(dir) then ! specifically the dir2dir
              call LUT_get_dir2dir(dz,kabs,ksca,g,angles(1),angles(2),C)
            else ! dir2diff
              call LUT_get_dir2diff(dz,kabs,ksca,g,angles(1),angles(2),C)
            endif
          else
            ! diff2diff
            call LUT_get_diff2diff(dz,kabs,ksca,g,C_diff)
          endif

        else if( coeff_mode.eq.i1) then
          if(present(inp_angles)) then ! obviously we want the direct coefficients
            if(dir) then ! specifically the dir2dir
              call ANN_get_dir2dir(dz,kabs,ksca,g,angles(1),angles(2),C)
            else ! dir2diff
              call ANN_get_dir2diff(dz,kabs,ksca,g,angles(1),angles(2),C)
            endif
          else
            ! diff2diff
            call ANN_get_diff2diff(dz,kabs,ksca,g,C_diff)
          endif

        endif

        if(.not.present(inp_angles)) then
          do isrc=1,diff_streams
            C( (isrc-1)*diff_streams+1:isrc*diff_streams ) = coeff_symmetry(isrc, C_diff )
          enddo
          deallocate(C_diff)
        endif

! This enables on-line calculations of coefficients with bmc code. This takes FOREVER! - use this only to check if LUT is working correctly!
        if(determine_coeff_error) then 
                call random_number(T_dir(1)) 
                if(optprop_debug.and.(T_dir(1).le.checking_limit)) then
                        if(present(inp_angles)) then 
                                if(dir) then !dir2dir
                                        do isrc=1,dir_streams
                                                call bmc_wrapper(isrc,dx,dy,dz,kabs,ksca,g,.True.,delta_scale,angles(1),angles(2),-1,S_diff,T_dir)
                                                print "('check ',i1,' dir2dir ',6e10.2,' :: RMSE ',2e13.4,' coeff ',8f7.4,' bmc ',8f7.4)",&
                                                        isrc,dz,kabs,ksca,g,angles(1),angles(2),RMSE(C((isrc-1)*dir_streams+1:isrc*dir_streams), T_dir),C((isrc-1)*dir_streams+1:isrc*dir_streams),T_dir
                                        enddo
                                else ! dir2diff
                                        do isrc=1,dir_streams
                                                call bmc_wrapper(isrc,dx,dy,dz,kabs,ksca,g,.True.,delta_scale,angles(1),angles(2),-1,S_diff,T_dir)
                                                print "('check ',i1,' dir2diff ',6e10.2,' :: RMSE ',2e13.4,' coeff ',10e10.2)",&
                                                        isrc,dz,kabs,ksca,g,angles(1),angles(2),RMSE(C((isrc-1)*diff_streams+1:isrc*diff_streams), S_diff),abs(C((isrc-1)*diff_streams+1:isrc*diff_streams)-S_diff)
                                                print "(a10,e13.4,10e13.4)",'C_interp',sum(C((isrc-1)*diff_streams+1:isrc*diff_streams)),C((isrc-1)*diff_streams+1:isrc*diff_streams)
                                                print "(a10,e13.4,10e13.4)",'C_bmc   ',sum(S_diff),                  S_diff
                                                print *,''
                                        enddo
                                endif
                        else
                                ! diff2diff
                                do isrc=1,diff_streams
                                        call bmc_wrapper(isrc,dx,dy,dz,kabs,ksca,g,.False.,delta_scale,zero,zero,-1,S_diff,T_dir)
                                        print "('check ',i1,' diff2diff ',4e10.2,' :: RMSE ',2e13.4,' coeff err',10e10.2)",&
                                                isrc,dz,kabs,ksca,g,RMSE(C((isrc-1)*diff_streams+1:isrc*diff_streams), S_diff),abs(C((isrc-1)*diff_streams+1:isrc*diff_streams)-S_diff)
                                        print "(a10,e13.4,10e13.4)",'C_interp',sum(C((isrc-1)*diff_streams+1:isrc*diff_streams)),C((isrc-1)*diff_streams+1:isrc*diff_streams)
                                        print "(a10,e13.4,10e13.4)",'C_bmc   ',sum(S_diff),S_diff
                                        print *,''
                                enddo
                        endif ! angles_present
                endif ! is in checking_limit
        endif ! want to check
end subroutine
function RMSE(a,b)
        real(ireals) :: RMSE(2) ,a(:),b(size(a))
        RMSE = [sqrt( sum( (a-b)**2 )/size(a) ), sqrt( sum( (a-b)**2 )/size(a) ) / max( sum(b)/size(b),tiny(RMSE) ) ]
end function
        function coeff_symmetry(isrc,coeff)
               real(ireals) :: coeff_symmetry(diff_streams)
               real(ireals),intent(in) :: coeff(:)
               integer(iintegers),intent(in) :: isrc
               integer(iintegers),parameter :: l=1
               integer(iintegers),parameter :: k=5
!               integer,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9
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
                        print *,'cant call coeff_symmetry with src being',isrc,'!'
                        call exit()
                end select
                ! turn coeffs by 90 deg. from boxmc coords to cosmo coords, where stream 2 for example is on x axis and now in cosmo is on north axes.
                ! coeff_symmetry = coeff_symmetry([1,2,7,8,9,10,3,4,5,6])
                if(sum(coeff_symmetry).ge.1._ireals) then
                        print *,'sum of diffuse coeff_symmetrys bigger one!',sum(coeff_symmetry),'for src=',isrc,'coeff_symmetry:',coeff_symmetry
                        call exit()
                endif
        end function

end module
