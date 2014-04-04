module tenstream_optprop_ANN_8_10
      use data_parameters, only : ireals, iintegers, zero,one,i1
      use boxmc_parameters_8_10, only : dir_streams,diff_streams
      use arrayio

      implicit none
      private
      public ANN_init,ANN_get_dir2dir,ANN_get_dir2diff,ANN_get_diff2diff

      logical,parameter :: ldebug=.True.,check_input=.True.

      integer,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9

      type ANN
        real(ireals),allocatable,dimension(:) :: weights, units
        integer(iintegers),allocatable,dimension(:) :: inno, outno
        real(ireals),allocatable,dimension(:,:) :: eni, deo, inlimits
        integer(iintegers),allocatable,dimension(:,:) :: conec
        real(ireals),allocatable,dimension(:) :: lastcall,lastresult
        integer(iintegers) :: in_size=-1, out_size=-1
        logical :: initialized=.False.
      end type

      type(ANN) :: diff2diff_network, dir2dir_network, dir2diff_network, direct_network

      contains

        subroutine ANN_init(dx,dy)
            real(ireals),intent(in) :: dx,dy
            character(300),parameter :: fname='/usr/users/jakub/data/tenstream/ANN/LUT2ANN.h5'
            character(300) :: netname
            integer(iintegers) :: ierr,k

            integer(iintegers) :: phi,theta,iphi,itheta
            if(dx.ne.dy) then
              print *,'dx ne dy,  this is probably wrong!.... exiting....'
              call exit()
            endif

            write(netname, FMT='("net_",I0,"/dir2diff")' ) nint(dx/10)*10
            call loadnet(fname, netname, dir2diff_network, ierr)

            write(netname, FMT='("net_",I0,"/dir2dir")' ) nint(dx/10)*10
            call loadnet(fname, netname, dir2dir_network, ierr)

            write(netname, FMT='("net_I0/diff2diff")' ) nint(dx/10)*10
            call loadnet(fname, netname, diff2diff_network, ierr)

!            write(netname, FMT='("net_I0/direct")' ) nint(dx/10)*10
!            call loadnet(fname, netname, direct_network, ierr)
        end subroutine

        subroutine loadnet(fname, netname,net,ierr)
            character(300) :: fname, netname
            type(ANN) :: net
            integer(iintegers),intent(out) :: ierr
            integer(iintegers) :: errcnt,k
            ierr=0
            if(.not.allocated(net%weights)) then
              call h5load([fname,netname,'weights' ],net%weights ,ierr) ; errcnt = ierr         ! ; print *,'loading weights ',ierr
              call h5load([fname,netname,'units'   ],net%units   ,ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading units   ',ierr
              call h5load([fname,netname,'inno'    ],net%inno    ,ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading inno    ',ierr
              call h5load([fname,netname,'outno'   ],net%outno   ,ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading outno   ',ierr
              call h5load([fname,netname,'conec'   ],net%conec   ,ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading conec   ',ierr
              call h5load([fname,netname,'deo'     ],net%deo     ,ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading deo     ',ierr
              call h5load([fname,netname,'eni'     ],net%eni     ,ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading eni     ',ierr
              call h5load([fname,netname,'inlimits'],net%inlimits,ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading inlimits',ierr
              if(ldebug) print *,'Loading ANN from:',fname,'name of Network: ',trim(netname),' resulted in errcnt',errcnt
              if(errcnt.ne.0) return

              net%in_size = size(net%inno)
              net%out_size= size(net%outno)

              if(ldebug) then
                print *,'shape eni',shape(net%eni),'weights',shape(net%weights),'conec',shape(net%conec)
                do k=1,ubound(net%inlimits,1)
                  print *,'input limits(',k,')',net%inlimits(k,:),'eni',net%eni(k,:)
                enddo
              endif

              allocate(net%lastcall(net%in_size) ) ; net%lastcall=-1
              allocate(net%lastresult(net%out_size) )

              net%initialized = .True.
            endif

        end subroutine

        !{{{ calc dir_coeffs
        subroutine ANN_get_dir2dir(dz, od, w, g, phi, theta, C)
            real(ireals),intent(in) :: dz,od,w,g,phi,theta
            real(ireals),intent(out) :: C(dir_streams**2)

            integer(iintegers) :: ierr,k,istart,iend
            real(ireals) :: coeff_sum

            call calc_net(C, [dz,od,w,g,phi,theta] , dir2dir_network,ierr )
            if(ierr.ne.0) then
              print *,'Error when calculating dir2dir_net coeffs',ierr
              call exit()
            endif

            ! make sure energy conservation is given
            C = min(one, max(C,zero) )

            do k=1,dir_streams
              istart = (k-i1)*dir_streams+i1
              iend   = k*dir_streams
              coeff_sum = sum( C(istart:iend) )
              if( coeff_sum.gt.one ) C(istart:iend) = C(istart:iend)/(coeff_sum+1e-8_ireals)
            enddo
        end subroutine
        subroutine ANN_get_dir2diff(dz, od, w, g, phi, theta, C)
            real(ireals),intent(in) :: dz,od,w,g,phi,theta
            real(ireals),intent(out) :: C(dir_streams*diff_streams)

            integer(iintegers) :: ierr,k,istart,iend
            real(ireals) :: coeff_sum

            call calc_net(C, [dz,od,w,g,phi,theta] , dir2diff_network,ierr )
            if(ierr.ne.0) then
              print *,'Error when calculating dir2diff_net coeffs',ierr
              call exit()
            endif

            ! make sure energy conservation is given
            C = min(one, max(C,zero) )

            do k=1,dir_streams
              istart = (k-i1)*diff_streams+i1
              iend   = k*diff_streams
              coeff_sum = sum( C(istart:iend) )
              if( coeff_sum.gt.one ) C(istart:iend) = C(istart:iend)/(coeff_sum+1e-8_ireals)
            enddo

        end subroutine
      !}}}

      !{{{ calc diff_coeffs
      subroutine ANN_get_diff2diff(dz, od, w, g, C)
        real(ireals),allocatable :: C(:)
        real(ireals),intent(in) :: dz,od,w,g
        integer(iintegers) :: ierr
        real(ireals) :: tmpsum

        if(.not.diff2diff_network%initialized) then
                print *,'network that is about to be used for coeffs is not loaded! diffuse:'
                call exit()
        endif
        allocate(C(diff2diff_network%out_size) )

        call calc_net(C, [dz,od,w,g],diff2diff_network,ierr )
        if(ierr.ne.0) then
                print *,'Error when calculating diff_net coeffs',ierr
                call exit()
        endif

        ! make sure energy conservation is given
        C = min(one, max(C,zero) )
      end subroutine
        !}}}

!{{{ calc_ANN output
subroutine calc_net(coeffs,inp,net,ierr)
        type(ANN),intent(inout) :: net
        real(ireals),intent(out):: coeffs(net%out_size )
        real(ireals),intent(in) :: inp(:)
        integer(iintegers),intent(out) :: ierr

        real(ireals) :: input(net%in_size),tmpsum

        integer(iintegers) :: k,srcnode,trgnode,ctrg,xn,phi,theta

        ierr=0

        input = inp

        if(all(approx(net%lastcall,inp) )) then
                coeffs = net%lastresult
                return
        endif
        net%lastcall = input

!        input([1,2]) = log(input([1,2]))
!        input([4,5]) = log(input([4,5]))
!        input(   7 ) = input(   7 ) *100._ireals

        ! Concerning the lower limits for optprops, just take the ANN limits. Theres not happening much anyway.
        input(1) = max(net%inlimits(1,1),input(1))
        input(2) = max(net%inlimits(2,1),input(2))
        input(4) = max(net%inlimits(4,1),input(4))
        input(5) = max(net%inlimits(5,1),input(5))

        ! Concerning the upper limits for optprops, take the ANN limits. This is huge TODO, as this might lead to super wrong results! This is a very dirty hack for the time being!
!        input(1) = min(net%inlimits(1,2),input(1))
!        input(2) = min(net%inlimits(2,2),input(2))
!        input(4) = min(net%inlimits(4,2),input(4))
!        input(5) = min(net%inlimits(5,2),input(5))

        ! for asymmetry parameter, this is odd, that the cosmo ANN input doesnt capture it.
        ! this should not happen for a real run but if this is a artificial case, it very well may be.
        ! TODO: make sure, this doesnt happen.
!        input(3) = max( net%inlimits(3,1), min( net%inlimits(3,2), input(3) ))
!        input(6) = max( net%inlimits(6,1), min( net%inlimits(6,2), input(6) ))

!{{{ check input
!       print *,'This is function: ANN_get_direct_Transmission'
        if(check_input) then
                if(inp(1).lt.net%inlimits(8,1).or.inp(1).gt.net%inlimits(8,2)) then
                        print *,'dz out of ANN range',inp(1),'limits:',net%inlimits(8,:)        ;ierr=8;! call exit()
                endif
                if(input(1).lt.net%inlimits(1,1).or.input(1).gt.net%inlimits(1,2)) then
                        print *,'bg kabs out of ANN range',input(1),'limits:',net%inlimits(1,:) ;ierr=1;! call exit()
                endif
                if(input(2).lt.net%inlimits(2,1).or.input(2).gt.net%inlimits(2,2)) then
                        print *,'bg ksca out of ANN range',input(2),'limits:',net%inlimits(2,:) ;ierr=2;! call exit()
                endif
                if(input(3).lt.net%inlimits(3,1).or.input(3).gt.net%inlimits(3,2)) then
                        print *,'bg g out of ANN range',input(3),'limits:',net%inlimits(3,:)    ;ierr=3;! call exit()
                endif
                if(input(4).lt.net%inlimits(4,1).or.input(4).gt.net%inlimits(4,2)) then
                        print *,'fg kabs out of ANN range',input(4),'limits:',net%inlimits(4,:) ;ierr=4;! call exit()
                endif
                if(input(5).lt.net%inlimits(5,1).or.input(5).gt.net%inlimits(5,2)) then
                        print *,'fg ksca out of ANN range',input(5),'limits:',net%inlimits(5,:) ;ierr=5;! call exit()
                endif
                if(input(6).lt.net%inlimits(6,1).or.input(6).gt.net%inlimits(6,2)) then
                        print *,'fg g out of ANN range',input(6),'limits:',net%inlimits(6,:)    ;ierr=6;! call exit()
                endif
                if(input(7).lt.net%inlimits(7,1).or.input(7).gt.net%inlimits(7,2)) then
                        print *,'fg frac out of ANN range',input(7),'limits:',net%inlimits(7,:) ;ierr=7;! call exit()
                endif
        endif
!}}}

        ! normalize input
        do k=1,ubound(net%inno,1)
          net%units(net%inno(k)) = net%eni(k,1) * input(k) + net%eni(k,2)
        enddo

        ! Propagate signal through conecs
        ctrg = net%conec(1,2)
        net%units(ctrg) = 0._ireals
        do xn=1,ubound(net%weights,1)
          srcnode = net%conec(xn,1)
          trgnode = net%conec(xn,2)
          !if next unit
          if (trgnode.ne.ctrg) then
             net%units(ctrg) = 1._ireals/(1._ireals+exp(-net%units(ctrg)))
             ctrg = trgnode
             if (srcnode.eq.0) then !handle bias
                 net%units(ctrg) = net%weights(xn)
             else
                 net%units(ctrg) = net%units(srcnode) * net%weights(xn)
             endif
          else
             if (srcnode.eq.0) then !handle bias
                net%units(ctrg) = net%units(ctrg) + net%weights(xn)
             else
                net%units(ctrg) = net%units(ctrg) + net%units(srcnode) * net%weights(xn)
             endif
          endif
        enddo
        net%units(ctrg) = 1._ireals/(1._ireals+exp(-net%units(ctrg))) !for last unit

        ! Denormalize output
        do k=1,ubound(net%outno,1)
          coeffs(k) = net%deo(k,1) * net%units(net%outno(k)) + net%deo(k,2)
        enddo

        net%lastcall = input
        net%lastresult = coeffs(:)
        contains 
          pure elemental logical function approx(a,b)
            real(ireals),intent(in) :: a,b
            real(ireals),parameter :: eps=1e-5
            if( a.le.b+eps .and. a.ge.b-eps ) then
              approx = .True.
            else
              approx = .False.
            endif
          end function

      end subroutine
      !}}}
!
!      !{{{ get dir2dir
!      subroutine ANN_get_dir2dir(src,netout,coeff)
!        real(ireals),intent(out) :: coeff(3)
!        real(ireals),intent(in) :: netout(:)
!        integer(iintegers),intent(in) :: src
!
!        select case (src)
!                case(1)
!                        coeff = netout(1:3)
!                case(2)
!                        coeff = netout(7:9)
!                case(3)
!                        coeff = netout(4:6)
!                case default
!                        print *,'cant call ANN_get_diffuse with src being',src,'!'
!                        call exit()
!        end select
!        if(sum(coeff).ge.1._ireals) then
!                print *,'sum of dir2dir coeffs bigger one!',sum(coeff),'for src=',src,'coeff:',coeff
!                call exit()
!        endif
!      end subroutine
!      !}}}
!      !{{{ get_dir2diff
!      subroutine ANN_get_dir2diff(src,netout,coeff)
!        real(ireals),intent(out) :: coeff(10)
!        real(ireals),intent(in) :: netout(:)
!        integer(iintegers),intent(in) :: src
!
!        real(ireals) :: dir(3),norm
!
!        select case (src)
!                case(1)
!                        coeff = netout(10:19)
!                case(2)
!                        coeff = netout(20:29)
!                case(3)
!                        coeff = netout(30:39)
!                case default
!                        print *,'cant call ANN_get_diffuse with src being',src,'!'
!                        call exit()
!        end select
!
!        ! turn coeffs by 90 deg. from boxmc coords to cosmo coords, where stream 2 for example is on x axis and now in cosmo is on north axes.
!        coeff = coeff([1,2,7,8,9,10,3,4,5,6])
!
!        call ANN_get_dir2dir(src,netout,dir)
!        norm = max(0._ireals, min( sum(coeff),1._ireals-sum(dir)-1e-8_ireals))
!        coeff = coeff/max(sum(coeff),tiny(coeff)) * norm
!        if(sum(coeff).gt.1._ireals-sum(dir) .or. sum(coeff).lt.0._ireals ) then
!                print *,'sum of dir2diff coeffs bigger one-direct!',sum(coeff),'for src=',src,'coeff:',coeff,'direct',dir,'sum(dir)',sum(dir),'norm',norm,'==>',coeff/sum(coeff) * norm
!                call exit()
!        endif
!      end subroutine
!      !}}}
!
!!{{{ diff2diff
!      subroutine ANN_get_diffuse(src,netout,coeff)
!        real(ireals),intent(out) :: coeff(10)
!        real(ireals) :: netout(:)
!        integer(iintegers),intent(in) :: src
!        integer(iintegers),parameter :: l=1 ! this is the first coeff for diffuse radiation.(depends on the partitioning of target data... if unclear: ask Fabian)
!        integer(iintegers),parameter :: k=5 ! this is the scnd part of coeffs for diffuse radiation.(depends on the partitioning of target data... if unclear: ask Fabian)
!        ! the idea is: if a photon walks into the cube and wants to go through each stream one after another, what is the symmetry to either the first(l::0) stream or the second(k::2) stream
!        if( k+7.ne.size(netout)) then
!                print *,'this is wondrous, I thought, the last diffuse coeff should be the end of the ANN output',k+7,'ne',size(netout)
!                call exit()
!        endif
!
!        select case (src)
!                case(E_up)
!                        coeff = netout([l+1, l+0, l+3, l+3, l+2, l+2, l+3, l+3, l+2, l+2])
!                case(E_dn)
!                        coeff = netout([l+0, l+1, l+2, l+2, l+3, l+3, l+2, l+2, l+3, l+3])
!                case(E_le_m)
!                        coeff = netout([k+0, k+1, k+3, k+2, k+5, k+4, k+7, k+7, k+6, k+6])
!                case(E_ri_m)
!                        coeff = netout([k+0, k+1, k+2, k+3, k+4, k+5, k+6, k+6, k+7, k+7])
!                case(E_le_p)
!                        coeff = netout([k+1, k+0, k+5, k+4, k+3, k+2, k+7, k+7, k+6, k+6])
!                case(E_ri_p)
!                        coeff = netout([k+1, k+0, k+4, k+5, k+2, k+3, k+6, k+6, k+7, k+7])
!                case(E_ba_m)
!                        coeff = netout([k+0, k+1, k+6, k+6, k+7, k+7, k+3, k+2, k+5, k+4])
!                case(E_fw_m)
!                        coeff = netout([k+0, k+1, k+6, k+6, k+7, k+7, k+2, k+3, k+4, k+5])
!                case(E_ba_p)
!                        coeff = netout([k+1, k+0, k+7, k+7, k+6, k+6, k+5, k+4, k+3, k+2])
!                case(E_fw_p)
!                        coeff = netout([k+1, k+0, k+7, k+7, k+6, k+6, k+4, k+5, k+2, k+3])
!                case default
!                        print *,'cant call ANN_get_diffuse with src being',src,'!'
!                        call exit()
!        end select
!        ! turn coeffs by 90 deg. from boxmc coords to cosmo coords, where stream 2 for example is on x axis and now in cosmo is on north axes.
!        coeff = coeff([1,2,7,8,9,10,3,4,5,6])
!        if(sum(coeff).ge.1._ireals) then
!                print *,'sum of diffuse coeffs bigger one!',sum(coeff),'for src=',src,'coeff:',coeff
!                call exit()
!        endif
!      end subroutine
!      !}}}

end module
