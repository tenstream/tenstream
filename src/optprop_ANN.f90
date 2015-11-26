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

module m_optprop_ANN
  USE m_data_parameters, ONLY : ireals, iintegers, zero,one,i1, mpiint
  use m_optprop_parameters, only: ldebug_optprop, lut_basename, &
      Ndz_8_10,Nkabs_8_10,Nksca_8_10,Ng_8_10,Nphi_8_10,Ntheta_8_10,Ndir_8_10,Ndiff_8_10, &
      ldelta_scale,delta_scale_truncate
  use m_netcdfio
  use mpi
  use m_helper_functions, only : imp_bcast

  implicit none
  private
  public ANN_init,ANN_get_dir2dir,ANN_get_dir2diff,ANN_get_diff2diff

  logical,parameter :: check_input=.True.

  integer(mpiint) :: myid,comm_size,mpierr

  type ANN
    real(ireals),allocatable,dimension(:) :: weights, units
    integer(iintegers),allocatable,dimension(:) :: inno, outno
    real(ireals),allocatable,dimension(:,:) :: eni, deo, inlimits
    integer(iintegers),allocatable,dimension(:,:) :: conec
    real(ireals),allocatable,dimension(:) :: lastcall,lastresult
    integer(iintegers) :: in_size=-1, out_size=-1
    logical :: initialized=.False.
  end type

  type(ANN),save :: diff2diff_network, dir2dir_network, dir2diff_network, direct_network

  real(ireals),parameter :: min_lim_coeff = zero
  logical,parameter :: lrenormalize=.True.
! logical,parameter :: lrenormalize=.False.

contains

  subroutine ANN_init(dx,dy,comm)
      real(ireals),intent(in) :: dx,dy
      integer(iintegers) :: idx,idy
      character(len=300) :: basename, netname, descr
      integer(iintegers) :: ierr

      integer(mpiint), intent(in) :: comm
      
      integer(iintegers),parameter :: horiz_rounding=1 ! round LUT for various horizontal distances: e.g. horiz_rounding=10 -> dx=66.7 ==> dx=70
!      integer(iintegers) :: phi,theta,iphi,itheta

      call MPI_Comm_rank(comm, myid, mpierr)
      call MPI_Comm_size(comm, comm_size, mpierr)

      if (myid.eq.0) then
        if(dx.ne.dy) then
          print *,'dx ne dy,  we probably dont have a network for asymmetric grid sizes!.... exiting....'
          call exit()
        endif

        idx = nint( dx/horiz_rounding  ) * horiz_rounding
        idy = nint( dy/horiz_rounding  ) * horiz_rounding

        basename = trim(lut_basename)//'_dstorder_8_10.'

        write(descr,FMT='("diffuse.dx",I0,".pspace.dz",I0,".kabs",I0,".ksca",I0,".g",I0,".delta_",L1,"_",F0.3)') &
          idx,Ndz_8_10,Nkabs_8_10,Nksca_8_10,Ng_8_10,ldelta_scale,delta_scale_truncate

        netname = trim(basename)//trim(descr)//'_diff2diff.ANN.nc'
        call loadnet(netname, diff2diff_network, ierr)

        write(descr,FMT='("direct.dx",I0,".pspace.dz",I0,".kabs",I0,".ksca",I0,".g",I0,".phi",I0,".theta",I0,".delta_",L1,"_",F0.3)') &
          idx,Ndz_8_10,Nkabs_8_10,Nksca_8_10,Ng_8_10,Nphi_8_10,Ntheta_8_10,ldelta_scale,delta_scale_truncate

        netname = trim(basename)//trim(descr)//'_dir2diff.ANN.nc'
        call loadnet(netname, dir2diff_network, ierr)

        netname = trim(basename)//trim(descr)//'_dir2dir.ANN.nc'
        call loadnet(netname, dir2dir_network, ierr)
      endif

      if (comm_size.gt.1) then 
        call scatter_ANN ( diff2diff_network )
        call scatter_ANN ( dir2diff_network  )
        call scatter_ANN ( dir2dir_network   )
      endif

  end subroutine

  subroutine scatter_ANN ( net )
      type(ANN) :: net

      call imp_bcast(net%weights    , 0_mpiint, myid )
      call imp_bcast(net%units      , 0_mpiint, myid )
      call imp_bcast(net%inno       , 0_mpiint, myid )
      call imp_bcast(net%outno      , 0_mpiint, myid )
      call imp_bcast(net%conec      , 0_mpiint, myid )
      call imp_bcast(net%deo        , 0_mpiint, myid )
      call imp_bcast(net%eni        , 0_mpiint, myid )
      call imp_bcast(net%inlimits   , 0_mpiint, myid )
      call imp_bcast(net%lastcall   , 0_mpiint, myid )
      call imp_bcast(net%lastresult , 0_mpiint, myid )
      call imp_bcast(net%in_size    , 0_mpiint, myid )
      call imp_bcast(net%out_size   , 0_mpiint, myid )
      call imp_bcast(net%initialized, 0_mpiint, myid )
  end subroutine

  subroutine loadnet(netname,net,ierr)
      character(300) :: netname
      type(ANN) :: net
      integer(iintegers),intent(out) :: ierr
      integer(iintegers) :: errcnt,k

      errcnt=0
      if(.not.allocated(net%weights)) then
        call ncload([netname,'weights' ],net%weights ,ierr) ; errcnt = ierr        !  ; print *,'loading weights ',ierr
        call ncload([netname,'units'   ],net%units   ,ierr) ; errcnt = errcnt+ierr !  ; print *,'loading units   ',ierr
        call ncload([netname,'inno'    ],net%inno    ,ierr) ; errcnt = errcnt+ierr !  ; print *,'loading inno    ',ierr
        call ncload([netname,'outno'   ],net%outno   ,ierr) ; errcnt = errcnt+ierr !  ; print *,'loading outno   ',ierr
        call ncload([netname,'conec'   ],net%conec   ,ierr) ; errcnt = errcnt+ierr !  ; print *,'loading conec   ',ierr
        call ncload([netname,'deo'     ],net%deo     ,ierr) ; errcnt = errcnt+ierr !  ; print *,'loading deo     ',ierr
        call ncload([netname,'eni'     ],net%eni     ,ierr) ; errcnt = errcnt+ierr !  ; print *,'loading eni     ',ierr
        call ncload([netname,'inlimits'],net%inlimits,ierr) ; errcnt = errcnt+ierr !  ; print *,'loading inlimits',ierr
        
        if(ldebug_optprop) &
          print *,'Loading ANN from: ',trim(netname),' resulted in errcnt',errcnt
        if(errcnt.ne.0) then
          ierr = errcnt
          print *,'ERROR loading ANN for netname :: ',trim(netname),' ::: ',ierr
          stop 'Could not load ANN'
          return
        endif

        net%in_size = size(net%inno)
        net%out_size= size(net%outno)

        if(ldebug_optprop) then
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

  subroutine ANN_get_dir2dir(dz, kabs, ksca, g, phi, theta, C)
      real(ireals),intent(in) :: dz,kabs,ksca,g,phi,theta
      real(ireals),intent(out) :: C(:)
      real(ireals) :: C2(dir2diff_network%out_size)

      integer(iintegers) :: ierr,isrc
      real(ireals) :: norm

      call calc_net(C, [dz,kabs,ksca,g,phi,theta] , dir2dir_network,ierr )
      if(ierr.ne.0) then
        print *,'Error when calculating dir2dir_net coeffs',ierr
        call exit()
      endif

      C = min(one, max(C,zero) )
      where( C.le.min_lim_coeff )
        C = zero
      endwhere

      !Check for energy conservation:
      if(lrenormalize) then
        do isrc=1,Ndir_8_10
          norm = sum( C( isrc:size(C):Ndir_8_10 ) )
          if(real(norm).gt.one) then
            C( isrc:size(C):Ndir_8_10 ) = C( isrc:size(C):Ndir_8_10 )/norm 
            !          print *,'dir2dir renormalization:',norm,' ::: ',sum( C( isrc:size(C):Ndir_8_10 ) )
          endif
        enddo
      endif

!     if(lrenormalize) then                                                               
!       call calc_net(C2, [dz,kabs,ksca,g,phi,theta], dir2diff_network,ierr)               
!       do isrc=1,Ndir_8_10                                                                
!         norm = sum( C(isrc:size(C):Ndir_8_10) ) + sum( C2(isrc:size(C2):Ndiff_8_10) )    
!         if(real(norm).gt.one) then                                                       
!           C(isrc:size(C):Ndir_8_10)=C( isrc:size(C):Ndir_8_10 )/norm                     
!         endif                                                                            
!       enddo                                                                              
!     endif                                                                                

   end subroutine

  subroutine ANN_get_dir2diff(dz, kabs, ksca, g, phi, theta, C)
      real(ireals),intent(in) :: dz,kabs,ksca,g,phi,theta
      real(ireals),intent(out) :: C(:)
      real(ireals) :: C2(dir2dir_network%out_size)

      integer(iintegers) :: ierr,isrc
      real(ireals) :: norm

      call calc_net(C, [dz,kabs,ksca,g,phi,theta] , dir2diff_network,ierr )
      C = C/1000.0
      if(ierr.ne.0) then
        print *,'Error when calculating dir2diff_net coeffs',ierr
        call exit()
      endif

      C = min(one, max(C,zero) )
      where( C.le.min_lim_coeff )
        C = zero
      endwhere

      !Check for energy conservation:
      if(lrenormalize) then
        do isrc=1,Ndiff_8_10
          norm = sum( C( isrc:size(C):Ndir_8_10 ) )
          if(real(norm).gt.one) then
            C( isrc:size(C):Ndir_8_10 ) = C( isrc:size(C):Ndir_8_10 )/ norm
            !          print *,'dir2diff renormalization:',norm,' ::: ',sum( C( isrc:size(C):Ndir_8_10 ) )
          endif
        enddo
      endif
!     if(lrenormalize) then                                                              
!       call calc_net(C2, [dz,kabs,ksca,g,phi,theta], dir2dir_network,ierr)              
!       do isrc=1,Ndir_8_10                                                              
!         norm = sum( C(isrc:size(C):Ndiff_8_10) ) + sum( C2(isrc:size(C2):Ndir_8_10) )  
!         if(real(norm).gt.one) then                                                     
!           C(isrc:size(C):Ndiff_8_10)=C( isrc:size(C):Ndiff_8_10 )/norm                 
!         endif                                                                          
!       enddo                                                                            
!     endif                                                                              
   end subroutine 

   subroutine ANN_get_diff2diff(dz, kabs, ksca, g, C)
      real(ireals),intent(out) :: C(:)
      real(ireals),intent(in) :: dz,kabs,ksca,g
      integer(iintegers) :: ierr,isrc
      real(ireals) :: norm

      if(.not.diff2diff_network%initialized) then
        print *,'network that is about to be used for coeffs is not loaded! diffuse:'
        call exit()
      endif

      call calc_net(C, [dz,kabs,ksca,g],diff2diff_network,ierr )
 !     C = C/1000.0
      if(ierr.ne.0) then
        print *,'Error when calculating diff_net coeffs',ierr
        call exit()
      endif

      C = min(one, max(C,zero) )
      where( C.le.min_lim_coeff )
        C = zero
      endwhere

      !Check for energy conservation:
      if(lrenormalize) then
        do isrc=1,Ndiff_8_10
          norm = sum( C( isrc:size(C):Ndiff_8_10 ) )
          if(real(norm).gt.one) then
            C( isrc:size(C):Ndiff_8_10 ) = C( isrc:size(C):Ndiff_8_10 )/ norm
            !          print *,'diffuse renormalization:',norm,' ::: ',sum( C( isrc:size(C):Ndiff_8_10 ) )
          endif
        enddo
      endif
    end subroutine
 
  subroutine calc_net(coeffs,inp,net,ierr)
      type(ANN),intent(inout) :: net
      real(ireals),intent(out):: coeffs(net%out_size )
      real(ireals),intent(in) :: inp(:)
      integer(iintegers),intent(out) :: ierr

      real(ireals) :: input(net%in_size)

      integer(iintegers) :: k,srcnode,trgnode,ctrg,xn


      ierr=0

      input = inp

      if(all(approx(net%lastcall,inp) )) then
        coeffs = net%lastresult
        return
      endif

      !        input([1,2]) = log(input([1,2]))
      !        input([4,5]) = log(input([4,5]))
      !        input(   7 ) = input(   7 ) *100._ireals

      ! Concerning the lower limits for optprops, just take the ANN limits. Theres not happening much anyway.
      input(2) = max(net%inlimits(2,1),input(2))
      input(3) = max(net%inlimits(3,1),input(3))
      input(4) = max(net%inlimits(4,1),input(4))

!      if(net%in_size.ge.5) then ! we should not fudge solar angles... this might confuse users...
!        input(5) = max(net%inlimits(5,1),input(5))
!        input(6) = max(net%inlimits(6,1),input(6))
!      endif

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
        if(inp(1).lt.net%inlimits(1,1).or.inp(1).gt.net%inlimits(1,2)) then
          print *,'dz out of ANN range',inp(1),'limits:',net%inlimits(1,:)        ;ierr=1;! call exit()
        endif
        if(input(2).lt.net%inlimits(2,1).or.input(2).gt.net%inlimits(2,2)) then
          print *,'kabs out of ANN range',input(2),'limits:',net%inlimits(2,:) ;ierr=2;! call exit()
        endif
        if(input(3).lt.net%inlimits(3,1).or.input(3).gt.net%inlimits(3,2)) then
          print *,'ksca out of ANN range',input(3),'limits:',net%inlimits(3,:)    ;ierr=3;! call exit()
        endif
        if(input(4).lt.net%inlimits(4,1).or.input(4).gt.net%inlimits(4,2)) then
          print *,'g out of ANN range',input(4),'limits:',net%inlimits(4,:) ;ierr=4;! call exit()
        endif
        if(net%in_size.ge.5) then
          if(input(5).lt.net%inlimits(5,1).or.input(5).gt.net%inlimits(5,2)) then
            print *,'phi out of ANN range',input(5),'limits:',net%inlimits(5,:) ;ierr=5;! call exit()
          endif
          if(input(6).lt.net%inlimits(6,1).or.input(6).gt.net%inlimits(6,2)) then
            print *,'theta out of ANN range',input(6),'limits:',net%inlimits(6,:)    ;ierr=6;! call exit()
          endif
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

end module
