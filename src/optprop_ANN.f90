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
  USE m_data_parameters, ONLY : ireals, iintegers, zero,one,i1, mpiint, default_str_len
  use m_optprop_parameters, only: ldebug_optprop, lut_basename, &
      Ntau, Nw0, Ng, Ndir_8_10, Ndiff_8_10, &
      ldelta_scale, delta_scale_truncate
  use m_netcdfio
  use mpi
  use m_helper_functions, only : imp_bcast, search_sorted_bisection


  implicit none
  private
  public ANN_init, ANN_destroy, ANN_get_dir2dir, ANN_get_dir2diff, ANN_get_diff2diff

  logical,parameter :: check_input=.True.

  integer(mpiint) :: myid,comm_size,mpierr

  type ANN
    real(ireals),allocatable,dimension(:) :: weights, units, tau, w0, g, phi, theta
    integer(iintegers),allocatable,dimension(:) :: inno, outno
    real(ireals),allocatable,dimension(:,:) :: eni, deo, inlimits
    integer(iintegers),allocatable,dimension(:,:) :: conec
    real(ireals),allocatable,dimension(:) :: lastcall,lastresult
    integer(iintegers) :: in_size=-1, out_size=-1
    logical :: initialized=.False.
  end type

  type(ANN),allocatable,save :: diff2diff_network, dir2dir_network, dir2diff_network

  real(ireals),parameter :: min_lim_coeff = 0.0001
! logical,parameter :: lrenormalize=.True.
  logical,parameter :: lrenormalize=.False.

contains
  subroutine ANN_destroy()
    if(allocated(diff2diff_network)) deallocate(diff2diff_network)
    if(allocated(dir2diff_network)) deallocate(dir2diff_network)
    if(allocated(dir2dir_network)) deallocate(dir2dir_network)
  end subroutine

  subroutine ANN_init( comm, ierr)
      character(len=default_str_len) :: basename, netname, descr
      integer(mpiint) :: ierr

      integer(mpiint), intent(in) :: comm

      call MPI_Comm_rank(comm, myid, mpierr)
      call MPI_Comm_size(comm, comm_size, mpierr)

      if (myid.eq.0) then

        basename = trim(lut_basename)//'_dstorder_8_10.'

        write(descr,FMT='("diffuse.tau",I0,".w0",I0,".g",I0,".delta_",L1,"_",F0.3)') &
          Ntau, Nw0, Ng, ldelta_scale, delta_scale_truncate

        allocate(diff2diff_network)
        netname = trim(basename)//trim(descr)//'_diff2diff.ANN.nc'
        call loadnet(netname, diff2diff_network, ierr)

        write(descr,FMT='("direct.tau",I0,".w0",I0,".g",I0,".delta_",L1,"_",F0.3)') &
          Ntau, Nw0, Ng, ldelta_scale, delta_scale_truncate

        allocate(dir2diff_network)
        netname = trim(basename)//trim(descr)//'_dir2diff.ANN.nc'
        call loadnet(netname, dir2diff_network, ierr)

        allocate(dir2dir_network)
        netname = trim(basename)//trim(descr)//'_dir2dir.ANN.nc'
        call loadnet(netname, dir2dir_network, ierr)
      endif
      call imp_bcast(comm, ierr, 0_mpiint)

      if (comm_size.gt.1 .and. ierr.eq.0) then
        if(myid.ne.0) then
          allocate( diff2diff_network )
          allocate( dir2diff_network  )
          allocate( dir2dir_network   )
        endif
        call scatter_ANN (comm, diff2diff_network)
        call scatter_ANN (comm, dir2diff_network)
        call scatter_ANN (comm, dir2dir_network)
      endif

  end subroutine

  subroutine scatter_ANN (comm, net)
    type(ANN) :: net
    integer(mpiint), intent(in) :: comm
    logical :: l_have_angles

    call imp_bcast(comm, net%weights    , 0_mpiint)
    call imp_bcast(comm, net%units      , 0_mpiint)
    call imp_bcast(comm, net%inno       , 0_mpiint)
    call imp_bcast(comm, net%outno      , 0_mpiint)
    call imp_bcast(comm, net%conec      , 0_mpiint)
    call imp_bcast(comm, net%deo        , 0_mpiint)
    call imp_bcast(comm, net%eni        , 0_mpiint)
    call imp_bcast(comm, net%inlimits   , 0_mpiint)
    call imp_bcast(comm, net%lastcall   , 0_mpiint)
    call imp_bcast(comm, net%lastresult , 0_mpiint)
    call imp_bcast(comm, net%in_size    , 0_mpiint)
    call imp_bcast(comm, net%out_size   , 0_mpiint)
    call imp_bcast(comm, net%initialized, 0_mpiint)
    call imp_bcast(comm, net%tau        , 0_mpiint)
    call imp_bcast(comm, net%w0         , 0_mpiint)
    call imp_bcast(comm, net%g          , 0_mpiint)
    if (myid.eq.0) then
      l_have_angles = allocated(net%phi) .and. allocated(net%theta)
    endif
    call imp_bcast(comm, l_have_angles  , 0_mpiint)

    if (l_have_angles) then
      call imp_bcast(comm, net%phi        , 0_mpiint)
      call imp_bcast(comm, net%theta      , 0_mpiint)
    endif
  end subroutine

  subroutine loadnet(netname,net,ierr)
      character(default_str_len) :: netname
      type(ANN) :: net
      integer(mpiint),intent(out) :: ierr
      integer(mpiint) :: errcnt,k

      errcnt=0
      if(.not.allocated(net%weights)) then
        call ncload([netname,'weights' ],net%weights , ierr) ; errcnt = ierr        ! ; print *,'loading weights ',ierr
        call ncload([netname,'units'   ],net%units   , ierr) ; errcnt = errcnt+ierr ! ; print *,'loading units   ',ierr
        call ncload([netname,'inno'    ],net%inno    , ierr) ; errcnt = errcnt+ierr ! ; print *,'loading inno    ',ierr
        call ncload([netname,'outno'   ],net%outno   , ierr) ; errcnt = errcnt+ierr ! ; print *,'loading outno   ',ierr
        call ncload([netname,'conec'   ],net%conec   , ierr) ; errcnt = errcnt+ierr ! ; print *,'loading conec   ',ierr
        call ncload([netname,'deo'     ],net%deo     , ierr) ; errcnt = errcnt+ierr ! ; print *,'loading deo     ',ierr
        call ncload([netname,'eni'     ],net%eni     , ierr) ; errcnt = errcnt+ierr ! ; print *,'loading eni     ',ierr
        call ncload([netname,'inlimits'],net%inlimits, ierr) ; errcnt = errcnt+ierr ! ; print *,'loading inlimits     ',ierr

        call ncload([netname,'pspace.tau'  ],net%tau  , ierr) ; errcnt = errcnt+ierr   !  ; print *,'loading inlimits',ierr
        call ncload([netname,'pspace.w0'   ],net%w0   , ierr) ; errcnt = errcnt+ierr   ! ; print *,'loading eni     ',ierr
        call ncload([netname,'pspace.g'    ],net%g    , ierr) ; errcnt = errcnt+ierr  ! ; print *,'loading inlimits',ierr
        call ncload([netname,'pspace.phi'  ],net%phi  , ierr) ; print *,'loading phi  ',ierr, allocated(net%phi  ), net%inlimits
        call ncload([netname,'pspace.theta'],net%theta, ierr) ; print *,'loading theta',ierr, allocated(net%theta), net%inlimits

        if(ldebug_optprop) &
          print *,'Loading ANN from: ',trim(netname),' resulted in errcnt',errcnt
        if(errcnt.ne.0) then
          ierr = errcnt
          print *,'ERROR loading ANN for netname :: ',trim(netname),' ::: ',ierr
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


  subroutine ANN_get_dir2dir(taux, tauz, w0, g, phi, theta, C)
      real(ireals),intent(in) :: taux, tauz, w0, g, phi, theta
      real(ireals) :: ind_taux, ind_tauz, ind_w0, ind_g, ind_phi, ind_theta
      real(ireals),intent(out) :: C(:)
      real(ireals) :: C2(dir2diff_network%out_size)

      integer(iintegers) :: ierr,isrc
      real(ireals) :: norm

      ind_taux  = search_sorted_bisection(dir2dir_network%tau  , taux)
      ind_tauz  = search_sorted_bisection(dir2dir_network%tau  , tauz)
      ind_w0    = search_sorted_bisection(dir2dir_network%w0   , w0  )
      ind_g     = search_sorted_bisection(dir2dir_network%g    ,g    )
      ind_phi   = search_sorted_bisection(dir2dir_network%phi  ,phi  )
      ind_theta = search_sorted_bisection(dir2dir_network%theta,theta)

      call calc_net(C, [ind_taux, ind_tauz, ind_w0, ind_g, ind_phi, ind_theta] , dir2dir_network,ierr )
      if(ierr.ne.0) then
        print *,'Error when calculating dir2dir_net coeffs',ierr
        call exit()
      endif

      C = min(one, max(C,zero) )
      where( C.le.min_lim_coeff )
        C = zero
      endwhere

      !Check for energy conservation:
!      if(lrenormalize) then
!        do isrc=1,Ndir_8_10
!          norm = sum( C( isrc:size(C):Ndir_8_10 ) )
!          if(real(norm).gt.one) then
!            C( isrc:size(C):Ndir_8_10 ) = C( isrc:size(C):Ndir_8_10 )/norm
!            !          print *,'dir2dir renormalization:',norm,' ::: ',sum( C( isrc:size(C):Ndir_8_10 ) )
!          endif
!        enddo
!      endif

      if(lrenormalize) then
        call calc_net(C2, [ind_taux, ind_tauz, ind_w0, ind_g, ind_phi, ind_theta], dir2diff_network,ierr)
        do isrc=1,Ndir_8_10
          norm = sum( C(isrc:size(C):Ndir_8_10) ) + sum( C2(isrc:size(C2):Ndiff_8_10) )
          if(real(norm).gt.one) then
            C(isrc:size(C):Ndir_8_10)=C( isrc:size(C):Ndir_8_10 )/norm
          endif
        enddo
      endif

   end subroutine

  subroutine ANN_get_dir2diff(taux, tauz, w0, g, phi, theta, C)
      real(ireals),intent(in) :: taux, tauz, w0, g,phi,theta
      real(ireals) :: ind_taux, ind_tauz, ind_w0, ind_g, ind_phi, ind_theta
      real(ireals),intent(out) :: C(:)
      real(ireals) :: C2(dir2dir_network%out_size)

      integer(mpiint) :: ierr,isrc
      real(ireals) :: norm

      ind_taux  = search_sorted_bisection(dir2diff_network%tau  , taux )
      ind_tauz  = search_sorted_bisection(dir2diff_network%tau  , tauz )
      ind_w0    = search_sorted_bisection(dir2diff_network%w0   , w0   )
      ind_g     = search_sorted_bisection(dir2diff_network%g    , g    )
      ind_phi   = search_sorted_bisection(dir2diff_network%phi  , phi  )
      ind_theta = search_sorted_bisection(dir2diff_network%theta, theta)

      call calc_net(C, [ind_taux, ind_tauz, ind_w0, ind_g, ind_phi, ind_theta] , dir2diff_network,ierr )
!      C = C/1000.0
      if(ierr.ne.0) then
        print *,'Error when calculating dir2diff_net coeffs',ierr
        call exit()
      endif

      C = min(one, max(C,zero) )
      where( C.le.min_lim_coeff )
        C = zero
      endwhere

      !Check for energy conservation:
!     if(lrenormalize) then
!       do isrc=1,Ndiff_8_10
!         norm = sum( C( isrc:size(C):Ndir_8_10 ) )
!         if(real(norm).gt.one) then
!           C( isrc:size(C):Ndir_8_10 ) = C( isrc:size(C):Ndir_8_10 )/ norm
!           !          print *,'dir2diff renormalization:',norm,' ::: ',sum( C( isrc:size(C):Ndir_8_10 ) )
!         endif
!       enddo
!     endif
      if(lrenormalize) then
        call calc_net(C2, [ind_taux, ind_tauz, ind_w0 ,ind_g, ind_phi, ind_theta], dir2dir_network, ierr)
        do isrc=1,Ndir_8_10
          norm = sum( C(isrc:size(C):Ndiff_8_10) ) + sum(C2(isrc:size(C2):Ndir_8_10))
          if(real(norm).gt.one) then
            C(isrc:size(C):Ndiff_8_10)=C( isrc:size(C):Ndiff_8_10 )/norm
          endif
        enddo
      endif
   end subroutine

   subroutine ANN_get_diff2diff(taux, tauz, w0, g, C)
      real(ireals),intent(out) :: C(:)
      real(ireals),intent(in) :: taux, tauz, w0, g
      real(ireals) :: ind_taux, ind_tauz, ind_w0, ind_g
      integer(mpiint) :: ierr,isrc
      real(ireals) :: norm

      if(.not.diff2diff_network%initialized) then
        print *,'network that is about to be used for coeffs is not loaded! diffuse:'
        call exit()
      endif

      ind_taux = search_sorted_bisection(diff2diff_network%tau , taux)
      ind_tauz = search_sorted_bisection(diff2diff_network%tau , tauz)
      ind_w0   = search_sorted_bisection(diff2diff_network%w0  , w0  )
      ind_g    = search_sorted_bisection(diff2diff_network%g   , g   )

      call calc_net(C, [ind_taux, ind_tauz, ind_w0, ind_g], diff2diff_network, ierr )
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
      integer(mpiint),intent(out) :: ierr

      real(ireals) :: input(net%in_size)

      integer(iintegers) :: k,srcnode,trgnode,ctrg,xn


      ierr=0

      input = inp-one

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
          print *,'taux out of ANN range',inp(1),'limits:',net%inlimits(1,:)        ;ierr=1;! call exit()
        endif
        if(input(2).lt.net%inlimits(2,1).or.input(2).gt.net%inlimits(2,2)) then
          print *,'tauz out of ANN range',input(2),'limits:',net%inlimits(2,:) ;ierr=2;! call exit()
        endif
        if(input(3).lt.net%inlimits(3,1).or.input(3).gt.net%inlimits(3,2)) then
          print *,'w0 out of ANN range',input(3),'limits:',net%inlimits(3,:)    ;ierr=3;! call exit()
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
