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

module m_optprop_LUT

  use mpi!, only: MPI_BCAST,MPI_LAND,MPI_LOR

  use m_helper_functions, only : approx,rel_approx,imp_bcast,mpi_logical_and,mpi_logical_or, search_sorted_bisection, CHKERR
  use m_data_parameters, only : ireals, iintegers,      &
    one, zero, i0, i1, i3, mpiint, nil, inil,           &
    imp_int, imp_real, imp_comm, imp_logical, numnodes, &
    default_str_len
  use m_optprop_parameters, only:          &
    ldebug_optprop, lut_basename,          &
    Naspect, Ntau, Nw0, Ng, Nphi, Ntheta,  &
    Ndir_1_2,Ndiff_1_2,interp_mode_1_2,    &
    Ndir_8_10,Ndiff_8_10,interp_mode_8_10, &
    ldelta_scale,delta_scale_truncate,     &
    stddev_atol, use_prescribed_LUT_dims,  &
    preset_aspect, preset_tau, preset_w0,  &
    preset_g, preset_theta
  use m_boxmc, only: t_boxmc,t_boxmc_8_10,t_boxmc_1_2
  use m_tenstream_interpolation, only: interp_4d
  use m_netcdfio

  implicit none

  private
  public :: t_optprop_LUT, t_optprop_LUT_8_10,t_optprop_LUT_1_2
  ! This module loads and generates the LUT-tables for Tenstream Radiation
  ! computations.
  ! It also holds functions for interpolation on the regular LUT grid.

  integer(iintegers) :: iierr
  integer(mpiint) :: myid,comm_size,mpierr

  type parameter_space
    real(ireals),allocatable ,dimension(:) :: aspect
    real(ireals),allocatable ,dimension(:) :: tau
    real(ireals),allocatable ,dimension(:) :: w0
    real(ireals),allocatable ,dimension(:) :: g
    real(ireals),allocatable ,dimension(:) :: phi
    real(ireals),allocatable ,dimension(:) :: theta
    real(ireals) , dimension(2)      :: range_aspect  = [ .01        , 5.0           ] ! is defined in set_parameter_space, so that min and max transmissions are met
    real(ireals) , dimension(2)      :: range_tau     = [ nil        , nil           ] ! is defined in set_parameter_space, so that min and max transmissions are met
    real(ireals) , dimension(2)      :: range_w0      = [ zero       , .999_ireals   ]
    real(ireals) , dimension(2)      :: range_g       = [ zero       , .999_ireals   ]
    real(ireals) , dimension(2)      :: range_phi     = [ zero       , 90._ireals    ]
    real(ireals) , dimension(2)      :: range_theta   = [ zero       , 90._ireals    ]
  end type

  type table
    real(ireals),allocatable :: c(:,:,:,:,:)
    real(ireals),allocatable :: stddev_tol(:,:,:,:,:)
    character(default_str_len),allocatable :: table_name_c(:)
    character(default_str_len),allocatable :: table_name_tol(:)
  end type

  type diffuseTable
    character(default_str_len) :: fname
    type(table) :: S
    type(parameter_space) :: pspace
  end type

  type directTable
    character(default_str_len),allocatable ::  fname(:,:) !dim=(phi,theta)
    type(table),allocatable :: S(:,:) !dim=(phi,theta)
    type(table),allocatable :: T(:,:) !dim=(phi,theta)
    type(parameter_space) :: pspace
  end type

  type,abstract :: t_optprop_LUT
    class(t_boxmc),allocatable :: bmc
    type(directTable),allocatable :: dirLUT
    type(diffuseTable),allocatable :: diffLUT

    integer(iintegers) :: Naspect, Ntau, Nw0, Ng, Nphi, Ntheta, interp_mode
    integer(iintegers) :: dir_streams=inil,diff_streams=inil
    logical :: LUT_initialized=.False.,optprop_LUT_debug=ldebug_optprop
    character(default_str_len) :: lutbasename

    contains
      procedure :: init
      procedure :: destroy
      procedure :: LUT_get_dir2dir
      procedure :: LUT_get_dir2diff
      procedure :: LUT_get_diff2diff
      procedure :: bmc_wrapper
      procedure :: scatter_LUTtables
      procedure :: createLUT_dir
      procedure :: createLUT_diff
      procedure :: loadLUT_dir
      procedure :: loadLUT_diff
      procedure :: set_parameter_space
  end type

  type,extends(t_optprop_LUT) :: t_optprop_LUT_1_2
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_8_10
  end type

contains

  subroutine destroy(OPP)
      class(t_optprop_LUT) :: OPP
      if(allocated(OPP%dirLUT )) deallocate(OPP%dirLUT )
      if(allocated(OPP%diffLUT)) deallocate(OPP%diffLUT)
      if(allocated(OPP%bmc    )) deallocate(OPP%bmc    )
      OPP%LUT_initialized=.False.
  end subroutine

  subroutine init(OPP, azis,szas, comm)
      class(t_optprop_LUT) :: OPP
      real(ireals),intent(in) :: szas(:),azis(:) ! all solar zenith angles that happen in this scene
      integer(mpiint) ,intent(in) :: comm

      character(default_str_len) :: descr

      call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
      call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

      if(.not.allocated(OPP%bmc)) then
        select type (OPP)
          class is (t_optprop_LUT_1_2)
            OPP%dir_streams  =  Ndir_1_2
            OPP%diff_streams =  Ndiff_1_2
            OPP%lutbasename=trim(lut_basename)//'_1_2.'
            allocate(t_boxmc_1_2::OPP%bmc)

          class is (t_optprop_LUT_8_10)
            OPP%dir_streams  = Ndir_8_10
            OPP%diff_streams = Ndiff_8_10
            OPP%lutbasename=trim(lut_basename)//'_8_10.'
            allocate(t_boxmc_8_10::OPP%bmc)

          class default
            stop 'initialize LUT: unexpected type for optprop_LUT object!'
        end select

        call OPP%bmc%init(comm)
      endif

      if(.not.allocated(OPP%diffLUT)) allocate(OPP%diffLUT)
      if(.not.allocated(OPP%dirLUT)) allocate(OPP%dirLUT)

      call OPP%set_parameter_space(OPP%diffLUT%pspace)
      call OPP%set_parameter_space(OPP%dirLUT%pspace )

      ! Name will be set in loadLUT_dir, different for each sun angle

      call OPP%loadLUT_dir(azis, szas, comm)

      ! Load diffuse LUT
      write(descr,FMT='("diffuse.aspect",I0,".tau",I0,".w0",I0,".g",I0,".delta_",L1,"_",F0.3)') &
                         OPP%Naspect, OPP%Ntau, OPP%Nw0, OPP%Ng, ldelta_scale, delta_scale_truncate

      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Loading diffuse LUT from ',trim(descr)
      OPP%diffLUT%fname = trim(OPP%lutbasename)//trim(descr)//'.nc'

      call OPP%loadLUT_diff(comm)

      ! Copy LUT to all ranks
      if(comm_size.gt.1) call OPP%scatter_LUTtables(comm, azis,szas)

      OPP%LUT_initialized=.True.
      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Done loading LUTs (shape diffLUT',shape(OPP%diffLUT%S%c),')'
  end subroutine

subroutine loadLUT_diff(OPP, comm)
    class(t_optprop_LUT) :: OPP
    integer(mpiint),intent(in) :: comm
    integer(iintegers) :: errcnt
    character(default_str_len) :: str(3)
    logical :: lstddev_inbounds

    associate ( LUT => OPP%diffLUT, S => OPP%diffLUT%S )

    if(allocated(S%c)) then
      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'diffLUT already loaded! is this a second call?'
      return
    endif


    write(str(1),FMT='(A)') 'diffuse'

    if(.not.allocated(S%table_name_c) ) allocate(S%table_name_c(3))
    write(str(2),FMT='(A)') "S" ; S%table_name_c = [LUT%fname,str(1),str(2)]

    if(.not.allocated(S%table_name_tol) ) allocate(S%table_name_tol(3))
    write(str(2),FMT='(A)') "S_tol"  ; S%table_name_tol = [LUT%fname,str(1),str(2)]

    errcnt=-1
    if(myid.eq.0) then
      lstddev_inbounds=.False.
      errcnt=0

      ! Load LUT's from file
      call ncload(S%table_name_c ,S%c,iierr) ; errcnt = errcnt+iierr

      ! check if all coefficients are in range between 0 and 1 and if we
      ! actually hold a lookuptable for the here specified parameter ranges.
      if(allocated(S%c) ) then
        if( any(S%c.gt.one) .or. any(S%c.lt.zero) ) errcnt=3
        call check_diffLUT_matches_pspace(LUT)
      endif

      call ncload(S%table_name_tol, S%stddev_tol,iierr) ; errcnt = errcnt+iierr

      if( allocated(S%stddev_tol) ) then
        lstddev_inbounds = all( S%stddev_tol.le.stddev_atol+epsilon(stddev_atol)*10 )

        if(.not.lstddev_inbounds) &
            print *,'Loading of diffuse tables S :: failed bc lstddev_inbounds',lstddev_inbounds,'::',maxval(S%stddev_tol),'(',stddev_atol+epsilon(stddev_atol)*10,')',errcnt

        if(OPP%optprop_LUT_debug .and. myid.eq.0) &
            print *,'... loading diffuse LUT',errcnt,lstddev_inbounds,'::',maxval(S%stddev_tol),'(',stddev_atol+epsilon(stddev_atol)*10,')'
      endif
    endif !master

    call mpi_bcast(errcnt           , 1_mpiint , imp_int     , 0_mpiint , comm , mpierr); call CHKERR(mpierr) ! inform all nodes if we were able to load the LUT
    call mpi_bcast(lstddev_inbounds , 1_mpiint , imp_logical , 0_mpiint , comm , mpierr); call CHKERR(mpierr) ! and if all coefficients are valid

    if(errcnt.ne.0 .or. .not.lstddev_inbounds ) then ! something is missing... lets try to recompute missing values in LUT
      if(myid.eq.0) then
        write(str(2),FMT='(A)') "pspace"

        if(OPP%optprop_LUT_debug) then
          print *,'Loading of diffuse tables failed for ',trim(LUT%fname),' :: ',trim(str(1)),'::',errcnt,'stddev required',lstddev_inbounds
          print *,''
          print *,'will write pspace with arguments:'
          print *,'0 ',trim(LUT%fname)
          print *,'1 ',trim(str(1))
          print *,'2 ',trim(str(2))
          print *,''
        endif

        write(str(3),FMT='(A)') "range_aspect"; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%range_aspect,iierr)
        write(str(3),FMT='(A)') "range_tau"   ; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%range_tau   ,iierr)
        write(str(3),FMT='(A)') "range_w0 "   ; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%range_w0    ,iierr)
        write(str(3),FMT='(A)') "range_g  "   ; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%range_g     ,iierr)

        write(str(3),FMT='(A)') "aspect   "   ; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%aspect,iierr)
        write(str(3),FMT='(A)') "tau      "   ; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%tau   ,iierr)
        write(str(3),FMT='(A)') "w0       "   ; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%w0    ,iierr)
        write(str(3),FMT='(A)') "g        "   ; call ncwrite([LUT%fname,str(1),str(2),str(3)],LUT%pspace%g     ,iierr)
      endif !master

      call OPP%createLUT_diff(LUT, comm)
    endif !error

    if(myid.eq.0) then
      deallocate(S%stddev_tol)
      if(OPP%optprop_LUT_debug) &
          print *,'Done loading diffuse LUTs',errcnt
    endif

    end associate
end subroutine
subroutine loadLUT_dir(OPP, azis,szas, comm)
    class(t_optprop_LUT) :: OPP
    real(ireals),intent(in) :: szas(:),azis(:) ! all solar zenith angles that happen in this scene
    integer(mpiint),intent(in) :: comm
    integer(iintegers) :: errcnt,iphi,itheta
    character(default_str_len) :: descr,str(5),varname(4)
    logical :: angle_mask(OPP%Nphi,OPP%Ntheta),lstddev_inbounds

    if(.not. allocated(OPP%dirLUT%fname) ) allocate( OPP%dirLUT%fname(OPP%Nphi,OPP%Ntheta) )
    if(.not. allocated(OPP%dirLUT%S    ) ) allocate( OPP%dirLUT%S    (OPP%Nphi,OPP%Ntheta) )
    if(.not. allocated(OPP%dirLUT%T    ) ) allocate( OPP%dirLUT%T    (OPP%Nphi,OPP%Ntheta) )

    write(str(1),FMT='(A)') 'direct'

    call determine_angles_to_load(comm, OPP%interp_mode, OPP%dirLUT, azis, szas, angle_mask)

    do itheta=1,OPP%Ntheta
      do iphi  =1,OPP%Nphi

        ! Set filename of LUT
        write(descr,FMT='("direct.aspect",I0,".tau",I0,".w0",I0,".g",I0,".phi",I0,".theta",I0,".delta_",L1,"_",F0.3)') &
            OPP%Naspect,OPP%Ntau, OPP%Nw0, OPP%Ng, &
            int(OPP%dirLUT%pspace%phi(iphi)),      &
            int(OPP%dirLUT%pspace%theta(itheta)),  &
            ldelta_scale,delta_scale_truncate

        OPP%dirLUT%fname(iphi,itheta) = trim(OPP%lutbasename)//trim(descr)//'.nc'

        associate ( phi   => int(OPP%dirLUT%pspace%phi(iphi)),      &
                    theta => int(OPP%dirLUT%pspace%theta(itheta)),  &
                    LUT   => OPP%dirLUT,                            &
                    fname => OPP%dirLUT%fname(iphi,itheta),         &
                    S     => OPP%dirLUT%S(iphi,itheta),             &
                    T     => OPP%dirLUT%T(iphi,itheta) )

        errcnt=0
        if( allocated(S%c) .or. .not.angle_mask(iphi,itheta) ) cycle

        if(OPP%optprop_LUT_debug .and. myid.eq.0) &
          print *,'Loading direct LUT from ',fname,' for azi',phi,': sza :',theta,':: ',descr

        write(str(2),FMT='("phi",I0)')   phi
        write(str(3),FMT='("theta",I0)') theta

        if(myid.eq.0) then
          if(OPP%optprop_LUT_debug) print *,'Trying to load the LUT from file... ',fname
            write(str(4),FMT='(A)') 'S' ; call ncload([fname,str(1),str(2),str(3),str(4)],S%c,iierr) ; errcnt = errcnt+iierr

!            if(OPP%optprop_LUT_debug) print *,'loaded the LUT from file...',[trim(fname),trim(str(1)),trim(str(2)),trim(str(3)),trim(str(4)),trim(str(5)),trim(str(6))]!,S%c
            if(iierr.eq.0) then
              if(any( S%c.gt.one ).or.any(S%c.lt.zero) ) errcnt=errcnt+100
            endif

            write(str(4),FMT='(A)') 'T' ; call ncload([fname,str(1),str(2),str(3),str(4)],T%c,iierr) ; errcnt = errcnt+iierr
!            if(OPP%optprop_LUT_debug) print *,'loaded the LUT from file...',[trim(fname),trim(str(1)),trim(str(2)),trim(str(3)),trim(str(4)),trim(str(5)),trim(str(6))]!,S%c
            if(iierr.eq.0) then
              if(any( T%c.gt.one ).or.any(T%c.lt.zero) ) errcnt=errcnt+200
              call check_dirLUT_matches_pspace(LUT,fname)
            endif


            ! Check if the precision requirements are all met and if we can load the %stddev_tol array
            write(str(4),FMT='(A)') 'S_tol' ; call ncload([fname,str(1),str(2),str(3),str(4)],S%stddev_tol,iierr)
            if(iierr.ne.0) then
              lstddev_inbounds = .False. ! if we could not load stddev...
            else
              lstddev_inbounds = .True. ! first assume that precision is met and then check if this is still the case...
              if(lstddev_inbounds) lstddev_inbounds = iierr.eq.i0
              if(lstddev_inbounds) lstddev_inbounds = all(S%stddev_tol.le.stddev_atol+epsilon(stddev_atol)*10)

              if(.not.lstddev_inbounds) &
                  print *,'Loading of direct tables S :: ',phi,theta,' failed bc lstddev_inbounds',lstddev_inbounds,'::',maxval(S%stddev_tol),'(',stddev_atol+epsilon(stddev_atol)*10,')',errcnt

              write(str(4),FMT='(A)') 'T_tol' ; call ncload([fname,str(1),str(2),str(3),str(4)],T%stddev_tol,iierr)

              if(lstddev_inbounds) lstddev_inbounds = iierr.eq.i0
              if(lstddev_inbounds) lstddev_inbounds = all(T%stddev_tol.le.stddev_atol+epsilon(stddev_atol)*10)

              if(.not.lstddev_inbounds) &
                  print *,'Loading of direct tables T :: ',phi,theta,' failed bc lstddev_inbounds',lstddev_inbounds,'::',maxval(T%stddev_tol),'(',stddev_atol+epsilon(stddev_atol)*10,')',errcnt

              if(OPP%optprop_LUT_debug) &
                  print *,'Tried to load the LUT from file... result is errcnt:',errcnt,'lstddev_inbounds',lstddev_inbounds,':',trim(str(1)),trim(str(2)),trim(str(3))
            endif
        endif !master

        call mpi_bcast(errcnt           , 1_mpiint , imp_int     , 0_mpiint , comm , mpierr); call CHKERR(mpierr)
        call mpi_bcast(lstddev_inbounds , 1_mpiint , imp_logical , 0_mpiint , comm , mpierr); call CHKERR(mpierr)

        if((errcnt.ne.0) .or. (.not.lstddev_inbounds) ) then
          if(myid.eq.0) then ! master -- setup netcdf files:
            write(str(4),FMT='(A)') 'pspace'
            write(str(5),FMT='(A)') "range_aspect" ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%range_aspect, iierr)
            write(str(5),FMT='(A)') "range_tau   " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%range_tau   , iierr)
            write(str(5),FMT='(A)') "range_w0    " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%range_w0    , iierr)
            write(str(5),FMT='(A)') "range_g     " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%range_g     , iierr)
            write(str(5),FMT='(A)') "range_phi   " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%range_phi   , iierr)
            write(str(5),FMT='(A)') "range_theta " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%range_theta , iierr)

            write(str(5),FMT='(A)') "aspect     " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%aspect      , iierr)
            write(str(5),FMT='(A)') "tau        " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%tau         , iierr)
            write(str(5),FMT='(A)') "w0         " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%w0          , iierr)
            write(str(5),FMT='(A)') "g          " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%g           , iierr)
            write(str(5),FMT='(A)') "phi        " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%phi         , iierr)
            write(str(5),FMT='(A)') "theta      " ; call ncwrite([fname , str(1),str(4),str(5) ] , LUT%pspace%theta       , iierr)

            write(varname(1),FMT='(A)') "S     "
            write(varname(2),FMT='(A)') "S_tol"
            write(varname(3),FMT='(A)') "T     "
            write(varname(4),FMT='(A)') "T_tol"

            if(.not.allocated(S%table_name_c  ) ) allocate(S%table_name_c  (5))
            if(.not.allocated(S%table_name_tol) ) allocate(S%table_name_tol(5))
            if(.not.allocated(T%table_name_c  ) ) allocate(T%table_name_c  (5))
            if(.not.allocated(T%table_name_tol) ) allocate(T%table_name_tol(5))

            S%table_name_c   = [fname,str(1),str(2),str(3),varname(1)]
            S%table_name_tol = [fname,str(1),str(2),str(3),varname(2)]
            T%table_name_c   = [fname,str(1),str(2),str(3),varname(3)]
            T%table_name_tol = [fname,str(1),str(2),str(3),varname(4)]
          endif !master

          ! all call createLUT
          call OPP%createLUT_dir(LUT, comm, iphi, itheta)

        endif

        if(myid.eq.0) then
          deallocate(S%stddev_tol)
          deallocate(T%stddev_tol)

          if(OPP%optprop_LUT_debug) print *,'Done loading direct LUT, for phi/theta:',phi,theta,':: shape(T)',shape(T%c),':: shape(S)',shape(S%c)
        endif !master

      end associate
      enddo !iphi
    enddo! itheta

end subroutine

subroutine createLUT_diff(OPP, LUT, comm)
    class(t_optprop_LUT) :: OPP
    type(diffuseTable) :: LUT
    integer(mpiint),intent(in) :: comm

    logical :: gotmsg
    integer(mpiint) :: status(MPI_STATUS_SIZE)
    integer(iintegers) :: workinput(5) !isrc, iaspect, itauz, iw0, ig
    integer(iintegers) :: idummy, workindex
    real(ireals) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(ireals) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)

    integer(mpiint), parameter :: READYMSG=1,HAVERESULTSMSG=2, WORKMSG=3, FINALIZEMSG=4, RESULTMSG=5

    if(myid.eq.0) then
      select type(OPP)
        class is (t_optprop_LUT_8_10)
          call prepare_table_space(LUT%S,OPP%diff_streams**2)

        class is (t_optprop_LUT_1_2)
          stop 'Twostream LUT generation is broken at the  moment -- need to write loadbalancer stuff'
          call prepare_table_space(LUT%S,2_iintegers)

        class default
          stop 'dont know optprop_class'
      end select
    endif

    if(myid.le.0 .and. numnodes.le.1) stop 'At the moment creation of diffuse Lookuptable needs at least two mpi-ranks to work... please run with more ranks.'
    if(myid.eq.0) then
      call master(LUT%S)
    else
      call worker()
    endif

    if(myid.eq.0) print *,'done calculating diffuse coefficients',shape(LUT%S%c)
    contains
      subroutine master(S)
        type(table),intent(inout) :: S
        integer(iintegers) :: total_size, cnt, finalizedworkers

        integer(iintegers),allocatable :: allwork(:,:) ! dimension (N,size(workinput)) ==> vector over work dimensions and 5 integers
        integer(iintegers) :: iaspect, itauz, iw0, ig
        integer(iintegers) :: isrc,idst,ind

        logical :: ldone(OPP%diff_streams)

        finalizedworkers=0
        total_size = OPP%Ng*OPP%Naspect*OPP%Ntau *OPP%Nw0 *OPP%diff_streams

        allocate( allwork(total_size, size(workinput) ) )
        cnt=1
        do ig = 1,OPP%Ng
          do iw0   = 1,OPP%Nw0
            do itauz = 1,OPP%Ntau
              do iaspect = 1,OPP%Naspect
                do isrc = 1,OPP%diff_streams
                  allwork(cnt, :) = [isrc,iaspect, itauz, iw0, ig]
                  cnt=cnt+1
                enddo
              enddo
            enddo
          enddo
        enddo

        cnt=1
        do

          ! Check if we already calculated the coefficients
          if(cnt.le.total_size) then
            isrc    = allwork(cnt, 1)
            iaspect = allwork(cnt, 2)
            itauz   = allwork(cnt, 3)
            iw0     = allwork(cnt, 4)
            ig      = allwork(cnt, 5)

            do idst = 1,OPP%diff_streams
                ind = (idst-1)*OPP%diff_streams + isrc
                ldone(idst) = ( ( S%c         ( ind, iaspect, itauz, iw0, ig ).ge.zero)            &
                          .and. ( S%c         ( ind, iaspect, itauz, iw0, ig ).le.one )            &
                          .and. ( S%stddev_tol( ind, iaspect, itauz, iw0, ig ).le.stddev_atol ) )
            enddo


            if(all(ldone)) then
              if( mod(cnt-1, total_size/100).eq.0 ) & !every 1 percent report status
                  print *,'Resuming from diffuse LUT...',cnt/(total_size/100),'%'
              cnt=cnt+1
              cycle
            endif
          endif

          ! Now that we know we got something to do, lets find a suitable worker
          gotmsg=.False.
          call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

          if (gotmsg) then

            select case (status(MPI_TAG))

            case(READYMSG)

              ! capture the READY MSG -- we should not leave messages hanging around.
              call mpi_recv(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), READYMSG, comm, status, mpierr) ; call CHKERR(mpierr)

              if(cnt.le.total_size) then ! we got something to do for a worker -- send him...
                call mpi_send(cnt, 1_mpiint, imp_int, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(allwork(cnt,:), size(allwork(cnt,:)), imp_int, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)

              else ! no more work to do... tell the worker to quit
                call mpi_send(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), FINALIZEMSG, comm, mpierr); call CHKERR(mpierr)
              endif
              cnt = cnt+1

            case(HAVERESULTSMSG)
              call mpi_recv(workindex, 1_mpiint, imp_int, status(MPI_SOURCE), HAVERESULTSMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_diff, size(S_diff), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_tol , size(S_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
!              call mpi_recv(T_dir , size(T_dir ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
!              call mpi_recv(T_tol , size(T_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)

              ! Sort coefficients into destination ordering and put em in LUT
              isrc    = allwork(workindex, 1)
              iaspect = allwork(workindex, 2)
              itauz   = allwork(workindex, 3)
              iw0     = allwork(workindex, 4)
              ig      = allwork(workindex, 5)

              do idst = 1, OPP%diff_streams
                ind = (idst-1)*OPP%diff_streams + isrc
                S%c         ( ind, iaspect, itauz, iw0, ig) = S_diff(idst)
                S%stddev_tol( ind, iaspect, itauz, iw0, ig) = S_tol (idst)
              enddo

              if( mod(workindex-1, total_size/100).eq.0 ) & !every 1 percent report status
                  print *,'Calculated diffuse LUT...',(100*(workindex-1)/total_size),'%'

              if( mod(workindex-1, total_size/100).eq.0 ) then !every 1 percent of LUT dump it.
                print *,'Writing diffuse table to file...',workindex,total_size
                call ncwrite(S%table_name_c  , S%c         ,iierr)
                call ncwrite(S%table_name_tol, S%stddev_tol,iierr)
                print *,'done writing!',iierr
              endif

            case(FINALIZEMSG)
              call mpi_recv(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
              finalizedworkers = finalizedworkers+1
              if(finalizedworkers.eq.numnodes-1) then

                gotmsg=.False.
                call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)
                if(gotmsg) then
                  print *,'Found message where I would not think there should be one',status
                  stop 'error'
                endif

                exit ! all work is done
              endif


            end select
          endif
        enddo

        print *,'Writing diffuse table to file...'
        call ncwrite(S%table_name_c  , S%c         ,iierr)
        call ncwrite(S%table_name_tol, S%stddev_tol,iierr)
        print *,'done writing!',iierr,':: max_atol',maxval(S%stddev_tol)
      end subroutine
      subroutine worker()
          ! workers send READY message to master
          call mpi_send(-i1, 1_mpiint, imp_int, 0_mpiint, READYMSG, comm, mpierr); call CHKERR(mpierr)

          do
            ! ask what to do
            gotmsg=.False.
            call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

            if (gotmsg) then

              select case (status(MPI_TAG))

              case(WORKMSG)
                ! wait for work to arrive
                call mpi_recv( workindex, 1_mpiint, imp_int, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_recv( workinput , size(workinput), imp_int, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)

                call OPP%bmc_wrapper(workinput(1),  &
                    LUT%pspace%aspect(workinput(2)),&
                    LUT%pspace%tau(workinput(3)),   &
                    LUT%pspace%w0(workinput(4)),    &
                    LUT%pspace%g(workinput(5)),     &
                    .False.,zero,zero,mpi_comm_self,&
                    S_diff,T_dir, S_tol,T_tol)

                call mpi_send(workindex, 1_mpiint, imp_int, status(MPI_SOURCE), HAVERESULTSMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(S_diff, size(S_diff), imp_real, status(MPI_SOURCE), RESULTMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(S_tol , size(S_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, mpierr); call CHKERR(mpierr)
!                call mpi_send(T_dir , size(T_dir ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, mpierr); call CHKERR(mpierr)
!                call mpi_send(T_tol , size(T_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, mpierr); call CHKERR(mpierr)

                call mpi_send(-i1, 1_mpiint, imp_int, 0_mpiint, READYMSG, comm, mpierr); call CHKERR(mpierr)

              case(FINALIZEMSG)
                call mpi_recv(idummy, 1_mpiint, imp_int, 0_mpiint, FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_send(-i1, 1_mpiint, imp_int, 0_mpiint, FINALIZEMSG, comm, mpierr); call CHKERR(mpierr)
                exit

              end select

            endif !gotmsg
          enddo
      end subroutine
      subroutine prepare_table_space(S,Ncoeff)
          type(table) :: S
          integer(iintegers) :: Ncoeff
          if(.not. allocated(S%stddev_tol) ) then
            allocate(S%stddev_tol(Ncoeff, OPP%Naspect,OPP%Ntau ,OPP%Nw0, OPP%Ng))
            S%stddev_tol = 1e8_ireals
            call ncwrite(S%table_name_tol, S%stddev_tol, iierr)
          endif

          if(.not.allocated(S%c) ) then
            allocate(S%c(Ncoeff, OPP%Naspect,OPP%Ntau ,OPP%Nw0, OPP%Ng))
            S%c = nil
            call ncwrite(S%table_name_c,S%c,iierr)
          endif
      end subroutine
end subroutine
subroutine createLUT_dir(OPP,LUT, comm, iphi,itheta)
    class(t_optprop_LUT) :: OPP
    type(directTable) :: LUT
    integer(mpiint),intent(in) :: comm
    integer(iintegers),intent(in) :: iphi,itheta

    logical :: gotmsg
    integer(mpiint) :: status(MPI_STATUS_SIZE)
    integer(iintegers) :: workinput(5)
    integer(iintegers) :: idummy, workindex
    real(ireals) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(ireals) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)

    integer(mpiint), parameter :: READYMSG=1,HAVERESULTSMSG=2, WORKMSG=3, FINALIZEMSG=4, RESULTMSG=5

    if(myid.eq.0) &
      call prepare_table_space(LUT%S(iphi,itheta), OPP%diff_streams, LUT%T(iphi,itheta), OPP%dir_streams)

    if(myid.le.0 .and. numnodes.le.1) &
      stop 'At the moment creation of direct Lookuptable needs at least two mpi-ranks to work... please run with more ranks.'

    if(myid.eq.0) then
      call master( LUT%S(iphi,itheta), LUT%T(iphi,itheta) )
      print *,'done calculating direct coefficients'
    else
      call worker()
    endif

    contains
      subroutine master(S,T)
        type(table),intent(inout) :: S,T

        integer(iintegers) :: total_size, cnt, finalizedworkers

        integer(iintegers),allocatable :: allwork(:,:) ! dimension (N,size(workinput)) ==> vector over work dimensions and 5 integers
        integer(iintegers) :: iaspect, itauz, iw0, ig
        integer(iintegers) :: isrc,idst,ind

        logical :: ldoneS(OPP%diff_streams), ldoneT(OPP%dir_streams)

        finalizedworkers=0
        total_size = OPP%Ng * OPP%Naspect * OPP%Ntau * OPP%Nw0 * OPP%dir_streams

        allocate(allwork(total_size, size(workinput)))
        cnt=1
        do ig = 1,OPP%Ng
          do iw0 = 1,OPP%Nw0
            do itauz = 1,OPP%Ntau
              do iaspect = 1,OPP%Naspect
                do isrc = 1,OPP%dir_streams
                  allwork(cnt, :) = [isrc, iaspect, itauz ,iw0, ig]
                  cnt=cnt+1
                enddo
              enddo
            enddo
          enddo
        enddo

        cnt=1
        do

          ! Check if we already calculated the coefficients
          if(cnt.le.total_size) then
            isrc  = allwork(cnt, 1)
            iaspect = allwork(cnt, 2)
            itauz = allwork(cnt, 3)
            iw0   = allwork(cnt, 4)
            ig    = allwork(cnt, 5)

            do idst = 1,OPP%diff_streams
                ind = (idst-1)*OPP%dir_streams + isrc
                ldoneS(idst) = ( ( S%c         ( ind, iaspect, itauz, iw0, ig ).ge.zero)            &
                           .and. ( S%c         ( ind, iaspect, itauz, iw0, ig ).le.one )            &
                           .and. ( S%stddev_tol( ind, iaspect, itauz, iw0, ig ).le.stddev_atol ) )
            enddo
            do idst = 1,OPP%dir_streams
                ind = (idst-1)*OPP%dir_streams + isrc
                ldoneT(idst) = ( ( T%c         ( ind, iaspect, itauz, iw0, ig ).ge.zero)            &
                           .and. ( T%c         ( ind, iaspect, itauz, iw0, ig ).le.one )            &
                           .and. ( T%stddev_tol( ind, iaspect, itauz, iw0, ig ).le.stddev_atol ) )
            enddo

            if( all(ldoneS) .and. all(ldoneT) ) then
              if( mod(cnt-1, total_size/100).eq.0 ) & !every 1 percent report status
                  print *,'Resuming from direct LUT(',int(LUT%pspace%phi(iphi)),int(LUT%pspace%theta(itheta)),')... ',cnt/(total_size/100),'%'
              cnt=cnt+1
              cycle
            endif
          endif

          ! Now that we know we got something to do, lets find a suitable worker
          gotmsg=.False.
          call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

          if (gotmsg) then

            select case (status(MPI_TAG))

            case(READYMSG)

              ! capture the READY MSG -- we should not leave messages hanging around.
              call mpi_recv(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), READYMSG, comm, status, mpierr) ; call CHKERR(mpierr)

              if(cnt.le.total_size) then ! we got something to do for a worker -- send him...
                call mpi_send(cnt, 1_mpiint, imp_int, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(allwork(cnt,:), size(allwork(cnt,:)), imp_int, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)

              else ! no more work to do... tell the worker to quit
                call mpi_send(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), FINALIZEMSG, comm, mpierr); call CHKERR(mpierr)
              endif
              cnt = cnt+1

            case(HAVERESULTSMSG)
              call mpi_recv(workindex, 1_mpiint, imp_int, status(MPI_SOURCE), HAVERESULTSMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_diff, size(S_diff), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_dir , size(T_dir ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_tol , size(S_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_tol , size(T_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)

              isrc    = allwork(workindex, 1)
              iaspect = allwork(workindex, 2)
              itauz   = allwork(workindex, 3)
              iw0     = allwork(workindex, 4)
              ig      = allwork(workindex, 5)

              !print *,myid,'Saving values for ',workindex,'::',isrc,iaspect,itauz,iw0,ig,'::',T_dir
              !print *,myid,'supposed optprop:', LUT%pspace%aspect(iaspect),LUT%pspace%tau(itauz),LUT%pspace%w0(iw0),LUT%pspace%g(ig),LUT%pspace%phi(iphi),LUT%pspace%theta(itheta)

              ! Sort coefficients into destination ordering and put em in LUT
              do idst = 1, OPP%diff_streams
                ind = (idst-1)*OPP%dir_streams + isrc
                S%c         (ind, iaspect, itauz, iw0, ig) = S_diff(idst)
                S%stddev_tol(ind, iaspect, itauz, iw0, ig) = S_tol (idst)
              enddo
              do idst = 1, OPP%dir_streams
                ind = (idst-1)*OPP%dir_streams + isrc
                T%c         (ind, iaspect, itauz, iw0, ig) = T_dir (idst)
                T%stddev_tol(ind, iaspect, itauz, iw0, ig) = T_tol (idst)
              enddo

              !do idst = 1, OPP%dir_streams
              !  print *,myid,'T%c for idst',idst,T%c((idst-1)*OPP%dir_streams+1:idst*OPP%dir_streams, iaspect, itauz, iw0, ig)
              !enddo

              if( mod(workindex-1, total_size/100).eq.0 ) & !every 1 percent report status
                  print *,'Calculated direct LUT(',int(LUT%pspace%phi(iphi)),int(LUT%pspace%theta(itheta)),')...',(100*(workindex-1))/total_size,'%'

              if( mod(workindex-1, total_size/100 ).eq.0 ) then !every 1 percent of LUT dump it.
                print *,'Writing direct table to file...'
                call ncwrite(S%table_name_c  , S%c         ,iierr)
                call ncwrite(S%table_name_tol, S%stddev_tol,iierr)
                call ncwrite(T%table_name_c  , T%c         ,iierr)
                call ncwrite(T%table_name_tol, T%stddev_tol,iierr)
                print *,'done writing!',iierr
              endif

            case(FINALIZEMSG)
              call mpi_recv(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
              finalizedworkers = finalizedworkers+1
              if(finalizedworkers.eq.numnodes-1) exit ! all work is done

            end select
          endif
        enddo

        print *,'Writing direct table to file...'
        call ncwrite(S%table_name_c  , S%c         ,iierr)
        call ncwrite(S%table_name_tol, S%stddev_tol,iierr)
        call ncwrite(T%table_name_c  , T%c         ,iierr)
        call ncwrite(T%table_name_tol, T%stddev_tol,iierr)
        print *,'done writing!',iierr,':: max_atol S',maxval(S%stddev_tol),'max_atol T',maxval(T%stddev_tol)
      end subroutine
      subroutine worker()
          ! workers send READY message to master
          call mpi_send(-i1, 1_mpiint, imp_int, 0_mpiint, READYMSG, comm, mpierr); call CHKERR(mpierr)

          do
            ! ask what to do
            gotmsg=.False.
            call mpi_iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, gotmsg, status, mpierr); call CHKERR(mpierr)

            if (gotmsg) then

              select case (status(MPI_TAG))

              case(WORKMSG)
                ! wait for work to arrive
                call mpi_recv( workindex, 1_mpiint, imp_int, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_recv( workinput , size(workinput), imp_int, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)

                call OPP%bmc_wrapper(workinput(1),  &
                    LUT%pspace%aspect(workinput(2)),&
                    LUT%pspace%tau(workinput(3)),   &
                    LUT%pspace%w0(workinput(4)),    &
                    LUT%pspace%g(workinput(5)),     &
                    .True. ,                        &
                    LUT%pspace%phi(iphi),           &
                    LUT%pspace%theta(itheta),       &
                    mpi_comm_self,                  &
                    S_diff,T_dir,S_tol,T_tol)

                  !print *,myid,'Computed values for ',workindex,'::',workinput,'::',T_dir
                  !print *,myid,'optprop:', LUT%pspace%aspect(workinput(2)),LUT%pspace%tau(workinput(3)),LUT%pspace%w0(workinput(4)),LUT%pspace%g(workinput(5)),LUT%pspace%phi(iphi),LUT%pspace%theta(itheta)

                call mpi_send(workindex , 1_mpiint     , imp_int  , status(MPI_SOURCE) , HAVERESULTSMSG , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(S_diff    , size(S_diff) , imp_real , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(T_dir     , size(T_dir ) , imp_real , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(S_tol     , size(S_tol ) , imp_real , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(T_tol     , size(T_tol ) , imp_real , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)

                call mpi_send(-i1       , 1_mpiint     , imp_int  , 0_mpiint           , READYMSG       , comm , mpierr); call CHKERR(mpierr)

              case(FINALIZEMSG)
                call mpi_recv(idummy, 1_mpiint, imp_int, 0_mpiint, FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_send(-i1, 1_mpiint, imp_int, 0_mpiint, FINALIZEMSG, comm, mpierr); call CHKERR(mpierr)
                exit

              end select

            endif !gotmsg
          enddo
      end subroutine
      subroutine prepare_table_space(S, NcoeffS, T, NcoeffT )
          type(table),intent(inout) :: S,T
          integer(iintegers),intent(in) :: NcoeffS,NcoeffT

          integer(iintegers) :: errcnt

          print *,'preparing space for LUT for direct coeffs at angles: ',iphi,itheta
          errcnt=0

          if(.not.allocated(S%stddev_tol) ) then
            allocate(S%stddev_tol(NcoeffS*NcoeffT, OPP%Naspect, OPP%Ntau ,OPP%Nw0, OPP%Ng))
            S%stddev_tol = 1e8_ireals
            call ncwrite(S%table_name_tol, S%stddev_tol, iierr); errcnt = errcnt +iierr
          endif

          if(.not.allocated(T%stddev_tol) ) then
            allocate(T%stddev_tol(NcoeffT**2, OPP%Naspect, OPP%Ntau, OPP%Nw0, OPP%Ng))
            T%stddev_tol = 1e8_ireals
            call ncwrite(T%table_name_tol, T%stddev_tol, iierr); errcnt = errcnt +iierr
          endif

          if(.not. allocated(S%c) ) then
            allocate(S%c(NcoeffS*NcoeffT, OPP%Naspect, OPP%Ntau, OPP%Nw0, OPP%Ng))
            S%c = nil
            call ncwrite(S%table_name_c, S%c,iierr); errcnt = errcnt +iierr
          endif

          if(.not. allocated(T%c) ) then
            allocate(T%c(NcoeffT**2, OPP%Naspect, OPP%Ntau, OPP%Nw0, OPP%Ng))
            T%c = nil
            call ncwrite(T%table_name_c, T%c,iierr); errcnt = errcnt +iierr
          endif

          if(errcnt.ne.0) stop 'createLUT_dir :: could somehow not write to file... exiting...'
      end subroutine
end subroutine

  subroutine scatter_LUTtables(OPP, comm, azis, szas)
      integer(mpiint) ,intent(in) :: comm
      class(t_optprop_LUT) :: OPP
      real(ireals),intent(in) :: szas(:),azis(:)
      integer(iintegers) :: iphi,itheta
      logical :: angle_mask(OPP%Nphi,OPP%Ntheta)

      call determine_angles_to_load(comm, OPP%interp_mode, OPP%dirLUT, azis, szas, angle_mask)

      do itheta=1,OPP%Ntheta
        do iphi  =1,OPP%Nphi
          if(.not.angle_mask(iphi,itheta) ) cycle                               ! this LUT is not needed, skip....
          if ( mpi_logical_and(comm, allocated(OPP%dirLUT%T(iphi,itheta)%c)) ) cycle ! if all nodes have the LUT already, we dont need to scatter it...

          call imp_bcast(comm, OPP%dirLUT%T(iphi,itheta)%c, 0_mpiint)  ! DIRECT 2 DIRECT
          call imp_bcast(comm, OPP%dirLUT%S(iphi,itheta)%c, 0_mpiint)  ! DIRECT 2 DIFFUSE
        enddo
      enddo

      ! DIFFUSE 2 DIFFUSE
      if( mpi_logical_or(comm, .not.allocated(OPP%diffLUT%S%c) )) & ! then ! if one or more nodes do not have it, guess we have to send it...
        call imp_bcast(comm, OPP%diffLUT%S%c, 0_mpiint)

  end subroutine
subroutine bmc_wrapper(OPP, src, aspect, tauz, w0, g, dir, phi, theta, comm, S_diff, T_dir, S_tol, T_tol)
    class(t_optprop_LUT) :: OPP
    integer(iintegers),intent(in) :: src
    integer(mpiint),intent(in) :: comm
    logical,intent(in) :: dir
    real(ireals),intent(in) :: aspect, tauz, w0, g, phi, theta
    real(ireals) :: dx,dy

    real(ireals),intent(out) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(ireals),intent(out) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)

    real(ireals) :: bg(3)
    real(ireals), parameter :: dz = 100

    dx = dz / aspect
    dy = dx

    bg(1) = tauz / dz * (one-w0)
    bg(2) = tauz / dz * w0
    bg(3) = g

    S_diff=nil
    T_dir=nil

    !print *,comm,'BMC :: calling bmc_get_coeff',bg,'src',src,'phi/theta',phi,theta,dz
    call OPP%bmc%get_coeff(comm,bg,src,dir,phi,theta,dx,dy,dz,S_diff,T_dir,S_tol,T_tol)
    !print *,'BMC :: dir',T_dir,'diff',S_diff
end subroutine

  subroutine check_diffLUT_matches_pspace(LUT)
      type(diffuseTable),intent(in) :: LUT
      real(ireals),allocatable :: buf(:)
      character(default_str_len) :: str(3)
      integer(iintegers) align(3);
      write(str(1),FMT='(A)') "diffuse"
      write(str(2),FMT='(A)') "pspace"
      align=0

      write(str(3),FMT='(A)') "aspect  " ; call ncload([LUT%fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf, LUT%pspace%aspect)) align(1 ) =1 ; if(allocated(buf))deallocate(buf )
      write(str(3),FMT='(A)') "tau     " ; call ncload([LUT%fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf, LUT%pspace%tau)   ) align(1 ) =1 ; if(allocated(buf))deallocate(buf )
      write(str(3),FMT='(A)') "w0      " ; call ncload([LUT%fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf, LUT%pspace%w0 )   ) align(2 ) =1 ; if(allocated(buf))deallocate(buf )
      write(str(3),FMT='(A)') "g       " ; call ncload([LUT%fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf, LUT%pspace%g  )   ) align(3 ) =1 ; if(allocated(buf))deallocate(buf )

      if(any(align.ne.0)) stop 'parameter space of direct LUT coefficients is not aligned!'
  end subroutine
  subroutine check_dirLUT_matches_pspace(LUT,fname)
      type(directTable),intent(in) :: LUT
      character(len=*),intent(in) :: fname
      real(ireals),allocatable :: buf(:)
      character(default_str_len) :: str(3)
      integer(iintegers) align(5);
      write(str(1),FMT='(A)') "direct"
      write(str(2),FMT='(A)') "pspace"
      align=0

      write(str(3),FMT='(A)') "aspect  "  ; call ncload([fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf,LUT%pspace%aspect)  ) align(1) =1 ; if(allocated(buf)) deallocate(buf )
      write(str(3),FMT='(A)') "tau     "  ; call ncload([fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf,LUT%pspace%tau   )  ) align(1) =1 ; if(allocated(buf)) deallocate(buf )
      write(str(3),FMT='(A)') "w0      "  ; call ncload([fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf,LUT%pspace%w0    )  ) align(2) =1 ; if(allocated(buf)) deallocate(buf )
      write(str(3),FMT='(A)') "g       "  ; call ncload([fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf,LUT%pspace%g     )  ) align(3) =1 ; if(allocated(buf)) deallocate(buf )
      write(str(3),FMT='(A)') "phi     "  ; call ncload([fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf,LUT%pspace%phi   )  ) align(4) =1 ; if(allocated(buf)) deallocate(buf )
      write(str(3),FMT='(A)') "theta   "  ; call ncload([fname,str(1),str(2),str(3)],buf,iierr ) ; if(.not.compare_same( buf,LUT%pspace%theta )  ) align(5) =1 ; if(allocated(buf)) deallocate(buf )

      if(any(align.ne.0)) stop 'parameter space of direct LUT coefficients is not aligned!'
    end subroutine
  function compare_same(a,b)
    !> @brief compare 2 arrays that they are approximatly the same and
    !> if not print them next to each other
    logical :: compare_same
    real(ireals),intent(in) :: a(:),b(:)
    integer(iintegers) :: k
    if( all( shape(a).ne.shape(b) ) ) then
      print *,'compare_same : was called with differing shapes!', shape(a),'::',shape(b)
      compare_same = .False.
      return
    endif
    compare_same = all( rel_approx( a,b,1e-1_ireals ) )
    if(.not. compare_same ) then
      print *,'Compare_Same :: Arrays do not have the same values:'
      do k=1,size(a)
        print *,'k=',k,'a',a(k),'b',b(k),'diff',a(k)-b(k)
      enddo
    endif
end function

subroutine determine_angles_to_load(comm, interp_mode, LUT, azis, szas, mask)
    integer(mpiint) ,intent(in) :: comm
    integer(iintegers),intent(in) :: interp_mode
    type(directTable) :: LUT
    real(ireals),intent(in) :: szas(:),azis(:) ! all solar zenith angles that happen in this scene
    logical,intent(out) :: mask(size(LUT%pspace%phi),size(LUT%pspace%theta)) ! boolean array, which LUT entries should be loaded

    integer(iintegers) :: itheta, iphi, itheta1, iphi1
    logical :: lneed_azi(2), lneed_sza(2)
    real(ireals) :: theta(2),phi(2) ! sza and azimuth angle

    mask = .False.
    ! Check if really need to load it... i.e. we only want to load angles which are necessary for this run.
    do itheta=1,size(LUT%pspace%theta)
      do iphi  =1,size(LUT%pspace%phi)

        iphi1   = min(size(LUT%pspace%phi),   iphi+1)
        itheta1 = min(size(LUT%pspace%theta), itheta+1)

        phi   = LUT%pspace%phi( [ iphi, iphi1 ] )
        theta = LUT%pspace%theta( [ itheta, itheta1 ]  )

        select case(interp_mode)

        case(1:3) ! only load nearest neighbors
            lneed_azi(1) = any( azis .ge. phi(1)       .and. azis .le. sum(phi)/2 )
            lneed_azi(2) = any( azis .ge. sum(phi)/2   .and. azis .le. phi(2) )

        case(4) ! interpolate azimuth -- therefore also need surrounding tables
            lneed_azi = any( azis .ge. phi(1)       .and. azis .le. phi(2) )

        case default
            stop 'interpolation mode not implemented yet! please choose something else! '
        end select

        select case(interp_mode)

        case(1:2) ! only load nearest neighbors
            lneed_sza(1) = any( szas .ge. theta(1)     .and. szas .le. sum(theta)/2 )
            lneed_sza(2) = any( szas .ge. sum(theta)/2 .and. szas .le. theta(2) )

        case(3:4) ! interpolate sza -- therefore also need surrounding tables
            lneed_sza = any( szas .ge. theta(1)     .and. szas .le. theta(2) )

        case default
            stop 'interpolation mode not implemented yet! please choose something else! '
        end select

        !print *,'determine_angles_to_load: occuring azimuths',': phi,theta',phi,theta,'need_azi',lneed_azi,'lneed_sza',lneed_sza

        if( lneed_azi(1) .and. lneed_sza(1) ) mask(iphi  , itheta ) = .True.
        if( lneed_azi(1) .and. lneed_sza(2) ) mask(iphi  , itheta1) = .True.
        if( lneed_azi(2) .and. lneed_sza(1) ) mask(iphi1 , itheta ) = .True.
        if( lneed_azi(2) .and. lneed_sza(2) ) mask(iphi1 , itheta1) = .True.

        !if (all( lneed_azi .and. lneed_sza )) mask([iphi,iphi+1],[itheta,itheta+1]) = .True. !todo breaks if we need either theta+1 or phi+1 i.e. uneven sza or phi=90
        !if (all( lneed_azi .and. lneed_sza )) mask([iphi],[itheta]) = .True. !todo breaks if we need either theta+1 or phi+1 i.e. uneven sza or phi=90
      enddo
    enddo
    if(ldebug_optprop .and. myid.eq.0) then
      print *,'       phis',LUT%pspace%range_phi
      do itheta=1,size(LUT%pspace%theta)
        print *,'theta=',LUT%pspace%theta(itheta),' :: ',mask(:,itheta)
      enddo
    endif

    ! in case ranks would require different angles, we should broadcast this here. in principal ranks may only load the LUT they need, this approach may however not be the easiest to implement?

    do itheta=1,size(LUT%pspace%theta)
      do iphi  =1,size(LUT%pspace%phi)
        mask(iphi  , itheta ) = mpi_logical_or(comm, mask(iphi, itheta )) ! todo this is terrible that we send many messages instead of one bigger one...
      enddo
    enddo
end subroutine


function exp_index_to_param(index,range,N,expn)
    real(ireals) :: exp_index_to_param
    real(ireals),intent(in) :: index,range(2),expn
    integer(iintegers),intent(in) :: N
    real(ireals) :: expn1
    expn1=one/expn
    exp_index_to_param = lin_index_to_param( index, range**expn1, N) ** expn
end function
function lin_index_to_param(index,range,N)
    real(ireals) :: lin_index_to_param
    real(ireals),intent(in) :: index,range(2)
    integer(iintegers),intent(in) :: N
    if(N.gt.i1) then
      lin_index_to_param = range(1) + (index-one) * ( range(2)-range(1) ) / (N-1)
    else
      lin_index_to_param = range(1)
    endif
end function

subroutine set_parameter_space(OPP,ps)
    class(t_optprop_LUT) :: OPP
    type(parameter_space),intent(inout) :: ps
    real(ireals),parameter :: maximum_transmission=one-1e-7_ireals !one-epsilon(maximum_transmission) ! this parameter defines the lambert beer transmission we want the LUT to have given a pathlength of the box diameter
    real(ireals),parameter :: minimum_transmission=1e-30_ireals
    real(ireals) :: transmission
    integer(iintegers) :: k

    OPP%Naspect= Naspect
    OPP%Ntau   = Ntau
    OPP%Nw0    = Nw0
    OPP%Ng     = Ng
    OPP%Nphi   = Nphi
    OPP%Ntheta = Ntheta

    select type(OPP)
      class is (t_optprop_LUT_1_2)
          OPP%Nphi   = 1 ! azimithally average in 1D
          OPP%interp_mode = interp_mode_1_2
      class is (t_optprop_LUT_8_10)
          OPP%interp_mode = interp_mode_8_10
      class default
        stop 'set_parameter space: unexpected type for optprop_LUT object!'
    end select

    if(.not. allocated(ps%aspect)) allocate(ps%aspect(OPP%Naspect))
    if(.not. allocated(ps%tau   )) allocate(ps%tau   (OPP%Ntau   ))
    if(.not. allocated(ps%w0    )) allocate(ps%w0    (OPP%Nw0    ))
    if(.not. allocated(ps%g     )) allocate(ps%g     (OPP%Ng     ))
    if(.not. allocated(ps%phi   )) allocate(ps%phi   (OPP%Nphi   ))
    if(.not. allocated(ps%theta )) allocate(ps%theta (OPP%Ntheta ))

    ! -------------- Setup aspect support points

    do k=1,OPP%Naspect
      ps%aspect(k) = lin_index_to_param(one*k, ps%range_aspect, OPP%Naspect)
    enddo

    ! determine support points over a range
    ! ATTENTION, this is currently not good for tau...
    ! better use preset dimensions (happens at the end of the routine
    ! see in optprop_parameters for details

    ! -------------- Setup tau support points

    ps%range_tau   = [ -log(maximum_transmission), -log(minimum_transmission) ]

    do k=1,OPP%Ntau
      transmission = lin_index_to_param(one*k, [minimum_transmission, maximum_transmission], OPP%Ntau)
      ps%tau(k)    = -log(transmission)
    enddo

    ! -------------- Setup w0 support points

    do k=1,OPP%Nw0
      ps%w0(k) = lin_index_to_param(one*k, ps%range_w0, OPP%Nw0)
    enddo

    ! -------------- Setup g support points

    if(ldelta_scale) ps%range_g=[zero,.5_ireals]

    do k=1,OPP%Ng
      ps%g(k)     = exp_index_to_param(one*k,ps%range_g,OPP%Ng, .25_ireals )
    enddo

    if(OPP%Ng.eq.1) then
      ps%g(1)=zero
      ps%range_g=zero
    endif


    ! -------------- Setup phi/theta support points

    do k=1,OPP%Nphi
      ps%phi(k)   = lin_index_to_param(one*k,ps%range_phi,OPP%Nphi)
    enddo
    do k=1,OPP%Ntheta
      ps%theta(k) = lin_index_to_param(one*k,ps%range_theta,OPP%Ntheta)
    enddo


    ! ------------- Overwrite dimensions with preset values

    if (use_prescribed_LUT_dims) then
      ps%aspect     = preset_aspect
      ps%tau        = preset_tau
      ps%w0         = preset_w0
      ps%g          = preset_g
      ps%theta      = preset_theta

      ps%range_aspect = [ps%aspect(1), ps%aspect(Naspect)]
      ps%range_tau    = [ps%tau   (1), ps%tau   (Ntau   )]
      ps%range_w0     = [ps%w0    (1), ps%w0    (Nw0    )]
      ps%range_g      = [ps%g     (1), ps%g     (Ng     )]
      ps%range_theta  = [ps%theta (1), ps%theta (Ntheta )]
    endif
end subroutine

subroutine LUT_get_dir2dir(OPP, in_aspect, in_tauz, in_w0, g, phi, theta, C)
    class(t_optprop_LUT) :: OPP
    real(ireals),intent(in) :: in_aspect, in_tauz, in_w0, g, phi, theta
    real(ireals),intent(out):: C(:) ! dimension(OPP%dir_streams**2)
    real(ireals) :: aspect, tauz, w0, norm
    integer(iintegers) :: src

    real(ireals) :: pti(6)
    !real(ireals) :: vals(OPP%dir_streams**2, 2)

    aspect = in_aspect; tauz = in_tauz; w0 = in_w0
    if(ldebug_optprop) call catch_limits(OPP%dirLUT%pspace,aspect, tauz, w0, g)

    pti = get_indices_6d(aspect, tauz, w0, g, phi, theta, OPP%dirLUT%pspace)

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      C = OPP%dirLUT%T(nint(pti(5)), nint(pti(6)) )%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) )
    case(2)
      call interp_4d(pti, OPP%dirLUT%T(nint(pti(5)), nint(pti(6)) )%c, C)
    case(3)
      call interp_4p1d(pti([1,2,3,4,6]), OPP%dirLUT%T(nint(pti(5)), :), C)
    case(4)
      call interp_4p2d(pti, OPP%dirLUT%T(:, :), C)

    case default
      stop 'interpolation mode not implemented yet! please choose something else! '
    end select

    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%dir_streams
        norm = sum(C( src:size(C):OPP%dir_streams))
        if(real(norm).gt.one+1e-5_ireals) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        print *,'Error in dir2dir coeffs :: ierr',iierr,size(C),OPP%dir_streams,'::',C
        do src=1,OPP%dir_streams
          print *,'SUM dir2dir coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%dir_streams)),' :: coeff',C( src:size(C):OPP%dir_streams )
        enddo
        call exit(1)
      endif
    endif
end subroutine

subroutine LUT_get_dir2diff(OPP, in_aspect, in_tauz, in_w0, g, phi, theta, C)
    class(t_optprop_LUT) :: OPP
    real(ireals),intent(in) :: in_aspect, in_tauz, in_w0, g, phi, theta
    real(ireals),intent(out):: C(:) ! dimension(OPP%dir_streams*OPP%diff_streams)

    real(ireals) :: aspect, tauz, w0
    real(ireals) :: pti(6),norm
    integer(iintegers) :: src

    aspect = in_aspect; tauz = in_tauz; w0 = in_w0
    if(ldebug_optprop) then
      call catch_limits(OPP%dirLUT%pspace,aspect, tauz, w0, g)
      if(size(C).ne.OPP%dir_streams*OPP%diff_streams) stop 'LUT_get_dir2diff called with wrong array shape'
    endif

    pti = get_indices_6d(aspect, tauz, w0, g, phi, theta, OPP%dirLUT%pspace)

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      C = OPP%dirLUT%S( nint(pti(5)), nint(pti(6)) )%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) )

    case(2)
      call interp_4d(pti, OPP%dirLUT%S( nint(pti(5)), nint(pti(6)) )%c, C)

    case(3)
      call interp_4p1d(pti([1,2,3,4,6]), OPP%dirLUT%S(nint(pti(5)), :), C)

    case(4)
      call interp_4p2d(pti, OPP%dirLUT%S(:, :), C)

    case default
      stop 'interpolation mode not implemented yet! please choose something else! '
    end select


    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%diff_streams
        norm = sum( C( src:size(C):OPP%dir_streams ) )
        if(real(norm).gt.one+1e-5_ireals) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        print *,'Error in dir2diff coeffs :: ierr',iierr,':',in_aspect, in_tauz, in_w0, g, phi, theta,'::',C,'::',shape(OPP%dirLUT%S( nint(pti(5)), nint(pti(6)) )%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) ))
        print *,'Error in dir2dir coeffs :: ierr',iierr,'::',OPP%dirLUT%T( nint(pti(5)), nint(pti(6)) )%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) ),'::',shape(OPP%dirLUT%T( nint(pti(5)), nint(pti(6)) )%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) ))
        do src=1,OPP%diff_streams
          print *,'SUM dir2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%dir_streams)),' :: coeff',C( src:size(C):OPP%dir_streams )
        enddo
        call exit(1)
      endif
    endif
end subroutine
subroutine LUT_get_diff2diff(OPP, in_aspect, in_tauz, in_w0, g, C)
    class(t_optprop_LUT) :: OPP
    real(ireals),intent(in) :: in_aspect, in_tauz, in_w0, g
    real(ireals),intent(out):: C(:) ! dimension(OPP%diff_streams**2)

    real(ireals) :: aspect, tauz, w0
    real(ireals) :: pti(4),norm
    integer(iintegers) :: src

    aspect = in_aspect; tauz = in_tauz; w0 = in_w0
    if(ldebug_optprop) call catch_limits(OPP%diffLUT%pspace, aspect, tauz, w0, g)

    pti = get_indices_4d(aspect, tauz, w0, g, OPP%diffLUT%pspace)

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      C = OPP%diffLUT%S%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) )
    case(2:4)
      ! Linear interpolation
      call interp_4d(pti, OPP%diffLUT%S%c, C)
    case default
      stop 'interpolation mode not implemented yet! please choose something else! '
    end select

    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%diff_streams
        norm = sum( C( src:size(C):OPP%diff_streams ) )
        if(norm.gt.one+1e-5_ireals) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        print *,'Error in diff2diff coeffs :: ierr',iierr, ':', in_aspect, in_tauz, in_w0, g, '::', C
        do src=1,OPP%diff_streams
          print *,'SUM diff2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%diff_streams)),' :: coeff',C(src:size(C):OPP%diff_streams)
        enddo
        call exit(1)
      endif
    endif
end subroutine

subroutine interp_4p2d(pti,ctable,C)
        integer,parameter :: Ndim=6
        real(ireals),intent(in) :: pti(Ndim)
        type(table),intent(in) :: ctable(:,:)  ! contains Nphi, Ntheta databases
        real(ireals),intent(out) :: C(:)

        real(ireals) :: weights(Ndim)
        integer :: indices(2,2),fpti(Ndim)

        ! Instead of doing a full interpolation in 6 dimension we start out with
        ! 4 dimensions only at the cornerstones of the 4d hypercube
        real(ireals) :: C4(size(C),6)

        !stop 'todo atm not advisable to use this function .. we dont load the neighboring LUTs'

        ! First determine the array indices, where to look.
        fpti = floor(pti)
        weights = modulo(pti, one)

        indices(:,1) = max(i1, min( ubound(ctable,1), [0,1] +fpti(5) ) )
        indices(:,2) = max(i1, min( ubound(ctable,2), [0,1] +fpti(6) ) )

        call interp_4d( pti(1:4), ctable(indices(1,1), indices(1,2) )%c, C4(:,1) ) ! differing azimuth
        call interp_4d( pti(1:4), ctable(indices(2,1), indices(1,2) )%c, C4(:,2) ) !        "
        call interp_4d( pti(1:4), ctable(indices(1,1), indices(2,2) )%c, C4(:,3) )
        call interp_4d( pti(1:4), ctable(indices(2,1), indices(2,2) )%c, C4(:,4) )

        C4(:,5) = C4(:,1) + weights(5) * ( C4(:,2) - C4(:,1) )
        C4(:,6) = C4(:,3) + weights(5) * ( C4(:,4) - C4(:,3) )
        C       = C4(:,5) + weights(6) * ( C4(:,6) - C4(:,5) )
end subroutine
subroutine interp_4p1d(pti,ctable,C)
        integer,parameter :: Ndim=5
        real(ireals),intent(in) :: pti(Ndim)
        type(table),intent(in) :: ctable(:)  ! contains N databases
        real(ireals),intent(out) :: C(:)

        real(ireals) :: weights(Ndim)
        integer :: indices(2),fpti(Ndim)

        ! Instead of doing a full interpolation in 6 dimension we start out with
        ! 4 dimensions only at the cornerstones of the 4d hypercube
        real(ireals) :: C4(size(C),2)

        ! First determine the array indices, where to look.
        fpti = floor(pti)
        weights = modulo(pti, one)

        indices(:) = max(i1, min( ubound(ctable,1), [0,1] +fpti(5) ) )

        call interp_4d( pti(1:4), ctable(indices(1))%c, C4(:,1) ) ! differing zenith
        call interp_4d( pti(1:4), ctable(indices(2))%c, C4(:,2) ) !        "

        C = C4(:,1) + weights(5) * ( C4(:,2) - C4(:,1) )
end subroutine


function get_indices_4d(aspect, tauz, w0, g, ps)
    real(ireals) :: get_indices_4d(4)
    real(ireals),intent(in) :: aspect, tauz, w0, g
    type(parameter_space),intent(in) :: ps

    get_indices_4d(1) = search_sorted_bisection(ps%aspect, aspect)
    get_indices_4d(2) = search_sorted_bisection(ps%tau   , tauz)
    get_indices_4d(3) = search_sorted_bisection(ps%w0    , w0)
    get_indices_4d(4) = search_sorted_bisection(ps%g     , g)
end function
function get_indices_6d(aspect, tauz, w0, g, phi, theta, ps)
    real(ireals) :: get_indices_6d(6)
    real(ireals),intent(in) :: aspect, tauz, w0, g, phi, theta
    type(parameter_space),intent(in) :: ps

    get_indices_6d(1:4) = get_indices_4d(aspect, tauz, w0, g, ps)

    get_indices_6d(5) = search_sorted_bisection(ps%phi  ,phi )
    get_indices_6d(6) = search_sorted_bisection(ps%theta,theta)

end function

logical function valid_input(val,range)
    real(ireals),intent(in) :: val,range(2)
    if(val.lt.range(1) .or. val.gt.range(2) ) then
      valid_input=.False.
      print *,'ohoh, this val is not in the optprop database range!',val,'not in',range
    else
      valid_input=.True.
    endif
end function
subroutine catch_limits(ps, aspect, tauz, w0, g)
    type(parameter_space),intent(in) :: ps
    real(ireals),intent(in) :: aspect, tauz, w0, g

    iierr=0

    if( aspect.lt.ps%range_aspect(1) .or. aspect.gt.ps%range_aspect(2) ) then
      print *,'aspect ratio is not in LookUpTable Range',aspect, 'LUT range',ps%range_aspect
      iierr=iierr+1
    endif
    if( tauz.lt.ps%range_tau(1) .or. tauz.gt.ps%range_tau(2) ) then
      print *,'tau is not in LookUpTable Range',tauz, 'LUT range',ps%range_tau
      iierr=iierr+1
    endif
    if( w0.lt.ps%range_w0(1) .or. w0.gt.ps%range_w0(2) ) then
      print *,'w0 is not in LookUpTable Range',w0, 'LUT range',ps%range_w0
      iierr=iierr+1
    endif
    if( g.lt.ps%range_g(1) .or. g.gt.ps%range_g(2) ) then
      print *,'g is not in LookUpTable Range',g, 'LUT range',ps%range_g
      iierr=iierr+1
    endif
    if(iierr.ne.0) print*, 'The LookUpTable was asked to give a coefficient, it was not defined for. Please specify a broader range.',iierr
end subroutine

end module
