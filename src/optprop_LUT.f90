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

  use m_helper_functions, only : approx,  &
    rel_approx, imp_bcast,                &
    mpi_logical_and, mpi_logical_or,      &
    search_sorted_bisection, CHKERR, itoa

  use m_data_parameters, only : ireals, iintegers, &
    one, zero, i0, i1, i3, mpiint, nil, inil,      &
    imp_int, imp_real, imp_logical,                &
    default_str_len

  use m_optprop_parameters, only:         &
    ldebug_optprop, lut_basename,         &
    Naspect, Ntau, Nw0, Ng, Nphi, Ntheta, &
    interp_mode_1_2,interp_mode_8_10,     &
    interp_mode_3_6,interp_mode_3_10,     &
    interp_mode_wedge_5_8,                &
    ldelta_scale,delta_scale_truncate,    &
    stddev_atol, stddev_rtol,             &
    use_prescribed_LUT_dims,              &
    preset_aspect, preset_tau, preset_w0, &
    preset_g,                             &
    OPP_LUT_ALL_ANGLES, luse_memory_map

  use m_boxmc, only: t_boxmc,t_boxmc_8_10,t_boxmc_1_2, t_boxmc_3_6, t_boxmc_3_10, &
    t_boxmc_wedge_5_8
  use m_tenstream_interpolation, only: interp_4d
  use m_netcdfio

  use m_mmap, only : arr_to_mmap, munmap_mmap_ptr

  implicit none

  private
  public :: t_optprop_LUT, t_optprop_LUT_8_10,t_optprop_LUT_1_2,t_optprop_LUT_3_6, t_optprop_LUT_3_10, &
    t_optprop_LUT_wedge_5_8
  ! This module loads and generates the LUT-tables for Tenstream Radiation
  ! computations.
  ! It also holds functions for interpolation on the regular LUT grid.

  integer(mpiint) :: iierr
  integer(mpiint) :: mpierr

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
    integer(iintegers) :: dir_streams = inil, diff_streams = inil
    integer(iintegers) :: Naspect, Ntau, Nw0, Ng, Nphi, Ntheta, interp_mode
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
  type,extends(t_optprop_LUT) :: t_optprop_LUT_3_10
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_3_6
  end type
  type,extends(t_optprop_LUT) :: t_optprop_LUT_wedge_5_8
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
      integer(mpiint) :: comm_size, myid

      call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
      call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

      if(.not.allocated(OPP%bmc)) then
        select type (OPP)
          class is (t_optprop_LUT_1_2)
            OPP%dir_streams  =  1
            OPP%diff_streams =  2
            OPP%lutbasename=trim(lut_basename)//'_1_2.'
            allocate(t_boxmc_1_2::OPP%bmc)

          class is (t_optprop_LUT_8_10)
            OPP%dir_streams  = 8
            OPP%diff_streams = 10
            OPP%lutbasename=trim(lut_basename)//'_8_10.'
            allocate(t_boxmc_8_10::OPP%bmc)

          class is (t_optprop_LUT_3_10)
            OPP%dir_streams  = 3
            OPP%diff_streams = 10
            OPP%lutbasename=trim(lut_basename)//'_3_10.'
            allocate(t_boxmc_3_10::OPP%bmc)

          class is (t_optprop_LUT_3_6)
            OPP%dir_streams  = 3
            OPP%diff_streams = 6
            OPP%lutbasename=trim(lut_basename)//'_3_6.'
            allocate(t_boxmc_3_6::OPP%bmc)

          class is (t_optprop_LUT_wedge_5_8)
            OPP%dir_streams  = 5
            OPP%diff_streams = 8
            OPP%lutbasename=trim(lut_basename)//'_wedge_5_8.'
            allocate(t_boxmc_wedge_5_8::OPP%bmc)

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
      !write(descr,FMT='("diffuse.outin.aspect",I0,".tau",I0,".w0",I0,".g",I0,".delta_",L1,"_",F0.3)') &
      !                   OPP%Naspect, OPP%Ntau, OPP%Nw0, OPP%Ng, ldelta_scale, delta_scale_truncate
      descr = gen_lut_basename('diffuse', OPP%Naspect, OPP%Ntau, OPP%Nw0, OPP%Ng)

      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Loading diffuse LUT from ',trim(descr)
      OPP%diffLUT%fname = trim(OPP%lutbasename)//trim(descr)//'.nc'

      call OPP%loadLUT_diff(comm)

      ! Copy LUT to all ranks
      if(comm_size.gt.1) call OPP%scatter_LUTtables(comm, azis,szas)

      OPP%LUT_initialized=.True.
      if(OPP%optprop_LUT_debug .and. myid.eq.0) print *,'Done loading LUTs (shape diffLUT',shape(OPP%diffLUT%S%c),')'
  end subroutine

  function gen_lut_basename(prefix, Naspect, Ntau, Nw0, Ng, phi, theta) result(lutname)
    character(len=default_str_len) :: lutname
    character(len=*), intent(in) :: prefix
    integer(iintegers), intent(in) :: Naspect, Ntau, Nw0, Ng
    integer(iintegers), intent(in), optional :: phi, theta
    lutname = prefix//'.aspect'//itoa(Naspect)//'.tau'//itoa(Ntau)//'.w0'//itoa(Nw0)//'.g'//itoa(Ng)
    if(present(phi))   lutname = trim(lutname)//'.phi'//itoa(phi)
    if(present(theta)) lutname = trim(lutname)//'.theta'//itoa(theta)
  end function

subroutine loadLUT_diff(OPP, comm)
    class(t_optprop_LUT) :: OPP
    integer(mpiint),intent(in) :: comm
    integer(iintegers) :: errcnt
    character(default_str_len) :: str(3)
    logical :: lstddev_inbounds

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

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
            print *,'Loading of diffuse tables S :: failed bc lstddev_inbounds',lstddev_inbounds,':: is ',maxval(S%stddev_tol),'(should be ',stddev_atol+epsilon(stddev_atol)*10,')',errcnt

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

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(.not. allocated(OPP%dirLUT%fname) ) allocate( OPP%dirLUT%fname(OPP%Nphi,OPP%Ntheta) )
    if(.not. allocated(OPP%dirLUT%S    ) ) allocate( OPP%dirLUT%S    (OPP%Nphi,OPP%Ntheta) )
    if(.not. allocated(OPP%dirLUT%T    ) ) allocate( OPP%dirLUT%T    (OPP%Nphi,OPP%Ntheta) )

    write(str(1),FMT='(A)') 'direct'

    call determine_angles_to_load(comm, OPP%interp_mode, OPP%dirLUT, azis, szas, angle_mask)

    do itheta=1,OPP%Ntheta
      do iphi  =1,OPP%Nphi

        ! Set filename of LUT
        descr = gen_lut_basename('direct', OPP%Naspect, OPP%Ntau, OPP%Nw0, OPP%Ng, &
          int(OPP%dirLUT%pspace%phi(iphi), iintegers), int(OPP%dirLUT%pspace%theta(itheta), iintegers))

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
                  print *,'Loading of direct tables S :: ',phi,theta,' failed bc lstddev_inbounds',lstddev_inbounds,':: is',maxval(S%stddev_tol),'( should be',stddev_atol+epsilon(stddev_atol)*10,')',errcnt

              write(str(4),FMT='(A)') 'T_tol' ; call ncload([fname,str(1),str(2),str(3),str(4)],T%stddev_tol,iierr)

              if(lstddev_inbounds) lstddev_inbounds = iierr.eq.i0
              if(lstddev_inbounds) lstddev_inbounds = all(T%stddev_tol.le.stddev_atol+epsilon(stddev_atol)*10)

              if(.not.lstddev_inbounds) &
                  print *,'Loading of direct tables T :: ',phi,theta,' failed bc lstddev_inbounds',lstddev_inbounds,':: is',maxval(T%stddev_tol),'(should be',stddev_atol+epsilon(stddev_atol)*10,')',errcnt

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

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(myid.eq.0) then
      select type(OPP)
        class is (t_optprop_LUT_8_10)
          call prepare_table_space(LUT%S,OPP%diff_streams**2)

        class is (t_optprop_LUT_3_10)
          call prepare_table_space(LUT%S,OPP%diff_streams**2)

        class is (t_optprop_LUT_1_2)
          stop 'Twostream LUT generation is broken at the  moment -- need to write loadbalancer stuff'
          call prepare_table_space(LUT%S,2_iintegers)

        class is (t_optprop_LUT_3_6)
          call prepare_table_space(LUT%S,OPP%diff_streams**2)

        class is (t_optprop_LUT_wedge_5_8)
          call prepare_table_space(LUT%S,OPP%diff_streams**2)

        class default
          stop 'dont know optprop_class'
      end select
    endif

    if(myid.le.0 .and. comm_size.le.1) stop 'At the moment creation of diffuse Lookuptable needs at least two mpi-ranks to work... please run with more ranks.'
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
                isrc = modulo(cnt-1, Nsrc) +1
                lutindex = (cnt-1) / Nsrc +1
                call mpi_send(lutindex, 1_mpiint, imp_int, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(isrc, 1_mpiint, imp_int, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)
                call mpi_send(present(T), 1_mpiint, imp_logical, status(MPI_SOURCE), WORKMSG, comm, mpierr); call CHKERR(mpierr)

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

              if( mod(workindex-1, total_size/10).eq.0 ) then !every 10 percent of LUT dump it.
                print *,'Writing diffuse table to file...',workindex,total_size
                call ncwrite(S%table_name_c  , S%c         ,iierr)
                call ncwrite(S%table_name_tol, S%stddev_tol,iierr)
                print *,'done writing!',iierr
              endif

            case(FINALIZEMSG)
              call mpi_recv(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
              finalizedworkers = finalizedworkers+1
              if(finalizedworkers.eq.comm_size-1) then

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

    integer(mpiint) :: comm_size, myid

    call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)
    call MPI_Comm_size(comm, comm_size, mpierr); call CHKERR(mpierr)

    if(myid.eq.0) &
      call prepare_table_space(LUT%S(iphi,itheta), OPP%diff_streams, LUT%T(iphi,itheta), OPP%dir_streams)

    if(myid.le.0 .and. comm_size.le.1) &
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
              call mpi_recv(lutindex, 1_mpiint, imp_int, status(MPI_SOURCE), HAVERESULTSMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(isrc, 1_mpiint, imp_int, status(MPI_SOURCE), HAVERESULTSMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_diff, size(S_diff), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(S_tol , size(S_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_dir , size(T_dir ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)
              call mpi_recv(T_tol , size(T_tol ), imp_real, status(MPI_SOURCE), RESULTMSG, comm, status, mpierr); call CHKERR(mpierr)

              ! Sort coefficients into destination ordering and put em in LUT
              do idst = 1, OPP%diff_streams
                ind = (idst-1) * Nsrc + isrc
                S%c         (ind, lutindex) = S_diff(idst)
                S%stddev_tol(ind, lutindex) = S_tol (idst)
              enddo
              if(present(T)) then
                do idst = 1, OPP%dir_streams
                  ind = (idst-1)*OPP%dir_streams + isrc
                  T%c         (ind, lutindex) = T_dir (idst)
                  T%stddev_tol(ind, lutindex) = T_tol (idst)
                enddo
              endif

              !do idst = 1, Nsrc
              !  print *, myid, 'S%c for isrc', isrc, 'idst', idst, S_diff(idst)
              !enddo

              if( mod(lutindex*(Nsrc-1)+isrc-1, total_size/100).eq.0 ) & !every 1 percent report status
                  print *,'Calculated LUT...',(100*(lutindex*(Nsrc-1)+isrc-1))/total_size,'%'

              if( mod(lutindex*(Nsrc-1)+isrc, total_size/3 ).eq.0 ) then !every 30 percent of LUT dump it.
                print *,'Writing table to file...', S%table_name_c
                call ncwrite(S%table_name_c  , S%c         ,iierr)
                print *,'Writing table to file...', S%table_name_tol
                call ncwrite(S%table_name_tol, S%stddev_tol,iierr)
                if(present(T)) then
                  print *,'Writing table to file...', T%table_name_c
                  call ncwrite(T%table_name_c  , T%c         ,iierr)
                  print *,'Writing table to file...', T%table_name_tol
                  call ncwrite(T%table_name_tol, T%stddev_tol,iierr)
                endif
                print *,'done writing!',iierr
              endif

            case(FINALIZEMSG)
              call mpi_recv(idummy, 1_mpiint, imp_int, status(MPI_SOURCE), FINALIZEMSG, comm, status, mpierr); call CHKERR(mpierr)
              finalizedworkers = finalizedworkers+1
              if(finalizedworkers.eq.comm_size-1) exit ! all work is done

            end select
          endif
        enddo

        print *,'Writing table to file...'
        call ncwrite(S%table_name_c  , S%c         ,iierr)
        call ncwrite(S%table_name_tol, S%stddev_tol,iierr)
        if(present(T)) then
          call ncwrite(T%table_name_c  , T%c         ,iierr)
          call ncwrite(T%table_name_tol, T%stddev_tol,iierr)
          print *,'done writing!',iierr,':: max_atol S',maxval(S%stddev_tol),'max_atol T',maxval(T%stddev_tol)
        else
          print *,'done writing!',iierr,':: max_atol S',maxval(S%stddev_tol)
        endif
      end subroutine
      subroutine worker(config)
          type(t_lut_config), intent(in) :: config
          integer(iintegers) :: isrc
          logical :: ldir

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
                call mpi_recv( lutindex, 1_mpiint, imp_int, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_recv( isrc, 1_mpiint, imp_int, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)
                call mpi_recv( ldir, 1_mpiint, imp_logical, 0_mpiint, WORKMSG, comm, status, mpierr); call CHKERR(mpierr)

                call OPP%LUT_bmc_wrapper(config, lutindex, isrc, ldir, &
                  mpi_comm_self, S_diff, T_dir, S_tol, T_tol)

                !print *,'Computed isrc',isrc,'aspect',aspect_zx,'tau',tau_z, w0, g,':', phi, theta
                !print *,myid,'Computed values for ',lutindex, isrc, ldir

                call mpi_send(lutindex , 1_mpiint     , imp_int  , status(MPI_SOURCE) , HAVERESULTSMSG , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(isrc , 1_mpiint     , imp_int  , status(MPI_SOURCE) , HAVERESULTSMSG , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(S_diff    , size(S_diff) , imp_real , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(S_tol     , size(S_tol ) , imp_real , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
                call mpi_send(T_dir     , size(T_dir ) , imp_real , status(MPI_SOURCE) , RESULTMSG      , comm , mpierr); call CHKERR(mpierr)
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
end subroutine createLUT

subroutine prepare_table_space(OPP, config, S, T)
  class(t_optprop_LUT) :: OPP
  type(t_lut_config), intent(in) :: config
  type(t_table),intent(inout) :: S
  type(t_table),intent(inout), optional :: T

  integer(iintegers) :: errcnt

  print *,'Allocating Space for LUTs'
  errcnt = 0
  if(present(T)) then
    if(.not.associated(S%c         )) allocate(S%c         (OPP%diff_streams*OPP%dir_streams, product(config%dims(:)%N)), source=nil)
    if(.not.allocated (S%stddev_tol)) allocate(S%stddev_tol(OPP%diff_streams*OPP%dir_streams, product(config%dims(:)%N)), source=1e8_ireals)
    if(.not.associated(T%c         )) allocate(T%c         (OPP%dir_streams**2, product(config%dims(:)%N)), source=nil)
    if(.not.allocated (T%stddev_tol)) allocate(T%stddev_tol(OPP%dir_streams**2, product(config%dims(:)%N)), source=1e8_ireals)
  else
    if(.not.associated(S%c         )) allocate(S%c         (OPP%diff_streams**2, product(config%dims(:)%N)), source=nil)
    if(.not.allocated (S%stddev_tol)) allocate(S%stddev_tol(OPP%diff_streams**2, product(config%dims(:)%N)), source=1e8_ireals)
  endif
end subroutine

! return the integer in config%dims that corresponds to the given dimension
function find_lut_dim_by_name(config, dimname) result(kdim)
  type(t_lut_config), intent(in) :: config
  character(len=*), intent(in) :: dimname
  integer(iintegers) :: kdim

  integer(iintegers) :: k
  do k=1,size(config%dims)
    if(trim(dimname).eq.trim(config%dims(k)%dimname)) then
      kdim = k
      return
    endif
  enddo
  kdim=-1
end function

subroutine get_sample_pnt_by_name_and_index(config, dimname, index_1d, sample_pnt, ierr)
  type(t_lut_config), intent(in) :: config
  character(len=*), intent(in) :: dimname
  integer(iintegers), intent(in) :: index_1d
  real(ireals), intent(out) :: sample_pnt
  integer(mpiint), intent(out) :: ierr

  integer(iintegers) :: kdim, nd_indices(size(config%dims))

  nd_indices = ind_1d_to_nd(config%offsets, index_1d)

  kdim = find_lut_dim_by_name(config, trim(dimname))
  if(kdim.lt.i1) then ! could not find the corresponding dimension
    ierr = 1
    sample_pnt = nil
    return
  endif
  if(nd_indices(kdim).gt.size(config%dims(kdim)%v)) then
    print *,index_1d,'nd_indices', nd_indices
    call CHKERR(1_mpiint, 'wrong indices in kdim')
  endif
  sample_pnt = config%dims(kdim)%v(nd_indices(kdim))
  ierr = 0
end subroutine

subroutine LUT_bmc_wrapper(OPP, config, index_1d, src, dir, comm, S_diff, T_dir, S_tol, T_tol)
    class(t_optprop_LUT) :: OPP
    type(t_lut_config), intent(in) :: config
    integer(iintegers), intent(in) :: index_1d
    integer(iintegers), intent(in) :: src
    logical, intent(in) :: dir
    integer(mpiint), intent(in) :: comm

    real(ireals),intent(out) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(ireals),intent(out) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)

    real(ireals) :: aspect_zx, aspect_zy, tauz, w0, g, phi, theta
    integer(mpiint) :: ierr

    call get_sample_pnt_by_name_and_index(config, 'aspect_zx', index_1d, aspect_zx, ierr); call CHKERR(ierr, 'aspect_zx has to be present')
    call get_sample_pnt_by_name_and_index(config, 'tau', index_1d, tauz, ierr); call CHKERR(ierr, 'tauz has to be present')
    call get_sample_pnt_by_name_and_index(config, 'w0', index_1d, w0, ierr); call CHKERR(ierr, 'w0 has to be present')
    call get_sample_pnt_by_name_and_index(config, 'g', index_1d, g, ierr); call CHKERR(ierr, 'g has to be present')
    if(dir) then
      call get_sample_pnt_by_name_and_index(config, 'phi', index_1d, phi, ierr); call CHKERR(ierr, 'phi has to be present for direct calculations')
      call get_sample_pnt_by_name_and_index(config, 'theta', index_1d, theta, ierr); call CHKERR(ierr, 'theta has to be present for direct calculations')
    endif

    call get_sample_pnt_by_name_and_index(config, 'aspect_zy', index_1d, aspect_zy, ierr)
    if(ierr.ne.0) then
      aspect_zy = aspect_zx ! set dy = dy
    endif

    call bmc_wrapper(OPP, src, aspect_zx, aspect_zy, tauz, w0, g, dir, phi, theta, comm, S_diff, T_dir, S_tol, T_tol)

end subroutine

subroutine bmc_wrapper(OPP, src, aspect_zx, aspect_zy, tauz, w0, g, dir, phi, theta, comm, S_diff, T_dir, S_tol, T_tol)
    class(t_optprop_LUT) :: OPP
    integer(iintegers),intent(in) :: src
    logical,intent(in) :: dir
    integer(mpiint),intent(in) :: comm
    real(ireals), intent(in) :: aspect_zx, aspect_zy, tauz, w0, g, phi, theta
    real(ireals) :: dx,dy

    real(ireals),intent(out) :: S_diff(OPP%diff_streams),T_dir(OPP%dir_streams)
    real(ireals),intent(out) :: S_tol (OPP%diff_streams),T_tol(OPP%dir_streams)

    real(ireals) :: bg(3)
    real(ireals), parameter :: dz = 100

    dx = dz / aspect_zx
    dy = dz / aspect_zy

    bg(1) = tauz / dz * (one-w0)
    bg(2) = tauz / dz * w0
    bg(3) = g

    S_diff=nil
    T_dir=nil

    !print *,comm,'BMC :: calling bmc_get_coeff',bg,'src',src,'phi/theta',phi,theta,dz
    call OPP%bmc%get_coeff(comm, bg, src, &
      dir, phi, theta, &
      dx, dy, dz, &
      S_diff, T_dir, S_tol, T_tol, &
      inp_atol=stddev_atol-epsilon(stddev_atol)*10, &
      inp_rtol=stddev_rtol-epsilon(stddev_rtol)*10 )
    !print *,'BMC :: dir',T_dir,'diff',S_diff
end subroutine
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

subroutine populate_LUT_dim(dimname, N, lut_dim, vrange, preset)
  character(len=*),intent(in) :: dimname
  integer(iintegers), intent(in) :: N
  type(t_LUT_dim),intent(out) :: lut_dim
  real(ireals), optional :: vrange(:), preset(:)
  integer(iintegers) :: k
  if(allocated(lut_dim%v)) return ! already done
  allocate(lut_dim%v(N))

  if(present(vrange)) then
    do k=1,N
      lut_dim%v(k) = lin_index_to_param(one*k, vrange, N)
    enddo
  elseif(present(preset)) then
    if(size(preset).ne.N) &
      call CHKERR(1_mpiint, 'Given preset size does not conform to proposed size N '//itoa(N)//' vs '//itoa(size(preset)))
    lut_dim%v = preset
  else
    call CHKERR(1_mpiint, 'Have to provide either a number of a preset for LUT dimension')
  endif

  lut_dim%vrange = [lut_dim%v(1), lut_dim%v(N)]
  lut_dim%dimname = trim(dimname)
  lut_dim%N = size(lut_dim%v)
end subroutine

subroutine set_parameter_space(OPP)
    class(t_optprop_LUT) :: OPP
    allocate(OPP%dirconfig)
    allocate(OPP%diffconfig)

    select type(OPP)
      !class is (t_optprop_LUT_1_2)
      !    OPP%Nphi   = 1 ! azimithally average in 1D
      !    OPP%interp_mode = interp_mode_1_2
      class is (t_optprop_LUT_8_10)
          OPP%interp_mode = interp_mode_8_10
          allocate(OPP%dirconfig%dims(6))
          call populate_LUT_dim('tau',       Ntau, OPP%dirconfig%dims(1), preset=preset_tau)
          call populate_LUT_dim('w0',        Nw0, OPP%dirconfig%dims(2), preset=preset_w0)
          call populate_LUT_dim('g',         Ng, OPP%dirconfig%dims(3), preset=preset_g)
          call populate_LUT_dim('aspect_zx', Naspect, OPP%dirconfig%dims(4), preset=preset_aspect)
          call populate_LUT_dim('phi',       Nphi, OPP%dirconfig%dims(5), vrange=real([0,90], ireals))
          call populate_LUT_dim('theta',     Ntheta, OPP%dirconfig%dims(6), vrange=real([0,90], ireals))
          !call populate_LUT_dim('theta', Ntheta, OPP%dirconfig%dims(6), preset=preset_theta)
          allocate(OPP%diffconfig%dims(4))
          call populate_LUT_dim('tau',       Ntau, OPP%diffconfig%dims(1), preset=preset_tau)
          call populate_LUT_dim('w0',        Nw0, OPP%diffconfig%dims(2), preset=preset_w0)
          call populate_LUT_dim('g',         Ng, OPP%diffconfig%dims(3), preset=preset_g)
          call populate_LUT_dim('aspect_zx', Naspect, OPP%diffconfig%dims(4), preset=preset_aspect)

      !class is (t_optprop_LUT_3_10)
      !    OPP%interp_mode = interp_mode_3_10
      !class is (t_optprop_LUT_3_6)
      !    OPP%interp_mode = interp_mode_3_6
      !class is (t_optprop_LUT_wedge_5_8)
      !    OPP%interp_mode = interp_mode_wedge_5_8
      !    ps%range_phi = [-70, 70]
      class default
        stop 'set_parameter space: unexpected type for optprop_LUT object!'
    end select

    ! Determine offsets
    allocate(OPP%dirconfig%offsets(size(OPP%dirconfig%dims)))
    OPP%dirconfig%offsets = ndarray_offsets(OPP%dirconfig%dims(:)%N)

    allocate(OPP%diffconfig%offsets(size(OPP%diffconfig%dims)))
    OPP%diffconfig%offsets = ndarray_offsets(OPP%diffconfig%dims(:)%N)

    !if(ldebug) then
    !  print *,'set_parameter space dims:', size(OPP%diffconfig%dims), size(OPP%dirconfig%dims)
    !  do k=1,size(OPP%dirconfig%dims)
    !    print *,'dim ',trim(OPP%dirconfig%dims(k)%dimname), OPP%dirconfig%offsets(k), OPP%dirconfig%dims(k)%vrange, ':', OPP%dirconfig%dims(k)%v
    !  enddo
    !endif
end subroutine

  subroutine scatter_LUTtables(OPP, comm)
      use m_optprop_parameters, only: luse_memory_map
      integer(mpiint) ,intent(in) :: comm
      class(t_optprop_LUT) :: OPP

      integer(mpiint) :: myid, ierr
      real(ireals), pointer :: mmap_ptr(:,:)=>NULL()


      call MPI_Comm_rank(comm, myid, mpierr); call CHKERR(mpierr)

      if (luse_memory_map) then
        call arr_to_mmap(comm, trim(OPP%Sdiff%table_name_c(1))//'.Sdiff.mmap', OPP%Sdiff%c, mmap_ptr, ierr)
        if(associated(OPP%Sdiff%c)) deallocate(OPP%Sdiff%c)
        OPP%Sdiff%c => mmap_ptr

        call arr_to_mmap(comm, trim(OPP%Sdir%table_name_c(1))//'.Sdir.mmap', OPP%Sdir%c, mmap_ptr, ierr)
        if(associated(OPP%Sdir%c)) deallocate(OPP%Sdir%c)
        OPP%Sdir%c => mmap_ptr

        call arr_to_mmap(comm, trim(OPP%Tdir%table_name_c(1))//'.Tdir.mmap', OPP%Tdir%c, mmap_ptr, ierr)
        if(associated(OPP%Tdir%c)) deallocate(OPP%Tdir%c)
        OPP%Tdir%c => mmap_ptr

      else
        if( mpi_logical_or(comm, .not.associated(OPP%Sdir%c) )) &
          call imp_bcast(comm, OPP%Sdir%c, 0_mpiint)  ! DIRECT 2 DIRECT

        if( mpi_logical_or(comm, .not.associated(OPP%Tdir%c) )) &
          call imp_bcast(comm, OPP%Tdir%c, 0_mpiint)  ! DIRECT 2 DIFFUSE

        if( mpi_logical_or(comm, .not.associated(OPP%Sdiff%c) )) &
          call imp_bcast(comm, OPP%Sdiff%c, 0_mpiint)

      endif
  end subroutine

subroutine LUT_get_dir2dir(OPP, sample_pts, C)
    class(t_optprop_LUT) :: OPP
    real(ireals),intent(in) :: sample_pts(:)
    real(ireals),intent(out):: C(:) ! dimension(OPP%dir_streams**2)

    integer(iintegers) :: src, kdim, ind1d
    real(ireals) :: pti(size(sample_pts)), norm

    do kdim = 1, size(sample_pts)
      pti(kdim) = search_sorted_bisection(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti))
      C = OPP%Tdir%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti, OPP%Tdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(OPP%interp_mode)//' not implemented yet! please choose something else!')
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
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      endif
    endif
    !call CHKERR(1_mpiint, 'DEBUG')
end subroutine

subroutine LUT_get_dir2diff(OPP, sample_pts, C)
    class(t_optprop_LUT) :: OPP
    real(ireals),intent(in) :: sample_pts(:)
    real(ireals),intent(out):: C(:) ! dimension(OPP%dir_streams*OPP%diff_streams)

    integer(iintegers) :: src, kdim, ind1d
    real(ireals) :: pti(size(sample_pts)), norm

    do kdim = 1, size(sample_pts)
      pti(kdim) = search_sorted_bisection(OPP%dirconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%dirconfig%offsets, nint(pti))
      C = OPP%Sdir%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti, OPP%Sdir%c, OPP%dirconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(OPP%interp_mode)//' not implemented yet! please choose something else!')
    end select

    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%diff_streams
        norm = sum( C( src:size(C):OPP%dir_streams ) )
        if(real(norm).gt.one+1e-5_ireals) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        do src=1,OPP%diff_streams
          print *,'SUM dir2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%dir_streams)),' :: coeff',C( src:size(C):OPP%dir_streams )
        enddo
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      endif
    endif
end subroutine

subroutine LUT_get_diff2diff(OPP, sample_pts, C)
    class(t_optprop_LUT) :: OPP
    real(ireals),intent(in) :: sample_pts(:)
    real(ireals),intent(out):: C(:) ! dimension(OPP%diff_streams**2)

    integer(iintegers) :: src, kdim, ind1d
    real(ireals) :: pti(size(sample_pts)), norm

    do kdim = 1, size(sample_pts)
      pti(kdim) = search_sorted_bisection(OPP%diffconfig%dims(kdim)%v, sample_pts(kdim))
    enddo

    select case(OPP%interp_mode)
    case(1)
      ! Nearest neighbour
      ind1d = ind_nd_to_1d(OPP%diffconfig%offsets, nint(pti))
      C = OPP%Sdiff%c(:, ind1d)
    case(2)
      call interp_vec_simplex_nd(pti, OPP%Sdiff%c, OPP%diffconfig%offsets, C)
    case default
      call CHKERR(1_mpiint, 'interpolation mode '//itoa(OPP%interp_mode)//' not implemented yet! please choose something else!')
    end select

    if(ldebug_optprop) then
      !Check for energy conservation:
      iierr=0
      do src=1,OPP%diff_streams
        norm = sum( C( src:size(C):OPP%diff_streams ) )
        if(norm.gt.one+1e-5_ireals) iierr=iierr+1
      enddo
      if(iierr.ne.0) then
        do src=1,OPP%diff_streams
          print *,'SUM diff2diff coeff for src ',src,' :: sum ',sum(C( src:size(C):OPP%diff_streams)),' :: coeff',C(src:size(C):OPP%diff_streams)
        enddo
        call CHKERR(1_mpiint, 'Check for energy conservation failed')
      endif
    endif
end subroutine

!subroutine interp_4p2d(pti,ctable,C)
!        integer,parameter :: Ndim=6
!        real(ireals),intent(in) :: pti(Ndim)
!        type(table),intent(in) :: ctable(:,:)  ! contains Nphi, Ntheta databases
!        real(ireals),intent(out) :: C(:)
!
!        real(ireals) :: weights(Ndim)
!        integer(iintegers) :: indices(2,2),fpti(Ndim)
!
!        ! Instead of doing a full interpolation in 6 dimension we start out with
!        ! 4 dimensions only at the cornerstones of the 4d hypercube
!        real(ireals) :: C4(size(C),6)
!
!        ! First determine the array indices, where to look.
!        fpti = floor(pti)
!        weights = modulo(pti, one)
!
!        indices(:,1) = max(i1, min( ubound(ctable,1, kind=iintegers), [i0,i1] +fpti(5) ) )
!        indices(:,2) = max(i1, min( ubound(ctable,2, kind=iintegers), [i0,i1] +fpti(6) ) )
!
!        call interp_4d( pti(1:4), ctable(indices(1,1), indices(1,2) )%c, C4(:,1) ) ! differing azimuth
!        call interp_4d( pti(1:4), ctable(indices(2,1), indices(1,2) )%c, C4(:,2) ) !        "
!        call interp_4d( pti(1:4), ctable(indices(1,1), indices(2,2) )%c, C4(:,3) )
!        call interp_4d( pti(1:4), ctable(indices(2,1), indices(2,2) )%c, C4(:,4) )
!
!        C4(:,5) = C4(:,1) + weights(5) * ( C4(:,2) - C4(:,1) )
!        C4(:,6) = C4(:,3) + weights(5) * ( C4(:,4) - C4(:,3) )
!        C       = C4(:,5) + weights(6) * ( C4(:,6) - C4(:,5) )
!end subroutine
!subroutine interp_4p1d(pti,ctable,C)
!        integer,parameter :: Ndim=5
!        real(ireals),intent(in) :: pti(Ndim)
!        type(table),intent(in) :: ctable(:)  ! contains N databases
!        real(ireals),intent(out) :: C(:)
!
!        real(ireals) :: weights(Ndim)
!        integer(iintegers) :: indices(2),fpti(Ndim)
!
!        ! Instead of doing a full interpolation in 6 dimension we start out with
!        ! 4 dimensions only at the cornerstones of the 4d hypercube
!        real(ireals) :: C4(size(C),2)
!
!        ! First determine the array indices, where to look.
!        fpti = floor(pti)
!        weights = modulo(pti, one)
!
!        indices(:) = max(i1, min( ubound(ctable,1, kind=iintegers), [i0,i1] +fpti(5)))
!
!        call interp_4d( pti(1:4), ctable(indices(1))%c, C4(:,1) ) ! differing zenith
!        call interp_4d( pti(1:4), ctable(indices(2))%c, C4(:,2) ) !        "
!
!        C = C4(:,1) + weights(5) * ( C4(:,2) - C4(:,1) )
!end subroutine
!
!function get_indices_4d(aspect, tauz, w0, g, ps)
!    real(ireals) :: get_indices_4d(4)
!    real(ireals),intent(in) :: aspect, tauz, w0, g
!    type(parameter_space),intent(in) :: ps
!
!    get_indices_4d(1) = search_sorted_bisection(ps%aspect, aspect)
!    get_indices_4d(2) = search_sorted_bisection(ps%tau   , tauz)
!    get_indices_4d(3) = search_sorted_bisection(ps%w0    , w0)
!    get_indices_4d(4) = search_sorted_bisection(ps%g     , g)
!end function
!function get_indices_6d(aspect, tauz, w0, g, phi, theta, ps)
!    real(ireals) :: get_indices_6d(6)
!    real(ireals),intent(in) :: aspect, tauz, w0, g, phi, theta
!    type(parameter_space),intent(in) :: ps
!
!    get_indices_6d(1:4) = get_indices_4d(aspect, tauz, w0, g, ps)
!
!    get_indices_6d(5) = search_sorted_bisection(ps%phi  ,phi )
!    get_indices_6d(6) = search_sorted_bisection(ps%theta,theta)
!end function
!
!logical function valid_input(val,range)
!    real(ireals),intent(in) :: val,range(2)
!    if(val.lt.range(1) .or. val.gt.range(2) ) then
!      valid_input=.False.
!      print *,'ohoh, this val is not in the optprop database range!',val,'not in',range
!    else
!      valid_input=.True.
!    endif
!end function

!subroutine catch_limits(ps, aspect, tauz, w0, g)
!    type(parameter_space),intent(in) :: ps
!    real(ireals),intent(in) :: aspect, tauz, w0, g
!
!    iierr=0
!
!    if( aspect.lt.ps%range_aspect(1) .or. aspect.gt.ps%range_aspect(2) ) then
!      print *,'aspect ratio is not in LookUpTable Range',aspect, 'LUT range',ps%range_aspect
!      iierr=iierr+1
!    endif
!    if( tauz.lt.ps%range_tau(1) .or. tauz.gt.ps%range_tau(2) ) then
!      print *,'tau is not in LookUpTable Range',tauz, 'LUT range',ps%range_tau
!      iierr=iierr+1
!    endif
!    if( w0.lt.ps%range_w0(1) .or. w0.gt.ps%range_w0(2) ) then
!      print *,'w0 is not in LookUpTable Range',w0, 'LUT range',ps%range_w0
!      iierr=iierr+1
!    endif
!    if( g.lt.ps%range_g(1) .or. g.gt.ps%range_g(2) ) then
!      print *,'g is not in LookUpTable Range',g, 'LUT range',ps%range_g
!      iierr=iierr+1
!    endif
!    if(iierr.ne.0) print*, 'The LookUpTable was asked to give a coefficient, it was not defined for. Please specify a broader range.',iierr
!end subroutine

end module
