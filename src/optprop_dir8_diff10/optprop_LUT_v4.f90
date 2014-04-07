module tenstream_optprop_LUT_8_10
  use helper_functions, only : approx
  use data_parameters, only : ireals, iintegers, one,zero,i0,i1,i3,mpiint,nil
  use boxmc_parameters_8_10, only: dir_streams,diff_streams, Ndz,Nkabs,Nksca,Ng,Nphi,Ntheta,interp_mode,delta_scale
  use boxmc_oop, only: boxmc_8_10
  use tenstream_interpolation, only: interp_4d,interp_6d,interp_6d_recursive,interp_4p2d
  use arrayio

  use mpi!, only: MPI_Comm_rank,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_Bcast


  implicit none

  private
  public :: init_LUT,LUT_get_dir2dir,LUT_get_dir2diff,LUT_get_diff2diff,delta_scale,bmc_wrapper
  ! This module loads and generates the LUT-tables for Tenstream Radiation
  ! computations.
  ! It also holds functions for interpolation on the regular LUT grid.

  logical :: LUT_initiliazed=.False.,optprop_debug=.True.

  integer(iintegers) :: iierr
  integer(mpiint) :: ierr

  type parameter_space
    real(ireals),dimension(Ndz)    :: dz
    real(ireals),dimension(Nkabs ) :: kabs
    real(ireals),dimension(Nksca)  :: ksca
    real(ireals),dimension(Ng)     :: g
    real(ireals),dimension(Nphi)   :: phi
    real(ireals),dimension(Ntheta) :: theta
    real(ireals) :: dz_exponent,kabs_exponent,ksca_exponent,g_exponent
    real(ireals),dimension(2)      :: range_dz      = [ 50._ireals   , 5001._ireals ]
    real(ireals),dimension(2)      :: range_kabs    = [ 1e-99_ireals, 10._ireals    ] !lower limit for kabs,ksca is set in set_parameter_space
    real(ireals),dimension(2)      :: range_ksca    = [ 1e-99_ireals, one           ] !lower limit for kabs,ksca is set in set_parameter_space
    real(ireals),dimension(2)      :: range_g       = [ zero        , .999_ireals   ]
    real(ireals),dimension(2)      :: range_phi     = [ zero        , 90._ireals    ]
    real(ireals),dimension(2)      :: range_theta   = [ zero        , 90._ireals    ]
  end type

  type table
    real(ireals),allocatable :: c(:,:,:,:,:)
  end type

  type diffuseTable
    real(ireals) :: dx,dy
    character(len=300) :: fname
    type(table) :: S
    type(parameter_space) :: pspace
  end type

  type directTable
    real(ireals) :: dx,dy
    character(len=300) :: fname
    type(table),allocatable :: S(:,:) !dim=(phi,theta)
    type(table),allocatable :: T(:,:) !dim=(phi,theta)
    type(parameter_space) :: pspace
  end type

  type(directTable) :: dirLUT
  type(diffuseTable) :: diffLUT

  type(boxmc_8_10) :: bmc

contains
  subroutine scatter_LUTtables(azis,szas,comm)
      real(ireals),intent(in) :: szas(:),azis(:) 
      integer(mpiint) ,intent(in) :: comm
      integer(mpiint) :: myid,Ntot,Ncoeff
      integer(iintegers) :: iphi,itheta
      real(ireals),allocatable,dimension(:) :: tmp
      logical :: angle_mask(Nphi,Ntheta)

      call MPI_Comm_rank(comm, myid, ierr)

      call determine_angles_to_load(dirLUT, azis, szas, angle_mask,comm)

      do itheta=1,Ntheta
        do iphi  =1,Nphi

!          if(dirLUT%pspace%theta(itheta).le.1e-3_ireals .and. dirLUT%pspace%phi(iphi).gt.1e-3_ireals ) cycle ! dont need to calculate different azimuth angles except the zero one... rest is symmetric
!          if(itheta.gt.2.or.iphi.gt.2) cycle !TODO just a shortcut, so that we can already calculate stuff until we wait for the LUT calculations
          if(.not.angle_mask(iphi,itheta) ) cycle

          ! DIRECT 2 DIRECT
          if(myid.eq.0) then
            Ncoeff = size(dirLUT%T(iphi,itheta)%c, 1)
            Ntot   = size(dirLUT%T(iphi,itheta)%c) 
            print *,myid,'Scattering LUT tables....',Ncoeff,Ntot,' iphi,itheta',iphi,itheta 
          endif
          call mpi_bcast(Ncoeff, 1, MPI_INTEGER, 0, comm, ierr)
          call mpi_bcast(Ntot  , 1, MPI_INTEGER, 0, comm, ierr)
          allocate(tmp(Ntot)) 
          if(myid.eq.0) tmp = reshape(dirLUT%T(iphi,itheta)%c,shape(tmp) )

          call mpi_bcast(tmp, Ntot, MPI_DOUBLE_PRECISION, 0, comm, ierr)
          if(myid.gt.0) then
            allocate(dirLUT%T(iphi,itheta)%c(Ncoeff, Ndz, Nkabs, Nksca, Ng) )
            dirLUT%T(iphi,itheta)%c = reshape(tmp, [i1*Ncoeff, Ndz, Nkabs, Nksca, Ng] )
          endif
          deallocate(tmp)

          ! DIRECT 2 DIFFUSE
          if(myid.eq.0) then
            Ncoeff = size(dirLUT%S(iphi,itheta)%c, 1)
            Ntot   = size(dirLUT%S(iphi,itheta)%c) 
            !          print *,'Scattering LUT tables....',Ncoeff,Ntot
          endif
          call mpi_bcast(Ncoeff, 1, MPI_INTEGER, 0, comm, ierr)
          call mpi_bcast(Ntot  , 1, MPI_INTEGER, 0, comm, ierr)
          allocate(tmp(Ntot)) 
          if(myid.eq.0) tmp = reshape(dirLUT%S(iphi,itheta)%c,[Ntot] )

          call mpi_bcast(tmp, Ntot, MPI_DOUBLE_PRECISION, 0, comm, ierr)
          if(myid.gt.0) then
            allocate(dirLUT%S(iphi,itheta)%c(Ncoeff, Ndz, Nkabs, Nksca, Ng) )
            dirLUT%S(iphi,itheta)%c = reshape(tmp, [i1*Ncoeff, Ndz, Nkabs, Nksca, Ng] )
          endif
          deallocate(tmp)

        enddo
      enddo

      ! DIFFUSE 2 DIFFUSE
      if(myid.eq.0) then
        Ncoeff = size(diffLUT%S%c, 1)
        Ntot   = size(diffLUT%S%c) 
        !          print *,'Scattering LUT tables....',Ncoeff,Ntot
      endif
      call mpi_bcast(Ncoeff, 1, MPI_INTEGER, 0, comm, ierr)
      call mpi_bcast(Ntot  , 1, MPI_INTEGER, 0, comm, ierr)
      allocate(tmp(Ntot)) 
      if(myid.eq.0) tmp = reshape(diffLUT%S%c,[Ntot] )

      call mpi_bcast(tmp, Ntot, MPI_DOUBLE_PRECISION, 0, comm, ierr)
      if(myid.gt.0) then
        allocate(diffLUT%S%c(Ncoeff, Ndz, Nkabs, Nksca, Ng) )
        diffLUT%S%c = reshape(tmp, [i1*Ncoeff, Ndz, Nkabs, Nksca, Ng] )
      endif
      deallocate(tmp)

      print *,myid,' :: ',diffLUT%S%c(1,1,1,1,1)
!      call mpi_barrier(comm,ierr)
  end subroutine
  !{{{ init LUT
  subroutine init_LUT(dx,dy, azis,szas, comm)
      real(ireals),intent(in) :: dx,dy
      real(ireals),intent(in) :: szas(:),azis(:) ! all solar zenith angles that happen in this scene
      integer(mpiint) ,intent(in) :: comm
      integer(mpiint) :: myid,comm_size

      character(len=*),parameter :: lutbasename='/home/opt/cosmo_tica_lib/tenstream/optpropLUT/LUT_8_10.'
      character(len=300) :: descr
      integer(iintegers) :: idx,idy
      call MPI_Comm_rank(comm, myid, ierr)
      call MPI_Comm_size(comm, comm_size, ierr)
      idx = nint( dx/10  ) * 10
      idy = nint( dy/10  ) * 10

      if(myid.eq.0) print *,'Loading optical properties look-up tables... Tables have diffuse/direct',&
          Ndz*Nkabs*Nksca*Ng,'/',Ndz*Nkabs*Nksca*Ng*Nphi*Ntheta,'entries --> approx.', &
          (Ndz*Nkabs*Nksca*Ng*12+2*Ndz*Nkabs*Nksca*Ng*Nphi*Ntheta*39)*sizeof(one)/1024/1024,'MB'

      ! Load diffuse LUT
      write(descr,FMT='("diffuse.dx",I0,".pspace.dz",I0,".kabs",I0,".ksca",I0,".g",I0)') idx,Ndz,Nkabs,Nksca,Ng
      if(myid.eq.0) print *,'Loading diffuse LUT from ',descr
      diffLUT%fname = trim(lutbasename)//trim(descr)//'.h5'
      diffLUT%dx    = idx
      diffLUT%dy    = idy

      call set_parameter_space(diffLUT%pspace,diffLUT%dx)
      call loadLUT_diff(diffLUT,comm)

      ! Load direct LUTS
      write(descr,FMT='("direct.dx",I0,".pspace.dz",I0,".kabs",I0,".ksca",I0,".g",I0,".phi",I0,".theta",I0,".delta",L1)') idx,Ndz,Nkabs,Nksca,Ng,Nphi,Ntheta,delta_scale
      if(myid.eq.0) print *,'Loading direct LUT from ',descr
      dirLUT%fname = trim(lutbasename)//trim(descr)//'.h5'
      dirLUT%dx    = idx
      dirLUT%dy    = idy

      call set_parameter_space(dirLUT%pspace,dirLUT%dx)
      call loadLUT_dir(dirLUT, azis, szas, comm)

      if(comm_size.gt.1) call scatter_LUTtables(azis,szas,comm)

      call check_diffLUT_matches_pspace(diffLUT)
      call check_dirLUT_matches_pspace(dirLUT)

      LUT_initiliazed=.True.
      if(myid.eq.0) print *,'Done loading LUTs'

  end subroutine
  !}}}
  !{{{ check for pspace alignment
  subroutine check_diffLUT_matches_pspace(LUT)
      type(diffuseTable),intent(in) :: LUT
      real(ireals),allocatable :: buf(:)
      character(300) :: str(2)
      integer(iintegers) align(4);
      write(str(1),FMT='("dx",I0)')   int(LUT%dx)
      write(str(2),FMT='("dy",I0)')   int(LUT%dy)
      align=0
      call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","dz      "],buf,iierr) ; if(.not.all(approx( buf, LUT%pspace%dz  ,1e-6_ireals ) )) align(1)=1 ; deallocate(buf)
      call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","kabs    "],buf,iierr) ; if(.not.all(approx( buf, LUT%pspace%kabs,1e-6_ireals ) )) align(2)=1 ; deallocate(buf)
      call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","ksca    "],buf,iierr) ; if(.not.all(approx( buf, LUT%pspace%ksca,1e-6_ireals ) )) align(3)=1 ; deallocate(buf)
      call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","g       "],buf,iierr) ; if(.not.all(approx( buf, LUT%pspace%g   ,1e-6_ireals ) )) align(4)=1 ; deallocate(buf)

      if(any(align.ne.0)) then
        print *,'ERROR, the LUT we found at ',LUT%fname,'diffuse',str(1),str(2),' does not match the parameter space!',align
        call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","dz      "],buf,iierr)
        print *,'LUT :: dz ',LUT%pspace%dz,'::',buf ; deallocate(buf)
        call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","kabs      "],buf,iierr)
        print *,'LUT :: kabs ',LUT%pspace%kabs,'::',buf ; deallocate(buf)
        call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","ksca      "],buf,iierr)
        print *,'LUT :: ksca ',LUT%pspace%ksca,'::',buf ; deallocate(buf)
        call h5load([LUT%fname,'diffuse',str(1),str(2),"pspace","g         "],buf,iierr)
        print *,'LUT :: g ',LUT%pspace%g,'::',buf ; deallocate(buf)
        call exit
      endif
  end subroutine                                   
  subroutine check_dirLUT_matches_pspace(LUT)
      type(directTable),intent(in) :: LUT
      real(ireals),allocatable :: buf(:)
      character(300) :: str(2)
      integer(iintegers) align(6);
      write(str(1),FMT='("dx",I0)')   int(LUT%dx)
      write(str(2),FMT='("dy",I0)')   int(LUT%dy)
      align=0
      call h5load([LUT%fname,'direct',str(1),str(2),"pspace","dz      "],buf,iierr) ; if(.not.all(approx( buf,LUT%pspace%dz   ,1e-6_ireals) )) align(1)=1 ; deallocate(buf)
      call h5load([LUT%fname,'direct',str(1),str(2),"pspace","kabs    "],buf,iierr) ; if(.not.all(approx( buf,LUT%pspace%kabs ,1e-6_ireals) )) align(2)=1 ; deallocate(buf)
      call h5load([LUT%fname,'direct',str(1),str(2),"pspace","ksca    "],buf,iierr) ; if(.not.all(approx( buf,LUT%pspace%ksca ,1e-6_ireals) )) align(3)=1 ; deallocate(buf)
      call h5load([LUT%fname,'direct',str(1),str(2),"pspace","g       "],buf,iierr) ; if(.not.all(approx( buf,LUT%pspace%g    ,1e-6_ireals) )) align(4)=1 ; deallocate(buf)
      call h5load([LUT%fname,'direct',str(1),str(2),"pspace","phi     "],buf,iierr) ; if(.not.all(approx( buf,LUT%pspace%phi  ,1e-6_ireals) )) align(5)=1 ; deallocate(buf)
      call h5load([LUT%fname,'direct',str(1),str(2),"pspace","theta   "],buf,iierr) ; if(.not.all(approx( buf,LUT%pspace%theta,1e-6_ireals) )) align(6)=1 ; deallocate(buf)

      if(any(align.ne.0)) then
        print *,'ERROR, the LUT we found at ',LUT%fname,'direct',str(1),str(2),' does not match the parameter space!',align
        call h5load([LUT%fname,'direct',str(1),str(2),"pspace","dz      "],buf,iierr)
        print *,'LUT :: dz ',LUT%pspace%dz,'::',buf ; deallocate(buf)
        call h5load([LUT%fname,'direct',str(1),str(2),"pspace","kabs      "],buf,iierr)
        print *,'LUT :: kabs ',LUT%pspace%kabs,'::',buf ; deallocate(buf)
        call h5load([LUT%fname,'direct',str(1),str(2),"pspace","ksca      "],buf,iierr)
        print *,'LUT :: ksca ',LUT%pspace%ksca,'::',buf ; deallocate(buf)
        call h5load([LUT%fname,'direct',str(1),str(2),"pspace","g         "],buf,iierr)
        print *,'LUT :: g ',LUT%pspace%g,'::',buf ; deallocate(buf)
        call exit
      endif
  end subroutine                                   
!}}}
!{{{ load LUT                                          
subroutine loadLUT_diff(LUT,comm)
    type(diffuseTable) :: LUT
    integer(mpiint),intent(in) :: comm
    integer(iintegers) :: errcnt
    integer(mpiint) :: myid
    character(300) :: str(2)

    call MPI_Comm_rank(comm, myid, ierr)

    if(myid.eq.0) print *,'... loading diffuse LUT',myid
    if(allocated(LUT%S%c)) then
      print *,'LUT already loaded! is this a second call?',myid
      return
    endif

    write(str(1),FMT='("dx",I0)')   int(LUT%dx)
    write(str(2),FMT='("dy",I0)')   int(LUT%dy)

    errcnt=0
    if(myid.eq.0) then
      call h5load([LUT%fname,'diffuse',str(1),str(2),"S"],LUT%S%c,iierr) ; errcnt = errcnt+iierr
      if(allocated(LUT%S%c) ) then
        if( any(LUT%S%c.gt.one) .or. any(LUT%S%c.lt.zero) ) errcnt=100
        call check_diffLUT_matches_pspace(LUT)
      endif
    endif
    call mpi_bcast(errcnt,1 , MPI_INTEGER, 0, comm, ierr)

    if(errcnt.ne.0) then
      if(myid.eq.0) then
        print *,'Loading of diffuse tables failed for',trim(LUT%fname),'  diffuse ',trim(str(1)),' ',trim(str(2)),'::',errcnt
        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","range_dz   "],LUT%pspace%range_dz   ,iierr)
        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","range_kabs   "],LUT%pspace%range_kabs   ,iierr)
        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","range_ksca    "],LUT%pspace%range_ksca    ,iierr)
        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","range_g    "],LUT%pspace%range_g    ,iierr)

        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","dz   "],LUT%pspace%dz   ,iierr)
        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","kabs   "],LUT%pspace%kabs   ,iierr)
        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","ksca    "],LUT%pspace%ksca    ,iierr)
        call h5write([LUT%fname,'diffuse',str(1),str(2),"pspace","g    "],LUT%pspace%g    ,iierr)
      endif

      call createLUT_diff(LUT,[LUT%fname,'diffuse',str(1),str(2),"S"],comm)

      if(myid.eq.0) call h5write([LUT%fname,'diffuse',str(1),str(2),"S"],LUT%S%c,iierr)

      call mpi_barrier(comm,ierr)
      call exit() ! TODO: We exit here in order to split the jobs for shorter runtime.
    endif

    if(myid.eq.0) print *,'Done loading diffuse LUTs'
end subroutine

subroutine determine_angles_to_load(LUT,azis,szas, mask,comm)
    type(directTable) :: LUT
    real(ireals),intent(in) :: szas(:),azis(:) ! all solar zenith angles that happen in this scene
    integer(mpiint),intent(in) :: comm
    logical,intent(out) :: mask(Nphi,Ntheta) ! boolean array, which LUT entries should be loaded

    integer(iintegers) :: itheta, iphi
    integer(mpiint) :: myid
    logical :: lneed_azi, lneed_sza
    real(ireals) :: theta(2),phi(2) ! sza and azimuth angle

    call MPI_Comm_rank(comm, myid, ierr)
    mask = .False.
    ! Check if really need to load it... i.e. we only want to load angles which are necessary for this run.
    do itheta=1,Ntheta-1
      do iphi  =1,Nphi-1
        phi   = LUT%pspace%phi( [ iphi, iphi+1 ] )
        theta = LUT%pspace%theta( [ itheta, itheta+1 ]  )
        lneed_azi = any( azis.ge.phi(1) .and. azis.lt.phi(2) )
        lneed_sza = any( szas.ge.theta(1) .and. szas.lt.theta(2) )
        if( lneed_azi .and. lneed_sza ) mask([iphi,iphi+1],[itheta,itheta+1]) = .True.
      enddo
    enddo
    if(myid.eq.0) then
      print *,'       phis',LUT%pspace%range_phi
      do itheta=1,Ntheta
        print *,'theta=',LUT%pspace%theta(itheta),' :: ',mask(:,itheta)
      enddo
    endif
!    call mpi_barrier(comm,ierr)
end subroutine
subroutine loadLUT_dir(LUT, azis,szas, comm)
    real(ireals),intent(in) :: szas(:),azis(:) ! all solar zenith angles that happen in this scene
    type(directTable) :: LUT
    integer(mpiint),intent(in) :: comm
    integer(mpiint) :: myid
    integer(iintegers) :: errcnt,iphi,itheta
    character(300) :: str(4)
    logical :: angle_mask(Nphi,Ntheta)

    call MPI_Comm_rank(comm, myid, ierr)

    if(myid.eq.0) print *,'... loading direct LUT'
    if(allocated(LUT%S).and.allocated(LUT%T)) then
      print *,'LUTs already loaded!',myid
      return
    endif
    allocate( LUT%S(Nphi,Ntheta) )
    allocate( LUT%T(Nphi,Ntheta) )

    write(str(1),FMT='("dx",I0)') nint(LUT%dx)
    write(str(2),FMT='("dy",I0)') nint(LUT%dy)

    call determine_angles_to_load(LUT, azis, szas, angle_mask,comm)
    errcnt=0
    do itheta=1,Ntheta
      do iphi  =1,Nphi
!        if(LUT%pspace%theta(itheta).le.1e-3_ireals .and. LUT%pspace%phi(iphi).gt.1e-3_ireals ) cycle ! dont need to calculate different azimuth angles except the zero one... rest is symmetric
        if(.not.angle_mask(iphi,itheta) ) cycle

        write(str(3),FMT='("phi",I0)')  int(LUT%pspace%phi(iphi)    )
        write(str(4),FMT='("theta",I0)')int(LUT%pspace%theta(itheta))

        if(myid.eq.0) then
            call h5load([LUT%fname,'direct',str(1),str(2),str(3),str(4),"S"],LUT%S(iphi,itheta)%c,iierr) ; errcnt = errcnt+iierr
            call h5load([LUT%fname,'direct',str(1),str(2),str(3),str(4),"T"],LUT%T(iphi,itheta)%c,iierr) ; errcnt = errcnt+iierr


            if(allocated(LUT%S(iphi,itheta)%c) ) then
              if(any( LUT%S(iphi,itheta)%c.gt.one ).or.any(LUT%S(iphi,itheta)%c.lt.zero) ) errcnt=errcnt+100
            endif
            if(allocated(LUT%T(iphi,itheta)%c) ) then
              if(any( LUT%T(iphi,itheta)%c.gt.one ).or.any(LUT%T(iphi,itheta)%c.lt.zero) ) errcnt=errcnt+100
              call check_dirLUT_matches_pspace(LUT)
            endif

        endif

        call mpi_bcast(errcnt,1 , MPI_INTEGER, 0, comm, ierr)

        if(errcnt.ne.0) then
          if(myid.eq.0) then
            print *,'Loading of direct tables failed for',trim(LUT%fname),'  direct ',trim(str(1)),' ',trim(str(2)),' ',trim(str(3)),' ',trim(str(4)),'::',errcnt
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","range_dz   "],LUT%pspace%range_dz   ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","range_kabs "],LUT%pspace%range_kabs   ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","range_ksca "],LUT%pspace%range_ksca    ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","range_g    "],LUT%pspace%range_g    ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","range_phi  "],LUT%pspace%range_phi  ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","range_theta"],LUT%pspace%range_theta,iierr)

            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","dz   "],LUT%pspace%dz   ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","kabs "],LUT%pspace%kabs   ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","ksca "],LUT%pspace%ksca    ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","g    "],LUT%pspace%g    ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","phi  "],LUT%pspace%phi  ,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),"pspace","theta"],LUT%pspace%theta,iierr)
          endif

          call createLUT_dir(LUT,[LUT%fname,'direct',str(1),str(2),str(3),str(4),"T"],[LUT%fname,'direct',str(1),str(2),str(3),str(4),"S"],comm,iphi,itheta)

          if(myid.eq.0) then
            call h5write([LUT%fname,'direct',str(1),str(2),str(3),str(4),"S"],LUT%S(iphi,itheta)%c,iierr)
            call h5write([LUT%fname,'direct',str(1),str(2),str(3),str(4),"T"],LUT%T(iphi,itheta)%c,iierr)
          endif

          call mpi_barrier(comm,ierr)
          call exit() ! TODO: We exit here in order to split the jobs for shorter runtime.
        endif

      enddo
    enddo

    if(myid.eq.0) print *,'Done loading direct LUTs'
end subroutine
!}}}
!{{{ create LUT
subroutine createLUT_diff(LUT,coeff_table_name,comm)
    type(diffuseTable) :: LUT
    integer(mpiint),intent(in) :: comm
    character(300),intent(in) :: coeff_table_name(:)

    integer(mpiint) :: myid
    integer(iintegers) :: idz,ikabs ,iksca,ig,total_size,cnt
    real(ireals) :: S_diff(diff_streams),T_dir(dir_streams)
    logical :: ldone

    call MPI_Comm_rank(comm, myid, ierr)

    if(myid.eq.0.and. .not. allocated(LUT%S%c) ) then
      allocate(LUT%S%c(12, Ndz,Nkabs ,Nksca,Ng))
      LUT%S%c = nil
      call h5write(coeff_table_name,LUT%S%c,iierr)
    endif

    total_size = Ng*Nksca*Nkabs *Ndz
    cnt=0
    do ig    =1,Ng
      do iksca    =1,Nksca    
        do ikabs    =1,Nkabs    
          do idz   =1,Ndz   
            cnt=cnt+1
            ! Check if we already calculated the coefficients and inform the other nodes
            ldone=.False.
            if(myid.eq.0) then
              if( all( LUT%S%c( :, idz,ikabs ,iksca,ig).ge.zero).and.all( LUT%S%c( :, idz,ikabs ,iksca,ig).le.one) ) ldone = .True.
            endif
            call mpi_bcast(ldone,1 , MPI_LOGICAL, 0, comm, ierr)
            if(ldone) cycle

            if(myid.eq.0) print *,'diff dx',LUT%dx,'dz',LUT%pspace%dz(idz),' :: ',LUT%pspace%kabs (ikabs ),LUT%pspace%ksca(iksca),LUT%pspace%g(ig),'(',100*cnt/total_size,'%)'
            ! src=1
            call bmc_wrapper(i1,LUT%dx,LUT%dy,LUT%pspace%dz(idz),LUT%pspace%kabs (ikabs ),LUT%pspace%ksca(iksca),LUT%pspace%g(ig),.False.,delta_scale,zero,zero,comm,S_diff,T_dir)
            if(myid.eq.0) LUT%S%c( 1, idz,ikabs ,iksca,ig) = S_diff(1)
            if(myid.eq.0) LUT%S%c( 2, idz,ikabs ,iksca,ig) = S_diff(2)
            if(myid.eq.0) LUT%S%c( 3, idz,ikabs ,iksca,ig) = sum(S_diff([3,4,7, 8]) )/4
            if(myid.eq.0) LUT%S%c( 4, idz,ikabs ,iksca,ig) = sum(S_diff([5,6,9,10]) )/4
            ! src=3
            call bmc_wrapper(i3,LUT%dx,LUT%dy,LUT%pspace%dz(idz),LUT%pspace%kabs (ikabs ),LUT%pspace%ksca(iksca),LUT%pspace%g(ig),.False.,delta_scale,zero,zero,comm,S_diff,T_dir)
            if(myid.eq.0) LUT%S%c( 5:10, idz,ikabs ,iksca,ig) = S_diff(1:6)
            if(myid.eq.0) LUT%S%c( 11  , idz,ikabs ,iksca,ig) = sum(S_diff(7: 8))/2
            if(myid.eq.0) LUT%S%c( 12  , idz,ikabs ,iksca,ig) = sum(S_diff(9:10))/2

          enddo !dz
        enddo !kabs

        if(myid.eq.0) then
          print *,'Writing diffuse table to file...'
          call h5write(coeff_table_name,LUT%S%c,iierr)
          print *,'done writing!',iierr
        endif

      enddo !ksca
    enddo !g
    if(myid.eq.0) print *,'done calculating diffuse coefficients'
end subroutine
subroutine createLUT_dir(LUT,dir_coeff_table_name,diff_coeff_table_name,comm,iphi,itheta)
    type(directTable) :: LUT
    integer(iintegers),intent(in) :: iphi,itheta
    integer(mpiint),intent(in) :: comm
    character(300),intent(in) :: dir_coeff_table_name(:),diff_coeff_table_name(:)

    integer(mpiint) :: myid
    integer(iintegers) :: idz,ikabs ,iksca,ig,src,total_size,cnt
    real(ireals) :: S_diff(diff_streams),T_dir(dir_streams)
    logical :: ldone

    call MPI_Comm_rank(comm, myid, ierr)

    if(myid.eq.0) then
      print *,'calculating direct coefficients for ',iphi,itheta

      if(.not. allocated(LUT%S(iphi,itheta)%c) ) then
        allocate(LUT%S(iphi,itheta)%c(dir_streams*diff_streams, Ndz,Nkabs ,Nksca,Ng))
        LUT%S(iphi,itheta)%c = nil
        call h5write(diff_coeff_table_name,LUT%S(iphi,itheta)%c,iierr) ; ierr = ierr+int(iierr)
      endif

      if(.not. allocated(LUT%T(iphi,itheta)%c) ) then
        allocate(LUT%T(iphi,itheta)%c(dir_streams*dir_streams, Ndz,Nkabs ,Nksca,Ng))
        LUT%T(iphi,itheta)%c = nil
        call h5write(dir_coeff_table_name ,LUT%T(iphi,itheta)%c,iierr) ; ierr = ierr+int(iierr)
      endif
    endif

    total_size = Ng*Nksca*Nkabs *Ndz
    cnt=0
    do ig    =1,Ng
      do iksca    =1,Nksca    
        do ikabs   =1,Nkabs    
          do idz   =1,Ndz   
            cnt=cnt+1
            ! Check if we already calculated the coefficients and inform the other nodes
            ldone=.False.
            if(myid.eq.0) then
              if( all( LUT%S(iphi,itheta)%c( :, idz,ikabs ,iksca,ig).ge.zero).and.all( LUT%S(iphi,itheta)%c( :, idz,ikabs ,iksca,ig).le.one) .and. &
                  all( LUT%T(iphi,itheta)%c( :, idz,ikabs ,iksca,ig).ge.zero).and.all( LUT%T(iphi,itheta)%c( :, idz,ikabs ,iksca,ig).le.one) ) ldone = .True.
            endif
            call mpi_bcast(ldone,1 , MPI_LOGICAL, 0, comm, ierr)
            if(ldone) cycle

            if(myid.eq.0) print *,'direct dx',LUT%dx,'dz',LUT%pspace%dz(idz),'phi0,theta0',LUT%pspace%phi(iphi),LUT%pspace%theta(itheta),' :: ',LUT%pspace%kabs (ikabs ),LUT%pspace%ksca(iksca),LUT%pspace%g(ig),'(',100*cnt/total_size,'%)'
            do src=1,dir_streams
              call bmc_wrapper(src,LUT%dx,LUT%dy,LUT%pspace%dz(idz),LUT%pspace%kabs (ikabs ),LUT%pspace%ksca(iksca),LUT%pspace%g(ig),.True.,delta_scale,LUT%pspace%phi(iphi),LUT%pspace%theta(itheta),comm,S_diff,T_dir)
              if(myid.eq.0) LUT%T(iphi,itheta)%c( (src-1)*dir_streams+1:src*dir_streams, idz,ikabs ,iksca,ig) = T_dir
              if(myid.eq.0) LUT%S(iphi,itheta)%c( (src-1)*diff_streams+1:(src)*diff_streams, idz,ikabs ,iksca,ig) = S_diff
            enddo
          enddo !dz
        enddo !kabs
        
        if(myid.eq.0) then
          print *,'Writing direct table to file...'
          call h5write(diff_coeff_table_name,LUT%S(iphi,itheta)%c,iierr) ; ierr = ierr+int(iierr)
          call h5write(dir_coeff_table_name ,LUT%T(iphi,itheta)%c,iierr) ; ierr = ierr+int(iierr)
          print *,'done writing!',ierr
        endif

      enddo !ksca
    enddo !g
    if(myid.eq.0) print *,'done calculating direct coefficients'
end subroutine
!}}} 
!{{{ bmc_wrapper 
subroutine bmc_wrapper(src,dx,dy,dz,kabs ,ksca,g,dir,delta_scale,phi,theta,comm,S_diff,T_dir)
    integer(iintegers),intent(in) :: src
    integer(mpiint),intent(in) :: comm
    logical,intent(in) :: dir,delta_scale
    real(ireals),intent(in) :: dx,dy,dz,kabs ,ksca,g,phi,theta

    real(ireals),intent(out) :: S_diff(diff_streams),T_dir(dir_streams)

    real(ireals) :: bg(3)

    bg(1) = kabs
    bg(2) = ksca
    bg(3) = g

!    print *,'BMC :: calling bmc_get_coeff',bg,'src',src,'phi/theta',phi,theta,dz
    if(.not. bmc%initialized) call bmc%init(comm)
    call bmc%get_coeff(comm,bg,src,S_diff,T_dir,dir,delta_scale,phi,theta,dx,dy,dz)
    !        print *,'BMC :: dir',T_dir,'diff',S_diff
end subroutine
!}}}

!{{{ set parameter space
 real(ireals) function exp_param_to_index(val,range,N,expn)
    real(ireals),intent(in) :: val,range(2),expn
    integer(iintegers),intent(in) :: N
    real(ireals) :: expn1,k
    expn1=one/expn
    k=range(1)**expn1
    if(.not.valid_input(val,range)) continue
    exp_param_to_index = min(one*N, max(   one, one+ (N-1)*(val**expn1-k)/(range(2)**expn1-k)   ))
end function
 real(ireals) function exp_index_to_param(index,range,N,expn)
    real(ireals),intent(in) :: index,range(2),expn
    integer(iintegers),intent(in) :: N
    real(ireals) :: expn1
    expn1=one/expn
    exp_index_to_param = lin_index_to_param( index, range**expn1, N) ** expn
end function
 real(ireals) function lin_param_to_index(val,range,N)
    real(ireals),intent(in) :: val,range(2)
    integer(iintegers),intent(in) :: N
    if(.not.valid_input(val,range)) continue
    lin_param_to_index = min(one*N, max(   one, one+ (N-one) * (val-range(1) )/( range(2)-range(1) )   ))
end function
 real(ireals) function lin_index_to_param(index,range,N)
    real(ireals),intent(in) :: index,range(2)
    integer(iintegers),intent(in) :: N
    lin_index_to_param = range(1) + (index-one) * ( range(2)-range(1) ) / (N-1)
end function

subroutine set_parameter_space(ps,dx)
    type(parameter_space),intent(inout) :: ps
    real(ireals),intent(in) :: dx
    real(ireals) :: diameter ! diameter of max. cube size
    real(ireals),parameter :: maximum_transmission=one-1e-2_ireals
    integer(iintegers) :: k
    ! LUT Extend is already set in type definition in header, however we could overwrite the range here.

    ps%range_dz      = [ min(ps%range_dz(1), dx/10_ireals )  , min( ps%range_dz(2), dx ) ]
    diameter = sqrt(2*dx**2 +  ps%range_dz(2)**2 )

!    ps%range_dz      = [ dx/100_ireals , 2*dx ]
!    diameter = sqrt(2*dx**2 +  (2*ps%range_dz(2))**2 )

    ps%range_kabs(1) = - log(maximum_transmission) / diameter
    ps%range_ksca(1) = - log(maximum_transmission) / diameter

    !        ps%range_g       = [ zero    , one-1e-2 ]
    !        ps%range_phi     = [ zero    , one*90.  ]
    !        ps%range_theta   = [ zero    , one*90.  ]

    if( any([Ndz,Nkabs ,Nksca,Ng,Nphi,Ntheta].le.1) ) then
      print *,'Need at least 2 support points for LUT!'
      call exit()
    endif

    ps%dz_exponent=6.
    ps%kabs_exponent=10.
    ps%ksca_exponent=10.
    ps%g_exponent=.25

    do k=1,Ndz
      ps%dz(k)    = exp_index_to_param(one*k,ps%range_dz,Ndz, ps%dz_exponent )
    enddo
    do k=1,Nkabs 
      ps%kabs (k) = exp_index_to_param(one*k,ps%range_kabs ,Nkabs,ps%kabs_exponent )
    enddo
    do k=1,Nksca
      ps%ksca (k) = exp_index_to_param(one*k,ps%range_ksca ,Nksca,ps%ksca_exponent )
    enddo
    do k=1,Ng
      ps%g(k)     = exp_index_to_param(one*k,ps%range_g,Ng, ps%g_exponent )
    enddo
    do k=1,Nphi
      ps%phi(k)   = lin_index_to_param(one*k,ps%range_phi,Nphi)
    enddo
    do k=1,Ntheta
      ps%theta(k) = lin_index_to_param(one*k,ps%range_theta,Ntheta)
    enddo
end subroutine
!}}}

!{{{ diffuse coeff symmetry
function coeff_symmetry(isrc,coeff)
    real(ireals) :: coeff_symmetry(diff_streams)
    real(ireals),intent(in) :: coeff(ubound(diffLUT%S%c,1) )
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
    if(sum(coeff_symmetry).gt.1._ireals) then
      print *,'sum of diffuse coeff_symmetrys bigger one!',sum(coeff_symmetry),'for src=',isrc,'coeff_symmetry:',coeff_symmetry
      call exit()
    endif
end function
!}}}

subroutine catch_upper_limit_kabs(ps,kabs,ksca)
    ! If we hit the upper limit of the LUT for kabs, we try to scale ksca down,
    ! so that single scatter albedo w stays constant ..... this is a hack for
    ! really big absorption optical depths, where we dampen the scattering
    ! strength
    type(parameter_space),intent(in) :: ps
    real(ireals),intent(inout) :: kabs,ksca
    real(ireals) :: w,scaled_kabs,scaled_ksca
    if(kabs.gt.ps%range_kabs(2) ) then
      w = ksca/(kabs+ksca)
      scaled_kabs = ps%range_kabs(2)
      scaled_ksca = w*scaled_kabs / (one-w)
      if(optprop_debug) print *,'rescaling kabs because it is too big kabs',kabs,'->',scaled_kabs,'ksca',ksca,'->',scaled_ksca
      ksca=scaled_ksca
      kabs=scaled_kabs
    endif

    kabs = max( ps%range_kabs(1), kabs ) ! Also use lower limit of LUT table....
    ksca = max( ps%range_ksca(1), ksca ) ! Lets hope that we have a meaningful lower bound, as we will not get a warning for this.
end subroutine

!{{{ get_coeff routines
subroutine LUT_get_dir2dir(dz,in_kabs ,in_ksca,g,phi,theta,C)
    real(ireals),intent(in) :: dz,in_kabs ,in_ksca,g,phi,theta
    real(ireals),intent(out):: C(dir_streams**2)
    real(ireals) :: kabs,ksca
    integer(iintegers) :: i

    real(ireals) :: pti(6),weights(6),imap(2) ! index of point

    kabs = in_kabs; ksca = in_ksca
    call catch_upper_limit_kabs(dirLUT%pspace,kabs,ksca)

    pti = get_indices_6d(dz,kabs ,ksca,g,phi,theta,dirLUT%pspace)

    select case(interp_mode)
    case(1)
      ! Nearest neighbour
      C = dirLUT%T(nint(pti(5)), nint(pti(6)) )%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) )
    case(2)
      weights = modulo(pti,one)
      imap = [ exp_index_to_param(dble(floor  (pti(1))), dirLUT%pspace%range_dz , Ndz, dirLUT%pspace%dz_exponent), &
               exp_index_to_param(dble(ceiling(pti(1))), dirLUT%pspace%range_dz , Ndz, dirLUT%pspace%dz_exponent)  ]
      if(approx(imap(2), imap(1)) ) then
        weights(1)=zero
      else
        weights(1) = min(one, (dz -imap(1))/(imap(2)-imap(1)) )
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(2))), dirLUT%pspace%range_kabs , Nkabs, dirLUT%pspace%kabs_exponent ), &
               exp_index_to_param(dble(ceiling(pti(2))), dirLUT%pspace%range_kabs , Nkabs, dirLUT%pspace%kabs_exponent )  ]
      if(approx(imap(2),imap(1)) ) then
        weights(2)=zero
      else
        weights(2) = (kabs -imap(1))/(imap(2)-imap(1))
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(3))), dirLUT%pspace%range_ksca , Nksca, dirLUT%pspace%ksca_exponent ), &
               exp_index_to_param(dble(ceiling(pti(3))), dirLUT%pspace%range_ksca , Nksca, dirLUT%pspace%ksca_exponent )  ]
      if(approx(imap(2),imap(1)) ) then
        weights(3)=zero
      else
        weights(3) = (ksca -imap(1))/(imap(2)-imap(1))
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(4))), dirLUT%pspace%range_g , Ng, dirLUT%pspace%g_exponent ), &
               exp_index_to_param(dble(ceiling(pti(4))), dirLUT%pspace%range_g , Ng, dirLUT%pspace%g_exponent )  ]
      if(approx(imap(2),imap(1)) ) then
        weights(4)=zero
      else
        weights(4) = (g -imap(1))/(imap(2)-imap(1))
      endif

      call interp_4d(pti, weights, dirLUT%T(nint(pti(5)), nint(pti(6)) )%c, C)
!                              print *,'C',C
!                              call interp_4p2d(pti, weights, dirLUT%T( nint(pti(5)), nint(pti(6)) )%c, C)
!                              print *,'C',C
!                              call exit
!                              print *,"lin interp dir :: pti",pti
!                              print *,' weights',weights
!                              print *,' imap',imap,'C',C
!                              print *,'C',C
    case default
      print *,'interpolation mode',interp_mode,' not implemented yet! please choose something else! '
    end select
    if(optprop_debug) then
      !Check for energy conservation:
      ierr=0
      do i=1,dir_streams
        if(sum(C( (i-1)*dir_streams+1:i*dir_streams)).gt.one) ierr=ierr+1
      enddo
      if(ierr.ne.0) then
        print *,'Error in dir2dir coeffs :: ierr',ierr
        do i=1,dir_streams
          print *,'SUM dir2dir coeff for src ',i,' :: sum ',sum(C( (i-1)*dir_streams+1:i*dir_streams)),' :: coeff',C( (i-1)*dir_streams+1:i*dir_streams)
        enddo
      endif
    endif
end subroutine 
subroutine LUT_get_dir2diff(dz,in_kabs ,in_ksca,g,phi,theta,C)
    real(ireals),intent(in) :: dz,in_kabs ,in_ksca,g,phi,theta
    real(ireals),intent(out):: C(dir_streams*diff_streams)

    real(ireals) :: kabs,ksca
    real(ireals) :: pti(6),weights(6),imap(2) ! index of point
    integer(iintegers) :: i

    kabs = in_kabs; ksca = in_ksca
    call catch_upper_limit_kabs(dirLUT%pspace,kabs,ksca)

    pti = get_indices_6d(dz,kabs ,ksca,g,phi,theta,dirLUT%pspace)

    select case(interp_mode)
    case(1)
      ! Nearest neighbour
      C = dirLUT%S( nint(pti(5)), nint(pti(6)) )%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) )
    case(2)
      !                        print *,'linear interpolation not implemented yet!'
      weights = modulo(pti,one)
      imap = [ exp_index_to_param(dble(floor  (pti(1))), diffLUT%pspace%range_dz , Ndz, dirLUT%pspace%dz_exponent), &
          exp_index_to_param(dble(ceiling(pti(1))), diffLUT%pspace%range_dz , Ndz, dirLUT%pspace%dz_exponent)  ]
      if(approx(imap(2), imap(1)) ) then
        weights(1)=zero
      else
        weights(1) = min(one, (dz -imap(1))/(imap(2)-imap(1)) )
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(2))), dirLUT%pspace%range_kabs , Nkabs, dirLUT%pspace%kabs_exponent ), &
          exp_index_to_param(dble(ceiling(pti(2))), dirLUT%pspace%range_kabs , Nkabs, dirLUT%pspace%kabs_exponent )  ]
      if(approx(imap(2), imap(1)) ) then
        weights(2)=zero
      else
        weights(2) = min(one, (kabs -imap(1))/(imap(2)-imap(1)) )
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(3))), dirLUT%pspace%range_ksca , Nksca, dirLUT%pspace%ksca_exponent ), &
          exp_index_to_param(dble(ceiling(pti(3))), dirLUT%pspace%range_ksca , Nksca, dirLUT%pspace%ksca_exponent )  ]
      if(approx(imap(2), imap(1)) ) then
        weights(3)=zero
      else
        weights(3) = min(one, (ksca -imap(1))/(imap(2)-imap(1)) )
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(4))), dirLUT%pspace%range_g , Ng, dirLUT%pspace%g_exponent ), &
               exp_index_to_param(dble(ceiling(pti(4))), dirLUT%pspace%range_g , Ng, dirLUT%pspace%g_exponent )  ]
      if(approx(imap(2),imap(1)) ) then
        weights(4)=zero
      else
        weights(4) = (g -imap(1))/(imap(2)-imap(1))
      endif

      call interp_4d(pti, weights, dirLUT%S( nint(pti(5)), nint(pti(6)) )%c, C)
      !                        call interp_4p2d(pti, weights, dirLUT%S, C)
      !                        print *,'lin interp dir2diff weights',weights,'imap',imap,'C',C
      !                        call interp_6d_recursive(pti, weights, dirLUT%S, C)
      !                        print *,'lin interp dir2diff recursive weights',weights,'imap',imap,'C',C
      !                        call exit
    case default
      print *,'interpolation mode',interp_mode,' not implemented yet! please choose something else! '
    end select

    if(optprop_debug) then
      !Check for energy conservation:
      ierr=0
      do i=1,dir_streams
        if(sum(C( (i-1)*diff_streams+1:i*diff_streams)).gt.one) ierr=ierr+1
      enddo
      if(ierr.ne.0) then
        print *,'Error in dir2diff coeffs :: ierr',ierr
        do i=1,dir_streams
          print *,'SUM dir2dir coeff for src ',i,' :: sum ',sum(C( (i-1)*diff_streams+1:i*diff_streams)),' :: coeff',C( (i-1)*diff_streams+1:i*diff_streams)
        enddo
      endif
    endif
end subroutine
subroutine LUT_get_diff2diff(dz,in_kabs ,in_ksca,g,C)
    real(ireals),intent(in) :: dz,in_kabs ,in_ksca,g
    real(ireals),allocatable,intent(out):: C(:)

    real(ireals) :: kabs,ksca
    real(ireals) :: pti(4),weights(4),imap(2) ! index of point

    kabs = in_kabs; ksca = in_ksca
    call catch_upper_limit_kabs(diffLUT%pspace,kabs,ksca)

    allocate( C(ubound(diffLUT%S%c,1) ) )
    pti = get_indices_4d(dz,kabs ,ksca,g,diffLUT%pspace)

    select case(interp_mode)
    case(1)
      ! Nearest neighbour
      C = diffLUT%S%c(:,nint(pti(1)), nint(pti(2)), nint(pti(3)), nint(pti(4)) )
    case(2)
      !                        print *,'linear interpolation not implemented yet!'
      weights = modulo(pti,one)
      imap = [ exp_index_to_param(dble(floor  (pti(1))), diffLUT%pspace%range_dz , Ndz, diffLUT%pspace%dz_exponent), &
          exp_index_to_param(dble(ceiling(pti(1))), diffLUT%pspace%range_dz , Ndz, diffLUT%pspace%dz_exponent)  ]
      if(approx(imap(2),imap(1)) ) then
        weights(1)=zero
      else
        weights(1) = min(one, (dz -imap(1))/(imap(2)-imap(1)) )
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(2))), diffLUT%pspace%range_kabs , Nkabs, diffLUT%pspace%kabs_exponent), &
          exp_index_to_param(dble(ceiling(pti(2))), diffLUT%pspace%range_kabs , Nkabs, diffLUT%pspace%kabs_exponent)  ]
      if(approx(imap(2), imap(1)) ) then
        weights(2)=zero
      else
        weights(2) = min(one, (kabs -imap(1))/(imap(2)-imap(1)) )
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(3))), diffLUT%pspace%range_ksca , Nksca, diffLUT%pspace%ksca_exponent), &
          exp_index_to_param(dble(ceiling(pti(3))), diffLUT%pspace%range_ksca , Nksca, diffLUT%pspace%ksca_exponent)  ]
      if(approx(imap(2), imap(1)) ) then
        weights(3)=zero
      else
        weights(3) = min(one, (ksca -imap(1))/(imap(2)-imap(1)) )
      endif

      imap = [ exp_index_to_param(dble(floor  (pti(4))), dirLUT%pspace%range_g , Ng, dirLUT%pspace%g_exponent ), &
               exp_index_to_param(dble(ceiling(pti(4))), dirLUT%pspace%range_g , Ng, dirLUT%pspace%g_exponent )  ]
      if(approx(imap(2),imap(1)) ) then
        weights(4)=zero
      else
        weights(4) = (g -imap(1))/(imap(2)-imap(1))
      endif

      call interp_4d(pti, weights, diffLUT%S%c, C)
      !                        print *,'lin interp diff weights',weights,'imap',imap,'C',C
      !                        call interp_4d_recursive(pti, weights, diffLUT%S, C)
      !                        print *,'recursive lin interp diff weights',weights,'imap',imap,'C',C
      !                        call exit
    case default
      print *,'interpolation mode',interp_mode,' not implemented yet! please choose something else! '
    end select
end subroutine 
!}}}

!{{{ index mappings
function get_indices_4d(dz,kabs ,ksca,g,ps)
    real(ireals) :: get_indices_4d(4)
    real(ireals),intent(in) :: dz,kabs ,ksca,g
    type(parameter_space),intent(in) :: ps

    get_indices_4d(1) = exp_param_to_index( dz    ,ps%range_dz   ,Ndz  , ps%dz_exponent )
    get_indices_4d(2) = exp_param_to_index( kabs  ,ps%range_kabs ,Nkabs, ps%kabs_exponent )
    get_indices_4d(3) = exp_param_to_index( ksca  ,ps%range_ksca ,Nksca, ps%ksca_exponent )
    get_indices_4d(4) = exp_param_to_index( g     ,ps%range_g    ,Ng   , ps%g_exponent )
end function
function get_indices_6d(dz,kabs ,ksca,g,phi,theta,ps)
    real(ireals) :: get_indices_6d(6)
    real(ireals),intent(in) :: dz,kabs ,ksca,g,phi,theta
    type(parameter_space),intent(in) :: ps

    get_indices_6d(1:4) = get_indices_4d(dz,kabs ,ksca,g,ps)
    get_indices_6d(5) = lin_param_to_index( phi   ,ps%range_phi  ,Nphi )
    get_indices_6d(6) = lin_param_to_index( theta ,ps%range_theta,Ntheta)
end function
!}}}
!{{{ check input values
logical function valid_input(val,range)
    real(ireals),intent(in) :: val,range(2)
    if(val.lt.range(1) .or. val.gt.range(2) ) then 
      valid_input=.False.
      print *,'ohoh, this val is not in the optprop database range!',val,'not in',range
    else
      valid_input=.True.
    endif
end function
!}}}



end module
