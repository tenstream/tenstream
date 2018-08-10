module test_tenstr_atm
  use iso_fortran_env, only: REAL32, REAL64
  use m_data_parameters, only : &
    init_mpi_data_parameters,   &
    iintegers, ireals, mpiint,  &
    zero, one, default_str_len

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm

  use pfunit_mod

  implicit none

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    ! Tidy up
    class (MpiTestMethod), intent(inout) :: this
    integer(mpiint) :: ierr
    call PetscFinalize(ierr)
  end subroutine teardown

  @test(npes =[1,2])
  subroutine concat_tenstr_atm(this)
    class (MpiTestMethod), intent(inout) :: this

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid

    integer(iintegers),parameter :: ncol=3, nzp=18
    real(ireals), dimension(nzp+1, ncol) :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1, ncol) :: tlev ! Temperature on layer interfaces [K]

    real(ireals), dimension(nzp, ncol) :: tlay !, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp, ncol) :: lwc, reliq

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(default_str_len),parameter :: atm_filename='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: k, kt, icld

    type(t_tenstr_atm) :: atm
    logical, parameter :: lverbose=.True.

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 500 hPa and a
    ! Temperature difference of 32.5K
    do k=1,nzp+1
      plev(k,:) = 1000_ireals - (k-one)*500._ireals/(nzp)
      tlev(k,:) = 288._ireals - (k-one)*(5*6.5_ireals)/(nzp)
    enddo

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 0

    icld = (nzp+1)/2
    lwc  (icld, :) = 1e-2_ireals
    reliq(icld, :) = 10._ireals

    tlev (icld  , :) = 288._ireals
    tlev (icld+1, :) = tlev (icld ,:)

    ! mean layer temperature is approx. arithmetic mean between layer interfaces
    tlay = (tlev(1:nzp,:) + tlev(2:nzp+1,:))/2

    call setup_tenstr_atm(comm, .False., atm_filename, plev, tlev, atm, &
      d_tlay=tlay, d_lwc=lwc, d_reliq=reliq)

    if(lverbose .and. myid.eq.0) then
      print *,'Shape of background profile:', shape(atm%bg_atm%plev)
      print *,'Shape of dynamics grid:', shape(plev)
      print *,'Shape of merged grid:', shape(atm%plev)
      print *,'atm_ke', atm%atm_ke
    endif

    if(lverbose .and. myid.eq.0) then
      do k = 1,size(atm%bg_atm%plev,1)
        print *,'Pressure on Background', k,':', atm%bg_atm%plev(k), 'Tlev', atm%bg_atm%tlev(k)
      enddo

      do k = size(plev,1),1,-1
        print *,'Pressure on Dynamics Grid', k,':', plev(k,1), 'Tlev', tlev(k,1)
      enddo

      do k = size(atm%plev,1), 1, -1
        print *,'Pressure on Merged Grid', k,':', atm%plev(k,1),'Tlev', atm%tlev(k,1)
      enddo
    endif

    @mpiassertEqual([113], shape(atm%bg_atm%plev), 'background pressure profile has wrong shape')
    @mpiassertEqual([54+nzp,3], shape(atm%plev), 'merged pressure grid has wrong shape')


    do k = 1, size(plev,1)
      @mpiassertEqual(plev(k,:), atm%plev(k,:), 1e-8_ireals)
      @mpiassertEqual(tlev(k,:), atm%tlev(k,:), 1e-8_ireals)
    enddo

    do k = 1, atm%atm_ke
      kt = size(atm%plev,1) - k + 1 ! index from top
      @mpiassertEqual(atm%bg_atm%plev(k), atm%plev(kt, 1))
      @mpiassertEqual(atm%bg_atm%plev(k), atm%plev(kt, ncol))
      @mpiassertEqual(atm%bg_atm%tlev(k), atm%tlev(kt, 1))
      @mpiassertEqual(atm%bg_atm%tlev(k), atm%tlev(kt, ncol))
    enddo

    ! Layer Quantities
    if(lverbose .and. myid.eq.0) then
      do k = 1, size(tlay,1)
        print *,'Tlay on Dynamics Grid', k,':', tlay(k,1), 'cld', lwc(k,1), reliq(k,1)
      enddo

      do k = size(atm%tlay,1), 1, -1
        print *,'Tlay on Merged Grid', k,':', atm%tlay(k,1), 'cld', atm%lwc(k,1), atm%reliq(k,1),':dz', atm%dz(k,1)
      enddo
    endif

    do k = 1, size(tlay,1)
      @mpiassertEqual(tlay(k,:), atm%tlay(k,:), 1e-8_ireals)
    enddo
    do k = 1, atm%atm_ke
      kt = size(atm%tlay,1) - k + 1 ! index from top
      @mpiassertEqual(atm%bg_atm%tlay(k), atm%tlay(kt, ncol))
    enddo

    @mpiassertEqual(lwc  (icld, :), atm%lwc  (icld,:), 1e-8_ireals)
    do k = 1, icld-1
      @mpiassertEqual(zero, atm%lwc  (k,:))
    enddo
    do k = icld+1, size(atm%lwc,1)
      @mpiassertEqual(zero, atm%lwc  (k,:))
    enddo

    @mpiassertEqual(reliq(icld, :), atm%reliq(icld,:), 1e-8_ireals)
    do k = 1, icld-1
      @mpiassertEqual(zero, atm%reliq  (k,:))
    enddo
    do k = icld+1, size(atm%reliq,1)
      @mpiassertEqual(zero, atm%reliq  (k,:))
    enddo
  end subroutine

  @test(npes =[1,2])
  subroutine test_dont_update_bg_entries(this)
    class (MpiTestMethod), intent(inout) :: this

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid

    integer(iintegers),parameter :: ncol=3, nzp=18
    real(ireals), dimension(nzp+1, ncol) :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1, ncol) :: tlev ! Temperature on layer interfaces [K]

    real(ireals), dimension(nzp, ncol) :: tlay !, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp, ncol) :: lwc, reliq

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(default_str_len),parameter :: atm_filename='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: k, kt, icld

    type(t_tenstr_atm) :: atm, atm2
    logical, parameter :: lverbose=.True.

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 500 hPa and a
    ! Temperature difference of 32.5K
    do k=1,nzp+1
      plev(k,:) = 1000_ireals - (k-one)*500._ireals/(nzp)
      tlev(k,:) = 288._ireals - (k-one)*(5*6.5_ireals)/(nzp)
    enddo

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 0

    icld = (nzp+1)/2
    lwc  (icld, :) = 1e-2_ireals
    reliq(icld, :) = 10._ireals

    tlev (icld  , :) = 288._ireals
    tlev (icld+1, :) = tlev (icld ,:)

    ! mean layer temperature is approx. arithmetic mean between layer interfaces
    tlay = (tlev(1:nzp,:) + tlev(2:nzp+1,:))/2

    call setup_tenstr_atm(comm, .False., atm_filename, plev, tlev, atm, &
      d_tlay=tlay, d_lwc=lwc, d_reliq=reliq)

    if(lverbose .and. myid.eq.0) then
      print *,'Shape of background profile:', shape(atm%bg_atm%plev)
      print *,'Shape of dynamics grid:', shape(plev)
      print *,'Shape of merged grid:', shape(atm%plev)
      print *,'atm_ke', atm%atm_ke
    endif

    if(lverbose .and. myid.eq.0) then
      do k = 1,size(atm%bg_atm%plev,1)
        print *,'Pressure on Background', k,':', atm%bg_atm%plev(k), 'Tlev', atm%bg_atm%tlev(k)
      enddo

      do k = size(plev,1),1,-1
        print *,'Pressure on Dynamics Grid', k,':', plev(k,1), 'Tlev', tlev(k,1)
      enddo

      do k = size(atm%plev,1), 1, -1
        print *,'Pressure on Merged Grid', k,':', atm%plev(k,1),'Tlev', atm%tlev(k,1)
      enddo
    endif

    ! Layer Quantities
    if(lverbose .and. myid.eq.0) then
      do k = 1, size(tlay,1)
        print *,'Tlay on Dynamics Grid', k,':', tlay(k,1), 'cld', lwc(k,1), reliq(k,1)
      enddo

      do k = size(atm%tlay,1), 1, -1
        print *,'Tlay on Merged Grid', k,':', atm%tlay(k,1), 'cld', atm%lwc(k,1), atm%reliq(k,1),':dz', atm%dz(k,1)
      enddo
    endif

    ! Setup a second atmoshpere
    call setup_tenstr_atm(comm, .False., atm_filename, plev, tlev, atm2, &
      d_tlay=tlay, d_lwc=lwc, d_reliq=reliq)

    ! Now override the background profile
    atm2%bg_atm%plev = -1
    atm2%bg_atm%tlev = -1
    atm2%bg_atm%tlay = -1

    ! Update atm2
    call setup_tenstr_atm(comm, .False., atm_filename, plev+one, tlev+one, atm2, &
      d_tlay=tlay+one, d_lwc=lwc*-one, d_reliq=reliq*-one)

    ! and make sure that the background profile values have not changed
    @mpiassertEqual(atm%plev(atm%d_ke1+1:ubound(atm%plev,1),:), atm2%plev(atm2%d_ke1+1:ubound(atm%plev,1),:))
    @mpiassertEqual(atm%tlev(atm%d_ke1+1:ubound(atm%tlev,1),:), atm2%tlev(atm2%d_ke1+1:ubound(atm%tlev,1),:))

    @mpiassertEqual(atm%tlay(atm%d_ke+1:ubound(atm%tlay,1),:), atm2%tlay(atm2%d_ke+1:ubound(atm%tlay,1),:))

    @mpiassertEqual(atm%lwc  (atm%d_ke+1:ubound(atm%lwc  ,1),:), atm2%lwc  (atm2%d_ke+1:ubound(atm%lwc  ,1),:))
    @mpiassertEqual(atm%reliq(atm%d_ke+1:ubound(atm%reliq,1),:), atm2%reliq(atm2%d_ke+1:ubound(atm%reliq,1),:))

    ! whereas the dynamics grid should be changed

    @mpiassertEqual(atm%plev(1:atm%d_ke1, :) + one, atm2%plev(1:atm%d_ke1, :))
    @mpiassertEqual(atm%tlev(1:atm%d_ke1, :) + one, atm2%tlev(1:atm%d_ke1, :))

    @mpiassertEqual(atm%tlay (1:atm%d_ke, :) + one, atm2%tlay (1:atm%d_ke, :))
    @mpiassertEqual(atm%lwc  (1:atm%d_ke, :) *-one, atm2%lwc  (1:atm%d_ke, :))
    @mpiassertEqual(atm%reliq(1:atm%d_ke, :) *-one, atm2%reliq(1:atm%d_ke, :))

    if(lverbose .and. myid.eq.0) then
      do k = 1,size(atm2%bg_atm%plev,1)
        print *,'Pressure on Background', k,':', atm2%bg_atm%plev(k), 'Tlev', atm2%bg_atm%tlev(k)
      enddo

      do k = size(plev,1),1,-1
        print *,'Pressure on Dynamics Grid', k,':', plev(k,1), 'Tlev', tlev(k,1)
      enddo

      do k = size(atm2%plev,1), 1, -1
        print *,'Pressure on Merged Grid', k,':', atm2%plev(k,1),'Tlev', atm2%tlev(k,1)
      enddo
    endif

    ! Layer Quantities
    if(lverbose .and. myid.eq.0) then
      do k = 1, size(tlay,1)
        print *,'Tlay on Dynamics Grid', k,':', tlay(k,1), 'cld', lwc(k,1), reliq(k,1)
      enddo

      do k = size(atm2%tlay,1), 1, -1
        print *,'Tlay on Merged Grid', k,':', atm2%tlay(k,1), 'cld', atm2%lwc(k,1), atm2%reliq(k,1),':dz', atm2%dz(k,1)
      enddo
    endif
  end subroutine
end module
