module test_pprts_rrtm_icollapse
  use iso_fortran_env, only: REAL32, REAL64
  use m_data_parameters, only : &
    init_mpi_data_parameters,   &
    iintegers, ireals, mpiint,  &
    i1, zero, one, default_str_len
  use m_helper_functions, only : linspace, itoa

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

  use m_pprts_base, only : t_solver_3_10

  use pfunit_mod

  implicit none

  type(t_solver_3_10) :: solver
  type(t_solver_3_10) :: solver_collapsed
  type(t_tenstr_atm) :: atm

  logical,parameter :: ldebug=.True.

  integer(iintegers),parameter :: nxp=3, nyp=3, nzp=20 ! local domain size for each rank
  real(ireals),parameter :: dx=100, dy=dx              ! horizontal grid spacing in [m]
  real(ireals),parameter :: phi0=180, theta0=60        ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
  real(ireals),parameter :: albedo_th=0, albedo_sol=.3 ! broadband ground albedo for solar and thermal spectrum
  real(ireals),parameter :: atolerance = 1             ! absolute tolerance when regression testing fluxes
  integer(iintegers),allocatable :: nxproc(:), nyproc(:)

contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y

    real(ireals), dimension(nzp+1,nxp,nyp), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1,nxp,nyp), target :: tlev ! Temperature on layer interfaces [K]

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `pprts_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp,nxp,nyp), target :: lwc, reliq

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(default_str_len),parameter :: atm_filename='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: k, icld

    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq

    call init_mpi_data_parameters(this%getMpiCommunicator())

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if(N_ranks_y*N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    endif
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate(nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 750 hPa and a
    ! Temperature difference of 20K
    do k=1,nzp+1
      plev(k,:,:) = linspace(k, [1e3_ireals, 750._ireals], nzp+1)
      tlev(k,:,:) = linspace(k, [290._ireals, 270._ireals], nzp+1)
    enddo

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 0

    icld = int(real(nzp+1)/2)
    lwc  (icld, :,:) = 1e-2_ireals
    reliq(icld, :,:) = 10._ireals

    tlev (icld  , :,:) = 288._ireals
    tlev (icld+1, :,:) = tlev (icld  , :,:)

    ! Setup Atmosphere
    pplev(1:size(plev,1),1:size(plev,2)*size(plev,3)) => plev
    ptlev(1:size(plev,1),1:size(tlev,2)*size(tlev,3)) => tlev
    plwc (1:size(lwc ,1),1:size(lwc ,2)*size(lwc ,3)) => lwc
    preliq(1:size(reliq,1),1:size(reliq,2)*size(reliq,3)) => reliq

    call setup_tenstr_atm(comm, .False., atm_filename, &
      pplev, ptlev, atm, &
      d_lwc=plwc, d_reliq=preliq)
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    ! Tidy up
    call destroy_pprts_rrtmg(solver, lfinalizepetsc=.False.)
    call destroy_pprts_rrtmg(solver_collapsed, lfinalizepetsc=.True.)
    call destroy_tenstr_atm(atm)
    if(allocated(nxproc)) deallocate(nxproc)
    if(allocated(nyproc)) deallocate(nyproc)
  end subroutine teardown

  @test(npes =[2,1])
  subroutine pprts_rrtm_icollapse_solar(this)
    class (MpiTestMethod), intent(inout) :: this

    logical :: lthermal, lsolar
    integer(mpiint) :: comm, myid

    integer(iintegers) :: k, k0

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals),allocatable, dimension(:,:,:) :: edir0, edn0, eup0, abso0 ! [nlev_merged(-1), nxp, nyp]
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    comm     = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    ! For comparison, compute lw and sw separately
    if(myid.eq.0 .and. ldebug) print *,'Computing Solar Radiation:'
    lthermal=.False.; lsolar=.True.

    solver%solvername = '3_10_base_solver'
    call pprts_rrtmg(comm, solver, atm, nxp, nyp, &
      dx, dy, phi0, theta0,   &
      albedo_th, albedo_sol,  &
      lthermal, lsolar,       &
      edir0, edn0, eup0, abso0, &
      nxproc=nxproc, nyproc=nyproc, &
      opt_time=zero)

    solver_collapsed%solvername = '3_10_collapsed_solver'
    call pprts_rrtmg(comm, solver_collapsed, atm, nxp, nyp, &
      dx, dy, phi0, theta0,   &
      albedo_th, albedo_sol,  &
      lthermal, lsolar,       &
      edir, edn, eup, abso, &
      icollapse=-i1, &
      nxproc=nxproc, nyproc=nyproc, &
      opt_time=zero)

    ! Determine number of actual output levels from returned flux arrays.
    ! We dont know the size before hand because we get the fluxes on the merged
    ! grid (dynamics + background profile)
    if(myid.eq.0) then
      if(ldebug) then
        print *,'---------------------------------- Base Solve'
        do k=1,ubound(edn0,1)-1
          print *,k,'edir', edir0(k,1,1), 'edn', edn0(k,1,1), 'eup', eup0(k,1,1), abso0(k,1,1)
        enddo
        print *,k,'edir', edir0(k,1,1), 'edn', edn0(k,1,1), 'eup', eup0(k,1,1)
      endif
    endif
    if(myid.eq.0) then
      if(ldebug) then
        print *,''
        print *,'---------------------------------- Collapsed'
        do k=1,ubound(edn,1)-1
          print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(k,1,1)
        enddo
        print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1)
        print *,''
      endif
    endif

    @assertEqual(edir0(1,1,1), edir(1,1,1), 1e-3_ireals, 'TOA solar incoming radiation')

    do k=2,ubound(edir,1)
      k0 = ubound(edir0,1) - ubound(edir,1) + k
      @assertEqual(edir0(k0,1,1), edir(k,1,1), 1e-0_ireals, 'solar direct  radiation @level '//itoa(k0)//' collapsed '//itoa(k))
      @assertEqual(edn0 (k0,1,1), edn (k,1,1), 1e-0_ireals, 'solar diff dn radiation @level '//itoa(k0)//' collapsed '//itoa(k))
      @assertEqual(eup0 (k0,1,1), eup (k,1,1), 1e-0_ireals, 'solar diff up radiation @level '//itoa(k0)//' collapsed '//itoa(k))
    enddo
    do k=2,ubound(abso,1)
      k0 = ubound(abso0,1) - ubound(abso,1) + k
      @assertEqual(abso0(k0,1,1), abso(k,1,1), 1e-4_ireals, 'solar abso @level '//itoa(k0)//' collapsed '//itoa(k))
    enddo
  end subroutine


  @test(npes =[2,1])
  subroutine pprts_rrtm_icollapse_thermal(this)
    class (MpiTestMethod), intent(inout) :: this

    logical :: lthermal, lsolar
    integer(mpiint) :: comm, myid

    integer(iintegers) :: k, k0

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals),allocatable, dimension(:,:,:) :: edir0, edn0, eup0, abso0 ! [nlev_merged(-1), nxp, nyp]
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    comm     = this%getMpiCommunicator()
    myid     = this%getProcessRank()

    ! For comparison, compute lw and sw separately
    if(myid.eq.0 .and. ldebug) print *,'Computing Thermal Radiation:'
    lthermal=.True.; lsolar=.False.

    solver%solvername = '3_10_base_solver'
    call pprts_rrtmg(comm, solver, atm, nxp, nyp, &
      dx, dy, phi0, theta0,   &
      albedo_th, albedo_sol,  &
      lthermal, lsolar,       &
      edir0, edn0, eup0, abso0, &
      nxproc=nxproc, nyproc=nyproc, &
      opt_time=zero)

    solver_collapsed%solvername = '3_10_collapsed_solver'
    call pprts_rrtmg(comm, solver_collapsed, atm, nxp, nyp, &
      dx, dy, phi0, theta0,   &
      albedo_th, albedo_sol,  &
      lthermal, lsolar,       &
      edir, edn, eup, abso, &
      icollapse=-i1, &
      nxproc=nxproc, nyproc=nyproc, &
      opt_time=zero)

    ! Determine number of actual output levels from returned flux arrays.
    ! We dont know the size before hand because we get the fluxes on the merged
    ! grid (dynamics + background profile)
    if(myid.eq.0) then
      if(ldebug) then
        print *,'---------------------------------- Base Solve'
        do k=1,ubound(edn0,1)-1
          print *,k, 'edn', edn0(k,1,1), 'eup', eup0(k,1,1), abso0(k,1,1)
        enddo
        print *,k, 'edn', edn0(k,1,1), 'eup', eup0(k,1,1)
      endif
    endif
    if(myid.eq.0) then
      if(ldebug) then
        print *,''
        print *,'---------------------------------- Collapsed'
        do k=1,ubound(edn,1)-1
          print *,k, 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(k,1,1)
        enddo
        print *,k, 'edn', edn(k,1,1), 'eup', eup(k,1,1)
        print *,''
      endif
    endif

    @assertEqual(eup0(1,1,1), eup(1,1,1), 2e0_ireals, 'TOA thermal up radiation')
    @assertEqual(edn0(1,1,1), edn(1,1,1), 1e-8_ireals, 'TOA thermal up radiation')

    do k=2,ubound(edn,1)
      k0 = ubound(edn0,1) - ubound(edn,1) + k
      @assertEqual(edn0(k0,1,1), edn(k,1,1), 3e-0_ireals, 'thermal down radiation @level '//itoa(k0)//' collapsed '//itoa(k))
      @assertEqual(eup0(k0,1,1), eup(k,1,1), 3e-0_ireals, 'thermal up radiation @level '//itoa(k0)//' collapsed '//itoa(k))
    enddo
    do k=2,ubound(abso,1)
      k0 = ubound(abso0,1) - ubound(abso,1) + k
      @assertEqual(abso0(k0,1,1), abso(k,1,1), 1e-2_ireals, 'thermal abso @level '//itoa(k0)//' collapsed '//itoa(k))
    enddo
  end subroutine
end module
