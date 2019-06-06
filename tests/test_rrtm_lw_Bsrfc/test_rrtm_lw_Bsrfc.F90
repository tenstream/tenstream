module test_rrtm_lw_Bsrfc
  use iso_fortran_env, only: REAL32, REAL64
  use m_data_parameters, only : &
    init_mpi_data_parameters,   &
    iintegers, ireals, mpiint,  &
    zero, one, default_str_len

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

  use m_pprts_base, only : t_solver_3_10

  use pfunit_mod

  implicit none

  type(t_solver_3_10) :: solver
  type(t_tenstr_atm) :: atm
contains

  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    call init_mpi_data_parameters(this%getMpiCommunicator())
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    ! Tidy up
    call destroy_pprts_rrtmg(solver, lfinalizepetsc=.True.)
    call destroy_tenstr_atm(atm)
  end subroutine teardown

  @test(npes =[1,2])
  subroutine rrtm_lw_Bsrfc(this)
    class (MpiTestMethod), intent(inout) :: this

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y

    integer(iintegers),parameter :: nxp=3, nyp=3, nzp=60 ! local domain size for each rank
    integer(iintegers),parameter :: icollapse=-1
    real(ireals),parameter :: dx=500, dy=dx              ! horizontal grid spacing in [m]
    real(ireals),parameter :: phi0=180, theta0=60        ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals),parameter :: albedo_th=0.1, albedo_sol=.3 ! broadband ground albedo for solar and thermal spectrum
    real(ireals),parameter :: atolerance = 1             ! absolute tolerance when regression testing fluxes

    real(ireals), dimension(nzp+1,nxp,nyp), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1,nxp,nyp), target :: tlev ! Temperature on layer interfaces [K]

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `pprts_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: tlay, h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp,nxp,nyp), target :: lwc, reliq
    real(ireals), dimension(nxp,nyp), target :: tskin

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(default_str_len),parameter :: atm_filename='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev, icld, iter
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq
    real(ireals), pointer, dimension(:)   :: ptskin
    real(ireals) :: target_Eup

    logical,parameter :: ldebug=.True.
    logical :: lthermal, lsolar

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

    ! Start with a dynamics grid ranging from 1000 hPa up to 500 hPa and a
    ! Temperature difference of 50K
    do k=1,nzp+1
      plev(k,:,:) = 1000_ireals - real(k-1, ireals)*500._ireals/real(nzp, ireals)
      tlev(k,:,:) = 288._ireals - real(k-1, ireals)* 50._ireals/real(nzp, ireals)
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
    ptskin(1:size(tskin)) => tskin

    do iter = 0, 30, 10
      tskin = 288 + iter
      call setup_tenstr_atm(comm, .False., atm_filename, &
        pplev, ptlev, atm, &
        d_lwc=plwc, d_reliq=preliq, &
        d_skin_temperature=ptskin )

      lsolar=.False.
      lthermal=.True.

      call pprts_rrtmg(comm, solver, atm, nxp, nyp, &
        dx, dy, phi0, theta0,   &
        albedo_th, albedo_sol,  &
        lthermal, lsolar,       &
        edir, edn, eup, abso,   &
        nxproc=nxproc, nyproc=nyproc, &
        opt_time=zero, icollapse=icollapse)

      ! Determine number of actual output levels from returned flux arrays.
      ! We dont know the size before hand because we get the fluxes on the merged
      ! grid (dynamics + background profile)
      nlev = ubound(edn,1)
      if(myid.eq.0) then
        if(ldebug) then
          do k=1,nlev
            print *,k, 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
          enddo
        endif
      endif
      target_Eup = 5.67e-8_ireals*tskin(1,1)**4 * (one - albedo_th) + edn(nlev,1,1) * albedo_th
      @assertEqual(target_Eup, eup(nlev,1,1), one, 'emission not as expected. it should use the skin temperature emission?!')
    enddo
  end subroutine
end module
