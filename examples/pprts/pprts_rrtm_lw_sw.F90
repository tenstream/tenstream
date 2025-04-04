module m_example_pprts_rrtm_lw_sw

#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, default_str_len

  use m_helper_functions, only: linspace, CHKERR, spherical_2_cartesian, meanval, get_petsc_opt

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only: pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

  implicit none

contains
  subroutine ex_pprts_rrtm_lw_sw(comm, nxp, nyp, nzp, dx, dy, &
      & phi0, theta0, albedo_th, albedo_sol, &
      & lthermal, lsolar, atm_filename, &
      & edir, edn, eup, abso, &
      & vlwc, viwc)
    integer(mpiint), intent(in) :: comm
    integer(iintegers), intent(in) :: nxp, nyp, nzp   ! local domain size for each rank
    real(ireals), intent(in) :: dx, dy                ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0          ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: albedo_th, albedo_sol ! broadband ground albedo for solar and thermal spectrum
    logical, intent(in) :: lthermal, lsolar                       ! switches to enable/disable spectral integration
    character(len=*), intent(in) :: atm_filename ! ='afglus_100m.dat'
    real(ireals), intent(in), optional :: vlwc, viwc            ! liquid/ice water content to be set in a layer

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, myid, N_ranks_x, N_ranks_y, ierr

    real(ireals), dimension(nzp + 1, nxp, nyp), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp + 1, nxp, nyp), target :: tlev ! Temperature on layer interfaces [K]

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp, nxp, nyp), target :: lwc, reliq, iwc, reice

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    logical :: lflg

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev, icld, iter
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev, plwc, preliq, piwc, preice

    real(ireals) :: sundir(3)

    logical, parameter :: ldebug = .true.

    class(t_solver), allocatable :: pprts_solver
    type(t_tenstr_atm) :: atm

    call MPI_COMM_SIZE(comm, numnodes, ierr)
    call MPI_COMM_RANK(comm, myid, ierr)

    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if (N_ranks_y * N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    end if
    if (ldebug .and. myid .eq. 0) print *, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate (nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate (nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 500 hPa and a
    ! Temperature difference of 50K
    do k = 1, nzp + 1
      plev(k, :, :) = linspace(k, [1e3_ireals, 500._ireals], nzp + 1)
      tlev(k, :, :) = linspace(k, [288._ireals, 250._ireals], nzp + 1)
    end do

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)
    ! h2ovmr = zero
    ! o3vmr  = zero
    ! co2vmr = zero
    ! ch4vmr = zero
    ! n2ovmr = zero
    ! o2vmr  = zero

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 0

    if (present(vlwc)) then
      icld = int(real(nzp + 1) / 2)
      lwc(icld, :, :) = vlwc
      reliq(icld, :, :) = 10
      tlev(icld + 1, :, :) = tlev(icld, :, :)
    end if

    iwc = 0
    reice = 0

    if (present(viwc)) then
      icld = nzp
      iwc(icld, :, :) = viwc
      reice(icld, :, :) = 60
      tlev(icld + 1, :, :) = tlev(icld, :, :)
    end if

    if (myid .eq. 0 .and. ldebug) print *, 'Setup Atmosphere...'

    pplev(1:size(plev, 1), 1:size(plev, 2) * size(plev, 3)) => plev
    ptlev(1:size(tlev, 1), 1:size(tlev, 2) * size(tlev, 3)) => tlev
    plwc(1:size(lwc, 1), 1:size(lwc, 2) * size(lwc, 3)) => lwc
    preliq(1:size(reliq, 1), 1:size(reliq, 2) * size(reliq, 3)) => reliq
    piwc(1:size(iwc, 1), 1:size(iwc, 2) * size(iwc, 3)) => iwc
    preice(1:size(reice, 1), 1:size(reice, 2) * size(reice, 3)) => reice

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          pplev, ptlev, atm, &
                          d_lwc=plwc, d_reliq=preliq, &
                          d_iwc=piwc, d_reice=preice)

    sundir = spherical_2_cartesian(phi0, theta0)

    call allocate_pprts_solver_from_commandline(pprts_solver, '3_10', ierr); call CHKERR(ierr)

    iter = 1
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-iter", iter, lflg, ierr)

    do k = 1, iter
      call pprts_rrtmg(comm, pprts_solver, atm, nxp, nyp, &
                       dx, dy, sundir, &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, &
                       opt_time=zero)
    end do

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          if (allocated(edir)) then
            print *, k, 'edir', meanval(edir(k, :, :)), 'edn', meanval(edn(k, :, :)), 'eup', meanval(eup(k, :, :)), &
              & 'abso', meanval(abso(min(nlev - 1, k), :, :))
          else
            print *, k, 'edn', meanval(edn(k, :, :)), 'eup', meanval(eup(k, :, :)), meanval(abso(min(nlev - 1, k), :, :))
          end if
        end do
      end if

      if (allocated(edir)) &
        print *, 'surface :: direct flux', meanval(edir(nlev, :, :))
      print *, 'surface :: downw flux ', meanval(edn(nlev, :, :))
      print *, 'surface :: upward fl  ', meanval(eup(nlev, :, :))
      print *, 'surface :: absorption ', meanval(abso(nlev - 1, :, :))

      if (allocated(edir)) &
        print *, 'TOA :: direct flux', meanval(edir(1, :, :))
      print *, 'TOA :: downw flux ', meanval(edn(1, :, :))
      print *, 'TOA :: upward fl  ', meanval(eup(1, :, :))
      print *, 'TOA :: absorption ', meanval(abso(1, :, :))
    end if

    ! Tidy up
    call destroy_pprts_rrtmg(pprts_solver, lfinalizepetsc=.true.)
    call destroy_tenstr_atm(atm)
  end subroutine

end module
