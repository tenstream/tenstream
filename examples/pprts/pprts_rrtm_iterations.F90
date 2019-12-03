module m_example_pprts_rrtm_iterations

#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, default_str_len

  use m_helper_functions, only : linspace, CHKERR, meanval
  use m_search, only: search_sorted_bisection

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only : t_solver, allocate_pprts_solver_from_commandline

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, &
    destroy_tenstr_atm, print_tenstr_atm

  implicit none

contains
  subroutine example_rrtm_lw_sw(nxp, nyp, nzp, dx, dy)

    implicit none
    integer(iintegers), intent(in) :: nxp, nyp, nzp  ! local domain size for each rank
    real(ireals), intent(in) :: dx, dy               ! horizontal grid spacing in [m]

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y, ierr

    real(ireals) :: phi0, theta0, theta              ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals) :: lwc0                             ! cloud lwc content
    real(ireals),parameter :: albedo_th=0, albedo_sol=.3 ! broadband ground albedo for solar and thermal spectrum

    real(ireals), dimension(nzp+1,nxp,nyp), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1,nxp,nyp), target :: tlev ! Temperature on layer interfaces [K]

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp,nxp,nyp), target :: lwc, reliq

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
    character(len=default_str_len) :: atm_filename ! ='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev, icld(2)
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq

    logical :: lverbose, lthermal, lsolar, lrandom_lwc

    class(t_solver), allocatable :: pprts_solver
    type(t_tenstr_atm) :: atm

    integer(iintegers) :: solve_iterations, iter
    real(ireals) :: solve_iterations_scale
    logical :: lflg

    comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(comm, numnodes, ierr)
    call MPI_COMM_RANK(comm, myid, ierr)

    lverbose = .True.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-v", lverbose, lflg,ierr) ; call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "--verbose", lverbose, lflg,ierr) ; call CHKERR(ierr)

    if (myid.eq.0.and.lverbose) print *,'Running rrtm_lw_sw example with grid size:', nxp, nyp, nzp


    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if(N_ranks_y*N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    endif
    if(lverbose .and. myid.eq.0) print *, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate(nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 250 hPa and a
    ! Temperature difference of 60K
    do k=1,nzp+1
      plev(k,:,:) = linspace(k, [1e3_ireals, 250._ireals], nzp+1)
      tlev(k,:,:) = linspace(k, [290._ireals, 230._ireals], nzp+1)
    enddo

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
    reliq = 10

    icld(1) = nint(search_sorted_bisection(plev(:,1,1), 800._ireals))
    icld(2) = nint(search_sorted_bisection(plev(:,1,1), 700._ireals))

    lwc0 = 1e-2
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-lwc", lwc0, lflg, ierr)
    lwc  (icld(1):icld(2), :,:) = lwc0

    lrandom_lwc = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-random_lwc", lrandom_lwc, lflg,ierr) ; call CHKERR(ierr)
    if(lrandom_lwc) then
      call random_number(lwc(icld(1):icld(2), :,:))
      lwc(icld(1):icld(2), :,:) = lwc(icld(1):icld(2), :,:) * lwc0
    endif

    !tlev (icld  , :,:) = 288
    !tlev (icld+1, :,:) = tlev (icld  , :,:)

    if(myid.eq.0 .and. lverbose) print *,'Setup Atmosphere...'

    pplev(1:size(plev,1),1:size(plev,2)*size(plev,3)) => plev
    ptlev(1:size(tlev,1),1:size(tlev,2)*size(tlev,3)) => tlev
    plwc (1:size(lwc ,1),1:size(lwc ,2)*size(lwc ,3)) => lwc
    preliq(1:size(reliq,1),1:size(reliq,2)*size(reliq,3)) => reliq

    atm_filename='afglus_100m.dat'
    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-atm_filename', &
      atm_filename, lflg, ierr); call CHKERR(ierr)

    phi0 = 180
    theta0 = 0
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi0", phi0, lflg, ierr); call CHKERR(ierr)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta0", theta0, lflg, ierr); call CHKERR(ierr)

    lsolar = .True.
    lthermal = .True.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr); call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(pprts_solver, '3_10')

    solve_iterations = 1
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solve_iterations", &
      solve_iterations, lflg,ierr) ; call CHKERR(ierr)
    solve_iterations_scale = 0
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solve_iterations_scale", &
      solve_iterations_scale, lflg,ierr) ; call CHKERR(ierr)

    do iter = 1, solve_iterations
      if(lthermal) then
        ! Perform a circular shift of cloud field by one column
        lwc   = cshift(lwc  , 1, dim=2)
        reliq = cshift(reliq, 1, dim=2)
      endif
      if(lsolar) then
        theta = theta0 + real(iter-1, ireals) * solve_iterations_scale
      else
        theta = 0
      endif

      call setup_tenstr_atm(comm, .False., atm_filename, &
        pplev, ptlev, atm, &
        d_lwc=plwc, d_reliq=preliq)

      ! Hack the background atmosphere so that there is no jump in temperatures,
      ! i.e. add the difference at the glue location to all background temps
      atm%tlev(atm%d_ke1+1:size(atm%tlev,1),:) = atm%tlev(atm%d_ke1+1:size(atm%tlev,1),:) &
        + meanval(atm%tlev(atm%d_ke1,:)-atm%tlev(atm%d_ke1+1,:))
      if(iter.eq.1.and.myid.eq.0) call print_tenstr_atm(atm)

      if(lsolar.and.myid.eq.0) print *,'theta0 =', theta
      call pprts_rrtmg(comm, pprts_solver, atm, nxp, nyp, &
        dx, dy, phi0, theta,   &
        albedo_th, albedo_sol,  &
        lthermal, lsolar,       &
        edir, edn, eup, abso,   &
        nxproc=nxproc, nyproc=nyproc )

      nlev = ubound(edn,1)
      if(myid.eq.0 .and. lverbose) then
        do k=1,nlev
          if(allocated(edir)) then
            print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
          else
            print *,k, 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
          endif
        enddo

        if(allocated(edir)) &
          print *,'surface :: direct flux', meanval(edir(nlev,:,:))
        print *,'surface :: downw flux ', meanval(edn (nlev,:,:))
        print *,'surface :: upward fl  ', meanval(eup (nlev,:,:))
        print *,'surface :: absorption ', meanval(abso(nlev-1,:,:))

        if(allocated(edir)) &
          print *,'TOA :: direct flux', meanval(edir(1,:,:))
        print *,'TOA :: downw flux ', meanval(edn (1,:,:))
        print *,'TOA :: upward fl  ', meanval(eup (1,:,:))
        print *,'TOA :: absorption ', meanval(abso(1,:,:))

        if(allocated(edir)) &
          print *,'icloud :: direct flux  ', meanval(edir(nlev-icld,:,:))
        if(allocated(edir)) &
          print *,'icloud+1 :: direct flux', meanval(edir(nlev-icld+1,:,:))
        print *,'icloud :: downw flux   ', meanval(edn (nlev-icld+1,:,:))
        print *,'icloud :: upward fl    ', meanval(eup (nlev-icld  ,:,:))
        print *,'icloud :: absorption   ', meanval(abso(nlev-icld  ,:,:))
      endif
    enddo

    ! Tidy up
    call destroy_pprts_rrtmg(pprts_solver, lfinalizepetsc=.True.)
    call destroy_tenstr_atm(atm)
  end subroutine

end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only : iintegers, mpiint, ireals
  use m_example_pprts_rrtm_iterations, only: example_rrtm_lw_sw

  implicit none

  integer(mpiint) :: ierr, myid
  integer(iintegers) :: Nx, Ny, Nz
  real(ireals)       :: dx, dy
  logical :: lflg

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER ,ierr)

  Nx=3; Ny=3; Nz=5
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg, ierr)

  dx = 500
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr)
  dy = dx
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr)

  call example_rrtm_lw_sw(Nx, Ny, Nz, dx, dy)

  call mpi_finalize(ierr)
end program
