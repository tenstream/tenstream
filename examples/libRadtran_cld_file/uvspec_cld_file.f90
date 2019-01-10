module m_example_pprts_rrtm_lw_sw
  use mpi
  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only : t_solver, t_solver_3_10, t_solver_3_16
  use m_pprts, only: gather_all_toZero

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, default_str_len

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, hydrostat_dp

  use m_helper_functions, only: CHKERR, itoa
  use m_netcdfio, only: ncload, ncwrite, get_global_attribute

  implicit none

contains
  subroutine example_uvspec_cld_file(comm, cldfile, outfile, albedo_th, albedo_sol, lsolar, lthermal, phi0, theta0, Tsrfc, dTdz)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, outfile
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: Tsrfc, dTdz

    real(ireals), dimension(:,:,:), allocatable, target :: lwc, reliq ! will have global shape Nz, Nx, Ny
    real(ireals), dimension(:,:,:), allocatable, target :: plev, tlev ! will have local shape nzp+1, nxp, nyp
    real(ireals), dimension(:), allocatable :: hhl ! dim Nz+1
    character(len=default_str_len) :: groups(2)

    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    real(ireals),allocatable, dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf

    type(t_solver_3_10) :: pprts_solver

    real(ireals) :: dx, dy
    integer(mpiint) :: myid, numnodes, ierr
    integer(iintegers) :: N_ranks_x, N_ranks_y, iproc, jproc, is, ie, js, je
    integer(iintegers) :: nxp, nyp, nzp ! local sizes of domain, nzp being number of layers

    integer(iintegers) :: k

    call mpi_comm_rank(comm, myid, ierr)

    ! Load LibRadtran Cloud File
    call get_global_attribute(cldfile, 'dx', dx)
    call get_global_attribute(cldfile, 'dy', dy)
    groups(1) = trim(cldfile)
    groups(2) = trim('lwc'); call ncload(groups, lwc, ierr); call CHKERR(ierr)
    groups(2) = trim('reff'); call ncload(groups, reliq, ierr); call CHKERR(ierr)
    groups(2) = trim('z'); call ncload(groups, hhl, ierr); call CHKERR(ierr)

    dx  = dx  * 1e+3_ireals
    dy  = dy  * 1e+3_ireals
    hhl = hhl * 1e+3_ireals

    if(myid.eq.0) then
      print *,'Loaded LibRadtran Cloud File with:'
      print *,'dx, dy:', dx, dy
      print *,'hhl', hhl
      print *,'shape lwc ', shape(lwc)
      print *,'shape reliq', shape(reliq)
    endif

    ! Determine Domain Decomposition
    call mpi_comm_size(comm, numnodes, ierr)
    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if(N_ranks_y*N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    endif
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    nxp = size(lwc, dim=2) / N_ranks_x
    nyp = size(lwc, dim=3) / N_ranks_y
    call CHKERR(modulo(size(lwc, dim=2), N_ranks_x), &
      'x-dimension is not evenly distributable on given communicator!'// &
      'cant put'//itoa(size(lwc, dim=2))//' pixels on '//itoa(N_ranks_x)//' ranks')
    call CHKERR(modulo(size(lwc, dim=3), N_ranks_y), &
      'y-dimension is not evenly distributable on given communicator!'// &
      'cant put'//itoa(size(lwc, dim=3))//' pixels on '//itoa(N_ranks_y)//' ranks')

    if(myid.eq.0) then
      print *,'Local Domain sizes are:',nxp,nyp
    endif
    jproc = myid / N_ranks_x
    iproc = modulo(myid, N_ranks_x)

    js = 1 + jproc * nyp
    je = js + nyp -1

    is = 1 + (myid - jproc*N_ranks_x) * nxp
    ie = is + nxp -1

    call mpi_barrier(comm, ierr)
    print *,myid,'i,j proc', iproc, jproc,' local portion: ', is, ie, 'and', js, je

    ! Start with a dynamics grid starting at 1000 hPa with a specified lapse rate and surface temperature
    nzp = size(hhl)-1
    allocate(plev(nzp+1, nxp, nyp), tlev(nzp+1, nxp, nyp))
    plev(1,:,:) = 1000_ireals
    tlev(1,:,:) = Tsrfc
    do k=2,nzp+1
      tlev(k,:,:) = tlev(k-1,:,:) + dTdz * (hhl(k)-hhl(k-1))
      plev(k,:,:) = plev(k-1,:,:) - hydrostat_dp(hhl(k)-hhl(k-1), plev(k-1,:,:), (tlev(k,:,:)+tlev(k-1,:,:))/2)
    enddo

    if(myid.eq.0) then
      do k=1,nzp+1
        print *, k, 'plev', plev(k,is,js), 'Tlev', tlev(k,is,js)
      enddo
    endif

    call run_rrtmg_lw_sw(pprts_solver, dx, dy, phi0, theta0, &
      plev, tlev, &
      lwc(:, is:ie, js:je), reliq(:, is:ie, js:je), &
      albedo_th, albedo_sol, lsolar, lthermal, &
      edir, edn, eup, abso)

    groups(1) = trim(outfile)

    if(allocated(edir)) then
      call gather_all_toZero(pprts_solver%C_one_atm1, edir, gedir)
      if(myid.eq.0) then
        print *,'dumping direct radiation with local and global shape', shape(edir), ':', shape(gedir)
        groups(2) = 'edir'; call ncwrite(groups, gedir, ierr); call CHKERR(ierr)
      endif
    endif
    call gather_all_toZero(pprts_solver%C_one_atm1, edn, gedn)
    call gather_all_toZero(pprts_solver%C_one_atm1, eup, geup)
    call gather_all_toZero(pprts_solver%C_one_atm, abso, gabso)
    if(myid.eq.0) then
      print *,'dumping edn radiation with local and global shape', shape(edn), ':', shape(gedn)
      groups(2) = 'edn' ; call ncwrite(groups, gedn , ierr); call CHKERR(ierr)
      print *,'dumping eup radiation with local and global shape', shape(eup), ':', shape(geup)
      groups(2) = 'eup' ; call ncwrite(groups, geup , ierr); call CHKERR(ierr)
      print *,'dumping abso radiation with local and global shape', shape(abso), ':', shape(gabso)
      groups(2) = 'abso'; call ncwrite(groups, gabso, ierr); call CHKERR(ierr)
    endif

    call destroy_pprts_rrtmg(pprts_solver, lfinalizepetsc=.True.)
  end subroutine
  subroutine run_rrtmg_lw_sw(pprts_solver, dx, dy, phi0, theta0, &
      plev, tlev, lwc, reliq, albedo_th, albedo_sol, lsolar, lthermal, &
      edir, edn, eup, abso)
    type(t_solver_3_10) :: pprts_solver
    real(ireals),intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in), dimension(:,:,:), contiguous, target :: lwc, reliq ! dim(Nz,Nx,Ny)
    real(ireals), intent(in) :: albedo_th, albedo_sol ! broadband ground albedo for solar and thermal spectrum
    logical, intent(in) :: lsolar, lthermal ! switches if solar or thermal computations should be done
    real(ireals), dimension(:,:,:), contiguous, target, intent(in) :: plev ! pressure on layer interfaces [hPa]   dim=nzp+1,nxp,nyp
    real(ireals), dimension(:,:,:), contiguous, target, intent(in) :: tlev ! Temperature on layer interfaces [K]  dim=nzp+1,nxp,nyp

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y, mpierr


    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    ! real(ireals), dimension(nzp,nxp,nyp), target :: lwc, reliq

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
    character(len=default_str_len),parameter :: atm_filename='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq

    logical,parameter :: ldebug=.True.

    type(t_tenstr_atm) :: atm

    comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(comm, numnodes, mpierr)
    call MPI_COMM_RANK(comm, myid, mpierr)

    N_ranks_y = int(sqrt(1.*numnodes))
    N_ranks_x = numnodes / N_ranks_y
    if(N_ranks_y*N_ranks_x .ne. numnodes) then
      N_ranks_x = numnodes
      N_ranks_y = 1
    endif
    if(ldebug .and. myid.eq.0) print *, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y, '::', numnodes

    allocate(nxproc(N_ranks_x), source=size(plev,2)) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=size(plev,3)) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Not much going on in the dynamics grid, we actually don't supply trace
    ! gases to the TenStream solver... this will then be interpolated from the
    ! background profile (read from `atm_filename`)
    ! h2ovmr = zero
    ! o3vmr  = zero
    ! co2vmr = zero
    ! ch4vmr = zero
    ! n2ovmr = zero
    ! o2vmr  = zero

    if(myid.eq.0 .and. ldebug) print *,'Setup Atmosphere...'

    pplev(1:size(plev,1),1:size(plev,2)*size(plev,3)) => plev
    ptlev(1:size(tlev,1),1:size(tlev,2)*size(tlev,3)) => tlev
    plwc (1:size(lwc ,1),1:size(lwc ,2)*size(lwc ,3)) => lwc
    preliq(1:size(reliq,1),1:size(reliq,2)*size(reliq,3)) => reliq

    call setup_tenstr_atm(comm, .False., atm_filename, &
      pplev, ptlev, atm, &
      d_lwc=plwc, d_reliq=preliq)

    call pprts_rrtmg(comm, pprts_solver, atm, size(plev,2), size(plev,3), &
      dx, dy, phi0, theta0,                    &
      albedo_th, albedo_sol,                   &
      lthermal, lsolar,                        &
      edir, edn, eup, abso,                    &
      nxproc=nxproc, nyproc=nyproc, opt_time=zero)

    nlev = ubound(edn,1)
    if(myid.eq.0) then
      if(ldebug) then
        do k=1,nlev
          if(allocated(edir)) then
          print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
        else
          print *,k, 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
        endif
        enddo
      endif

      if(allocated(edir)) &
        print *,'surface :: direct flux', edir(nlev,1,1)
      print *,'surface :: downw flux ', edn (nlev,1,1)
      print *,'surface :: upward fl  ', eup (nlev,1,1)
      print *,'surface :: absorption ', abso(nlev-1,1,1)

      if(allocated(edir)) &
        print *,'TOA :: direct flux', edir(1,1,1)
      print *,'TOA :: downw flux ', edn (1,1,1)
      print *,'TOA :: upward fl  ', eup (1,1,1)
      print *,'TOA :: absorption ', abso(1,1,1)

    endif

    ! Tidy up
    call destroy_tenstr_atm(atm)
  end subroutine

end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only : iintegers, mpiint
  use m_example_pprts_rrtm_lw_sw

  implicit none

  integer(mpiint) :: ierr, myid
  logical :: lthermal, lsolar, lflg
  character(len=10*default_str_len) :: cldfile, outfile
  real(ireals) :: Ag, phi0, theta0, Tsrfc, dTdz

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)
  call init_mpi_data_parameters(mpi_comm_world)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-cld', cldfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a cloud filename... please call with -cld <libRadtran_cloud_file.nc>'

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a output filename... please call with -out <output.nc>'

  Ag = .1
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)

  lsolar = .True.
  lthermal = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr) ; call CHKERR(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr) ; call CHKERR(ierr)


  phi0 = 0
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg,ierr) ; call CHKERR(ierr)
  theta0 = 20
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg,ierr) ; call CHKERR(ierr)

  Tsrfc = 288
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Tsrfc", Tsrfc, lflg,ierr) ; call CHKERR(ierr)
  dTdz = -6.5_ireals * 1e-3_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dTdz", dTdz, lflg,ierr) ; call CHKERR(ierr)

  call example_uvspec_cld_file(mpi_comm_world, cldfile, outfile, zero, Ag, lsolar, lthermal, phi0, theta0, Tsrfc, dTdz)

  call mpi_finalize(ierr)
end program

