module m_example_uvspec_cld_file
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_pprts_base, only : t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: gather_all_toZero

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, &
    i0, i1, i2, zero, one, default_str_len

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, hydrostat_dp

  use m_helper_functions, only: CHKERR, itoa, imp_bcast, reverse, spherical_2_cartesian, resize_arr, &
    domain_decompose_2d
  use m_netcdfio, only: ncload, ncwrite, get_global_attribute

  use m_icon_plex_utils, only: create_2d_regular_plex, dmplex_2D_to_3D, &
    rank0_f90vec_to_plex, dmplex_gvec_from_f90_array, plex_gvec_tozero

  use m_plex_grid, only: t_plexgrid, &
    setup_plexgrid, create_plex_section

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

  use m_plex_rt, only: init_plex_rt_solver

  use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

  implicit none

contains
  subroutine example_uvspec_cld_file(comm, &
      cldfile, atm_filename, outfile, &
      albedo_th, albedo_sol, &
      lsolar, lthermal, &
      phi0, theta0, &
      Tsrfc, dTdz)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, atm_filename, outfile
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

    class(t_solver), allocatable :: pprts_solver

    real(ireals) :: dx, dy
    integer(mpiint) :: myid, N_ranks_x, N_ranks_y, ierr
    integer(iintegers) :: iproc, jproc, is, ie, js, je
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

    if(size(lwc  ,dim=2).eq.1) call resize_arr(3_iintegers, lwc  , dim=2, lrepeat=.True.)
    if(size(reliq,dim=2).eq.1) call resize_arr(3_iintegers, reliq, dim=2, lrepeat=.True.)

    if(size(lwc  ,dim=3).eq.1) call resize_arr(3_iintegers, lwc  , dim=3, lrepeat=.True.)
    if(size(reliq,dim=3).eq.1) call resize_arr(3_iintegers, reliq, dim=3, lrepeat=.True.)

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
    call domain_decompose_2d(comm, N_ranks_x, N_ranks_y, ierr); call CHKERR(ierr)
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y

    nxp = size(lwc, dim=2) / N_ranks_x
    nyp = size(lwc, dim=3) / N_ranks_y
    call CHKERR(modulo(size(lwc, dim=2,kind=mpiint), int(N_ranks_x,mpiint)), &
      'x-dimension is not evenly distributable on given communicator!'// &
      'cant put'//itoa(size(lwc, dim=2))//' pixels on '//itoa(N_ranks_x)//' ranks')
    call CHKERR(modulo(size(lwc, dim=3, kind=mpiint), int(N_ranks_y, mpiint)), &
      'y-dimension is not evenly distributable on given communicator!'// &
      'cant put'//itoa(size(lwc, dim=3))//' pixels on '//itoa(N_ranks_y)//' ranks')

    if(myid.eq.0) then
      print *,'Local Domain sizes are:',nxp,nyp
    endif
    jproc = myid / N_ranks_x
    iproc = modulo(myid, int(N_ranks_x, mpiint))

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

    call allocate_pprts_solver_from_commandline(pprts_solver, default_solver='8_10')

    call run_rrtmg_lw_sw(pprts_solver, atm_filename, dx, dy, phi0, theta0, &
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

  subroutine run_rrtmg_lw_sw(pprts_solver, atm_filename, dx, dy, phi0, theta0, &
      plev, tlev, lwc, reliq, albedo_th, albedo_sol, lsolar, lthermal, &
      edir, edn, eup, abso)
    class(t_solver) :: pprts_solver
    real(ireals),intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in), dimension(:,:,:), contiguous, target :: lwc, reliq ! dim(Nz,Nx,Ny)
    real(ireals), intent(in) :: albedo_th, albedo_sol ! broadband ground albedo for solar and thermal spectrum
    logical, intent(in) :: lsolar, lthermal ! switches if solar or thermal computations should be done
    real(ireals), dimension(:,:,:), contiguous, target, intent(in) :: plev ! pressure on layer interfaces [hPa]   dim=nzp+1,nxp,nyp
    real(ireals), dimension(:,:,:), contiguous, target, intent(in) :: tlev ! Temperature on layer interfaces [K]  dim=nzp+1,nxp,nyp

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: comm, myid, N_ranks_x, N_ranks_y, ierr


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
    character(len=*), intent(in) :: atm_filename

    !------------ Local vars ------------------
    integer(iintegers) :: k, nlev
    integer(iintegers),allocatable :: nxproc(:), nyproc(:)

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq

    logical,parameter :: ldebug=.True.

    type(t_tenstr_atm) :: atm

    comm = mpi_comm_world
    call mpi_comm_rank(comm, myid, ierr)

    ! Determine Domain Decomposition
    call domain_decompose_2d(comm, N_ranks_x, N_ranks_y, ierr); call CHKERR(ierr)
    if(myid.eq.0) print *, myid, 'Domain Decomposition will be', N_ranks_x, 'and', N_ranks_y

    allocate(nxproc(N_ranks_x), source=size(plev,2, kind=iintegers)) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=size(plev,3, kind=iintegers)) ! values have to define the local domain sizes on each rank (here constant on all processes)

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

    call pprts_rrtmg(comm, pprts_solver, atm, &
      size(plev,2, kind=iintegers), size(plev,3, kind=iintegers), &
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

  subroutine example_uvspec_cld_file_with_plexrt(comm, &
      cldfile, atm_filename, outfile, &
      albedo_th, albedo_sol, &
      lsolar, lthermal, &
      phi0, theta0, &
      Tsrfc, dTdz)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, atm_filename, outfile
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: Tsrfc, dTdz

    real(ireals), dimension(:,:,:), allocatable, target :: glob_lwc, glob_reliq ! will have global shape Nz, Nx, Ny
    real(ireals), dimension(:,:), allocatable, target :: plev, tlev ! will have local shape nzp+1, Ncol
    real(ireals), dimension(:), allocatable :: hhl ! dim Nz+1
    character(len=default_str_len) :: groups(2)

    real(ireals),allocatable, dimension(:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    !real(ireals),allocatable, dimension(:,:) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf

    real(ireals) :: dx, dy
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: k

    integer(iintegers) :: Nx_global, Ny_global, Nz
    integer(iintegers) :: fStart, fEnd, Ncol, Nlev

    type(tDM) :: dm2d, dm2d_dist, dm3d
    type(tPetscSF) :: migration_sf
    type(tPetscSection) :: parCellSection

    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)
    class(t_plex_solver), allocatable :: solver
    type(t_tenstr_atm) :: atm
    real(ireals), pointer :: xlwc(:), xreff(:), col_lwc(:,:), col_reff(:,:)
    type(tVec) :: vlwc, vreff

    real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system

    call mpi_comm_rank(comm, myid, ierr)

    if(myid.eq.0) then
      ! Load LibRadtran Cloud File
      call get_global_attribute(cldfile, 'dx', dx)
      call get_global_attribute(cldfile, 'dy', dy)
      groups(1) = trim(cldfile)
      groups(2) = trim('lwc') ; call ncload(groups, glob_lwc,   ierr); call CHKERR(ierr)
      groups(2) = trim('reff'); call ncload(groups, glob_reliq, ierr); call CHKERR(ierr)
      groups(2) = trim('z')   ; call ncload(groups, hhl,        ierr); call CHKERR(ierr)

      if(size(glob_lwc  ,dim=2).eq.1) call resize_arr(3_iintegers, glob_lwc  , dim=2, lrepeat=.True.)
      print *,'shape lwc', shape(glob_lwc)
      if(size(glob_reliq,dim=2).eq.1) call resize_arr(3_iintegers, glob_reliq, dim=2, lrepeat=.True.)
      print *,'shape lwc', shape(glob_lwc)

      if(size(glob_lwc  ,dim=3).eq.1) call resize_arr(3_iintegers, glob_lwc  , dim=3, lrepeat=.True.)
      print *,'shape lwc', shape(glob_lwc)
      if(size(glob_reliq,dim=3).eq.1) call resize_arr(3_iintegers, glob_reliq, dim=3, lrepeat=.True.)
      print *,'shape lwc', shape(glob_lwc)

      dx  = dx  * 1e+3_ireals
      dy  = dy  * 1e+3_ireals
      hhl = hhl * 1e+3_ireals

      Nz        = size(glob_lwc,dim=1)
      Nx_global = size(glob_lwc,dim=2)
      Ny_global = size(glob_lwc,dim=3)

      print *,'Loaded LibRadtran Cloud File with:'
      print *,'dx, dy:', dx, dy
      print *,'global Nx/Ny', Nx_global, Ny_global
      print *,'hhl', hhl
      print *,'shape lwc ', shape(glob_lwc)
      print *,'shape reliq', shape(glob_reliq)
    endif

    call imp_bcast(comm, Nx_global, 0_mpiint)
    call imp_bcast(comm, Ny_global, 0_mpiint)
    call imp_bcast(comm, Nz, 0_mpiint)
    call imp_bcast(comm, dx, 0_mpiint)
    call imp_bcast(comm, dy, 0_mpiint)
    call imp_bcast(comm, hhl, 0_mpiint)

    call create_2d_regular_plex(comm, Nx_global+1, Ny_global+1, dm2d, dm2d_dist, &
      opt_migration_sf=migration_sf, opt_dx=dx)

    call DMPlexGetHeightStratum(dm2d_dist, i0, fStart, fEnd, ierr); call CHKERR(ierr)
    Ncol = fEnd - fStart

    !if(myid.eq.0) then
      print *,myid,'Local Domain sizes are:',Ncol,Nz
    !endif

    ! Start with a dynamics grid starting at 1000 hPa with a specified lapse rate and surface temperature
    allocate(plev(Nz+1, fStart:fEnd-1), tlev(Nz+1, fStart:fEnd-1))
    plev(1,:) = 1000_ireals
    tlev(1,:) = Tsrfc
    do k=2,Nz+1
      tlev(k,:) = tlev(k-1,:) + dTdz * (hhl(k)-hhl(k-1))
      plev(k,:) = plev(k-1,:) - hydrostat_dp(hhl(k)-hhl(k-1), plev(k-1,:), (tlev(k,:)+tlev(k-1,:))/2)
    enddo

    if(myid.eq.0) then
      do k=1,Nz+1
        print *, k, 'plev', plev(k,fStart), 'Tlev', tlev(k,fStart)
      enddo
    endif
    hhl = reverse(hhl)

    call setup_tenstr_atm(comm, .False., atm_filename, &
      plev, tlev, atm)

    Nlev = size(atm%plev,1,kind=iintegers)
    call dmplex_2D_to_3D(dm2d_dist, Nlev, reverse(atm%zt(:, i1)), dm3d, zindex, lpolar_coords=.False.)

    call setup_plexgrid(dm2d_dist, dm3d, Nlev-1, zindex, plex, reverse(atm%zt(:, i1)))
    deallocate(zindex)

    if(myid.eq.0) then
      ! increase the size by two in the x direction,
      ! i.e. from boxes to wedges with 2 elements per box
      call resize_arr(size(glob_lwc  ,2,kind=iintegers)*2, glob_lwc  , dim=2, fillVal=-1._ireals)
      call resize_arr(size(glob_reliq,2,kind=iintegers)*2, glob_reliq, dim=2, fillVal=-1._ireals)
      glob_lwc  (:, 1:size(glob_lwc  ,dim=2):2, :) = glob_lwc  (:, 1:size(glob_lwc  ,dim=2)/2, :)
      glob_lwc  (:, 2:size(glob_lwc  ,dim=2):2, :) = glob_lwc  (:, 1:size(glob_lwc  ,dim=2):2, :)
      glob_reliq(:, 1:size(glob_reliq,dim=2):2, :) = glob_reliq(:, 1:size(glob_reliq,dim=2)/2, :)
      glob_reliq(:, 2:size(glob_reliq,dim=2):2, :) = glob_reliq(:, 1:size(glob_reliq,dim=2):2, :)

      col_lwc (1:size(glob_lwc ,1),1:size(glob_lwc ,2)*size(glob_lwc ,3)) => glob_lwc
      col_reff(1:size(glob_lwc ,1),1:size(glob_lwc ,2)*size(glob_lwc ,3)) => glob_reliq
      print *,myid,'global shape col lwc', shape(col_lwc)
    endif

    call rank0_f90vec_to_plex(dm2d, dm2d_dist, migration_sf, col_lwc , parCellSection, vlwc)
    call rank0_f90vec_to_plex(dm2d, dm2d_dist, migration_sf, col_reff, parCellSection, vreff)
    nullify(col_lwc)
    nullify(col_reff)

    call VecGetArrayReadF90(vlwc , xlwc , ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(vreff, xreff, ierr); call CHKERR(ierr)

    col_lwc (1:Nz, fstart:fEnd-1) => xlwc
    col_reff(1:Nz, fstart:fEnd-1) => xreff

    call setup_tenstr_atm(comm, .False., atm_filename, &
      plev, tlev, atm, &
      d_lwc=col_lwc, d_reliq=col_reff)

    nullify(col_lwc)
    nullify(col_reff)
    call VecRestoreArrayReadF90(vlwc , xlwc , ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(vreff, xreff, ierr); call CHKERR(ierr)

    ! Finished preparing input, lets do the computations

    call allocate_plexrt_solver_from_commandline(solver, '5_8')
    call init_plex_rt_solver(plex, solver)

    sundir = spherical_2_cartesian(phi0,theta0)
    print *,'sundir', sundir

    call plexrt_rrtmg(solver, atm, sundir, &
      albedo_thermal=albedo_th, albedo_solar=albedo_sol, &
      lthermal=lthermal, lsolar=lsolar, &
      edir=edir, edn=edn, eup=eup, abso=abso)

    call dump_results()

    call DMDestroy(dm2d, ierr); call CHKERR(ierr)
    call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)

    call destroy_plexrt_rrtmg(solver, lfinalizepetsc=.False.)
    contains
      subroutine transfer_arr(col_inp, arr3d, idof_start, idof_end)
        real(ireals), intent(in) :: col_inp(:,:)
        real(ireals), intent(out) :: arr3d(:,:,:)
        integer(iintegers), intent(in) :: idof_start, idof_end
        integer(iintegers) :: i, j, icol, idof
        do j=1,Ny_global
          do i=1,Nx_global
            icol = Nx_global*(j-1) + i
            arr3d(:,i,j) = 0
            do idof = idof_start, idof_end
              arr3d(:,i,j) = arr3d(:,i,j) + col_inp(:,2*(icol-1)+idof)
            enddo
            arr3d(:,i,j) = arr3d(:,i,j) / real(idof_end-idof_start+1, ireals)
          enddo
        enddo
      end subroutine
      subroutine dump_var(var, varname)
        real(ireals), allocatable, intent(in) :: var(:,:)
        character(len=*), intent(in) :: varname
        type(tPetscSection) :: flxSection, r0flxSection
        type(tVec) :: v_var
        type(tVec) :: r0var

        real(ireals), pointer :: xarr(:), xxarr(:,:)
        real(ireals), dimension(:,:,:), allocatable :: oarr

        groups(1) = trim(outfile)
        if(myid.eq.0) allocate(oarr (size(var,dim=1), Nx_global, Ny_global))

        if(.not. allocated(var)) return

        call create_plex_section(dm2d_dist, 'face_section', i1, &
          [i0], [size(var,dim=1,kind=iintegers)], [i0], [i0], flxSection)

        call dmplex_gVec_from_f90_array(comm, var, v_var)
        call plex_gVec_toZero(dm2d_dist, migration_sf, flxSection, v_var, &
          r0flxSection, r0var)

        call VecGetArrayF90(r0var, xarr,ierr); call CHKERR(ierr)
        if(myid.eq.0) then
          xxarr(1:size(var,dim=1), 1:Nx_global*Ny_global*2) => xarr

          call transfer_arr(xxarr, oarr, i1, i1)
          groups(2) = trim(varname)//'_leftvertices'; call ncwrite(groups, oarr, ierr); call CHKERR(ierr)
          call transfer_arr(xxarr, oarr, i2, i2)
          groups(2) = trim(varname)//'_rightvertices'; call ncwrite(groups, oarr, ierr); call chkerr(ierr)
          call transfer_arr(xxarr, oarr, i1, i2)
          groups(2) = trim(varname) ; call ncwrite(groups, oarr, ierr); call chkerr(ierr)

          nullify(xxarr)
        endif
        call VecRestoreArrayF90(r0var, xarr, ierr); call CHKERR(ierr)

        call VecDestroy(v_var, ierr); call CHKERR(ierr)
        call VecDestroy(r0var, ierr); call CHKERR(ierr)
      end subroutine
      subroutine dump_results()

        call dump_var(edir, 'edir')
        call dump_var(edn , 'edn')
        call dump_var(eup , 'eup')
        call dump_var(abso, 'abso')

        if(myid.eq.0) then
          groups(1) = trim(outfile)
          groups(2) = 'lwc'  ; call ncwrite(groups, glob_lwc  , ierr); call CHKERR(ierr)
          groups(2) = 'reff' ; call ncwrite(groups, glob_reliq, ierr); call CHKERR(ierr)
        endif
      end subroutine
  end subroutine

end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only : mpiint
  use m_example_uvspec_cld_file

  implicit none

  integer(mpiint) :: ierr, myid
  logical :: lthermal, lsolar, lflg
  character(len=10*default_str_len) :: cldfile, outfile
  real(ireals) :: Ag, phi0, theta0, Tsrfc, dTdz
  character(len=default_str_len) :: atm_filename

  logical :: luse_plexrt

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-cld', cldfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply a cloud filename... please call with -cld <libRadtran_cloud_file.nc>')

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply a output filename... please call with -out <output.nc>')

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-atm_filename', atm_filename, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply an atmosphere filename... please call with -atm_filename <atm.dat>')

  Ag = .1
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg,ierr) ; call CHKERR(ierr)

  lsolar = .True.
  lthermal = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr) ; call CHKERR(ierr)
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr) ; call CHKERR(ierr)


  phi0 = 270
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg,ierr) ; call CHKERR(ierr)
  theta0 = 60
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg,ierr) ; call CHKERR(ierr)

  Tsrfc = 288
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Tsrfc", Tsrfc, lflg,ierr) ; call CHKERR(ierr)
  dTdz = -6.5_ireals * 1e-3_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dTdz", dTdz, lflg,ierr) ; call CHKERR(ierr)

  luse_plexrt = .False.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-use_plexrt", luse_plexrt, lflg,ierr) ; call CHKERR(ierr)

  if(luse_plexrt) then
    call example_uvspec_cld_file_with_plexrt(mpi_comm_world, cldfile, atm_filename, outfile, &
      zero, Ag, lsolar, lthermal, phi0, theta0, Tsrfc, dTdz)
  else
    call example_uvspec_cld_file(mpi_comm_world, cldfile, atm_filename, outfile, &
      zero, Ag, lsolar, lthermal, phi0, theta0, Tsrfc, dTdz)
  endif
  call PetscFinalize(ierr)
  call mpi_finalize(ierr)
end program
