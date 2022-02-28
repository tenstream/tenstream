module m_wetterstein
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, default_str_len

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only: &
    & CHKERR, &
    & domain_decompose_2d_petsc, &
    & imp_bcast, &
    & reverse, &
    & spherical_2_cartesian

  use m_petsc_helpers, only: &
    & f90VecToPetsc, &
    & getvecpointer, &
    & petscGlobalVecToZero, &
    & petscVecToF90, &
    & restorevecpointer

  use m_netcdfIO, only: ncwrite, ncload, set_global_attribute

  use m_pprts_rrtmg, only: pprts_rrtmg, destroy_pprts_rrtmg
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm
  use m_pprts, only: gather_all_toZero

  implicit none

contains
  subroutine ex_wetterstein(comm, &
      & input_filename, atm_filename, output_filename, &
      & phi0, theta0, albedo_sol, albedo_th)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: input_filename
    character(len=*), intent(in) :: atm_filename
    character(len=*), intent(in) :: output_filename
    real(ireals), intent(in) :: phi0, theta0
    real(ireals), intent(in) :: albedo_sol, albedo_th

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, myid

    integer(iintegers) :: Nx, Ny, Nlay   ! global domain size
    integer(iintegers) :: nxp, nyp, xs, ys ! local domain size in x and y aswell as start indices
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    class(t_solver), allocatable :: solver
    type(t_tenstr_atm), target :: atm

    real(ireals), parameter :: dx = 100, dy = dx

    real(ireals), allocatable, dimension(:, :, :) :: glob_plev, glob_tlev  ! nlay+1, Nx, Ny
    real(ireals), allocatable, dimension(:, :, :) :: glob_lwc              ! nlay, Nx, Ny
    real(ireals), allocatable, dimension(:, :, :), target :: plev, tlev   ! nlay+1, nxp, nyp
    real(ireals), allocatable, dimension(:, :, :), target :: lwc, reliq    ! nlay  , nxp, nyp
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev, plwc, preliq ! reshape pointers to convert to column vecs
    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! nlyr(+1), global_nx, global_ny
    !real(ireals),allocatable, dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf
    real(ireals), pointer, dimension(:, :, :) :: patm
    real(ireals), pointer :: hhl(:, :, :, :) => null(), hhl1d(:) => null()

    character(len=default_str_len) :: nc_path(2) ! [ filename, varname ]

    real(ireals) :: sundir(3)
    integer(iintegers) :: k
    integer(mpiint) :: ncerr, ierr

    logical, parameter :: lthermal = .false., lsolar = .true.
    logical, parameter :: ldebug = .true.

    call init_mpi_data_parameters(comm)
    call mpi_comm_size(comm, numnodes, ierr)
    call mpi_comm_rank(comm, myid, ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      nc_path(1) = trim(input_filename)
      nc_path(2) = 'plev'; call ncload(nc_path, glob_plev, ncerr); call CHKERR(ncerr)
      nc_path(2) = 'tlev'; call ncload(nc_path, glob_tlev, ncerr); call CHKERR(ncerr)
      nc_path(2) = 'lwc'; call ncload(nc_path, glob_lwc, ncerr); call CHKERR(ncerr)

      if (myid .eq. 0) print *, 'plev shape', shape(glob_plev)
    end if
    call imp_bcast(comm, glob_plev, 0_mpiint)
    call imp_bcast(comm, glob_tlev, 0_mpiint)
    call imp_bcast(comm, glob_lwc, 0_mpiint)

    Nlay = ubound(glob_plev, 1) - 1
    Nx = ubound(glob_plev, 2)
    Ny = ubound(glob_plev, 3)

    call domain_decompose_2d_petsc(comm, Nx, Ny, &
      & nxp, nyp, xs, ys, nxproc, nyproc, ierr); call CHKERR(ierr)
    if (myid .eq. 0) print *, myid, 'Domain Decomposition on', numnodes, 'ranks will be', nxproc, 'and', nyproc

    allocate (plev(nlay + 1, nxp, nyp))
    allocate (tlev(nlay + 1, nxp, nyp))
    allocate (lwc(nlay, nxp, nyp))
    plev = glob_plev(:, 1 + xs:xs + nxp, 1 + ys:ys + nyp)
    tlev = glob_tlev(:, 1 + xs:xs + nxp, 1 + ys:ys + nyp)
    lwc = glob_lwc(:, 1 + xs:xs + nxp, 1 + ys:ys + nyp)

    allocate (reliq(nlay, nxp, nyp))
    reliq = 10

    if (myid .eq. 0 .and. ldebug) then
      do k = 1, nlay + 1
        print *, k, 'plev', plev(k, 1, 1), 'Temp', tlev(k, 1, 1)
      end do
    end if

    sundir = spherical_2_cartesian(phi0, theta0)

    pplev(1:Nlay + 1, 1:nxp * nyp) => plev
    ptlev(1:Nlay + 1, 1:nxp * nyp) => tlev
    plwc(1:Nlay, 1:nxp * nyp) => lwc
    preliq(1:Nlay, 1:nxp * nyp) => reliq

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          pplev, ptlev, atm, &
                          d_lwc=plwc, d_reliq=preliq)

    call pprts_rrtmg(comm, solver, atm, nxp, nyp, &
                     dx, dy, sundir, &
                     albedo_th, albedo_sol, &
                     lthermal, lsolar, &
                     edir, edn, eup, abso, &
                     nxproc=nxproc, nyproc=nyproc)

    !if(allocated(edir)) then
    !  call gather_all_toZero(solver%C_one_atm1, edir, gedir)
    !  if(myid.eq.0) then
    !    print *,'dumping direct radiation with local and global shape', shape(edir), ':', shape(gedir)
    !    nc_path(2) = 'edir'; call ncwrite(nc_path, gedir, ierr); call CHKERR(ierr)
    !  endif
    !endif
    !call gather_all_toZero(solver%C_one_atm1, edn , gedn)
    !call gather_all_toZero(solver%C_one_atm1, eup , geup)
    !call gather_all_toZero(solver%C_one_atm , abso, gabso)

    !if(myid.eq.0) then
    !  nc_path(2) = 'edn' ; call ncwrite(nc_path, gedn , ierr); call CHKERR(ierr)
    !  nc_path(2) = 'eup' ; call ncwrite(nc_path, geup , ierr); call CHKERR(ierr)
    !  nc_path(2) = 'abso'; call ncwrite(nc_path, gabso, ierr); call CHKERR(ierr)
    !endif

    nc_path(1) = trim(output_filename)

    associate (&
        & C => solver%C_one,     &
        & C1 => solver%C_one1,    &
        & Ca => solver%C_one_atm, &
        & Ca1 => solver%C_one_atm1_box, &
        & Cs => solver%Csrfc_one)

      call dump_vec(Ca%da, solver%atm%dz, 'dz')

      patm(Ca1%zs:Ca1%ze, Ca1%xs:Ca1%xe, Ca1%ys:Ca1%ye) => atm%plev
      call dump_vec(Ca1%da, reverse(patm), 'p_lev')

      patm(Ca1%zs:Ca1%ze, Ca1%xs:Ca1%xe, Ca1%ys:Ca1%ye) => atm%tlev
      call dump_vec(Ca1%da, reverse(patm), 't_lev')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%tlay
      call dump_vec(Ca%da, reverse(patm), 't_lay')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%o3_lay
      call dump_vec(Ca%da, reverse(patm), 'o3_lay')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%o2_lay
      call dump_vec(Ca%da, reverse(patm), 'o2_lay')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%h2o_lay
      call dump_vec(Ca%da, reverse(patm), 'h2o_lay')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%co2_lay
      call dump_vec(Ca%da, reverse(patm), 'co2_lay')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%n2o_lay
      call dump_vec(Ca%da, reverse(patm), 'no2_lay')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%lwc
      call dump_vec(Ca%da, reverse(patm), 'lwc')

      patm(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%reliq
      call dump_vec(Ca%da, reverse(patm), 'reliq')

      if (allocated(edir)) &
        call dump_vec(C1%da, edir, 'edir')
      call dump_vec(C1%da, edn, 'edn')
      call dump_vec(C1%da, eup, 'eup')
      call dump_vec(C%da, abso, 'abso')

      call getVecPointer(Ca1%da, solver%atm%hhl, hhl1d, hhl)
      call dump_vec(Ca1%da, hhl(0, Ca1%zs:Ca1%ze, Ca1%xs:Ca1%xe, Ca1%ys:Ca1%ye), 'hhl')
      call dump_vec_2d(Cs%da, hhl(0, Ca1%ze, Ca1%xs:Ca1%xe, Ca1%ys:Ca1%ye), 'h_srfc')
      call restoreVecPointer(Ca1%da, solver%atm%hhl, hhl1d, hhl)

      call dump_vec_2d(Cs%da, edir(size(edir, 1), :, :), 'edir_srfc')
      call dump_vec_2d(Cs%da, edn(size(edn, 1), :, :), 'edn_srfc')
      call dump_vec_2d(Cs%da, eup(size(eup, 1), :, :), 'eup_srfc')

      if (myid .eq. 0) then
        call set_global_attribute(nc_path(1), 'Nx', C%glob_xm, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'Ny', C%glob_ym, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'Nlay', Nlay, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'dx', dx, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'dy', dy, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'phi0', phi0, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'theta0', theta0, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'Ag_solar', albedo_sol, ierr); call CHKERR(ierr)
        call set_global_attribute(nc_path(1), 'Ag_thermal', albedo_th, ierr); call CHKERR(ierr)
      end if
    end associate

    call destroy_pprts_rrtmg(solver, lfinalizepetsc=.true.)

  contains
    subroutine dump_vec_2d(dm, arr, varname)
      type(tDM), intent(in) :: dm
      real(ireals), intent(in) :: arr(:, :)
      character(len=*), intent(in) :: varname
      type(tVec) :: gvec, lVec
      real(ireals), allocatable :: larr(:, :)
      integer(mpiint) :: ierr

      call DMGetGlobalVector(dm, gvec, ierr); call CHKERR(ierr)
      call f90VecToPetsc(arr, dm, gvec)
      call petscGlobalVecToZero(gvec, dm, lVec)
      if (myid .eq. 0) then
        call petscVecToF90(lVec, dm, larr, only_on_rank0=.true.)

        nc_path(2) = trim(varname)
        call ncwrite(nc_path, larr, ierr); call CHKERR(ierr)
      end if
      call VecDestroy(lVec, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(dm, gvec, ierr); call CHKERR(ierr)
    end subroutine
    subroutine dump_vec(dm, arr, varname)
      type(tDM), intent(in) :: dm
      real(ireals), intent(in) :: arr(:, :, :)
      character(len=*), intent(in) :: varname
      type(tVec) :: gvec, lVec
      real(ireals), allocatable :: larr(:, :, :)
      integer(mpiint) :: ierr

      call DMGetGlobalVector(dm, gvec, ierr); call CHKERR(ierr)
      call f90VecToPetsc(arr, dm, gvec)
      call petscGlobalVecToZero(gvec, dm, lVec)
      if (myid .eq. 0) then
        call petscVecToF90(lVec, dm, larr, only_on_rank0=.true.)

        nc_path(2) = trim(varname)
        call ncwrite(nc_path, larr, ierr); call CHKERR(ierr)
      end if
      call VecDestroy(lVec, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(dm, gvec, ierr); call CHKERR(ierr)
    end subroutine

  end subroutine
end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi, only: mpi_init, mpi_finalize, MPI_COMM_WORLD
  use m_data_parameters, only: &
    & default_str_len, &
    & init_mpi_data_parameters, &
    & ireals, &
    & mpiint
  use m_helper_functions, only: CHKERR, get_petsc_opt
  use m_wetterstein, only: ex_wetterstein

  implicit none

  character(len=default_str_len) :: atm_filename
  character(len=default_str_len) :: input_filename
  character(len=default_str_len) :: output_filename

  real(ireals) :: phi0, theta0
  real(ireals) :: albedo_sol, albedo_th

  logical :: lflg
  integer(mpiint) :: ierr

  call mpi_init(ierr)
  call init_mpi_data_parameters(MPI_COMM_WORLD)

  input_filename = 'input.nc'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-inp', &
                     input_filename, lflg, ierr); call CHKERR(ierr)

  output_filename = 'output.nc'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-out', &
                     output_filename, lflg, ierr); call CHKERR(ierr)

  atm_filename = 'afglus_100m.dat'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', &
                     atm_filename, lflg, ierr); call CHKERR(ierr)

  phi0 = 180
  theta0 = 40
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-phi', &
    & phi0, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-theta', &
    & theta0, lflg, ierr); call CHKERR(ierr)

  albedo_sol = 0.2
  albedo_th = 0.05
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-albedo_sol', &
    & albedo_sol, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-albedo_th', &
    & albedo_th, lflg, ierr); call CHKERR(ierr)

  call ex_wetterstein( &
    & MPI_COMM_WORLD, &
    & input_filename, &
    & atm_filename, &
    & output_filename, &
    & phi0, theta0, &
    & albedo_sol, &
    & albedo_th)

  call mpi_finalize(ierr)
end program
