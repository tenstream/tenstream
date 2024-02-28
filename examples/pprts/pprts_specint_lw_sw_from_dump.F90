module m_example_pprts_specint_lw_sw_from_dump

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
  use m_pprts, only: gather_all_toZero

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy, load_input_dump

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, abso2hr, print_tenstr_atm

  use m_petsc_helpers, only: getvecpointer, restorevecpointer
  use m_netcdfio, only: ncwrite, set_attribute

  implicit none

contains
  subroutine ex_pprts_specint_lw_sw_from_dump(specint, comm, inpfile, outfile)
    character(len=*), intent(in) :: specint           ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: inpfile
    character(len=*), intent(in) :: outfile

    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso, hr ! [nlev_merged(-1), nxp, nyp]

    integer(mpiint) :: myid, ierr

    class(t_solver), allocatable :: pprts_solver

    integer(iintegers) :: Nx_local, Ny_local
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)
    real(ireals) :: dx, dy
    real(ireals), allocatable :: sundir(:)
    real(ireals) :: albedo_solar, albedo_thermal
    logical :: lsolar, lthermal
    type(t_tenstr_atm) :: atm
    real(ireals), allocatable :: opt_time, solar_albedo_2d(:, :), thermal_albedo_2d(:, :), opt_solar_constant

    integer(iintegers) :: icollapse, k, nlev
    logical :: lflg

    call init_mpi_data_parameters(comm)

    call MPI_COMM_RANK(comm, myid, ierr)

    call load_input_dump(&
      & comm, &
      & inpfile, &
      & Nx_local, Ny_local, &
      & nxproc, nyproc, &
      & dx, dy, &
      & sundir, &
      & albedo_thermal, albedo_solar, &
      & lsolar, lthermal, &
      & atm, &
      & opt_time, &
      & solar_albedo_2d, thermal_albedo_2d, &
      & opt_solar_constant, &
      & ierr)
    call CHKERR(ierr)

    call print_tenstr_atm(atm)

!    sundir = spherical_2_cartesian(phi0, theta0)
!
    call allocate_pprts_solver_from_commandline(pprts_solver, '3_10', ierr); call CHKERR(ierr)

    icollapse = 1
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-icollapse", icollapse, lflg, ierr)

    if (.not. allocated(opt_time)) then
      allocate (opt_time)
      opt_time = 0._ireals
    end if
    if (.not. allocated(solar_albedo_2d)) then
      allocate (solar_albedo_2d(Nx_local, Ny_local), source=albedo_solar)
    end if
    if (.not. allocated(thermal_albedo_2d)) then
      allocate (thermal_albedo_2d(Nx_local, Ny_local), source=albedo_thermal)
    end if
    if (allocated(opt_solar_constant)) call CHKERR(1_mpiint, 'use of opt_solar_constant not implemented') ! Easy to add here by hand    but cumbersome to implement generically in Fortran

    call specint_pprts(specint, comm, pprts_solver, atm, &
                       Nx_local, Ny_local, &
                       dx, dy, sundir, &
                       albedo_thermal, albedo_solar, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, &
                       icollapse=icollapse, &
                       opt_time=opt_time, &
                       solar_albedo_2d=solar_albedo_2d, &
                       thermal_albedo_2d=thermal_albedo_2d)

    allocate (hr(size(abso, 1), size(abso, 2), size(abso, 3)))
    call abso2hr(atm, abso, hr, ierr); call CHKERR(ierr)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      do k = 1, nlev
        if (allocated(edir)) then
          print *, k, 'edir', meanval(edir(k, :, :)), 'edn', meanval(edn(k, :, :)), 'eup', meanval(eup(k, :, :)), &
            & 'abso', meanval(abso(min(nlev - 1, k), :, :)), 'hr', meanval(hr(min(nlev - 1, k), :, :)) * 3600 * 24
        else
          print *, k, 'edir', 0, 'edn', meanval(edn(k, :, :)), 'eup', meanval(eup(k, :, :)), &
            & 'abso', meanval(abso(min(nlev - 1, k), :, :)), 'hr', meanval(hr(min(nlev - 1, k), :, :)) * 3600 * 24
        end if
      end do

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

    call write_results_to_file()

    !if (allocated(edir)) &
    !  & call gather_all_toZero(pprts_solver%C_one1, edir, gedir)
    !call gather_all_toZero(pprts_solver%C_one1, edn, gedn)
    !call gather_all_toZero(pprts_solver%C_one1, eup, geup)
    !call gather_all_toZero(pprts_solver%C_one, abso, gabso)

    !  dimnames(1) = 'zlev'
    !  dimnames(2) = 'nx'
    !  dimnames(3) = 'ny'
    !  groups(1) = trim(outfile)
    !  if (lsolar) then
    !    groups(2) = 'edir'; call ncwrite(groups, gedir, ierr, dimnames=dimnames); call CHKERR(ierr)
    !  end if
    !  groups(2) = 'edn'; call ncwrite(groups, gedn, ierr, dimnames=dimnames); call CHKERR(ierr)
    !  groups(2) = 'eup'; call ncwrite(groups, geup, ierr, dimnames=dimnames); call CHKERR(ierr)
    !  dimnames(1) = 'zlay'
    !  groups(2) = 'abso'; call ncwrite(groups, gabso, ierr, dimnames=dimnames); call CHKERR(ierr)

    !  print *, 'dumping z coords'
    !  associate (Ca1 => pprts_solver%C_one_atm1_box)
    !    call getVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
    !    dimnames(1) = 'nlev'
    !    groups(2) = 'zlev'
    !    call ncwrite(groups, z(0, Ca1%zs:Ca1%ze, Ca1%xs, Ca1%ys), ierr, dimnames=dimnames(1:1))
    !    call CHKERR(ierr)
    !    dimnames(1) = 'nlay'
    !    groups(2) = 'zlay'
    !    call ncwrite(groups, &
    !                 & (z(0, Ca1%zs:Ca1%ze - 1, Ca1%xs, Ca1%ys) &
    !                 & + z(0, Ca1%zs + 1:Ca1%ze, Ca1%xs, Ca1%ys) &
    !                 & )*.5_ireals, &
    !                 & ierr, dimnames=dimnames(1:1))
    !    call CHKERR(ierr)
    !    call restoreVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
    !  end associate
    !end if

    ! Tidy up
    call specint_pprts_destroy(specint, pprts_solver, lfinalizepetsc=.true., ierr=ierr)
    call destroy_tenstr_atm(atm)
  contains
    subroutine write_results_to_file()
      if (len_trim(outfile) .gt. 0) then
        if (myid .eq. 0_mpiint) print *, 'Dumping results to out file: ', trim(outfile)
      else
        return
      end if

      associate (C => pprts_solver%C_diff)
        if (allocated(edir)) then
          call ncwrite(&
            & comm, &
            & groups=[character(len=default_str_len) :: outfile, 'edir'], &
            & arr=edir, &
            & ierr=ierr, &
            & arr_shape=[integer :: C%glob_zm, C%glob_xm, C%glob_ym], &
            & dimnames=[character(len=default_str_len) :: 'z', 'x', 'y'], &
            & startp=[integer :: C%zs, C%xs, C%ys] + 1, &
            & countp=shape(edir), &
            & verbose=.True., &
          & deflate_lvl=0)
          call CHKERR(ierr)
          call set_attribute(outfile, 'edir', 'units', 'W/m2', ierr); call CHKERR(ierr)
        end if
        call ncwrite(&
          & comm, &
          & groups=[character(len=default_str_len) :: outfile, 'edn'], &
          & arr=edn, &
          & ierr=ierr, &
          & arr_shape=[integer :: C%glob_zm, C%glob_xm, C%glob_ym], &
          & dimnames=[character(len=default_str_len) :: 'z', 'x', 'y'], &
          & startp=[integer :: C%zs, C%xs, C%ys] + 1, &
          & countp=shape(edn), &
          & verbose=.True., &
          & deflate_lvl=0)
        call CHKERR(ierr)
        call set_attribute(outfile, 'edn', 'units', 'W/m2', ierr); call CHKERR(ierr)
        call ncwrite(&
          & comm, &
          & groups=[character(len=default_str_len) :: outfile, 'eup'], &
          & arr=eup, &
          & ierr=ierr, &
          & arr_shape=[integer :: C%glob_zm, C%glob_xm, C%glob_ym], &
          & dimnames=[character(len=default_str_len) :: 'z', 'x', 'y'], &
          & startp=[integer :: C%zs, C%xs, C%ys] + 1, &
          & countp=shape(eup), &
          & deflate_lvl=0)
        call CHKERR(ierr)
        call set_attribute(outfile, 'eup', 'units', 'W/m2', ierr); call CHKERR(ierr)
      end associate
      associate (C => pprts_solver%C_one)
        call ncwrite(&
          & comm, &
          & groups=[character(len=default_str_len) :: outfile, 'abso'], &
          & arr=abso, &
          & ierr=ierr, &
          & arr_shape=[integer :: C%glob_zm, C%glob_xm, C%glob_ym], &
          & dimnames=[character(len=default_str_len) :: 'zlay', 'x', 'y'], &
          & startp=[integer :: C%zs, C%xs, C%ys] + 1, &
          & countp=shape(abso))
        call CHKERR(ierr)
        call set_attribute(outfile, 'abso', 'units', 'W/m3', ierr); call CHKERR(ierr)
        call ncwrite(&
          & comm, &
          & groups=[character(len=default_str_len) :: outfile, 'hr'], &
          & arr=hr*3600*24, &
          & ierr=ierr, &
          & arr_shape=[integer :: C%glob_zm, C%glob_xm, C%glob_ym], &
          & dimnames=[character(len=default_str_len) :: 'zlay', 'x', 'y'], &
          & startp=[integer :: C%zs, C%xs, C%ys] + 1, &
          & countp=shape(hr))
        call CHKERR(ierr)
        call set_attribute(outfile, 'hr', 'units', 'K/d', ierr); call CHKERR(ierr)
      end associate

    end subroutine
  end subroutine

end module
