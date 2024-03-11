module m_example_pprts_specint_lw_sw_from_dump

#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, default_str_len, pi

  use m_helper_functions, only: &
    & CHKERR, &
    & get_petsc_opt, &
    & linspace, &
    & meanval, &
    & rotate_angle_x, &
    & rotate_angle_z, &
    & spherical_2_cartesian, &
    & toStr

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: gather_all_toZero

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy, load_input_dump

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, abso2hr, print_tenstr_atm

  use m_petsc_helpers, only: getvecpointer, restorevecpointer
  use m_netcdfio, only: ncwrite, set_attribute

  use m_buildings, only: t_pprts_buildings

  implicit none

contains
  subroutine ex_pprts_specint_lw_sw_from_dump(specint, comm, inpfile, outfile)
    character(len=*), intent(in) :: specint           ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: inpfile
    character(len=*), intent(in) :: outfile

    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    real(ireals), allocatable, dimension(:, :, :), target :: hr ! [nlev_merged(-1), nxp, nyp]

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
    real(ireals), allocatable, dimension(:, :, :) :: opt_tau_solar, opt_w0_solar, opt_g_solar, opt_tau_thermal

    type(t_pprts_buildings), allocatable :: opt_buildings_solar
    type(t_pprts_buildings), allocatable :: opt_buildings_thermal

    real(ireals) :: dt
    real(ireals), pointer, dimension(:, :) :: phr

    integer(iintegers) :: icollapse, k, nlev, Niter, kiter
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
      & opt_buildings_solar, opt_buildings_thermal, &
      & opt_tau_solar, opt_w0_solar, opt_g_solar, opt_tau_thermal, &
      & ierr)
    call CHKERR(ierr)

    if (myid .eq. 0) call print_tenstr_atm(atm)

!    sundir = spherical_2_cartesian(phi0, theta0)
!
    call allocate_pprts_solver_from_commandline(pprts_solver, '3_10', ierr); call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, '-solar', lsolar, lflg, ierr); call CHKERR(ierr)
    call get_petsc_opt(PETSC_NULL_CHARACTER, '-thermal', lthermal, lflg, ierr); call CHKERR(ierr)

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

    Niter = 1
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-iter", Niter, lflg, ierr)

    dt = 60 ! sec
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-dt", dt, lflg, ierr)

    do kiter = 1, Niter
      call specint_pprts( &
        specint, comm, pprts_solver, atm, &
        Nx_local, Ny_local, &
        dx, dy, sundir, &
        albedo_thermal, albedo_solar, &
        lthermal, lsolar, &
        edir, edn, eup, abso, &
        nxproc=nxproc, nyproc=nyproc, &
        icollapse=icollapse, &
        opt_time=opt_time, &
        solar_albedo_2d=solar_albedo_2d, &
        thermal_albedo_2d=thermal_albedo_2d, &
        opt_solar_constant=opt_solar_constant, &
        opt_buildings_solar=opt_buildings_solar, &
        opt_buildings_thermal=opt_buildings_thermal, &
        opt_tau_solar=opt_tau_solar, &
        opt_w0_solar=opt_w0_solar, &
        opt_g_solar=opt_g_solar, &
        opt_tau_thermal=opt_tau_thermal)

      if (.not. allocated(hr)) allocate (hr(size(abso, 1), size(abso, 2), size(abso, 3)))

      call abso2hr(atm, abso, hr, ierr); call CHKERR(ierr)
      phr(1:size(hr, 1), 1:size(hr, 2) * size(hr, 3)) => hr
      atm%tlay(:, :) = atm%tlay(:, :) + phr * dt
      if (myid .eq. 0) call print_tenstr_atm(atm)
      if (allocated(sundir)) then
        sundir = rotate_angle_x(sundir, dt * 90._ireals / (12._ireals * 3600))
        sundir = rotate_angle_z(sundir, dt * 180._ireals / (12._ireals * 3600))
      end if
      opt_time = opt_time + dt
      if (myid .eq. 0) print *, 'iter ('//toStr(kiter)//'/'//toStr(Niter)//') ', opt_time, 'sundir after', sundir

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

    end do

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
      integer(mpiint) :: global_shape(3), startp(3)
      if (len_trim(outfile) .gt. 0) then
        if (myid .eq. 0_mpiint) print *, 'Dumping results to out file: ', trim(outfile)
      else
        return
      end if

      associate (C => pprts_solver%C_diff)
        global_shape = [integer(mpiint) :: int(C%glob_zm, mpiint), int(C%glob_xm, mpiint), int(C%glob_ym, mpiint)]
        startp = [integer(mpiint) :: int(C%zs, mpiint), int(C%xs, mpiint), int(C%ys, mpiint)] + 1_mpiint

        if (allocated(edir)) then
          call ncwrite(&
            & comm=comm, &
            & groups=[character(len=default_str_len) :: outfile, 'edir'], &
            & arr=edir, &
            & ierr=ierr, &
            & arr_shape=global_shape, &
            & dimnames=[character(len=default_str_len) :: 'z', 'x', 'y'], &
            & startp=startp, &
            & countp=shape(edir), &
          & deflate_lvl=0)
          call CHKERR(ierr)
          if (myid .eq. 0_mpiint) then
            call set_attribute(outfile, 'edir', 'units', 'W/m2', ierr); call CHKERR(ierr)
          end if
        end if
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: outfile, 'edn'], &
          & arr=edn, &
          & ierr=ierr, &
          & arr_shape=global_shape, &
          & dimnames=[character(len=default_str_len) :: 'z', 'x', 'y'], &
          & startp=startp, &
          & countp=shape(edn), &
          & deflate_lvl=0)
        call CHKERR(ierr)
        if (myid .eq. 0_mpiint) then
          call set_attribute(outfile, 'edn', 'units', 'W/m2', ierr); call CHKERR(ierr)
        end if
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: outfile, 'eup'], &
          & arr=eup, &
          & ierr=ierr, &
          & arr_shape=global_shape, &
          & dimnames=[character(len=default_str_len) :: 'z', 'x', 'y'], &
          & startp=startp, &
          & countp=shape(eup), &
          & deflate_lvl=0)
        call CHKERR(ierr)
        if (myid .eq. 0_mpiint) then
          call set_attribute(outfile, 'eup', 'units', 'W/m2', ierr); call CHKERR(ierr)
        end if
      end associate

      associate (C => pprts_solver%C_one)
        global_shape = [integer(mpiint) :: int(C%glob_zm, mpiint), int(C%glob_xm, mpiint), int(C%glob_ym, mpiint)]
        startp = [integer(mpiint) :: int(C%zs, mpiint), int(C%xs, mpiint), int(C%ys, mpiint)] + 1_mpiint

        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: outfile, 'abso'], &
          & arr=abso, &
          & ierr=ierr, &
          & arr_shape=global_shape, &
          & dimnames=[character(len=default_str_len) :: 'zlay', 'x', 'y'], &
          & startp=startp, &
          & countp=shape(abso))
        call CHKERR(ierr)
        if (myid .eq. 0_mpiint) then
          call set_attribute(outfile, 'abso', 'units', 'W/m3', ierr); call CHKERR(ierr)
        end if
        call ncwrite(&
          & comm=comm, &
          & groups=[character(len=default_str_len) :: outfile, 'hr'], &
          & arr=hr*3600*24, &
          & ierr=ierr, &
          & arr_shape=global_shape, &
          & dimnames=[character(len=default_str_len) :: 'zlay', 'x', 'y'], &
          & startp=startp, &
          & countp=shape(hr))
        call CHKERR(ierr)
        if (myid .eq. 0_mpiint) then
          call set_attribute(outfile, 'hr', 'units', 'K/d', ierr); call CHKERR(ierr)
        end if
      end associate

    end subroutine
  end subroutine

end module
