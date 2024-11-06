module m_examples_plex_specint_fish

#include "petsc/finclude/petsc.h"
  use petsc

  use m_helper_functions, only: &
    & CHKERR, CHKWARN, &
    & cstr, &
    & imp_reduce_mean, &
    & meanval, &
    & reverse, &
    & toStr

  use m_data_parameters, only: ireals, iintegers, mpiint, imp_ireals, &
                               default_str_len, &
                               i0, i1, i2, i3, i4, i5, &
                               zero, one, &
                               init_mpi_data_parameters

  use m_plex_grid, only: t_plexgrid, setup_plexgrid

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

  use m_plex_rt, only: init_plex_rt_solver

  use m_specint_plexrt, only: specint_plexrt, specint_plexrt_destroy

  use m_dyn_atm_to_rrtmg, only: &
    & t_tenstr_atm, t_bg_atm, &
    & load_atmfile, &
    & setup_tenstr_atm

  use m_icon_plex_utils, only: &
    & create_2d_fish_plex, &
    & create_2d_regular_plex, &
    & dmplex_2D_to_3D

  use m_tenstream_interpolation, only: interp_1d
  use m_search, only: find_real_location

  implicit none

contains

  subroutine ex_plex_specint_fish(&
      & specint, &
      & comm, &
      & lverbose, &
      & lregular_mesh, &
      & lthermal, lsolar, &
      & atm_filename, &
      & Nx, Ny, Nz, &
      & dx, dz, &
      & Ag_thermal, &
      & Ag_solar, &
      & sundir, &
      & lwc, &
      & edir, edn, eup, abso)

    character(len=*), intent(in) :: specint
    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lverbose, lregular_mesh, lthermal, lsolar
    character(len=*), intent(in) :: atm_filename
    integer(iintegers), intent(in) :: Nx, Ny, Nz
    real(ireals), intent(in) :: dx, dz, Ag_solar, Ag_thermal
    real(ireals), intent(in) :: sundir(3)
    real(ireals), intent(in) :: lwc
    real(ireals), allocatable, dimension(:, :), intent(out) :: edir, edn, eup, abso

    type(tDM) :: dm2d, dm2d_dist, dm3d
    real(ireals) :: hhl(Nz), z_location

    integer(mpiint) :: myid, numnodes, ierr

    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)
    class(t_plex_solver), allocatable :: solver

    type(t_bg_atm), allocatable :: bg_atm
    type(t_tenstr_atm) :: atm
    integer(iintegers) :: k, fStart, fEnd, Ncol, dNlev, dNlay, Nlev
    real(ireals), allocatable :: col_plev(:, :), col_tlev(:, :)
    real(ireals), allocatable :: col_lwc(:, :), col_reff(:, :)
    real(ireals), allocatable :: col_tskin(:)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

    if (lregular_mesh) then
      if (lverbose .and. myid .eq. 0) print *, 'Initialize rectilinear mesh...'
      call create_2d_regular_plex(comm, Nx, Ny, dm2d, dm2d_dist, opt_dx=dx)
    else
      if (lverbose .and. myid .eq. 0) print *, 'Initialize fish mesh...'
      call create_2d_fish_plex(comm, Nx, Ny, dm2d, dm2d_dist, opt_dx=dx)
    end if

    hhl(1) = 0
    do k = 2, Nz
      hhl(k) = hhl(k - 1) + dz
    end do

    call DMPlexGetHeightStratum(dm2d_dist, i0, fStart, fEnd, ierr); call CHKERR(ierr)
    Ncol = fEnd - fStart
    dNlev = Nz; dNlay = dNlev - 1

    if (lverbose .and. myid .eq. 0) print *, 'Dynamics Grid has Size Nlev, Ncol:', dNlev, Ncol
    if (Ncol .eq. 0) call CHKERR(1_mpiint, 'We have a process that has nothing to do? '// &
                                 'Maybe decrease the number of processes or increase the problem size')

    ! prepare atmosphere
    allocate (col_tlev(dNlev, Ncol))
    allocate (col_plev(dNlev, Ncol))

    ! Load the background atmosphere file and interpolate pressure and temperature from that
    call load_atmfile(comm, atm_filename, bg_atm)
    do k = 1, Nz
      z_location = find_real_location(bg_atm%zt, hhl(k))
      col_plev(k, :) = interp_1d(z_location, bg_atm%plev)
      col_tlev(k, :) = interp_1d(z_location, bg_atm%tlev)
    end do
    deallocate (bg_atm)

    if (lverbose .and. myid .eq. 0) then
      print *, 'Dynamics Grid Pressure and Temperature'
      do k = 1, dNlev
        print *, k, col_plev(k, i1), col_tlev(k, i1)
      end do
    end if

    do k = 2, Ncol
      col_plev(:, k) = col_plev(:, i1)
      col_tlev(:, k) = col_tlev(:, i1)
    end do

    allocate (col_lwc(dNlay, Ncol), col_reff(dNlay, Ncol))
    col_lwc = 0
    col_reff = 10

    allocate (col_tskin(Ncol))
    col_tskin = 300

    if (lwc .gt. 0) then
      if (lverbose .and. myid .eq. 0) print *, "Adding cloud at top of dynamics grid k = "//toStr(dNlay)
      col_lwc(dNlay, :) = lwc
    end if

    call setup_tenstr_atm(&
      & comm, .false., atm_filename, &
      & col_plev, col_tlev, atm, &
      & d_skin_temperature=col_tskin, &
      & d_lwc=col_lwc, d_reliq=col_reff)

    Nlev = size(atm%plev, 1, kind=iintegers)
    call dmplex_2D_to_3D(dm2d_dist, Nlev, reverse(atm%zt(:, i1)), [zero, zero, -huge(zero) * 1e-1_ireals], dm3d, zindex)

    call setup_plexgrid(dm2d_dist, dm3d, Nlev - 1, zindex, plex, hhl=reverse(atm%zt(:, i1)))

    call DMDestroy(dm2d, ierr); call CHKERR(ierr)
    call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)
    deallocate (zindex)

    if (lregular_mesh) then
      call allocate_plexrt_solver_from_commandline(solver, 'rectilinear_5_8')
    else
      call allocate_plexrt_solver_from_commandline(solver, '5_8')
    end if
    call init_plex_rt_solver(plex, solver)

    call specint_plexrt(&
      & specint, &
      & solver, &
      & atm, &
      & sundir, &
      & albedo_thermal=Ag_thermal, &
      & albedo_solar=Ag_solar, &
      & lthermal=lthermal, &
      & lsolar=lsolar, &
      & edir=edir, edn=edn, eup=eup, abso=abso)

    if (lverbose) call print_overview()

    if (meanval(edn) .lt. epsilon(edn)) &
      call CHKERR(1_mpiint, 'Mean Edn is really small,'// &
                  'the solver probably had a problem but did not fail.'// &
                  'Caution! Your results may be garbage!')

    call specint_plexrt_destroy(specint, solver, lfinalizepetsc=.false., ierr=ierr); call CHKERR(ierr)

  contains
    subroutine print_overview()
      real(ireals) :: medir, medn, meup, mabso
      integer(iintegers) :: k
      if (lsolar) then
        if (myid .eq. 0) then
          print *, ''
          print *, cstr('Solar Radiation', 'green'), lsolar, cstr('Thermal Radiation', 'green'), lthermal
          print *, cstr('Avg.horiz', 'blue')// &
            cstr(' k   Edir          Edn            Eup            abso', 'blue')
        end if
        do k = 1, ubound(abso, 1)
          call imp_reduce_mean(comm, edir(k, :), medir)
          call imp_reduce_mean(comm, edn(k, :), medn)
          call imp_reduce_mean(comm, eup(k, :), meup)
          call imp_reduce_mean(comm, abso(k, :), mabso)
          if (myid .eq. 0) print *, k, medir, medn, meup, mabso
        end do
        call imp_reduce_mean(comm, edir(ubound(edn, 1), :), medir)
        call imp_reduce_mean(comm, edn(ubound(edn, 1), :), medn)
        call imp_reduce_mean(comm, eup(ubound(edn, 1), :), meup)
        if (myid .eq. 0) print *, k, medir, medn, meup
      else
        if (myid .eq. 0) then
          print *, ''
          print *, cstr('Solar Radiation', 'green'), lsolar, cstr('Thermal Radiation', 'green'), lthermal
          print *, cstr('Avg.horiz', 'blue')// &
            cstr(' k  Edn             Eup           abso', 'blue')
        end if
        do k = 1, ubound(abso, 1)
          call imp_reduce_mean(comm, edn(k, :), medn)
          call imp_reduce_mean(comm, eup(k, :), meup)
          call imp_reduce_mean(comm, abso(k, :), mabso)
          if (myid .eq. 0) print *, k, medn, meup, mabso
        end do
        call imp_reduce_mean(comm, edn(ubound(edn, 1), :), medn)
        call imp_reduce_mean(comm, eup(ubound(edn, 1), :), meup)
        if (myid .eq. 0) print *, k, medn, meup
      end if
    end subroutine

  end subroutine

end module
