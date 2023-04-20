module m_examples_pprts_specint_buildings_from_file

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & iintegers, ireals, mpiint, &
    & default_str_len

  use m_helper_functions, only: &
    & CHKERR, &
    & cstr, &
    & domain_decompose_2d_petsc, &
    & get_arg, &
    & get_petsc_opt, &
    & is_inrange, &
    & imp_bcast, &
    & linspace, &
    & spherical_2_cartesian, &
    & toStr

  use m_netcdfio, only: ncload

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline

  use m_pprts, only: gather_all_to_all

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy

  ! tenstr_atm holds info about tracer and merges dynamics grid vars with background grids
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

  use m_buildings, only: &
    & t_pprts_buildings, &
    & init_buildings, &
    & clone_buildings, &
    & check_buildings_consistency, &
    & faceidx_by_cell_plus_offset

  use m_boxmc_geometry, only: &
    & PPRTS_TOP_FACE

  implicit none

contains
  subroutine ex_pprts_specint_buildings_from_file(&
      & specint,                            &
      & comm, lverbose,                     &
      & lthermal, lsolar,                   &
      & buildings_filename,                 &
      & dx, dy,                             &
      & atm_filename,                       &
      & phi0, theta0,                       &
      & Ag_solar, Ag_thermal,               &
      & skin_temp,                          &
      & gedir, gedn, geup, gabso,           &
      & buildings_solar, buildings_thermal, &
      & local_dims, icollapse)

    character(len=*), intent(in) :: specint           ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lverbose, lthermal, lsolar
    character(len=*), intent(in) :: buildings_filename
    real(ireals), intent(in) :: dx, dy               ! grid spacing in [m]
    character(len=*), intent(in) :: atm_filename     ! e.g. 'afglus_100m.dat'
    real(ireals), intent(in) :: phi0, theta0         ! sun azimuth(phi) and zenith(theta) angle
    real(ireals), intent(in) :: Ag_solar, Ag_thermal ! surface albedo
    real(ireals), intent(in) :: skin_temp            ! surface skin temperature
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: gedir, gedn, geup, gabso
    type(t_pprts_buildings), allocatable, intent(inout), optional :: buildings_solar, buildings_thermal
    integer(iintegers), intent(out), optional :: local_dims(:) ! local domain indices (zs, zm, xs, xm, ys, ym), dim(6)
    integer(iintegers), intent(in), optional :: icollapse

    integer(iintegers) :: Nx, Ny, Nlay ! global domain size from netcdf buildings file
    integer(iintegers) :: nxp, nyp, xs, ys ! local domain size in x and y aswell as start indices
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    real(ireals) :: sundir(3) ! cartesian vector pointing from sun to earth

    real(ireals), allocatable, dimension(:, :, :), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), allocatable, dimension(:, :, :), target :: tlev ! Temperature on layer interfaces [K]
    real(ireals), allocatable, dimension(:, :), target :: tskin ! Temperature on ground surface [K]
    real(ireals), allocatable, dimension(:, :, :, :) :: glob_buildings_temperature, glob_buildings_albedo ! dim [6,Nz,Ny,Nx]

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev
    real(ireals), pointer, dimension(:) :: ptskin

    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    class(t_solver), allocatable :: solver
    type(t_tenstr_atm) :: atm

    integer(iintegers) :: k, iface
    integer(mpiint) :: myid, ierr

    call init_mpi_data_parameters(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      call ncload([character(len=default_str_len) :: buildings_filename, 'albedo'], &
        & glob_buildings_albedo, ierr); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: buildings_filename, 'temperature'], &
        & glob_buildings_temperature, ierr); call CHKERR(ierr)
    end if
    call imp_bcast(comm, glob_buildings_albedo, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, glob_buildings_temperature, 0_mpiint, ierr); call CHKERR(ierr)
    print *, 'shape glob_buildings_albedo', shape(glob_buildings_albedo)
    Nx = size(glob_buildings_albedo, dim=4)
    Ny = size(glob_buildings_albedo, dim=3)
    Nlay = size(glob_buildings_albedo, dim=2)

    call domain_decompose_2d_petsc(comm, Nx, Ny, &
      & nxp, nyp, xs, ys, nxproc, nyproc, ierr); call CHKERR(ierr)

    sundir = spherical_2_cartesian(phi0, theta0)
    if (lverbose) print *, 'sundir:', sundir

    allocate (plev(Nlay + 1, nxp, nyp))
    allocate (tlev(Nlay + 1, nxp, nyp))
    allocate (tskin(nxp, nyp))

    ! Start with a dynamics grid ranging from 1000 hPa up to 900 hPa and a
    ! Temperature difference of 10K
    do k = 1, Nlay + 1
      plev(k, :, :) = linspace(k, [1e3_ireals, 900._ireals], Nlay + 1)
      tlev(k, :, :) = linspace(k, [288._ireals, 278._ireals], Nlay + 1)
    end do
    tskin(:, :) = skin_temp

    pplev(1:size(plev, 1), 1:size(plev, 2) * size(plev, 3)) => plev
    ptlev(1:size(tlev, 1), 1:size(tlev, 2) * size(tlev, 3)) => tlev
    ptskin(1:size(tskin)) => tskin

    call setup_tenstr_atm( &
      & comm, .false.,     &
      & atm_filename,      &
      & pplev, ptlev,      &
      & atm,               &
      & d_skin_temperature=ptskin)

    ! only init grid structures
    call specint_pprts(specint,&
      & comm,                  &
      & solver, atm,           &
      & nxp, nyp, dx, dy,      &
      & sundir,                &
      & Ag_thermal, Ag_solar,  &
      & lthermal, lsolar,      &
      & edir, edn, eup, abso,  &
      & icollapse=get_arg(1_iintegers, icollapse), &
      & nxproc=nxproc,         &
      & nyproc=nyproc,         &
      & lonly_initialize=.true.)

    if (present(local_dims)) then
      local_dims(:) = [ &
        & solver%C_one%zs, solver%C_one%zm, &
        & solver%C_one%xs, solver%C_one%xm, &
        & solver%C_one%ys, solver%C_one%ym]
    end if

    if (present(buildings_solar) .and. present(buildings_thermal)) then
      call setup_buildings(glob_buildings_albedo, glob_buildings_temperature)

      call specint_pprts(specint,&
        & comm,                  &
        & solver, atm,           &
        & nxp, nyp, dx, dy,      &
        & sundir,                &
        & Ag_thermal, Ag_solar,  &
        & lthermal, lsolar,      &
        & edir, edn, eup, abso,  &
        & icollapse=get_arg(1_iintegers, icollapse), &
        & nxproc=nxproc,         &
        & nyproc=nyproc,         &
        & opt_buildings_solar=buildings_solar, &
        & opt_buildings_thermal=buildings_thermal)

      if (lverbose) then
        if (lsolar) then
          print *, 'Flux on Top Face of buildings'//new_line('')// &
            & '   Building cell'// &
            & cstr('     edir', 'red')// &
            & cstr('     incoming', 'green')// &
            & cstr('     outgoing', 'blue')
          associate (B => buildings_solar)
            do k = 0, size(B%iface) / 6 - 1
              iface = k * 6 + PPRTS_TOP_FACE
              print *, k, &
                & ' '//cstr(toStr(B%edir(iface)), 'red')// &
                & ' '//cstr(toStr(B%incoming(iface)), 'green')// &
                & ' '//cstr(toStr(B%outgoing(iface)), 'blue')
            end do
          end associate
        else
          print *, 'Flux on Top Face of buildings'//new_line('')// &
            & '   Building cell'// &
            & cstr('     incoming', 'green')// &
            & cstr('     outgoing', 'blue')
          associate (B => buildings_thermal)
            do k = 0, size(B%iface) / 6 - 1
              iface = k * 6 + PPRTS_TOP_FACE
              print *, k, &
                & ' '//cstr(toStr(B%incoming(iface)), 'green')// &
                & ' '//cstr(toStr(B%outgoing(iface)), 'blue')
            end do
          end associate
        end if
      end if

    else ! without buildings
      call specint_pprts(specint,&
        & comm,                  &
        & solver, atm,           &
        & nxp, nyp, dx, dy,      &
        & sundir,                &
        & Ag_thermal, Ag_solar,  &
        & lthermal, lsolar,      &
        & edir, edn, eup, abso,  &
        & icollapse=get_arg(1_iintegers, icollapse), &
        & nxproc=nxproc,         &
        & nyproc=nyproc)
    end if

    if (allocated(edir)) &
      & call gather_all_to_all(solver%C_one1, edir, gedir)
    call gather_all_to_all(solver%C_one1, edn, gedn)
    call gather_all_to_all(solver%C_one1, eup, geup)
    call gather_all_to_all(solver%C_one, abso, gabso)

    ! Tidy up
    call specint_pprts_destroy(specint, solver, lfinalizepetsc=.true., ierr=ierr); call CHKERR(ierr)
    call destroy_tenstr_atm(atm)

  contains
    subroutine setup_buildings(glob_buildings_albedo, glob_buildings_temperature)
      real(ireals), dimension(:, :, :, :), intent(in) :: glob_buildings_temperature, glob_buildings_albedo
      integer(iintegers) :: i, j, k
      integer(iintegers) :: Nbuildings ! local count cells that we have to set
      integer(iintegers) :: Nfaces     ! local count faces that we have to set

      logical :: lflg
      real(ireals) :: v

      ! First count how many boxes, i.e. faces we have to allocate on this rank
      Nbuildings = 0

      associate (C1 => solver%C_one)
        do k = 1, Nlay
          do j = 1, Ny
            do i = 1, Nx
              if (any(glob_buildings_albedo(:, k, j, i) .ge. 0)) then
                if (have_box(C1%glob_zm - k + 1, i, j)) then
                  print *, myid, 'have building block at i,j,k', i, j, k, &
                    & 'albedo', glob_buildings_albedo(1, k, j, i), 'T', glob_buildings_temperature(1, k, j, i)
                  Nbuildings = Nbuildings + 1
                end if
              end if
            end do
          end do
        end do

        Nfaces = Nbuildings * 6
        if (lverbose) then
          print *, 'rank '//toStr(myid)//' has '//toStr(Nbuildings)//' buildings cells with '//toStr(Nfaces)//' faces'
        end if

        call init_buildings(buildings_solar, &
          & [integer(iintegers) :: 6, C1%zm, C1%xm, C1%ym], &
          & Nfaces, &
          & ierr); call CHKERR(ierr)

        Nbuildings = 0
        do k = 1, Nlay
          do j = 1, Ny
            do i = 1, Nx
              if (any(glob_buildings_albedo(:, k, j, i) .ge. 0)) then
                if (have_box(C1%glob_zm - k + 1, i, j)) then
                  Nbuildings = Nbuildings + 1
                  call fill_cell_with_building(buildings_solar, Nbuildings, C1%glob_zm - k + 1, i, j)
                  buildings_solar%albedo((Nbuildings - 1) * 6 + 1:(Nbuildings - 1) * 6 + 6) = glob_buildings_albedo(:, k, j, i)
                end if
              end if
            end do
          end do
        end do

        call get_petsc_opt(PETSC_NULL_CHARACTER, "-override_buildings_albedo", v, lflg, ierr); call CHKERR(ierr)
        if (lflg) buildings_solar%albedo(:) = v

        call clone_buildings(&
          & buildings_solar, &
          & buildings_thermal, &
          & l_copy_data=.true., &
          & ierr=ierr); call CHKERR(ierr)
        if (.not. allocated(buildings_thermal%temp)) allocate (buildings_thermal%temp(Nfaces))

        Nbuildings = 0
        do k = 1, Nlay
          do j = 1, Ny
            do i = 1, Nx
              if (any(glob_buildings_albedo(:, k, j, i) .ge. 0)) then
                if (have_box(C1%glob_zm - k + 1, i, j)) then
                  Nbuildings = Nbuildings + 1
                  buildings_thermal%temp((Nbuildings - 1) * 6 + 1:(Nbuildings - 1) * 6 + 6) = glob_buildings_temperature(:, k, j, i)
                end if
              end if
            end do
          end do
        end do
        call get_petsc_opt(PETSC_NULL_CHARACTER, "-override_buildings_temperature", v, lflg, ierr); call CHKERR(ierr)
        if (lflg) buildings_thermal%temp(:) = v

        call check_buildings_consistency(buildings_solar, C1%zm, C1%xm, C1%ym, ierr); call CHKERR(ierr)
        call check_buildings_consistency(buildings_thermal, C1%zm, C1%xm, C1%ym, ierr); call CHKERR(ierr)
      end associate
    end subroutine

    subroutine fill_cell_with_building(B, m, gk, gi, gj)
      type(t_pprts_buildings), intent(inout) :: B
      integer(iintegers), intent(in) :: m ! building cell nr
      integer(iintegers), intent(in) :: gk, gi, gj ! global inidces
      integer(iintegers) :: k, i, j ! local subdomain indices
      integer(iintegers) :: iface

      k = gk - solver%C_one%zs
      i = gi - solver%C_one%xs
      j = gj - solver%C_one%ys

      do iface = 1, 6
        B%iface((m - 1) * 6 + iface) = faceidx_by_cell_plus_offset( &
          & B%da_offsets, &
          & k, &
          & i, &
          & j, iface)
      end do
    end subroutine

    logical function have_box(gk, gi, gj)
      integer(iintegers), intent(in) :: gk, gi, gj ! global inidces
      have_box = all([ &
        & is_inrange(gk, solver%C_one%zs + 1, solver%C_one%ze + 1), &
        & is_inrange(gi, solver%C_one%xs + 1, solver%C_one%xe + 1), &
        & is_inrange(gj, solver%C_one%ys + 1, solver%C_one%ye + 1)])
    end function
  end subroutine
end module
