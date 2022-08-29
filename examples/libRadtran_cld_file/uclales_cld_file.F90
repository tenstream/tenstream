module m_example_uclales_cld_file
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: gather_all_toZero

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only: init_mpi_data_parameters, iintegers, ireals, mpiint, &
                               i0, i1, i2, zero, one, default_str_len, &
                               R_DRY_AIR, CP_DRY_AIR

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy

  use m_dyn_atm_to_rrtmg, only: &
    & destroy_tenstr_atm, &
    & print_tenstr_atm, &
    & reff_from_lwc_and_N, &
    & setup_tenstr_atm, &
    & t_tenstr_atm

  use m_tenstream_options, only: read_commandline_options

  use m_helper_functions, only: &
    & CHKERR, &
    & domain_decompose_2d_petsc, &
    & deg2rad, &
    & get_petsc_opt, &
    & imp_bcast, &
    & resize_arr, &
    & reverse, &
    & spherical_2_cartesian, &
    & toStr

  use m_netcdfio, only: ncload, ncwrite, get_global_attribute, list_global_attributes, get_dim_info

  use m_petsc_helpers, only: getvecpointer, restorevecpointer

  use m_icon_plex_utils, only: create_2d_regular_plex, dmplex_2D_to_3D, &
                               rank0_f90vec_to_plex, dmplex_gvec_from_f90_array, plex_gvec_tozero

  use m_plex_grid, only: t_plexgrid, &
                         setup_plexgrid, create_plex_section

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

  use m_plex_rt, only: init_plex_rt_solver

  use m_plexrt_rrtmg, only: plexrt_rrtmg, destroy_plexrt_rrtmg

  use m_tenstream_interpolation, only: interp_1d
  use m_search, only: find_real_location
  implicit none

contains

  ! (liquid water) potential temperature to temperature
  elemental function Tpot2T(Tpot, p, p0, l) result(T)
    real(ireals), intent(in) :: Tpot, p, p0
    real(ireals), intent(in), optional :: l
    real(ireals) :: T
    real(ireals), parameter :: k = R_DRY_AIR / CP_DRY_AIR, Lv = 2.5e6
    if (present(l)) then
      ! liquid water potential temperature to potential temperature
      T = Tpot + Lv / CP_DRY_AIR * l
      T = T * (p / p0)**k
    else
      T = Tpot * (p / p0)**k
    end if
  end function

  subroutine load_meta_data(comm, cldfile, time, zlev, Nx, Ny, dx, dy, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile
    real(ireals), allocatable, intent(out) :: time(:), zlev(:)
    integer(iintegers) :: Nx, Ny
    real(ireals), intent(out) :: dx, dy
    integer(mpiint), intent(out) :: ierr

    character(len=default_str_len), parameter :: dimnametime = 'time', dimnamez = 'zm', dimnamex = 'xm', dimnamey = 'ym'

    logical :: lflg
    integer(mpiint) :: myid
    real(ireals), allocatable :: dimx(:), dimy(:)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      call list_global_attributes(cldfile, ierr); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, dimnametime], time, ierr, lock_waittime=10_mpiint); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, dimnamez], zlev, ierr, lock_waittime=10_mpiint); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, dimnamex], dimx, ierr, lock_waittime=10_mpiint); call CHKERR(ierr)
      call ncload([character(len=default_str_len) :: cldfile, dimnamey], dimy, ierr, lock_waittime=10_mpiint); call CHKERR(ierr)
      Nx = size(dimx)
      Ny = size(dimy)
      dx = dimx(2) - dimx(1)
      dy = dimy(2) - dimy(1)
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
      call get_petsc_opt(PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr); call CHKERR(ierr)
      print *, 'Timesteps: ', size(time), 'Nx', size(dimx), 'Ny', size(dimy), 'Nz', size(zlev), 'dx', dx, 'dy', dy
    end if
    call imp_bcast(comm, time, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, zlev, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Nx, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, Ny, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dx, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, dy, 0_mpiint, ierr); call CHKERR(ierr)

    ierr = 0
  end subroutine

  subroutine load_timestep_data(comm, cldfile, it, plev, tlev, qv, ql, reliq, ierr)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile
    integer(iintegers), intent(in) :: it
    real(ireals), dimension(:, :, :, :), allocatable, intent(inout) :: plev, tlev, qv, ql, reliq
    integer(mpiint), intent(out) :: ierr

    integer(mpiint) :: myid
    integer :: ostart(4)

    ostart = [integer(iintegers) :: 1, 1, 1, it]

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      call ncload([character(len=default_str_len) :: cldfile, 'p'], plev, ierr, &
        & ostart=ostart, ocount=shape(plev)); call CHKERR(ierr)
      plev = plev * 1e-2 ! Pa to hPa

      call ncload([character(len=default_str_len) :: cldfile, 'q'], qv, ierr, &
        & ostart=ostart, ocount=shape(qv)); call CHKERR(ierr)

      call ncload([character(len=default_str_len) :: cldfile, 'l'], ql, ierr, &
        & ostart=ostart, ocount=shape(ql)); call CHKERR(ierr)
      ql = ql * 1e3

      call ncload([character(len=default_str_len) :: cldfile, 't'], tlev, ierr, &
        & ostart=ostart, ocount=shape(tlev)); call CHKERR(ierr)
      tlev = Tpot2T(tlev, plev, 1013._ireals) ! ignore liquid water potential temperature part for now - layer and lvl handling here is a bit off TODO:

      reliq = reff_from_lwc_and_N(ql, 500._ireals) ! TODO: this is not exact, should give it g/m3 but we have g/kg
      reliq = min(max(reliq, 2.5_ireals), 60._ireals) ! bounds for rrtmg

      print *, 'Min/Max plev ', minval(plev), maxval(plev)
      print *, 'Min/Max tlev ', minval(tlev), maxval(tlev)
      print *, 'Min/Max qv   ', minval(qv), maxval(qv)
      print *, 'Min/Max ql   ', minval(ql), maxval(ql)
      print *, 'Min/Max reliq', minval(reliq), maxval(reliq)
    end if
    call imp_bcast(comm, plev, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, tlev, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, qv, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, ql, 0_mpiint, ierr); call CHKERR(ierr)
    call imp_bcast(comm, reliq, 0_mpiint, ierr); call CHKERR(ierr)

    ierr = 0
  end subroutine

  subroutine run_lw_sw(&
      & specint, &
      & comm, &
      & pprts_solver, &
      & atm, &
      & atm_filename, &
      & nxproc, nyproc, &
      & dx, dy, &
      & phi0, theta0, &
      & albedo_th, albedo_sol, &
      & lsolar, lthermal, &
      & time, &
      & plev, tlev, qv, ql, reliq, &
      & edir, edn, eup, abso)
    character(len=*), intent(in) :: specint                  ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    class(t_solver) :: pprts_solver
    type(t_tenstr_atm), intent(inout) :: atm
    character(len=*), intent(in) :: atm_filename
    integer(iintegers), intent(in) :: nxproc(:), nyproc(:)
    real(ireals), intent(in) :: dx, dy       ! horizontal grid spacing in [m]
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals), intent(in) :: albedo_th, albedo_sol ! broadband ground albedo for solar and thermal spectrum
    logical, intent(in) :: lsolar, lthermal ! switches if solar or thermal computations should be done
    real(ireals), intent(in) :: time

    real(ireals), intent(in), dimension(:, :, :), contiguous, target :: plev, tlev, qv, ql, reliq ! dim(Nz(+1),Nx,Ny)

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:, :, :), intent(out) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    !------------ Local vars ------------------
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: k, nlev

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev, pqv, pql, preliq
    real(ireals) :: sundir(3)

    logical, parameter :: ldebug = .true.

    real(ireals) :: timeofday, tod_offset, tod_theta, tod_phi
    logical :: lflg

    call mpi_comm_rank(comm, myid, ierr)

    if (myid .eq. 0 .and. ldebug) print *, 'Setup Atmosphere...'

    pplev(1:size(plev, 1), 1:size(plev, 2) * size(plev, 3)) => plev
    ptlev(1:size(tlev, 1), 1:size(tlev, 2) * size(tlev, 3)) => tlev
    pqv(1:size(qv, 1), 1:size(qv, 2) * size(qv, 3)) => qv
    pql(1:size(ql, 1), 1:size(ql, 2) * size(ql, 3)) => ql
    preliq(1:size(reliq, 1), 1:size(reliq, 2) * size(reliq, 3)) => reliq

    sundir = spherical_2_cartesian(phi0, theta0)

    tod_offset = 0
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-tod_offset", tod_offset, lflg, ierr); call CHKERR(ierr)

    if (lflg) then
      timeofday = modulo(time / 86400._ireals + tod_offset, 1._ireals)
      tod_phi = modulo(timeofday * 360._ireals, 360._ireals)
      tod_theta = 1._ireals - max(0._ireals, sin(deg2rad(timeofday * 360._ireals - 90))) ! range [0,1]
      sundir = spherical_2_cartesian(tod_phi, tod_theta * 90)
      if (myid .eq. 0) then
        print *, timeofday, 'new phi0', tod_phi
        print *, timeofday, 'new theta0', tod_theta * 90
        print *, timeofday, 'new sundir', sundir
      endif
    end if

    call setup_tenstr_atm(comm, .false., atm_filename, &
                          pplev, ptlev, atm, &
                          d_h2ovmr=pqv, &
                          d_lwc=pql, d_reliq=preliq)

    !if(myid.eq.0) call print_tenstr_atm(atm)

    call specint_pprts(specint, comm, pprts_solver, atm, &
                       size(plev, 2, kind=iintegers), size(plev, 3, kind=iintegers), &
                       dx, dy, sundir, &
                       albedo_th, albedo_sol, &
                       lthermal, lsolar, &
                       edir, edn, eup, abso, &
                       nxproc=nxproc, nyproc=nyproc, opt_time=time)

    nlev = ubound(edn, 1)
    if (myid .eq. 0) then
      if (ldebug) then
        do k = 1, nlev
          if (allocated(edir)) then
            print *, k, 'edir', edir(k, 1, 1), 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
          else
            print *, k, 'edn', edn(k, 1, 1), 'eup', eup(k, 1, 1), abso(min(nlev - 1, k), 1, 1)
          end if
        end do
      end if

      if (allocated(edir)) &
        print *, 'surface :: direct flux', edir(nlev, 1, 1)
      print *, 'surface :: downw flux ', edn(nlev, 1, 1)
      print *, 'surface :: upward fl  ', eup(nlev, 1, 1)
      print *, 'surface :: absorption ', abso(nlev - 1, 1, 1)

      if (allocated(edir)) &
        print *, 'TOA :: direct flux', edir(1, 1, 1)
      print *, 'TOA :: downw flux ', edn(1, 1, 1)
      print *, 'TOA :: upward fl  ', eup(1, 1, 1)
      print *, 'TOA :: absorption ', abso(1, 1, 1)
    end if
  end subroutine

  subroutine example_uclales_cld_file_pprts(specint, comm, &
                                            cldfile, atm_filename, outfile, &
                                            albedo_th, albedo_sol, &
                                            lsolar, lthermal, &
                                            phi0, theta0, &
                                            tstart, tend, tinc)

    character(len=*), intent(in) :: specint ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: cldfile, atm_filename, outfile
    real(ireals), intent(in) :: albedo_th, albedo_sol
    logical, intent(in) :: lsolar, lthermal
    real(ireals), intent(in) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    integer(iintegers), intent(in) :: tstart, tend, tinc

    real(ireals), pointer :: z(:, :, :, :) => null(), z1d(:) => null() ! dim Nz+1
    character(len=default_str_len) :: groups(3), dimnames(3)

    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    real(ireals), allocatable, dimension(:, :, :) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf

    class(t_solver), allocatable :: pprts_solver
    type(t_tenstr_atm) :: atm

    real(ireals) :: dx, dy
    real(ireals), allocatable :: time(:), hhl(:)

    integer(mpiint) :: myid, ierr
    integer(iintegers) :: Nx_global, Ny_global, is, ie, js, je
    integer(iintegers) :: nxp, nyp ! local sizes of domain, nzp being number of layers
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    real(ireals), allocatable, dimension(:, :, :, :) :: plev, tlev, qv, ql, reliq ! dim(Nz(+1),Nx,Ny)

    integer(iintegers) :: it

    call mpi_comm_rank(comm, myid, ierr)

    call load_meta_data(comm, cldfile, time, hhl, Nx_global, Ny_global, dx, dy, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      print *, 'Loaded UCLALES File with:'
      print *, 'Nx, Ny:', Nx_global, Ny_global
      print *, 'dx, dy:', dx, dy
      print *, 'hhl', hhl
    end if
    if (tstart .lt. lbound(time, 1)) &
      & call CHKERR(1_mpiint, 'tstart '//toStr(tstart)//' invalid because smaller than '//toStr(lbound(time, 1)))
    if (tstart .gt. ubound(time, 1)) &
      & call CHKERR(1_mpiint, 'tstart '//toStr(tstart)//' invalid because larger  than '//toStr(ubound(time, 1)))
    if (tend .lt. lbound(time, 1)) &
      & call CHKERR(1_mpiint, 'tend '//toStr(tend)//' invalid because smaller than '//toStr(lbound(time, 1)))
    if (tend .gt. ubound(time, 1)) &
      & call CHKERR(1_mpiint, 'tend '//toStr(tend)//' invalid because larger  than '//toStr(ubound(time, 1)))

    ! Determine Domain Decomposition
    call domain_decompose_2d_petsc(comm, &
      & Nx_global=Nx_global, &
      & Ny_global=Ny_global, &
      & Nx_local=nxp, &
      & Ny_local=nyp, &
      & xStart=is, &
      & yStart=js, &
      & nxproc=nxproc, &
      & nyproc=nyproc, &
      & ierr=ierr); call CHKERR(ierr)
    is = is + 1 ! fortran based indices
    js = js + 1 ! fortran based indices
    ie = is + nxp - 1
    je = js + nyp - 1

    call allocate_pprts_solver_from_commandline(pprts_solver, default_solver='3_10', ierr=ierr); call CHKERR(ierr)

    allocate (plev(size(hhl), Nx_global, Ny_global, 1))
    allocate (tlev(size(hhl), Nx_global, Ny_global, 1))
    allocate (qv(size(hhl) - 1, Nx_global, Ny_global, 1))
    allocate (ql(size(hhl) - 1, Nx_global, Ny_global, 1))
    allocate (reliq(size(hhl) - 1, Nx_global, Ny_global, 1))

    do it = tstart, tend, tinc
      if (myid .eq. 0) print *, 'Computing timestep', it, '(time='//toStr(time(it))
      call load_timestep_data(comm, cldfile, it, plev, tlev, qv, ql, reliq, ierr); call CHKERR(ierr)

      call run_lw_sw(&
        & specint, &
        & comm, &
        & pprts_solver, &
        & atm, &
        & atm_filename, &
        & nxproc, nyproc, &
        & dx, dy, phi0, theta0, &
        & albedo_th, albedo_sol, &
        & lsolar, lthermal, &
        & time(it), &
        & plev(:, is:ie, js:je, 1), &
        & tlev(:, is:ie, js:je, 1), &
        & qv(:, is:ie, js:je, 1), &
        & ql(:, is:ie, js:je, 1), &
        & reliq(:, is:ie, js:je, 1), &
        & edir, edn, eup, abso)

      if (len_trim(outfile) .gt. 0) then
        groups(1) = trim(outfile)
        groups(3) = trim(toStr(it))

        associate (&
            & C => pprts_solver%C_one,     &
            & C1 => pprts_solver%C_one1,    &
            & Ca => pprts_solver%C_one_atm, &
            & Ca1 => pprts_solver%C_one_atm1_box)

          dimnames(1) = 'zlev'
          dimnames(2) = 'nx'
          dimnames(3) = 'ny'
          if (allocated(edir)) then
            call gather_all_toZero(C1, edir, gedir)
            if (myid .eq. 0) then
              print *, 'dumping direct radiation with local and global shape', shape(edir), ':', shape(gedir)
              groups(2) = 'edir'; call ncwrite(groups, gedir, ierr, dimnames=dimnames); call CHKERR(ierr)
            end if
          end if
          call gather_all_toZero(C1, edn, gedn)
          call gather_all_toZero(C1, eup, geup)
          call gather_all_toZero(C, abso, gabso)
          if (myid .eq. 0) then
            print *, 'dumping edn radiation with local and global shape', shape(edn), ':', shape(gedn)
            groups(2) = 'edn'; call ncwrite(groups, gedn, ierr, dimnames=dimnames); call CHKERR(ierr)

            print *, 'dumping eup radiation with local and global shape', shape(eup), ':', shape(geup)
            groups(2) = 'eup'; call ncwrite(groups, geup, ierr, dimnames=dimnames); call CHKERR(ierr)

            print *, 'dumping abso radiation with local and global shape', shape(abso), ':', shape(gabso)
            dimnames(1) = 'zlay'
            groups(2) = 'abso'; call ncwrite(groups, gabso, ierr, dimnames=dimnames); call CHKERR(ierr)

            print *, 'dumping z coords'
            call getVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
            dimnames(1) = 'nlev'
            groups(2) = 'zlev'; call ncwrite(groups, z(0, Ca1%zs:Ca1%ze, Ca1%xs, Ca1%ys), ierr, dimnames=dimnames(1:1))
            call CHKERR(ierr)
            dimnames(1) = 'nlay'
            groups(2) = 'zlay'; call ncwrite(groups, &
 & (z(0, Ca1%zs:Ca1%ze - 1, Ca1%xs, Ca1%ys) + z(0, Ca1%zs + 1:Ca1%ze, Ca1%xs, Ca1%ys))*.5_ireals, &
 & ierr, dimnames=dimnames(1:1))
            call CHKERR(ierr)
            call restoreVecPointer(Ca1%da, pprts_solver%atm%hhl, z1d, z)
          end if

        end associate
      end if
    end do

    call specint_pprts_destroy(specint, pprts_solver, lfinalizepetsc=.true., ierr=ierr)
    call destroy_tenstr_atm(atm)
  end subroutine
end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only: mpiint, share_dir
  use m_helper_functions, only: get_petsc_opt
  use m_example_uclales_cld_file

  implicit none

  integer(mpiint) :: ierr, myid
  logical :: lthermal, lsolar, lflg
  character(len=10*default_str_len) :: cldfile, outfile
  real(ireals) :: Ag, phi0, theta0, Tsrfc, dTdz
  character(len=default_str_len) :: atm_filename, specint
  integer(iintegers) :: tstart, tend, tinc

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call read_commandline_options(mpi_comm_world)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  specint = 'no default set'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-specint', specint, lflg, ierr); call CHKERR(ierr)

  call get_petsc_opt(PETSC_NULL_CHARACTER, '-cld', cldfile, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) call CHKERR(1_mpiint, 'need to supply a cloud filename... please call with -cld <libRadtran_cloud_file.nc>')

  outfile = ''
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)

  atm_filename = share_dir//'tenstream_default.atm'
  call get_petsc_opt(PETSC_NULL_CHARACTER, '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)

  Ag = .1
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Ag", Ag, lflg, ierr); call CHKERR(ierr)

  lsolar = .true.
  lthermal = .true.
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-solar", lsolar, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg, ierr); call CHKERR(ierr)

  phi0 = 270
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  theta0 = 60
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr); call CHKERR(ierr)

  Tsrfc = 288
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-Tsrfc", Tsrfc, lflg, ierr); call CHKERR(ierr)
  dTdz = -6.5_ireals * 1e-3_ireals
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-dTdz", dTdz, lflg, ierr); call CHKERR(ierr)

  tstart = 800
  tend = 900
  tinc = 1
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-tstart", tstart, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-tend", tend, lflg, ierr); call CHKERR(ierr)
  call get_petsc_opt(PETSC_NULL_CHARACTER, "-tinc", tinc, lflg, ierr); call CHKERR(ierr)

  call example_uclales_cld_file_pprts(&
    & specint, &
    & mpi_comm_world, cldfile, atm_filename, outfile, &
    & zero, Ag, lsolar, lthermal, phi0, theta0, &
    & tstart, tend, tinc)
  call mpi_finalize(ierr)
end program
