!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

module m_compute_training_data
  use m_data_parameters, only: &
    & AVOGADRO, &
    & EARTHACCEL, &
    & MOLMASSAIR, &
    & K_BOLTZMANN, &
    & R_DRY_AIR, &
    & default_str_len, &
    & iintegers, &
    & ireals, &
    & mpiint

  use m_helper_functions, only: &
    & CHKERR, &
    & CHKWARN, &
    & get_arg, &
    & imp_min_mean_max, &
    & resize_arr, &
    & reverse, &
    & spherical_2_cartesian, &
    & toStr

  use m_search, only: find_real_location
  use m_netcdfIO, only: ncload, ncwrite
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, print_tenstr_atm, planck
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, pprts_get_result, set_angles

  use m_repwvl_base, only: t_repwvl_data, repwvl_init, repwvl_dtau
  use m_repwvl_pprts, only: repwvl_optprop
  use m_mie_tables, only: mie_tables_init, mie_water_table, mie_optprop
  use m_fu_ice, only: fu_ice_init
  use m_rayleigh, only: rayleigh

  implicit none

  logical, parameter :: lboundary_flx_only=.True.

contains

  subroutine load_atmospheres(comm, bg_atm_filename, atm_filename, Npert_atmo, atm, ierr, lverbose)
    ! Garand profiles (9 columns)
    ! p, T, z, H2O, O3, CO2, N2O, CO, CH4
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: bg_atm_filename, atm_filename
    integer(iintegers), intent(in) :: Npert_atmo
    type(t_tenstr_atm), intent(out) :: atm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose

    integer(iintegers) :: ke, ke1, Natm
    real(ireals), allocatable :: atmdat(:, :, :) ! dim (cols, nlev, Natm)
    character(len=default_str_len) :: groups(2)

    groups(1) = trim(atm_filename)
    groups(2) = 'Atmosphere'
    if (get_arg(.false., lverbose)) then
      print *, 'Reading Atmospheres from '//trim(groups(1))//' : '//trim(groups(2))
    end if
    call ncload(groups, atmdat, ierr, lverbose); call CHKERR(ierr)
    if (get_arg(.false., lverbose)) then
      print *, 'Shape Atmospheres: (tracer, nlev, Natm)', shape(atmdat)
    end if
    ke1 = size(atmdat, dim=2)
    ke = ke1 - 1_iintegers

    Natm = size(atmdat, dim=3)

    call resize_arr(Natm * Npert_atmo, atmdat, dim=3_mpiint, lrepeat=.true.)

    select case (size(atmdat, dim=1))
    case (9) ! Garand
      call setup_tenstr_atm(                      &
        & comm,                                   &
        & lTOA_to_srfc=.false.,                 &
        & atm_filename=bg_atm_filename,         &
        & d_plev=atmdat(1, :, :) * 1e-2_ireals, &
        & d_tlev=atmdat(2, :, :),               &
        & atm=atm,                         &
        & d_h2ovmr=(atmdat(4, 1:ke, :) + atmdat(4, 2:ke1, :))*.5_ireals, &
        & d_o3vmr=(atmdat(5, 1:ke, :) + atmdat(5, 2:ke1, :))*.5_ireals, &
        & d_co2vmr=(atmdat(6, 1:ke, :) + atmdat(6, 2:ke1, :))*.5_ireals, &
        & d_n2ovmr=(atmdat(7, 1:ke, :) + atmdat(7, 2:ke1, :))*.5_ireals, &
        & d_ch4vmr=(atmdat(9, 1:ke, :) + atmdat(9, 2:ke1, :))*.5_ireals)
    case default
      call CHKERR(1_mpiint, 'No clue what kind of atmosphere file this could be')
    end select
    if (get_arg(.false., lverbose)) call print_tenstr_atm(atm, 1_iintegers)
  end subroutine

  subroutine perturb_atmospheres(comm, Npert_atmo, atm, ierr, lverbose)
    integer(mpiint), intent(in) :: comm
    integer(iintegers), intent(in) :: Npert_atmo
    type(t_tenstr_atm), target, intent(inout) :: atm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose

    real(ireals), pointer :: ptr(:, :, :)
    real(ireals), allocatable :: rnd(:, :)
    integer(iintegers) :: Nx, Ny, kmin, kmax, ke, ke1
    integer(iintegers) :: icld, ipert

    integer(iintegers) :: h2o, co2, o3, ch4
    real(ireals), parameter :: H2O_pert(*) = [real(ireals) :: .01, .1, 1., 2., 10]
    real(ireals), parameter :: CO2_pert(*) = [real(ireals) :: .1, 1., 10.]
    real(ireals), parameter :: O3_pert(*) = [real(ireals) :: .1, 1., 10.]
    real(ireals), parameter :: CH4_pert(*) = [real(ireals) :: .1, 1., 10.]

    ierr = 0

    ke1 = size(atm%plev, dim=1)
    ke = ke1 - 1
    Nx = size(atm%plev, dim=2) / Npert_atmo
    Ny = Npert_atmo
    ipert = 1

    do h2o = 1, size(H2O_pert)
      do co2 = 1, size(CO2_pert)
        do o3 = 1, size(O3_pert)
          do ch4 = 1, size(CH4_pert)
            ipert = ipert + 1
            call do_gas_pert(ipert, h2o, co2, o3, ch4)

            do icld = 0, 1 ! ipert from 8 to 8 + 4*3 -1 = 19
              ! water cloud perturbations
              ipert = ipert + 1
              call do_gas_pert(ipert, h2o, co2, o3, ch4)
              kmin = int(find_real_location(atm%plev(:, 1), 850._ireals + real(icld, ireals) * 50._ireals), iintegers)
              kmax = int(find_real_location(atm%plev(:, 1), 800._ireals + real(icld, ireals) * 50._ireals), iintegers)
              ptr(1:ke, 1:Nx, 1:Ny) => atm%lwc; ptr(kmin:kmax, :, ipert) = .1
              ptr(1:ke, 1:Nx, 1:Ny) => atm%reliq; ptr(kmin:kmax, :, ipert) = 10

              ! ice cloud perturbations
              ipert = ipert + 1
              call do_gas_pert(ipert, h2o, co2, o3, ch4)
              kmin = int(find_real_location(atm%plev(:, 1), 300._ireals + real(icld, ireals) * 50._ireals), iintegers)
              kmax = int(find_real_location(atm%plev(:, 1), 200._ireals + real(icld, ireals) * 50._ireals), iintegers)
              ptr(1:ke, 1:Nx, 1:Ny) => atm%iwc; ptr(kmin:kmax, :, ipert) = .01
              ptr(1:ke, 1:Nx, 1:Ny) => atm%reice; ptr(kmin:kmax, :, ipert) = 20

              ! water + ice cloud perturbations
              ipert = ipert + 1
              call do_gas_pert(ipert, h2o, co2, o3, ch4)
              kmin = int(find_real_location(atm%plev(:, 1), 850._ireals + real(icld, ireals) * 50._ireals), iintegers)
              kmax = int(find_real_location(atm%plev(:, 1), 800._ireals + real(icld, ireals) * 50._ireals), iintegers)
              ptr(1:ke, 1:Nx, 1:Ny) => atm%lwc; ptr(kmin:kmax, :, ipert) = .1
              ptr(1:ke, 1:Nx, 1:Ny) => atm%reliq; ptr(kmin:kmax, :, ipert) = 10
              kmin = int(find_real_location(atm%plev(:, 1), 300._ireals + real(icld, ireals) * 50._ireals), iintegers)
              kmax = int(find_real_location(atm%plev(:, 1), 200._ireals + real(icld, ireals) * 50._ireals), iintegers)
              ptr(1:ke, 1:Nx, 1:Ny) => atm%iwc; ptr(kmin:kmax, :, ipert) = .01
              ptr(1:ke, 1:Nx, 1:Ny) => atm%reice; ptr(kmin:kmax, :, ipert) = 20
            end do
          end do
        end do
      end do
    end do

    ! Random H2O perturbations
    allocate (rnd(1:ke, 1:Nx))
    call random_number(rnd)
    ipert = ipert + 1; ptr(1:ke, 1:Nx, 1:Ny) => atm%h2o_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * (rnd * 5)

    print *, 'Nr. of atmosphere perturbation entries:', ipert, '/', Npert_atmo
    call CHKERR(int(ipert - Npert_atmo, mpiint), 'havent used all perturbation slots')

  contains
    subroutine do_gas_pert(ipert, h2o, co2, o3, ch4)
      integer(iintegers), intent(in) :: ipert, h2o, co2, o3, ch4
      ptr(1:ke, 1:Nx, 1:Ny) => atm%h2o_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * H2O_pert(h2o)
      ptr(1:ke, 1:Nx, 1:Ny) => atm%co2_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * CO2_pert(co2)
      ptr(1:ke, 1:Nx, 1:Ny) => atm%o3_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * O3_pert(o3)
      ptr(1:ke, 1:Nx, 1:Ny) => atm%ch4_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * CH4_pert(ch4)
    end subroutine
  end subroutine

  subroutine monochrome_solve(comm, solver, atm, lsolar, repwvl_data, iwvl, edn, eup, abso, ierr, edir)
    integer(mpiint), intent(in) :: comm
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in) :: atm
    logical, intent(in) :: lsolar
    type(t_repwvl_data), intent(in) :: repwvl_data
    integer(iintegers), intent(in) :: iwvl
    real(ireals), allocatable, dimension(:, :, :), intent(inout) :: edn, eup, abso ! [nlyr(+1), local_nx, local_ny ]
    integer(mpiint), intent(out) :: ierr
    real(ireals), allocatable, dimension(:, :, :), intent(inout), optional :: edir

    real(ireals), parameter :: edirTOA = 1._ireals
    real(ireals), allocatable, dimension(:, :, :) :: kabs, ksca, kg ! [nlyr, local_nx, local_ny]
    real(ireals), allocatable :: Blev(:, :, :), Bsrfc(:, :)

    real(ireals) :: albedo
    integer(iintegers) :: k, i, j, icol

    ierr = 0

    allocate (kabs(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym), source=-999._ireals)
    allocate (ksca(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym), source=-999._ireals)
    allocate (kg(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym), source=-999._ireals)
    if (.not. lsolar) then
      allocate (Blev(solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym), source=-999._ireals)
      allocate (Bsrfc(solver%C_one1%xm, solver%C_one1%ym), source=-999._ireals)
    end if

    call repwvl_optprop(repwvl_data, atm, lsolar, iwvl, kabs, ksca, kg, ierr); call CHKERR(ierr)

    if (lsolar) then
      albedo = .15
      call set_optical_properties( &
        & solver,                  &
        & albedo,                  &
        & kabs,                    &
        & ksca,                    &
        & kg)
    else
      albedo = .03

      do j = 1, solver%C_one%ym
        do i = 1, solver%C_one%xm
          icol = i + (j - 1_iintegers) * solver%C_one%xm
          do k = 1, solver%C_one1%zm
            Blev(solver%C_one1%zm + 1 - k, i, j) = planck(repwvl_data%wvls(iwvl) * 1e-9_ireals, atm%tlev(k, icol)) * 1e-9_ireals
          end do
        end do
      end do

      call set_optical_properties( &
        & solver,                  &
        & albedo,                  &
        & kabs,                    &
        & ksca,                    &
        & kg,                      &
        & planck=Blev)
    end if

    call solve_pprts(solver, lthermal=lsolar .eqv. .false., lsolar=lsolar, edirTOA=edirTOA)

    call pprts_get_result(solver, edn, eup, abso, edir)
  end subroutine

  subroutine solve_scene(comm, solver, atm, repwvl_data, edn, eup, abso, ierr, edir)
    integer(mpiint), intent(in) :: comm
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in) :: atm
    type(t_repwvl_data), intent(in) :: repwvl_data
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl ]
    integer(mpiint), intent(out) :: ierr
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout), optional :: edir

    real(ireals), allocatable, dimension(:, :, :) :: spec_edir, spec_abso ! [nlyr(+1), nx, ny ]
    real(ireals), allocatable, dimension(:, :, :) :: spec_edn, spec_eup   ! [nlyr(+1), nx, ny ]

    integer(iintegers) :: iwvl
    integer(mpiint) :: myid

    logical :: lsolar
    lsolar = present(edir)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call set_angles(solver, sundir=spherical_2_cartesian(0._ireals, 60._ireals))

    ! solar part
    if(lboundary_flx_only) then
      if (lsolar) allocate (edir(2, solver%C_dir%xm, solver%C_dir%ym, repwvl_data%Nwvl))
      allocate (edn(2, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
      allocate (eup(2, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
    else
      if (lsolar) allocate (edir(solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym, repwvl_data%Nwvl))
      allocate (edn(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
      allocate (eup(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
    endif
    allocate (abso(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym, repwvl_data%Nwvl), source=0._ireals)

    do iwvl = 1, size(repwvl_data%wvls)

      if (myid .eq. 0) print *, 'Computing wavelengths '//toStr(iwvl)//' / '//toStr(size(repwvl_data%wvls))//&
        & ' -- '//toStr(100._ireals * real(iwvl, ireals) / real(size(repwvl_data%wvls), ireals))//' %'// &
        & ' ('//toStr(repwvl_data%wvls(iwvl))//' nm)'

      if (lsolar) then
        call monochrome_solve(&
          & comm, &
          & solver, &
          & atm, &
          & lsolar, &
          & repwvl_data, &
          & iwvl, &
          & spec_edn, &
          & spec_eup, &
          & spec_abso, &
          & ierr, &
          & spec_edir); call CHKERR(ierr)
      else
        call monochrome_solve(&
          & comm, &
          & solver, &
          & atm, &
          & lsolar, &
          & repwvl_data, &
          & iwvl, &
          & spec_edn, &
          & spec_eup, &
          & spec_abso, &
          & ierr); call CHKERR(ierr)
      end if

      if(lboundary_flx_only) then
        if (lsolar) then
          edir(1, :, :, iwvl) = edir(1, :, :, iwvl) + spec_edir(lbound(spec_edir,1),:,:) * repwvl_data%wgts(iwvl)
          edir(2, :, :, iwvl) = edir(2, :, :, iwvl) + spec_edir(ubound(spec_edir,1),:,:) * repwvl_data%wgts(iwvl)
        endif
        edn(1, :, :, iwvl) = edn(1, :, :, iwvl) + spec_edn(lbound(spec_edn,1),:,:) * repwvl_data%wgts(iwvl)
        edn(2, :, :, iwvl) = edn(2, :, :, iwvl) + spec_edn(ubound(spec_edn,1),:,:) * repwvl_data%wgts(iwvl)
        eup(1, :, :, iwvl) = eup(1, :, :, iwvl) + spec_eup(lbound(spec_eup,1),:,:) * repwvl_data%wgts(iwvl)
        eup(2, :, :, iwvl) = eup(2, :, :, iwvl) + spec_eup(ubound(spec_eup,1),:,:) * repwvl_data%wgts(iwvl)
      else
        if (lsolar) edir(:, :, :, iwvl) = edir(:, :, :, iwvl) + spec_edir * repwvl_data%wgts(iwvl)
        edn(:, :, :, iwvl) = edn(:, :, :, iwvl) + spec_edn * repwvl_data%wgts(iwvl)
        eup(:, :, :, iwvl) = eup(:, :, :, iwvl) + spec_eup * repwvl_data%wgts(iwvl)
      endif
      abso(:, :, :, iwvl) = abso(:, :, :, iwvl) + spec_abso * repwvl_data%wgts(iwvl)

    end do
  end subroutine

  subroutine init_solver(comm, atm, Npert_atmo, solver, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_tenstr_atm), intent(in) :: atm
    integer(iintegers), intent(in) :: Npert_atmo
    class(t_solver), intent(inout) :: solver
    integer(mpiint), intent(out) :: ierr

    integer(iintegers) :: ke1, ke, xm, ym
    integer(iintegers) :: i, j, icol

    real(ireals), parameter :: dx = 1, dy = 1
    real(ireals), parameter :: sundir(3) = [0, 0, -1]

    real(ireals), allocatable :: dz_t2b(:, :, :) ! dz (t2b := top 2 bottom)

    ierr = 0

    ke1 = size(atm%plev, dim=1)
    ke = ke1 - 1_iintegers

    xm = size(atm%plev, dim=2) / Npert_atmo
    ym = Npert_atmo

    allocate (dz_t2b(size(atm%dz, dim=1), xm, ym))
    do j = 1_iintegers, ym
      do i = 1_iintegers, xm
        icol = i + (j - 1_iintegers) * xm
        dz_t2b(:, i, j) = reverse(atm%dz(:, icol))
      end do
    end do

    call init_pprts(comm, ke, xm, ym, dx, dy, sundir, solver, dz3d=dz_t2b)

  end subroutine

  subroutine write_output_data(out_filename, prefix, rdata, edn, eup, abso, ierr, edir, lverbose)
    character(len=*), intent(in) :: out_filename, prefix
    type(t_repwvl_data), intent(in) :: rdata
    real(ireals), dimension(:, :, :, :), intent(in) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl]
    integer(mpiint), intent(out) :: ierr
    real(ireals), dimension(:, :, :, :), intent(in), optional :: edir
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2), dimnames(4)

    dimnames(1) = trim(prefix)//'nlev'
    dimnames(2) = trim(prefix)//'nx'
    dimnames(3) = trim(prefix)//'ny'
    dimnames(4) = trim(prefix)//'wvl'

    groups(1) = trim(out_filename)
    groups(2) = trim(prefix)//'edn'
    call ncwrite(groups, edn, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    groups(2) = trim(prefix)//'eup'
    call ncwrite(groups, eup, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    if (present(edir)) then
      groups(2) = trim(prefix)//'edir'
      call ncwrite(groups, edir, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if

    dimnames(1) = trim(prefix)//'nlay'
    groups(2) = trim(prefix)//'abso'
    call ncwrite(groups, abso, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

    dimnames(1) = trim(prefix)//'wvl'
    groups(2) = trim(prefix)//'wvl'
    call ncwrite(groups, rdata%wvls, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
  end subroutine

  subroutine compute_training_data(comm, bg_atm_filename, atm_filename, out_filename, lsolar, edn, eup, abso, ierr, lverbose, edir)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: bg_atm_filename, atm_filename, out_filename
    logical, intent(in) :: lsolar
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl]
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout), optional :: edir

    type(t_tenstr_atm) :: atm
    class(t_solver), allocatable :: solver
    type(t_repwvl_data), allocatable :: repwvl_data_solar, repwvl_data_thermal

    integer(iintegers), parameter :: Npert_atmo = 947 ! we use the y axis to perturb the loaded atmospheres
    integer(iintegers) :: ipert

    ierr = 0

    call load_atmospheres(comm, bg_atm_filename, atm_filename, Npert_atmo, atm, ierr, lverbose)
    call perturb_atmospheres(comm, Npert_atmo, atm, ierr, lverbose)

    call allocate_pprts_solver_from_commandline(solver, '2str', ierr); call CHKERR(ierr)

    call repwvl_init(repwvl_data_solar, repwvl_data_thermal, ierr, lverbose=.false.); call CHKERR(ierr)
    if (lverbose) print *, 'repwvl_data_thermal Nwvl', size(repwvl_data_thermal%wvls)
    if (lverbose) print *, 'repwvl_data_solar Nwvl', size(repwvl_data_solar%wvls)
    call mie_tables_init(comm, ierr, lverbose=.false.)
    call fu_ice_init(comm, ierr, lverbose=.false.); call CHKERR(ierr)

    call init_solver(comm, atm, Npert_atmo, solver, ierr)

    if (lsolar) then
      call solve_scene(comm, solver, atm, repwvl_data_solar, edn, eup, abso, ierr, edir); call CHKERR(ierr)
      call write_output_data(out_filename, 'solar_', repwvl_data_solar, &
        & edn, eup, abso, ierr, edir, lverbose); call CHKERR(ierr)
      do ipert = 1, Npert_atmo
        print *, 'edir @ srfc ( perturbation='//toStr(ipert)//')', sum(edir(solver%C_dir%ze, :, ipert, :), dim=2)
        print *, 'edn  @ srfc ( perturbation='//toStr(ipert)//')', sum(edn(solver%C_diff%ze, :, ipert, :), dim=2)
      end do
    else
      call solve_scene(comm, solver, atm, repwvl_data_thermal, edn, eup, abso, ierr); call CHKERR(ierr)
      call write_output_data(out_filename, 'thermal_', repwvl_data_thermal, &
        & edn, eup, abso, ierr, lverbose=lverbose); call CHKERR(ierr)
      do ipert = 1, Npert_atmo
        print *, 'edn @ srfc ( perturbation='//toStr(ipert)//')', sum(edn(solver%C_diff%ze, :, ipert, :), dim=2)
        print *, 'eup @ srfc ( perturbation='//toStr(ipert)//')', sum(eup(solver%C_diff%ze, :, ipert, :), dim=2)
      end do
    end if
  end subroutine
end module

program main
  use mpi, only: MPI_COMM_WORLD

  use m_data_parameters, only: &
    & default_str_len, &
    & finalize_mpi, &
    & init_mpi_data_parameters, &
    & mpiint, &
    & ireals

  use m_helper_functions, only: CHKERR, get_petsc_opt
  use m_tenstream_options, only: read_commandline_options

  use m_compute_training_data, only: compute_training_data

  integer(mpiint) :: comm, ierr
  character(len=default_str_len) :: bg_atm_filename, atm_filename, out_filename
  logical :: lflg, lverbose, lsolar, lthermal

  real(ireals), allocatable, dimension(:, :, :, :) :: edir, edn, eup, abso

  comm = MPI_COMM_WORLD
  call init_mpi_data_parameters(comm)
  call read_commandline_options(comm)

  bg_atm_filename = 'atm.dat'
  call get_petsc_opt('', '-bg_atm', bg_atm_filename, lflg, ierr); call CHKERR(ierr)

  atm_filename = 'garand_profiles.nc'
  call get_petsc_opt('', '-atm', atm_filename, lflg, ierr); call CHKERR(ierr)

  out_filename = 'repwvl_training_flx.nc'
  call get_petsc_opt('', '-out', out_filename, lflg, ierr); call CHKERR(ierr)

  lverbose = .true.
  call get_petsc_opt('', '-verbose', lverbose, lflg, ierr); call CHKERR(ierr)

  lsolar = .false.
  call get_petsc_opt('', '-solar', lsolar, lflg, ierr); call CHKERR(ierr)

  lthermal = .false.
  call get_petsc_opt('', '-thermal', lthermal, lflg, ierr); call CHKERR(ierr)

  if (lsolar) then
    call compute_training_data(&
      & comm,                  &
      & bg_atm_filename,       &
      & atm_filename,          &
      & out_filename,          &
      & .true.,                &
      & edn, eup, abso,        &
      & ierr,                  &
      & lverbose=lverbose,     &
      & edir=edir); call CHKERR(ierr)
  else
    print *, 'Not computing solar part: enable with -solar'
  end if

  if (lthermal) then
    call compute_training_data(&
      & comm,                  &
      & bg_atm_filename,       &
      & atm_filename,          &
      & out_filename,          &
      & .false.,               &
      & edn, eup, abso,        &
      & ierr,                  &
      & lverbose); call CHKERR(ierr)
  else
    print *, 'Not computing thermal part: enable with -thermal'
  end if

  call finalize_mpi(comm, .true., .true.)
end program
