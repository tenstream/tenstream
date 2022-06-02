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

module m_compute_repwvl_training_data
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
    & imp_bcast, &
    & domain_decompose_2d_petsc, &
    & imp_min_mean_max, &
    & resize_arr, &
    & reverse, &
    & spherical_2_cartesian, &
    & toStr

  use m_search, only: find_real_location
  use m_netcdfIO, only: ncload, ncwrite
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, print_tenstr_atm, planck
  use m_pprts_base, only: t_coord, t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, pprts_get_result, set_angles, gather_all_toZero

  use m_repwvl_base, only: t_repwvl_data, repwvl_init, repwvl_dtau
  use m_repwvl_pprts, only: repwvl_optprop
  use m_mie_tables, only: mie_tables_init, t_mie_table, mie_optprop, destroy_mie_table
  use m_fu_ice, only: fu_ice_init
  use m_rayleigh, only: rayleigh

  implicit none

  logical, parameter :: lboundary_flx_only = .true.

contains

  subroutine load_atmospheres(comm, bg_atm_filename, atm_filename, atm, Nx, Ny, nxproc, nyproc, ierr, lverbose)
    ! Garand profiles (9 columns)
    ! p, T, z, H2O, O3, CO2, N2O, CO, CH4
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: bg_atm_filename, atm_filename
    type(t_tenstr_atm), intent(out) :: atm
    integer(iintegers), intent(out) :: Nx, Ny
    integer(iintegers), allocatable, intent(out) :: nxproc(:), nyproc(:)
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose

    integer(mpiint) :: myid
    integer(iintegers) :: ke, ke1, Natm
    real(ireals), allocatable :: atmdat(:, :, :) ! dim (cols, nlev, Natm)
    character(len=default_str_len) :: groups(2)

    integer(iintegers) :: xs, ys

    integer(iintegers) :: i, j, icol, icolatm
    real(ireals), allocatable, dimension(:, :) :: plev, tlev, h2ovmr, o3vmr, co2vmr, n2ovmr, ch4vmr

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if (myid .eq. 0) then
      groups(1) = trim(atm_filename)
      groups(2) = 'Atmosphere'
      if (get_arg(.false., lverbose)) then
        print *, 'Reading Atmospheres from '//trim(groups(1))//' : '//trim(groups(2))
      end if
      call ncload(groups, atmdat, ierr, lverbose); call CHKERR(ierr)
      if (get_arg(.false., lverbose)) then
        print *, 'Shape Atmospheres: (tracer, nlev, Natm)', shape(atmdat)
      end if
    end if

    call imp_bcast(comm, atmdat, ierr); call CHKERR(ierr)

    ke1 = size(atmdat, dim=2)
    ke = ke1 - 1_iintegers
    Natm = size(atmdat, dim=3)

    call domain_decompose_2d_petsc(comm, Natm, 1_iintegers, Nx, Ny, xs, ys, nxproc, nyproc, ierr); call CHKERR(ierr)
    if (get_arg(.false., lverbose)) &
      & print *, 'domain decomp Nx', Nx, 'Ny', Ny, 'nxproc', nxproc, 'nyproc', nyproc, 'xs/ys', xs, ys

    select case (size(atmdat, dim=1))
    case (9) ! Garand

      allocate (plev(1:ke1, Nx * Ny))
      allocate (tlev(1:ke1, Nx * Ny))
      allocate (h2ovmr(1:ke, Nx * Ny))
      allocate (o3vmr(1:ke, Nx * Ny))
      allocate (co2vmr(1:ke, Nx * Ny))
      allocate (n2ovmr(1:ke, Nx * Ny))
      allocate (ch4vmr(1:ke, Nx * Ny))

      do j = ys + 1, ys + Ny
        do i = xs + 1, xs + Nx
          icolatm = i + (j - 1_iintegers) * Natm
          icol = (i - xs) + ((j - ys) - 1_iintegers) * Nx

          plev(:, icol) = atmdat(1, :, icolatm) * 1e-2_ireals
          tlev(:, icol) = atmdat(2, :, icolatm)

          h2ovmr(:, icol) = (atmdat(4, 1:ke, icolatm) + atmdat(4, 2:ke1, icolatm))*.5_ireals
          o3vmr(:, icol) = (atmdat(5, 1:ke, icolatm) + atmdat(5, 2:ke1, icolatm))*.5_ireals
          co2vmr(:, icol) = (atmdat(6, 1:ke, icolatm) + atmdat(6, 2:ke1, icolatm))*.5_ireals
          n2ovmr(:, icol) = (atmdat(7, 1:ke, icolatm) + atmdat(7, 2:ke1, icolatm))*.5_ireals
          ch4vmr(:, icol) = (atmdat(9, 1:ke, icolatm) + atmdat(9, 2:ke1, icolatm))*.5_ireals
        end do
      end do

      call setup_tenstr_atm(            &
        & comm,                         &
        & lTOA_to_srfc=.false.,         &
        & atm_filename=bg_atm_filename, &
        & d_plev=plev,                  &
        & d_tlev=tlev,                  &
        & atm=atm,                      &
        & d_h2ovmr=h2ovmr,              &
        & d_o3vmr=o3vmr,                &
        & d_co2vmr=co2vmr,              &
        & d_n2ovmr=n2ovmr,              &
        & d_ch4vmr=ch4vmr)
    case default
      call CHKERR(1_mpiint, 'No clue what kind of atmosphere file this could be')
    end select
    if (myid .eq. 0 .and. get_arg(.false., lverbose)) call print_tenstr_atm(atm, 1_iintegers)
  end subroutine

  subroutine perturb_atmospheres(comm, Npert_atmo, atm, ierr, lverbose, ldryrun)
    integer(mpiint), intent(in) :: comm
    integer(iintegers), intent(inout) :: Npert_atmo
    type(t_tenstr_atm), target, intent(inout) :: atm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose, ldryrun

    real(ireals), pointer :: ptr(:, :, :)
    real(ireals), allocatable :: rnd(:, :)
    integer(iintegers) :: Nx, Ny, kmin, kmax, ke, ke1
    integer(iintegers) :: icld, ipert
    integer(mpiint) :: myid
    logical :: run

    integer(iintegers) :: h2o, co2, o3, ch4
    real(ireals), parameter :: H2O_pert(*) = [real(ireals) :: 0.0001, .01, .1, .5, 1., 2., 10, 20]
    real(ireals), parameter :: CO2_pert(*) = [real(ireals) :: .1, 1., 10.]
    real(ireals), parameter :: O3_pert(*) = [real(ireals) :: .1, 1., 10.]
    real(ireals), parameter :: CH4_pert(*) = [real(ireals) :: .1, 1., 10.]
    !real(ireals), parameter :: H2O_pert(*) = [real(ireals) :: 1.]
    !real(ireals), parameter :: CO2_pert(*) = [real(ireals) :: 1.]
    !real(ireals), parameter :: O3_pert(*) = [real(ireals) :: 1.]
    !real(ireals), parameter :: CH4_pert(*) = [real(ireals) :: 1.]

    ierr = 0

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    ipert = 1
    run = .not. get_arg(.false., ldryrun)
    if (run) then
      ke1 = size(atm%plev, dim=1)
      ke = ke1 - 1
      Nx = size(atm%plev, dim=2)
      Ny = Npert_atmo

      if (allocated(atm%plev)) call resize_arr(Nx * Npert_atmo, atm%plev, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%tlev)) call resize_arr(Nx * Npert_atmo, atm%tlev, dim=2_mpiint, lrepeat=.true.)

      if (allocated(atm%zm)) call resize_arr(Nx * Npert_atmo, atm%zm, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%dz)) call resize_arr(Nx * Npert_atmo, atm%dz, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%tlay)) call resize_arr(Nx * Npert_atmo, atm%tlay, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%h2o_lay)) call resize_arr(Nx * Npert_atmo, atm%h2o_lay, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%o3_lay)) call resize_arr(Nx * Npert_atmo, atm%o3_lay, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%co2_lay)) call resize_arr(Nx * Npert_atmo, atm%co2_lay, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%ch4_lay)) call resize_arr(Nx * Npert_atmo, atm%ch4_lay, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%n2o_lay)) call resize_arr(Nx * Npert_atmo, atm%n2o_lay, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%o2_lay)) call resize_arr(Nx * Npert_atmo, atm%o2_lay, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%lwc)) call resize_arr(Nx * Npert_atmo, atm%lwc, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%reliq)) call resize_arr(Nx * Npert_atmo, atm%reliq, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%iwc)) call resize_arr(Nx * Npert_atmo, atm%iwc, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%reice)) call resize_arr(Nx * Npert_atmo, atm%reice, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%opt_tau)) call resize_arr(Nx * Npert_atmo, atm%opt_tau, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%opt_w0)) call resize_arr(Nx * Npert_atmo, atm%opt_w0, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%opt_g)) call resize_arr(Nx * Npert_atmo, atm%opt_g, dim=2_mpiint, lrepeat=.true.)
      if (allocated(atm%tskin)) call resize_arr(Nx * Npert_atmo, atm%tskin, lrepeat=.true.)
    end if

    do h2o = 1, size(H2O_pert)
      do co2 = 1, size(CO2_pert)
        do o3 = 1, size(O3_pert)
          do ch4 = 1, size(CH4_pert)

            ipert = ipert + 1
            if (run) call do_gas_pert(ipert, h2o, co2, o3, ch4)

            do icld = 0, 2 ! ipert from 8 to 8 + 4*3 -1 = 19
              ! water cloud perturbations
              ipert = ipert + 1
              if (run) then
                call do_gas_pert(ipert, h2o, co2, o3, ch4)
                kmin = int(find_real_location(atm%plev(:, 1), 850._ireals + real(icld, ireals) * 50._ireals), iintegers)
                kmax = int(find_real_location(atm%plev(:, 1), 800._ireals + real(icld, ireals) * 50._ireals), iintegers)
                ptr(1:ke, 1:Nx, 1:Ny) => atm%lwc; ptr(kmin:kmax, :, ipert) = .1
                ptr(1:ke, 1:Nx, 1:Ny) => atm%reliq; ptr(kmin:kmax, :, ipert) = 10
              end if

              ! ice cloud perturbations
              ipert = ipert + 1
              if (run) then
                call do_gas_pert(ipert, h2o, co2, o3, ch4)
                kmin = int(find_real_location(atm%plev(:, 1), 300._ireals + real(icld, ireals) * 50._ireals), iintegers)
                kmax = int(find_real_location(atm%plev(:, 1), 200._ireals + real(icld, ireals) * 50._ireals), iintegers)
                ptr(1:ke, 1:Nx, 1:Ny) => atm%iwc; ptr(kmin:kmax, :, ipert) = .01
                ptr(1:ke, 1:Nx, 1:Ny) => atm%reice; ptr(kmin:kmax, :, ipert) = 20
              end if

              ! water + ice cloud perturbations
              ipert = ipert + 1
              if (run) then
                call do_gas_pert(ipert, h2o, co2, o3, ch4)
                kmin = int(find_real_location(atm%plev(:, 1), 850._ireals + real(icld, ireals) * 50._ireals), iintegers)
                kmax = int(find_real_location(atm%plev(:, 1), 800._ireals + real(icld, ireals) * 50._ireals), iintegers)
                ptr(1:ke, 1:Nx, 1:Ny) => atm%lwc; ptr(kmin:kmax, :, ipert) = .1
                ptr(1:ke, 1:Nx, 1:Ny) => atm%reliq; ptr(kmin:kmax, :, ipert) = 10
                kmin = int(find_real_location(atm%plev(:, 1), 300._ireals + real(icld, ireals) * 50._ireals), iintegers)
                kmax = int(find_real_location(atm%plev(:, 1), 200._ireals + real(icld, ireals) * 50._ireals), iintegers)
                ptr(1:ke, 1:Nx, 1:Ny) => atm%iwc; ptr(kmin:kmax, :, ipert) = .01
                ptr(1:ke, 1:Nx, 1:Ny) => atm%reice; ptr(kmin:kmax, :, ipert) = 20
              end if
            end do
          end do
        end do
      end do
    end do

    ! Random H2O perturbations
    ipert = ipert + 1

    if (run) then
      allocate (rnd(1:ke, 1:Nx))
      call random_number(rnd)
      ptr(1:ke, 1:Nx, 1:Ny) => atm%h2o_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * (rnd * 5)

      if (myid .eq. 0 .and. get_arg(.false., lverbose)) print *, 'Nr. of atmosphere perturbation entries:', ipert, '/', Npert_atmo
      call CHKERR(int(ipert - Npert_atmo, mpiint), 'havent used all perturbation slots '//toStr(ipert)//' / '//toStr(Npert_atmo))
    else
      Npert_atmo = ipert
    end if

  contains
    subroutine do_gas_pert(ipert, h2o, co2, o3, ch4)
      integer(iintegers), intent(in) :: ipert, h2o, co2, o3, ch4
      ptr(1:ke, 1:Nx, 1:Ny) => atm%h2o_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * H2O_pert(h2o)
      ptr(1:ke, 1:Nx, 1:Ny) => atm%co2_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * CO2_pert(co2)
      ptr(1:ke, 1:Nx, 1:Ny) => atm%o3_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * O3_pert(o3)
      ptr(1:ke, 1:Nx, 1:Ny) => atm%ch4_lay; ptr(:, :, ipert) = ptr(:, :, ipert) * CH4_pert(ch4)
    end subroutine
  end subroutine

  subroutine monochrome_solve(solver, atm, lsolar, repwvl_data, mie_table, iwvl, edn, eup, abso, ierr, edir)
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in) :: atm
    logical, intent(in) :: lsolar
    type(t_repwvl_data), intent(in) :: repwvl_data
    type(t_mie_table), intent(in) :: mie_table
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

    call repwvl_optprop(repwvl_data, atm, mie_table, lsolar, iwvl, kabs, ksca, kg, ierr); call CHKERR(ierr)

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

  subroutine solve_scene(comm, solver, atm, repwvl_data, mie_table, edn, eup, abso, ierr, edir, lverbose)
    integer(mpiint), intent(in) :: comm
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in) :: atm
    type(t_repwvl_data), intent(in) :: repwvl_data
    type(t_mie_table), intent(in) :: mie_table
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl ]
    integer(mpiint), intent(out) :: ierr
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout), optional :: edir
    logical, intent(in), optional :: lverbose

    real(ireals), allocatable, dimension(:, :, :) :: spec_edir, spec_abso ! [nlyr(+1), nx, ny ]
    real(ireals), allocatable, dimension(:, :, :) :: spec_edn, spec_eup   ! [nlyr(+1), nx, ny ]

    integer(iintegers) :: iwvl
    integer(mpiint) :: myid

    logical :: lsolar
    lsolar = present(edir)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call set_angles(solver, sundir=spherical_2_cartesian(0._ireals, 60._ireals))

    ! solar part
    if (lboundary_flx_only) then
      if (lsolar) allocate (edir(2, solver%C_dir%xm, solver%C_dir%ym, repwvl_data%Nwvl))
      allocate (edn(2, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
      allocate (eup(2, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
    else
      if (lsolar) allocate (edir(solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym, repwvl_data%Nwvl))
      allocate (edn(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
      allocate (eup(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
    end if
    allocate (abso(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym, repwvl_data%Nwvl), source=0._ireals)

    do iwvl = 1, size(repwvl_data%wvls)

      if (myid .eq. 0 .and. get_arg(.false., lverbose)) &
        & print *, 'Computing wavelengths '//toStr(iwvl)//' / '//toStr(size(repwvl_data%wvls))//&
        & ' -- '//toStr(100._ireals * real(iwvl, ireals) / real(size(repwvl_data%wvls), ireals))//' %'// &
        & ' ('//toStr(repwvl_data%wvls(iwvl))//' nm)'

      if (lsolar) then
        call monochrome_solve(&
          & solver, &
          & atm, &
          & lsolar, &
          & repwvl_data, &
          & mie_table, &
          & iwvl, &
          & spec_edn, &
          & spec_eup, &
          & spec_abso, &
          & ierr, &
          & spec_edir); call CHKERR(ierr)
      else
        call monochrome_solve(&
          & solver, &
          & atm, &
          & lsolar, &
          & repwvl_data, &
          & mie_table, &
          & iwvl, &
          & spec_edn, &
          & spec_eup, &
          & spec_abso, &
          & ierr); call CHKERR(ierr)
      end if

      if (lboundary_flx_only) then
        if (lsolar) then
          edir(1, :, :, iwvl) = edir(1, :, :, iwvl) + spec_edir(lbound(spec_edir, 1), :, :) * repwvl_data%wgts(iwvl)
          edir(2, :, :, iwvl) = edir(2, :, :, iwvl) + spec_edir(ubound(spec_edir, 1), :, :) * repwvl_data%wgts(iwvl)
        end if
        edn(1, :, :, iwvl) = edn(1, :, :, iwvl) + spec_edn(lbound(spec_edn, 1), :, :) * repwvl_data%wgts(iwvl)
        edn(2, :, :, iwvl) = edn(2, :, :, iwvl) + spec_edn(ubound(spec_edn, 1), :, :) * repwvl_data%wgts(iwvl)
        eup(1, :, :, iwvl) = eup(1, :, :, iwvl) + spec_eup(lbound(spec_eup, 1), :, :) * repwvl_data%wgts(iwvl)
        eup(2, :, :, iwvl) = eup(2, :, :, iwvl) + spec_eup(ubound(spec_eup, 1), :, :) * repwvl_data%wgts(iwvl)
      else
        if (lsolar) edir(:, :, :, iwvl) = edir(:, :, :, iwvl) + spec_edir * repwvl_data%wgts(iwvl)
        edn(:, :, :, iwvl) = edn(:, :, :, iwvl) + spec_edn * repwvl_data%wgts(iwvl)
        eup(:, :, :, iwvl) = eup(:, :, :, iwvl) + spec_eup * repwvl_data%wgts(iwvl)
      end if
      abso(:, :, :, iwvl) = abso(:, :, :, iwvl) + spec_abso * repwvl_data%wgts(iwvl)

    end do
  end subroutine

  subroutine init_solver(comm, atm, Nx, Ny, nxproc, nyproc, solver, ierr)
    integer(mpiint), intent(in) :: comm
    type(t_tenstr_atm), intent(in) :: atm
    integer(iintegers), intent(in) :: Nx, Ny, nxproc(:), nyproc(:)
    class(t_solver), intent(inout) :: solver
    integer(mpiint), intent(out) :: ierr

    integer(iintegers) :: ke1, ke
    integer(iintegers) :: i, j, icol

    real(ireals), parameter :: dx = 1, dy = 1
    real(ireals), parameter :: sundir(3) = [0, 0, -1]

    real(ireals), allocatable :: dz_t2b(:, :, :) ! dz (t2b := top 2 bottom)

    ierr = 0

    ke1 = size(atm%plev, dim=1)
    ke = ke1 - 1_iintegers

    allocate (dz_t2b(size(atm%dz, dim=1), Nx, Ny))
    do j = 1_iintegers, Ny
      do i = 1_iintegers, Nx
        icol = i + (j - 1_iintegers) * Nx
        dz_t2b(:, i, j) = reverse(atm%dz(:, icol))
      end do
    end do

    call init_pprts(comm, &
      & ke, Nx, Ny, &
      & dx, dy, &
      & sundir, &
      & solver, &
      & dz3d=dz_t2b, &
      & nxproc=nxproc, &
      & nyproc=nyproc)

  end subroutine

  subroutine write_output_data(comm, solver, out_filename, prefix, rdata, edn, eup, abso, ierr, edir, lverbose)
    integer(mpiint), intent(in) :: comm
    class(t_solver), intent(in) :: solver
    character(len=*), intent(in) :: out_filename, prefix
    type(t_repwvl_data), intent(in) :: rdata
    real(ireals), dimension(:, :, :, :), intent(in) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl]
    integer(mpiint), intent(out) :: ierr
    real(ireals), dimension(:, :, :, :), intent(in), optional :: edir
    logical, intent(in), optional :: lverbose

    character(len=default_str_len) :: groups(2), dimnames(4)
    real(ireals), allocatable, dimension(:, :, :, :) :: gflx ! global domain [nlyr(+1), nx, ny, Nwvl]

    integer(mpiint) :: myid

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    dimnames(1) = trim(prefix)//'nlev'
    dimnames(2) = trim(prefix)//'nx'
    dimnames(3) = trim(prefix)//'ny'
    dimnames(4) = trim(prefix)//'wvl'

    groups(1) = trim(out_filename)
    groups(2) = trim(prefix)//'edn'
    if (lboundary_flx_only) then
      call local2global_bounds_only(solver%Csrfc_one, edn, gflx)
    else
      call local2global(solver%C_one_atm1, edn, gflx)
    end if
    if (myid .eq. 0) then
      call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if
    if (allocated(gflx)) deallocate (gflx)

    if (lboundary_flx_only) then
      call local2global_bounds_only(solver%Csrfc_one, eup, gflx)
    else
      call local2global(solver%C_one_atm1, eup, gflx)
    end if
    groups(2) = trim(prefix)//'eup'
    if (myid .eq. 0) then
      call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if
    if (allocated(gflx)) deallocate (gflx)

    if (present(edir)) then
      groups(2) = trim(prefix)//'edir'
      if (lboundary_flx_only) then
        call local2global_bounds_only(solver%Csrfc_one, edir, gflx)
      else
        call local2global(solver%C_one_atm1, edir, gflx)
      end if
      if (myid .eq. 0) then
        call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      if (allocated(gflx)) deallocate (gflx)
    end if

    dimnames(1) = trim(prefix)//'nlay'
    groups(2) = trim(prefix)//'abso'
    call local2global(solver%C_one_atm, abso, gflx)
    if (myid .eq. 0) then
      call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if

    dimnames(1) = trim(prefix)//'wvl'
    groups(2) = trim(prefix)//'wvl'
    if (myid .eq. 0) then
      call ncwrite(groups, rdata%wvls, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if

  contains

    subroutine local2global_bounds_only(C, local, gflx)
      type(t_coord), intent(in) :: C
      real(ireals), intent(in) :: local(:, :, :, :)
      real(ireals), allocatable, intent(out) :: gflx(:, :, :, :) ! global array on rank0
      real(ireals), allocatable :: tmp(:, :, :) ! global array on rank0
      integer(iintegers) :: iwvl, k
      if (myid .eq. 0) then
        allocate (gflx(size(local, dim=1), C%glob_xm, C%glob_ym, size(local, dim=4)))
      end if
      do iwvl = 1, size(local, dim=4)
        do k = 1, 2
          call gather_all_toZero(C, local(k:k, :, :, iwvl), tmp)
          if (myid .eq. 0) gflx(k:k, :, :, iwvl) = tmp
        end do
      end do
    end subroutine
    subroutine local2global(C, local, gflx)
      type(t_coord), intent(in) :: C
      real(ireals), intent(in) :: local(:, :, :, :)
      real(ireals), allocatable, intent(out) :: gflx(:, :, :, :) ! global array on rank0
      real(ireals), allocatable :: tmp(:, :, :) ! global array on rank0
      integer(iintegers) :: iwvl
      if (myid .eq. 0) then
        allocate (gflx(C%glob_zm, C%glob_xm, C%glob_ym, size(local, dim=4)))
      end if
      do iwvl = 1, size(local, dim=4)
        call gather_all_toZero(C, local(:, :, :, iwvl), tmp)
        if (myid .eq. 0) gflx(:, :, :, iwvl) = tmp
      end do
    end subroutine
  end subroutine

  subroutine compute_repwvl_training_data(&
      & comm,            &
      & bg_atm_filename, &
      & atm_filename,    &
      & out_filename,    &
      & lsolar,          &
      & edn,             &
      & eup,             &
      & abso,            &
      & ierr,            &
      & lverbose,        &
      & edir)
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

    integer(iintegers) :: Npert_atmo, ipert
    integer(mpiint) :: myid

    integer(iintegers) :: Nx, Ny ! local domain sizes, we use the y axis to perturb the loaded atmospheres
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    type(t_mie_table), allocatable :: mie_table

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call perturb_atmospheres(comm, Npert_atmo, atm, ierr, lverbose, ldryrun=.true.); call CHKERR(ierr)
    if (myid .eq. 0 .and. lverbose) print *, 'Number of perturbations', Npert_atmo

    call load_atmospheres(comm, bg_atm_filename, atm_filename, atm, Nx, Ny, nxproc, nyproc, ierr, lverbose); call CHKERR(ierr)
    Ny = Npert_atmo
    nyproc(:) = Ny
    call perturb_atmospheres(comm, Ny, atm, ierr, lverbose)

    call allocate_pprts_solver_from_commandline(solver, '2str', ierr); call CHKERR(ierr)

    call repwvl_init(comm, repwvl_data_solar, repwvl_data_thermal, ierr, lverbose=.false.); call CHKERR(ierr)
    if (myid .eq. 0 .and. lverbose) print *, 'repwvl_data_thermal Nwvl', size(repwvl_data_thermal%wvls)
    if (myid .eq. 0 .and. lverbose) print *, 'repwvl_data_solar Nwvl', size(repwvl_data_solar%wvls)
    call mie_tables_init(comm, mie_table, ierr, lverbose=.false.)
    call fu_ice_init(comm, ierr, lverbose=.false.); call CHKERR(ierr)

    call init_solver(comm, atm, Nx, Ny, nxproc, nyproc, solver, ierr)

    if (lsolar) then
      call solve_scene(&
        & comm, &
        & solver, &
        & atm, &
        & repwvl_data_solar, &
        & mie_table, &
        & edn, &
        & eup, &
        & abso, &
        & ierr, &
        & edir, &
        & lverbose=lverbose); call CHKERR(ierr)

      call write_output_data(comm, solver, out_filename, 'solar_', repwvl_data_solar, &
        & edn, eup, abso, ierr, edir, lverbose); call CHKERR(ierr)

      if (myid .eq. 0) then
        do ipert = 1, Ny
          print *, 'edir @ srfc ( perturbation='//toStr(ipert)//')', sum(edir(ubound(edir, 1), :, ipert, :), dim=2)
          print *, 'edn  @ srfc ( perturbation='//toStr(ipert)//')', sum(edn(ubound(edn, 1), :, ipert, :), dim=2)
        end do
      end if
    else

      call solve_scene(&
        & comm, &
        & solver, &
        & atm, &
        & repwvl_data_thermal, &
        & mie_table, &
        & edn, &
        & eup, &
        & abso, &
        & ierr, &
        & lverbose=lverbose); call CHKERR(ierr)

      call write_output_data(comm, solver, out_filename, 'thermal_', repwvl_data_thermal, &
        & edn, eup, abso, ierr, lverbose=lverbose); call CHKERR(ierr)

      if (myid .eq. 0) then
        do ipert = 1, Ny
          print *, 'edn @ srfc ( perturbation='//toStr(ipert)//')', sum(edn(ubound(edn, 1), :, ipert, :), dim=2)
          print *, 'eup @ srfc ( perturbation='//toStr(ipert)//')', sum(eup(ubound(eup, 1), :, ipert, :), dim=2)
        end do
      end if
    end if
    call destroy_mie_table(mie_table, ierr); call CHKERR(ierr)
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

  use m_compute_repwvl_training_data, only: compute_repwvl_training_data

  integer(mpiint) :: comm, ierr
  character(len=default_str_len) :: bg_atm_filename, atm_filename, out_filename
  logical :: lflg, lverbose, lsolar, lthermal

  real(ireals), allocatable, dimension(:, :, :, :) :: edir, edn, eup, abso
  integer(mpiint) :: myid

  comm = MPI_COMM_WORLD
  call init_mpi_data_parameters(comm)
  call read_commandline_options(comm)

  call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

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
    call compute_repwvl_training_data(&
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
    if (myid .eq. 0) print *, 'Not computing solar part: enable with -solar'
  end if

  if (lthermal) then
    call compute_repwvl_training_data(&
      & comm,                  &
      & bg_atm_filename,       &
      & atm_filename,          &
      & out_filename,          &
      & .false.,               &
      & edn, eup, abso,        &
      & ierr,                  &
      & lverbose); call CHKERR(ierr)
  else
    if (myid .eq. 0) print *, 'Not computing thermal part: enable with -thermal'
  end if

  call finalize_mpi(comm, .true., .true.)
end program
