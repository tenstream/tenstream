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
    & domain_decompose_2d_petsc, &
    & get_arg, &
    & get_petsc_opt, &
    & imp_bcast, &
    & imp_min_mean_max, &
    & imp_min_mean_max, &
    & resize_arr, &
    & reverse, &
    & spherical_2_cartesian, &
    & toStr

  use m_search, only: find_real_location
  use m_netcdfIO, only: ncload, ncwrite, get_global_attribute, list_global_attributes
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm, print_tenstr_atm, planck
  use m_pprts_base, only: t_coord, t_solver, allocate_pprts_solver_from_commandline
  use m_pprts, only: init_pprts, set_optical_properties, solve_pprts, pprts_get_result, set_angles, gather_all_toZero

  use m_repwvl_base, only: t_repwvl_data, repwvl_init
  use m_repwvl_optprop, only: repwvl_optprop
  use m_mie_tables, only: mie_tables_init, t_mie_table, mie_optprop, destroy_mie_table
  use m_fu_ice, only: fu_ice_init
  use m_rayleigh, only: rayleigh

  implicit none

  logical :: lboundary_flx_only

contains

  subroutine load_atmospheres(comm, bg_atm_filename, atm_filename, atm, Nx, Ny, nxproc, nyproc, ierr, lverbose)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: bg_atm_filename, atm_filename
    type(t_tenstr_atm), intent(out) :: atm
    integer(iintegers), intent(out) :: Nx, Ny
    integer(iintegers), allocatable, intent(out) :: nxproc(:), nyproc(:)
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose

    integer(mpiint) :: myid
    integer(iintegers) :: input_mode
    character(len=default_str_len) :: attr

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    input_mode = -1

    if (myid .eq. 0) then
      if (get_arg(.false., lverbose)) then
        call list_global_attributes(trim(atm_filename), ierr); call CHKERR(ierr)
      end if
      call get_global_attribute(trim(atm_filename), 'project', attr, ierr)
      if (ierr .ne. 0) input_mode = 1
      if (trim(attr) .eq. 'CKDMIP') input_mode = 2
    end if
    call imp_bcast(comm, input_mode, 0_mpiint, ierr); call CHKERR(ierr)

    select case (input_mode)
    case (1)
      call garand_input()
    case (2)
      call CKDMIP_input()
    case default
      call CHKERR(-1_mpiint, "dont know what atmosphere input file that is supposed to be")
    end select

  contains
    subroutine CKDMIP_input()
      ! CKDMIP profiles
      ! p, T, H2O, O3, CO2, N2O, O2, CH4
      character(len=default_str_len) :: groups(2)
      real(ireals), allocatable, dimension(:, :) :: plev, play, tlev, tlay, h2ovmr, o3vmr, co2vmr, n2ovmr, ch4vmr, o2vmr
      integer(iintegers) :: ke, ke1, Natm
      integer(iintegers) :: xs, ys

      if (myid .eq. 0) then
        groups(1) = trim(atm_filename)
        if (get_arg(.false., lverbose)) then
          print *, 'Reading Atmospheres from '//trim(groups(1))
        end if
        groups(2) = 'pressure_hl'; call ncload(groups, plev, ierr, lverbose); call CHKERR(ierr)
        if (get_arg(.false., lverbose)) then
          print *, 'Shape Atmospheres: (nlev, Natm)', shape(plev)
        end if
        groups(2) = 'pressure_fl'; call ncload(groups, play, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'temperature_hl'; call ncload(groups, tlev, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'temperature_fl'; call ncload(groups, tlay, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'h2o_mole_fraction_fl'; call ncload(groups, h2ovmr, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'o3_mole_fraction_fl'; call ncload(groups, o3vmr, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'co2_mole_fraction_fl'; call ncload(groups, co2vmr, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'ch4_mole_fraction_fl'; call ncload(groups, ch4vmr, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'n2o_mole_fraction_fl'; call ncload(groups, n2ovmr, ierr, lverbose); call CHKERR(ierr)
        groups(2) = 'o2_mole_fraction_fl'; call ncload(groups, o2vmr, ierr, lverbose); call CHKERR(ierr)

        plev = reverse(plev)
        play = reverse(play)
        tlev = reverse(tlev)
        tlay = reverse(tlay)
        h2ovmr = reverse(h2ovmr)
        o3vmr = reverse(o3vmr)
        co2vmr = reverse(co2vmr)
        ch4vmr = reverse(ch4vmr)
        n2ovmr = reverse(n2ovmr)
        o2vmr = reverse(o2vmr)
      end if

      call imp_bcast(comm, plev, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, play, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, tlev, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, tlay, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, h2ovmr, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, o3vmr, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, co2vmr, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, ch4vmr, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, n2ovmr, 0_mpiint, ierr); call CHKERR(ierr)
      call imp_bcast(comm, o2vmr, 0_mpiint, ierr); call CHKERR(ierr)

      ke1 = size(plev, dim=1)
      ke = ke1 - 1_iintegers
      Natm = size(plev, dim=2)

      call domain_decompose_2d_petsc(comm, Natm, 1_iintegers, Nx, Ny, xs, ys, nxproc, nyproc, ierr); call CHKERR(ierr)
      if (get_arg(.false., lverbose)) &
        & print *, 'domain decomp Nx', Nx, 'Ny', Ny, 'nxproc', nxproc, 'nyproc', nyproc, 'xs/ys', xs, ys

      call setup_tenstr_atm(              &
        & comm,                           &
        & lTOA_to_srfc=.false.,           &
        & atm_filename=bg_atm_filename,   &
        & d_plev=plev(:, xs + 1:xs + Nx) * 1e-2, &
        & d_tlev=tlev(:, xs + 1:xs + Nx),      &
        & atm=atm,                        &
        & d_h2ovmr=h2ovmr(:, xs + 1:xs + Nx),  &
        & d_o3vmr=o3vmr(:, xs + 1:xs + Nx),  &
        & d_co2vmr=co2vmr(:, xs + 1:xs + Nx),  &
        & d_ch4vmr=ch4vmr(:, xs + 1:xs + Nx),  &
        & d_n2ovmr=n2ovmr(:, xs + 1:xs + Nx),  &
        & d_o2vmr=o2vmr(:, xs + 1:xs + Nx))

      if (myid .eq. 0 .and. get_arg(.false., lverbose)) call print_tenstr_atm(atm, 1_iintegers)
    end subroutine

    subroutine garand_input()
      ! Garand profiles (9 columns)
      ! p, T, z, H2O, O3, CO2, N2O, CO, CH4
      integer(iintegers) :: ke, ke1, Natm
      real(ireals), allocatable :: atmdat(:, :, :) ! dim (cols, nlev, Natm)
      character(len=default_str_len) :: groups(2)

      integer(iintegers) :: xs, ys

      integer(iintegers) :: i, j, icol, icolatm
      real(ireals), allocatable, dimension(:, :) :: plev, tlev, h2ovmr, o3vmr, co2vmr, n2ovmr, ch4vmr
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

      call imp_bcast(comm, atmdat, 0_mpiint, ierr); call CHKERR(ierr)

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
  end subroutine

  subroutine perturb_atmospheres(comm, lperturb, Npert_atmo, atm, ierr, lverbose, ldryrun)
    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lperturb
    integer(iintegers), intent(inout) :: Npert_atmo
    type(t_tenstr_atm), target, intent(inout) :: atm
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose, ldryrun

    real(ireals), pointer :: ptr(:, :, :)
    integer(iintegers) :: Nx, Ny, kmin, kmax, ke, ke1
    integer(iintegers) :: icld, ipert
    integer(mpiint) :: myid
    logical :: run

    integer(iintegers) :: h2o, co2, o3, ch4
    !real(ireals), parameter :: H2O_pert(*) = [real(ireals) :: .1, .5, 1., 2., 10.]
    !real(ireals), parameter :: CO2_pert(*) = [real(ireals) :: .5, 1., 2.]
    !real(ireals), parameter :: O3_pert(*) = [real(ireals) :: .5, 1., 2.]
    !real(ireals), parameter :: CH4_pert(*) = [real(ireals) :: .5, 1., 2.]
    real(ireals), parameter :: H2O_pert(*) = [real(ireals) :: .5, 1., 2.]
    real(ireals), parameter :: CO2_pert(*) = [real(ireals) :: 1.]
    real(ireals), parameter :: O3_pert(*) = [real(ireals) :: 1.]
    real(ireals), parameter :: CH4_pert(*) = [real(ireals) :: 1.]

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

    if (lperturb) then
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
                  print *, 'kmin/max', kmin, kmax, 850._ireals + real(icld, ireals) * 50._ireals
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%lwc; ptr(kmin:kmax, :, ipert) = 0.1 * (icld + 1)
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%reliq; ptr(kmin:kmax, :, ipert) = 10
                end if

                ! ice cloud perturbations
                ipert = ipert + 1
                if (run) then
                  call do_gas_pert(ipert, h2o, co2, o3, ch4)
                  kmin = int(find_real_location(atm%plev(:, 1), 300._ireals + real(icld, ireals) * 50._ireals), iintegers)
                  kmax = int(find_real_location(atm%plev(:, 1), 200._ireals + real(icld, ireals) * 50._ireals), iintegers)
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%iwc; ptr(kmin:kmax, :, ipert) = 0.05 * (icld + 1)
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%reice; ptr(kmin:kmax, :, ipert) = 20
                end if

                ! water + ice cloud perturbations
                ipert = ipert + 1
                if (run) then
                  call do_gas_pert(ipert, h2o, co2, o3, ch4)
                  kmin = int(find_real_location(atm%plev(:, 1), 750._ireals + real(icld, ireals) * 50._ireals), iintegers)
                  kmax = int(find_real_location(atm%plev(:, 1), 700._ireals + real(icld, ireals) * 50._ireals), iintegers)
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%lwc; ptr(kmin:kmax, :, ipert) = 0.2 * (icld + 1)
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%reliq; ptr(kmin:kmax, :, ipert) = 20

                  kmin = int(find_real_location(atm%plev(:, 1), 200._ireals + real(icld, ireals) * 50._ireals), iintegers)
                  kmax = int(find_real_location(atm%plev(:, 1), 100._ireals + real(icld, ireals) * 50._ireals), iintegers)
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%iwc; ptr(kmin:kmax, :, ipert) = 0.1 * (icld + 1)
                  ptr(1:ke, 1:Nx, 1:Ny) => atm%reice; ptr(kmin:kmax, :, ipert) = 40
                end if
              end do
            end do
          end do
        end do
      end do

      if (run) then
        if (myid .eq. 0 .and. get_arg(.false., lverbose)) print *, 'Nr. of atmosphere perturbation entries:', ipert, '/', Npert_atmo
        call CHKERR(int(ipert - Npert_atmo, mpiint), 'havent used all perturbation slots '//toStr(ipert)//' / '//toStr(Npert_atmo))
      end if
    end if

    if (.not. run) then
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

    do j = 1, solver%C_one%ym
      do i = 1, solver%C_one%xm
        icol = i + (j - 1_iintegers) * solver%C_one%xm

        do k = 1, solver%C_one%zm
          call repwvl_optprop(&
            & repwvl_data, atm, mie_table, &
            & lsolar, k, icol, iwvl, &
            & kabs(size(kabs, dim=1) + 1 - k, i, j), &
            & ksca(size(ksca, dim=1) + 1 - k, i, j), &
            & kg(size(kg, dim=1) + 1 - k, i, j), &
            & ierr); call CHKERR(ierr)
        end do
      end do
    end do

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

  subroutine solve_scene(comm, solver, atm, repwvl_data, mie_table, lspectral_output, &
      & edn, eup, abso, edn_int, eup_int, abso_int, ierr, edir, edir_int, lverbose)
    integer(mpiint), intent(in) :: comm
    class(t_solver), intent(inout) :: solver
    type(t_tenstr_atm), intent(in) :: atm
    type(t_repwvl_data), intent(in) :: repwvl_data
    type(t_mie_table), intent(in) :: mie_table
    logical, intent(in) :: lspectral_output
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl ]
    real(ireals), allocatable, dimension(:, :, :), intent(inout) :: edn_int, eup_int, abso_int ! [nlyr(+1), nx, ny ]
    integer(mpiint), intent(out) :: ierr
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout), optional :: edir
    real(ireals), allocatable, dimension(:, :, :), intent(inout), optional :: edir_int
    logical, intent(in), optional :: lverbose

    real(ireals), allocatable, dimension(:, :, :) :: spec_edir, spec_abso ! [nlyr(+1), nx, ny ]
    real(ireals), allocatable, dimension(:, :, :) :: spec_edn, spec_eup   ! [nlyr(+1), nx, ny ]

    integer(iintegers) :: iwvl
    integer(mpiint) :: myid

    logical :: lsolar
    lsolar = present(edir)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call set_angles(solver, sundir=spherical_2_cartesian(0._ireals, 60._ireals))

    if (lspectral_output) then
      if (lboundary_flx_only) then
        if (lsolar) allocate (edir(2, solver%C_dir%xm, solver%C_dir%ym, repwvl_data%Nwvl))
        allocate (edn(2, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
        allocate (eup(2, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
      else
        if (lsolar) allocate (edir(solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym, repwvl_data%Nwvl))
        allocate (edn(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
        allocate (eup(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym, repwvl_data%Nwvl))
      end if
      allocate (abso(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym, repwvl_data%Nwvl))

      if (lsolar) edir = 0
      edn = 0
      eup = 0
      abso = 0
    end if

    if (lboundary_flx_only) then
      if (lsolar) allocate (edir_int(2, solver%C_dir%xm, solver%C_dir%ym))
      allocate (edn_int(2, solver%C_diff%xm, solver%C_diff%ym))
      allocate (eup_int(2, solver%C_diff%xm, solver%C_diff%ym))
    else
      if (lsolar) allocate (edir_int(solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym))
      allocate (edn_int(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
      allocate (eup_int(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    end if
    allocate (abso_int(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))

    if (lsolar) edir_int = 0
    edn_int = 0
    eup_int = 0
    abso_int = 0

    do iwvl = 1, size(repwvl_data%wvls)

      if (myid .eq. 0 .and. get_arg(.false., lverbose)) &
        & print *, 'Computing wavelengths '//toStr(iwvl)//' / '//toStr(size(repwvl_data%wvls))//&
        & ' -- '//toStr(100._ireals * real(iwvl, ireals) / real(size(repwvl_data%wvls), ireals))//' %'// &
        & ' ('//toStr(repwvl_data%wvls(iwvl))//' nm)'//' wgt= '//toStr(repwvl_data%wgts(iwvl))

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

      if (lspectral_output) then
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
      end if

      if (lboundary_flx_only) then
        if (lsolar) then
          edir_int(1, :, :) = edir_int(1, :, :) + spec_edir(lbound(spec_edir, 1), :, :) * repwvl_data%wgts(iwvl)
          edir_int(2, :, :) = edir_int(2, :, :) + spec_edir(ubound(spec_edir, 1), :, :) * repwvl_data%wgts(iwvl)
        end if
        edn_int(1, :, :) = edn_int(1, :, :) + spec_edn(lbound(spec_edn, 1), :, :) * repwvl_data%wgts(iwvl)
        edn_int(2, :, :) = edn_int(2, :, :) + spec_edn(ubound(spec_edn, 1), :, :) * repwvl_data%wgts(iwvl)
        eup_int(1, :, :) = eup_int(1, :, :) + spec_eup(lbound(spec_eup, 1), :, :) * repwvl_data%wgts(iwvl)
        eup_int(2, :, :) = eup_int(2, :, :) + spec_eup(ubound(spec_eup, 1), :, :) * repwvl_data%wgts(iwvl)
      else
        if (lsolar) edir_int(:, :, :) = edir_int(:, :, :) + spec_edir * repwvl_data%wgts(iwvl)
        edn_int(:, :, :) = edn_int(:, :, :) + spec_edn * repwvl_data%wgts(iwvl)
        eup_int(:, :, :) = eup_int(:, :, :) + spec_eup * repwvl_data%wgts(iwvl)
      end if
      abso_int(:, :, :) = abso_int(:, :, :) + spec_abso * repwvl_data%wgts(iwvl)

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

  subroutine write_output_data_local(comm, solver, out_filename, prefix, rdata, atm, &
      & edn, eup, abso, edn_int, eup_int, abso_int, ierr, edir, edir_int, lverbose)
    integer(mpiint), intent(in) :: comm
    class(t_solver), intent(in) :: solver
    character(len=*), intent(in) :: out_filename, prefix
    type(t_repwvl_data), intent(in) :: rdata
    type(t_tenstr_atm), target, intent(in) :: atm
    real(ireals), allocatable, dimension(:, :, :, :), intent(in) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl]
    real(ireals), allocatable, dimension(:, :, :), intent(in) :: edn_int, eup_int, abso_int ! [nlyr(+1), nx, ny]
    integer(mpiint), intent(out) :: ierr
    real(ireals), allocatable, dimension(:, :, :, :), intent(in), optional :: edir
    real(ireals), allocatable, dimension(:, :, :), intent(in), optional :: edir_int
    logical, intent(in), optional :: lverbose

    real(ireals), pointer, dimension(:, :, :) :: pplev
    character(len=default_str_len) :: groups(2), dimnames(4)
    real(ireals), allocatable :: xcoord(:), ycoord(:)
    integer(iintegers) :: i

    character(len=4) :: rankstr
    integer(mpiint) :: myid

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    write (rankstr, '(i0.4)') myid
    groups(1) = trim(out_filename)//'.'//rankstr

    dimnames(1) = trim(prefix)//'wvl'
    groups(2) = trim(prefix)//'wvl'
    call ncwrite(groups, rdata%wvls, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

    allocate (xcoord(solver%C_diff%xs:solver%C_diff%xe))
    allocate (ycoord(solver%C_diff%ys:solver%C_diff%ye))

    do i = solver%C_diff%xs, solver%C_diff%xe
      xcoord(i) = (i + .5_ireals) * solver%atm%dx
    end do
    dimnames(1) = trim(prefix)//'nx'
    groups(2) = trim(prefix)//'nx'
    call ncwrite(groups, xcoord, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

    do i = solver%C_diff%ys, solver%C_diff%ye
      ycoord(i) = (i + .5_ireals) * solver%atm%dy
    end do
    dimnames(1) = trim(prefix)//'ny'
    groups(2) = trim(prefix)//'ny'
    call ncwrite(groups, ycoord, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

    pplev(1:solver%C_one_atm1%zm, 1:solver%C_one_atm1%xm, 1:solver%C_one_atm1%ym) => atm%plev
    dimnames(1) = trim(prefix)//'nlev_pressure'
    dimnames(2) = trim(prefix)//'nx'
    dimnames(3) = trim(prefix)//'ny'
    groups(2) = trim(prefix)//'plev'
    call ncwrite(groups, reverse(pplev), ierr, dimnames=dimnames(1:3), deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

    dimnames(1) = trim(prefix)//'nlev'
    dimnames(2) = trim(prefix)//'nx'
    dimnames(3) = trim(prefix)//'ny'
    dimnames(4) = trim(prefix)//'wvl'
    if (allocated(edn)) then
      groups(2) = trim(prefix)//'edn'
      call ncwrite(groups, edn, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if
    groups(2) = trim(prefix)//'edn_int'
    call ncwrite(groups, edn_int, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

    if (allocated(eup)) then
      groups(2) = trim(prefix)//'eup'
      call ncwrite(groups, eup, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if
    groups(2) = trim(prefix)//'eup_int'
    call ncwrite(groups, eup_int, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

    if (present(edir)) then
      if (allocated(edir)) then
        groups(2) = trim(prefix)//'edir'
        call ncwrite(groups, edir, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      groups(2) = trim(prefix)//'edir_int'
      call ncwrite(groups, edir_int, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if

    dimnames(1) = trim(prefix)//'nlay'
    if (allocated(abso)) then
      groups(2) = trim(prefix)//'abso'
      call ncwrite(groups, abso, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if
    groups(2) = trim(prefix)//'abso_int'
    call ncwrite(groups, abso_int, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)

  end subroutine

  subroutine write_output_data_global(comm, solver, out_filename, prefix, rdata, atm, &
      & edn, eup, abso, edn_int, eup_int, abso_int, ierr, edir, edir_int, lverbose)
    integer(mpiint), intent(in) :: comm
    class(t_solver), intent(in) :: solver
    character(len=*), intent(in) :: out_filename, prefix
    type(t_repwvl_data), intent(in) :: rdata
    type(t_tenstr_atm), target, intent(in) :: atm
    real(ireals), allocatable, dimension(:, :, :, :), intent(in) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl]
    real(ireals), allocatable, dimension(:, :, :), intent(in) :: edn_int, eup_int, abso_int ! [nlyr(+1), nx, ny, Nwvl]
    integer(mpiint), intent(out) :: ierr
    real(ireals), allocatable, dimension(:, :, :, :), intent(in), optional :: edir
    real(ireals), allocatable, dimension(:, :, :), intent(in), optional :: edir_int
    logical, intent(in), optional :: lverbose

    real(ireals), pointer, dimension(:, :, :) :: pplev
    character(len=default_str_len) :: groups(2), dimnames(4)
    real(ireals), allocatable, dimension(:, :, :, :) :: gflx ! global domain [nlyr(+1), nx, ny, Nwvl]
    real(ireals), allocatable, dimension(:, :, :) :: gflx_int ! global domain [nlyr(+1), nx, ny]
    real(ireals), allocatable, dimension(:, :, :) :: pglobal ! global domain [nlyr(+1), nx, ny]

    integer(mpiint) :: myid

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    groups(1) = trim(out_filename)

    dimnames(1) = trim(prefix)//'wvl'
    groups(2) = trim(prefix)//'wvl'
    if (myid .eq. 0) then
      call ncwrite(groups, rdata%wvls, ierr, dimnames=dimnames(1:1), deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
    end if

    dimnames(1) = trim(prefix)//'nlev_pressure'
    dimnames(2) = trim(prefix)//'nx'
    dimnames(3) = trim(prefix)//'ny'
    groups(2) = trim(prefix)//'plev'
    pplev(1:solver%C_one_atm1%zm, 1:solver%C_one_atm1%xm, 1:solver%C_one_atm1%ym) => atm%plev
    call gather_all_toZero(solver%C_one_atm1, pplev, pglobal)
    if (myid .eq. 0) then
      call ncwrite(groups, reverse(pglobal), ierr, dimnames=dimnames(1:3), deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      deallocate (pglobal)
    end if

    dimnames(1) = trim(prefix)//'nlev'
    dimnames(2) = trim(prefix)//'nx'
    dimnames(3) = trim(prefix)//'ny'
    dimnames(4) = trim(prefix)//'wvl'

    ! Edn
    if (allocated(edn)) then
      if (lboundary_flx_only) then
        call local2global_bounds_only(solver%Csrfc_one, edn, gflx)
      else
        call local2global(solver%C_one_atm1, edn, gflx)
      end if
      if (myid .eq. 0) then
        groups(2) = trim(prefix)//'edn'
        call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      if (allocated(gflx)) deallocate (gflx)
    end if

    ! Edn_int
    if (allocated(edn_int)) then
      if (lboundary_flx_only) then
        call local2global_bounds_only_int(solver%Csrfc_one, edn_int, gflx_int)
      else
        call local2global_int(solver%C_one_atm1, edn_int, gflx_int)
      end if
      if (myid .eq. 0) then
        groups(2) = trim(prefix)//'edn_int'
        call ncwrite(groups, gflx_int, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      if (allocated(gflx_int)) deallocate (gflx_int)
    end if

    ! Eup
    if (allocated(eup)) then
      if (lboundary_flx_only) then
        call local2global_bounds_only(solver%Csrfc_one, eup, gflx)
      else
        call local2global(solver%C_one_atm1, eup, gflx)
      end if
      if (myid .eq. 0) then
        groups(2) = trim(prefix)//'eup'
        call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      if (allocated(gflx)) deallocate (gflx)
    end if

    ! Eup_int
    if (allocated(eup_int)) then
      if (lboundary_flx_only) then
        call local2global_bounds_only_int(solver%Csrfc_one, eup_int, gflx_int)
      else
        call local2global_int(solver%C_one_atm1, eup_int, gflx_int)
      end if
      if (myid .eq. 0) then
        groups(2) = trim(prefix)//'eup_int'
        call ncwrite(groups, gflx_int, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      if (allocated(gflx_int)) deallocate (gflx_int)
    end if

    ! Edir
    if (present(edir)) then
      if (allocated(edir)) then
        if (lboundary_flx_only) then
          call local2global_bounds_only(solver%Csrfc_one, edir, gflx)
        else
          call local2global(solver%C_one_atm1, edir, gflx)
        end if
        if (myid .eq. 0) then
          groups(2) = trim(prefix)//'edir'
          call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
        end if
        if (allocated(gflx)) deallocate (gflx)
      end if
    end if

    ! Edir_int
    if (present(edir_int)) then
      if (allocated(edir_int)) then
        if (lboundary_flx_only) then
          call local2global_bounds_only_int(solver%Csrfc_one, edir_int, gflx_int)
        else
          call local2global_int(solver%C_one_atm1, edir_int, gflx_int)
        end if
        if (myid .eq. 0) then
          groups(2) = trim(prefix)//'edir_int'
          call ncwrite(groups, gflx_int, ierr, dimnames=dimnames(1:3), deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
        end if
        if (allocated(gflx_int)) deallocate (gflx_int)
      end if
    end if

    ! Abso
    dimnames(1) = trim(prefix)//'nlay'
    if (allocated(abso)) then
      call local2global(solver%C_one_atm, abso, gflx)
      if (myid .eq. 0) then
        groups(2) = trim(prefix)//'abso'
        call ncwrite(groups, gflx, ierr, dimnames=dimnames, deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
        groups(2) = trim(prefix)//'abso_int'
        call ncwrite(groups, sum(gflx, dim=4), ierr, dimnames=dimnames(1:3), deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      if (allocated(gflx)) deallocate (gflx)
    end if

    ! Abso_int
    dimnames(1) = trim(prefix)//'nlay'
    if (allocated(abso_int)) then
      call local2global_int(solver%C_one_atm, abso_int, gflx_int)
      if (myid .eq. 0) then
        groups(2) = trim(prefix)//'abso_int'
        call ncwrite(groups, gflx_int, ierr, dimnames=dimnames(1:3), deflate_lvl=0, verbose=lverbose); call CHKERR(ierr)
      end if
      if (allocated(gflx_int)) deallocate (gflx_int)
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

    subroutine local2global_bounds_only_int(C, local, gflx)
      type(t_coord), intent(in) :: C
      real(ireals), intent(in) :: local(:, :, :)
      real(ireals), allocatable, intent(out) :: gflx(:, :, :) ! global array on rank0
      real(ireals), allocatable :: tmp(:, :, :) ! global array on rank0
      integer(iintegers) :: k
      if (myid .eq. 0) then
        allocate (gflx(size(local, dim=1), C%glob_xm, C%glob_ym))
      end if
      do k = 1, 2
        call gather_all_toZero(C, local(k:k, :, :), tmp)
        if (myid .eq. 0) gflx(k:k, :, :) = tmp
      end do
    end subroutine
    subroutine local2global_int(C, local, gflx)
      type(t_coord), intent(in) :: C
      real(ireals), intent(in) :: local(:, :, :)
      real(ireals), allocatable, intent(out) :: gflx(:, :, :) ! global array on rank0
      if (myid .eq. 0) then
        allocate (gflx(C%glob_zm, C%glob_xm, C%glob_ym))
      end if
      call gather_all_toZero(C, local, gflx)
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
      & edn_int,         &
      & eup_int,         &
      & abso_int,        &
      & ierr,            &
      & lverbose,        &
      & edir, edir_int)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: bg_atm_filename, atm_filename, out_filename
    logical, intent(in) :: lsolar
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout) :: edn, eup, abso ! [nlyr(+1), nx, ny, Nwvl]
    real(ireals), allocatable, dimension(:, :, :), intent(inout) :: edn_int, eup_int, abso_int ! [nlyr(+1), nx, ny]
    integer(mpiint), intent(out) :: ierr
    logical, intent(in), optional :: lverbose
    real(ireals), allocatable, dimension(:, :, :, :), intent(inout), optional :: edir
    real(ireals), allocatable, dimension(:, :, :), intent(inout), optional :: edir_int

    type(t_tenstr_atm) :: atm
    class(t_solver), allocatable :: solver
    type(t_repwvl_data), allocatable :: repwvl_data_solar, repwvl_data_thermal

    integer(iintegers) :: Npert_atmo, ipert
    integer(mpiint) :: myid

    integer(iintegers) :: Nx, Ny ! local domain sizes, we use the y axis to perturb the loaded atmospheres
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    type(t_mie_table), allocatable :: mie_table

    real(ireals) :: mmm_edir(3), mmm_edn(3), mmm_eup(3), mmm_abso(3)
    logical :: lglobal_output, lspectral_output, lflg
    logical :: lperturb

    ierr = 0
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    lperturb = .false.
    call get_petsc_opt('', '-perturb', lperturb, lflg, ierr); call CHKERR(ierr)

    call perturb_atmospheres(comm, lperturb, Npert_atmo, atm, ierr, lverbose, ldryrun=.true.); call CHKERR(ierr)
    if (myid .eq. 0 .and. lverbose) print *, 'Number of perturbations', Npert_atmo

    call load_atmospheres(comm, bg_atm_filename, atm_filename, atm, Nx, Ny, nxproc, nyproc, ierr, lverbose); call CHKERR(ierr)
    Ny = Npert_atmo

    nyproc(:) = Ny
    call perturb_atmospheres(comm, lperturb, Ny, atm, ierr, lverbose)

    call allocate_pprts_solver_from_commandline(solver, '2str', ierr); call CHKERR(ierr)

    call repwvl_init(comm, repwvl_data_solar, repwvl_data_thermal, ierr, lverbose=.false.); call CHKERR(ierr)
    if (myid .eq. 0 .and. lverbose) print *, 'repwvl_data_thermal Nwvl', size(repwvl_data_thermal%wvls)
    if (myid .eq. 0 .and. lverbose) print *, 'repwvl_data_solar Nwvl', size(repwvl_data_solar%wvls)
    call mie_tables_init(comm, mie_table, ierr, lverbose=.false.)
    call fu_ice_init(comm, ierr, lverbose=.false.); call CHKERR(ierr)

    call init_solver(comm, atm, Nx, Ny, nxproc, nyproc, solver, ierr)

    lspectral_output = .true.
    call get_petsc_opt('', '-spectral_output', lspectral_output, lflg, ierr); call CHKERR(ierr)

    lglobal_output = .true.
    call get_petsc_opt('', '-globalIO', lglobal_output, lflg, ierr); call CHKERR(ierr)

    lboundary_flx_only = .true.
    call get_petsc_opt('', '-boundary_flx_only', lboundary_flx_only, lflg, ierr); call CHKERR(ierr)

    if (lsolar) then
      call solve_scene(&
        & comm, &
        & solver, &
        & atm, &
        & repwvl_data_solar, &
        & mie_table, &
        & lspectral_output, &
        & edn, &
        & eup, &
        & abso, &
        & edn_int, &
        & eup_int, &
        & abso_int, &
        & ierr, &
        & edir, &
        & edir_int, &
        & lverbose=lverbose); call CHKERR(ierr)

      if (lglobal_output) then
        call write_output_data_global(comm, solver, out_filename, 'solar_', repwvl_data_solar, &
          & atm, edn, eup, abso, edn_int, eup_int, abso_int, ierr, edir, edir_int, lverbose); call CHKERR(ierr)
      else
        call write_output_data_local(comm, solver, out_filename, 'solar_', repwvl_data_solar, &
          & atm, edn, eup, abso, edn_int, eup_int, abso_int, ierr, edir, edir_int, lverbose); call CHKERR(ierr)
      end if

      if (myid .eq. 0) then
        do ipert = 1, Ny
          print *, 'edir @ srfc ( perturbation='//toStr(ipert)//')', edir_int(ubound(edir_int, 1), :, ipert)
          print *, 'edn  @ srfc ( perturbation='//toStr(ipert)//')', edn_int(ubound(edn_int, 1), :, ipert)
        end do
      end if
      if (lspectral_output) then
        call imp_min_mean_max(comm, edir, mmm_edir)
        call imp_min_mean_max(comm, edn, mmm_edn)
        call imp_min_mean_max(comm, eup, mmm_eup)
        call imp_min_mean_max(comm, abso, mmm_abso)
        if (myid .eq. 0) then
          print *, 'Min/Mean/Max edir', mmm_edir
          print *, 'Min/Mean/Max edn ', mmm_edn
          print *, 'Min/Mean/Max eup ', mmm_eup
          print *, 'Min/Mean/Max abso', mmm_abso
        end if
      end if

    else ! lthermal

      call solve_scene(&
        & comm, &
        & solver, &
        & atm, &
        & repwvl_data_thermal, &
        & mie_table, &
        & lspectral_output, &
        & edn, &
        & eup, &
        & abso, &
        & edn_int, &
        & eup_int, &
        & abso_int, &
        & ierr, &
        & lverbose=lverbose); call CHKERR(ierr)

      if (lglobal_output) then
        call write_output_data_global(comm, solver, out_filename, 'thermal_', repwvl_data_thermal, &
          & atm, edn, eup, abso, edn_int, eup_int, abso_int, ierr, lverbose=lverbose); call CHKERR(ierr)
      else
        call write_output_data_local(comm, solver, out_filename, 'thermal_', repwvl_data_thermal, &
          & atm, edn, eup, abso, edn_int, eup_int, abso_int, ierr, lverbose=lverbose); call CHKERR(ierr)
      end if

      if (myid .eq. 0) then
        do ipert = 1, Ny
          print *, 'edn @ srfc ( perturbation='//toStr(ipert)//')', edn_int(ubound(edn_int, 1), :, ipert)
          print *, 'eup @ srfc ( perturbation='//toStr(ipert)//')', eup_int(ubound(eup_int, 1), :, ipert)
        end do
      end if
      if (lspectral_output) then
        call imp_min_mean_max(comm, edn, mmm_edn)
        call imp_min_mean_max(comm, eup, mmm_eup)
        call imp_min_mean_max(comm, abso, mmm_abso)
        if (myid .eq. 0) then
          print *, 'Min/Mean/Max edn ', mmm_edn
          print *, 'Min/Mean/Max eup ', mmm_eup
          print *, 'Min/Mean/Max abso', mmm_abso
        end if
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
  real(ireals), allocatable, dimension(:, :, :) :: edir_int, edn_int, eup_int, abso_int
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
      & edn_int, eup_int, abso_int,&
      & ierr,                  &
      & lverbose=lverbose,     &
      & edir=edir, edir_int=edir_int); call CHKERR(ierr)
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
      & edn_int, eup_int, abso_int,&
      & ierr,                  &
      & lverbose); call CHKERR(ierr)
  else
    if (myid .eq. 0) print *, 'Not computing thermal part: enable with -thermal'
  end if

  call finalize_mpi(comm, .true., .true.)
end program
