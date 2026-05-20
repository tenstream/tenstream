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

module m_pprts_1D_solvers

  use m_data_parameters, only: ireals, iintegers, mpiint, i0, i1, zero

  use m_helper_functions, only: &
    & CHKERR, &
    & get_petsc_opt, &
    & ind_1d_to_nd, &
    & toStr

  use m_pprts_base, only: &
    & atmk, &
    & t_solver, &
    & t_state_container

  use m_schwarzschild, only: schwarzschild

  use m_twostream, only: delta_eddington_twostream, adding_delta_eddington_twostream

  use m_tenstr_disort, only: default_flx_computation

  use m_buildings, only: t_pprts_buildings

  use m_boxmc_geometry, only: PPRTS_TOP_FACE

  implicit none
  private
  public :: disort, twostream, schwarz

  logical, parameter :: ldebug = .false.

contains

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine twostream(solver, edirTOA, solution, opt_buildings)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), target :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    real(ireals), pointer, dimension(:, :, :, :) :: xv_dir, xv_diff, xv_abso
    integer(iintegers) :: i, j, k, src

    real(ireals), allocatable :: dtau(:), kext(:), w0(:), g(:), S(:), Edn(:), Eup(:)
    real(ireals) :: mu0, incSolar, fac, Ag, Bsrfc

    integer(iintegers), allocatable :: buildings_mink(:, :)
    integer(iintegers), allocatable :: buildings_face(:, :)
    integer(iintegers) :: m, bk, idx(4)

    xv_dir => null()
    xv_diff => null()
    xv_abso => null()

    associate (atm => solver%atm, &
               C_diff => solver%C_diff, &
               C_dir => solver%C_dir, &
               C_one => solver%C_one, &
               C_one_atm => solver%C_one_atm, &
               C_one_atm1 => solver%C_one_atm1)

      if (present(opt_buildings)) then
        associate (B => opt_buildings)
          ! on each pixel, find the top most face
          allocate (buildings_mink(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          allocate (buildings_face(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          buildings_mink(:, :) = C_dir%ze
          buildings_face(:, :) = -1
          do m = 1, size(opt_buildings%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) - 1 + [C_dir%zs, C_dir%xs, C_dir%ys]

            associate (d => idx(1), k => idx(2), i => idx(3), j => idx(4))
              if (d .eq. PPRTS_TOP_FACE) then
                if (k .lt. buildings_mink(i, j)) then
                  buildings_mink(i, j) = k
                  buildings_face(i, j) = m
                end if
              end if
            end associate
          end do
        end associate
      end if

      if (solution%lsolar_rad) then
        solution%edir = zero
      end if
      solution%ediff = zero

      allocate (dtau(C_one_atm%zs:C_one_atm%ze))
      allocate (kext(C_one_atm%zs:C_one_atm%ze))
      allocate (w0(C_one_atm%zs:C_one_atm%ze))
      allocate (g(C_one_atm%zs:C_one_atm%ze))

      if (solution%lsolar_rad) then
        xv_dir => solution%edir
        mu0 = solver%sun%mu
        incSolar = edirTOA
      else
        mu0 = 0
        incSolar = 0
      end if

      xv_diff => solution%ediff
      xv_abso => solution%abso

      allocate (S(C_one_atm1%zs:C_one_atm1%ze))
      allocate (Eup(C_one_atm1%zs:C_one_atm1%ze))
      allocate (Edn(C_one_atm1%zs:C_one_atm1%ze))

      do j = C_one_atm%ys, C_one_atm%ye
        do i = C_one_atm%xs, C_one_atm%xe

          kext = atm%kabs(:, i, j) + atm%ksca(:, i, j)
          dtau = atm%dz(:, i, j) * kext
          w0 = atm%ksca(:, i, j) / max(kext, epsilon(kext))
          g = atm%g(:, i, j)

          Ag = atm%albedo(i, j)

          if (allocated(atm%planck)) then
            if (allocated(atm%Bsrfc)) then
              Bsrfc = atm%Bsrfc(i, j)
            else
              Bsrfc = atm%planck(C_one_atm1%ze, i, j)
            end if
          end if

          if (present(opt_buildings)) then
            m = buildings_face(i, j)
            if (m .gt. i0) then
              Ag = opt_buildings%albedo(m)
              bk = atmk(atm, buildings_mink(i, j) - 1)
              if (allocated(atm%planck)) then
                Bsrfc = opt_buildings%planck(m)
              end if
            else
              bk = C_one_atm%ze
            end if

            S(:) = zero
            Edn(:) = zero
            Eup(:) = zero

            if (allocated(atm%planck)) then
              call delta_eddington_twostream(&
                & dtau(C_one_atm%zs:bk), &
                & w0(C_one_atm%zs:bk), &
                & g(C_one_atm%zs:bk), &
                & mu0, incSolar, Ag, &
                & S(C_one_atm%zs:bk + 1), &
                & Edn(C_one_atm%zs:bk + 1), &
                & Eup(C_one_atm%zs:bk + 1), &
                & planck=atm%planck(C_one_atm%zs:bk + 1, i, j), &
                & planck_srfc=Bsrfc)
            else

              call delta_eddington_twostream(&
                & dtau(C_one_atm%zs:bk), &
                & w0(C_one_atm%zs:bk), &
                & g(C_one_atm%zs:bk), &
                & mu0, incSolar, Ag, &
                & S(C_one_atm%zs:bk + 1), &
                & Edn(C_one_atm%zs:bk + 1), &
                & Eup(C_one_atm%zs:bk + 1))
            end if

          else ! not buildings

            if (allocated(atm%planck)) then
              call delta_eddington_twostream(dtau, w0, g, &
                                             mu0, incSolar, Ag, &
                                             S, Edn, Eup, &
                                             planck=atm%planck(:, i, j), &
                                             planck_srfc=Bsrfc)
            else
              call adding_delta_eddington_twostream(dtau, w0, g, mu0, incSolar, atm%albedo(i, j), S, Edn, Eup)
            end if
          end if

          if (solution%lsolar_rad) then
            fac = real(solver%dirtop%area_divider, ireals) / real(solver%dirtop%streams, ireals)
            do src = i0, solver%dirtop%dof - 1
              xv_dir(src, C_dir%zs + 1:C_dir%ze, i, j) = S(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_dir(src, C_dir%zs, i, j) = S(C_one_atm1%zs) * fac
            end do
          end if

          fac = real(solver%difftop%area_divider, ireals) / real(solver%difftop%streams, ireals)
          do src = 1, solver%difftop%dof
            if (solver%difftop%is_inward(src)) then ! E_dn
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = Edn(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = Edn(C_one_atm1%zs) * fac
            else ! E_up
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = Eup(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = Eup(C_one_atm1%zs) * fac
            end if
          end do

          xv_abso(i0, :, i, j) = &
            & +Edn(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & - Edn(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) &
            & - Eup(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & + Eup(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)

          if (solution%lsolar_rad) then
            xv_abso(i0, :, i, j) = xv_abso(i0, :, i, j) &
              & + S(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
              & - S(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)
          end if
          do k = C_one%zs, C_one%ze
            xv_abso(i0, k, i, j) = xv_abso(i0, k, i, j) / atm%dz(atmk(atm, k), i, j)
          end do
        end do
      end do

      xv_dir => null()
      xv_diff => null()
      xv_abso => null()

      !Twostream solver returns fluxes as [W]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is not up to date
      solution%lchanged = .false.

      deallocate (S)
      deallocate (Edn)
      deallocate (Eup)

    end associate
  end subroutine

  !> @brief wrapper for the disort solver
  subroutine disort(solver, edirTOA, solution, opt_buildings)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), target :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    real(ireals), pointer, dimension(:, :, :, :) :: xv_dir, xv_diff, xv_abso
    integer(iintegers) :: i, j, k, src, nstreams

    real, allocatable :: &
      & Bfrac(:), &
      & Blev(:), &
      & tlev(:), &
      & kext(:), &
      & dtau(:), &
      & w0(:), &
      & g(:), &
      & FLDIR(:), &
      & FLDN(:), &
      & FLUP(:), &
      & DFDT(:), &
      & UAVG(:)
    real :: mu0, Ag, Bskin
    real(ireals) :: fac
    integer(mpiint) :: ierr
    logical :: lflg

    xv_dir => null()
    xv_diff => null()
    xv_abso => null()

    if (present(opt_buildings)) call CHKERR(1_mpiint, "buildings not implemented for pprts_disort")

    associate (atm => solver%atm, &
               C_diff => solver%C_diff, &
               C_dir => solver%C_dir, &
               C_one => solver%C_one, &
               C_one_atm => solver%C_one_atm, &
               C_one_atm1 => solver%C_one_atm1)

      nstreams = 16
      call get_petsc_opt(solver%prefix, &
                         "-disort_streams", nstreams, lflg, ierr); call CHKERR(ierr)

      if (solution%lsolar_rad) then
        solution%edir = zero
      end if
      solution%ediff = zero

      allocate (dtau(C_one_atm%zs:C_one_atm%ze))
      allocate (kext(C_one_atm%zs:C_one_atm%ze))
      allocate (w0(C_one_atm%zs:C_one_atm%ze))
      allocate (g(C_one_atm%zs:C_one_atm%ze))
      allocate (Bfrac(C_one_atm%zs:C_one_atm%ze))
      allocate (tlev(C_one_atm1%zs:C_one_atm1%ze))
      Bfrac = 1
      tlev = 300

      if (solution%lsolar_rad) then
        xv_dir => solution%edir
        mu0 = real(solver%sun%mu)
      else
        mu0 = 1
        allocate (Blev(C_one_atm%zs:C_one_atm%ze))
      end if

      xv_abso => solution%abso
      xv_diff => solution%ediff

      allocate (FLDIR(C_one_atm1%zs:C_one_atm1%ze))
      allocate (FLDN(C_one_atm1%zs:C_one_atm1%ze))
      allocate (FLUP(C_one_atm1%zs:C_one_atm1%ze))
      allocate (DFDT(C_one_atm1%zs:C_one_atm1%ze))
      allocate (UAVG(C_one_atm1%zs:C_one_atm1%ze))

      do j = C_one_atm%ys, C_one_atm%ye
        do i = C_one_atm%xs, C_one_atm%xe

          kext = real(atm%kabs(:, i, j) + atm%ksca(:, i, j))
          dtau = real(atm%dz(:, i, j) * kext)
          w0 = real(atm%ksca(:, i, j) / max(kext, epsilon(kext)))
          g = real(atm%g(:, i, j))
          Ag = real(atm%albedo(i, j))
          if (.not. solution%lsolar_rad) then
            Blev = real(atm%planck(:, i, j))
            if (allocated(atm%Bsrfc)) then
              Bskin = real(atm%Bsrfc(i, j))
            else
              Bskin = real(atm%planck(ubound(atm%planck, 1), i, j))
            end if
          end if

          call default_flx_computation(       &
            & mu0,                            &
            & real(max(0._ireals, edirTOA / solver%sun%mu)),  &
            & Ag,                             &
            & 300.,                           & ! tskin (ignored because we provide planck values directly)
            & .not. solution%lsolar_rad,      & ! lthermal
            & [0., 1.],                       & ! wavenumbers (ignored because we provide planck values directly)
            & Bfrac,                          &
            & dtau,                           &
            & w0,                             &
            & g,                              &
            & tlev,                           & ! (ignored because we provide planck values directly)
            & FLDIR, FLDN, FLUP, DFDT, UAVG,  &
            & int(nstreams), lverbose=.false., &
            & Blev=Blev, Bskin=Bskin)

          if (solution%lsolar_rad) then
            fac = real(solver%dirtop%area_divider, ireals) / real(solver%dirtop%streams, ireals)
            do src = i0, solver%dirtop%dof - 1
              xv_dir(src, C_dir%zs + 1:C_dir%ze, i, j) = FLDIR(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_dir(src, C_dir%zs, i, j) = FLDIR(C_one_atm1%zs) * fac
            end do
          end if

          fac = real(solver%difftop%area_divider, ireals) / real(solver%difftop%streams, ireals)
          do src = 1, solver%difftop%dof
            if (solver%difftop%is_inward(src)) then
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = FLDN(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = FLDN(C_one_atm1%zs) * fac
            else
              xv_diff(src - 1, C_diff%zs + 1:C_diff%ze, i, j) = FLUP(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) * fac
              xv_diff(src - 1, C_diff%zs, i, j) = FLUP(C_one_atm1%zs) * fac
            end if
          end do

          xv_abso(i0, :, i, j) = &
            & +FLDN(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & - FLDN(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) &
            & - FLUP(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & + FLUP(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)

          if (solution%lsolar_rad) then
            xv_abso(i0, :, i, j) = xv_abso(i0, :, i, j) &
              & + FLDIR(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
              & - FLDIR(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)
          end if

          do k = C_one%zs, C_one%ze
            xv_abso(i0, k, i, j) = xv_abso(i0, k, i, j) / atm%dz(atmk(atm, k), i, j)
          end do

        end do
      end do

      xv_dir => null()
      xv_diff => null()
      xv_abso => null()

      !Disort returns fluxes as [W]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is not up to date
      solution%lchanged = .false.

    end associate
  end subroutine

  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine schwarz(solver, solution, opt_buildings)
    class(t_solver) :: solver
    type(t_state_container), target :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    real(ireals), pointer, dimension(:, :, :, :) :: xv_diff, xv_abso
    integer(iintegers) :: i, j, k, idof
    integer(iintegers) :: Nmu, ak
    logical :: lflg

    real(ireals), allocatable :: dtau(:), Edn(:), Eup(:)
    real(ireals) :: Bsrfc, Ag
    integer(mpiint) :: ierr

    integer(iintegers), allocatable :: buildings_mink(:, :)
    integer(iintegers), allocatable :: buildings_face(:, :)
    integer(iintegers) :: m, bk, idx(4)

    xv_diff => null()
    xv_abso => null()

    associate ( &
      atm => solver%atm, &
      C_diff => solver%C_diff, &
      C_one => solver%C_one, &
      C_one1 => solver%C_one1, &
      C_one_atm => solver%C_one_atm, &
      C_one_atm1 => solver%C_one_atm1)

      if (present(opt_buildings)) then
        associate (B => opt_buildings)
          ! on each pixel, find the top most face
          allocate (buildings_mink(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          allocate (buildings_face(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
          buildings_mink(:, :) = C_diff%ze
          buildings_face(:, :) = -1
          do m = 1, size(B%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) - 1 + [C_diff%zs, C_diff%xs, C_diff%ys]

            associate (d => idx(1), k => idx(2), i => idx(3), j => idx(4))
              if (d .eq. PPRTS_TOP_FACE) then
                if (k .lt. buildings_mink(i, j)) then
                  buildings_mink(i, j) = k
                  buildings_face(i, j) = m
                end if
              end if
            end associate
          end do
        end associate
      end if

      if (solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling schwarschild solver for solar calculation -- stopping!')
      if (.not. allocated(atm%planck)) &
        & call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

      solution%ediff = zero

      allocate (dtau(C_one_atm%zs:C_one_atm%ze))

      xv_abso => solution%abso
      xv_diff => solution%ediff

      allocate (Eup(C_one_atm1%zs:C_one_atm1%ze))
      allocate (Edn(C_one_atm1%zs:C_one_atm1%ze))

      Nmu = 2
      call get_petsc_opt(solver%prefix, &
                         "-schwarzschild_Nmu", Nmu, lflg, ierr); call CHKERR(ierr)

      if (solver%myid .eq. 0 .and. ldebug) print *, ' CALCULATING schwarzschild ::'

      do j = C_diff%ys, C_diff%ye
        do i = C_diff%xs, C_diff%xe

          dtau = atm%dz(:, i, j) * atm%kabs(:, i, j)
          Ag = atm%albedo(i, j)

          if (allocated(atm%Bsrfc)) then
            Bsrfc = atm%Bsrfc(i, j)
          else
            Bsrfc = atm%planck(C_one_atm1%ze, i, j)
          end if

          if (present(opt_buildings)) then
            m = buildings_face(i, j)
            if (m .gt. i0) then
              Ag = opt_buildings%albedo(m)
              bk = atmk(atm, buildings_mink(i, j) - 1)
              if (allocated(atm%planck)) then
                Bsrfc = opt_buildings%planck(m)
              end if
            else
              bk = C_one_atm%ze
            end if

            Edn(:) = zero
            Eup(:) = zero

            call schwarzschild(          &
              & Nmu,                     &
              & dtau(C_one_atm%zs:bk),   &
              & Ag,                      &
              & Edn(C_one_atm%zs:bk + 1),  &
              & Eup(C_one_atm%zs:bk + 1),  &
              & atm%planck(C_one_atm%zs:bk + 1, i, j), &
              & opt_srfc_emission=Bsrfc)
          else ! not buildings
            call schwarzschild(          &
              & Nmu,                     &
              & dtau,                    &
              & Ag,                      &
              & Edn, Eup,                &
              & atm%planck(:, i, j),     &
              & opt_srfc_emission=Bsrfc)
          end if

          ! icollapse needs special case for TOA flx's
          do idof = 0, solver%difftop%dof - 1
            if (solver%difftop%is_inward(i1 + idof)) then ! Edn
              xv_diff(idof, C_diff%zs, i, j) = Edn(0) / real(solver%difftop%streams, ireals)
            else ! Eup
              xv_diff(idof, C_diff%zs, i, j) = Eup(0) / real(solver%difftop%streams, ireals)
            end if
          end do

          ! rest of the atmosphere
          do k = C_diff%zs + 1, C_diff%ze
            ak = atmk(atm, k)
            do idof = 0, solver%difftop%dof - 1
              if (solver%difftop%is_inward(i1 + idof)) then ! Edn
                xv_diff(idof, k, i, j) = Edn(ak) / real(solver%difftop%streams, ireals)
              else ! Eup
                xv_diff(idof, k, i, j) = Eup(ak) / real(solver%difftop%streams, ireals)
              end if
            end do
          end do

          xv_abso(i0, :, i, j) = &
            & +Edn(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & - Edn(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze) &
            & - Eup(atmk(atm, C_one_atm1%zs):C_one_atm1%ze - 1) &
            & + Eup(atmk(atm, C_one_atm1%zs) + 1:C_one_atm1%ze)

          do k = C_one%zs, C_one%ze
            xv_abso(i0, k, i, j) = xv_abso(i0, k, i, j) / atm%dz(atmk(atm, k), i, j)
          end do

        end do
      end do

      xv_diff => null()
      xv_abso => null()

      !Schwarzschild solver returns fluxes as [W/m^2]
      solution%lWm2_dir = .true.
      solution%lWm2_diff = .true.
      ! and mark solution that it is up to date
      solution%lchanged = .false.

      deallocate (Edn)
      deallocate (Eup)

    end associate
  end subroutine

end module m_pprts_1D_solvers
