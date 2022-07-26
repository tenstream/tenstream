!-------------------------------------------------------------------------
! This file is part of the TenStream solver.
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

module m_repwvl_optprop

  use m_data_parameters, only: &
    & AVOGADRO, &
    & EARTHACCEL, &
    & iintegers, &
    & ireals, &
    & K_BOLTZMANN, &
    & MOLMASSAIR, &
    & mpiint, &
    & R_DRY_AIR

  use m_helper_functions, only: &
    & CHKERR, &
    & toStr

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm
  use m_fu_ice, only: fu_ice_optprop, fu_ice_data_solar, fu_ice_data_thermal
  use m_mie_tables, only: t_mie_table, mie_optprop
  use m_rayleigh, only: rayleigh
  use m_repwvl_base, only: t_repwvl_data, repwvl_log_events
  use m_search, only: find_real_location

  implicit none

  private
  public :: &
    & check_fu_table_consistency, &
    & repwvl_optprop

  ! constant trace gas volume mixing ratios
  real(ireals), parameter :: CO = 1e-9_ireals
  real(ireals), parameter :: HNO3 = 1e-9_ireals
  real(ireals), parameter :: N2 = 0.78102_ireals

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .true.
#else
  logical, parameter :: ldebug = .true.
#endif

contains

  subroutine repwvl_optprop(repwvl_data, atm, mie_table, lsolar, k, icol, iwvl, kabs, ksca, kg, ierr)
    type(t_repwvl_data), intent(in) :: repwvl_data
    type(t_tenstr_atm), intent(in), target :: atm
    type(t_mie_table), intent(in) :: mie_table
    logical, intent(in) :: lsolar
    integer(iintegers), intent(in) :: k, icol, iwvl
    real(ireals), intent(out) :: kabs, ksca, kg
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: VMRS(repwvl_data%Ntracer)
    real(ireals) :: tabs, tsca, g
    real(ireals) :: P, dP, dtau, rayleigh_xsec, N, lwc_vmr, qext_cld, w0_cld, g_cld, iwp

    logical, parameter :: lprofile = ldebug

    ierr = 0

    !do j = lbound(kabs, 3), ubound(kabs, 3)
    !do i = lbound(kabs, 2), ubound(kabs, 2)
    !icol = i + (j - 1) * size(kabs, dim=2)
    !do k = lbound(kabs, 1), ubound(kabs, 1)

    P = (atm%plev(k, icol) + atm%plev(k + 1, icol))*.5_ireals * 1e2_ireals
    dP = (atm%plev(k, icol) - atm%plev(k + 1, icol)) * 1e2_ireals

    ! Repwvl molecular absorption cross section
    if (lprofile) then
      call PetscLogEventBegin(repwvl_log_events%repwvl_optprop_dtau, ierr); call CHKERR(ierr)
    end if
    VMRS(:) = [ &
      & atm%h2o_lay(k, icol), &
      & atm%h2o_lay(k, icol), &
      & atm%co2_lay(k, icol), &
      & atm%o3_lay(k, icol), &
      & atm%n2o_lay(k, icol), &
      & CO, &
      & atm%ch4_lay(k, icol), &
      & atm%o2_lay(k, icol), &
      & HNO3, &
      & N2]

    call repwvl_dtau(&
      & repwvl_data, &
      & iwvl, &
      & P, &
      & dP, &
      & atm%tlay(k, icol), &
      & VMRS, &
      & dtau, &
      & ierr); call CHKERR(ierr)

    tabs = dtau / atm%dz(k, icol)
    if (lprofile) then
      call PetscLogEventEnd(repwvl_log_events%repwvl_optprop_dtau, ierr); call CHKERR(ierr)
    end if
    if (ldebug) then
      if (tabs .lt. 0) call CHKERR(1_mpiint, 'kabs from repwvl negative!'//toStr(tabs))
    end if

    ! Repwvl molecular scattering cross section
    if (lprofile) then
      call PetscLogEventBegin(repwvl_log_events%repwvl_optprop_rayleigh, ierr); call CHKERR(ierr)
    end if
    call rayleigh(&
      & repwvl_data%wvls(iwvl) * 1e-3_ireals, &
      & atm%co2_lay(k, icol), &
      & rayleigh_xsec, &
      & ierr); call CHKERR(ierr)
    if (ldebug) then
      if (rayleigh_xsec .lt. 0) call CHKERR(1_mpiint, 'rayleigh xsec negative!'//toStr(rayleigh_xsec))
    end if

    N = dP * AVOGADRO / EARTHACCEL / MOLMASSAIR
    tsca = N * rayleigh_xsec * 1e-4 / atm%dz(k, icol) ! [1e-4 from cm2 to m2]
    if (ldebug) then
      if (tsca .lt. 0) call CHKERR(1_mpiint, 'rayleigh scattering coeff negative!'//toStr(tsca))
    end if
    g = 0                                             ! rayleigh has symmetric asymmetry parameter

    if (lprofile) then
      call PetscLogEventEnd(repwvl_log_events%repwvl_optprop_rayleigh, ierr); call CHKERR(ierr)
    end if

    ! Repwvl water cloud
    if (atm%lwc(k, icol) > 0) then
      if (lprofile) then
        call PetscLogEventBegin(repwvl_log_events%repwvl_optprop_mie, ierr); call CHKERR(ierr)
      end if
      call mie_optprop(&
        & mie_table, &
        & repwvl_data%wvls(iwvl) * 1e-3_ireals, &
        & atm%reliq(k, icol), &
        & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)

      lwc_vmr = atm%lwc(k, icol) * dP / (EARTHACCEL * atm%dz(k, icol)) ! have lwc in [ g / kg ], lwc_vmr in [ g / m3 ]
      qext_cld = qext_cld * 1e-3 * lwc_vmr                             ! from [km^-1 / (g / m^3)] to [1/m]

      g = (g * tsca + g_cld * qext_cld) / (tsca + qext_cld)
      tabs = tabs + qext_cld * max(0._ireals, (1._ireals - w0_cld))
      tsca = tsca + qext_cld * w0_cld
      if (lprofile) then
        call PetscLogEventEnd(repwvl_log_events%repwvl_optprop_mie, ierr); call CHKERR(ierr)
      end if
    end if

    ! Repwvl ice cloud
    if (atm%iwc(k, icol) > 0) then
      if (lprofile) then
        call PetscLogEventBegin(repwvl_log_events%repwvl_optprop_fu_ice, ierr); call CHKERR(ierr)
      end if

      call get_fu_ice_optprop()
      iwp = atm%iwc(k, icol) * dP / (EARTHACCEL * atm%dz(k, icol)) ! have iwc in [ g / kg ], iwp in [ g / m3 ]
      qext_cld = qext_cld * iwp                                    ! from [m^-1 / (g / m^3)] to [1/m]

      g = (g * tsca + g_cld * qext_cld) / (tsca + qext_cld)
      tabs = tabs + qext_cld * max(0._ireals, (1._ireals - w0_cld))
      tsca = tsca + qext_cld * w0_cld
      if (lprofile) then
        call PetscLogEventEnd(repwvl_log_events%repwvl_optprop_fu_ice, ierr); call CHKERR(ierr)
      end if
    end if

    !kabs(size(kabs, dim=1) + 1 - k, i, j) = tabs
    !ksca(size(ksca, dim=1) + 1 - k, i, j) = tsca
    !kg(size(kg, dim=1) + 1 - k, i, j) = g
    kabs = tabs
    ksca = tsca
    kg = g

    !end do
    !end do
    !end do

  contains

    subroutine get_fu_ice_optprop()
      if (lsolar) then
        if (fu_ice_data_solar%is_repwvl) then
          call fu_ice_optprop(&
            & fu_ice_data_solar, &
            & iwvl, &
            & atm%reice(k, icol), &
            & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
        else
          call fu_ice_optprop(&
            & fu_ice_data_solar, &
            & repwvl_data%wvls(iwvl) * 1e-3_ireals, &
            & atm%reice(k, icol), &
            & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
        end if
      else
        if (fu_ice_data_thermal%is_repwvl) then
          call fu_ice_optprop(&
            & fu_ice_data_thermal, &
            & iwvl, &
            & atm%reice(k, icol), &
            & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
        else
          call fu_ice_optprop(&
            & fu_ice_data_thermal, &
            & repwvl_data%wvls(iwvl) * 1e-3_ireals, &
            & atm%reice(k, icol), &
            & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
        end if
      end if
    end subroutine

  end subroutine

  subroutine check_fu_table_consistency(repwvl_data_solar, repwvl_data_thermal)
    type(t_repwvl_data) :: repwvl_data_solar, repwvl_data_thermal
    if (fu_ice_data_solar%is_repwvl) then
      if (size(fu_ice_data_solar%wvl) .ne. size(repwvl_data_solar%wvls)) then
        call CHKERR(1_mpiint, 'loaded a repwvl fu table that does not fit the solar repwvl file'// &
          & ' Nwvl_fu '//toStr(size(fu_ice_data_solar%wvl))//&
          & ' Nwvl_repwvl = '//toStr(size(repwvl_data_solar%wvls)))
      end if
    end if
    if (fu_ice_data_thermal%is_repwvl) then
      if (size(fu_ice_data_thermal%wvl) .ne. size(repwvl_data_thermal%wvls)) then
        call CHKERR(1_mpiint, 'loaded a repwvl fu table that does not fit the thermal repwvl file'// &
          & ' Nwvl_fu = '//toStr(size(fu_ice_data_thermal%wvl))//&
          & ' Nwvl_repwvl = '//toStr(size(repwvl_data_thermal%wvls)))
      end if
    end if
  end subroutine

  subroutine repwvl_dtau(    &
      & repwvl_data,         & !> data tables
      & iwvl,                & !> wavelenght index
      & P,                   & !> mean pressure of layer [Pa]
      & dP,                  & !> delta pressure of layer [Pa]
      & T,                   & !> mean Temperature of layer [K]
      & VMRS,                & !> vol mix ratios as in the repwvl data, e.g. (H2O, H2O, CO2, O3, N2O, CO, CH4, O2, HNO3, N2)
      & dtau,                & !> optical thickness of layer
      & ierr)

    type(t_repwvl_data), intent(in) :: repwvl_data
    integer(iintegers), intent(in) :: iwvl
    real(ireals), intent(in) :: P, dP, T
    real(ireals), intent(in) :: VMRS(:)
    real(ireals), intent(out) :: dtau
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: numDens, wP, wT
    real(ireals) :: wgt_P, wgt_T, xsecP0, xsecP1
    integer(iintegers) :: ip0, ip1, iT0, iT1, spec
    ierr = 0

    if (any([T, P, dP] .lt. 0)) then
      ierr = 1
      call CHKERR(ierr, "Found bad input, one of variables T,p,dP < 0:"//new_line('')// &
        & " T = "//toStr(T)//new_line('')// &
        & " p = "//toStr(P)//new_line('')// &
        & " dP= "//toStr(dP)//new_line('')// &
        & "")
    end if

    if (any(VMRS .lt. 0) .or. any(VMRS .gt. 1)) then
      ierr = 1
      call CHKERR(ierr, "Found bad input, one of VMRS < 0 or > 1:"//new_line('')// &
        & " VMRS = "//toStr(VMRS)//new_line('')// &
        & "")
    end if

    numDens = dP * AVOGADRO / MOLMASSAIR / EARTHACCEL

    wP = find_real_location(repwvl_data%p_ref, P)
    ip0 = int(floor(wP), iintegers)
    ip1 = ip0 + 1

    wT = find_real_location(repwvl_data%t_pert + repwvl_data%t_ref(ip0), T)
    iT0 = int(floor(wT), iintegers)
    iT1 = iT0 + 1

    wgt_P = wP - real(ip0, ireals)
    wgt_T = wT - real(iT0, ireals)

    dtau = 0

    ! H2O continuum is quadratic
    spec = 1
    xsecP0 = (repwvl_data%xsec(ip0, iwvl, spec, iT1) * wgt_T &
           & + repwvl_data%xsec(ip0, iwvl, spec, iT0) * (1._ireals - wgt_T)) &
           & / repwvl_data%vmrs_ref(ip0, spec)
    xsecP1 = (repwvl_data%xsec(ip1, iwvl, spec, iT1) * wgt_T &
           & + repwvl_data%xsec(ip1, iwvl, spec, iT0) * (1._ireals - wgt_T)) &
           & / repwvl_data%vmrs_ref(ip1, spec)
    dtau = dtau + (xsecP1 * wgt_P + xsecP0 * (1._ireals - wgt_P)) * numDens * (VMRS(spec) * VMRS(spec))

    do spec = 2, size(VMRS)
      xsecP0 = (repwvl_data%xsec(ip0, iwvl, spec, iT1) * wgt_T &
             & + repwvl_data%xsec(ip0, iwvl, spec, iT0) * (1._ireals - wgt_T))
      xsecP1 = (repwvl_data%xsec(ip1, iwvl, spec, iT1) * wgt_T &
             & + repwvl_data%xsec(ip1, iwvl, spec, iT0) * (1._ireals - wgt_T))
      dtau = dtau + (xsecP1 * wgt_P + xsecP0 * (1._ireals - wgt_P)) * numDens * VMRS(spec)
    end do

    if (allocated(repwvl_data%crs_o3)) &
      & call additional_tracer_from_bremen_data(repwvl_data%crs_o3(:, iwvl), VMRS(4), dtau)
    if (allocated(repwvl_data%crs_n2o)) &
      & call additional_tracer_from_bremen_data(repwvl_data%crs_n2o(:, iwvl), VMRS(5), dtau)

    if (dtau .lt. 0) then
      ierr = 2
      call CHKERR(ierr, "Resulted in negative dtau ?! dtau = "//toStr(dtau))
    end if

  contains
    subroutine additional_tracer_from_bremen_data(coeff, vmr, tau)
      real(ireals), intent(in) :: coeff(:), vmr
      real(ireals), intent(inout) :: tau

      real(ireals), parameter :: T0 = 273.15
      real(ireals) :: rho, dz, N, sigma, dtau

      rho = P / (R_DRY_AIR * T)
      dz = dP / (rho * EARTHACCEL)
      N = P / (K_BOLTZMANN * T) * 1e-4_ireals * dz
      sigma = max(0._ireals, (coeff(1) + coeff(2) * (T - T0) + coeff(3) * (T - T0)**2) * 1e-20_ireals)
      dtau = vmr * N * sigma
      tau = tau + dtau
    end subroutine
  end subroutine

end module
