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

module m_ecckd_optprop

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
  use m_ecckd_base, only: t_ecckd_data, ecckd_log_events
  use m_search, only: find_real_location

  implicit none

  private
  public :: &
    & check_fu_table_consistency, &
    & ecckd_optprop

  ! constant trace gas volume mixing ratios
  real(ireals), parameter :: CFC11 = 0e-9_ireals
  real(ireals), parameter :: CFC12 = 0e-9_ireals
  real(ireals), parameter :: N2 = 0.78102_ireals

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .true.
#else
  logical, parameter :: ldebug = .true.
#endif

! Concentration dependence of individual gases
enum, bind(c)
  enumerator :: IConcDependenceNone = 0, &
    &        IConcDependenceLinear, &
    &        IConcDependenceLUT, &
    &        IConcDependenceRelativeLinear
end enum


contains

  subroutine ecckd_optprop(ecckd_data, atm, mie_table, lsolar, k, icol, iband, kabs, ksca, kg, ierr)
    type(t_ecckd_data), intent(in) :: ecckd_data
    type(t_tenstr_atm), intent(in), target :: atm
    type(t_mie_table), intent(in) :: mie_table
    logical, intent(in) :: lsolar
    integer(iintegers), intent(in) :: k, icol, iband
    real(ireals), intent(out) :: kabs, ksca, kg
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: VMRS(ecckd_data%n_gases)
    real(ireals) :: tabs, tsca, g
    real(ireals) :: P, dP, dtau, rayleigh_xsec, N, lwc_vmr, qext_cld, w0_cld, g_cld, iwp

    logical, parameter :: lprofile = ldebug

    ierr = 0

    P = (atm%plev(k, icol) + atm%plev(k + 1, icol)) * .5_ireals * 1e2_ireals
    dP = (atm%plev(k, icol) - atm%plev(k + 1, icol)) * 1e2_ireals

    ! ecckd molecular absorption cross section
    if (lprofile) then
      call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_dtau, ierr); call CHKERR(ierr)
    end if
    VMRS(:) = [ &
      & atm%o2_lay(k, icol) + N2 + atm%n2o_lay(k, icol) + atm%ch4_lay(k, icol), & ! composite
      & atm%h2o_lay(k, icol), &
      & atm%o3_lay(k, icol), &
      & atm%co2_lay(k, icol), &
      & atm%ch4_lay(k, icol), &
      & atm%n2o_lay(k, icol), &
      & CFC11, CFC12]

    call ecckd_dtau(&
      & ecckd_data, &
      & iband, &
      & P, &
      & dP, &
      & atm%tlay(k, icol), &
      & VMRS, &
      & dtau, &
      & ierr); call CHKERR(ierr)

    tabs = dtau / atm%dz(k, icol)
    if (lprofile) then
      call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_dtau, ierr); call CHKERR(ierr)
    end if
    if (ldebug) then
      if (tabs .lt. 0) call CHKERR(1_mpiint, 'kabs from ecckd negative!'//toStr(tabs))
    end if

    ! ecckd molecular scattering cross section
    if (lprofile) then
      call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_rayleigh, ierr); call CHKERR(ierr)
    end if

    tsca = 0
    if(allocated(ecckd_data%rayleigh_molar_scattering_coeff)) then
      tsca = dP * ecckd_data%rayleigh_molar_scattering_coeff(iband) / &
        & (EARTHACCEL * MOLMASSAIR * atm%dz(k, icol))
    endif

    !N = dP * AVOGADRO / EARTHACCEL / MOLMASSAIR
    !tsca = N * rayleigh_xsec * 1e-4 / atm%dz(k, icol) ! [1e-4 from cm2 to m2]
    if (ldebug) then
      if (tsca .lt. 0) call CHKERR(1_mpiint, 'rayleigh scattering coeff negative!'//toStr(tsca))
    end if
    g = 0                                             ! rayleigh has symmetric asymmetry parameter

    if (lprofile) then
      call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_rayleigh, ierr); call CHKERR(ierr)
    end if

    ! ecckd water cloud
!    if (atm%lwc(k, icol) > 0) then
!      if (lprofile) then
!        call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_mie, ierr); call CHKERR(ierr)
!      end if
!      call mie_optprop(&
!        & mie_table, &
!        & ecckd_data%wvls(iwvl) * 1e-3_ireals, &
!        & atm%reliq(k, icol), &
!        & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
!
!      lwc_vmr = atm%lwc(k, icol) * dP / (EARTHACCEL * atm%dz(k, icol)) ! have lwc in [ g / kg ], lwc_vmr in [ g / m3 ]
!      qext_cld = qext_cld * 1e-3 * lwc_vmr                             ! from [km^-1 / (g / m^3)] to [1/m]
!
!      g = (g * tsca + g_cld * qext_cld) / (tsca + qext_cld)
!      tabs = tabs + qext_cld * max(0._ireals, (1._ireals - w0_cld))
!      tsca = tsca + qext_cld * w0_cld
!      if (lprofile) then
!        call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_mie, ierr); call CHKERR(ierr)
!      end if
!    end if
!
!    ! ecckd ice cloud
!    if (atm%iwc(k, icol) > 0) then
!      if (lprofile) then
!        call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_fu_ice, ierr); call CHKERR(ierr)
!      end if
!
!      call get_fu_ice_optprop()
!      iwp = atm%iwc(k, icol) * dP / (EARTHACCEL * atm%dz(k, icol)) ! have iwc in [ g / kg ], iwp in [ g / m3 ]
!      qext_cld = qext_cld * iwp                                    ! from [m^-1 / (g / m^3)] to [1/m]
!
!      g = (g * tsca + g_cld * qext_cld) / (tsca + qext_cld)
!      tabs = tabs + qext_cld * max(0._ireals, (1._ireals - w0_cld))
!      tsca = tsca + qext_cld * w0_cld
!      if (lprofile) then
!        call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_fu_ice, ierr); call CHKERR(ierr)
!      end if
!    end if

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
!TODO          call fu_ice_optprop(&
!TODO            & fu_ice_data_solar, &
!TODO            & ecckd_data%wvls(iwvl) * 1e-3_ireals, &
!TODO            & atm%reice(k, icol), &
!TODO            & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
      else
!TODO          call fu_ice_optprop(&
!TODO            & fu_ice_data_thermal, &
!TODO            & ecckd_data%wvls(iwvl) * 1e-3_ireals, &
!TODO            & atm%reice(k, icol), &
!TODO            & qext_cld, w0_cld, g_cld, ierr); call CHKERR(ierr)
      end if
    end subroutine

  end subroutine

  subroutine check_fu_table_consistency(ecckd_data_solar, ecckd_data_thermal)
    type(t_ecckd_data) :: ecckd_data_solar, ecckd_data_thermal
  end subroutine

  ! Solar ckd
  !	"composite h2o o3 co2 ch4 n2o" ;
  !	:composite_constituent_id = "o2 n2 n2o ch4" ;
  ! Thermal ckd
  !	:constituent_id = "composite h2o o3 co2 ch4 n2o cfc11 cfc12" ;
  ! :composite_constituent_id = "o2 n2 n2o ch4" ;

  subroutine ecckd_dtau(     &
      & ecckd_data,          & !> data tables
      & igpt,                & !> wavelenght index
      & P,                   & !> mean pressure of layer [Pa]
      & dP,                  & !> delta pressure of layer [Pa]
      & T,                   & !> mean Temperature of layer [K]
      & VMRS,                & !> vol mix ratios as in the ecckd data
      & dtau,                & !> optical thickness of layer
      & ierr)

    type(t_ecckd_data), intent(in) :: ecckd_data
    integer(iintegers), intent(in) :: igpt
    real(ireals), intent(in) :: P, dP, T
    real(ireals), intent(in) :: VMRS(:)
    real(ireals), intent(out) :: dtau
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: numDens, wP, wT, simple_multiplier
    real(ireals) :: wgt_P0, wgt_P1, wgt_T0, wgt_T1
    integer(iintegers) :: ip0, ip1, iT0, iT1
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

    wP = find_real_location(ecckd_data%pressure, P)
    ip0 = int(floor(wP), iintegers)
    ip1 = ip0 + 1

    wT = find_real_location(ecckd_data%temperature(ip0,:), T)
    iT0 = int(floor(wT), iintegers)
    iT1 = iT0 + 1

    wgt_P1 = wP - real(ip0, ireals)
    wgt_T1 = wT - real(iT0, ireals)
    wgt_P0 = (1._ireals - wgt_P1)
    wgt_T0 = (1._ireals - wgt_T1)

    dtau = 0

    simple_multiplier = dP / (MOLMASSAIR / EARTHACCEL)
    !call additional_tracer(simple_multiplier, VMRS[0], ecckd_data%composite_conc_dependence_code, ecckd_data%composite_molar_absorption_coeff, dtau)

    !call additional_tracer(simple_multiplier, VMRS(2), &
    !  & ecckd_data%h2o_conc_dependence_code, &
    !  & ecckd_data%h2o_molar_absorption_coeff, dtau, &
    !  & h2o_mole_fraction(:), &
    !  & )

    call additional_tracer( &
      & simple_multiplier, VMRS(3), &
      & ecckd_data%o3_conc_dependence_code, &
      & ecckd_data%o3_molar_absorption_coeff, &
      & dtau)
    !print *, 'dtau1', dtau

    call additional_tracer( &
      & simple_multiplier, VMRS(4), &
      & ecckd_data%co2_conc_dependence_code, &
      & ecckd_data%co2_molar_absorption_coeff, &
      & dtau)
    !print *, 'dtau2', dtau

    call additional_tracer( &
      & simple_multiplier, VMRS(5), &
      & ecckd_data%ch4_conc_dependence_code, &
      & ecckd_data%ch4_molar_absorption_coeff, &
      & dtau, &
      & ref_vmr=ecckd_data%ch4_reference_mole_fraction)
    !print *, 'dtau3', dtau

    dtau = max(0._ireals, dtau)

!    ! H2O continuum is quadratic
!    spec = 1
!    xsecP0 = (ecckd_data%xsec(ip0, iwvl, spec, iT1) * wgt_T &
!           & + ecckd_data%xsec(ip0, iwvl, spec, iT0) * (1._ireals - wgt_T)) &
!           & / ecckd_data%vmrs_ref(ip0, spec)
!    xsecP1 = (ecckd_data%xsec(ip1, iwvl, spec, iT1) * wgt_T &
!           & + ecckd_data%xsec(ip1, iwvl, spec, iT0) * (1._ireals - wgt_T)) &
!           & / ecckd_data%vmrs_ref(ip1, spec)
!    dtau = dtau + (xsecP1 * wgt_P + xsecP0 * (1._ireals - wgt_P)) * numDens * (VMRS(spec) * VMRS(spec))
!
!    do spec = 2, size(VMRS)
!      xsecP0 = (ecckd_data%xsec(ip0, iwvl, spec, iT1) * wgt_T &
!             & + ecckd_data%xsec(ip0, iwvl, spec, iT0) * (1._ireals - wgt_T))
!      xsecP1 = (ecckd_data%xsec(ip1, iwvl, spec, iT1) * wgt_T &
!             & + ecckd_data%xsec(ip1, iwvl, spec, iT0) * (1._ireals - wgt_T))
!      dtau = dtau + (xsecP1 * wgt_P + xsecP0 * (1._ireals - wgt_P)) * numDens * VMRS(spec)
!    end do
!
!    if (allocated(ecckd_data%crs_o3)) &
!      & call additional_tracer_from_bremen_data(ecckd_data%crs_o3(:, iwvl), VMRS(4), dtau)
!    if (allocated(ecckd_data%crs_n2o)) &
!      & call additional_tracer_from_bremen_data(ecckd_data%crs_n2o(:, iwvl), VMRS(5), dtau)
!
!    if (dtau .lt. 0) then
!      ierr = 2
!      call CHKERR(ierr, "Resulted in negative dtau ?! dtau = "//toStr(dtau))
!    end if
!
  contains

    subroutine additional_tracer(simple_multiplier, vmr, dep_code, molar_abs, dtau, ref_vmr)
      real(ireals), intent(in) :: simple_multiplier, vmr
      integer(iintegers), intent(in) :: dep_code
      real(ireals), intent(in) :: molar_abs(:,:,:)
      real(ireals), intent(inout) :: dtau
      real(ireals), intent(in), optional :: ref_vmr

      select case (dep_code)
      case (IConcDependenceNone)
        dtau = dtau + simple_multiplier * ( &
          &   wgt_T0 * (wgt_P0 * molar_abs(igpt, ip0, iT0) + wgt_P1 * molar_abs(igpt, ip1, iT0)) &
          & + wgt_T1 * (wgt_P0 * molar_abs(igpt, ip0, iT1) + wgt_P1 * molar_abs(igpt, ip1, iT1)) &
          & )
      case (IConcDependenceLinear)
        dtau = dtau + simple_multiplier * vmr * ( &
          &   wgt_T0 * (wgt_P0 * molar_abs(igpt, ip0, iT0) + wgt_P1 * molar_abs(igpt, ip1, iT0)) &
          & + wgt_T1 * (wgt_P0 * molar_abs(igpt, ip0, iT1) + wgt_P1 * molar_abs(igpt, ip1, iT1)) &
          & )
      case (IConcDependenceRelativeLinear)
        !print *, 'IConcDependenceRelativeLinear fac', simple_multiplier, 'vmr', vmr, 'ref', ref_vmr
        dtau = dtau + simple_multiplier * (vmr - ref_vmr) * ( &
          &   wgt_T0 * (wgt_P0 * molar_abs(igpt, ip0, iT0) + wgt_P1 * molar_abs(igpt, ip1, iT0)) &
          & + wgt_T1 * (wgt_P0 * molar_abs(igpt, ip0, iT1) + wgt_P1 * molar_abs(igpt, ip1, iT1)) &
          & )

      case default
        call CHKERR(1_mpiint, 'dep_code '//toStr(dep_code)//' not implemented?!')
      end select

    end subroutine
  end subroutine

end module
