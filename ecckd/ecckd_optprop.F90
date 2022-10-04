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
    & PI_inv, &
    & R_DRY_AIR

  use m_helper_functions, only: &
    & CHKERR, &
    & toStr

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm
  use m_fu_ice, only: fu_ice_optprop, fu_ice_data_solar, fu_ice_data_thermal
  use m_mie_tables, only: t_mie_table, mie_optprop
  use m_rayleigh, only: rayleigh

  use m_ecckd_base, only: &
    & ecckd_log_events, &
    & IConcDependenceLinear, &
    & IConcDependenceLUT, &
    & IConcDependenceNone, &
    & IConcDependenceRelativeLinear, &
    & t_ecckd_atm_gas, &
    & t_ecckd_data

  use m_search, only: find_real_location

  implicit none

  private
  public :: &
    & check_fu_table_consistency, &
    & ecckd_optprop, &
    & ecckd_planck

#ifdef __RELEASE_BUILD__
  logical, parameter :: ldebug = .false.
#else
  logical, parameter :: ldebug = .true.
#endif

contains

  subroutine ecckd_optprop(ecckd_data, atm, mie_table, lsolar, k, icol, iband, kabs, ksca, kg, ierr)
    type(t_ecckd_data), intent(in) :: ecckd_data
    type(t_tenstr_atm), intent(in), target :: atm
    type(t_mie_table), intent(in) :: mie_table
    logical, intent(in) :: lsolar
    integer(iintegers), intent(in) :: k, icol, iband
    real(ireals), intent(out) :: kabs, ksca, kg
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: tabs, tsca, g
    real(ireals) :: P, dP, dtau, rayleigh_xsec, N, lwc_vmr, qext_cld, w0_cld, g_cld, iwp

    logical, parameter :: lprofile = ldebug

    ierr = 0

    P = (atm%plev(k, icol) + atm%plev(k + 1, icol))*.5_ireals * 1e2_ireals
    dP = (atm%plev(k, icol) - atm%plev(k + 1, icol)) * 1e2_ireals

    ! ecckd molecular absorption cross section
    if (lprofile) then
      call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_dtau, ierr); call CHKERR(ierr)
    end if

    call ecckd_dtau(&
      & k, icol, &
      & ecckd_data, &
      & iband, &
      & P, &
      & dP, &
      & atm%tlay(k, icol), &
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
    if (allocated(ecckd_data%rayleigh_molar_scattering_coeff)) then
      tsca = dP * ecckd_data%rayleigh_molar_scattering_coeff(iband) / &
        & (EARTHACCEL * MOLMASSAIR * atm%dz(k, icol))
    end if

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
  !        "composite h2o o3 co2 ch4 n2o" ;
  !        :composite_constituent_id = "o2 n2 n2o ch4" ;
  ! Thermal ckd
  !        :constituent_id = "composite h2o o3 co2 ch4 n2o cfc11 cfc12" ;
  ! :composite_constituent_id = "o2 n2 n2o ch4" ;

  subroutine ecckd_dtau(     &
      & atm_k,               & !> vertical index in atm
      & atm_icol,            & !> column index in atm
      & ecckd_data,          & !> data tables
      & igpt,                & !> wavelenght index
      & P,                   & !> mean pressure of layer [Pa]
      & dP,                  & !> delta pressure of layer [Pa]
      & T,                   & !> mean Temperature of layer [K]
      & dtau,                & !> optical thickness of layer
      & ierr)

    integer(iintegers), intent(in) :: atm_k, atm_icol
    type(t_ecckd_data), intent(in) :: ecckd_data
    integer(iintegers), intent(in) :: igpt
    real(ireals), intent(in) :: P, dP, T
    real(ireals), intent(out) :: dtau
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: numDens, wP, wT, simple_multiplier, taugas
    real(ireals) :: wgt_P0, wgt_P1, wgt_T0, wgt_T1
    integer(iintegers) :: ip0, ip1, iT0, iT1, igas
    ierr = 0

    if (any([T, P, dP] .lt. 0)) then
      ierr = 1
      call CHKERR(ierr, "Found bad input, one of variables T,p,dP < 0:"//new_line('')// &
        & " T = "//toStr(T)//new_line('')// &
        & " p = "//toStr(P)//new_line('')// &
        & " dP= "//toStr(dP)//new_line('')// &
        & "")
    end if

    numDens = dP * AVOGADRO / MOLMASSAIR / EARTHACCEL

    wP = find_real_location(ecckd_data%pressure, P)
    ip0 = int(floor(wP), iintegers)
    ip1 = ip0 + 1

    wT = find_real_location(ecckd_data%temperature(ip0, :), T)
    iT0 = int(floor(wT), iintegers)
    iT1 = iT0 + 1

    wgt_P1 = wP - real(ip0, ireals)
    wgt_T1 = wT - real(iT0, ireals)
    wgt_P0 = (1._ireals - wgt_P1)
    wgt_T0 = (1._ireals - wgt_T1)

    dtau = 0

    simple_multiplier = dP / (MOLMASSAIR * EARTHACCEL)

    do igas = 1, size(ecckd_data%gases)
      call additional_tracer( &
        & simple_multiplier, &
        & ecckd_data%gases(igas), &
        & taugas, ierr)
      !print *, igas, trim(ecckd_data%gases(igas)%id), ' => tau', taugas
      call CHKERR(ierr, 'dep_code '//toStr(ecckd_data%gases(igas)%conc_dependence_code)// &
        & ' not implemented for gas: '//trim(ecckd_data%gases(igas)%id))
      dtau = dtau + taugas
    end do

    dtau = max(0._ireals, dtau)

  contains

    !pure &
    subroutine additional_tracer(simple_multiplier, gas, taugas, ierr)
      real(ireals), intent(in) :: simple_multiplier
      type(t_ecckd_atm_gas), intent(in) :: gas
      real(ireals), intent(out) :: taugas
      integer(mpiint), intent(out) :: ierr

      real(ireals) :: wC, wgt_C0, wgt_C1
      integer(iintegers) :: iC0, iC1

      ierr = 0

      associate ( &
          & code => gas%conc_dependence_code, &
          & ref_vmr => gas%reference_mole_fraction, &
          & mabs3 => gas%molar_absorption_coeff3, &
          & mabs4 => gas%molar_absorption_coeff4, &
          & mfrac1 => gas%mole_fraction1)

        select case (code)
        case (IConcDependenceNone)
          if (ldebug .and. .not. associated(gas%molar_absorption_coeff3)) &
            & call CHKERR(1_mpiint, 'gas%molar_absorption_coeff3 not associated for gas '//trim(gas%id))

          taugas = simple_multiplier * ( &
            &   wgt_T0 * (wgt_P0 * mabs3(igpt, ip0, iT0) + wgt_P1 * mabs3(igpt, ip1, iT0)) &
            & + wgt_T1 * (wgt_P0 * mabs3(igpt, ip0, iT1) + wgt_P1 * mabs3(igpt, ip1, iT1)) &
            & )
        case (IConcDependenceLinear)
          if (ldebug .and. .not. associated(gas%molar_absorption_coeff3)) &
            & call CHKERR(1_mpiint, 'gas%molar_absorption_coeff3 not associated for gas '//trim(gas%id))

          if (associated(gas%vmr)) then
            taugas = simple_multiplier * gas%vmr(atm_k, atm_icol) * ( &
              &   wgt_T0 * (wgt_P0 * mabs3(igpt, ip0, iT0) + wgt_P1 * mabs3(igpt, ip1, iT0)) &
              & + wgt_T1 * (wgt_P0 * mabs3(igpt, ip0, iT1) + wgt_P1 * mabs3(igpt, ip1, iT1)) &
              & )
          else
            if (ldebug .and. gas%vmr_const .lt. 0) &
              & call CHKERR(1_mpiint, 'gas%vmr_const has a bad value'//toStr(gas%vmr_const)//' '//trim(gas%id))
            taugas = simple_multiplier * gas%vmr_const * ( &
              &   wgt_T0 * (wgt_P0 * mabs3(igpt, ip0, iT0) + wgt_P1 * mabs3(igpt, ip1, iT0)) &
              & + wgt_T1 * (wgt_P0 * mabs3(igpt, ip0, iT1) + wgt_P1 * mabs3(igpt, ip1, iT1)) &
              & )
          end if

        case (IConcDependenceRelativeLinear)
          if (ldebug .and. .not. associated(gas%molar_absorption_coeff3)) &
            & call CHKERR(1_mpiint, 'gas%molar_absorption_coeff3 not associated for gas '//trim(gas%id))
          if (ldebug .and. ref_vmr .lt. 0) call CHKERR(1_mpiint, 'bad ref_vmr value'//toStr(ref_vmr)//' for gas '//trim(gas%id))
          if (ldebug .and. .not. associated(gas%vmr)) call CHKERR(1_mpiint, 'gas%vmr not associated for gas '//trim(gas%id))

          taugas = simple_multiplier * (gas%vmr(atm_k, atm_icol) - ref_vmr) * ( &
            &   wgt_T0 * (wgt_P0 * mabs3(igpt, ip0, iT0) + wgt_P1 * mabs3(igpt, ip1, iT0)) &
            & + wgt_T1 * (wgt_P0 * mabs3(igpt, ip0, iT1) + wgt_P1 * mabs3(igpt, ip1, iT1)) &
            & )

        case (IConcDependenceLUT)
          if (ldebug .and. .not. associated(gas%molar_absorption_coeff4)) &
            & call CHKERR(1_mpiint, 'gas%molar_absorption_coeff4 not associated for gas '//trim(gas%id))
          if (ldebug .and. .not. associated(gas%mole_fraction1)) &
            & call CHKERR(1_mpiint, 'gas%mole_fraction1 not associated for gas '//trim(gas%id))
          if (ldebug .and. .not. associated(gas%vmr)) call CHKERR(1_mpiint, 'gas%vmr not associated for gas '//trim(gas%id))

          wC = find_real_location(mfrac1, gas%vmr(atm_k, atm_icol))
          iC0 = int(floor(wC), iintegers)
          iC1 = iC0 + 1
          wgt_C1 = wC - real(iC0, ireals)
          wgt_C0 = (1._ireals - wgt_C1)

          !print *,trim(gas%id)//' interpolation vmr', vmr(atm_k, atm_icol), 'in', mfrac1, 'idx:', wC

          taugas = simple_multiplier * gas%vmr(atm_k, atm_icol) * ( &
            & wgt_C0 * ( &
            &   wgt_T0 * (wgt_P0 * mabs4(igpt, ip0, iT0, iC0) + wgt_P1 * mabs4(igpt, ip1, iT0, iC0)) &
            & + wgt_T1 * (wgt_P0 * mabs4(igpt, ip0, iT1, iC0) + wgt_P1 * mabs4(igpt, ip1, iT1, iC0)) &
            & ) + &
            & wgt_C1 * ( &
            &   wgt_T0 * (wgt_P0 * mabs4(igpt, ip0, iT0, iC1) + wgt_P1 * mabs4(igpt, ip1, iT0, iC1)) &
            & + wgt_T1 * (wgt_P0 * mabs4(igpt, ip0, iT1, iC1) + wgt_P1 * mabs4(igpt, ip1, iT1, iC1)) &
            & ))

        end select
      end associate

    end subroutine
  end subroutine

  subroutine ecckd_planck(   &
      & ecckd_data,          & !> data tables
      & igpt,                & !> wavelenght index
      & T,                   & !> Temperature [K]
      & B,                   & !> emitted blackbody flux [W/m2]
      & ierr)

    type(t_ecckd_data), intent(in) :: ecckd_data
    integer(iintegers), intent(in) :: igpt
    real(ireals), intent(in) :: T
    real(ireals), intent(out) :: B
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: wT
    real(ireals) :: wgt_T0, wgt_T1
    integer(iintegers) :: iT0, iT1
    ierr = 0

    if (T .lt. 0) then
      ierr = 1
      call CHKERR(ierr, "Found bad input, for Temperature < 0:"//new_line('')// &
        & " T = "//toStr(T)//new_line('')// &
        & "")
    end if

    wT = find_real_location(ecckd_data%temperature_planck, T)
    iT0 = int(floor(wT), iintegers)
    iT1 = iT0 + 1

    wgt_T1 = wT - real(iT0, ireals)
    wgt_T0 = (1._ireals - wgt_T1)

    B = wgt_T0 * ecckd_data%planck_function(igpt, iT0) + &
      & wgt_T1 * ecckd_data%planck_function(igpt, iT1)
    B = B * PI_inv
  end subroutine

end module
