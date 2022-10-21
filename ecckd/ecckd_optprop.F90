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
    & irealLUT, &
    & K_BOLTZMANN, &
    & MOLMASSAIR, &
    & mpiint, &
    & PI_inv, &
    & R_DRY_AIR

  use m_helper_functions, only: &
    & CHKERR, &
    & delta_scale_optprop, &
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

  logical, parameter :: lgeometric_optics = .false.
  real(ireals) :: geometric_sw_single_scattering_albedo = 0.999999_ireals
  real(ireals) :: geometric_sw_asymmetry_factor = 0.86_ireals
  real(ireals) :: geometric_lw_single_scattering_albedo = 0.538_ireals
  real(ireals) :: geometric_lw_asymmetry_factor = 0.925_ireals

contains

  subroutine ecckd_optprop(ecckd_data, atm, lsolar, k, icol, igpt, kabs, ksca, kg, ierr)
    type(t_ecckd_data), intent(in) :: ecckd_data
    type(t_tenstr_atm), intent(in) :: atm
    logical, intent(in) :: lsolar
    integer(iintegers), intent(in) :: k, icol, igpt
    real(ireals), intent(out) :: kabs, ksca, kg
    integer(mpiint), intent(out) :: ierr

    real(ireals) :: P, dP, dtau
    real(ireals) :: qext_cld_l, w0_cld_l, g_cld_l
    real(ireals) :: qext_cld_i, w0_cld_i, g_cld_i

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
      & igpt, &
      & P, &
      & dP, &
      & atm%tlay(k, icol), &
      & dtau, &
      & ierr); call CHKERR(ierr)

    kabs = dtau / atm%dz(k, icol)
    if (lprofile) then
      call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_dtau, ierr); call CHKERR(ierr)
    end if
    if (ldebug) then
      if (kabs .lt. 0) call CHKERR(1_mpiint, 'kabs from ecckd negative!'//toStr(kabs))
    end if

    ! ecckd molecular scattering cross section
    if (lprofile) then
      call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_rayleigh, ierr); call CHKERR(ierr)
    end if

    ksca = 0
    if (allocated(ecckd_data%rayleigh_molar_scattering_coeff)) then
      ksca = dP * ecckd_data%rayleigh_molar_scattering_coeff(igpt) / &
        & (EARTHACCEL * MOLMASSAIR * atm%dz(k, icol))
    end if

    if (ldebug) then
      if (ksca .lt. 0) call CHKERR(1_mpiint, 'rayleigh scattering coeff negative!'//toStr(ksca))
    end if
    kg = 0                                             ! rayleigh has symmetric asymmetry parameter

    if (lprofile) then
      call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_rayleigh, ierr); call CHKERR(ierr)
    end if

    if (atm%lwc(k, icol) > 0) then
      if (lprofile) then
        call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_mie, ierr); call CHKERR(ierr)
      end if

      call get_liq_cld_optprop(ecckd_data, lsolar, igpt, &
        & atm%lwc(k, icol), atm%reliq(k, icol), atm%dz(k, icol), dP, &
        & qext_cld_l, w0_cld_l, g_cld_l, ierr); call CHKERR(ierr)
      !print *, 'liq cld optprop single', qext_cld_l, w0_cld_l, g_cld_l

      kg = (kg * ksca + g_cld_l * qext_cld_l * w0_cld_l) / (ksca + qext_cld_l * w0_cld_l)
      kabs = kabs + qext_cld_l * max(0._ireals, (1._ireals - w0_cld_l))
      ksca = ksca + qext_cld_l * w0_cld_l

      if (lprofile) then
        call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_mie, ierr); call CHKERR(ierr)
      end if
    end if

    ! ecckd ice cloud
    if (atm%iwc(k, icol) > 0) then
      if (lprofile) then
        call PetscLogEventBegin(ecckd_log_events%ecckd_optprop_fu_ice, ierr); call CHKERR(ierr)
      end if

      call get_ice_cld_optprop(ecckd_data, lsolar, igpt, &
        & atm%iwc(k, icol), atm%reice(k, icol), atm%dz(k, icol), dP, &
        & qext_cld_i, w0_cld_i, g_cld_i, ierr); call CHKERR(ierr)
      !print *, 'ice cld optprop single', qext_cld_i, w0_cld_i, g_cld_i

      kg = (kg * ksca + g_cld_i * qext_cld_i * w0_cld_i) / (ksca + qext_cld_i * w0_cld_i)
      kabs = kabs + qext_cld_i * max(0._ireals, (1._ireals - w0_cld_i))
      ksca = ksca + qext_cld_i * w0_cld_i

      if (lprofile) then
        call PetscLogEventEnd(ecckd_log_events%ecckd_optprop_fu_ice, ierr); call CHKERR(ierr)
      end if
    end if

  end subroutine

  subroutine get_liq_cld_optprop(ecckd_data, lsolar, igpt, lwc, reliq, dz, dP, qext_cld_l, w0_cld_l, g_cld_l, ierr)
    type(t_ecckd_data), intent(in) :: ecckd_data
    logical, intent(in) :: lsolar
    integer(iintegers), intent(in) :: igpt
    real(ireals), intent(in) :: lwc, reliq, dz, dP
    real(ireals), intent(out) :: qext_cld_l, w0_cld_l, g_cld_l
    integer(mpiint), intent(out) :: ierr

    real(ireals), parameter :: DensityLiquidWater = 1000.0_ireals ! kg m-3

    real(ireals) :: lwp
    integer(iintegers) :: iR0, iR1
    real(ireals) :: wR, wR0, wR1

    ierr = 0

    if (lgeometric_optics) then

      if (lsolar) then
        qext_cld_l = (3.0_ireals / (2.0_ireals * DensityLiquidWater) * 1e3_ireals / EARTHACCEL) * &
          & lwc * dP / (dz * reliq)
        w0_cld_l = geometric_sw_single_scattering_albedo
        g_cld_l = geometric_sw_asymmetry_factor
      else
        qext_cld_l = lwc * dP / (EARTHACCEL * dz) * 137.22_ireals
        w0_cld_l = geometric_lw_single_scattering_albedo
        g_cld_l = geometric_lw_asymmetry_factor
      end if

      call delta_scale_optprop(qext_cld_l, w0_cld_l, g_cld_l, g_cld_l**2)

    else

      qext_cld_l = 0
      w0_cld_l = 0
      g_cld_l = 0

      wR = find_real_location(ecckd_data%mie_table%reff, real(reliq, irealLUT))

      iR0 = int(floor(wR), iintegers)
      iR1 = iR0 + 1

      wR1 = wR - real(iR0, ireals)
      wR0 = (1._ireals - wR1)

      if (wR1 .lt. sqrt(epsilon(wR0))) then
        qext_cld_l = ecckd_data%mie_table%qext(iR0, igpt)
        w0_cld_l = ecckd_data%mie_table%w0(iR0, igpt)
        g_cld_l = ecckd_data%mie_table%g(iR0, igpt)
      else
        qext_cld_l = ecckd_data%mie_table%qext(iR0, igpt) * wR0 + ecckd_data%mie_table%qext(iR1, igpt) * wR1
        w0_cld_l = ecckd_data%mie_table%w0(iR0, igpt) * wR0 + ecckd_data%mie_table%w0(iR1, igpt) * wR1
        g_cld_l = ecckd_data%mie_table%g(iR0, igpt) * wR0 + ecckd_data%mie_table%g(iR1, igpt) * wR1
      end if

      lwp = lwc * dP / (EARTHACCEL * dz) ! have lwc in [ g / kg ] -> lwc_vmr in [ g / m3 ]
      qext_cld_l = qext_cld_l * 1e-3_ireals * lwp ! from [km^-1 / (g / m^3)] to [1/m]

    end if
  end subroutine

  subroutine get_ice_cld_optprop(ecckd_data, lsolar, igpt, iwc, reice, dz, dP, qext_cld_i, w0_cld_i, g_cld_i, ierr)
    type(t_ecckd_data), intent(in) :: ecckd_data
    logical, intent(in) :: lsolar
    integer(iintegers), intent(in) :: igpt
    real(ireals), intent(in) :: iwc, reice, dz, dP
    real(ireals), intent(out) :: qext_cld_i, w0_cld_i, g_cld_i
    integer(mpiint), intent(out) :: ierr

    real(ireals), parameter :: DensitySolidIce = 916.7_ireals  ! kg m-3

    real(ireals) :: iwp
    integer(iintegers) :: iwvnr
    real(ireals) :: qext, w0, g, wgt, wvl_lo, wvl_hi, wvl

    ierr = 0

    if (lgeometric_optics) then

      qext_cld_i = (3.0_ireals / (2.0_ireals * DensitySolidIce) * 1e3_ireals / EARTHACCEL) * &
        & iwc * dP / (dz * reice)

      if (lsolar) then
        w0_cld_i = geometric_sw_single_scattering_albedo
        g_cld_i = geometric_sw_asymmetry_factor
      else
        w0_cld_i = geometric_lw_single_scattering_albedo
        g_cld_i = geometric_lw_asymmetry_factor
      end if

    else

      qext_cld_i = 0
      w0_cld_i = 0
      g_cld_i = 0

      do iwvnr = 1, size(ecckd_data%gpoint_fraction, dim=1)
        wgt = ecckd_data%gpoint_fraction(iwvnr, igpt)
        if (wgt .gt. 0) then
          wvl_lo = 1e7_ireals / ecckd_data%wavenumber2(iwvnr)
          wvl_hi = 1e7_ireals / ecckd_data%wavenumber1(iwvnr)
          wvl = (wvl_lo + wvl_hi)*.5

          if (lsolar) then
            call fu_ice_optprop(&
              & fu_ice_data_solar, &
              & wvl * 1e-3_ireals, &
              & reice, &
              & qext, w0, g, ierr); call CHKERR(ierr)
          else
            call fu_ice_optprop(&
              & fu_ice_data_thermal, &
              & wvl * 1e-3_ireals, &
              & reice, &
              & qext, w0, g, ierr); call CHKERR(ierr)
          end if

          qext_cld_i = qext_cld_i + wgt * qext
          w0_cld_i = w0_cld_i + wgt * w0
          g_cld_i = g_cld_i + wgt * g
        end if
      end do

      iwp = iwc * dP / (EARTHACCEL * dz) ! have iwc in [ g / kg ] -> iwp in [ g / m3 ]
      qext_cld_i = qext_cld_i * iwp      ! from [m^-1 / (g / m^3)] to [1/m]
    end if

    !call delta_scale_optprop(qext_cld_i, w0_cld_i, g_cld_i, g_cld_i**2)
  end subroutine

  subroutine check_fu_table_consistency()
    !type(t_ecckd_data) :: ecckd_data_solar, ecckd_data_thermal
    if (fu_ice_data_solar%is_repwvl) then
      call CHKERR(1_mpiint, 'solar fu table is repwvl but this is ecckd')
    end if
    if (fu_ice_data_thermal%is_repwvl) then
      call CHKERR(1_mpiint, 'thermal fu table is repwvl but this is ecckd')
    end if
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

    real(ireals) :: wP, wT, simple_multiplier, taugas
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

    !wP = find_real_location(ecckd_data%pressure, P)
    wP = find_real_location(ecckd_data%log_pressure, log(P))

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
      taugas = 0

      associate ( &
          & code => gas%conc_dependence_code, &
          & ref_vmr => gas%reference_mole_fraction, &
          & mabs3 => gas%molar_absorption_coeff3, &
          & mabs4 => gas%molar_absorption_coeff4, &
          & mfrac1 => gas%mole_fraction1, &
          & log_mfrac1 => gas%log_mole_fraction1)

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
          if (ldebug .and. .not. allocated(gas%log_mole_fraction1)) &
            & call CHKERR(1_mpiint, 'gas%log_mole_fraction1 not allocated for gas '//trim(gas%id))
          if (ldebug .and. .not. associated(gas%vmr)) call CHKERR(1_mpiint, 'gas%vmr not associated for gas '//trim(gas%id))

          !wC = find_real_location(mfrac1, gas%vmr(atm_k, atm_icol))
          wC = find_real_location(log_mfrac1, log(gas%vmr(atm_k, atm_icol)))
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
