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

module m_rayleigh
  use m_data_parameters, only: &
    & iintegers, &
    & ireals, &
    & ireal_dp, &
    & mpiint, &
    & pi

  implicit none

  private
  public :: rayleigh

  interface rayleigh
    module procedure rayleigh_bodhaine
    module procedure rayleigh_bodhaine29
  end interface
contains

  !> @brief Calculates the Rayleigh scattering cross section with co2 correction
  !> @details according to the formula by Bodhaine et al,
  !> \n `On Rayleigh optical depth calculations', J. Atm. Ocean Technol., 16, 1854-1861, 1999.
  elemental subroutine rayleigh_bodhaine(lambda, mixing_ratio_co2, ksca, ierr)
    real(ireals), intent(in) :: lambda ! wvl in mu
    real(ireals), intent(in) :: mixing_ratio_co2 ! CO2 mixing ratio
    real(ireals), intent(out) :: ksca ! rayleigh cross section
    integer(mpiint), intent(out) :: ierr

    real(ireals), parameter :: N_s = 2.546899e+19_ireals
    real(ireals), parameter :: ray_const = (24._ireals * pi**3 / N_s) / N_s

    real(ireals) :: l2, lm2, n2
    real(ireals) :: co2, lambda_cm
    real(ireals) :: F_N_2, F_O_2, F_air
    real(ireals) :: ref_ratio, n_300, n ! n = refractive index of air at desired co2 concentration

    ierr = 0
    co2 = mixing_ratio_co2 * 1.e-4 ! This is co2 in parts per volume by percent, top of page, 1858 of Bodhaine et al paper.
    lambda_cm = lambda * 1e-4_ireals

    l2 = lambda**2
    lm2 = 1._ireals / l2

    n_300 = (8060.51_ireals + (2480990._ireals / (132.274_ireals - lm2)) + (17455.7_ireals / (39.32957_ireals - lm2))) &
           & * 1e-08_ireals ! Eq. (18)
    n = (1._ireals + 0.54_ireals * (mixing_ratio_co2 * 1e-6_ireals - 0.0003_ireals)) * n_300 + 1._ireals ! Eq. (19)
    n2 = n**2
    ref_ratio = (n2 - 1)**2 / (n2 + 2)**2; 
    F_N_2 = 1.034_ireals + 3.17e-4_ireals / l2 ! Eq. (5)
    F_O_2 = 1.096_ireals + 1.385e-3_ireals / l2 + 1.448e-4_ireals / l2 / l2 ! Eq. (6)
    F_air = (78.084_ireals * F_N_2 + 20.946_ireals * F_O_2 + 0.934_ireals * 1.0_ireals + co2 * 1.15_ireals) / &
          & (78.084_ireals + 20.946_ireals + 0.934_ireals + co2) ! Eq. (23)

    ksca = (ray_const / lambda_cm**4) * ref_ratio * F_air
  end subroutine

  !> @brief Calculates the Rayleigh scattering cross section
  !> @details according to the formula by Bodhaine et al, Eq. 29
  !> \n `On Rayleigh optical depth calculations', J. Atm. Ocean Technol., 16, 1854-1861, 1999.
  elemental subroutine rayleigh_bodhaine29(lambda, ksca, ierr)
    real(ireals), intent(in) :: lambda ! wvl in mu
    real(ireals), intent(out) :: ksca ! rayleigh cross section
    integer(mpiint), intent(out) :: ierr
    real(ireals) :: l2
    ierr = 0
    if (lambda .lt. .125_ireals) ierr = 1
    l2 = lambda**2
    ksca = 1.e-28_ireals * (1.0455996_ireals - 341.29061_ireals / l2 - 0.90230850_ireals * l2) / &
         & (1._ireals + 0.0027059889_ireals / l2 - 85.968563_ireals * l2)
  end subroutine

end module
