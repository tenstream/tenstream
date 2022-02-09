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

module m_tenstream_options

  use m_data_parameters, only : init_mpi_data_parameters, &
    & ireals, irealLUT, iintegers, mpiint, &
    & zero, one, i0, default_str_len
  use m_optprop_parameters, only: lut_basename, stddev_atol, stddev_rtol
  use m_helper_functions, only: CHKERR, CHKWARN, get_petsc_opt

#include "petsc/finclude/petsc.h"
  use petsc

  implicit none

  logical :: &
    luse_eddington    =.True. , & ! use delta eddington coefficients for upper atmosphere, if False , we use boxmc 2-str coeffs
    lcalc_nca         =.False., & ! calculate twostream and modify absorption with NCA algorithm
    lschwarzschild    =.False., & ! use schwarzschild solver instead of twostream for thermal calculations
    lmcrts            =.False., & ! use monte carlo solver
    lskip_thermal     =.False., & ! Skip thermal calculations and just return zero for fluxes and absorption
    ltopography       =.False., & ! use raybending to include surface topography
    lLUT_mockup       =.False.

  real(ireals) :: twostr_ratio, &
    options_max_solution_err, options_max_solution_time

  integer(iintegers) :: pert_xshift, pert_yshift, mcrts_photons_per_pixel

contains
  subroutine show_options()
    print *,'------------------------------------------------------------------------------------------------------------------'
    print *,'------------------------------------------------------------------------------------------------------------------'
    print *,'Tenstream options:'
    print *,'-show_options         :: show this text                                                                           '
    print *,'-schwarzschild        :: use schwarzschild solver instead of twostream for thermal calculations                   '
    print *,'-mcrts                :: use a montecarlo solver'
    print *,'-mcrts_photons_per_px :: number of photons per pixel'
    print *,'-twostr_ratio <limit> :: when aspect ratio (dz/dx) is larger than <limit> then we use twostr_coeffs(default = 2.)'
    print *,'-calc_nca             :: calculate twostream and modify absorption with NCA algorithm (Klinger)                   '
    print *,'-skip_thermal         :: skip thermal calculations and just return zero for flux and absorption                   '
    print *,'-topography           :: use raybending to include surface topography, needs a 3D dz information                  '
    print *,'-pert_xshift <i>      :: shift optical properties in x direction by <i> pixels                                    '
    print *,'-pert_yshift <j>      :: shift optical properties in Y direction by <j> pixels                                    '
    print *,'-max_solution_err [W] :: if max error of solution is estimated below this value, skip calculation                 '
    print *,'-max_solution_time[s] :: if last update of solution is older, update irrespective of estimated error              '
    print *,'-lut_basename         :: path to LUT table files -- default is local dir                                          '
    print *,'------------------------------------------------------------------------------------------------------------------'
    print *,'------------------------------------------------------------------------------------------------------------------'
  end subroutine
  subroutine read_commandline_options(comm)
    integer(mpiint), intent(in) :: comm
    logical :: lflg=.False.
    integer(mpiint) :: ierr
    logical :: lshow_options=.False.
    logical :: ltenstr_view=.False.
    logical :: file_exists

    integer(mpiint) :: myid, numnodes
    character(len=default_str_len) :: env_lut_basename

    call init_mpi_data_parameters(comm)

    call MPI_COMM_RANK( comm, myid, ierr); call CHKERR(ierr)
    call MPI_Comm_size( comm, numnodes, ierr); call CHKERR(ierr)

    inquire(file='tenstream.options', exist=file_exists)
    if(file_exists) then
      call PetscOptionsInsertFile(comm, PETSC_NULL_OPTIONS, 'tenstream.options', PETSC_FALSE, ierr); call CHKERR(ierr)
    endif

    call get_petsc_opt(PETSC_NULL_CHARACTER,"-show_options",lshow_options,lflg,ierr) ;call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) then
      if(lshow_options) then
        if(myid.eq.0) call show_options()
        call mpi_barrier(comm, ierr)
        call CHKERR(1_mpiint, 'Exiting after show_options')
      endif
    endif

    options_max_solution_err = 5e3_ireals/real(3600*24, ireals)
    call get_petsc_opt(PETSC_NULL_CHARACTER,"-max_solution_err",&
      options_max_solution_err, lflg,ierr)  ; call CHKERR(ierr)

    options_max_solution_time = 0
    call get_petsc_opt(PETSC_NULL_CHARACTER,"-max_solution_time",&
      options_max_solution_time, lflg,ierr)  ; call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-eddington",luse_eddington,lflg,ierr); call CHKERR(ierr)

    twostr_ratio = 2._ireals
    call get_petsc_opt(PETSC_NULL_CHARACTER,"-twostr_ratio",twostr_ratio, lflg,ierr); call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER,"-pert_xshift",pert_xshift, lflg,ierr); call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) pert_xshift=0
    call get_petsc_opt(PETSC_NULL_CHARACTER,"-pert_yshift",pert_yshift, lflg,ierr); call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) pert_yshift=0

    call get_environment_variable("LUT_BASENAME", env_lut_basename, status=ierr)
    if(ierr.eq.0) lut_basename = trim(env_lut_basename)

    call get_petsc_opt(PETSC_NULL_CHARACTER,'-lut_basename', &
      lut_basename, lflg, ierr); call CHKERR(ierr)

    lLUT_mockup=.False.
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-LUT_mockup", &
      lLUT_mockup , lflg , ierr) ;call CHKERR(ierr)
    if(lLUT_mockup) then
      call CHKWARN(1_mpiint, 'Using LUT_mockup, setting the LUT constraints to zero. Your results will be wrong!')
      stddev_atol = .1_irealLUT
      stddev_rtol = .5_irealLUT
    endif

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-calc_nca", &
      lcalc_nca , lflg , ierr) ;call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-topography", &
      ltopography, lflg, ierr); call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-skip_thermal", &
      lskip_thermal, lflg, ierr); call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-schwarzschild", &
      lschwarzschild, lflg, ierr); call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-mcrts", &
      lmcrts, lflg, ierr); call CHKERR(ierr)

    mcrts_photons_per_pixel = 1000
    call get_petsc_opt(PETSC_NULL_CHARACTER, "-mcrts_photons_per_px", &
      mcrts_photons_per_pixel, lflg,ierr); call CHKERR(ierr)

    call get_petsc_opt(PETSC_NULL_CHARACTER, "-tenstr_view", &
      ltenstr_view, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-twostr_only", lflg, ierr) ; call CHKERR(ierr)
    if(lflg) call CHKERR(1_mpiint, 'Option -twostr_only is deprecated in favor of a distinct solver option, e.g. -solver 2str')

    if(myid.eq.0.and.ltenstr_view) then
      print *,'********************************************************************'
      print *,'***   nr. of Nodes:',numnodes
      print *,'***   eddington    ',luse_eddington
      print *,'***   calc_nca     ',lcalc_nca
      print *,'***   schwarzschild',lschwarzschild
      print *,'***   mcrts        ',lmcrts
      print *,'***   skip_thermal ',lskip_thermal
      print *,'***   topography   ',ltopography
      print *,'***   twostr_ratio ',twostr_ratio
      print *,'***   size_of ireal/iintegers',sizeof(one),sizeof(i0)
      print *,'***   max_solution_err       ',options_max_solution_err
      print *,'***   max_solution_time      ',options_max_solution_time
      print *,'***   lut_basename           ',trim(lut_basename)
      print *,'********************************************************************'
      print *,''
    endif
  end subroutine

end module
