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

  use m_data_parameters, only : init_mpi_data_parameters, ireals, iintegers, mpiint, &
    zero, one, i0, default_str_len
  use m_optprop_parameters, only: lut_basename, coeff_mode, stddev_atol, stddev_rtol
  use m_helper_functions, only: CHKERR, CHKWARN

#include "petsc/finclude/petsc.h"
  use petsc

  implicit none

  logical :: ltwostr  =.False., & ! additionally calculate delta eddington twostream solution
    ltwostr_only      =.False., & ! only calculate twostream
    luse_eddington    =.True. , & ! use delta eddington coefficients for upper atmosphere            , if False , we use boxmc 2-str coeffs
    luse_hdf5_guess   =.False., & ! try loading initial guess from file
    luse_twostr_guess =.False., & ! use twostream solution as first guess
    lcalc_nca         =.False., & ! calculate twostream and modify absorption with NCA algorithm
    lschwarzschild    =.False., & ! use schwarzschild solver instead of twostream for thermal calculations
    lmcrts            =.False., & ! use monte carlo solver
    lskip_thermal     =.False., & ! Skip thermal calculations and just return zero for fluxes and absorption
    ltopography       =.False., & ! use raybending to include surface topography
    lforce_phi        =.False., & ! Force to use the phi given in options entries(overrides values given to tenstream calls
    lforce_theta      =.False., & !
    lLUT_mockup       =.False.

  real(ireals) :: twostr_ratio, &
    ident_dx,               &
    ident_dy,               &
    options_phi,            &
    options_theta,          &
    options_max_solution_err, options_max_solution_time

  integer(iintegers) :: pert_xshift, pert_yshift, mcrts_photons_per_pixel

  character(len=default_str_len) :: ident,output_prefix
  character(len=default_str_len) :: basepath

contains
  subroutine show_options()
    print *,'------------------------------------------------------------------------------------------------------------------'
    print *,'------------------------------------------------------------------------------------------------------------------'
    print *,'Tenstream options:'
    print *,'-show_options         :: show this text                                                                           '
    print *,'-ident <run_*>        :: load optical properties from hdf5 -read petsc_solver::load_optprop (default = run_test)  '
    print *,'-ident run_test       :: load optical properties from function in petsc_solver::load_test_optprop                 '
    print *,'-out                  :: output prefix (default = ts)                                                             '
    print *,'-basepath             :: output directory (default = ./)                                                          '
    print *,'-dx -dy               :: domain size in [m] (mandatory if running with -ident <run_*> )                           '
    print *,'-phi -theta           :: solar azimuth and zenith angle (default = (180,0) == south,overhead sun)                 '
    print *,'-force_phi/theta      :: force using options provided phi/theta                                                   '
    print *,'-writeall             :: dump intermediate results                                                                '
    print *,'-twostr_only          :: only calculate twostream solution -- dont bother calculating 3D Radiation                '
    print *,'-twostr               :: calculate delta eddington twostream solution                                             '
    print *,'-schwarzschild        :: use schwarzschild solver instead of twostream for thermal calculations                   '
    print *,'-mcrts                :: use a montecarlo solver'
    print *,'-mcrts_photons_per_px :: number of photons per pixel'
    print *,'-hdf5_guess           :: if run earlier with -writeall can now use dumped solutions as initial guess              '
    print *,'-twostr_guess         :: use delta eddington twostream solution as first guess                                    '
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
    logical :: lflg=.False.,lflg_ident=.False.
    integer(mpiint) :: ierr
    logical :: lshow_options=.False.
    logical :: ltenstr_view=.False.

    integer(mpiint) :: myid, numnodes
    character(len=default_str_len) :: env_lut_basename

    call init_mpi_data_parameters(comm)

    call MPI_COMM_RANK( comm, myid, ierr); call CHKERR(ierr)
    call MPI_Comm_size( comm, numnodes, ierr); call CHKERR(ierr)

    call PetscOptionsInsertFile(comm, PETSC_NULL_OPTIONS, 'tenstream.options', PETSC_FALSE, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-show_options",lshow_options,lflg,ierr) ;call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) then
      if(lshow_options) then
        if(myid.eq.0) call show_options()
        call mpi_barrier(comm, ierr)
        call CHKERR(1_mpiint, 'Exiting after show_options')
      endif
    endif

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,'-ident',ident,lflg_ident,ierr) ; call CHKERR(ierr)
    if(lflg_ident.eqv.PETSC_FALSE) ident = 'run_test'

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,'-out',output_prefix,lflg,ierr) ; call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) output_prefix = 'ts'

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,'-basepath',basepath,lflg,ierr) ; call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) basepath = './'

    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-dx",ident_dx, lflg,ierr)  ; call CHKERR(ierr)
    if( (lflg.eqv.PETSC_FALSE) .and. (lflg_ident.eqv.PETSC_TRUE) ) then
      print *,'If we run with -ident, you need to specify "-dx" commandline option e.g. -dx 70'
      call CHKERR(1_mpiint, 'option -ident '//trim(ident)//' requires also -dx option')
    endif

    options_max_solution_err = 5e3_ireals/real(3600*24, ireals)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-max_solution_err",&
      options_max_solution_err, lflg,ierr)  ; call CHKERR(ierr)

    options_max_solution_time = 0
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-max_solution_time",&
      options_max_solution_time, lflg,ierr)  ; call CHKERR(ierr)

    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-dy",ident_dy, lflg,ierr)  ; call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) ident_dy = ident_dx

    options_phi=180; options_theta=0
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi"  , options_phi, lflg,ierr); call CHKERR(ierr)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", options_theta, lflg,ierr); call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-force_phi", lforce_phi, lflg , ierr); call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-force_theta", lforce_theta, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-eddington",luse_eddington,lflg,ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-twostr", ltwostr, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-hdf5_guess", &
      luse_hdf5_guess, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-twostr_guess", &
      luse_twostr_guess, lflg, ierr); call CHKERR(ierr)
    if(luse_twostr_guess) ltwostr = .True.

    if(luse_twostr_guess.and.luse_hdf5_guess) then
      print *,'cant use twostr_guess .AND. hdf5_guess at the same time'
      call CHKERR(1_mpiint)
    endif

    twostr_ratio = 2._ireals
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-twostr_ratio",twostr_ratio, lflg,ierr); call CHKERR(ierr)

    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-pert_xshift",pert_xshift, lflg,ierr); call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) pert_xshift=0
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,"-pert_yshift",pert_yshift, lflg,ierr); call CHKERR(ierr)
    if(lflg.eqv.PETSC_FALSE) pert_yshift=0

    call get_environment_variable("LUT_BASENAME", env_lut_basename, status=ierr)
    if(ierr.eq.0) lut_basename = trim(env_lut_basename)

    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,'-lut_basename', &
      lut_basename, lflg, ierr); call CHKERR(ierr)

    lLUT_mockup=.False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-LUT_mockup", &
      lLUT_mockup , lflg , ierr) ;call CHKERR(ierr)
    if(lLUT_mockup) then
      call CHKWARN(1_mpiint, 'Using LUT_mockup, setting the LUT constraints to zero. Your results will be wrong!')
      stddev_atol = 1._ireals
      stddev_rtol = 1._ireals
    endif

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-calc_nca", &
      lcalc_nca , lflg , ierr) ;call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-twostr_only", &
      ltwostr_only, lflg, ierr); call CHKERR(ierr)
    if(ltwostr_only) then
      twostr_ratio=zero
      ltwostr=.True.
      luse_twostr_guess=.True.
    endif

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-topography", &
      ltopography, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-skip_thermal", &
      lskip_thermal, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-schwarzschild", &
      lschwarzschild, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-mcrts", &
      lmcrts, lflg, ierr); call CHKERR(ierr)

    mcrts_photons_per_pixel = 1000
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-mcrts_photons_per_px", &
      mcrts_photons_per_pixel, lflg,ierr); call CHKERR(ierr)

    coeff_mode=0 ! use LUT by default
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-coeff_mode", &
      coeff_mode, lflg, ierr); call CHKERR(ierr)

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-tenstr_view", &
      ltenstr_view, lflg, ierr); call CHKERR(ierr)
    if(myid.eq.0.and.ltenstr_view) then
      print *,'********************************************************************'
      print *,'***   nr. of Nodes:',numnodes
      print *,'***   eddington    ',luse_eddington
      print *,'***   coeff_mode   ',coeff_mode
      print *,'***   twostr_only  ',ltwostr_only
      print *,'***   twostr       ',ltwostr
      print *,'***   twostr_guess ',luse_twostr_guess
      print *,'***   calc_nca     ',lcalc_nca
      print *,'***   schwarzschild',lschwarzschild
      print *,'***   mcrts        ',lmcrts
      print *,'***   skip_thermal ',lskip_thermal
      print *,'***   topography   ',ltopography
      print *,'***   hdf5_guess   ',luse_hdf5_guess
      print *,'***   twostr_ratio ',twostr_ratio
      print *,'***   out          ',trim(output_prefix)
      print *,'***   solar azimuth',options_phi, lforce_phi
      print *,'***   solar zenith ',options_theta, lforce_theta
      print *,'***   size_of ireal/iintegers',sizeof(one),sizeof(i0)
      print *,'***   max_solution_err       ',options_max_solution_err
      print *,'***   max_solution_time      ',options_max_solution_time
      print *,'***   lut_basename           ',trim(lut_basename)
      print *,'********************************************************************'
      print *,''
    endif
  end subroutine

end module
