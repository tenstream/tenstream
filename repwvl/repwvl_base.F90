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

!> \page Routines to call tenstream with optical properties from a representative wavelength approach

module m_repwvl_base

  use m_data_parameters, only: &
    & default_str_len, &
    & iintegers, &
    & ireals, &
    & mpiint

  use m_helper_functions, only: &
    & CHKERR, &
    & get_arg

  use m_netcdfIO, only: ncload

  implicit none

  private
  public :: t_repwvl_data, repwvl_init, calc_dtau

  type t_repwvl_data
    ! dims xsec(xsec_nbooks, xsec_npages, xsec_nrows, xsec_ncols)
    !           temp         species      wvl         pres
    real(ireals), allocatable :: xsec(:,:,:,:)
    real(ireals), allocatable :: wvls(:)
    real(ireals), allocatable :: wgts(:)
    real(ireals), allocatable :: p_grid(:)
    real(ireals), allocatable :: t_ref(:)
    real(ireals), allocatable :: t_pert(:)
    real(ireals), allocatable :: vmrs_ref(:,:)
  end type

  contains

    subroutine load_data(fname, repwvl_data, ierr, lverbose)
      character(len=*), intent(in) :: fname
      type(t_repwvl_data), intent(inout) :: repwvl_data
      integer(mpiint), intent(out) :: ierr
      logical, intent(in), optional :: lverbose

      character(len=default_str_len) :: groups(2)

      groups(1) = trim(fname)
      groups(2) = 'xsec'; call ncload(groups, repwvl_data%xsec, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'ChosenWvls'; call ncload(groups, repwvl_data%wvls, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'ChosenWeights'; call ncload(groups, repwvl_data%wgts, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'p_grid'; call ncload(groups, repwvl_data%p_grid, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 't_ref'; call ncload(groups, repwvl_data%t_ref, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 't_pert'; call ncload(groups, repwvl_data%t_pert, ierr, lverbose); call CHKERR(ierr)
      groups(2) = 'vmrs_ref'; call ncload(groups, repwvl_data%vmrs_ref, ierr, lverbose); call CHKERR(ierr)
    end subroutine

    subroutine repwvl_init(repwvl_data_thermal, ierr, fname_repwvl_thermal, fname_repwvl_solar)
      type(t_repwvl_data), intent(inout) :: repwvl_data_thermal
      integer(mpiint), intent(out) :: ierr
      character(len=*), intent(in), optional :: fname_repwvl_thermal, fname_repwvl_solar
      character(len=default_str_len) :: fname_thermal, fname_solar

      ierr = 0

      fname_thermal = get_arg('repwvl_thermal.lut', fname_repwvl_thermal)
      fname_solar = get_arg('repwvl_solar.lut', fname_repwvl_solar)

      print *,'Reading representative wavelength data from '//trim(fname_thermal)
      print *,'Reading representative wavelength data from '//trim(fname_solar)

      call load_data(fname_thermal, repwvl_data_thermal, ierr, lverbose=.True.); call CHKERR(ierr)
    end subroutine

    subroutine calc_dtau(&
        & repwvl_data_thermal, &
        & iwvl, &
        & temp, &
        & pres, &
        & VMRS, &
        & dtau  )
      type(t_repwvl_data), intent(in) :: repwvl_data_thermal
      integer(iintegers) :: iwvl
      real(ireals), intent(in) :: temp
      real(ireals), intent(in) :: pres
      real(ireals), intent(in) :: VMRS(:)
      real(ireals), intent(out) :: dtau

      print *, iwvl, temp, pres, VMRS
      dtau = -1 + repwvl_data_thermal%wvls(1)

    end subroutine
end module

!module m_pprts_rrtmg
!  use, intrinsic :: iso_c_binding
!
!#include "petsc/finclude/petsc.h"
!  use petsc
!
!  use mpi, only : mpi_comm_rank
!  use m_tenstr_parkind_sw, only: im => kind_im, rb => kind_rb
!  use m_tenstream_options, only: read_commandline_options
!  use m_data_parameters, only : init_mpi_data_parameters, &
!      iintegers, ireals, zero, one, i0, i1, i2, i9,         &
!      mpiint, pi, default_str_len
!
!  use m_pprts_base, only : t_solver, compute_gradient, destroy_pprts, atmk
!
!  use m_pprts, only : init_pprts, set_angles, set_optical_properties, solve_pprts,&
!      pprts_get_result, pprts_get_result_toZero
!
!  use m_buildings, only: &
!    & t_pprts_buildings, &
!    & PPRTS_BOT_FACE
!
!  use m_xdmf_export, only: &
!    & xdmf_pprts_buildings, &
!    & xdmf_pprts_srfc_flux
!
!  use m_adaptive_spectral_integration, only: need_new_solution
!
!  use m_helper_functions, only : &
!    & CHKERR, &
!    & approx, &
!    & cross_3d, &
!    & deg2rad, &
!    & get_arg, &
!    & gradient, &
!    & imp_allreduce_max, &
!    & imp_allreduce_mean, &
!    & imp_allreduce_min, &
!    & imp_bcast, &
!    & ind_1d_to_nd, &
!    & meanvec, &
!    & mpi_logical_all_same, &
!    & read_ascii_file_2d, &
!    & reverse, &
!    & spherical_2_cartesian, &
!    & toStr
!
!  use m_search, only: find_real_location
!  use m_petsc_helpers, only: dmda_convolve_ediff_srfc, &
!    getvecpointer, restorevecpointer, f90vectopetsc
!
!  use m_netcdfIO, only : ncwrite
!
!  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, plkint, print_tenstr_atm, vert_integral_coeff
!
!  use m_optprop_rrtmg, only: optprop_rrtm_lw, optprop_rrtm_sw, get_spectral_bands
!
!  use m_tenstr_disort, only: default_flx_computation
!  use m_tenstr_rrtmg_base, only: t_rrtmg_log_events, setup_log_events
!  use m_buildings, only: t_pprts_buildings, clone_buildings, destroy_buildings
!
!  use m_pprts_external_solvers, only: destroy_rayli_info
!
!  use m_tenstr_rrtmg_lw_init, only: rrtmg_lw_ini
!  use m_tenstr_rrtmg_sw_init, only: rrtmg_sw_ini
!
!  implicit none
!
!  private
!  public :: pprts_rrtmg, destroy_pprts_rrtmg
!
!!  logical,parameter :: ldebug=.True.
!  logical,parameter :: ldebug=.False.
!
!  type(t_rrtmg_log_events) :: log_events
!contains
!
!  subroutine init_pprts_rrtmg(comm, solver, dx, dy, dz, &
!                  sundir, &
!                  xm, ym, zm, &
!                  nxproc, nyproc, &
!                  pprts_icollapse)
!
!    integer(mpiint), intent(in) :: comm
!
!    real(ireals), intent(in)      :: dx, dy, sundir(3)
!    real(ireals), intent(in)      :: dz(:,:) ! bot to top, e.g. from atm%dz
!    integer(iintegers),intent(in) :: xm, ym, zm
!    class(t_solver),intent(inout) :: solver
!
!    ! arrays containing xm and ym for all nodes :: dim[x-ranks, y-ranks]
!    integer(iintegers),intent(in), optional :: nxproc(:), nyproc(:), pprts_icollapse
!
!    integer(iintegers) :: i, j, icol
!
!    ! vertical thickness in [m]
!    real(ireals),allocatable :: dz_t2b(:,:,:) ! dz (t2b := top 2 bottom)
!
!    call setup_log_events(log_events, 'pprts_')
!
!    if(present(nxproc) .neqv. present(nyproc)) then
!      print *,'Wrong call to init_tenstream_rrtm_lw --    &
!            & in order to work, we need both arrays for &
!            & the domain decomposition, call with nxproc AND nyproc'
!      call CHKERR(1_mpiint, 'init_tenstream_rrtm_lw -- missing arguments nxproc or nyproc')
!    endif
!
!    allocate(dz_t2b(zm, xm, ym))
!    do j=i1,ym
!      do i=i1,xm
!        icol =  i+(j-1)*xm
!        dz_t2b(:,i,j) = reverse(dz(:,icol))
!      enddo
!    enddo
!
!    if(present(nxproc) .and. present(nyproc)) then
!      call init_pprts(comm, zm, xm, ym, dx, dy, sundir, solver, nxproc=nxproc, nyproc=nyproc, dz3d=dz_t2b, &
!        collapseindex=pprts_icollapse)
!    else ! we let petsc decide where to put stuff
!      call init_pprts(comm, zm, xm, ym, dx, dy, sundir, solver, dz3d=dz_t2b, &
!        collapseindex=pprts_icollapse)
!    endif
!
!    call rrtmg_lw_ini(1006._rb)
!    call rrtmg_sw_ini(1006._rb)
!  end subroutine
!
!  subroutine smooth_surface_fluxes(solver, edn, eup)
!    class(t_solver), intent(inout)  :: solver                       ! solver type (e.g. t_solver_8_10)
!    real(ireals),allocatable, dimension(:,:,:), intent(inout) :: edn, eup  ! [nlyr+1, local_nx, local_ny ]
!    integer(iintegers) :: i, kernel_width, Niter
!    real(ireals) :: radius, mflx_up, mflx_dn
!    logical :: lflg
!    integer(mpiint) :: myid, ierr
!
!    call PetscLogEventBegin(log_events%smooth_surface_fluxes, ierr); call CHKERR(ierr)
!    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
!
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!      "-pprts_smooth_srfc_flx", radius , lflg , ierr) ;call CHKERR(ierr)
!
!    if (lflg) then
!
!      if(radius.lt.zero) then
!        call imp_allreduce_mean(solver%comm,  edn(ubound(edn,1),:,:), mflx_dn); edn(ubound(edn,1),:,:) = mflx_dn
!        call imp_allreduce_mean(solver%comm,  eup(ubound(eup,1),:,:), mflx_up); eup(ubound(eup,1),:,:) = mflx_up
!        if(ldebug.and.myid.eq.0) &
!          print *,'Smoothing diffuse srfc fluxes over the entire domain '// &
!                  'mean downward flx', mflx_dn, 'up', mflx_up
!
!      else
!
!        call find_iter_and_kernelwidth(Niter, kernel_width)
!        do i = 1, Niter
!        call dmda_convolve_ediff_srfc(solver%C_diff%da, kernel_width, edn(ubound(edn,1):ubound(edn,1),:,:))
!        call dmda_convolve_ediff_srfc(solver%C_diff%da, kernel_width, eup(ubound(eup,1):ubound(eup,1),:,:))
!        enddo
!      endif
!    endif
!    call PetscLogEventEnd(log_events%smooth_surface_fluxes, ierr); call CHKERR(ierr)
!    contains
!      subroutine find_iter_and_kernelwidth(Niter, kernel_width)
!        integer(iintegers), intent(out) :: kernel_width, Niter
!        integer(iintegers), parameter :: test_Ni=500
!        integer(iintegers), parameter :: Ni_limit=250
!        integer(iintegers) :: i, min_iter
!        real(ireals) :: test_k(test_Ni), residuals(test_Ni)
!        real(ireals) :: radius_in_pixel
!
!        radius_in_pixel = radius / ((solver%atm%dx+solver%atm%dy)/2)
!        ! Try a couple of number of iterations and determine optimal kernel_width
!        min_iter = 1
!        do i = 1, test_Ni
!          test_k(i) = (sqrt( (12*radius_in_pixel**2 + real(i, ireals)) / real(i, ireals) +1 ) -1 )/2
!          if(nint(test_k(i)).ge.min(solver%C_diff%xm,solver%C_diff%ym)) min_iter = i+1
!        enddo
!        if(min_iter.gt.Ni_limit) &
!          call CHKERR(int(min_iter, mpiint), 'the smoothing iteration count would be larger than '// &
!            toStr(Ni_limit)//'... this could become expensive, you can set the value higher but I'// &
!            'suspect that you are trying something weird')
!        ! We want it
!        ! as close as possible to an integer
!        ! and as big as possible
!        ! but not bigger than the local domain size
!        ! or a certain max value
!        residuals = abs(test_k - nint(test_k))
!
!        Niter = min_iter -1 + minloc(residuals(min_iter:test_Ni),dim=1)
!        kernel_width = nint(test_k(Niter))
!
!        call imp_allreduce_max(solver%comm, Niter, i); Niter = i
!        call imp_allreduce_min(solver%comm, kernel_width, i); kernel_width = i
!
!        if(kernel_width.eq.0) Niter=0
!
!        if(ldebug.and.myid.eq.0) then
!          do i = 1, test_Ni
!          print *, 'iter', i, 'test_k', test_k(i), residuals(i)
!          enddo
!          print *,'Smoothing diffuse srfc fluxes with radius', solver%atm%dx, radius, radius_in_pixel, &
!            'Niter', Niter, 'kwidth', kernel_width
!        endif
!
!      end subroutine
!  end subroutine
!
!  subroutine slope_correction_fluxes(solver, edir)
!    class(t_solver), intent(inout)  :: solver                         ! solver type (e.g. t_solver_8_10)
!    real(ireals),allocatable, dimension(:,:,:), intent(inout) :: edir ! [nlyr+1, local_nx, local_ny ]
!    logical :: lslope_correction, latm_correction, lflg
!    integer(mpiint) :: myid, ierr
!
!    real(ireals),pointer :: grad   (:,:,:,:) =>null()
!    real(ireals),pointer :: grad_1d(:)       =>null()
!    real(ireals) :: fac, n(3)
!    integer(iintegers) :: i,j,k
!
!    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
!
!    lslope_correction = .False.
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!      "-pprts_slope_correction", lslope_correction, lflg, ierr) ;call CHKERR(ierr)
!
!    latm_correction = .False.
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!      "-pprts_atm_correction", latm_correction, lflg, ierr) ;call CHKERR(ierr)
!
!    associate(&
!        atm     => solver%atm,  &
!        sun     => solver%sun,  &
!        C_dir   => solver%C_dir,&
!        C_two1  => solver%C_two1)
!
!      call getVecPointer(C_two1%da, atm%hgrad, grad_1d, grad)
!
!      if(lslope_correction) then
!        k = C_two1%ze
!        do j=C_two1%ys,C_two1%ye
!          do i=C_two1%xs,C_two1%xe
!
!            n = cross_3d([one, zero, grad(i0, k, i, j)], [zero, one, grad(i1, k, i, j)])
!            n = n/norm2(n)
!
!            fac = dot_product(solver%sun%sundir, n) / dot_product(solver%sun%sundir, [zero, zero, one])
!            edir(k-C_two1%zs+1, i-C_two1%xs+1, j-C_two1%ys+1) = edir(k-C_two1%zs+1, i-C_two1%xs+1, j-C_two1%ys+1) * fac
!          enddo
!        enddo
!      endif
!
!      if(latm_correction) then
!        do j=C_two1%ys,C_two1%ye
!          do i=C_two1%xs,C_two1%xe
!            do k=C_two1%zs,C_two1%ze
!
!              n = cross_3d([one, zero, grad(i0, k, i, j)], [zero, one, grad(i1, k, i, j)])
!              n = n/norm2(n)
!
!              fac = dot_product(solver%sun%sundir, [zero, zero, one]) / dot_product(solver%sun%sundir, n)
!              edir(k-C_two1%zs+1, i-C_two1%xs+1, j-C_two1%ys+1) = edir(k-C_two1%zs+1, i-C_two1%xs+1, j-C_two1%ys+1) * fac
!            enddo
!          enddo
!        enddo
!      endif
!
!      call restoreVecPointer(C_two1%da, atm%hgrad, grad_1d, grad)
!    end associate
!  end subroutine
!
!  subroutine pprts_rrtmg(comm, solver, atm, ie, je,   &
!      & dx, dy, sundir,                               &
!      & albedo_thermal, albedo_solar,                 &
!      & lthermal, lsolar,                             &
!      & edir,edn,eup,abso,                            &
!      & nxproc, nyproc, icollapse,                    &
!      & opt_time, solar_albedo_2d, thermal_albedo_2d, &
!      & opt_solar_constant,                           &
!      & opt_buildings_solar, opt_buildings_thermal,   &
!      & opt_tau_solar,                                &
!      & opt_w0_solar,                                 &
!      & opt_g_solar,                                  &
!      & opt_tau_thermal,                              &
!      & lonly_initialize)
!
!    integer(mpiint), intent(in)     :: comm ! MPI Communicator
!
!    class(t_solver), intent(inout)  :: solver                       ! solver type (e.g. t_solver_8_10)
!    type(t_tenstr_atm), intent(in)  :: atm                          ! contains info on atmospheric constituents
!    integer(iintegers), intent(in)  :: ie, je                       ! local domain size in x and y direction
!    real(ireals), intent(in)        :: dx, dy                       ! horizontal grid spacing in [m]
!    real(ireals), intent(in)        :: sundir(:)                    ! cartesian sun direction, pointing away from the sun, dim(3)
!    real(ireals), intent(in)        :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum
!
!    ! Compute solar or thermal radiative transfer. Or compute both at once.
!    logical, intent(in) :: lsolar, lthermal
!
!    ! nxproc dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
!    ! nyproc dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
!    ! if not present, we let petsc decide how to decompose the fields(probably does not fit the decomposition of a host model)
!    integer(iintegers),intent(in),optional :: nxproc(:), nyproc(:)
!
!    integer(iintegers),intent(in),optional :: icollapse ! experimental, dont use it if you dont know what you are doing.
!
!    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions
!    ! and compute new solutions only after threshold estimate is exceeded.
!    ! If solar_albedo_2d is present, we use a 2D surface albedo
!    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:,:), thermal_albedo_2d(:,:), opt_solar_constant
!
!    ! buildings information, setup broadband thermal and solar albedo on faces inside the domain, not just on the surface
!    ! see definition for details on how to set it up
!    type(t_pprts_buildings), intent(inout), optional :: opt_buildings_solar
!    type(t_pprts_buildings), intent(inout), optional :: opt_buildings_thermal
!
!    ! optional optical properties with dims(nlyr(srfc to TOA), local_nx, local_ny, nr-g-points)
!    ! used e.g. for aerosol or vegetation, you can provide only tau or (tau and w0) defaults to w0=0 and g=0
!    ! note that first dimension can also be smaller than nlyr, we will fill it up from the ground,
!    ! i.e. if you provide only two layers, the two lowermost layers near the surface will be filled with addition optprops
!    real(ireals), intent(in), optional, dimension(:,:,:,:) :: opt_tau_solar, opt_w0_solar, opt_g_solar
!    real(ireals), intent(in), optional, dimension(:,:,:,:) :: opt_tau_thermal
!
!    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
!    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
!    ! the size of the merged grid. If you only want to use heating rates on the
!    ! dynamics grid, use the lower layers, i.e.,
!    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
!    ! or:
!    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
!    real(ireals),allocatable, dimension(:,:,:), intent(inout) :: edir, edn, eup  ! [nlyr+1, local_nx, local_ny ]
!    real(ireals),allocatable, dimension(:,:,:), intent(inout) :: abso            ! [nlyr  , local_nx, local_ny ]
!
!    ! if only_initialize we dont compute any radiation, merely setup the grid structures
!    logical, intent(in), optional :: lonly_initialize
!
!    ! ---------- end of API ----------------
!
!    ! Counters
!    integer(iintegers) :: ke, ke1
!
!    ! for debug purposes, can output variables into netcdf files
!    !character(default_str_len) :: output_path(2) ! [ filename, varname ]
!    !logical :: lfile_exists
!
!    integer(mpiint) :: myid, ierr
!    logical :: lrrtmg_only, lskip_thermal, lskip_solar, ldisort_only, lprint_atm, lflg
!
!    integer(iintegers) :: pprts_icollapse
!
!    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
!
!    if(ldebug) then ! make sure that all ranks give the same option for lsolar and lthermal
!      if(.not. mpi_logical_all_same(comm, lsolar)) &
!        call CHKERR(1_mpiint, 'all ranks have to give the same value for lsolar')
!      if(.not. mpi_logical_all_same(comm, lthermal)) &
!        call CHKERR(1_mpiint, 'all ranks have to give the same value for lthermal')
!    endif
!
!    if(.not.solver%linitialized) call read_commandline_options(comm) ! so that tenstream.options file are read in
!
!    lrrtmg_only=.False. ! by default use normal tenstream solver
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!      "-rrtmg_only" , lrrtmg_only , lflg , ierr) ;call CHKERR(ierr)
!
!    ldisort_only = .False.
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!      "-disort_only" , ldisort_only , lflg , ierr) ;call CHKERR(ierr)
!
!    pprts_icollapse = get_arg(i1, icollapse)
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!                             "-pprts_collapse" , pprts_icollapse, lflg , ierr) ;call CHKERR(ierr)
!    if(pprts_icollapse.eq.-1) then
!      if(ldebug.and.myid.eq.0) print *,'Collapsing background atmosphere', atm%atm_ke
!      pprts_icollapse = atm%atm_ke ! collapse the complete background atmosphere
!    endif
!    if(pprts_icollapse.ne.i1) then
!      if(ldebug.and.myid.eq.0) print *,'Collapsing atmosphere', pprts_icollapse
!      if(ldisort_only) call CHKERR(1_mpiint, 'for disort_only, pprts_collapse has to be set to 1')
!      if(lrrtmg_only) call CHKERR(1_mpiint, 'for rrtmg_only, pprts_collapse has to be set to 1')
!    endif
!    if(myid.eq.0.and.ldebug) print *,'icollapse', pprts_icollapse, 'disort_only', ldisort_only, 'rrtmg_only', lrrtmg_only
!
!    ke1 = ubound(atm%plev,1)
!    ke = ubound(atm%tlay,1)
!
!    lprint_atm = ldebug
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!      "-pprts_rrtmg_atm_view" , lprint_atm , lflg , ierr) ;call CHKERR(ierr)
!    if(lprint_atm .and. myid.eq.0) then
!      call print_tenstr_atm(atm)
!    endif
!
!    if(.not.solver%linitialized) then
!      call init_pprts_rrtmg(comm, solver, &
!        dx, dy, atm%dz, &
!        sundir, &
!        ie, je, ke, &
!        nxproc, nyproc, pprts_icollapse)
!    endif
!
!    ! Allocate space for results -- for integrated values...
!    if(.not.allocated(edn )) allocate(edn (solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    if(.not.allocated(eup )) allocate(eup (solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    if(.not.allocated(abso)) allocate(abso(solver%C_one%zm , solver%C_one%xm , solver%C_one%ym ))
!    edn = zero
!    eup = zero
!    abso= zero
!
!    if(get_arg(.False., lonly_initialize)) return
!
!    lskip_thermal = .False.
!    call get_petsc_opt(PETSC_NULL_CHARACTER , &
!      "-skip_thermal" , lskip_thermal, lflg , ierr) ;call CHKERR(ierr)
!    if(lthermal.and..not.lskip_thermal)then
!      call PetscLogStagePush(log_events%stage_rrtmg_thermal, ierr); call CHKERR(ierr)
!      call compute_thermal(                           &
!        & solver,                                     &
!        & atm,                                        &
!        & ie, je, ke, ke1,                            &
!        & albedo_thermal,                             &
!        & edn, eup, abso,                             &
!        & opt_time          = opt_time,               &
!        & lrrtmg_only       = lrrtmg_only,            &
!        & thermal_albedo_2d = thermal_albedo_2d,      &
!        & opt_buildings     = opt_buildings_thermal,  &
!        & opt_tau           = opt_tau_thermal         &
!        & )
!      call PetscLogStagePop(ierr); call CHKERR(ierr) ! pop solver%logs%stage_rrtmg_thermal
!    endif
!
!    if(lsolar.and..not.allocated(edir)) allocate(edir (solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    if(allocated(edir)) edir = zero
!
!    if(lsolar) then
!      lskip_solar = .False.
!      call get_petsc_opt(PETSC_NULL_CHARACTER , &
!        "-skip_solar" , lskip_solar, lflg , ierr) ;call CHKERR(ierr)
!      if(.not.lskip_solar) then
!        call PetscLogStagePush(log_events%stage_rrtmg_solar, ierr); call CHKERR(ierr)
!        call compute_solar(                          &
!          & solver,                                  &
!          & atm,                                     &
!          & ie, je, ke,                              &
!          & sundir, albedo_solar,                    &
!          & edir, edn, eup, abso,                    &
!          & opt_time          = opt_time,            &
!          & solar_albedo_2d   = solar_albedo_2d,     &
!          & lrrtmg_only       = lrrtmg_only,         &
!          & opt_solar_constant= opt_solar_constant,  &
!          & opt_buildings     = opt_buildings_solar, &
!          & opt_tau           = opt_tau_solar,       &
!          & opt_w0            = opt_w0_solar,        &
!          & opt_g             = opt_g_solar          &
!          & )
!        call PetscLogStagePop(ierr); call CHKERR(ierr) ! pop solver%logs%stage_rrtmg_solar
!      endif
!    endif
!
!    call smooth_surface_fluxes(solver, edn, eup)
!    if(lsolar) call slope_correction_fluxes(solver, edir)
!
!    call dump_results()
!
!    !if(myid.eq.0 .and. ldebug) then
!    !  if(present(opt_time)) then
!    !    write (output_path(1), "(A,I6.6,L1,L1,A3)") "3dout_",int(opt_time),lthermal,lsolar, '.nc'
!    !  else
!    !    write (output_path(1), "(A,L1,L1,A3)") "3dout_",lthermal,lsolar, '.nc'
!    !  endif
!    !  inquire( file=trim(output_path(1)), exist=lfile_exists )
!    !  if(.not. lfile_exists) then
!    !    output_path(2) = 'edir' ; call ncwrite(output_path, edir, i)
!    !    output_path(2) = 'edn'  ; call ncwrite(output_path, edn , i)
!    !    output_path(2) = 'eup'  ; call ncwrite(output_path, eup , i)
!    !    output_path(2) = 'abso' ; call ncwrite(output_path, abso, i)
!    !    if(present(d_lwc)) then
!    !      output_path(2) = 'lwc'  ; call ncwrite(output_path, d_lwc, i)
!    !    endif
!    !    if(present(d_iwc)) then
!    !      output_path(2) = 'iwc'  ; call ncwrite(output_path, d_iwc, i)
!    !    endif
!    !  endif
!    !endif
!    contains
!
!      subroutine dump_variable(var, dm, dumpstring, varname)
!        real(ireals), intent(in) :: var(:,:,:)
!        type(tDM), intent(in) :: dm
!        character(len=*), intent(in) :: dumpstring, varname
!        character(len=default_str_len) :: vname
!        logical :: lflg
!        type(tVec) :: dumpvec
!
!        call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
!          trim(dumpstring), lflg , ierr) ;call CHKERR(ierr)
!        if(lflg) then
!          vname = trim(varname)
!          if(present(opt_time)) vname = trim(vname)//'.t'//trim(adjustl(toStr(opt_time)))
!          if(.not.lsolar.or..not.lthermal) vname = trim(vname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
!          call DMGetGlobalVector(dm ,dumpvec ,ierr); call CHKERR(ierr)
!          call PetscObjectSetName(dumpvec, trim(vname), ierr); call CHKERR(ierr)
!          call f90VecToPetsc(var, dm, dumpvec)
!
!          call PetscObjectViewFromOptions(dumpvec, PETSC_NULL_VEC, &
!            trim(dumpstring), ierr); call CHKERR(ierr)
!          call DMRestoreGlobalVector(dm ,dumpvec ,ierr); call CHKERR(ierr)
!        endif
!      end subroutine
!
!      subroutine dump_results()
!        logical :: lflg
!        character(len=default_str_len) :: fname
!        integer(mpiint) :: ierr
!
!        if(lsolar) then
!          call dump_variable(edir, solver%C_one1%da, "-pprts_rrtmg_dump_edir", "edir")
!        endif
!        call dump_variable(edn , solver%C_one1%da, "-pprts_rrtmg_dump_edn", "edn")
!        call dump_variable(eup , solver%C_one1%da, "-pprts_rrtmg_dump_eup", "eup")
!        call dump_variable(abso, solver%C_one%da,  "-pprts_rrtmg_dump_abso", "abso")
!
!        if(lsolar.and.present(opt_buildings_solar)) then
!          call get_petsc_opt(PETSC_NULL_CHARACTER, &
!            & '-pprts_rrtmg_xdmf_buildings_solar', fname, lflg, ierr); call CHKERR(ierr)
!          if(lflg) then
!            if(present(opt_time)) fname = trim(fname)//'.t'//trim(adjustl(toStr(opt_time)))
!            if(.not.lsolar.or..not.lthermal) fname = trim(fname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
!            call xdmf_pprts_buildings(solver, opt_buildings_solar, fname, ierr, verbose=.True.); call CHKERR(ierr)
!          endif
!        endif
!
!        if(lthermal.and.present(opt_buildings_thermal)) then
!          call get_petsc_opt(PETSC_NULL_CHARACTER, &
!            & '-pprts_rrtmg_xdmf_buildings_thermal', fname, lflg, ierr); call CHKERR(ierr)
!          if(lflg) then
!            if(present(opt_time)) fname = trim(fname)//'.t'//trim(adjustl(toStr(opt_time)))
!            if(.not.lsolar.or..not.lthermal) fname = trim(fname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
!            call xdmf_pprts_buildings(solver, opt_buildings_thermal, fname, ierr, verbose=.True.); call CHKERR(ierr)
!          endif
!        endif
!
!        call get_petsc_opt(PETSC_NULL_CHARACTER, &
!          & '-pprts_rrtmg_xdmf', fname, lflg, ierr); call CHKERR(ierr)
!        if(lflg) then
!          if(present(opt_time)) fname = trim(fname)//'.t'//trim(adjustl(toStr(opt_time)))
!          if(.not.lsolar.or..not.lthermal) fname = trim(fname)//'.sol'//toStr(lsolar)//'.th'//toStr(lthermal)
!          if(lsolar) then
!            call xdmf_pprts_srfc_flux(solver, fname, edn, eup, ierr, edir=edir, verbose=.True.); call CHKERR(ierr)
!          else
!            call xdmf_pprts_srfc_flux(solver, fname, edn, eup, ierr, verbose=.True.); call CHKERR(ierr)
!          endif
!        endif
!      end subroutine
!  end subroutine
!
!  subroutine compute_thermal( &
!      & solver,               &
!      & atm,                  &
!      & ie, je, ke, ke1,      &
!      & albedo,               &
!      & edn,                  &
!      & eup,                  &
!      & abso,                 &
!      & opt_time,             &
!      & lrrtmg_only,          &
!      & thermal_albedo_2d,    &
!      & opt_buildings,        &
!      & opt_tau)
!
!    use m_tenstr_rrlw_wvn, only : ngb, wavenum1, wavenum2
!    use m_tenstr_parrrtm, only: ngptlw
!
!    class(t_solver),    intent(inout) :: solver
!    type(t_tenstr_atm), intent(in), target :: atm
!    integer(iintegers), intent(in)    :: ie, je, ke,ke1
!
!    real(ireals),intent(in) :: albedo
!
!    real(ireals),intent(inout),dimension(:,:,:) :: edn, eup, abso
!
!    real(ireals), optional, intent(in) :: opt_time, thermal_albedo_2d(:,:)
!    logical, optional, intent(in) :: lrrtmg_only
!    type(t_pprts_buildings), intent(inout), optional :: opt_buildings
!    real(ireals), intent(in), optional, dimension(:,:,:,:) :: opt_tau
!
!    real(ireals),allocatable, target, dimension(:,:,:,:) :: tau, Bfrac  ! [nlyr, ie, je, ngptlw]
!    real(ireals),allocatable, dimension(:,:,:) :: kabs, ksca, g, Blev   ! [nlyr(+1), local_nx, local_ny]
!    real(ireals),allocatable, dimension(:,:)   :: Bsrfc                 ! [local_nx, local_ny]
!    real(ireals),allocatable, dimension(:,:,:), target :: spec_edn,spec_eup,spec_abso  ! [nlyr(+1), local_nx, local_ny ]
!    real(ireals),allocatable, dimension(:) :: integral_coeff            ! [nlyr]
!
!    real(ireals), allocatable, dimension(:,:,:) :: ptau, pBfrac
!    real(ireals), pointer, dimension(:,:,:) :: patm_dz
!    real(ireals), pointer, dimension(:,:) :: pedn, peup, pabso
!
!    real(ireals) :: col_albedo, col_tskin(1)
!
!    integer(iintegers) :: i, j, k, icol, ib, current_ibnd, spectral_bands(2)
!    integer(iintegers) :: ak, idx(4)
!    logical :: need_any_new_solution, lflg
!
!    type(t_pprts_buildings), allocatable :: spec_buildings
!
!    integer(mpiint) :: myid, ierr
!
!    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
!
!    allocate(spec_edn (solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    allocate(spec_eup (solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    allocate(spec_abso(solver%C_one%zm , solver%C_one%xm , solver%C_one%ym ))
!    if(present(opt_buildings)) then
!
!      call clone_buildings(opt_buildings, spec_buildings, l_copy_data=.False., ierr=ierr); call CHKERR(ierr)
!      if(.not.allocated(spec_buildings%albedo )) allocate(spec_buildings%albedo (size(opt_buildings%albedo)))
!      if(.not.allocated(opt_buildings%incoming)) allocate(opt_buildings%incoming(size(opt_buildings%iface)))
!      if(.not.allocated(opt_buildings%outgoing)) allocate(opt_buildings%outgoing(size(opt_buildings%iface)))
!      spec_buildings%albedo  = opt_buildings%albedo
!      opt_buildings%incoming = zero
!      opt_buildings%outgoing = zero
!
!      if(.not.(allocated(opt_buildings%temp))) call CHKERR(1_mpiint, 'Thermal computation but opt_buildings%temp is not allocated')
!      if((allocated(opt_buildings%planck))) call CHKERR(1_mpiint, 'Thermal computation but opt_buildings%planck is allocated... '//&
!        & 'you should only provide temperatures. I`ll take care of planck emissison')
!
!      allocate(spec_buildings%temp  (size(opt_buildings%temp)))
!      allocate(spec_buildings%planck(size(opt_buildings%temp)))
!    endif
!
!
!    need_any_new_solution=.False.
!    do ib=1,ngptlw
!      if(need_new_solution(solver%comm, solver%solutions(500+ib), opt_time, solver%lenable_solutions_err_estimates)) &
!        need_any_new_solution=.True.
!    enddo
!    if(.not.need_any_new_solution) then
!      do ib=1,ngptlw
!        if(present(opt_buildings)) then
!          call pprts_get_result(solver, &
!            & spec_edn, spec_eup, spec_abso, &
!            & opt_solution_uid=500+ib, &
!            & opt_buildings=spec_buildings)
!          opt_buildings%incoming = opt_buildings%incoming + spec_buildings%incoming
!          opt_buildings%outgoing = opt_buildings%outgoing + spec_buildings%outgoing
!        else
!          call pprts_get_result(solver, &
!            & spec_edn, spec_eup, spec_abso, &
!            & opt_solution_uid=500+ib)
!        endif
!        edn  = edn  + spec_edn
!        eup  = eup  + spec_eup
!        abso = abso + spec_abso
!      enddo
!
!      if(present(opt_buildings)) then
!        call destroy_buildings(spec_buildings, ierr)
!      endif
!      return
!    endif
!
!    ! Compute optical properties with RRTMG
!    allocate(tau  (ke, i1:ie, i1:je, ngptlw))
!    allocate(Bfrac(ke1, i1:ie, i1:je, ngptlw))
!    allocate(ptau  (ke, i1, ngptlw))
!    allocate(pBfrac(ke, i1, ngptlw))
!    allocate(integral_coeff(ke))
!
!    col_albedo = albedo
!
!    if(lrrtmg_only) then
!      do j=i1,je
!        do i=i1,ie
!          icol =  i+(j-1)*ie
!
!          pedn (1:ke1, 1:1) => spec_edn (:,i,j)
!          peup (1:ke1, 1:1) => spec_eup (:,i,j)
!          pabso(1:ke , 1:1) => spec_abso(:,i,j)
!
!          do k=1,ke
!            integral_coeff(k) = vert_integral_coeff(atm%plev(k,icol), atm%plev(k+1,icol))
!          enddo
!
!          if(present(thermal_albedo_2d)) col_albedo = thermal_albedo_2d(i,j)
!          if(allocated(atm%tskin)) then
!            col_tskin = atm%tskin(icol)
!          else
!            col_tskin = atm%tlev(1,icol)
!          endif
!
!          call optprop_rrtm_lw(i1, ke, col_albedo,                          &
!            atm%plev(:,icol), atm%tlev(:,icol),                             &
!            atm%tlay(:, icol), col_tskin,                                   &
!            atm%h2o_lay(:,icol), atm%o3_lay (:,icol), atm%co2_lay(:,icol),  &
!            atm%ch4_lay(:,icol), atm%n2o_lay(:,icol), atm%o2_lay (:,icol) , &
!            atm%lwc(:,icol)*integral_coeff, atm%reliq(:, icol),             &
!            atm%iwc(:,icol)*integral_coeff, atm%reice(:, icol),             &
!            tau=ptau, Bfrac=pBfrac, opt_lwuflx=peup, opt_lwdflx=pedn, opt_lwhr=pabso, &
!            log_event=log_events%rrtmg_optprop_lw)
!
!          tau  (:,i,j,:) = ptau(:,i1,:)
!          Bfrac(2:ke1,i,j,:) = pBfrac(:,i1,:)
!
!          eup(:,i,j)  = eup(:,i,j)  + reverse(spec_eup (:,i,j))
!          edn(:,i,j)  = edn(:,i,j)  + reverse(spec_edn (:,i,j))
!          abso(:,i,j) = abso(:,i,j) + reverse( &
!            (spec_edn(2:ke1,i,j)-spec_edn(1:ke,i,j)+spec_eup(1:ke,i,j)-spec_eup(2:ke1,i,j)) / &
!             atm%dz(:,icol))
!
!        enddo
!      enddo
!      return
!    else
!      do j=i1,je
!        do i=i1,ie
!          icol =  i+(j-1)*ie
!          do k=1,ke
!            integral_coeff(k) = vert_integral_coeff(atm%plev(k,icol), atm%plev(k+1,icol))
!          enddo
!
!          if (present(thermal_albedo_2d)) col_albedo = thermal_albedo_2d(i,j)
!          if(allocated(atm%tskin)) then
!            col_tskin = atm%tskin(icol)
!          else
!            col_tskin = atm%tlev(1,icol)
!          endif
!
!          call optprop_rrtm_lw(i1, ke, col_albedo,                          &
!            atm%plev(:,icol), atm%tlev(:, icol),                            &
!            atm%tlay(:, icol), col_tskin,                                   &
!            atm%h2o_lay(:,icol), atm%o3_lay (:,icol), atm%co2_lay(:,icol),  &
!            atm%ch4_lay(:,icol), atm%n2o_lay(:,icol), atm%o2_lay (:,icol) , &
!            atm%lwc(:,icol)*integral_coeff, atm%reliq(:, icol),             &
!            atm%iwc(:,icol)*integral_coeff, atm%reice(:, icol),             &
!            tau=ptau, Bfrac=pBfrac,                                         &
!            log_event=log_events%rrtmg_optprop_lw)
!
!          tau  (:,i,j,:) = ptau(:,i1,:)
!          Bfrac(2:ke1,i,j,:) = pBfrac(:,i1,:)
!        enddo
!      enddo
!    endif
!    Bfrac(1,:,:,:) = Bfrac(2,:,:,:)
!
!    call add_optional_optprop(tau=tau, opt_tau=opt_tau)
!
!    allocate(kabs (ke , i1:ie, i1:je))
!    allocate(Blev (ke1, i1:ie, i1:je))
!    allocate(Bsrfc(i1:ie, i1:je))
!
!    ! rrtmg_lw does not support thermal scattering... set to zero
!    allocate(ksca (ke , i1:ie, i1:je), source=zero)
!    allocate(g    (ke , i1:ie, i1:je), source=zero)
!
!    current_ibnd = -1 ! current lw band
!
!    spectral_bands = get_spectral_bands(solver%comm, i1, int(ngptlw, iintegers))
!
!    if(compute_thermal_disort()) return
!
!    do ib=spectral_bands(1), spectral_bands(2)
!
!      if(need_new_solution(solver%comm, solver%solutions(500+ib), opt_time, solver%lenable_solutions_err_estimates)) then
!        ! divide by thickness to convert from tau to coefficients per meter
!        patm_dz(1:ke, i1:ie, i1:je) => atm%dz
!        kabs = max(zero, tau(:,:,:,ib)) / patm_dz
!        kabs = reverse(kabs)
!
!        !Compute Plank Emission for nbndlw
!        if(current_ibnd.eq.ngb(ib)) then ! still the same band, dont need to upgrade the plank emission
!          continue
!        else
!          do j=i1,je
!            do i=i1,ie
!              icol = i+(j-1)*ie
!              do k=i1,ke1
!                Blev(k,i,j) = plkint(real(wavenum1(ngb(ib))), real(wavenum2(ngb(ib))), real(atm%tlev(k, icol)))
!              enddo
!              if(allocated(atm%tskin)) then
!                Bsrfc(i,j) = plkint(real(wavenum1(ngb(ib))), real(wavenum2(ngb(ib))), real(atm%tskin(icol)))
!              else
!                Bsrfc(i,j) = Blev(1,i,j)
!              endif
!            enddo ! i
!          enddo ! j
!          current_ibnd = ngb(ib)
!
!          if(ldebug)then
!            if(any(Blev.lt.zero)) then
!              print *,'Min Max Planck:', minval(Blev), maxval(Blev), 'location min', minloc(Blev)
!              call CHKERR(1_mpiint, 'Found a negative Planck emission, this is not physical! Aborting...')
!            endif
!          endif
!
!          if(present(opt_buildings)) then
!            ! Use spec_buildings%temperature array to hold raw planck values ...
!            do k = 1, size(opt_buildings%temp)
!              spec_buildings%temp(k) = &
!                & plkint(real(wavenum1(ngb(ib))), real(wavenum2(ngb(ib))), real(opt_buildings%temp(k)))
!            enddo
!          endif
!        endif
!
!        if(present(opt_buildings)) then
!          ! ... and use spec_buildings planck array for raw planck values multiplied with the planck fractions
!          do icol = 1, size(spec_buildings%temp)
!            call ind_1d_to_nd(spec_buildings%da_offsets, spec_buildings%iface(icol), idx)
!            associate(d => idx(1), k => idx(2), i => idx(3), j => idx(4))
!              if(d.eq.PPRTS_BOT_FACE) then
!                ak = atmk(solver%atm, k+1)
!              else
!                ak = atmk(solver%atm, k)
!              endif
!              spec_buildings%planck(icol) = spec_buildings%temp(icol) * &
!                & Bfrac(1+size(Bfrac,dim=1)-ak,i,j,ib)
!            end associate
!          enddo
!        endif
!
!        call set_optical_properties(solver, albedo, kabs, ksca, g, &
!            reverse(Blev*Bfrac(:,:,:,ib)), planck_srfc=Bsrfc*Bfrac(1,:,:,ib), &
!            albedo_2d=thermal_albedo_2d)
!        call solve_pprts(solver, &
!          & lthermal=.True., &
!          & lsolar=.False., &
!          & edirTOA=zero, &
!          & opt_solution_uid=500+ib, &
!          & opt_solution_time=opt_time, &
!          & opt_buildings=spec_buildings)
!      endif
!
!      if(present(opt_buildings)) then
!        call pprts_get_result(solver, &
!          & spec_edn, spec_eup, spec_abso, &
!          & opt_solution_uid=500+ib, &
!          & opt_buildings=spec_buildings)
!        opt_buildings%incoming = opt_buildings%incoming + spec_buildings%incoming
!        opt_buildings%outgoing = opt_buildings%outgoing + spec_buildings%outgoing
!      else
!        call pprts_get_result(solver, &
!          & spec_edn, spec_eup, spec_abso, &
!          & opt_solution_uid=500+ib)
!      endif
!
!      edn  = edn  + spec_edn
!      eup  = eup  + spec_eup
!      abso = abso + spec_abso
!    enddo ! ib 1 -> nbndlw , i.e. spectral integration
!
!    if(present(opt_buildings)) then
!      call destroy_buildings(spec_buildings, ierr)
!    endif
!    contains
!      function compute_thermal_disort() result(ldisort_only)
!        logical :: ldisort_only
!        integer(iintegers) :: nstreams
!        real :: mu0, S0, col_albedo, wvnms(2), col_tskin
!        real, dimension(size(tau,1))   :: col_Bfrac, col_dtau, col_w0, col_g
!        real, dimension(size(tau,1)+1) :: col_temper
!        real, dimension(size(edn,1))   :: RFLDIR, RFLDN, FLUP, DFDT, UAVG
!
!        ldisort_only = .False.
!        call get_petsc_opt(PETSC_NULL_CHARACTER , &
!          "-disort_only" , ldisort_only , lflg , ierr) ;call CHKERR(ierr)
!
!        if(ldisort_only) then
!          nstreams = 16
!          call get_petsc_opt(PETSC_NULL_CHARACTER , &
!            "-disort_streams" , nstreams , lflg , ierr) ;call CHKERR(ierr)
!
!          mu0 = 0
!          S0  = 0
!          col_w0 = 0
!          col_g  = 0
!
!          col_albedo = real(albedo)
!          do j=1,je
!            do i=1,ie
!              icol =  i+(j-1)*ie
!
!              if(present(thermal_albedo_2d)) col_albedo = real(thermal_albedo_2d(i,j))
!
!              col_temper = real(reverse(atm%tlev(:, icol)))
!              if(allocated(atm%tskin)) then
!                col_tskin = real(atm%tskin(icol))
!              else
!                col_tskin = real(atm%tlev(1,icol))
!              endif
!
!
!              do ib=spectral_bands(1), spectral_bands(2)
!                col_Bfrac  = real(reverse(Bfrac(1:ke,i,j,ib)))
!                col_dtau   = max(tiny(col_dtau), real(reverse(tau(:,i,j,ib))))
!                wvnms = [real(wavenum1(ngb(ib))), real(wavenum2(ngb(ib)))]
!
!                call default_flx_computation(&
!                  mu0, &
!                  S0, &
!                  col_albedo, &
!                  col_tskin, &
!                  .True., wvnms, col_Bfrac, &
!                  col_dtau, &
!                  col_w0,   &
!                  col_g,    &
!                  col_temper, &
!                  RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
!                  int(nstreams), lverbose=.False.)
!
!                eup (:,i,j) = eup (:,i,j) + FLUP
!                edn (:,i,j) = edn (:,i,j) + RFLDN
!              enddo ! ib 1 -> nbndsw , i.e. spectral integration
!            enddo
!          enddo
!
!          patm_dz(1:ke, i1:ie, i1:je) => atm%dz
!          abso =( edn (1:ke,:,:) - edn (2:ke+1,:,:)   &
!                - eup (1:ke,:,:) + eup (2:ke+1,:,:) ) / reverse(patm_dz)
!        endif
!      end function
!  end subroutine compute_thermal
!
!  subroutine compute_solar( &
!      & solver,             &
!      & atm,                &
!      & ie,                 &
!      & je,                 &
!      & ke,                 &
!      & sundir,             &
!      & albedo,             &
!      & edir,               &
!      & edn,                &
!      & eup,                &
!      & abso,               &
!      & opt_time,           &
!      & solar_albedo_2d,    &
!      & lrrtmg_only,        &
!      & opt_solar_constant, &
!      & opt_buildings,      &
!      & opt_tau,            &
!      & opt_w0,             &
!      & opt_g)
!
!      use m_tenstr_parrrsw, only: ngptsw
!      use m_tenstr_rrtmg_sw_spcvrt, only: tenstr_solsrc
!
!    class(t_solver), intent(inout)  :: solver
!    type(t_tenstr_atm), intent(in), target :: atm
!    integer(iintegers),intent(in)   :: ie, je, ke
!
!    real(ireals),intent(in) :: albedo
!    real(ireals),intent(in) :: sundir(3)
!
!    real(ireals),intent(inout),dimension(:,:,:) :: edir, edn, eup, abso
!
!    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:,:)
!    logical, optional, intent(in) :: lrrtmg_only
!    real(ireals), intent(in), optional :: opt_solar_constant
!    type(t_pprts_buildings), intent(inout), optional :: opt_buildings
!    real(ireals), intent(in), optional, dimension(:,:,:,:) :: opt_tau, opt_w0, opt_g
!
!    real(ireals) :: edirTOA
!
!    real(ireals),allocatable, dimension(:,:,:,:) :: tau, w0, g          ! [nlyr, ie, je, ngptsw]
!    real(ireals),allocatable, dimension(:,:,:)   :: kabs, ksca, kg      ! [nlyr, local_nx, local_ny]
!    real(ireals),allocatable, dimension(:,:,:), target   :: spec_edir,spec_abso ! [nlyr(+1), local_nx, local_ny ]
!    real(ireals),allocatable, dimension(:,:,:), target   :: spec_edn, spec_eup  ! [nlyr(+1), local_nx, local_ny ]
!    real(ireals),allocatable, dimension(:) :: integral_coeff            ! [nlyr]
!
!    real(ireals), allocatable, dimension(:,:,:) :: ptau, pw0, pg
!    real(ireals), pointer, dimension(:,:,:) :: patm_dz
!    real(ireals), pointer, dimension(:,:) :: pedir, pedn, peup, pabso
!
!    real(ireals) :: col_albedo
!
!    integer(iintegers) :: i, j, k, icol, ib
!    logical :: need_any_new_solution
!
!    type(t_pprts_buildings), allocatable :: spec_buildings
!
!    logical :: lflg
!    integer(iintegers) :: spectral_bands(2)
!    integer(mpiint) :: ierr
!
!    allocate(spec_edir(solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    allocate(spec_edn (solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    allocate(spec_eup (solver%C_one1%zm, solver%C_one1%xm, solver%C_one1%ym))
!    allocate(spec_abso(solver%C_one%zm , solver%C_one%xm , solver%C_one%ym ))
!
!    if(present(opt_buildings)) then
!      call clone_buildings(opt_buildings, spec_buildings, l_copy_data=.False., ierr=ierr); call CHKERR(ierr)
!      if(.not.allocated(spec_buildings%albedo )) allocate(spec_buildings%albedo (size(opt_buildings%albedo)))
!      if(.not.allocated(opt_buildings%edir    )) allocate(opt_buildings%edir    (size(opt_buildings%iface)))
!      if(.not.allocated(opt_buildings%incoming)) allocate(opt_buildings%incoming(size(opt_buildings%iface)))
!      if(.not.allocated(opt_buildings%outgoing)) allocate(opt_buildings%outgoing(size(opt_buildings%iface)))
!      spec_buildings%albedo  = opt_buildings%albedo
!      opt_buildings%edir     = zero
!      opt_buildings%incoming = zero
!      opt_buildings%outgoing = zero
!    endif
!
!    need_any_new_solution=.False.
!    do ib=1,ngptsw
!      if(need_new_solution(solver%comm, solver%solutions(ib), opt_time, solver%lenable_solutions_err_estimates)) &
!        need_any_new_solution=.True.
!    enddo
!    if(.not.need_any_new_solution) then
!      do ib=1,ngptsw
!        if(present(opt_buildings)) then
!          call pprts_get_result(solver, &
!            & spec_edn, spec_eup, spec_abso, spec_edir, &
!            & opt_solution_uid=ib, &
!            & opt_buildings=spec_buildings)
!          opt_buildings%edir = opt_buildings%edir + spec_buildings%edir
!          opt_buildings%incoming = opt_buildings%incoming + spec_buildings%incoming
!          opt_buildings%outgoing = opt_buildings%outgoing + spec_buildings%outgoing
!        else
!          call pprts_get_result(solver, &
!            & spec_edn, spec_eup, spec_abso, spec_edir, &
!            & opt_solution_uid=ib)
!        endif
!        edir = edir + spec_edir
!        edn  = edn  + spec_edn
!        eup  = eup  + spec_eup
!        abso = abso + spec_abso
!      enddo
!
!      if(present(opt_buildings)) then
!        call destroy_buildings(spec_buildings, ierr)
!      endif
!      return
!    endif
!
!    ! Compute optical properties with RRTMG
!    allocate(tau(ke, i1:ie, i1:je, ngptsw))
!    allocate(w0 (ke, i1:ie, i1:je, ngptsw))
!    allocate(g  (ke, i1:ie, i1:je, ngptsw))
!    allocate(ptau(ke, i1, ngptsw))
!    allocate(pw0 (ke, i1, ngptsw))
!    allocate(pg  (ke, i1, ngptsw))
!
!    allocate(integral_coeff(ke))
!
!    call set_angles(solver, sundir)
!
!    if(lrrtmg_only) then
!      do j=1,je
!        do i=1,ie
!          icol =  i+(j-1)*ie
!
!          pEdir(1:size(edir,1), 1:1) => spec_edir(:,i,j)
!          pEdn (1:size(edn ,1), 1:1) => spec_edn (:,i,j)
!          pEup (1:size(eup ,1), 1:1) => spec_eup (:,i,j)
!          pabso(1:size(abso,1), 1:1) => spec_abso(:,i,j)
!
!          do k=1,ke
!            integral_coeff(k) = vert_integral_coeff(atm%plev(k,icol), atm%plev(k+1,icol))
!          enddo
!
!          if(present(solar_albedo_2d)) then
!            col_albedo = solar_albedo_2d(i,j)
!          else
!            col_albedo = albedo
!          endif
!          call optprop_rrtm_sw(i1, ke, &
!            solver%sun%theta, col_albedo, &
!            atm%plev(:,icol), atm%tlev(:,icol), atm%tlay(:,icol), &
!            atm%h2o_lay(:,icol), atm%o3_lay(:,icol), atm%co2_lay(:,icol), &
!            atm%ch4_lay(:,icol), atm%n2o_lay(:,icol), atm%o2_lay(:,icol), &
!            atm%lwc(:,icol)*integral_coeff, atm%reliq(:,icol), &
!            atm%iwc(:,icol)*integral_coeff, atm%reice(:,icol), &
!            ptau, pw0, pg, &
!            opt_swdirflx=pEdir, opt_swuflx=pEup, &
!            opt_swdflx=pEdn, opt_swhr=pabso, &
!            opt_solar_constant=opt_solar_constant, &
!            log_event=log_events%rrtmg_optprop_sw)
!
!          tau(:,i,j,:) = ptau(:,i1,:)
!          w0 (:,i,j,:) = pw0(:,i1,:)
!          g  (:,i,j,:) = pg(:,i1,:)
!
!          edir(:,i,j) = edir(:,i,j) + reverse(spec_edir(:, i, j))
!          eup (:,i,j) = eup (:,i,j) + reverse(spec_eup (:, i, j))
!          edn (:,i,j) = edn (:,i,j) + reverse(spec_edn (:, i, j))
!          abso(:,i,j) = abso(:,i,j) + reverse( ( &
!            - spec_edir(1:ke,i,j) + spec_edir(2:ke+1,i,j) &
!            - spec_edn (1:ke,i,j) + spec_edn (2:ke+1,i,j) &
!            + spec_eup (1:ke,i,j) - spec_eup (2:ke+1,i,j) ) / atm%dz(:,icol) )
!        enddo
!      enddo
!      return
!    else
!      do j=1,je
!        do i=1,ie
!          icol =  i+(j-1)*ie
!          do k=1,ke
!            integral_coeff(k) = vert_integral_coeff(atm%plev(k,icol), atm%plev(k+1,icol))
!          enddo
!
!          if(present(solar_albedo_2d)) then
!            col_albedo = solar_albedo_2d(i,j)
!          else
!            col_albedo = albedo
!          endif
!
!          call optprop_rrtm_sw(i1, ke, &
!            solver%sun%theta, col_albedo, &
!            atm%plev(:,icol), atm%tlev(:,icol), atm%tlay(:,icol), &
!            atm%h2o_lay(:,icol), atm%o3_lay(:,icol), atm%co2_lay(:,icol), &
!            atm%ch4_lay(:,icol), atm%n2o_lay(:,icol), atm%o2_lay(:,icol), &
!            atm%lwc(:,icol)*integral_coeff, atm%reliq(:,icol), &
!            atm%iwc(:,icol)*integral_coeff, atm%reice(:,icol), &
!            ptau, pw0, pg, &
!            opt_solar_constant=opt_solar_constant, &
!            log_event=log_events%rrtmg_optprop_sw)
!
!          tau(:,i,j,:) = ptau(:,i1,:)
!          w0 (:,i,j,:) = pw0(:,i1,:)
!          g  (:,i,j,:) = pg(:,i1,:)
!        enddo
!      enddo
!    endif
!
!    call add_optional_optprop(tau, w0, g, opt_tau, opt_w0, opt_g)
!
!    w0 = min(one, max(zero, w0))
!
!    spectral_bands = get_spectral_bands(solver%comm, i1, int(ngptsw, iintegers))
!
!    if(compute_solar_disort()) return
!
!    allocate(kabs(ke , i1:ie, i1:je))
!    allocate(ksca(ke , i1:ie, i1:je))
!    allocate(kg  (ke , i1:ie, i1:je))
!
!    do ib=spectral_bands(1), spectral_bands(2)
!
!      if(need_new_solution(solver%comm, solver%solutions(ib), opt_time, solver%lenable_solutions_err_estimates)) then
!        patm_dz(1:ke, i1:ie, i1:je) => atm%dz
!        kabs = max(zero, tau(:,:,:,ib)) * (one - w0(:,:,:,ib))
!        ksca = max(zero, tau(:,:,:,ib)) * w0(:,:,:,ib)
!        kg   = min(one, max(zero, g(:,:,:,ib)))
!        kabs = reverse(kabs / patm_dz)
!        ksca = reverse(ksca / patm_dz)
!        kg   = reverse(kg)
!
!        if(present(opt_solar_constant)) then
!          edirTOA = tenstr_solsrc(ib) / sum(tenstr_solsrc) * opt_solar_constant
!        else
!          edirTOA = tenstr_solsrc(ib)
!        endif
!
!        call set_optical_properties( &
!          & solver,                  &
!          & albedo,                  &
!          & kabs,                    &
!          & ksca,                    &
!          & kg,                      &
!          & albedo_2d=solar_albedo_2d)
!
!        call solve_pprts(solver, &
!          & lthermal=.False., &
!          & lsolar=.True., &
!          & edirTOA=edirTOA, &
!          & opt_solution_uid=ib,          &
!          & opt_solution_time=opt_time,   &
!          & opt_buildings=spec_buildings)
!
!      endif
!
!      if(present(opt_buildings)) then
!        call pprts_get_result(solver, &
!          & spec_edn, spec_eup, spec_abso, spec_edir, &
!          & opt_solution_uid=ib, &
!          & opt_buildings=spec_buildings)
!        opt_buildings%edir = opt_buildings%edir + spec_buildings%edir
!        opt_buildings%incoming = opt_buildings%incoming + spec_buildings%incoming
!        opt_buildings%outgoing = opt_buildings%outgoing + spec_buildings%outgoing
!      else
!        call pprts_get_result(solver, &
!          & spec_edn, spec_eup, spec_abso, spec_edir, &
!          & opt_solution_uid=ib)
!      endif
!
!      edir = edir + spec_edir
!      edn  = edn  + spec_edn
!      eup  = eup  + spec_eup
!      abso = abso + spec_abso
!
!    enddo ! ib 1 -> nbndsw , i.e. spectral integration
!
!    if(present(opt_buildings)) then
!      call destroy_buildings(spec_buildings, ierr)
!    endif
!    contains
!      function compute_solar_disort() result(ldisort_only)
!        logical :: ldisort_only
!        integer(iintegers) :: nstreams
!        real :: mu0, col_tskin
!        real, dimension(size(tau,1))   :: col_Bfrac, col_dtau, col_w0, col_g
!        real, dimension(size(tau,1)+1) :: col_temper
!        real, dimension(size(edn,1))   :: RFLDIR, RFLDN, FLUP, DFDT, UAVG
!
!        ldisort_only = .False.
!        call get_petsc_opt(PETSC_NULL_CHARACTER , &
!          "-disort_only" , ldisort_only , lflg , ierr) ;call CHKERR(ierr)
!
!        if(ldisort_only) then
!          nstreams = 16
!          call get_petsc_opt(PETSC_NULL_CHARACTER , &
!            "-disort_streams" , nstreams , lflg , ierr) ;call CHKERR(ierr)
!
!          col_temper = 0
!          col_tskin = 0
!          col_Bfrac = 1
!          col_albedo = albedo
!          do j=1,je
!            do i=1,ie
!              icol =  i+(j-1)*ie
!
!              if(present(solar_albedo_2d)) col_albedo = solar_albedo_2d(i,j)
!
!              do ib=spectral_bands(1), spectral_bands(2)
!                if(present(opt_solar_constant)) then
!                  edirTOA = tenstr_solsrc(ib) /sum(tenstr_solsrc) * opt_solar_constant
!                else
!                  edirTOA = tenstr_solsrc(ib)
!                endif
!
!                col_dtau   = max(tiny(col_dtau), real(reverse(tau(:,i,j,ib))))
!                col_w0     = max(tiny(col_w0  ), real(reverse(w0 (:,i,j,ib))))
!                col_g      = max(tiny(col_g   ), real(reverse(g  (:,i,j,ib))))
!
!                mu0 = real(cos(deg2rad(solver%sun%theta)))
!                call default_flx_computation(&
!                  mu0, &
!                  real(edirTOA), &
!                  real(col_albedo), &
!                  col_tskin, &
!                  .False., [0., 0.], col_Bfrac, &
!                  col_dtau, &
!                  col_w0,   &
!                  col_g,    &
!                  col_temper, &
!                  RFLDIR, RFLDN, FLUP, DFDT, UAVG, &
!                  int(nstreams), lverbose=.False.)
!
!                edir(:,i,j) = edir(:,i,j) + RFLDIR
!                eup (:,i,j) = eup (:,i,j) + FLUP
!                edn (:,i,j) = edn (:,i,j) + RFLDN
!              enddo ! ib 1 -> nbndsw , i.e. spectral integration
!            enddo
!          enddo
!
!          patm_dz(1:ke, i1:ie, i1:je) => atm%dz
!          abso =( edir(1:ke,:,:) - edir(2:ke+1,:,:)   &
!                + edn (1:ke,:,:) - edn (2:ke+1,:,:)   &
!                - eup (1:ke,:,:) + eup (2:ke+1,:,:) ) / reverse(patm_dz)
!        endif
!      end function
!  end subroutine compute_solar
!
!  subroutine destroy_pprts_rrtmg(solver, lfinalizepetsc)
!    class(t_solver)     :: solver
!    logical, intent(in) :: lfinalizepetsc
!    !TODO this should need to happen here. the reason we need this is in case that someone destroys petsc
!    !     and then initializes tenstream again, this is not valid anymore. Anyway, it should really reside in the pprts_solver object but
!    !     this currently would infer some serious refactoring to avoid cyclic dependencies because the part in pprt2plex2rayli brings in
!    !     the plexrt deps. For now, live with it but one day, we need to sort this out.
!    call destroy_rayli_info()
!
!    ! Tidy up the solver
!    call destroy_pprts(solver, lfinalizepetsc=lfinalizepetsc)
!
!  end subroutine
!
!  subroutine add_optional_optprop(tau, w0, g, opt_tau, opt_w0, opt_g)
!    real(ireals), intent(inout), optional, dimension(:,:,:,:) :: tau, w0, g
!    real(ireals), intent(in), optional, dimension(:,:,:,:) :: opt_tau, opt_w0, opt_g
!
!    real(ireals) :: tausca_old, tausca_new
!    integer(iintegers) :: k,i,j,l
!
!    call check_shape('tau', tau, opt_tau)
!    call check_shape('w0', w0, opt_w0)
!    call check_shape('g', g, opt_g)
!
!    if(present(opt_g) .and. .not.present(opt_w0)) call CHKERR(1_mpiint, 'if opt_g is provided, need also opt_w0')
!    if((present(opt_w0).or.present(opt_g)) .and. .not.present(opt_tau)) &
!      & call CHKERR(1_mpiint, 'if opt_g or opt_w0 is provided, need also opt_tau')
!
!    if(present(opt_g)) then
!      do l = 1, size(opt_tau,4)
!        do j = 1, size(opt_tau,3)
!          do i = 1, size(opt_tau,2)
!            do k = 1, size(opt_tau,1)
!              tausca_old = tau(k,i,j,l) * w0(k,i,j,l)
!              tausca_new = opt_tau(k,i,j,l) * opt_w0(k,i,j,l)
!              g(k,i,j,l) = (g(k,i,j,l) * tausca_old + opt_g (k,i,j,l) * tausca_new) / (tausca_old+tausca_new)
!              w0(k,i,j,l) = (tausca_old+tausca_new) / (tau(k,i,j,l) + opt_tau(k,i,j,l))
!              tau(k,i,j,l) = tau(k,i,j,l) + opt_tau(k,i,j,l)
!            enddo
!          enddo
!        enddo
!      enddo
!    elseif(present(opt_w0)) then
!      do l = 1, size(opt_tau,4)
!        do j = 1, size(opt_tau,3)
!          do i = 1, size(opt_tau,2)
!            do k = 1, size(opt_tau,1)
!              tausca_old = tau(k,i,j,l) * w0(k,i,j,l)
!              tausca_new = opt_tau(k,i,j,l) * opt_w0(k,i,j,l)
!              g(k,i,j,l) = g(k,i,j,l) * tausca_old / (tausca_old+tausca_new)
!              w0(k,i,j,l) = (tausca_old+tausca_new) / (tau(k,i,j,l) + opt_tau(k,i,j,l))
!              tau(k,i,j,l) = tau(k,i,j,l) + opt_tau(k,i,j,l)
!            enddo
!          enddo
!        enddo
!      enddo
!    elseif(present(opt_tau)) then
!      do l = 1, size(opt_tau,4)
!        do j = 1, size(opt_tau,3)
!          do i = 1, size(opt_tau,2)
!            do k = 1, size(opt_tau,1)
!              tau(k,i,j,l) = tau(k,i,j,l) + opt_tau(k,i,j,l)
!            enddo
!          enddo
!        enddo
!      enddo
!    endif
!
!    contains
!      subroutine check_shape(varname, var, opt_var)
!        character(len=*), intent(in) :: varname
!        real(ireals), intent(in), optional, dimension(:,:,:,:) :: var
!        real(ireals), intent(in), optional, dimension(:,:,:,:) :: opt_var
!        integer :: i
!
!        if(.not.present(opt_var)) return
!
!        if(present(var).neqv.present(opt_var)) &
!          & call CHKERR(-1_mpiint, trim(varname)//' / opt_'//trim(varname)//&
!            & ' :: cannot have one argument, need both or none '// &
!            & '('//toStr(present(var))//'/'//toStr(present(opt_var))//')')
!
!        if (size(opt_var,1) .gt. size(var,1)) then
!          call CHKERR(1_mpiint, trim(varname)//" :: first dimension (nlyr) of opt_var "// &
!            & " cannot be larger than the target array. "//new_line('')// &
!            & " shape(var)    "//toStr(shape(var))//new_line('')// &
!            & " shape(opt_var)"//toStr(shape(opt_var)))
!        endif
!        do i=2,4
!          if (size(opt_var,i) .gt. size(var,i)) then
!            call CHKERR(1_mpiint, trim(varname)//" :: dimension "//toStr(i)//" of opt_var "// &
!              & " has to be same as the target array. "//new_line('')// &
!              & " shape(var)    "//toStr(shape(var))//new_line('')// &
!              & " shape(opt_var)"//toStr(shape(opt_var)))
!          endif
!        enddo
!      end subroutine
!  end subroutine
!end module
