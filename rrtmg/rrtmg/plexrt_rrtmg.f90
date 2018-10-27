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

!> \page Routines to call tenstream with optical properties from RRTM
!! The function `tenstream_rrtmg` provides an easy interface to
!! couple the TenStream solvers to a host model.
!!
!! * Tasks that have to be performed:
!!   - Radiative transfer needs to read a background profile that goes up till Top of the Atmosphere
!!   - Interpolate or merge the provided dynamics grid variables onto the the background profile
!!   - Compute optical properties with RRTMG
!!   - Solve the radiative transfer equation
!!   - And will return the fluxes on layer interfaces as well as the mean absorption in the layers
!!
!! Please carefully study the input and output parameters of subroutine ``tenstream_rrtmg()``
!!
!!

module m_plexrt_rrtmg

#include "petsc/finclude/petsc.h"
  use petsc

  use mpi, only : mpi_comm_rank
  use m_tenstr_parkind_sw, only: im => kind_im, rb => kind_rb
  use m_data_parameters, only : init_mpi_data_parameters, &
      iintegers, ireals, zero, one, i0, i1, i2, i9,         &
      mpiint, pi, default_str_len
  use m_adaptive_spectral_integration, only: need_new_solution
  use m_helper_functions, only : read_ascii_file_2d, gradient, meanvec, imp_bcast, &
      imp_allreduce_min, imp_allreduce_max, search_sorted_bisection, CHKERR, deg2rad, &
      reverse, itoa, angle_between_two_vec, norm, rad2deg, get_arg
  use m_tenstream_interpolation, only : interp_1d

  use m_plex_grid, only: TOAFACE, get_inward_face_normal, compute_face_geometry
  use m_plex_rt, only: t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, destroy_plexrt_solver, &
    plexrt_get_result

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, plkint, print_tenstr_atm, vert_integral_coeff
  use m_optprop_rrtmg, only: optprop_rrtm_lw, optprop_rrtm_sw
  use m_icon_plex_utils, only: Nz_Ncol_vec_to_celldm1

  use m_netcdfIO, only : ncwrite

  implicit none

  private
  public :: plexrt_rrtmg, destroy_plexrt_rrtmg

  logical,parameter :: ldebug=.True.
!  logical,parameter :: ldebug=.False.

contains

  subroutine plexrt_rrtmg(solver, atm, sundir,     &
      albedo_thermal, albedo_solar,                &
      lthermal, lsolar,                            &
      edir,edn,eup,abso,                           &
      opt_time, solar_albedo_2d, thermal_albedo_2d, &
      opt_solar_constant)

    class(t_plex_solver), allocatable, intent(inout)  :: solver ! solver type -- includes a dmplex and more
    type(t_tenstr_atm), intent(in)  :: atm                      ! atmosphere construct, which describes physical properties on layers and levels up till TOA
    real(ireals), intent(in)        :: sundir(:)                ! cartesian direction of sun rays
    real(ireals), intent(in)        :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions and compute new solutions only after threshold estimate is exceeded.
    ! If solar_albedo_2d is present, we use a 2D surface albedo, dim=(ncol,)
    ! TODO: introduce solar diffuse albedo.
    ! Or hack something together ala Ritter/Geleyn -> see:
    ! icon-lem/src/atm_phy_nwp/mo_nwp_rrtm_interface.f90 l.1371
    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:), thermal_albedo_2d(:), opt_solar_constant

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals), allocatable, dimension(:,:), intent(out) :: edir,edn,eup,abso          ! [nlyr(+1), ncol]

    ! ---------- end of API ----------------

    type(tIS) :: toa_ids

    integer(iintegers) :: k, Ncol, ke1

    integer(mpiint) :: myid, comm, ierr
    logical :: lrrtmg_only, lflg

    if(.not.allocated(solver)) call CHKERR(1_mpiint, 'solver has to be setup beforehand')
    if(.not.allocated(solver%plex)) call CHKERR(1_mpiint, 'Solver has to have a ready to go Plexgrid')

    call PetscObjectGetComm(solver%plex%dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    lrrtmg_only=.False. ! by default use normal tenstream solver
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
                             "-rrtmg_only" , lrrtmg_only , lflg , ierr) ;call CHKERR(ierr)

    ke1 = solver%plex%Nlay+1
    call CHKERR(int(ke1 - size(atm%plev,dim=1),mpiint), 'Vertical Size of atm and plex solver dont match')

    call DMGetStratumIS(solver%plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
    call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) then
      call print_tenstr_atm(atm,Ncol)
      print *,'m_plexrt_rrtmg sundir:', sundir, 'albedo th,sol',albedo_thermal, albedo_solar,'lth/lsol', lthermal, lsolar
    endif
    if(present(opt_time)) print *,'time', opt_time

    if(.not.allocated(edn )) allocate(edn (ke1, Ncol))
    if(.not.allocated(eup )) allocate(eup (ke1, Ncol))
    if(.not.allocated(abso)) allocate(abso(ke1-1, Ncol))
    edn  = zero
    eup  = zero
    abso = zero

    ! make sure that optprop vecs are allocated
    call allocate_optprop_vec(solver%plex%cell1_dm, solver%kabs)
    call allocate_optprop_vec(solver%plex%cell1_dm, solver%ksca)
    call allocate_optprop_vec(solver%plex%cell1_dm, solver%g   )
    call allocate_optprop_vec(solver%plex%srfc_boundary_dm, solver%albedo)

    if(lthermal) then
      call allocate_optprop_vec(solver%plex%cell1_dm, solver%plck)
      call allocate_optprop_vec(solver%plex%srfc_boundary_dm, solver%srfc_emission)

      call compute_thermal(comm, solver, atm, &
        Ncol, ke1, &
        albedo_thermal, &
        edn, eup, abso, &
        opt_time=opt_time, &
        thermal_albedo_2d=thermal_albedo_2d, &
        lrrtmg_only=lrrtmg_only)

      call dump_vec(edn(1:ke1-1,:), '-plexrt_dump_thermal_Edn_1_ke')
      call dump_vec(edn(2:ke1,:) ,  '-plexrt_dump_thermal_Edn_2_ke1')
      call dump_vec(eup(1:ke1-1,:), '-plexrt_dump_thermal_Eup_1_ke')
      call dump_vec(eup(2:ke1,:),   '-plexrt_dump_thermal_Eup_2_ke1')
      call dump_vec(abso, '-plexrt_dump_thermal_abso')
    endif

    if(lsolar) then
      if(.not.allocated(edir)) allocate(edir(ke1, Ncol))
      edir = zero
      call compute_solar(comm, solver, atm, &
        Ncol, ke1, &
        sundir, albedo_solar, &
        edir, edn, eup, abso, &
        opt_time=opt_time, &
        solar_albedo_2d=solar_albedo_2d, &
        lrrtmg_only=lrrtmg_only, &
        opt_solar_constant=opt_solar_constant)

      call dump_vec(edir(1:ke1-1,:),'-plexrt_dump_Edir_1_ke')
      call dump_vec(edir(2:ke1,:),  '-plexrt_dump_Edir_2_ke1')
    endif

    if(ldebug.and.myid.eq.0) then
      if(lsolar) then
        print *,'vert level    edir           edn              eup          abso'
        do k = 1, ke1-1
          print *,k, edir(k,i1), edn(k,i1), eup(k,i1), abso(k,i1)
        enddo
        print *,k, edir(ke1,i1), edn(ke1,i1), eup(ke1,i1)

        if(present(solar_albedo_2d)) then
          print *,'MinMax Solar Albedo', minval(solar_albedo_2d), maxval(solar_albedo_2d)
        endif
      else
        print *,'vert level    edn              eup          abso'
        do k = 1, ke1-1
          print *,k, edn(k,i1), eup(k,i1), abso(k,i1)
        enddo
        print *,k, edn(k,i1), eup(k,i1)
        if(present(thermal_albedo_2d)) then
          print *,'MinMax Thermal Albedo', minval(thermal_albedo_2d), maxval(thermal_albedo_2d)
        endif
      endif
    endif

    call dump_vec(edn(1:ke1-1,:), '-plexrt_dump_Edn_1_ke')
    call dump_vec(edn(2:ke1,:) ,  '-plexrt_dump_Edn_2_ke1')
    call dump_vec(eup(1:ke1-1,:), '-plexrt_dump_Eup_1_ke')
    call dump_vec(eup(2:ke1,:),   '-plexrt_dump_Eup_2_ke1')

    call dump_vec(abso, '-plexrt_dump_abso')

    if(allocated(atm%lwc )) call dump_vec(atm%lwc, '-plexrt_dump_lwc', lreverse=.True.)
    if(allocated(atm%iwc )) call dump_vec(atm%iwc, '-plexrt_dump_iwc', lreverse=.True.)
    if(allocated(atm%tlay)) call dump_vec(atm%tlay, '-plexrt_dump_temp', lreverse=.True.)

    contains
      subroutine dump_vec(inp, vecshow_string, lreverse)
        real(ireals), intent(in) :: inp(:,:) ! Nlay, Ncol
        character(len=*), intent(in) :: vecshow_string
        logical, intent(in), optional :: lreverse
        type(tVec) :: vec
        logical :: lflg, lrev

        lrev = get_arg(.False., lreverse)

        lflg = .False.
        call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, vecshow_string, lflg, ierr); call CHKERR(ierr)
        if(lflg) then
          if(.not.allocated(solver%plex%cell1_dm)) return
          call DMGetGlobalVector(solver%plex%cell1_dm, vec, ierr); call CHKERR(ierr)
          if(lrev) then
            call Nz_Ncol_vec_to_celldm1(solver%plex, reverse(inp), vec)
          else
            call Nz_Ncol_vec_to_celldm1(solver%plex, inp, vec)
          endif
          call PetscObjectSetName(vec, 'dump_vec'//trim(vecshow_string), ierr);call CHKERR(ierr)
          call PetscObjectViewFromOptions(vec, PETSC_NULL_VEC, trim(vecshow_string), ierr); call CHKERR(ierr)
          call DMRestoreGlobalVector(solver%plex%cell1_dm, vec, ierr); call CHKERR(ierr)
        endif
      end subroutine
  end subroutine

  subroutine compute_thermal(comm, solver, atm, Ncol, ke1, albedo, &
      edn, eup, abso, opt_time, thermal_albedo_2d, lrrtmg_only)

    use m_tenstr_rrlw_wvn, only : ngb, wavenum1, wavenum2
    use m_tenstr_parrrtm, only: ngptlw, nbndlw

    integer(mpiint), intent(in) :: comm
    class(t_plex_solver), allocatable, intent(inout)  :: solver
    type(t_tenstr_atm), intent(in), target :: atm
    integer(iintegers),intent(in)   :: Ncol, ke1

    real(ireals),intent(in) :: albedo

    real(ireals),intent(inout),dimension(:,:) :: edn, eup, abso

    real(ireals), optional, intent(in) :: opt_time, thermal_albedo_2d(:)
    logical, optional, intent(in) :: lrrtmg_only

    real(ireals),allocatable, dimension(:,:,:) :: tau, Bfrac          ! [nlyr, ncol, ngptlw]
    real(ireals),allocatable, dimension(:,:)   :: spec_abso           ! [nlyr(+1), ncol ]
    real(ireals),allocatable, dimension(:,:)   :: spec_edn, spec_eup  ! [nlyr(+1), ncol ]
    real(ireals),allocatable, dimension(:,:)   :: tmp, plck           ! [nlyr, ncol ]
    real(ireals),allocatable, dimension(:)     :: integral_coeff      ! [nlyr]

    real(ireals), pointer :: xalbedo(:), xsrfc_emission(:)

    integer(iintegers) :: i, ib, icol, k, current_ibnd
    logical :: need_any_new_solution

    integer(mpiint) :: ierr

    allocate(spec_edn (ke1, Ncol))
    allocate(spec_eup (ke1, Ncol))
    allocate(spec_abso(ke1-1, Ncol))

    need_any_new_solution=.False.
    do ib=1,ngptlw
      if(need_new_solution(comm, solver%solutions(500+ib), opt_time, solver%lenable_solutions_err_estimates)) &
        need_any_new_solution=.True.
    enddo
    if(.not.need_any_new_solution) then
      do ib=1,ngptlw
        call plexrt_get_result(solver, spec_edn, spec_eup, spec_abso, opt_solution_uid=500+ib)
        edn  = edn  + spec_edn
        eup  = eup  + spec_eup
        abso = abso + spec_abso
      enddo
      return
    endif

    ! Compute optical properties with RRTMG
    allocate(tau  (ke1-i1, Ncol, ngptlw))
    allocate(Bfrac(ke1-i1, Ncol, ngptlw))
    allocate(integral_coeff(ke1-i1))

    if(lrrtmg_only) then
        do i = 1, Ncol
          integral_coeff = vert_integral_coeff(atm%plev(1:ke1-1,i), atm%plev(2:ke1,i))

          call optprop_rrtm_lw(i1, ke1-i1, albedo, &
            atm%plev(:,i), atm%tlev(:,i), atm%tlay(:,i), &
            atm%h2o_lay(:,i), atm%o3_lay(:,i), atm%co2_lay(:,i), &
            atm%ch4_lay(:,i), atm%n2o_lay(:,i), atm%o2_lay(:,i), &
            atm%lwc(:,i)*integral_coeff, atm%reliq(:,i), &
            atm%iwc(:,i)*integral_coeff, atm%reice(:,i), &
            tau(:,i:i,:), Bfrac(:,i:i,:), &
            spec_eup(:,i:i), spec_edn(:,i:i), spec_abso(:,i:i))

          eup (:,i) = eup (:,i) + reverse(spec_eup (:,i))
          edn (:,i) = edn (:,i) + reverse(spec_edn (:,i))
          !abso(:,i) = abso(:,i) + reverse(spec_abso(:,i)) ! This would be in K/day
          abso(:,i) = abso(:,i) + reverse( ( &
              - spec_edn(1:ke1-1,i) + spec_edn(2:ke1,i) &
              + spec_eup(1:ke1-1,i) - spec_eup(2:ke1,i) ) / atm%dz(:,i) )
        enddo
      return
    else
        do i = 1, Ncol
          integral_coeff = vert_integral_coeff(atm%plev(1:ke1-1,i), atm%plev(2:ke1,i))

          call optprop_rrtm_lw(i1, ke1-i1, albedo, &
            atm%plev(:,i), atm%tlev(:,i), atm%tlay(:,i), &
            atm%h2o_lay(:,i), atm%o3_lay(:,i), atm%co2_lay(:,i), &
            atm%ch4_lay(:,i), atm%n2o_lay(:,i), atm%o2_lay(:,i), &
            atm%lwc(:,i)*integral_coeff, atm%reliq(:,i), &
            atm%iwc(:,i)*integral_coeff, atm%reice(:,i), &
            tau(:,i:i,:), Bfrac(:,i:i,:))
        enddo
    endif

    if(allocated(atm%opt_tau)) then
      if(.not.all( shape(atm%opt_tau) .eq. shape(tau(1:atm%d_ke,:,:)) )) then
        print *,'shape atm%opt_tau', shape(atm%opt_tau)
        print *,'shape tau', shape(tau)
        print *,'shape(tau(1:atm%d_ke)', shape(tau(1:atm%d_ke,:,:))
        call CHKERR(1_mpiint, 'Bad shape of atm%opt_tau; is: '// &
          itoa(shape(atm%opt_tau))//' should be '//itoa(shape(tau(1:atm%d_ke,:,:))))
      endif
        tau(1:atm%d_ke,:,:) = tau(1:atm%d_ke,:,:) + atm%opt_tau
    endif

    ! then fill in data, first the spectrally invariant properties
    call VecSet(solver%ksca, zero, ierr); call CHKERR(ierr)
    call VecSet(solver%g   , zero, ierr); call CHKERR(ierr)

    if(present(thermal_albedo_2d)) then
      call VecGetArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
      xalbedo = thermal_albedo_2d
      call VecRestoreArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
    else
      call VecSet(solver%albedo, albedo, ierr); call CHKERR(ierr)
    endif

    ! tmp space for transformations of cell properties
    allocate(tmp(ke1-1, Ncol))
    allocate(plck(ke1, Ncol)) ! actually have ke entries for layers plus one for surface. the srfc plck val is put at the end of the vert. axis

    current_ibnd = -1 ! current lw band
    do ib=1, ngptlw

      if(need_new_solution(comm, solver%solutions(500+ib), opt_time, solver%lenable_solutions_err_estimates)) then

        tmp = reverse(max(zero, tau(:,:,ib)) / atm%dz)
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%kabs)

        !Compute Plank Emission for nbndlw
        if(current_ibnd.eq.ngb(ib)) then ! still the same band, dont need to upgrade the plank emission
          continue
        else
          do icol=i1,Ncol
            do k=i1,ke1-1
              plck(k,icol) = plkint(real(wavenum1(ngb(ib))), real(wavenum2(ngb(ib))), real(atm%tlay(k,icol)))
            enddo
            plck(k,icol) = plkint(real(wavenum1(ngb(ib))), real(wavenum2(ngb(ib))), real(atm%tlev(1,icol)))
          enddo
          current_ibnd = ngb(ib)
        endif

        tmp = reverse(plck(1:ke1-1,:) * Bfrac(:,:,ib))
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%plck)

        call VecGetArrayF90(solver%srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)
        do icol = i1, Ncol
          xsrfc_emission(icol) = plck(ke1,icol) * Bfrac(1,icol,ib)
        enddo
        call VecRestoreArrayF90(solver%srfc_emission, xsrfc_emission, ierr); call CHKERR(ierr)

        call run_plex_rt_solver(solver, lthermal=.True., lsolar=.False., sundir=[zero, zero, one], &
          opt_solution_uid=500+ib, opt_solution_time=opt_time)

      endif
      call plexrt_get_result(solver, spec_edn, spec_eup, spec_abso, opt_solution_uid=500+ib)

      edn  = edn  + spec_edn
      eup  = eup  + spec_eup
      abso = abso + spec_abso

    enddo ! ib 1 -> nbndlw , i.e. spectral integration
  end subroutine compute_thermal

  subroutine compute_solar(comm, solver, atm, Ncol, ke1, &
      sundir, albedo, &
      edir, edn, eup, abso, opt_time, solar_albedo_2d, &
      lrrtmg_only, opt_solar_constant)

      use m_tenstr_parrrsw, only: ngptsw
      use m_tenstr_rrtmg_sw_spcvrt, only: tenstr_solsrc

    integer(mpiint), intent(in) :: comm
    class(t_plex_solver), allocatable, intent(inout)  :: solver
    type(t_tenstr_atm), intent(in), target :: atm
    integer(iintegers),intent(in)   :: Ncol, ke1

    real(ireals),intent(in) :: albedo
    real(ireals),intent(in) :: sundir(:)

    real(ireals),intent(inout),dimension(:,:) :: edir, edn, eup, abso

    real(ireals), intent(in), optional :: opt_time, solar_albedo_2d(:)
    logical, intent(in), optional :: lrrtmg_only
    real(ireals), intent(in), optional :: opt_solar_constant

    real(ireals),allocatable, dimension(:,:,:) :: tau, w0, g          ! [nlyr, ncol, ngptsw]
    real(ireals),allocatable, dimension(:,:)   :: spec_edir,spec_abso ! [nlyr(+1), ncol ]
    real(ireals),allocatable, dimension(:,:)   :: spec_edn, spec_eup  ! [nlyr(+1), ncol ]
    real(ireals),allocatable, dimension(:,:)   :: tmp  ! [nlyr, ncol ]
    real(ireals),allocatable, dimension(:)     :: integral_coeff  ! [nlyr]

    real(ireals), pointer :: xalbedo(:)

    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    type(tPetscSection) :: geomSection
    real(ireals) :: face_normal(3), theta0, Ag, rescaled_sundir(3)

    type(tIS) :: toa_ids
    integer(iintegers) :: i, ib, iface
    integer(iintegers), pointer :: cell_support(:), xitoa_faces(:)
    logical :: need_any_new_solution

    integer(mpiint) :: myid,ierr

    allocate(spec_edir(ke1, Ncol))
    allocate(spec_edn (ke1, Ncol))
    allocate(spec_eup (ke1, Ncol))
    allocate(spec_abso(ke1-1, Ncol))

    need_any_new_solution=.False.
    do ib=1,ngptsw
      if(need_new_solution(comm, solver%solutions(ib), opt_time, solver%lenable_solutions_err_estimates)) &
        need_any_new_solution=.True.
    enddo
    if(.not.need_any_new_solution) then
      do ib=1,ngptsw
        call plexrt_get_result(solver, spec_edn, spec_eup, spec_abso, spec_edir, opt_solution_uid=ib)
        edir = edir + spec_edir
        edn  = edn  + spec_edn
        eup  = eup  + spec_eup
        abso = abso + spec_abso
      enddo
      return
    endif

    call DMGetSection(solver%plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
    call VecGetArrayReadF90(solver%plex%geomVec, geoms, ierr); call CHKERR(ierr)
    call DMGetStratumIS(solver%plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
    call ISGetIndicesF90(toa_ids, xitoa_faces, ierr); call CHKERR(ierr)

    ! Compute optical properties with RRTMG
    allocate(tau(ke1-i1, Ncol, ngptsw))
    allocate(w0 (ke1-i1, Ncol, ngptsw))
    allocate(g  (ke1-i1, Ncol, ngptsw))
    allocate(integral_coeff(ke1-i1))

    if(lrrtmg_only) then
        do i = 1, Ncol
          iface = xitoa_faces(i)
          call DMPlexGetSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          call get_inward_face_normal(iface, cell_support(1), geomSection, geoms, face_normal)
          call DMPlexRestoreSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          theta0 = rad2deg(angle_between_two_vec(face_normal, sundir))
          Ag = albedo
          if(present(solar_albedo_2d)) Ag = solar_albedo_2d(i)


          if(theta0.lt.90._ireals) then
            integral_coeff = vert_integral_coeff(atm%plev(1:ke1-1,i), atm%plev(2:ke1,i))

            call optprop_rrtm_sw(i1, ke1-i1, &
              theta0, Ag, &
              atm%plev(:,i), atm%tlev(:,i), atm%tlay(:,i), &
              atm%h2o_lay(:,i), atm%o3_lay(:,i), atm%co2_lay(:,i), &
              atm%ch4_lay(:,i), atm%n2o_lay(:,i), atm%o2_lay(:,i), &
              atm%lwc(:,i)*integral_coeff, atm%reliq(:,i), &
              atm%iwc(:,i)*integral_coeff, atm%reice(:,i), &
              tau(:,i:i,:), w0(:,i:i,:), g(:,i:i,:), &
              spec_eup(:,i:i), spec_edn(:,i:i), spec_abso(:,i:i), &
              opt_solar_constant=opt_solar_constant)

            edir(:,i) = edir(:,i) + zero
            eup (:,i) = eup (:,i) + reverse(spec_eup (:,i))
            edn (:,i) = edn (:,i) + reverse(spec_edn (:,i))
            !abso(:,i) = abso(:,i) + reverse(spec_abso(:,i)) ! This would be in K/day
            abso(:,i) = abso(:,i) + reverse( ( &
              - spec_edn(1:ke1-1,i) + spec_edn(2:ke1,i) &
              + spec_eup(1:ke1-1,i) - spec_eup(2:ke1,i) ) / atm%dz(:,i) )
          else
            edir(:,i) = edir(:,i) + zero
            eup (:,i) = eup (:,i) + zero
            edn (:,i) = edn (:,i) + zero
            abso(:,i) = abso(:,i) + zero
          endif
        enddo
      return
    else
        do i = 1, Ncol
          iface = xitoa_faces(i)
          call DMPlexGetSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          call get_inward_face_normal(iface, cell_support(1), geomSection, geoms, face_normal)
          call DMPlexRestoreSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          theta0 = rad2deg(angle_between_two_vec(face_normal, sundir))
          Ag = albedo
          if(present(solar_albedo_2d)) Ag = solar_albedo_2d(i)

          integral_coeff = vert_integral_coeff(atm%plev(1:ke1-1,i), atm%plev(2:ke1,i))

          call optprop_rrtm_sw(i1, ke1-i1, &
            theta0, Ag, &
            atm%plev(:,i), atm%tlev(:,i), atm%tlay(:,i), &
            atm%h2o_lay(:,i), atm%o3_lay(:,i), atm%co2_lay(:,i), &
            atm%ch4_lay(:,i), atm%n2o_lay(:,i), atm%o2_lay(:,i), &
            atm%lwc(:,i)*integral_coeff, atm%reliq(:,i), &
            atm%iwc(:,i)*integral_coeff, atm%reice(:,i), &
            tau(:,i:i,:), w0(:,i:i,:), g(:,i:i,:), &
            opt_solar_constant=opt_solar_constant)
        enddo
    endif
    if(ldebug) then
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      if(myid.eq.0) print *,'DEBUG theta0', theta0, 'deg; 2d albedo?', present(solar_albedo_2d)
    endif
    w0 = min(one, max(zero, w0))
    call ISRestoreIndicesF90(toa_ids, xitoa_faces, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call VecGetArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
    if(present(solar_albedo_2d)) then
      xalbedo = solar_albedo_2d
    else
      xalbedo = albedo
    endif
    call VecRestoreArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)

    allocate(tmp(ke1-1, Ncol))

    do ib=1, ngptsw

      if(need_new_solution(comm, solver%solutions(ib), opt_time, solver%lenable_solutions_err_estimates)) then

        tmp = reverse(max(zero, tau(:,:,ib)) * (one - w0(:,:,ib)) / atm%dz)
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%kabs)

        tmp = reverse(max(zero, tau(:,:,ib)) * w0(:,:,ib) / atm%dz)
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%ksca)

        tmp    = reverse(min(one, max(zero, g(:,:,ib))))
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%g)

        if(present(opt_solar_constant)) then
          rescaled_sundir = sundir/norm(sundir) * tenstr_solsrc(ib) / sum(tenstr_solsrc) * opt_solar_constant
        else
          rescaled_sundir = sundir/norm(sundir) * tenstr_solsrc(ib)
        endif

        call run_plex_rt_solver(solver, lthermal=.False., lsolar=.True., &
          sundir=rescaled_sundir, &
          opt_solution_uid=ib, opt_solution_time=opt_time)
      endif
      call plexrt_get_result(solver, spec_edn, spec_eup, spec_abso, spec_edir, opt_solution_uid=ib)

      edir = edir + spec_edir
      edn  = edn  + spec_edn
      eup  = eup  + spec_eup
      abso = abso + spec_abso

    enddo ! ib 1 -> nbndsw , i.e. spectral integration
  end subroutine compute_solar

  subroutine allocate_optprop_vec(dm, vec, val)
    type(tDM), intent(in) :: dm
    type(tVec), allocatable, intent(inout) :: vec
    real(ireals), intent(in), optional :: val
    integer(mpiint) :: ierr

    if(.not.allocated(vec)) then
      allocate(vec)
      call DMCreateGlobalVector(dm, vec, ierr); call CHKERR(ierr)
    endif
    if(present(val)) then
      call VecSet(vec, val, ierr); call CHKERR(ierr)
    endif
  end subroutine

  subroutine destroy_plexrt_rrtmg(solver, lfinalizepetsc)
    class(t_plex_solver) :: solver
    logical, intent(in) :: lfinalizepetsc
    ! Tidy up the solver
    call destroy_plexrt_solver(solver, lfinalizepetsc=lfinalizepetsc)
  end subroutine
end module
