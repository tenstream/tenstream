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
      reverse, itoa, angle_between_two_vec, norm
  use m_tenstream_interpolation, only : interp_1d

  use m_plex_grid, only: TOAFACE, get_inward_face_normal
  use m_plex_rt, only: t_plex_solver, init_plex_rt_solver, run_plex_rt_solver, destroy_plexrt_solver, &
    plexrt_get_result

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, plkint, print_tenstr_atm
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
      opt_time, solar_albedo_2d)

    class(t_plex_solver), allocatable, intent(inout)  :: solver ! solver type -- includes a dmplex and more
    class(t_tenstr_atm), intent(in) :: atm                      ! atmosphere construct, which describes physical properties on layers and levels up till TOA
    real(ireals), intent(in)        :: sundir(:)                ! cartesian direction of sun rays
    real(ireals), intent(in)        :: albedo_solar, albedo_thermal ! broadband ground albedo for solar and thermal spectrum

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    logical, intent(in) :: lsolar, lthermal

    ! opt_time is the model time in seconds. If provided we will track the error growth of the solutions and compute new solutions only after threshold estimate is exceeded.
    ! If solar_albedo_2d is present, we use a 2D surface albedo, dim=(ncol,)
    ! TODO: introduce thermal 2d albedo and solar diffuse albedo.
    ! Or hack something together ala Ritter/Geleyn -> see:
    ! icon-lem/src/atm_phy_nwp/mo_nwp_rrtm_interface.f90 l.1371
    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:)

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

    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
                             "-rrtmg_only" , lrrtmg_only , lflg , ierr) ;call CHKERR(ierr)
    if(.not.lflg) lrrtmg_only=.False. ! by default use normal tenstream solver

    if(ldebug.and.myid.eq.0) &
      call print_tenstr_atm(atm)
    print *,'debug', sundir, albedo_thermal, albedo_solar, lthermal, lsolar
    if(present(opt_time)) print *,'time', opt_time

    ke1 = solver%plex%Nz
    call CHKERR(int(ke1 - size(atm%plev,dim=1),mpiint), 'Vertical Size of atm and plex solver dont match')

    call DMGetStratumIS(solver%plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
    call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)

    print *,'Ncol of toa_ids', Ncol

    if(.not.allocated(edn )) allocate(edn (ke1, Ncol))
    if(.not.allocated(eup )) allocate(eup (ke1, Ncol))
    if(.not.allocated(abso)) allocate(abso(ke1-1, Ncol))
    edn  = zero
    eup  = zero
    abso = zero

    if(lthermal) then
    !  call compute_thermal(solver, atm, ie, je, ke, ke1, &
    !    albedo_thermal, &
    !    edn, eup, abso, opt_time=opt_time, lrrtmg_only=lrrtmg_only)
    endif

    if(lsolar) then
      if(.not.allocated(edir)) allocate(edir(ke1, Ncol))
      edir = zero
      call compute_solar(solver, atm, Ncol, ke1, &
        sundir, albedo_solar, &
        edir, edn, eup, abso, opt_time=opt_time, solar_albedo_2d=solar_albedo_2d, &
        lrrtmg_only=lrrtmg_only)
    endif

    if(ldebug.and.myid.eq.0) then
      if(lsolar) then
        print *,'vert level    edir           edn              eup          abso'
        do k = 1, ke1-1
          print *,k, edir(k,i1), edn(k,i1), eup(k,i1), abso(k,i1)
        enddo
        print *,k, edir(ke1,i1), edn(ke1,i1), eup(ke1,i1)
      else
        print *,'vert level    edn              eup          abso'
        do k = 1, ke1-1
          print *,k, edn(k,i1), eup(k,i1), abso(k,i1)
        enddo
        print *,k, edn(ke1,i1), eup(ke1,i1)
      endif
    endif

  end subroutine

  subroutine compute_solar(solver, atm, Ncol, ke1, &
      sundir, albedo, &
      edir, edn, eup, abso, opt_time, solar_albedo_2d, lrrtmg_only)

      use m_tenstr_parrrsw, only: ngptsw
      use m_tenstr_rrtmg_sw_spcvrt, only: tenstr_solsrc

    type(t_plex_solver), allocatable, intent(inout)  :: solver
    type(t_tenstr_atm), intent(in), target :: atm
    integer(iintegers),intent(in)   :: Ncol, ke1

    real(ireals),intent(in) :: albedo
    real(ireals),intent(in) :: sundir(:)

    real(ireals),intent(inout),dimension(:,:) :: edir, edn, eup, abso

    real(ireals), optional, intent(in) :: opt_time, solar_albedo_2d(:)
    logical, optional, intent(in) :: lrrtmg_only

    real(ireals),allocatable, dimension(:,:,:) :: tau, w0, g          ! [nlyr, ncol, ngptsw]
    real(ireals),allocatable, dimension(:,:)   :: spec_edir,spec_abso ! [nlyr(+1), ncol ]
    real(ireals),allocatable, dimension(:,:)   :: spec_edn, spec_eup  ! [nlyr(+1), ncol ]
    real(ireals),allocatable, dimension(:,:)   :: tmp  ! [nlyr, ncol ]

    real(ireals), pointer :: xalbedo(:)

    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    type(tPetscSection) :: geomSection
    real(ireals) :: face_normal(3), theta0

    type(tIS) :: toa_ids
    integer(iintegers) :: i, ib, iface
    integer(iintegers), pointer :: cell_support(:), xitoa_faces(:)
    logical :: need_any_new_solution


    integer(mpiint) :: ierr

    allocate(spec_edir(ke1, Ncol))
    allocate(spec_edn (ke1, Ncol))
    allocate(spec_eup (ke1, Ncol))
    allocate(spec_abso(ke1-1, Ncol))

    need_any_new_solution=.False.
    do ib=1,ngptsw
      if(need_new_solution(solver%solutions(ib), opt_time, solver%lenable_solutions_err_estimates)) &
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

    if(lrrtmg_only) then
        do i = 1, Ncol
          iface = xitoa_faces(i)
          call DMPlexGetSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          call get_inward_face_normal(iface, cell_support(1), geomSection, geoms, face_normal)
          call DMPlexRestoreSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          theta0 = angle_between_two_vec(face_normal, sundir)

          call optprop_rrtm_sw(i1, ke1-i1, &
            theta0, albedo, &
            atm%plev(:,i), atm%tlev(:,i), atm%tlay(:,i), &
            atm%h2o_lay(:,i), atm%o3_lay(:,i), atm%co2_lay(:,i), &
            atm%ch4_lay(:,i), atm%n2o_lay(:,i), atm%o2_lay(:,i), &
            atm%lwc(:,i)*atm%dz(:,i), atm%reliq(:,i), &
            atm%iwc(:,i)*atm%dz(:,i), atm%reice(:,i), &
            tau(:,i:i,:), w0(:,i:i,:), g(:,i:i,:), &
            spec_eup(:,i:i), spec_edn(:,i:i), spec_abso(:,i:i))

          edir(:,i) = edir(:,i) + zero
          eup (:,i) = eup (:,i) + reverse(spec_eup (:,i))
          edn (:,i) = edn (:,i) + reverse(spec_edn (:,i))
          abso(:,i) = abso(:,i) + reverse(spec_abso(:,i))
        enddo
      return
    else
        do i = 1, Ncol
          iface = xitoa_faces(i)
          call DMPlexGetSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          call get_inward_face_normal(iface, cell_support(1), geomSection, geoms, face_normal)
          call DMPlexRestoreSupport(solver%plex%geom_dm, iface, cell_support, ierr); call CHKERR(ierr)
          theta0 = angle_between_two_vec(face_normal, sundir)

          call optprop_rrtm_sw(i1, ke1-i1, &
            theta0, albedo, &
            atm%plev(:,i), atm%tlev(:,i), atm%tlay(:,i), &
            atm%h2o_lay(:,i), atm%o3_lay(:,i), atm%co2_lay(:,i), &
            atm%ch4_lay(:,i), atm%n2o_lay(:,i), atm%o2_lay(:,i), &
            atm%lwc(:,i)*atm%dz(:,i), atm%reliq(:,i), &
            atm%iwc(:,i)*atm%dz(:,i), atm%reice(:,i), &
            tau(:,i:i,:), w0(:,i:i,:), g(:,i:i,:))
        enddo
    endif
    w0 = min(one, max(zero, w0))
    call ISRestoreIndicesF90(toa_ids, xitoa_faces, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(solver%plex%geomVec, geoms, ierr); call CHKERR(ierr)

    if(.not.allocated(solver%kabs)) then
      allocate(solver%kabs)
      call DMCreateGlobalVector(solver%plex%cell1_dm, solver%kabs, ierr); call CHKERR(ierr)
    endif
    if(.not.allocated(solver%ksca)) then
      allocate(solver%ksca)
      call DMCreateGlobalVector(solver%plex%cell1_dm, solver%ksca, ierr); call CHKERR(ierr)
    endif
    if(.not.allocated(solver%g)) then
      allocate(solver%g)
      call DMCreateGlobalVector(solver%plex%cell1_dm, solver%g, ierr); call CHKERR(ierr)
    endif

    call VecGetArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)
    if(present(solar_albedo_2d)) then
      xalbedo = solar_albedo_2d
    else
      xalbedo = albedo
    endif
    call VecRestoreArrayF90(solver%albedo, xalbedo, ierr); call CHKERR(ierr)

    allocate(tmp(ke1-1, Ncol))

    do ib=1, ngptsw

      if(need_new_solution(solver%solutions(ib), opt_time, solver%lenable_solutions_err_estimates)) then

        tmp = reverse(max(zero, tau(:,:,ib)) * (one - w0(:,:,ib)) / atm%dz)
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%kabs)

        tmp = reverse(max(zero, tau(:,:,ib)) * w0(:,:,ib) / atm%dz)
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%ksca)

        tmp    = reverse(min(one, max(zero, g(:,:,ib))))
        call Nz_Ncol_vec_to_celldm1(solver%plex, tmp, solver%g)

        call run_plex_rt_solver(solver, sundir/norm(sundir)*tenstr_solsrc(ib), &
          opt_solution_uid=ib, opt_solution_time=opt_time)
      endif
      call plexrt_get_result(solver, spec_edn, spec_eup, spec_abso, spec_edir, opt_solution_uid=ib)

      edir = edir + spec_edir
      edn  = edn  + spec_edn
      eup  = eup  + spec_eup
      abso = abso + spec_abso

    enddo ! ib 1 -> nbndsw , i.e. spectral integration
  end subroutine compute_solar

  subroutine destroy_plexrt_rrtmg(solver, lfinalizepetsc)
    class(t_plex_solver) :: solver
    logical, intent(in) :: lfinalizepetsc
    ! Tidy up the solver
    call destroy_plexrt_solver(solver, lfinalizepetsc=lfinalizepetsc)
  end subroutine
end module
