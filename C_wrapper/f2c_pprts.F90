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

module m_f2c_pprts
#include "petsc/finclude/petsc.h"
      use petsc

      use iso_c_binding

      use m_data_parameters, only : init_mpi_data_parameters, &
        iintegers, ireals, mpiint, default_str_len, &
        i0, i1, zero, one

      use m_tenstream_options, only: read_commandline_options

      use m_helper_functions, only: imp_bcast, meanval, CHKERR, &
        resize_arr, spherical_2_cartesian, itoa

      use m_pprts_base, only : t_solver, t_solver_3_10, t_solver_3_6, t_solver_1_2, t_coord, &
        t_solver_8_10, t_solver_3_16, t_solver_8_16, t_solver_8_18

      use m_pprts, only : init_pprts, set_global_optical_properties, solve_pprts, destroy_pprts, &
        pprts_get_result_toZero

      use m_petsc_helpers, only : petscVecToF90, petscGlobalVecToZero, f90VecToPetsc

      use m_icon_plex_utils, only: create_2d_regular_plex, dmplex_2D_to_3D, &
        rank0_f90vec_to_plex, dmplex_gvec_from_f90_array, plex_gvec_tozero
      use m_plex_grid, only: t_plexgrid, setup_plexgrid, create_plex_section

      use m_plex_rt_base, only: allocate_plexrt_solver_from_commandline, &
        t_plex_solver, t_plex_solver_rectilinear_5_8

      use m_plex_rt, only: &
        init_plex_rt_solver, run_plex_rt_solver, plexrt_get_result, &
        destroy_plexrt_solver

      use m_netcdfio, only: ncwrite

      implicit none

      private
      public :: pprts_f2c_init, pprts_f2c_set_global_optical_properties, &
        pprts_f2c_solve, pprts_f2c_get_result, pprts_f2c_destroy

      class(t_solver), allocatable :: pprts_solver
      class(t_plex_solver), allocatable :: plex_solver
      type(tDM) :: dm2d, dm2d_dist
      type(tPetscSF) :: migration_sf
      real(ireals) :: sundir(3)

      logical, parameter :: ldebug=.False.

#include "f2c_solver_ids.h"

contains

  ! initialize pprts environment
  ! all nodes in communicator have to call this
  ! but only the zeroth node has to have meaningful values for the arguments except the communicator
  ! all but hhl is overwritten on nonzero nodes
  subroutine pprts_f2c_init(comm, solver_id, Nz,Nx,Ny,dx,dy,hhl, phi0, theta0, collapseindex) bind(C)
    integer(c_int), value :: comm
    integer(c_int),intent(inout) :: solver_id, Nx, Ny, Nz
    real(c_double),intent(inout) :: dx,dy
    real(c_float), intent(in),dimension(Nz+1) :: hhl
    real(c_float), intent(inout) :: phi0,theta0
    integer(c_int),intent(inout) :: collapseindex

    integer(iintegers) :: osolver_id, oNx, oNy, oNz, ocollapseindex
    real(ireals) :: odx,ody,ophi0,otheta0
    real(ireals),allocatable :: ohhl(:)

    real(ireals),allocatable :: odz(:)
    integer(iintegers) :: k
    integer(mpiint) :: myid, ierr

    logical,save :: initialized=.False.

    if(initialized) then
      ierr = 0
      if(allocated(pprts_solver)) then
        select type(pprts_solver)
        class is (t_solver_1_2)
          if(solver_id.ne.SOLVER_ID_PPRTS_1_2) ierr = 1
        class is (t_solver_3_6)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_6) ierr = 1
        class is (t_solver_3_10)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_10) ierr = 1
        class is (t_solver_8_10)
          if(solver_id.ne.SOLVER_ID_PPRTS_8_10) ierr = 1
        class is (t_solver_3_16)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_16) ierr = 1
        class is (t_solver_8_16)
          if(solver_id.ne.SOLVER_ID_PPRTS_8_16) ierr = 1
        class is (t_solver_8_18)
          if(solver_id.ne.SOLVER_ID_PPRTS_8_18) ierr = 1
        end select
      endif
      if(allocated(plex_solver)) then
        select type(plex_solver)
        class is (t_plex_solver_rectilinear_5_8)
          if(solver_id.ne.SOLVER_ID_PLEXRT_RECTILINEAR_5_8) ierr = 1
        end select
      endif
      if(ierr.ne.0) call CHKERR(ierr, 'seems you changed the solver type id in between calls...'// &
        'you must destroy the solver first before you use a different kind of solver')
      return
    endif

    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(myid.eq.0) then
      osolver_id = int(solver_id, kind=iintegers)
      oNx     = int(Nx, kind=iintegers)
      oNy     = int(Ny, kind=iintegers)
      oNz     = int(Nz, kind=iintegers)
      odx     = real(dx    ,kind=ireals)
      ody     = real(dy    ,kind=ireals)
      ophi0   = real(phi0  ,kind=ireals)
      otheta0 = real(theta0,kind=ireals)
      ocollapseindex = int(collapseindex, kind=iintegers)

      allocate( ohhl(size(hhl)) )
      ohhl = hhl
    endif

    call imp_bcast(comm, osolver_id, 0_mpiint)
    call imp_bcast(comm, oNx    ,0_mpiint)
    call imp_bcast(comm, oNy    ,0_mpiint)
    call imp_bcast(comm, oNz    ,0_mpiint)
    call imp_bcast(comm, odx    ,0_mpiint)
    call imp_bcast(comm, ody    ,0_mpiint)
    call imp_bcast(comm, ophi0  ,0_mpiint)
    call imp_bcast(comm, otheta0,0_mpiint)
    call imp_bcast(comm, ocollapseindex,0_mpiint)

    call imp_bcast(comm, ohhl,0_mpiint)

    ! and overwrite input values to propagate values back to caller...
    solver_id = int(osolver_id, kind=c_int)
    Nx     = int(oNx, kind=c_int)
    Ny     = int(oNy, kind=c_int)
    Nz     = int(oNz, kind=c_int)
    dx     = real(odx    , kind=c_double)
    dy     = real(ody    , kind=c_double)
    phi0   = real(ophi0  , kind=c_float)
    theta0 = real(otheta0, kind=c_float)
    collapseindex = int(ocollapseindex, kind=c_int)

    ! Now every process has the correct values
    if(ldebug) print *,myid,'Initializing pprts environment from C Language :: solver_id', solver_id
    if(ldebug) print *,myid,'Initializing pprts environment from C Language :: domainshape',oNx,oNy,oNz,'::',shape(ohhl)

    allocate(odz(oNz))
    do k=1,Nz
      odz(k) = ohhl(k) - ohhl(k+1)
    enddo

    ierr=1
    select case(solver_id)
    case(SOLVER_ID_PPRTS_1_2)
      allocate(t_solver_1_2::pprts_solver); ierr=0
      call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, pprts_solver, dz1d=odz)
    case(SOLVER_ID_PPRTS_3_6)
      allocate(t_solver_3_6::pprts_solver); ierr=0
      call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, pprts_solver, dz1d=odz)
    case(SOLVER_ID_PPRTS_3_10)
      allocate(t_solver_3_10::pprts_solver); ierr=0
      call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, pprts_solver, dz1d=odz)
    case(SOLVER_ID_PPRTS_8_10)
      allocate(t_solver_8_10::pprts_solver); ierr=0
      call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, pprts_solver, dz1d=odz)
    case(SOLVER_ID_PPRTS_3_16)
      allocate(t_solver_3_16::pprts_solver); ierr=0
      call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, pprts_solver, dz1d=odz)
    case(SOLVER_ID_PPRTS_8_16)
      allocate(t_solver_8_16::pprts_solver); ierr=0
      call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, pprts_solver, dz1d=odz)
    case(SOLVER_ID_PPRTS_8_18)
      allocate(t_solver_8_18::pprts_solver); ierr=0
      call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, pprts_solver, dz1d=odz)
    case(SOLVER_ID_PLEXRT_RECTILINEAR_5_8); ierr=0
      call plexrt_f2c_init(comm, osolver_id, oNx, oNy, oNz+1, odx, ohhl, ophi0, otheta0)
    end select
    call CHKERR(ierr, 'Could not find a suitable solver to allocate for solver_id: '//itoa(solver_id))

    initialized=.True.
  end subroutine

  subroutine pprts_f2c_set_global_optical_properties(Nz, Nx, Ny, albedo, kabs, ksca, g, planck) bind(c)
    integer(c_int), value :: Nx,Ny,Nz
    real(c_float),intent(in) :: albedo
    real(c_float),intent(in),dimension(Nz  ,Nx,Ny) :: kabs, ksca, g
    real(c_float),intent(in),dimension(Nz+1,Nx,Ny) :: planck

    real(ireals) :: oalbedo
    real(ireals),allocatable,dimension(:,:,:) :: okabs, oksca, og, oplanck

    if(allocated(pprts_solver)) then
      if(pprts_solver%myid.eq.0) then
        oalbedo = real(albedo, kind=ireals)
        allocate( okabs  (Nz  ,Nx,Ny) ); okabs   = real(kabs, ireals)
        allocate( oksca  (Nz  ,Nx,Ny) ); oksca   = real(ksca, ireals)
        allocate( og     (Nz  ,Nx,Ny) ); og      = real(g, ireals)
        allocate( oplanck(Nz+1,Nx,Ny) ); oplanck = real(planck, ireals)

        if(any(oplanck.gt.zero)) then
          call set_global_optical_properties(pprts_solver, oalbedo, okabs, oksca, og, oplanck)
        else
          call set_global_optical_properties(pprts_solver, oalbedo, okabs, oksca, og)
        endif

        if(ldebug) print *,'mean kabs  ',meanval(okabs)
        if(ldebug) print *,'mean ksca  ',meanval(oksca)
        if(ldebug) print *,'mean g     ',meanval(og)
        if(ldebug) print *,'mean planck',meanval(oplanck)
      else !slave
        call set_global_optical_properties(pprts_solver)
      endif
    endif

    if(allocated(plex_solver)) then
      call plexrt_f2c_set_global_optprop(Nz, Nx, Ny, albedo, kabs, ksca, g, planck)
    endif
  end subroutine

  subroutine pprts_f2c_solve(comm, edirTOA) bind(c)
    ! solve pprts equations
    ! optical properties have had to be set and environment had to be initialized
    ! incoming solar radiation need only be set by zeroth node
    integer(c_int), value :: comm
    real(c_float), value :: edirTOA
    real(ireals) :: oedirTOA
    logical :: lthermal
    integer(mpiint) :: myid, ierr
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(myid.eq.0) oedirTOA = real(edirTOA, ireals)
    call imp_bcast(comm, oedirTOA, 0_mpiint)

    if(allocated(pprts_solver)) then
      call solve_pprts(pprts_solver, oedirTOA)
    endif
    if(allocated(plex_solver)) then
      lthermal = allocated(plex_solver%plck)

      if(ldebug) &
        print *,'pprts_f2c_solve; '// &
                'lthermal', lthermal, &
                'lsolar', oedirTOA.gt.0, &
                'sundir=', sundir
      call run_plex_rt_solver(plex_solver, lthermal=lthermal, lsolar=oedirTOA.gt.0, sundir=sundir)
    endif
  end subroutine

  subroutine pprts_f2c_get_result(Nz,Nx,Ny, res_edn, res_eup, res_abso, res_edir) bind(c)
    ! after solving equations -- retrieve the results for edir,edn,eup and absorption
    ! only zeroth node gets the results back.

    integer(c_int), value :: Nx,Ny,Nz
    real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edn
    real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_eup
    real(c_float),intent(out),dimension(Nz  ,Nx,Ny) :: res_abso
    real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edir

    real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso
    integer(mpiint) :: comm, myid, ierr
    logical :: lflg
    character(len=default_str_len) :: dump_fname(2)

    if(allocated(pprts_solver)) then
      call pprts_get_result_toZero(pprts_solver,redn,reup,rabso,redir)
      comm = pprts_solver%comm
      if(pprts_solver%myid.eq.0) then
        res_edn  = real(redn (:, 1:Nx, 1:Ny), kind=c_float)
        res_eup  = real(reup (:, 1:Nx, 1:Ny), kind=c_float)
        res_abso = real(rabso(:, 1:Nx, 1:Ny), kind=c_float)
        res_edir = real(redir(:, 1:Nx, 1:Ny), kind=c_float)
      endif
    endif

    if(allocated(plex_solver)) then
      call PetscObjectGetComm(plex_solver%plex%dm, comm, ierr); call CHKERR(ierr)
      call plex_solver_retrieve_and_average_results(comm, res_edn, res_eup, res_abso, res_edir)
    endif

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(myid.eq.0) then
      if(ldebug) print *,'pprts_f2c_get_result result_edir first column', res_edir(:,1,1)
      if(ldebug) print *,'pprts_f2c_get_result result_edn  first column', res_edn (:,1,1)
      if(ldebug) print *,'pprts_f2c_get_result result_eup  first column', res_eup (:,1,1)
      if(ldebug) print *,'pprts_f2c_get_result rabso first column', res_abso(:,1,1)
      call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-f2c_dump_result", &
        dump_fname(1), lflg , ierr) ;call CHKERR(ierr)
      if(lflg) then
        dump_fname(2) = "edir"; call ncwrite(dump_fname, res_edir, ierr); call CHKERR(ierr)
        dump_fname(2) = "edn" ; call ncwrite(dump_fname, res_edn, ierr); call CHKERR(ierr)
        dump_fname(2) = "eup" ; call ncwrite(dump_fname, res_eup, ierr); call CHKERR(ierr)
        dump_fname(2) = "abso"; call ncwrite(dump_fname, res_abso, ierr); call CHKERR(ierr)
      endif
    endif

    contains
      subroutine plex_solver_retrieve_and_average_results(comm, res_edn, res_eup, res_abso, res_edir)
        integer(mpiint), intent(in) :: comm
        real(c_float), intent(out), dimension(:,:,:) :: res_edn, res_eup, res_abso, res_edir
        real(ireals),allocatable,dimension(:,:) :: redir,redn,reup,rabso

        call plexrt_get_result(plex_solver, redn, reup, rabso, redir)

        call transfer_one_var(comm, redn, res_edn)
        call transfer_one_var(comm, reup, res_eup)
        call transfer_one_var(comm, rabso, res_abso)
        call transfer_one_var(comm, redir, res_edir)
      end subroutine
      subroutine transfer_one_var(comm, var, var0)
        integer(mpiint), intent(in) :: comm
        real(ireals), intent(in) :: var(:,:) ! output from plexrt dim (Nz, Ncol), i.e. Ncol==2*Nx*Ny
        real(c_float), intent(out) :: var0(:,:,:) ! output on rank0 dim (Nz, Nx, Ny)
        type(tPetscSection) :: flxSection, r0flxSection
        type(tVec) :: v_var
        type(tVec) :: r0var
        real(ireals), pointer :: xv(:), xxv(:,:,:)
        integer(mpiint) :: myid, ierr

        call create_plex_section(dm2d_dist, 'face_section', i1, &
          [i0], [size(var,dim=1,kind=iintegers)], [i0], [i0], flxSection)

        call dmplex_gVec_from_f90_array(comm, var, v_var)
        call plex_gVec_toZero(dm2d_dist, migration_sf, flxSection, v_var, &
          r0flxSection, r0var)

        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
        if(myid.eq.0) then
          call VecGetArrayF90(r0var, xv,ierr); call CHKERR(ierr)
          call CHKERR(int(size(xv)/2-size(var0), mpiint), 'Global array sizes do not match, expected input size ('// &
            itoa(shape(var0))//') to be half the size of the plexrt mesh result: ('//itoa(shape(xv))//')')
          xxv(1:size(var,dim=1), 1:2*Nx, 1:Ny) => xv
          var0 = real( xxv(:, 1:size(xxv,2):2, :) + xxv(:, 2:size(xxv,2):2, :), c_float) / 2
          nullify(xxv)
          call VecRestoreArrayF90(r0var, xv,ierr); call CHKERR(ierr)
        endif
        call PetscSectionDestroy(flxSection, ierr); call CHKERR(ierr)
        call PetscSectionDestroy(r0flxSection, ierr); call CHKERR(ierr)
        call VecDestroy(v_var, ierr); call CHKERR(ierr)
        call VecDestroy(r0var, ierr); call CHKERR(ierr)
      end subroutine

  end subroutine

  subroutine pprts_f2c_destroy(ifinalizepetsc) bind(c)
    integer(c_int), value, intent(in) :: ifinalizepetsc
    logical :: lfinalizepetsc
    integer(mpiint) :: ierr
    lfinalizepetsc = ifinalizepetsc.ne.0
    if(allocated(pprts_solver)) then
      call destroy_pprts(pprts_solver, lfinalizepetsc=lfinalizepetsc)
      deallocate(pprts_solver)
    endif
    if(allocated(plex_solver)) then
      call PetscSFDestroy(migration_sf, ierr); call CHKERR(ierr)
      call DMDestroy(dm2d, ierr); call CHKERR(ierr)
      call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)
      call destroy_plexrt_solver(plex_solver, lfinalizepetsc=lfinalizepetsc)
      deallocate(plex_solver)
    endif
  end subroutine


  ! --------------------- Start of PLEXRT Routines -------------------
  subroutine plexrt_f2c_init(comm, solver_id, &
      Nx_global, Ny_global, Nlev, &
      dx, hhl, phi0, theta0 )
    integer(mpiint), intent(in) :: comm
    integer(iintegers), intent(in) :: solver_id, Nx_global, Ny_global, Nlev
    real(ireals), intent(in) :: dx, hhl(:), phi0, theta0

    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)
    integer(iintegers) :: fStart, fEnd, Ncol
    type(tDM) :: dm3d
    integer(mpiint) :: myid, ierr

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call create_2d_regular_plex(comm, Nx_global+1, Ny_global+1, dm2d, dm2d_dist, &
      opt_migration_sf=migration_sf, opt_dx=dx)

    call DMPlexGetHeightStratum(dm2d, i0, fStart, fEnd, ierr); call CHKERR(ierr)
    Ncol = fEnd - fStart

    if(ldebug) print *,myid,'plexrt_f2c_init: Global Domain sizes are:', Ncol, Nlev

    call DMPlexGetHeightStratum(dm2d_dist, i0, fStart, fEnd, ierr); call CHKERR(ierr)
    Ncol = fEnd - fStart

    if(ldebug) print *,myid,'plexrt_f2c_init: Local Domain sizes are:', Ncol, Nlev

    call dmplex_2D_to_3D(dm2d_dist, Nlev, hhl, dm3d, zindex, lpolar_coords=.False.)

    call setup_plexgrid(dm2d_dist, dm3d, Nlev-1,zindex, plex, hhl)
    call DMDestroy(dm3d, ierr); call CHKERR(ierr)
    deallocate(zindex)

    select case(solver_id)
    case(SOLVER_ID_PLEXRT_RECTILINEAR_5_8)
      call allocate_plexrt_solver_from_commandline(plex_solver, 'rectilinear_5_8')
    case default
      call CHKERR(1_mpiint, 'unknown solver id::'//itoa(solver_id))
    end select

    call init_plex_rt_solver(plex, plex_solver)

    sundir = spherical_2_cartesian(phi0,theta0)
    if(ldebug) print *,'f2c_pprts::sundir', sundir, '(', phi0,theta0, ')'
  end subroutine

  subroutine plexrt_f2c_set_global_optprop(Nz, Nx, Ny, albedo, kabs, ksca, g, planck) bind(c)
    integer(c_int), value :: Nx,Ny,Nz
    real(c_float),intent(in) :: albedo
    real(c_float),intent(in),dimension(Nz  ,Nx,Ny) :: kabs, ksca, g
    real(c_float),intent(in),dimension(Nz+1,Nx,Ny) :: planck

    real(ireals) :: oalbedo
    real(ireals),allocatable,target,dimension(:,:,:) :: work
    integer(mpiint) :: comm, myid, ierr
    logical :: lthermal

    call PetscObjectGetComm(plex_solver%plex%dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    if(myid.eq.0) oalbedo = real(albedo, ireals)
    call imp_bcast(comm, oalbedo, 0_mpiint); call CHKERR(ierr)

    if(.not.allocated(plex_solver%albedo)) then
      allocate(plex_solver%albedo)
      call DMCreateGlobalVector(plex_solver%plex%srfc_boundary_dm, plex_solver%albedo, ierr); call CHKERR(ierr)
    endif
    call VecSet(plex_solver%albedo, oalbedo, ierr); call CHKERR(ierr)

    if(myid.eq.0) allocate( work(Nz,2*Nx,Ny) )
    call propagate_from_Zero_to_solver_optprop(kabs, work, plex_solver%kabs)
    call propagate_from_Zero_to_solver_optprop(ksca, work, plex_solver%ksca)
    call propagate_from_Zero_to_solver_optprop(g   , work, plex_solver%g   )

    if(myid.eq.0) lthermal = any(planck.gt.zero)
    if(ldebug.and.myid.eq.0) print *,'plexrt_f2c_set_global_optprop lthermal', lthermal
    call imp_bcast(comm, lthermal, 0_mpiint)

    if(lthermal) then
      if(myid.eq.0) then
        deallocate(work)
        allocate( work(Nz+1,2*Nx,Ny) )
      endif
      call propagate_from_Zero_to_solver_optprop(planck, work, plex_solver%plck)
    endif

  contains
    subroutine propagate_from_Zero_to_solver_optprop(arr, work, solvervec)
      real(c_float), intent(in) :: arr(:,:,:) ! has global cartesian mesh dimensions
      real(ireals), intent(inout), contiguous, target :: work(:,:,:) ! has global rectilinear mesh dimensions
      type(tVec), allocatable, intent(inout) :: solvervec ! e.g. solver%kabs
      real(ireals), pointer :: col_arr(:,:)
      type(tPetscSection) :: parCellSection

      if(myid.eq.0) then
        work(:, 1:size(work,2):2, :) = arr
        work(:, 2:size(work,2):2, :) = arr
        col_arr(1:size(arr,1), 1:size(work,2)*size(work,3)) => work
      endif

      if(.not.allocated(solvervec)) allocate(solvervec)
      call rank0_f90vec_to_plex(dm2d, dm2d_dist, migration_sf, &
        col_arr, parCellSection, solvervec)
      call PetscSectionDestroy(parCellSection, ierr); call CHKERR(ierr)
    end subroutine
  end subroutine

end module
