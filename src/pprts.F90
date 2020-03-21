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

module m_pprts

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, irealLUT,     &
    init_mpi_data_parameters, mpiint,                            &
    zero, one, nil, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi, &
    default_str_len

  use m_helper_functions, only : CHKERR, CHKWARN, deg2rad, rad2deg, imp_allreduce_min, &
    imp_bcast, imp_allreduce_max, delta_scale, mpi_logical_and, meanval, get_arg, approx, &
    inc, ltoa, cstr, itoa, ftoa, imp_allreduce_mean

  use m_twostream, only: delta_eddington_twostream, adding_delta_eddington_twostream
  use m_schwarzschild, only: schwarzschild, B_eff
  use m_optprop, only: t_optprop, t_optprop_1_2, t_optprop_3_6, t_optprop_3_10, &
    t_optprop_8_10, t_optprop_3_16, t_optprop_8_16, t_optprop_8_18
  use m_eddington, only : eddington_coeff_zdun

  use m_tenstream_options, only : read_commandline_options, ltwostr, luse_eddington, twostr_ratio, &
    options_max_solution_err, options_max_solution_time, ltwostr_only, luse_twostr_guess,        &
    options_phi, lforce_phi, options_theta, lforce_theta, &
    lcalc_nca, lskip_thermal, lschwarzschild, ltopography, &
    lmcrts

  use m_petsc_helpers, only : petscGlobalVecToZero, scatterZerotoPetscGlobal, &
    petscVecToF90, f90VecToPetsc, getVecPointer, restoreVecPointer, hegedus_trick

  use m_mcrts_dmda, only : solve_mcrts

  use m_pprts_base, only : t_solver, t_solver_1_2, t_solver_3_6, t_solver_3_10, &
    t_solver_8_10, t_solver_3_16, t_solver_8_16, t_solver_8_18, &
    t_coord, t_suninfo, t_atmosphere, compute_gradient, atmk, &
    t_state_container, prepare_solution, destroy_solution, &
    t_dof, t_solver_log_events, setup_log_events

  implicit none
  private

  public :: init_pprts, &
            set_optical_properties, set_global_optical_properties, &
            solve_pprts, set_angles, destroy_pprts, pprts_get_result, &
            pprts_get_result_toZero, gather_all_toZero, scale_flx

  logical,parameter :: ldebug=.False.
  logical,parameter :: lcycle_dir=.True.
  logical,parameter :: lprealloc=.True.

  integer(iintegers),parameter :: minimal_dimension=3 ! this is the minimum number of gridpoints in x or y direction

  integer(mpiint) :: ierr

  contains

  !> @brief Main routine to setup PPRTS solver
  !> @details This will setup the PETSc DMDA grid and set other grid information, needed for the pprts
  !> \n Nx, Ny Nz are either global domain size or have to be local sizes if present(nxproc,nyproc)
  !> \n where nxproc and nyproc then are the number of pixel per rank for all ranks -- i.e. sum(nxproc) != Nx_global
  subroutine init_pprts(icomm, Nz,Nx,Ny, dx,dy, phi0, theta0, solver, dz1d, dz3d, nxproc, nyproc, collapseindex, solvername)
    MPI_Comm, intent(in)          :: icomm         !< @param MPI_Communicator for this solver
    integer(iintegers),intent(in) :: Nz            !< @param[in] Nz     Nz is the number of layers and Nz+1 would be the number of levels
    integer(iintegers),intent(in) :: Nx            !< @param[in] Nx     number of boxes in x-direction
    integer(iintegers),intent(in) :: Ny            !< @param[in] Ny     number of boxes in y-direction
    real(ireals),intent(in)       :: dx            !< @param[in] dx     physical size of grid in [m]
    real(ireals),intent(in)       :: dy            !< @param[in] dy     physical size of grid in [m]
    real(ireals),intent(in)       :: phi0          !< @param[in] phi0   solar azmiuth and zenith angle
    real(ireals),intent(in)       :: theta0        !< @param[in] theta0 solar azmiuth and zenith angle

    class(t_solver), intent(inout)         :: solver         !< @param[inout] solver
    real(ireals),optional,intent(in)       :: dz1d(:)        !< @param[in]    dz1d    if given, dz1d is used everywhere on the rank
    real(ireals),optional,intent(in)       :: dz3d(:,:,:)    !< @param[in]    dz3d    if given, dz3d has to be local domain size, cannot have global shape
    integer(iintegers),optional,intent(in) :: nxproc(:)      !< @param[in]    nxproc  if given, Nx has to be the local size, dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    integer(iintegers),optional,intent(in) :: nyproc(:)      !< @param[in]    nyproc  if given, Ny has to be the local size, dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    integer(iintegers),optional,intent(in) :: collapseindex  !< @param[in]    collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them
    character(len=*), optional, intent(in) :: solvername     !< @param[in] primarily for logging purposes, name will be prefix to logging stages

    integer(iintegers) :: k,i,j
    logical :: lpetsc_is_initialized

    if(.not.solver%linitialized) then

      solver%difftop%area_divider = 1
      solver%diffside%area_divider = 1

      solver%dirtop%area_divider = 1
      solver%dirside%area_divider = 1

      select type(solver)
        class is (t_solver_1_2)
          allocate(solver%difftop%is_inward(2))
          solver%difftop%is_inward = [.False.,.True.]

          allocate(solver%diffside%is_inward(0))

          allocate(solver%dirtop%is_inward(1))
          solver%dirtop%is_inward = [.True.]

          allocate(solver%dirside%is_inward(0))

        class is (t_solver_3_6)

          allocate(solver%difftop%is_inward(2))
          solver%difftop%is_inward = [.False.,.True.]

          allocate(solver%diffside%is_inward(2))
          solver%diffside%is_inward = [.False.,.True.]

          allocate(solver%dirtop%is_inward(1))
          solver%dirtop%is_inward = [.True.]

          allocate(solver%dirside%is_inward(1))
          solver%dirside%is_inward = [.True.]

        class is (t_solver_3_10)

          allocate(solver%difftop%is_inward(2))
          solver%difftop%is_inward = [.False.,.True.]

          allocate(solver%diffside%is_inward(4))
          solver%diffside%is_inward = [.False.,.True.,.False.,.True.]

          allocate(solver%dirtop%is_inward(1))
          solver%dirtop%is_inward = [.True.]

          allocate(solver%dirside%is_inward(1))
          solver%dirside%is_inward = [.True.]

        class is (t_solver_8_10)

          allocate(solver%difftop%is_inward(2))
          solver%difftop%is_inward = [.False.,.True.]

          allocate(solver%diffside%is_inward(4))
          solver%diffside%is_inward = [.False.,.True.,.False.,.True.]

          allocate(solver%dirtop%is_inward(4))
          solver%dirtop%is_inward = .True.
          solver%dirtop%area_divider = 4

          allocate(solver%dirside%is_inward(2))
          solver%dirside%is_inward = .True.
          solver%dirside%area_divider = 2

        class is (t_solver_3_16)

          allocate(solver%difftop%is_inward(8), source= &
            [.False.,.True.,.False.,.True.,.False.,.True.,.False.,.True.])

          allocate(solver%diffside%is_inward(4))
          solver%diffside%is_inward = [.False.,.True.,.False.,.True.]

          allocate(solver%dirtop%is_inward(1))
          solver%dirtop%is_inward = .True.

          allocate(solver%dirside%is_inward(1))
          solver%dirside%is_inward = .True.

        class is (t_solver_8_16)

          allocate(solver%difftop%is_inward(8), source= &
            [.False.,.True.,.False.,.True.,.False.,.True.,.False.,.True.])

          allocate(solver%diffside%is_inward(4), source=[.False.,.True.,.False.,.True.])

          allocate(solver%dirtop%is_inward(4), source=.True.)
          solver%dirtop%area_divider = 4

          allocate(solver%dirside%is_inward(2), source=.True.)
          solver%dirside%area_divider = 2

        class is (t_solver_8_18)

          allocate(solver%difftop%is_inward(10), source= &
            [.False.,.True.,.False.,.True.,.False.,.True.,.False.,.True.,.False.,.True.])

          allocate(solver%diffside%is_inward(4), source=[.False.,.True.,.False.,.True.])

          allocate(solver%dirtop%is_inward(4), source=.True.)
          solver%dirtop%area_divider = 4

          allocate(solver%dirside%is_inward(2), source=.True.)
          solver%dirside%area_divider = 2

        class default
          call CHKERR(1_mpiint, 'unexpected type for solver')
      end select

      solver%difftop%dof= size(solver%difftop%is_inward)
      solver%diffside%dof= size(solver%diffside%is_inward)
      solver%dirtop%dof= size(solver%dirtop%is_inward)
      solver%dirside%dof= size(solver%dirside%is_inward)

      solver%difftop%streams  = solver%difftop%dof/2
      solver%diffside%streams = solver%diffside%dof/2
      solver%dirtop%streams   = solver%dirtop%dof
      solver%dirside%streams  = solver%dirside%dof

      call init_mpi_data_parameters(icomm)

      solver%comm = icomm
      call mpi_comm_rank(solver%comm, solver%myid, ierr)     ; call CHKERR(ierr)
      call mpi_comm_size(solver%comm, solver%numnodes, ierr) ; call CHKERR(ierr)

      call read_commandline_options(solver%comm)

      solver%lenable_solutions_err_estimates = &
        options_max_solution_err.gt.zero .and. options_max_solution_time.gt.zero

      if(.not.approx(dx,dy)) &
        call CHKERR(1_mpiint, 'dx and dy currently have to be the same '//ftoa(dx)//' vs '//ftoa(dy))

      if(ldebug.and.solver%myid.eq.0) then
        print *,'atm dx/dy '//ftoa(dx)//' , '//ftoa(dy)
        print *,'Solver dirtop:', solver%dirtop%is_inward, ':', solver%dirtop%dof
        print *,'Solver dirside:', solver%dirside%is_inward, ':', solver%dirside%dof
        print *,'Solver difftop:', solver%difftop%is_inward, ':', solver%difftop%dof
        print *,'Solver diffside:', solver%diffside%is_inward, ':', solver%diffside%dof
      endif

      call PetscInitialized(lpetsc_is_initialized, ierr); call CHKERR(ierr)
      if(.not.lpetsc_is_initialized) call PetscInitialize(PETSC_NULL_CHARACTER, ierr); call CHKERR(ierr)
#ifdef _XLF
      call PetscPopSignalHandler(ierr); call CHKERR(ierr) ! in case of xlf ibm compilers, remove petsc signal handler -- otherwise we dont get fancy signal traps from boundschecking or FPE's
#endif

      if(present(nxproc) .and. present(nyproc) ) then
        if(ldebug.and.solver%myid.eq.0) print *,'nxproc',shape(nxproc),'::',nxproc
        if(ldebug.and.solver%myid.eq.0) print *,'nyproc',shape(nyproc),'::',nyproc
        call setup_grid(solver, Nz, Nx, Ny, nxproc,nyproc, collapseindex=collapseindex)
      else
        call setup_grid(solver, Nz, max(minimal_dimension, Nx), max(minimal_dimension, Ny), collapseindex=collapseindex)
      endif

      call setup_atm()

      ! init work vectors
      call init_memory(solver%C_dir, solver%C_diff, solver%incSolar, solver%b)

      if(present(solvername)) solver%solvername = trim(solver%solvername)//trim(solvername)
      ! init petsc logging facilities
      call setup_log_events(solver%logs, solver%solvername)

      solver%linitialized = .True.
    else
      print *,solver%myid,'You tried to initialize already initialized PPRTS     &
        &solver. This should not be done. If you need to reinitialize the grids, &
        &call destroy_pprts() first.'
      ierr=1; call CHKERR(ierr)
    endif

    !Todo: this is just here so that we do not break the API. User could also
    !       directly use the set_angle routines?
    call set_angles(solver, phi0, theta0)

  contains
    subroutine setup_atm()
      if(.not.allocated(solver%atm)) allocate(solver%atm)

      solver%atm%dx  = dx
      solver%atm%dy  = dy

      if(.not.allocated(solver%atm%dz) ) then
        allocate(solver%atm%dz(solver%C_one_atm%zs:solver%C_one_atm%ze, &
                               solver%C_one_atm%xs:solver%C_one_atm%xe, &
                               solver%C_one_atm%ys:solver%C_one_atm%ye ), stat=ierr)
      endif

      if(present(dz1d)) then
        do j=solver%C_one_atm%ys,solver%C_one_atm%ye
        do i=solver%C_one_atm%xs,solver%C_one_atm%xe
        solver%atm%dz(:,i,j) = dz1d
        enddo
        enddo
      endif
      if(present(dz3d)) then
        if( any( shape(dz3d).ne.shape(solver%atm%dz) ) ) then
          print *,'Whoops I got a 3D dz definition but this does not correspond to the grid definition :: shapes: ', &
            shape(dz3d), ' vs. ',shape(solver%atm%dz)
          print *,'please know that providing a 3D dz profile has to be of local shape, it can not have global size'
          call MPI_Abort(icomm,1*ierr,ierr)
        endif
        solver%atm%dz = dz3d

      endif
      if(.not.present(dz1d) .and. .not.present(dz3d)) then
        print *,'have to give either dz1d or dz3d in routine call....'
        ierr=1; call CHKERR(ierr)
      endif

      solver%sun%luse_topography = present(dz3d) .and. ltopography  ! if the user supplies 3d height levels and has set the topography option

      if(.not.allocated(solver%atm%l1d)) then
        allocate(solver%atm%l1d( &
          solver%C_one_atm%zs:solver%C_one_atm%ze, &
          solver%C_one_atm%xs:solver%C_one_atm%xe, &
          solver%C_one_atm%ys:solver%C_one_atm%ye ) )
      endif

      !TODO if we have a horiz. staggered grid, this may lead to the point where one 3d box has a outgoing sideward flux but the adjacent
      !1d box does not send anything back --> therefore huge absorption :( -- would need to introduce mirror boundary conditions for
      !sideward fluxes in 1d boxes
      do j=solver%C_one_atm%ys,solver%C_one_atm%ye
        do i=solver%C_one_atm%xs,solver%C_one_atm%xe
          solver%atm%l1d(solver%C_one_atm%ze,i,j) = (solver%atm%dz(solver%C_one_atm%ze,i,j) / solver%atm%dx) .gt. twostr_ratio
          do k=solver%C_one_atm%ze-1,solver%C_one_atm%zs,-1
            solver%atm%l1d(k,i,j) = (solver%atm%dz(k,i,j) / solver%atm%dx) .gt. twostr_ratio
          enddo
        enddo
      enddo
      if(ltwostr_only) solver%atm%l1d = .True.

      if(present(collapseindex)) then
        solver%atm%lcollapse=collapseindex.gt.i1
        solver%atm%icollapse=collapseindex
        solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs),:,:) = .True. ! if need to be collapsed, they have to be 1D.
        if(ldebug) print *,'Using icollapse:',collapseindex, solver%atm%lcollapse
      endif
    end subroutine
  end subroutine

  !> @brief Construct PETSC grid information for regular DMDA
  !> @details setup DMDA grid containers for direct, diffuse and absorption grid
  !>  \n and fill user context containers(t_coord) which include local as well as global array sizes
  !>  \n every mpi rank has to call this
  subroutine setup_grid(solver,Nz_in,Nx,Ny,nxproc,nyproc, collapseindex)
      class(t_solver), intent(inout) :: solver
      integer(iintegers), intent(in) :: Nz_in,Nx,Ny !< @param[in] local number of grid boxes -- in the vertical we have Nz boxes and Nz+1 levels
      integer(iintegers), optional :: nxproc(:), nyproc(:) ! size of local domains on each node
      integer(iintegers), optional, intent(in) :: collapseindex  !< @param[in] collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them

      DMBoundaryType :: bp=DM_BOUNDARY_PERIODIC, bn=DM_BOUNDARY_NONE, bg=DM_BOUNDARY_GHOSTED
      integer(iintegers) :: Nz

      Nz = Nz_in
      if(present(collapseindex)) Nz = Nz_in-collapseindex+i1

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Setting up the DMDA grid for',Nz,Nx,Ny,'using',solver%numnodes,'nodes'

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Configuring DMDA C_diff'
      call setup_dmda(solver%comm, solver%C_diff, Nz+1,Nx,Ny, bp, solver%difftop%dof + 2* solver%diffside%dof)

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Configuring DMDA C_dir'
      if(lcycle_dir) then
        call setup_dmda(solver%comm, solver%C_dir, Nz+1,Nx,Ny, bp, solver%dirtop%dof + 2* solver%dirside%dof)
      else
        call setup_dmda(solver%comm, solver%C_dir, Nz+1,Nx,Ny, bg, solver%dirtop%dof + 2* solver%dirside%dof)
      endif

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Configuring DMDA C_one'
      call setup_dmda(solver%comm, solver%C_one , Nz  , Nx,Ny,  bp, i1)
      call setup_dmda(solver%comm, solver%C_one1, Nz+1, Nx,Ny,  bp, i1)

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Configuring DMDA C_two'
      call setup_dmda(solver%comm, solver%C_two1, Nz+1, Nx,Ny,  bp, i2)

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Configuring DMDA atm'
      call setup_dmda(solver%comm, solver%C_one_atm , Nz_in  , Nx,Ny,  bp, i1)
      call setup_dmda(solver%comm, solver%C_one_atm1, Nz_in+1, Nx,Ny,  bp, i1)

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'DMDA grid ready'
    contains
      subroutine setup_dmda(icomm, C, Nz, Nx, Ny, boundary, dof)
        integer(mpiint), intent(in) :: icomm
        type(t_coord), allocatable :: C
        integer(iintegers), intent(in) :: Nz,Nx,Ny,dof
        DMBoundaryType, intent(in) :: boundary

        integer(iintegers), parameter :: stencil_size=1
        if(ldebug.and.present(nxproc).and.present(nyproc)) print *, solver%myid, 'setup_dmda', nxproc, nyproc

        allocate(C)

        C%comm = icomm

        C%dof = i1*dof
        if(present(nxproc) .and. present(nyproc) ) then
          call DMDACreate3d( C%comm ,                                                     &
            bn                      , boundary                   , boundary             , &
            DMDA_STENCIL_STAR       ,                                                     &
            Nz                      , sum(nxproc)                , sum(nyproc)          , &
            i1                      , size(nxproc,kind=iintegers), size(nyproc,kind=iintegers), &
            C%dof                   , stencil_size               ,                        &
            [Nz]                    , nxproc                     , nyproc               , &
            C%da                    , ierr)
          call CHKERR(ierr)
        else
          call DMDACreate3d( C%comm ,                                                   &
            bn                      , boundary                 , boundary             , &
            DMDA_STENCIL_STAR       ,                                                   &
            i1*Nz                   , Nx                       , Ny                   , &
            i1                      , PETSC_DECIDE             , PETSC_DECIDE         , &
            C%dof                   , stencil_size             ,                        &
            [Nz]                    , PETSC_NULL_INTEGER       , PETSC_NULL_INTEGER   , &
            C%da                    , ierr) ;call CHKERR(ierr)
        endif

        call DMSetMatType(C%da, MATAIJ, ierr)                                ; call CHKERR(ierr)
        if(lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ; call CHKERR(ierr)
        call DMSetFromOptions(C%da, ierr)                                    ; call CHKERR(ierr)
        call DMSetup(C%da,ierr)                                              ; call CHKERR(ierr)

        if(ldebug) call DMView(C%da, PETSC_VIEWER_STDOUT_WORLD ,ierr)        ; call CHKERR(ierr)
        call setup_coords(C)
      end subroutine
      subroutine setup_coords(C)
        type(t_coord) :: C
        DMBoundaryType :: bx, by, bz
        DMDAStencilType :: st
        integer(iintegers) :: stencil_width, nproc_x, nproc_y, nproc_z, Ndof
        integer(mpiint) :: numnodes

        call DMDAGetInfo(C%da, C%dim,                             &
          C%glob_zm, C%glob_xm, C%glob_ym,                        &
          nproc_z, nproc_x, nproc_y, Ndof, stencil_width,         &
          bx, by, bz, st, ierr) ;call CHKERR(ierr)

        call DMDAGetCorners(C%da, C%zs, C%xs,C%ys, C%zm, C%xm,C%ym, ierr) ;call CHKERR(ierr)
        C%xe = C%xs+C%xm-1
        C%ye = C%ys+C%ym-1
        C%ze = C%zs+C%zm-1

        call DMDAGetGhostCorners(C%da,C%gzs,C%gxs,C%gys,C%gzm,C%gxm,C%gym,ierr) ;call CHKERR(ierr)
        C%gxe = C%gxs+C%gxm-1
        C%gye = C%gys+C%gym-1
        C%gze = C%gzs+C%gzm-1

        if(ldebug) then
          print *,solver%myid,'Domain Corners z:: ',C%zs,':',C%ze,' (',C%zm,' entries)','global size',C%glob_zm
          print *,solver%myid,'Domain Corners x:: ',C%xs,':',C%xe,' (',C%xm,' entries)','global size',C%glob_xm
          print *,solver%myid,'Domain Corners y:: ',C%ys,':',C%ye,' (',C%ym,' entries)','global size',C%glob_ym
        endif

        allocate(C%neighbors(0:3**C%dim-1) )
        call DMDAGetNeighbors(C%da,C%neighbors,ierr) ;call CHKERR(ierr)
        call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)
        if(numnodes.gt.i1) then
          if(ldebug.and.(C%dim.eq.3)) print *,'PETSC id',solver%myid,C%dim,'Neighbors are',C%neighbors(10),  &
                                                                                           C%neighbors(4) ,  &
                                                                                           C%neighbors(16),  &
                                                                                           C%neighbors(22),  &
                                                                             'while I am ',C%neighbors(13)

          if(ldebug.and.(C%dim.eq.2)) print *,'PETSC id',solver%myid,C%dim,'Neighbors are',C%neighbors(1), &
                                                                                           C%neighbors(3), &
                                                                                           C%neighbors(7), &
                                                                                           C%neighbors(5), &
                                                                             'while I am ',C%neighbors(4)
        endif
        if(C%glob_xm.lt.i3) call CHKERR(1_mpiint, 'Global domain is too small in x-direction (Nx='//itoa(C%glob_xm)// &
          '). However, need at least 3 because of horizontal ghost cells')
        if(C%glob_ym.lt.i3) call CHKERR(1_mpiint, 'Global domain is too small in y-direction (Ny='//itoa(C%glob_ym)// &
          '). However, need at least 3 because of horizontal ghost cells')
      end subroutine
    end subroutine

    !> @brief initialize basic memory structs like PETSc vectors and matrices
    subroutine init_memory(C_dir, C_diff, incSolar,b)
      type(t_coord), intent(in) :: C_dir, C_diff
      type(tVec), intent(inout), allocatable :: b,incSolar

      if(ltwostr_only) return

      if(.not.allocated(incSolar)) allocate(incSolar)
      if(.not.allocated(b)) allocate(b)

      call DMCreateGlobalVector(C_dir%da,incSolar,ierr) ; call CHKERR(ierr)
      call DMCreateGlobalVector(C_diff%da,b,ierr)       ; call CHKERR(ierr)

      call VecSet(incSolar,zero,ierr) ; call CHKERR(ierr)
      call VecSet(b,zero,ierr)        ; call CHKERR(ierr)
    end subroutine

    subroutine set_angles(solver, phi0, theta0, phi2d, theta2d)
      class(t_solver), intent(inout)   :: solver
      real(ireals),intent(in)          :: phi0           !< @param[in] phi0   solar azmiuth and zenith angle
      real(ireals),intent(in)          :: theta0         !< @param[in] theta0 solar azmiuth and zenith angle
      real(ireals),optional,intent(in) :: phi2d  (:,:)   !< @param[in] phi2d   if given, horizontally varying azimuth
      real(ireals),optional,intent(in) :: theta2d(:,:)   !< @param[in] theta2d if given, and zenith angle

      logical :: lchanged_theta, lchanged_phi

      if(.not.solver%linitialized) then
          print *,solver%myid,'You tried to set angles in the PPRTS solver.  &
              & This should be called right after init_pprts'
          ierr=1; call CHKERR(ierr)
      endif

      lchanged_theta = .True.
      lchanged_phi   = .True.
      if(allocated(solver%sun%theta)) then ! was initialized
        if(present(theta2d)) then
          lchanged_theta = any(.not.approx(theta2d, solver%sun%theta(1,:,:)))
        else
          lchanged_theta = any(.not.approx(theta0 , solver%sun%theta(1,:,:)))
        endif
      endif
      if(allocated(solver%sun%phi)) then
        if(present(phi2d)) then
          lchanged_phi = any(.not.approx(phi2d, solver%sun%phi(1,:,:)))
        else
          lchanged_phi = any(.not.approx(phi0 , solver%sun%phi(1,:,:)))
        endif
      endif
      if(solver%myid.eq.0 .and. ldebug) print *,'pprts set_angles -- changed angles?',lchanged_theta, lchanged_phi
      if(.not. lchanged_theta .and. .not. lchanged_phi) then
        return
      endif

      call setup_suninfo(solver, phi0, theta0, solver%sun, solver%C_one, phi2d=phi2d, theta2d=theta2d)

      if(ltwostr_only .or. lmcrts) return ! dont need anything here, we just compute Twostream anyway

      ! init box montecarlo model
      select type(solver)
        class is (t_solver_1_2)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_1_2::solver%OPP)

        class is (t_solver_3_6)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_3_6::solver%OPP)

        class is (t_solver_3_10)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_3_10::solver%OPP)

        class is (t_solver_8_10)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_8_10::solver%OPP)

        class is (t_solver_3_16)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_3_16::solver%OPP)

        class is (t_solver_8_16)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_8_16::solver%OPP)

        class is (t_solver_8_18)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_8_18::solver%OPP)

        class default
           call CHKERR(1_mpiint, 'init pprts: unexpected type for solver')
      end select

      call solver%OPP%init(solver%comm)

      call init_Matrix(solver, solver%C_dir , solver%Mdir , setup_direct_preallocation)
      call init_Matrix(solver, solver%C_diff, solver%Mdiff, setup_diffuse_preallocation)
  end subroutine


  !> @brief set direction where sun stands
  !> @details save sun azimuth and zenith angle
  !>   \n sun azimuth is reduced to the range of [0,90] and the transmission of direct radiation is contributed for by a integer increment,
  !>   \n determining which neighbouring box is used in the horizontal direction
  subroutine setup_suninfo(solver, phi0, theta0, sun, C_one, phi2d, theta2d)
    class(t_solver), intent(in) :: solver
    real(ireals),intent(in) :: phi0, theta0
    type(t_suninfo),intent(inout) :: sun
    type(t_coord), intent(in) :: C_one
    real(ireals), optional, intent(in) :: phi2d(:,:), theta2d(:,:)

    real(ireals) :: avgphi
    logical :: use_avg_phi, lflg
    integer(mpiint) :: ierr
    integer(iintegers) :: i,j

    call alloc_sun_rfield(sun%symmetry_phi)
    call alloc_sun_rfield(sun%theta)
    call alloc_sun_rfield(sun%phi)
    call alloc_sun_rfield(sun%costheta)
    call alloc_sun_rfield(sun%sintheta)
    call alloc_sun_ifield(sun%xinc)
    call alloc_sun_ifield(sun%yinc)

    if(lforce_phi) then
      sun%phi(:,:,:) = options_phi
    else
      if(present(phi2d)) then
        do j = C_one%ys, C_one%ye
          do i = C_one%xs, C_one%xe
            sun%phi(:,i,j) = phi2d(i-C_one%xs+1,j-C_one%ys+1)
          enddo
        enddo
      else
        sun%phi(:,:,:) = phi0
      endif
      use_avg_phi=.False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-use_avg_phi", &
        use_avg_phi, lflg , ierr) ;call CHKERR(ierr)
      if(use_avg_phi) then
        call imp_allreduce_mean(solver%comm, sun%phi, avgphi)
        sun%phi(:,:,:) = avgphi
      endif
    endif

    if(lforce_theta) then
      sun%theta(:,:,:) = options_theta
    else
      if(present(theta2d)) then
        do j = C_one%ys, C_one%ye
          do i = C_one%xs, C_one%xe
            sun%theta(:,i,j) = theta2d(i-C_one%xs+1,j-C_one%ys+1)
          enddo
        enddo
      else
        sun%theta(:,:,:) = theta0
      endif
    endif

    if(sun%luse_topography) then
      call setup_topography(solver%atm, solver%C_one, solver%C_one1, solver%C_two1, sun)
    endif

    sun%costheta = max( cos(deg2rad(sun%theta)), zero)
    sun%sintheta = max( sin(deg2rad(sun%theta)), zero)

    ! use symmetry for direct beam: always use azimuth [0,90] an just reverse the order where we insert the coeffs
    sun%symmetry_phi = sym_rot_phi(sun%phi)

    where(sin(deg2rad(sun%phi)).gt.zero ) ! phi between 0 and 180 degreee
        sun%xinc=i0
    else where
        sun%xinc=i1
    end where

    where(cos(deg2rad(sun%phi)).lt.zero) ! phi between 90 and 270 degree
        sun%yinc=i1
    else where
        sun%yinc=i0
    end where

    if(ldebug) print *,solver%myid,'setup_dir_inc done', &
      count(sun%xinc.eq.0),count(sun%xinc.eq.i1), &
      count(sun%yinc.eq.0),count(sun%yinc.eq.i1), &
      '::', minval(sun%xinc), maxval(sun%xinc), minval(sun%yinc), maxval(sun%yinc)

    contains
        pure elemental function sym_rot_phi(phi)
            real(ireals) :: sym_rot_phi
            real(ireals),intent(in) :: phi
            ! ''swap'' phi axis down to the range of [0,180]
            sym_rot_phi = acos(cos(deg2rad(phi)))
            !print *,'1st phi swap',phi,' :: ',sym_rot_phi,'=',phi*pi/180,cos(phi*pi/180),acos(cos(phi*pi/180))
            ! and then mirror it onto range [0,90]
            sym_rot_phi = min(90._ireals, max(0._ireals, rad2deg( asin(sin(sym_rot_phi)) )))
            !print *,'2nd phi swap',phi,' :: ',sym_rot_phi,'=',&
            ! sin(sym_rot_phi),asin(sin(sym_rot_phi)),asin(sin(sym_rot_phi)) /pi * 180,int(asin(sin(sym_rot_phi)) /pi * 180)
        end function

        subroutine alloc_sun_rfield(f)
          real(ireals), allocatable, dimension(:,:,:) :: f
          if(.not.allocated(f)) allocate(f(C_one%zs:C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye))
        end subroutine
        subroutine alloc_sun_ifield(f)
          integer(iintegers), allocatable, dimension(:,:,:) :: f
          if(.not.allocated(f)) allocate(f(C_one%zs:C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye))
        end subroutine
  end subroutine

  !> @brief setup topography information
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information and
  !>   \n tweak the local sun angles to bend the rays.
  subroutine setup_topography(atm, C_one, C_one1, C_two1, sun)
    type(t_atmosphere), intent(in) :: atm
    type(t_coord), intent(in) :: C_one, C_one1, C_two1
    type(t_suninfo),intent(inout) :: sun

    real(ireals),Pointer :: grad(:,:,:,:)=>null(), grad_1d(:)=>null()

    integer(iintegers) :: i,j,k
    real(ireals) :: newtheta, newphi, xsun(3)
    real(ireals) :: rotmat(3, 3), newxsun(3)
    integer(mpiint) :: comm, myid, ierr

    if(.not.allocated(sun%theta)) call CHKERR(1_mpiint, 'You called  setup_topography() &
        &but the sun struct is not yet up, make sure setup_suninfo is called before')
    if(.not.allocated(atm%dz)) call CHKERR(1_mpiint, 'You called  setup_topography() &
        &but the atm struct is not yet up, make sure we have atm%dz before')

    call compute_gradient(atm, C_one1, C_two1, atm%hgrad)

    call PetscObjectGetComm(C_one%da, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call getVecPointer(atm%hgrad , C_two1%da, grad_1d, grad)

    rotmat = reshape((/ one , zero, zero,  &
                        zero, one , zero,  &
                        nil , nil , one/), &
                     (/3, 3/), order=(/2, 1/) )

    do j=C_one%ys,C_one%ye
      do i=C_one%xs,C_one%xe
        do k=C_one%zs,C_one%ze

          ! if we are at the global boundary we have to take care that the gradient does not get too
          ! steep. That wouldnt make sense for cyclic boundaries... ! todo: check what we should do?
          !if(i.eq.i0 .or. i.eq.C_one%glob_xm-i1 .or. j.eq.i0 .or. j.eq.C_one%glob_ym-i1) then
          !   !print *,solver%myid,'global edge at:',i,j
          !    cycle
          !endif

          ! Vector of sun direction
          xsun(1) = sin(deg2rad(sun%theta(k,i,j)))*sin(deg2rad(sun%phi(k,i,j)))
          xsun(2) = sin(deg2rad(sun%theta(k,i,j)))*cos(deg2rad(sun%phi(k,i,j)))
          xsun(3) = cos(deg2rad(sun%theta(k,i,j)))

          xsun = xsun / norm2(xsun)

          rotmat(3, 1) = grad(i0, k, i, j)
          rotmat(3, 2) = grad(i1, k, i, j)

          newxsun = matmul(rotmat, xsun)
          newxsun = newxsun / norm2(newxsun)

          newtheta = rad2deg(acos(newxsun(3)))

          !newphi in meteorologiecal definitions: clockwise from y-axis
          newphi = rad2deg(atan2(newxsun(1), newxsun(2)))

          !if(i.eq.C_one1%xs.and.k.eq.C_one%ze) print *,myid,i,j,k, &
          !  '::', sun%theta(k,i,j), newtheta, &
          !  '::', sun%phi  (k,i,j), newphi

          sun%theta(k,i,j) = max(zero, min( 90._ireals, newtheta ))
          sun%phi  (k,i,j) = newphi
        enddo
      enddo
    enddo

    call restoreVecPointer(atm%hgrad, grad_1d, grad)
  end subroutine

  !> @brief create PETSc matrix and inserts diagonal elements
  !> @details one important step for performance is to set preallocation of matrix structure.
  !>  \n  i.e. determine the number of local vs. remote number of entries in each row.
  !>  \n  DMDA actually provides a preallocation but this assumes that degrees of freedom on neighbouring boxes are fully connected to local ones.
  !>  \n  this does of course drastically overestimate non-zeros as we need only the streams that actually send radiation in the respective direction.
  !>  \n  at the moment preallocation routines determine nonzeros by manually checking bounadries --
  !>  \n  !todo we should really use some form of iterating through the entries as it is done in the matrix assembly routines and just flag the rows
  subroutine init_Matrix(solver, C, A, prealloc_subroutine)
    interface
      subroutine preallocation_sub(solver, C, d_nnz, o_nnz)
        import :: t_solver, t_coord, iintegers
        class(t_solver), intent(in) :: solver
        type(t_coord), intent(in) :: C
        integer(iintegers),allocatable :: d_nnz(:)
        integer(iintegers),allocatable :: o_nnz(:)
      end subroutine
    end interface

    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    type(tMat), allocatable, intent(inout) :: A
    procedure(preallocation_sub), optional :: prealloc_subroutine

    integer(iintegers),dimension(:),allocatable :: o_nnz,d_nnz
    integer(mpiint) :: numnodes

    if(.not.allocated(A)) then
      allocate(A)
      call DMCreateMatrix(C%da, A, ierr) ;call CHKERR(ierr)
    endif

    call mpi_comm_size(C%comm, numnodes, ierr) ; call CHKERR(ierr)

    if(numnodes.gt.1) then
      if(lprealloc.and.present(prealloc_subroutine)) then
        call prealloc_subroutine(solver, C, d_nnz, o_nnz)
        call MatMPIAIJSetPreallocation(A, C%dof+1, d_nnz, C%dof, o_nnz, ierr) ;call CHKERR(ierr)
      else ! poor mans perallocation uses way more memory...
        call CHKERR(1_mpiint, 'init_Matrix::setPreallocation : poor mans preallocation should really not be used...')
        call MatMPIAIJSetPreallocation(A, C%dof+1, PETSC_NULL_INTEGER, C%dof, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)
      endif
    else
      call MatSeqAIJSetPreallocation(A, C%dof+i1, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)
    endif

    ! If matrix is resetted, keep nonzero pattern and allow to non-zero allocations -- those should not be many
    ! call MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr) ;call CHKERR(ierr)

    ! pressure mesh  may wiggle a bit and change atm%l1d -- keep the nonzeros flexible
    !call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr) ;call CHKERR(ierr)

    ! call MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr) ;call CHKERR(ierr) ! dont throw away the zero -- this completely destroys preallocation performance

    call MatSetUp(A,ierr) ;call CHKERR(ierr)

    call mat_set_diagonal(A)
  contains
    subroutine mat_set_diagonal(A,vdiag)
      type(tMat) :: A
      real(ireals),intent(in),optional :: vdiag

      integer(iintegers) :: is, ie, irow
      real(ireals) :: v

      if(present(vdiag)) then
        v = vdiag
      else
        v = one
      endif

      call MatGetOwnershipRange(A, is, ie, ierr); call CHKERR(ierr)
      do irow = is, ie-1
        call MatSetValue(A, irow, irow, v, INSERT_VALUES, ierr); call CHKERR(ierr)
      enddo
    end subroutine

  end subroutine

  subroutine setup_direct_preallocation(solver, C, d_nnz, o_nnz)
    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    integer(iintegers),allocatable :: d_nnz(:)
    integer(iintegers),allocatable :: o_nnz(:)
    type(tVec) :: v_o_nnz,v_d_nnz
    type(tVec) :: g_o_nnz,g_d_nnz
    real(ireals),Pointer :: xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
    real(ireals),Pointer :: xo1d(:)=>null(),xd1d(:)=>null()

    integer(iintegers) :: vsize, i, j, k, isrc, idst, xinc, yinc, src, dst, icnt

    logical :: llocal_src, llocal_dst
    integer(mpiint) :: myid, ierr
    MatStencil :: row(4,C%dof), col(4,C%dof)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)


    if(myid.eq.0.and.ldebug) print *,myid,'building direct o_nnz for mat with',C%dof,'dof'
    call DMGetLocalVector(C%da,v_o_nnz,ierr) ;call CHKERR(ierr)
    call DMGetLocalVector(C%da,v_d_nnz,ierr) ;call CHKERR(ierr)

    call getVecPointer(v_o_nnz, C%da, xo1d, xo)
    call getVecPointer(v_d_nnz, C%da, xd1d, xd)

    xd = zero
    xo = zero

    icnt = -1
    do j=C%ys,C%ye
      do i=C%xs,C%xe
        do k=C%zs,C%ze-1

          if( solver%atm%l1d(atmk(solver%atm,k), i, j) ) then
            do idst = i0, solver%dirtop%dof-1
              call inc( xd(idst, k+1, i,j), one )
            enddo
          else
            xinc = solver%sun%xinc(k,i,j)
            yinc = solver%sun%yinc(k,i,j)

            do idst = 1, solver%dirtop%dof
              dst = idst
              row(MatStencil_j,dst) = i
              row(MatStencil_k,dst) = j
              row(MatStencil_i,dst) = k+1
              row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
            enddo

            do idst = 1, solver%dirside%dof
              dst = idst + solver%dirtop%dof
              row(MatStencil_j,dst) = i+xinc
              row(MatStencil_k,dst) = j
              row(MatStencil_i,dst) = k
              row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
            enddo

            do idst = 1, solver%dirside%dof
              dst = idst + solver%dirtop%dof + solver%dirside%dof
              row(MatStencil_j,dst) = i
              row(MatStencil_k,dst) = j+yinc
              row(MatStencil_i,dst) = k
              row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
            enddo

            do isrc = 1, solver%dirtop%dof
              src = isrc
              col(MatStencil_j,src) = i
              col(MatStencil_k,src) = j
              col(MatStencil_i,src) = k
              col(MatStencil_c,src) = src-i1 ! Define transmission towards the lower/upper lid
            enddo

            do isrc = 1, solver%dirside%dof
              src = isrc + solver%dirtop%dof
              col(MatStencil_j,src) = i+1-xinc
              col(MatStencil_k,src) = j
              col(MatStencil_i,src) = k
              col(MatStencil_c,src) = src-i1 ! Define transmission towards the left/right lid
            enddo

            do isrc = 1, solver%dirside%dof
              src = isrc + solver%dirtop%dof + solver%dirside%dof
              col(MatStencil_j,src) = i
              col(MatStencil_k,src) = j+1-yinc
              col(MatStencil_i,src) = k
              col(MatStencil_c,src) = src-i1 ! Define transmission towards the front/back lid
            enddo

            do idst = 1,C%dof
              icnt = icnt+1
              llocal_dst = .True.
              if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 .and. row(MatStencil_j,idst).lt.C%xs ) llocal_dst = .False. ! have real neighbor west  and is not local entry
              if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 .and. row(MatStencil_j,idst).gt.C%xe ) llocal_dst = .False. ! have real neighbor east  and is not local entry
              if( C%neighbors( 4).ne.myid .and. C%neighbors( 4).ge.i0 .and. row(MatStencil_k,idst).lt.C%ys ) llocal_dst = .False. ! have real neighbor south and is not local entry
              if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 .and. row(MatStencil_k,idst).gt.C%ye ) llocal_dst = .False. ! have real neighbor north and is not local entry

              do isrc = 1,C%dof
                llocal_src = .True.

                if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 .and. col(MatStencil_j,isrc).lt.C%xs ) llocal_src = .False. ! have real neighbor west  and is not local entry
                if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 .and. col(MatStencil_j,isrc).gt.C%xe ) llocal_src = .False. ! have real neighbor east  and is not local entry
                if( C%neighbors( 4).ne.myid .and. C%neighbors( 4).ge.i0 .and. col(MatStencil_k,isrc).lt.C%ys ) llocal_src = .False. ! have real neighbor south and is not local entry
                if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 .and. col(MatStencil_k,isrc).gt.C%ye ) llocal_src = .False. ! have real neighbor north and is not local entry

                !if(myid.eq.0) print *,myid,icnt,k,i,j,'::',idst,isrc,'::',llocal_dst,llocal_src
                if(llocal_dst .and. llocal_src) then
                  call inc(xd(row(4,idst),row(3,idst),row(2,idst),row(1,idst)), one)
                else
                  call inc(xo(row(4,idst),row(3,idst),row(2,idst),row(1,idst)), one)
                endif
              enddo
            enddo
          endif
        enddo
      enddo
    enddo

    call restoreVecPointer(v_o_nnz, xo1d, xo)
    call restoreVecPointer(v_d_nnz, xd1d, xd)

    call DMGetGlobalVector(C%da,g_o_nnz,ierr) ;call CHKERR(ierr)
    call DMGetGlobalVector(C%da,g_d_nnz,ierr) ;call CHKERR(ierr)
    call VecSet(g_o_nnz, zero, ierr); call CHKERR(ierr)
    call VecSet(g_d_nnz, zero, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da,v_o_nnz,ADD_VALUES,g_o_nnz,ierr) ;call CHKERR(ierr)
    call DMLocalToGlobalEnd  (C%da,v_o_nnz,ADD_VALUES,g_o_nnz,ierr) ;call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da,v_d_nnz,ADD_VALUES,g_d_nnz,ierr) ;call CHKERR(ierr)
    call DMLocalToGlobalEnd  (C%da,v_d_nnz,ADD_VALUES,g_d_nnz,ierr) ;call CHKERR(ierr)

    call DMRestoreLocalVector(C%da,v_o_nnz,ierr) ;call CHKERR(ierr)
    call DMRestoreLocalVector(C%da,v_d_nnz,ierr) ;call CHKERR(ierr)

    call getVecPointer(g_o_nnz, C%da, xo1d, xo)
    call getVecPointer(g_d_nnz, C%da, xd1d, xd)

    call VecGetLocalSize(g_d_nnz,vsize,ierr) ;call CHKERR(ierr)
    allocate(o_nnz(0:vsize-1))
    allocate(d_nnz(0:vsize-1))

    o_nnz=int(xo1d, kind=iintegers)
    d_nnz=int(xd1d, kind=iintegers) + i1  ! +1 for diagonal entries

    call restoreVecPointer(g_o_nnz, xo1d, xo)
    call restoreVecPointer(g_d_nnz, xd1d, xd)

    call DMRestoreGlobalVector(C%da,g_o_nnz,ierr) ;call CHKERR(ierr)
    call DMRestoreGlobalVector(C%da,g_d_nnz,ierr) ;call CHKERR(ierr)

    if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz), &
      'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
  end subroutine

  subroutine setup_diffuse_preallocation(solver, C, d_nnz, o_nnz)
    class(t_solver), intent(in) :: solver
    type(t_coord), intent(in) :: C
    integer(iintegers),allocatable :: d_nnz(:)
    integer(iintegers),allocatable :: o_nnz(:)
    type(tVec) :: v_o_nnz,v_d_nnz
    type(tVec) :: g_o_nnz,g_d_nnz
    real(ireals),Pointer :: xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
    real(ireals),Pointer :: xo1d(:)=>null(),xd1d(:)=>null()

    integer(iintegers) :: vsize, i, j, k, isrc, idst, src, dst, icnt
    integer(iintegers) :: dst_id, src_id, idof
    integer(mpiint) :: myid, ierr
    MatStencil :: row(4,0:C%dof-1), col(4,0:C%dof-1)

    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)

    if(myid.eq.0.and.ldebug) print *,myid,'building diffuse o_nnz for mat with',C%dof,'dof'
    call DMGetLocalVector(C%da, v_o_nnz, ierr) ;call CHKERR(ierr)
    call DMGetLocalVector(C%da, v_d_nnz, ierr) ;call CHKERR(ierr)

    call getVecPointer(v_o_nnz, C%da, xo1d, xo)
    call getVecPointer(v_d_nnz, C%da, xd1d, xd)


    xd = zero
    xo = zero

    icnt = -1
    do j=C%ys,C%ye
      do i=C%xs,C%xe
        do k=C%zs,C%ze-1

          if( solver%atm%l1d(atmk(solver%atm, k),i,j) ) then
            do idof=1, solver%difftop%dof
              if (solver%difftop%is_inward(idof)) then
                call inc( xd(idof-1, k+1, i, j), real(solver%difftop%dof, ireals) )
              else
                call inc( xd(idof-1, k  , i, j), real(solver%difftop%dof, ireals) )
              endif
            enddo
          else

            ! Diffuse Coefficients Code Begin
            do idof=1, solver%difftop%dof
              src = idof-1
              if (solver%difftop%is_inward(idof)) then
                col(MatStencil_j,src) = i
                col(MatStencil_k,src) = j
                col(MatStencil_i,src) = k
                col(MatStencil_c,src) = src
              else
                col(MatStencil_j,src) = i
                col(MatStencil_k,src) = j
                col(MatStencil_i,src) = k+1
                col(MatStencil_c,src) = src
              endif
            enddo

            do idof=1, solver%diffside%dof
              src = solver%difftop%dof + idof -1
              if (solver%diffside%is_inward(idof)) then
                col(MatStencil_j,src) = i
                col(MatStencil_k,src) = j
                col(MatStencil_i,src) = k
                col(MatStencil_c,src) = src
              else
                col(MatStencil_j,src) = i+1
                col(MatStencil_k,src) = j
                col(MatStencil_i,src) = k
                col(MatStencil_c,src) = src
              endif
            enddo

            do idof=1, solver%diffside%dof
              src = solver%difftop%dof + solver%diffside%dof + idof -1
              if (solver%diffside%is_inward(idof)) then
                col(MatStencil_j,src) = i
                col(MatStencil_k,src) = j
                col(MatStencil_i,src) = k
                col(MatStencil_c,src) = src
              else
                col(MatStencil_j,src) = i
                col(MatStencil_k,src) = j+1
                col(MatStencil_i,src) = k
                col(MatStencil_c,src) = src
              endif
            enddo

            do idof=1, solver%difftop%dof
              dst = idof-1
              if (solver%difftop%is_inward(idof)) then
                row(MatStencil_j,dst) = i
                row(MatStencil_k,dst) = j
                row(MatStencil_i,dst) = k+1
                row(MatStencil_c,dst) = dst
              else
                row(MatStencil_j,dst) = i
                row(MatStencil_k,dst) = j
                row(MatStencil_i,dst) = k
                row(MatStencil_c,dst) = dst
              endif
            enddo

            do idof=1, solver%diffside%dof
              dst = solver%difftop%dof + idof-1
              if (solver%diffside%is_inward(idof)) then
                row(MatStencil_j,dst) = i+1
                row(MatStencil_k,dst) = j
                row(MatStencil_i,dst) = k
                row(MatStencil_c,dst) = dst
              else
                row(MatStencil_j,dst) = i
                row(MatStencil_k,dst) = j
                row(MatStencil_i,dst) = k
                row(MatStencil_c,dst) = dst
              endif
            enddo

            do idof=1, solver%diffside%dof
              dst = solver%difftop%dof + solver%diffside%dof + idof-1
              if (solver%diffside%is_inward(idof)) then
                row(MatStencil_j,dst) = i
                row(MatStencil_k,dst) = j+1
                row(MatStencil_i,dst) = k
                row(MatStencil_c,dst) = dst
              else
                row(MatStencil_j,dst) = i
                row(MatStencil_k,dst) = j
                row(MatStencil_i,dst) = k
                row(MatStencil_c,dst) = dst
              endif
            enddo
            ! Diffuse Coefficients Code End

            do idst = 0,C%dof-1
              icnt = icnt+1
              dst_id = myid
              if( C%neighbors(10).ne.myid.and.C%neighbors(10).ge.i0.and.row(MatStencil_j,idst).lt.C%xs ) dst_id = C%neighbors(10) ! have real neighbor west  and is not local entry
              if( C%neighbors(16).ne.myid.and.C%neighbors(16).ge.i0.and.row(MatStencil_j,idst).gt.C%xe ) dst_id = C%neighbors(16) ! have real neighbor east  and is not local entry
              if( C%neighbors( 4).ne.myid.and.C%neighbors( 4).ge.i0.and.row(MatStencil_k,idst).lt.C%ys ) dst_id = C%neighbors( 4) ! have real neighbor south and is not local entry
              if( C%neighbors(22).ne.myid.and.C%neighbors(22).ge.i0.and.row(MatStencil_k,idst).gt.C%ye ) dst_id = C%neighbors(22) ! have real neighbor north and is not local entry

              do isrc = 0,C%dof-1
                src_id = myid

                if( C%neighbors(10).ne.myid.and.C%neighbors(10).ge.i0.and.col(MatStencil_j,isrc).lt.C%xs ) src_id = C%neighbors(10) ! have real neighbor west  and is not local entry
                if( C%neighbors(16).ne.myid.and.C%neighbors(16).ge.i0.and.col(MatStencil_j,isrc).gt.C%xe ) src_id = C%neighbors(16) ! have real neighbor east  and is not local entry
                if( C%neighbors( 4).ne.myid.and.C%neighbors( 4).ge.i0.and.col(MatStencil_k,isrc).lt.C%ys ) src_id = C%neighbors( 4) ! have real neighbor south and is not local entry
                if( C%neighbors(22).ne.myid.and.C%neighbors(22).ge.i0.and.col(MatStencil_k,isrc).gt.C%ye ) src_id = C%neighbors(22) ! have real neighbor north and is not local entry

                if(src_id .eq. dst_id) then
                  call inc(xd(row(4,idst),row(3,idst),row(2,idst),row(1,idst)), one)
                else
                  call inc(xo(row(4,idst),row(3,idst),row(2,idst),row(1,idst)), one)
                endif

              enddo
            enddo

          endif ! atm_1d / 3d?
        enddo ! k

        ! Surface entries E_up
        do idof=1, solver%difftop%dof
          src = idof-1
          if (.not.solver%difftop%is_inward(idof)) then
            call inc( xd(src, C%ze, i,j), real(solver%difftop%streams, ireals) )
          endif
        enddo

      enddo
    enddo

    call restoreVecPointer(v_o_nnz, xo1d, xo)
    call restoreVecPointer(v_d_nnz, xd1d, xd)

    call DMGetGlobalVector(C%da,g_o_nnz,ierr) ;call CHKERR(ierr)
    call DMGetGlobalVector(C%da,g_d_nnz,ierr) ;call CHKERR(ierr)
    call VecSet(g_o_nnz, zero, ierr); call CHKERR(ierr)
    call VecSet(g_d_nnz, zero, ierr); call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da,v_o_nnz,ADD_VALUES,g_o_nnz,ierr) ;call CHKERR(ierr)
    call DMLocalToGlobalEnd  (C%da,v_o_nnz,ADD_VALUES,g_o_nnz,ierr) ;call CHKERR(ierr)

    call DMLocalToGlobalBegin(C%da,v_d_nnz,ADD_VALUES,g_d_nnz,ierr) ;call CHKERR(ierr)
    call DMLocalToGlobalEnd  (C%da,v_d_nnz,ADD_VALUES,g_d_nnz,ierr) ;call CHKERR(ierr)

    call DMRestoreLocalVector(C%da,v_o_nnz,ierr) ;call CHKERR(ierr)
    call DMRestoreLocalVector(C%da,v_d_nnz,ierr) ;call CHKERR(ierr)

    call getVecPointer(g_o_nnz, C%da, xo1d, xo)
    call getVecPointer(g_d_nnz, C%da, xd1d, xd)

    call VecGetLocalSize(g_d_nnz, vsize, ierr) ;call CHKERR(ierr)
    allocate(o_nnz(0:vsize-1))
    allocate(d_nnz(0:vsize-1))

    o_nnz=int(xo1d, kind=iintegers)
    d_nnz=int(xd1d, kind=iintegers) + i1  ! +1 for diagonal entries

    call restoreVecPointer(g_o_nnz, xo1d, xo)
    call restoreVecPointer(g_d_nnz, xd1d, xd)

    call DMRestoreGlobalVector(C%da, g_o_nnz, ierr) ;call CHKERR(ierr)
    call DMRestoreGlobalVector(C%da, g_d_nnz, ierr) ;call CHKERR(ierr)

    if(myid.eq.0 .and. ldebug) print *,myid,'diffuse d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),&
      'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)

  end subroutine

  subroutine mat_info(comm, A)
    MPI_Comm, intent(in) :: comm
    type(tMat) :: A
    double precision :: info(MAT_INFO_SIZE)
    double precision :: mal, nz_allocated, nz_used, nz_unneeded
    integer(iintegers) :: m, n
    integer(mpiint) :: myid, ierr

    call mpi_comm_rank(comm, myid, ierr)     ; call CHKERR(ierr)

    call MatGetInfo(A,MAT_LOCAL,info,ierr) ;call CHKERR(ierr)
    mal          = info(MAT_INFO_MALLOCS)
    nz_allocated = info(MAT_INFO_NZ_ALLOCATED)
    nz_used      = info(MAT_INFO_NZ_USED)
    nz_unneeded  = info(MAT_INFO_NZ_UNNEEDED)

    call MatGetOwnershipRange(A, m,n, ierr)

    if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_MALLOCS',mal,'MAT_INFO_NZ_ALLOCATED',nz_allocated
    if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_USED',nz_used,'MAT_INFO_NZ_unneded',nz_unneeded
    if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: Ownership range',m,n
  end subroutine
  subroutine set_optical_properties(solver, albedo, &
      kabs, ksca, g, &
      planck, planck_srfc, &
      albedo_2d, &
      ldelta_scaling)
    class(t_solver)                                   :: solver
    real(ireals), intent(in)                          :: albedo
    real(ireals),intent(in),dimension(:,:,:),optional :: kabs, ksca, g ! dimensions (Nz  , Nx, Ny)
    real(ireals),intent(in),dimension(:,:,:),optional :: planck                    ! dimensions (Nz+1, Nx, Ny) planck radiation on levels
    real(ireals),intent(in),dimension(:,:),optional   :: planck_srfc               ! dimensions (Nx, Ny)
    real(ireals),intent(in),dimension(:,:),optional   :: albedo_2d                 ! dimensions (Nx, Ny)
    logical, intent(in), optional :: ldelta_scaling ! determines if we should try to delta scale these optprops

    real(ireals)        :: pprts_delta_scale_max_g
    integer(iintegers)  :: k, i, j
    logical :: lpprts_delta_scale, lflg

    call PetscLogEventBegin(solver%logs%set_optprop, ierr); call CHKERR(ierr)

    associate( atm => solver%atm, &
        C_one_atm => solver%C_one_atm, &
        C_one_atm1 => solver%C_one_atm1, &
        sun => solver%sun, &
        C_one => solver%C_one)

    if(.not.allocated(atm%kabs) )  &
      allocate( atm%kabs(C_one_atm%zs :C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye) )
    if(.not.allocated(atm%ksca) )  &
      allocate( atm%ksca(C_one_atm%zs :C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye) )
    if(.not.allocated(atm%g) )     &
      allocate( atm%g   (C_one_atm%zs :C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye) )

    if(.not.allocated(atm%albedo)) allocate(atm%albedo(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
    atm%albedo = albedo
    if(present(albedo_2d)) atm%albedo = albedo_2d

    if(present(kabs) ) atm%kabs = kabs
    if(present(ksca) ) atm%ksca = ksca
    if(present(g   ) ) atm%g    = g

    if(present(planck) ) then
      if(.not.allocated(atm%planck)) &
        allocate(atm%planck(C_one_atm1%zs:C_one_atm1%ze, C_one_atm1%xs:C_one_atm1%xe, C_one_atm1%ys:C_one_atm1%ye))
      atm%planck = planck
    else
      if(allocated(atm%planck)) deallocate(atm%planck)
    endif

    if(present(planck_srfc) ) then
      if(.not.allocated(atm%Bsrfc)) &
        allocate(atm%Bsrfc(C_one_atm1%xs:C_one_atm1%xe, C_one_atm1%ys:C_one_atm1%ye))
      atm%Bsrfc = planck_srfc
    else
      if(allocated(atm%Bsrfc)) deallocate(atm%Bsrfc)
    endif

    if(ldebug) then
      if( any([kabs,ksca,g].lt.zero) ) then
        print *,solver%myid,'set_optical_properties :: found illegal value in local_optical properties! abort!'
        do k=C_one_atm%zs,C_one_atm%ze
          print *,solver%myid,k,'kabs',kabs(k,:,:)
          print *,solver%myid,k,'ksca',ksca(k,:,:)
        enddo
        call CHKERR(1_mpiint, 'set_optical_properties :: found illegal value in local_optical properties! '//&
                              ' '//itoa(solver%myid)//&
                              ' kabs min '//ftoa(minval(kabs))//' max '//ftoa(maxval(kabs))//&
                              ' ksca min '//ftoa(minval(ksca))//' max '//ftoa(maxval(ksca))//&
                              ' g    min '//ftoa(minval(g   ))//' max '//ftoa(maxval(g   )))
      endif
      if( any(isnan([kabs,ksca,g]))) then
        call CHKERR(1_mpiint, 'set_optical_properties :: found NaN value in optical properties!'//&
                              ' NaN in kabs? '//ltoa(any(isnan(kabs)))// &
                              ' NaN in ksca? '//ltoa(any(isnan(ksca)))// &
                              ' NaN in g   ? '//ltoa(any(isnan(g   ))))
      endif
    endif
    if(ldebug) then
      if( (any([atm%kabs,atm%ksca,atm%g].lt.zero)) .or. (any(isnan([atm%kabs,atm%ksca,atm%g]))) ) then
        call CHKERR(1_mpiint, 'set_optical_properties :: found illegal value in optical properties! '//&
                              ' '//itoa(solver%myid)//&
                              ' kabs min '//ftoa(minval(atm%kabs))//' max '//ftoa(maxval(atm%kabs))//&
                              ' ksca min '//ftoa(minval(atm%ksca))//' max '//ftoa(maxval(atm%ksca))//&
                              ' g    min '//ftoa(minval(atm%g   ))//' max '//ftoa(maxval(atm%g   )))
      endif
    endif

    if(ldebug.and.solver%myid.eq.0) then
      if(present(kabs) ) then
        print *,'atm_kabs     ',maxval(atm%kabs  )  ,shape(atm%kabs  )
      endif
      if(present(ksca) ) then
        print *,'atm_ksca     ',maxval(atm%ksca  )  ,shape(atm%ksca  )
      endif
      if(present(g) ) then
        print *,'atm_g        ',maxval(atm%g     )  ,shape(atm%g     )
      endif
      if(present(planck) ) then
        print *,'atm_planck   ',maxval(atm%planck   )  ,shape(atm%planck   )
      endif

      print *,'Number of 1D layers: ', count(atm%l1d) , size(atm%l1d),'(',(100._ireals* count(atm%l1d) )/size(atm%l1d),'%)'
      if(present(kabs)) print *,'init local optprop:', shape(kabs), '::', shape(atm%kabs)
    endif

    lpprts_delta_scale = get_arg(.True., ldelta_scaling)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_delta_scale", &
      lpprts_delta_scale, lflg , ierr) ;call CHKERR(ierr)

    if(lpprts_delta_scale) then
      pprts_delta_scale_max_g=.649_ireals
      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_delta_scale_max_g", &
        pprts_delta_scale_max_g, lflg , ierr) ;call CHKERR(ierr)

      call delta_scale(atm%kabs, atm%ksca, atm%g, max_g=pprts_delta_scale_max_g)
    else
      if(solver%myid.eq.0.and.lflg) print *,"Skipping Delta scaling of optprops"
      if(any(atm%g.ge.0.649_ireals)) &
        call CHKWARN(1_mpiint, 'Skipping delta scaling but now we have values of '// &
        'g > '//ftoa(pprts_delta_scale_max_g)// &
        ' ('//ftoa(maxval(atm%g))//')')
    endif

    if(ltwostr_only) then
      if(ldebug .and. solver%myid.eq.0) then
        do k=C_one_atm%zs,C_one_atm%ze
          if(present(planck)) then
            print *,solver%myid,'Optical Properties:',k, &
              'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
              minval(atm%kabs(k,:,:)), minval(atm%ksca(k,:,:)), minval(atm%g(k,:,:)),&
              maxval(atm%kabs(k,:,:)), maxval(atm%ksca(k,:,:)), maxval(atm%g(k,:,:)),&
              '::',minval(atm%planck (k,:,:)),maxval(atm%planck        (k, :,:))
          else
            print *,solver%myid,'Optical Properties:',k, &
              'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
              minval(atm%kabs(k,:,:)), minval(atm%ksca(k,:,:)), minval(atm%g(k,:,:)),&
              maxval(atm%kabs(k,:,:)), maxval(atm%ksca(k,:,:)), maxval(atm%g(k,:,:))
          endif
        enddo
      endif
      return ! twostream should not depend on eddington coeffs... it will have to calculate it on its own.
    endif

    if(ldebug) then
      if( (any([atm%kabs,atm%ksca,atm%g].lt.zero)) .or. (any(isnan([atm%kabs,atm%ksca,atm%g]))) ) then
        print *,solver%myid,'set_optical_properties :: found illegal value in delta_scaled optical properties! abort!'
      endif
    endif

    if(luse_eddington) then
      if(.not.allocated(atm%a11)) allocate(atm%a11(C_one_atm%zs:C_one_atm%ze,C_one_atm%xs:C_one_atm%xe,C_one_atm%ys:C_one_atm%ye))   ! allocate space for twostream coefficients
      if(.not.allocated(atm%a12)) allocate(atm%a12(C_one_atm%zs:C_one_atm%ze,C_one_atm%xs:C_one_atm%xe,C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a21)) allocate(atm%a21(C_one_atm%zs:C_one_atm%ze,C_one_atm%xs:C_one_atm%xe,C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a22)) allocate(atm%a22(C_one_atm%zs:C_one_atm%ze,C_one_atm%xs:C_one_atm%xe,C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a13)) allocate(atm%a13(C_one_atm%zs:C_one_atm%ze,C_one_atm%xs:C_one_atm%xe,C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a23)) allocate(atm%a23(C_one_atm%zs:C_one_atm%ze,C_one_atm%xs:C_one_atm%xe,C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a33)) allocate(atm%a33(C_one_atm%zs:C_one_atm%ze,C_one_atm%xs:C_one_atm%xe,C_one_atm%ys:C_one_atm%ye))
    endif

    if(luse_eddington) then
      !DIR$ IVDEP
      do j=C_one_atm%ys,C_one_atm%ye
        do i=C_one_atm%xs,C_one_atm%xe
          do k=C_one_atm%zs,C_one_atm%ze
            if( atm%l1d(k,i,j) ) then
              call eddington_coeff_zdun ( &
                atm%dz(k,i,j) * max(tiny(one), atm%kabs(k,i,j) + atm%ksca(k,i,j)), & ! tau
                atm%ksca(k,i,j) / max(tiny(one), atm%kabs(k,i,j) + atm%ksca(k,i,j)), & ! w0
                atm%g(k,i,j), &
                sun%costheta(i1,i,j), &
                atm%a11(k,i,j),          &
                atm%a12(k,i,j),          &
                atm%a13(k,i,j),          &
                atm%a23(k,i,j),          &
                atm%a33(k,i,j))
            else
              !TODO :: we should really not have this memeory accessible at all....
              !     :: the fix would be trivial at the moment, as long as all 1d layers start at same 'k',
              !     :: however if that is to change i.e. on staggered grid, we would need staggered array constructs....
              if(ldebug) then
                atm%a11(k,i,j) = nil
                atm%a12(k,i,j) = nil
                atm%a13(k,i,j) = nil
                atm%a23(k,i,j) = nil
                atm%a33(k,i,j) = nil
              endif !ldebug
            endif !l1d
          enddo !k
        enddo !i
      enddo !j

      ! Set symmetric up and down transport coefficients. If only for one
      ! layer, they would be indeed symmetric. If we collapse the atmosphere
      ! however, we have different transmission if we come from top or from
      ! bottom.
      atm%a21 = atm%a12
      atm%a22 = atm%a11
      call handle_atm_collapse()
    endif

    if(ldebug .and. solver%myid.eq.0) then
      do k=C_one_atm%zs,C_one_atm%ze
        if(present(planck)) then
          print *,solver%myid,'Optical Properties:',k,&
            'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
            minval(atm%kabs(k,:,:)), minval(atm%ksca(k,:,:)), minval(atm%g(k,:,:)),&
            maxval(atm%kabs(k,:,:)), maxval(atm%ksca(k,:,:)), maxval(atm%g(k,:,:)),&
            '::',minval(atm%planck (k,:,:)),maxval(atm%planck        (k, :,:))
        else
          print *,solver%myid,'Optical Properties:',k,&
            'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
            minval(atm%kabs(k,:,:)), minval(atm%ksca(k,:,:)), minval(atm%g(k,:,:)),&
            maxval(atm%kabs(k,:,:)), maxval(atm%ksca(k,:,:)), maxval(atm%g(k,:,:)),&
            '::',minval(atm%a33 (k,:,:)),maxval(atm%a33(k,:,:))
        endif
      enddo
    endif
    end associate
    call PetscLogEventEnd(solver%logs%set_optprop, ierr); call CHKERR(ierr)

    contains
      subroutine handle_atm_collapse()
        real(ireals), allocatable :: Eup(:), Edn(:)
        integer(iintegers) :: i,j,ak
        associate( atm => solver%atm, &
            C_one     => solver%C_one, &
            C_one_atm => solver%C_one_atm )

          if(atm%lcollapse) then
            ak = atmk(atm, C_one%zs)
            if(ak.ne.i1) then
              if(present(planck)) then
                allocate(Edn(C_one_atm%zs:ak+1), Eup(C_one_atm%zs:ak+1))
                if(.not.allocated(atm%Bbot)) &
                  allocate(atm%Bbot(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
                if(.not.allocated(atm%Btop)) &
                  allocate(atm%Btop(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
              endif
              do j=C_one_atm%ys,C_one_atm%ye
                do i=C_one_atm%xs,C_one_atm%xe
                  if(present(planck)) then
                    call adding(&
                      atm%a11(C_one_atm%zs:ak, i, j), &
                      atm%a12(C_one_atm%zs:ak, i, j), &
                      atm%a21(C_one_atm%zs:ak, i, j), &
                      atm%a22(C_one_atm%zs:ak, i, j), &
                      atm%a13(C_one_atm%zs:ak, i, j), &
                      atm%a23(C_one_atm%zs:ak, i, j), &
                      atm%a33(C_one_atm%zs:ak, i, j), &
                      atm%dz(C_one_atm%zs:ak, i, j) * atm%kabs(C_one_atm%zs:ak, i, j), &
                      atm%planck(C_one_atm%zs:ak+1, i, j), &
                      Eup, Edn, atm%Btop(i, j), atm%Bbot(i, j))
                  else
                    call adding(&
                      atm%a11(C_one_atm%zs:ak, i, j), &
                      atm%a12(C_one_atm%zs:ak, i, j), &
                      atm%a21(C_one_atm%zs:ak, i, j), &
                      atm%a22(C_one_atm%zs:ak, i, j), &
                      atm%a13(C_one_atm%zs:ak, i, j), &
                      atm%a23(C_one_atm%zs:ak, i, j), &
                      atm%a33(C_one_atm%zs:ak, i, j))
                  endif
                enddo !i
              enddo !j
            endif
          endif !lcollapse
        end associate
      end subroutine
      subroutine adding(a11,a12,a21,a22,a13,a23,a33,dtau,planck,Eup,Edn,Btop,Bbot)
        real(ireals),intent(inout),dimension(:) :: a11,a12,a21,a22,a13,a23,a33
        real(ireals),intent(in)   ,dimension(:),optional :: dtau, planck
        real(ireals),intent(out)  ,dimension(:),optional :: Eup, Edn
        real(ireals),intent(out)  ,optional :: Btop, Bbot
        real(ireals) :: t, r, rdir, sdir, tdir

        integer(iintegers) :: N ! = size(a11)
        integer(iintegers) :: k
        real(ireals) :: rl, tl, Tbot, Ttop, Rbot, Rtop
        N = size(a11)

        t = a11(1)
        r = a12(1)

        tdir = a33(1)
        rdir = a13(1)
        sdir = a23(1)

        ! Reflectivity as seen from top
        do k=2,N
          rl = r
          tl = t

          r = r + (a12(k) * t**2) / (one - r*a12(k))
          t = t * a11(k) / (one - rl * a12(k))

          sdir = (a11(k) * sdir + tdir * a13(k) * rl * a11(k)) / (one - rl * a12(k)) + tdir * a23(k)
          rdir = rdir + ( tdir * a13(k) + sdir * a12(k) ) * tl
          tdir = tdir*a33(k)
        enddo

        Ttop = t
        Rtop = r

        a13(N) = rdir
        a23(N) = sdir
        a33(N) = tdir

        t = a22(N)
        r = a21(N)
        ! Reflectivity as seen from bottom
        do k=N-1,1,-1
          rl = r
          tl = t

          r = a12(k) + (r * a11(k)**2) / (one - r * a12(k))
          t = t * a11(k) / (one - rl * a21(k))
        enddo

        Tbot = t
        Rbot = r

        a12(N) = Rtop
        a22(N) = Ttop
        a11(N) = Tbot
        a21(N) = Rbot

        a11(1:N-1) = nil
        a12(1:N-1) = nil
        a21(1:N-1) = nil
        a22(1:N-1) = nil

        a13(1:N-1) = nil
        a23(1:N-1) = nil
        a33(1:N-1) = nil

        if(present(planck)) then
          call schwarzschild(3_iintegers, dtau, 0._ireals, Edn, Eup, planck, &
            opt_srfc_emission=0._ireals)
          Bbot = Edn(N+1) / pi
          Btop = Eup(1) / pi
        endif
      end subroutine
    end subroutine

  subroutine set_global_optical_properties(solver, global_albedo, global_kabs, global_ksca, global_g, global_planck)
    class(t_solver),intent(in) :: solver
    real(ireals),intent(inout),optional :: global_albedo
    real(ireals),intent(inout),dimension(:,:,:),allocatable,optional :: global_kabs, global_ksca, global_g
    real(ireals),intent(inout),dimension(:,:,:),allocatable,optional :: global_planck
    real(ireals),dimension(:,:,:),allocatable :: local_kabs, local_ksca, local_g
    real(ireals),dimension(:,:,:),allocatable :: local_planck
    real(ireals) :: local_albedo
    logical :: lhave_planck, lhave_kabs, lhave_ksca, lhave_g

    call PetscLogEventBegin(solver%logs%set_optprop, ierr); call CHKERR(ierr)

    if(.not.solver%linitialized) then
      call CHKERR(1_mpiint, 'You tried to set global optical properties but pprts environment seems not to be initialized....'// &
        'please call init first!')
    endif

    lhave_kabs   = present(global_kabs  ); call imp_bcast(solver%comm, lhave_kabs  , 0_mpiint)
    lhave_ksca   = present(global_ksca  ); call imp_bcast(solver%comm, lhave_ksca  , 0_mpiint)
    lhave_g      = present(global_g     ); call imp_bcast(solver%comm, lhave_g     , 0_mpiint)
    lhave_planck = present(global_planck); call imp_bcast(solver%comm, lhave_planck, 0_mpiint)

    ! Make sure that our domain has at least 3 entries in each dimension.... otherwise violates boundary conditions
    if(solver%myid.eq.0) then
      if(present(global_albedo)) local_albedo = global_albedo
      if( lhave_kabs   ) call extend_arr(global_kabs)
      if( lhave_ksca   ) call extend_arr(global_ksca)
      if( lhave_g      ) call extend_arr(global_g)
      if( lhave_planck ) call extend_arr(global_planck)
    endif

    ! Scatter global optical properties to MPI nodes
    call local_optprop()
    ! Now global_fields are local to mpi subdomain.
    call PetscLogEventEnd(solver%logs%set_optprop, ierr); call CHKERR(ierr)

    if(lhave_planck) then
      call set_optical_properties(solver, local_albedo, local_kabs, local_ksca, local_g, local_planck)
    else
      call set_optical_properties(solver, local_albedo, local_kabs, local_ksca, local_g)
    endif
  contains
    subroutine local_optprop()
      type(tVec) :: local_vec

      if(solver%myid.eq.0.and.ldebug .and. lhave_kabs) &
        print *,solver%myid,'copying optprop: global to local :: shape kabs',shape(global_kabs),&
        'xstart/end',solver%C_one_atm%xs,solver%C_one_atm%xe,&
        'ystart/end',solver%C_one_atm%ys,solver%C_one_atm%ye

      call imp_bcast(solver%comm, local_albedo, 0_mpiint)

      call DMGetGlobalVector(solver%C_one_atm%da, local_vec, ierr) ; call CHKERR(ierr)

      if(lhave_kabs) then
        call scatterZerotoPetscGlobal(global_kabs, solver%C_one_atm%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm%da, local_kabs)
      endif

      if(lhave_ksca) then
        call scatterZerotoPetscGlobal(global_ksca, solver%C_one_atm%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm%da, local_ksca)
      endif

      if(lhave_g) then
        call scatterZerotoPetscGlobal(global_g, solver%C_one_atm%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm%da, local_g)
      endif

      call DMRestoreGlobalVector(solver%C_one_atm%da, local_vec, ierr) ; call CHKERR(ierr)

      if(lhave_planck) then
        call DMGetGlobalVector(solver%C_one_atm1%da, local_vec, ierr) ; call CHKERR(ierr)
        call scatterZerotoPetscGlobal(global_planck, solver%C_one_atm1%da, local_vec)
        call petscVecToF90(local_vec, solver%C_one_atm1%da, local_planck)
        call DMRestoreGlobalVector(solver%C_one_atm1%da, local_vec, ierr) ; call CHKERR(ierr)
      endif
    end subroutine
    subroutine extend_arr(arr)
      real(ireals),intent(inout),allocatable :: arr(:,:,:)
      real(ireals),allocatable :: tmp(:,:)
      integer(iintegers) :: dims(3),i

      if(.not. allocated(arr) ) print *,solver%myid,'ERROR in SUBROUTINE extend_arr :: Cannot extend non allocated array!'

      dims = shape(arr)
      if( dims(2) .eq. 1 ) then
        allocate( tmp(dims(1),dims(3) ), SOURCE=arr(:, 1, :) )
        deallocate(arr)
        allocate( arr(dims(1), minimal_dimension, dims(3) ) )
        do i=1,minimal_dimension
          arr(:, i, :) = tmp
        enddo
        deallocate(tmp)
      endif

      dims = shape(arr)
      if( dims(3) .eq. 1 ) then
        allocate( tmp(dims(1),dims(2) ), SOURCE=arr(:,:,1) )
        deallocate(arr)
        allocate( arr(dims(1), dims(2), minimal_dimension ) )
        do i=1,minimal_dimension
          arr(:, :, i) = tmp
        enddo
        deallocate(tmp)
      endif
      if(any(shape(arr).lt.minimal_dimension) ) &
        call CHKERR(1_mpiint, 'set_optprop -> extend_arr :: dimension is smaller than we support.. please think of something here')
    end subroutine

  end subroutine

  subroutine solve_pprts(solver, edirTOA, opt_solution_uid, opt_solution_time)
    class(t_solver), intent(inout)          :: solver
    real(ireals),intent(in)                 :: edirTOA
    integer(iintegers),optional,intent(in)  :: opt_solution_uid
    real(ireals),      optional,intent(in)  :: opt_solution_time

    integer(iintegers) :: uid
    logical            :: lsolar
    logical            :: lskip_diffuse_solve, lflg

    associate(  solutions => solver%solutions, &
                C_dir     => solver%C_dir,     &
                C_diff    => solver%C_diff,    &
                Mdir      => solver%Mdir,      &
                Mdiff     => solver%Mdiff     )

    if(.not.allocated(solver%atm)) call CHKERR(1_mpiint, 'atmosphere is not allocated?!')
    if(.not.allocated(solver%atm%kabs)) &
      call CHKERR(1_mpiint, 'atmosphere%kabs is not allocated! - maybe you need to call set_optical_properties() first')
    if(.not.allocated(solver%atm%ksca)) &
      call CHKERR(1_mpiint, 'atmosphere%ksca is not allocated! - maybe you need to call set_optical_properties() first')

    uid = get_arg(0_iintegers, opt_solution_uid)

    lsolar = mpi_logical_and(solver%comm, edirTOA.gt.zero .and. any(solver%sun%theta.ge.zero))
    if(ldebug) print *,'uid', uid, 'lsolar', lsolar, 'edirTOA', edirTOA, ':', any(solver%sun%theta.ge.zero)

    if(.not.solutions(uid)%lset) then
      call prepare_solution(solver%C_dir%da, solver%C_diff%da, solver%C_one%da, &
        lsolar=lsolar, solution=solutions(uid), uid=uid)
    endif

    ! --------- Skip Thermal Computation (-lskip_thermal) --
    if(lskip_thermal .and. (solutions(uid)%lsolar_rad.eqv..False.) ) then !
      if(ldebug .and. solver%myid.eq.0) print *,'skipping thermal calculation -- returning zero flux'
      call VecSet(solutions(uid)%ediff, zero, ierr); call CHKERR(ierr)
      solutions(uid)%lchanged=.True.
      goto 99 ! quick exit
    endif

    ! --------- Calculate 1D Radiative Transfer ------------
    if(  ltwostr &
      .or. all(solver%atm%l1d.eqv..True.) &
      .or. ((solutions(uid)%lsolar_rad.eqv..False.) .and. lcalc_nca) &
      .or. ((solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild) ) then


      if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild ) then
        call PetscLogEventBegin(solver%logs%solve_schwarzschild, ierr)
        call schwarz(solver, solutions(uid))
        call PetscLogEventEnd(solver%logs%solve_schwarzschild, ierr)
      else
        call PetscLogEventBegin(solver%logs%solve_twostream, ierr)
        call twostream(solver, edirTOA,  solutions(uid) )
        call PetscLogEventEnd(solver%logs%solve_twostream, ierr)
      endif

      if(ldebug .and. solver%myid.eq.0) print *,'1D calculation done'


      if( ltwostr_only ) goto 99
      if( all(solver%atm%l1d.eqv..True.) ) goto 99
      if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lcalc_nca ) goto 99
      if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild ) goto 99
    endif

    if( lmcrts ) then
      call PetscLogEventBegin(solver%logs%solve_mcrts, ierr)
      call solve_mcrts(solver, edirTOA, solutions(uid))
      call PetscLogEventEnd(solver%logs%solve_mcrts, ierr)
      goto 99
    endif

    ! --------- scale from [W/m**2] to [W] -----------------
    call scale_flx(solver, solutions(uid), lWm2=.False. )

    ! ---------------------------- Edir  -------------------
    if( solutions(uid)%lsolar_rad ) then
      call PetscLogEventBegin(solver%logs%compute_Edir, ierr)

      call setup_incSolar(solver, solver%incSolar, edirTOA)

      call PetscLogEventBegin(solver%logs%setup_Mdir, ierr)
      call set_dir_coeff(solver, solver%sun, solver%Mdir, C_dir)
      call PetscLogEventEnd(solver%logs%setup_Mdir, ierr)

      if(ldebug) call mat_info(solver%comm, solver%Mdir)

      call PetscLogEventBegin(solver%logs%solve_Mdir, ierr)
      call setup_ksp(solver%atm, solver%ksp_solar_dir, C_dir, solver%Mdir, "solar_dir_")
      call solve(solver, &
        solver%ksp_solar_dir, &
        solver%incSolar, &
        solutions(uid)%edir, &
        solutions(uid)%Niter_dir, &
        solutions(uid)%dir_ksp_residual_history)

      call PetscLogEventEnd(solver%logs%solve_Mdir, ierr)

      solutions(uid)%lchanged=.True.
      solutions(uid)%lWm2_dir=.False.
      call PetscObjectSetName(solutions(uid)%edir,'debug_edir',ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solutions(uid)%edir, PETSC_NULL_VEC, "-show_debug_edir", ierr); call CHKERR(ierr)
      call PetscLogEventEnd(solver%logs%compute_Edir, ierr)
    endif


    ! ---------------------------- Source Term -------------
    call PetscLogEventBegin(solver%logs%compute_Ediff, ierr)
    call setup_b(solver, solutions(uid),solver%b)

    lskip_diffuse_solve = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-skip_diffuse_solve", lskip_diffuse_solve, lflg , ierr) ;call CHKERR(ierr)
    if(lskip_diffuse_solve) then
      call VecCopy(solver%b, solutions(uid)%ediff, ierr)
    else
      ! ---------------------------- Ediff -------------------
      call PetscLogEventBegin(solver%logs%setup_Mdiff, ierr)
      call set_diff_coeff(solver, Mdiff,C_diff)
      call PetscLogEventEnd(solver%logs%setup_Mdiff, ierr)

      if(ldebug) call mat_info(solver%comm, solver%Mdiff)

      call PetscLogEventBegin(solver%logs%solve_Mdiff, ierr)
      if( solutions(uid)%lsolar_rad ) then
        call setup_ksp(solver%atm, solver%ksp_solar_diff, C_diff, Mdiff, "solar_diff_")
        call solve(solver, &
          solver%ksp_solar_diff, &
          solver%b, &
          solutions(uid)%ediff, &
          solutions(uid)%Niter_diff, &
          solutions(uid)%diff_ksp_residual_history)
      else
        call setup_ksp(solver%atm, solver%ksp_thermal_diff, C_diff, Mdiff, "thermal_diff_")
        call solve(solver, &
          solver%ksp_thermal_diff, &
          solver%b, &
          solutions(uid)%ediff, &
          solutions(uid)%Niter_diff, &
          solutions(uid)%diff_ksp_residual_history)
      endif
      call PetscLogEventEnd(solver%logs%solve_Mdiff, ierr)
    endif

    solutions(uid)%lchanged=.True.
    solutions(uid)%lWm2_diff=.False. !Tenstream solver returns fluxes as [W]

    call PetscLogEventEnd(solver%logs%compute_Ediff, ierr)

    99 continue ! this is the quick exit final call where we clean up before the end of the routine

    call restore_solution(solver, solutions(uid), opt_solution_time)

    end associate
  end subroutine

  !> @brief renormalize fluxes with the size of a face(sides or lid)
  subroutine scale_flx(solver, solution, lWm2)
    class(t_solver), intent(inout)        :: solver
    type(t_state_container),intent(inout) :: solution   !< @param solution container with computed fluxes
    logical,intent(in)                    :: lWm2  !< @param determines direction of scaling, if true, scale to W/m**2

    call PetscLogEventBegin(solver%logs%scale_flx, ierr); call CHKERR(ierr)
    if(solution%lsolar_rad) then
      if(.not.allocated(solver%dir_scalevec_Wm2_to_W)) then
        allocate(solver%dir_scalevec_Wm2_to_W)
        call VecDuplicate(solution%edir, solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        call gen_scale_dir_flx_vec(solver, solver%dir_scalevec_Wm2_to_W, solver%C_dir)
      endif
      if(.not.allocated(solver%dir_scalevec_W_to_Wm2)) then
        allocate(solver%dir_scalevec_W_to_Wm2)
        call VecDuplicate(solver%dir_scalevec_Wm2_to_W, solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        call VecSet(solver%dir_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
        call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
          solver%dir_scalevec_W_to_Wm2, &
          solver%dir_scalevec_W_to_Wm2, &
          solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      endif
      if(solution%lWm2_dir .neqv. lWm2) then
        if(lWm2) then
          call VecPointwiseMult(solution%edir, solution%edir, solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        else
          call VecPointwiseMult(solution%edir, solution%edir, solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        endif
        solution%lWm2_dir = lWm2
      endif
    endif

    if(.not.allocated(solver%diff_scalevec_Wm2_to_W)) then
      allocate(solver%diff_scalevec_Wm2_to_W)
      call VecDuplicate(solution%ediff, solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      call gen_scale_diff_flx_vec(solver, solver%diff_scalevec_Wm2_to_W, solver%C_diff)
    endif
    if(.not.allocated(solver%diff_scalevec_W_to_Wm2)) then
      allocate(solver%diff_scalevec_W_to_Wm2)
      call VecDuplicate(solver%diff_scalevec_Wm2_to_W, solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
      call VecSet(solver%diff_scalevec_W_to_Wm2, one, ierr); call CHKERR(ierr)
      call VecPointwiseDivide( &  ! Computes 1./scalevec_Wm2_to_W
        solver%diff_scalevec_W_to_Wm2, &
        solver%diff_scalevec_W_to_Wm2, &
        solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
    endif

    if(solution%lWm2_diff .neqv. lWm2) then
      if(lWm2) then
        call VecPointwiseMult(solution%ediff, solution%ediff, solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
      else
        call VecPointwiseMult(solution%ediff, solution%ediff, solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
      endif
      solution%lWm2_diff = lWm2
    endif
    call PetscLogEventEnd(solver%logs%scale_flx, ierr); call CHKERR(ierr)

  contains
    subroutine gen_scale_dir_flx_vec(solver, v, C)
      class(t_solver)      :: solver
      type(tVec)           :: v
      type(t_coord)        :: C
      real(ireals),pointer :: xv  (:,:,:,:) =>null()
      real(ireals),pointer :: xv1d(:)       =>null()
      integer(iintegers)   :: i, j, k, d, iside
      real(ireals)         :: Ax, Ay, Az, fac

      associate( atm     => solver%atm )

      if(solver%myid.eq.0.and.ldebug) print *,'rescaling direct fluxes',C%zm,C%xm,C%ym
      call getVecPointer(v ,C%da ,xv1d, xv)

      Az  = solver%atm%dx*solver%atm%dy / real(solver%dirtop%area_divider, ireals)  ! size of a direct stream in m**2
      fac = Az

      ! Scaling top faces
      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze
            do iside=1,solver%dirtop%dof
              d = iside-1
              xv(d,k,i,j) = fac
            enddo
          enddo
        enddo
      enddo

      ! Scaling side faces
      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-1
            ! First the faces in x-direction
            Ax = solver%atm%dy*solver%atm%dz(k,i,j) / real(solver%dirside%area_divider, ireals)
            fac = Ax
            do iside=1,solver%dirside%dof
              d = solver%dirtop%dof + iside-1
              xv(d,k,i,j) = fac
            enddo

            ! Then the rest of the faces in y-direction
            Ay = atm%dy*atm%dz(k,i,j) / real(solver%dirside%area_divider, ireals)
            fac = Ay
            do iside=1,solver%dirside%dof
              d = solver%dirtop%dof + solver%dirside%dof + iside-1
              xv(d,k,i,j) = fac
            enddo
          enddo
          ! the side faces underneath the surface are always scaled by unity
          xv(solver%dirtop%dof:ubound(xv,dim=1),k,i,j) = one
        enddo
      enddo
      call restoreVecPointer(v, xv1d, xv )
      end associate
    end subroutine
    subroutine gen_scale_diff_flx_vec(solver, v, C)
      class(t_solver)      :: solver
      type(tVec)           :: v
      type(t_coord)        :: C
      real(ireals),pointer :: xv(:,:,:,:)=>null()
      real(ireals),pointer :: xv1d(:)=>null()

      integer(iintegers)  :: iside, src, i, j, k
      real(ireals)        :: Az, Ax, Ay, fac

      if(solver%myid.eq.0.and.ldebug) print *,'rescaling fluxes',C%zm,C%xm,C%ym
      call getVecPointer(v ,C%da ,xv1d, xv)

      if(C%dof.eq.i3 .or. C%dof.eq.i8) then
        print *,'scale_flx_vec is just for diffuse radia'
      endif

      ! Scaling top faces
      Az = solver%atm%dx*solver%atm%dy / real(solver%difftop%area_divider, ireals)
      fac = Az

      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze
            do iside=1,solver%difftop%dof
              src = iside -1
              xv(src ,k,i,j) = fac                  ! diffuse radiation
            enddo
          enddo
        enddo
      enddo

      ! Scaling side faces
      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-1
              ! faces in x-direction
              Ax = solver%atm%dy*solver%atm%dz(k,i,j) / real(solver%diffside%area_divider, ireals)
              fac = Ax

              do iside=1,solver%diffside%dof
                src = solver%difftop%dof + iside -1
                xv(src ,k,i,j) = fac
              enddo

              ! faces in y-direction
              Ay = solver%atm%dx*solver%atm%dz(k,i,j) / real(solver%difftop%area_divider, ireals)
              fac = Ay

              do iside=1,solver%diffside%dof
                src = solver%difftop%dof + solver%diffside%dof + iside -1
                xv(src ,k,i,j) = fac
              enddo
          enddo
          ! the side faces underneath the surface are always scaled by unity
          xv(solver%difftop%dof:ubound(xv,dim=1),k,i,j) = one
        enddo
      enddo
      call restoreVecPointer(v, xv1d, xv )

    end subroutine
  end subroutine

  subroutine restore_solution(solver, solution, time)
    ! restore_solution:: if flux have changed, we need to update absorption, save the residual history
    class(t_solver)         :: solver
    type(t_state_container) :: solution
    real(ireals),intent(in),optional :: time

    character(default_str_len) :: vecname
    real(ireals) :: inf_norm
    type(tVec) :: abso_old

    if( .not. solution%lset ) call CHKERR(1_mpiint, 'cant restore solution that was not initialized')
    if( .not. allocated(solution%abso) ) call CHKERR(1_mpiint, 'cant restore solution that was not initialized')

    if( .not. solution%lchanged ) call CHKERR(1_mpiint, 'cant restore solution which was not changed')

    if(present(time) .and. solver%lenable_solutions_err_estimates) then ! Create working vec to determine difference between old and new absorption vec
      call DMGetGlobalVector(solver%C_one%da, abso_old, ierr) ; call CHKERR(ierr)
      call VecCopy( solution%abso, abso_old, ierr)     ; call CHKERR(ierr)
    endif

    ! update absorption
    call calc_flx_div(solver, solution)

    if(ldebug .and. solver%myid.eq.0) &
      print *,'Saving Solution ',solution%uid

    ! make sure to bring the fluxes into [W/m**2]
    call scale_flx(solver, solution, lWm2=.True. )


    if(ldebug .and. solver%myid.eq.0) &
      print *,'Saving Solution done'
    solution%lchanged=.False.

    if(present(time) .and. solver%lenable_solutions_err_estimates) then ! Compute norm between old absorption and new one
      call VecAXPY(abso_old, -one, solution%abso, ierr); call CHKERR(ierr) ! overwrite abso_old with difference to new one
      call VecNorm(abso_old, NORM_INFINITY, inf_norm, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(solver%C_one%da, abso_old, ierr)   ; call CHKERR(ierr)

      ! Save norm for later analysis
      solution%maxnorm = eoshift ( solution%maxnorm, shift = -1) !shift all values by 1 to the right
      solution%time    = eoshift ( solution%time   , shift = -1) !shift all values by 1 to the right

      solution%maxnorm( 1 ) = inf_norm
      solution%time( 1 )    = time

      if(ldebug .and. solver%myid.eq.0) &
        print *,'Updating error statistics for solution ',solution%uid,'at time ',time,'::',solution%time(1), &
                ':: norm', inf_norm,'[W] :: hr_norm approx:',inf_norm*86.1,'[K/d]'

    endif !present(time) .and. solver%lenable_solutions_err_estimates


    if(allocated(solution%edir)) then
      write(vecname,FMT='("edir",I0)') solution%uid
      call PetscObjectSetName(solution%edir,vecname,ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_edir", ierr); call CHKERR(ierr)
    endif

    if(allocated(solution%ediff)) then
      write(vecname,FMT='("ediff",I0)') solution%uid
      call PetscObjectSetName(solution%ediff,vecname,ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-show_ediff", ierr); call CHKERR(ierr)
    endif

    if(allocated(solution%abso)) then
      write(vecname,FMT='("abso",I0)') solution%uid
      call PetscObjectSetName(solution%abso,vecname,ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, "-show_abso", ierr); call CHKERR(ierr)
    endif

    if(allocated(solver%b)) then
      write(vecname,FMT='("b",I0)') solution%uid
      call PetscObjectSetName(solver%b,vecname,ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solver%b, PETSC_NULL_VEC, "-show_b", ierr); call CHKERR(ierr)
    endif

    if(allocated(solver%incSolar)) then
      write(vecname,FMT='("incSolar",I0)') solution%uid
      call PetscObjectSetName(solver%incSolar,vecname,ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solver%incSolar, PETSC_NULL_VEC, "-show_incSolar", ierr); call CHKERR(ierr)
    endif
  end subroutine


  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine twostream(solver, edirTOA, solution)
    class(t_solver), intent(inout) :: solver
    real(ireals),intent(in)       :: edirTOA
    type(t_state_container)       :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: xv_dir=>null(),xv_diff=>null()
    real(ireals),pointer,dimension(:) :: xv_dir1d=>null(),xv_diff1d=>null()
    integer(iintegers) :: i,j,src

    real(ireals),allocatable :: dtau(:),kext(:),w0(:),g(:),S(:),Edn(:),Eup(:)
    real(ireals) :: mu0,incSolar,fac

    associate(atm         => solver%atm, &
              C_diff      => solver%C_diff, &
              C_dir       => solver%C_dir, &
              C_one_atm   => solver%C_one_atm, &
              C_one_atm1  => solver%C_one_atm1)

    if(solution%lsolar_rad) then
      call PetscObjectSetName(solution%edir,'twostream_edir_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
      call VecSet(solution%edir ,zero,ierr); call CHKERR(ierr)
    endif

    call PetscObjectSetName(solution%ediff,'twostream_ediff_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
    call VecSet(solution%ediff,zero,ierr); call CHKERR(ierr)

    allocate( dtau(C_one_atm%zm) )
    allocate( kext(C_one_atm%zm) )
    allocate(   w0(C_one_atm%zm) )
    allocate(    g(C_one_atm%zm) )

    if(solution%lsolar_rad) &
      call getVecPointer(solution%edir  ,C_dir%da  ,xv_dir1d , xv_dir)
    call getVecPointer(solution%ediff ,C_diff%da ,xv_diff1d, xv_diff)

    allocate( S  (C_one_atm1%zs:C_one_atm1%ze) )
    allocate( Eup(C_one_atm1%zs:C_one_atm1%ze) )
    allocate( Edn(C_one_atm1%zs:C_one_atm1%ze) )

    do j=C_one_atm%ys,C_one_atm%ye
      do i=C_one_atm%xs,C_one_atm%xe

        mu0 = solver%sun%costheta(C_one_atm1%zs,i,j)
        incSolar = edirTOA* mu0

        kext = atm%kabs(:,i,j) + atm%ksca(:,i,j)
        dtau = atm%dz(:,i,j)* kext
        w0   = atm%ksca(:,i,j) / max(kext, epsilon(kext))
        g    = atm%g(:,i,j)

        if(allocated(atm%planck) ) then
          if(allocated(atm%Bsrfc)) then
            call delta_eddington_twostream(dtau, w0, g, &
              mu0, incSolar, atm%albedo(i,j), &
              S, Edn, Eup, &
              planck=atm%planck(:,i,j), &
              planck_srfc=atm%Bsrfc(i,j))
          else
            call delta_eddington_twostream(dtau, w0, g, &
              mu0, incSolar, atm%albedo(i,j), &
              S, Edn, Eup, &
              planck=atm%planck(:,i,j) )
          endif
        else
          !
          ! call adding_delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup )
          !
          !TODO investigate if this one is really ok...
          ! I recently had valgrind errors in VecNorm after calling this:
          ! make -j ex_pprts_rrtm_lw_sw &&
          ! mpirun -np 1 -wdir ../examples/rrtm_lw_sw/ valgrind $(pwd)/bin/ex_pprts_rrtm_lw_sw -twostr_only
          call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup )
        endif

        if(solution%lsolar_rad) then
          fac = real(solver%dirtop%area_divider, ireals) / real(solver%dirtop%streams, ireals)
          do src=i0,solver%dirtop%dof-1
            xv_dir(src,C_dir%zs+1:C_dir%ze,i,j) = S(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
            xv_dir(src,C_dir%zs           ,i,j) = S(C_one_atm1%zs) * fac
          enddo
        endif

        fac = real(solver%difftop%area_divider, ireals) / real(solver%difftop%streams, ireals)
        do src = 1, solver%difftop%dof
          if(solver%difftop%is_inward(src)) then
            xv_diff(src-1,C_diff%zs+1:C_diff%ze,i,j) = Edn(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
            xv_diff(src-1,C_diff%zs            ,i,j) = Edn(C_one_atm1%zs) * fac
          else
            xv_diff(src-1,C_diff%zs+1:C_diff%ze,i,j) = Eup(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
            xv_diff(src-1,C_diff%zs            ,i,j) = Eup(C_one_atm1%zs) * fac
          endif
        enddo
      enddo
    enddo

    if(solution%lsolar_rad) &
      call restoreVecPointer(solution%edir, xv_dir1d, xv_dir  )
    call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

    !Twostream solver returns fluxes as [W]
    solution%lWm2_dir  = .True.
    solution%lWm2_diff = .True.
    ! and mark solution that it is not up to date
    solution%lchanged  = .True.

    deallocate(S)
    deallocate(Edn)
    deallocate(Eup)

    end associate
  end subroutine

  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine schwarz(solver, solution)
    class(t_solver)         :: solver
    type(t_state_container) :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: xv_diff=>null()
    real(ireals),pointer,dimension(:)       :: xv_diff1d=>null()
    integer(iintegers) :: i,j,k,idof
    integer(iintegers) :: Nmu, ak
    logical :: lflg

    real(ireals),allocatable :: dtau(:),Edn(:),Eup(:)


    associate( C_diff   => solver%C_diff, &
                atm     => solver%atm, &
                C_one   => solver%C_one, &
                C_one1  => solver%C_one1)

    if(solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling schwarschild solver for solar calculation -- stopping!')
    if(.not.allocated(atm%planck)) call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

    call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

    allocate(dtau(size(atm%dz,dim=1)))

    call getVecPointer(solution%ediff, C_diff%da, xv_diff1d, xv_diff)

    allocate( Eup(0:size(atm%dz,dim=1)) )
    allocate( Edn(0:size(atm%dz,dim=1)) )

    Nmu = 10
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
      "-schwarzschild_Nmu" , Nmu, lflg , ierr) ;call CHKERR(ierr)

    if(solver%myid.eq.0 .and. ldebug) print *,' CALCULATING schwarzschild ::'

    do j = C_diff%ys, C_diff%ye
      do i = C_diff%xs, C_diff%xe

        dtau = atm%dz(:, i, j) * atm%kabs(:, i, j)

        if(allocated(atm%Bsrfc)) then
          call schwarzschild(Nmu, dtau, atm%albedo(i,j), Edn, Eup, &
            atm%planck(:, i, j), opt_srfc_emission=atm%Bsrfc(i,j))
        else
          call schwarzschild(Nmu, dtau, atm%albedo(i,j), Edn, Eup, &
            atm%planck(:, i, j))
        endif

        ! icollapse needs special case for TOA flx's
        do idof = 0, solver%difftop%dof-1
          if (solver%difftop%is_inward(i1+idof)) then ! Edn
            xv_diff(idof,C_diff%zs,i,j) = Edn(0)
          else ! Eup
            xv_diff(idof,C_diff%zs,i,j) = Eup(0)
          endif
        enddo

        ! rest of the atmosphere
        do k=C_diff%zs+1,C_diff%ze
          ak = atmk(atm,k)
          do idof = 0, solver%difftop%dof-1
            if (solver%difftop%is_inward(i1+idof)) then ! Edn
              xv_diff(idof,k,i,j) = Edn(ak)
            else ! Eup
              xv_diff(idof,k,i,j) = Eup(ak)
            endif
          enddo
        enddo
      enddo
    enddo
    xv_diff = xv_diff / real(solver%difftop%streams, ireals)

    call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

    !Schwarzschild solver returns fluxes as [W/m^2]
    solution%lWm2_dir  = .True.
    solution%lWm2_diff = .True.
    ! and mark solution that it is not up to date
    solution%lchanged         = .True.

    deallocate(Edn)
    deallocate(Eup)

    end associate
  end subroutine

  !> @brief Compute flux divergence, i.e. absorption
  !> @details from gauss divergence theorem, the divergence in the volume is the integral of the flux through the surface
  !> \n we therefore sum up the incoming and outgoing fluxes to compute the divergence
  subroutine calc_flx_div(solver, solution)
    class(t_solver)         :: solver
    type(t_state_container) :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: xediff=>null(),xedir=>null(),xabso=>null()
    real(ireals),pointer,dimension(:)       :: xediff1d=>null(),xedir1d=>null(),xabso1d=>null()

    integer(iintegers) :: offset, isrc, src
    integer(iintegers) :: i, j, k, xinc, yinc
    type(tVec)         :: ledir,lediff ! local copies of vectors, including ghosts
    real(ireals)       :: Volume,Az

    if(allocated(solution%edir)) then
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_flxdiv_edir", ierr); call CHKERR(ierr)
    endif
    call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-show_flxdiv_ediff", ierr); call CHKERR(ierr)

    associate(  atm     => solver%atm, &
                C_dir   => solver%C_dir, &
                C_diff  => solver%C_diff, &
                C_one   => solver%C_one, &
                C_one1  => solver%C_one1)

    if(.not.allocated(solver%abso_scalevec)) then
      allocate(solver%abso_scalevec)
      call VecDuplicate(solution%abso, solver%abso_scalevec, ierr); call CHKERR(ierr)
      call getVecPointer(solver%abso_scalevec, C_one%da, xabso1d, xabso)
      Az = atm%dx * atm%dy

      do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe
          do k=C_one%zs,C_one%ze
            Volume = Az * atm%dz(atmk(atm, k),i,j)
            xabso(i0,k,i,j) = one / Volume
          enddo
        enddo
      enddo

      ! here a special case for icollapse, take the dz of all layers above dynamical grid
      do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe
          Volume = Az * sum(atm%dz(C_one%zs:atmk(atm, C_one%zs),i,j))
          xabso(i0,C_one%zs,i,j) = one / Volume
        enddo
      enddo

      call restoreVecPointer(solver%abso_scalevec, xabso1d ,xabso)
    endif

    if( (solution%lsolar_rad.eqv..False.) .and. lcalc_nca ) then ! if we should calculate NCA (Klinger), we can just return afterwards
      call scale_flx(solver, solution, lWm2=.True.)
      call nca_wrapper(solver, solution%ediff, solution%abso)
      return
    endif

    if(solver%myid.eq.0.and.ldebug) print *,'Calculating flux divergence',solution%lsolar_rad,lcalc_nca
    call VecSet(solution%abso,zero,ierr) ;call CHKERR(ierr)

    ! make sure to bring the fluxes into [W] before the absorption calculation
    call scale_flx(solver, solution, lWm2=.False.)

    ! if there are no 3D layers globally, we should skip the ghost value copying....
    !lhave_no_3d_layer = mpi_logical_and(solver%comm, all(atm%l1d.eqv..True.))
    if(ltwostr_only) then
      if(ldebug.and.solver%myid.eq.0) print *,'lhave_no_3d_layer => will use 1D absorption computation'

      if(solution%lsolar_rad) call getVecPointer(solution%edir, C_dir%da ,xedir1d ,xedir )

      call getVecPointer(solution%ediff, C_diff%da, xediff1d, xediff)
      call getVecPointer(solution%abso, C_one%da, xabso1d, xabso)

      !! calculate absorption by flux divergence

      do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe
          do k=C_one%zs,C_one%ze
            if(solution%lsolar_rad) then
              do isrc = i0, solver%dirtop%dof-1
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + xedir(isrc, k, i, j ) - xedir(isrc , k+i1 , i, j )
              enddo
            endif

            do isrc = 0, solver%difftop%dof-1
              if (solver%difftop%is_inward(i1+isrc)) then
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(isrc, k, i, j) - xediff(isrc, k+1, i, j))
              else
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(isrc, k+1, i, j) - xediff(isrc, k, i, j))
              endif
            enddo
          enddo
        enddo
      enddo

      if(solution%lsolar_rad) call restoreVecPointer(solution%edir, xedir1d, xedir )

      call restoreVecPointer(solution%ediff, xediff1d, xediff)
      call restoreVecPointer(solution%abso, xabso1d ,xabso)

      call VecPointwiseMult(solution%abso, solution%abso, solver%abso_scalevec, ierr); call CHKERR(ierr)
      return
    endif

    if(solution%lsolar_rad) then
      ! Copy ghosted values for direct vec
      call DMGetLocalVector(C_dir%da ,ledir ,ierr)                      ; call CHKERR(ierr)
      call VecSet(ledir ,zero,ierr)                                     ; call CHKERR(ierr)
      call DMGlobalToLocalBegin(C_dir%da ,solution%edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C_dir%da ,solution%edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)
      call getVecPointer(ledir, C_dir%da ,xedir1d ,xedir )
    endif

    ! Copy ghosted values for diffuse vec
    call DMGetLocalVector(C_diff%da,lediff,ierr)                      ; call CHKERR(ierr)
    call VecSet(lediff,zero,ierr)                                     ; call CHKERR(ierr)
    call DMGlobalToLocalBegin(C_diff%da,solution%ediff,ADD_VALUES,lediff,ierr) ; call CHKERR(ierr)
    call DMGlobalToLocalEnd  (C_diff%da,solution%ediff,ADD_VALUES,lediff,ierr) ; call CHKERR(ierr)
    call getVecPointer(lediff, C_diff%da, xediff1d, xediff)

    call getVecPointer(solution%abso, C_one%da, xabso1d, xabso)

    ! calculate absorption by flux divergence
    !DIR$ IVDEP
    do j=C_one%ys,C_one%ye
      do i=C_one%xs,C_one%xe
        do k=C_one%zs,C_one%ze

          ! Divergence = Incoming - Outgoing

          if(atm%l1d(atmk(atm, k),i,j)) then ! one dimensional i.e. twostream
            if(solution%lsolar_rad) then
              do isrc = 0, solver%dirtop%dof-1
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xedir(isrc, k, i, j )  - xedir(isrc , k+i1 , i, j ))
              enddo
            endif

            !xabso(i0,k,i,j) = xabso(i0,k,i,j) + ( xediff(E_up  ,k+1,i  ,j  )  - xediff(E_up  ,k  ,i  ,j  )  )
            !xabso(i0,k,i,j) = xabso(i0,k,i,j) + ( xediff(E_dn  ,k  ,i  ,j  )  - xediff(E_dn  ,k+1,i  ,j  )  )
            do isrc = 0, solver%difftop%dof-1
              if (solver%difftop%is_inward(i1+isrc)) then
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(isrc, k, i, j) - xediff(isrc, k+1, i, j))
              else
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(isrc, k+1, i, j) - xediff(isrc, k, i, j))
              endif
            enddo

          else ! 3D-radiation
            offset = solver%dirtop%dof + solver%dirside%dof*2

            ! direct part of absorption
            if(solution%lsolar_rad) then
              xinc = solver%sun%xinc(k,i,j)
              yinc = solver%sun%yinc(k,i,j)

              src = 0
              do isrc = 0,solver%dirtop%dof-1
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xedir(src, k, i, j) - xedir(src, k+i1, i, j))
                src = src+1
              enddo

              do isrc = 0, solver%dirside%dof-1
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xedir(src, k, i+1-xinc, j) - xedir(src, k, i+xinc, j))
                src = src+1
              enddo

              do isrc = 0, solver%dirside%dof-1
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xedir(src, k, i, j+i1-yinc) - xedir(src, k, i, j+yinc))
                src = src+1
              enddo
            endif

            ! diffuse part of absorption
            src = 0
            do isrc = 0, solver%difftop%dof-1
              if (solver%difftop%is_inward(i1+isrc)) then
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(src, k  , i, j) - xediff(src, k+1, i, j))
              else
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(src, k+1, i, j) - xediff(src, k  , i, j))
              endif
              src = src+1
            enddo

            do isrc = 0, solver%diffside%dof-1
              if (solver%diffside%is_inward(i1+isrc)) then
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(src, k, i, j) - xediff(src, k, i+1, j))
              else
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(src, k, i+1, j) - xediff(src, k, i, j))
              endif
              src = src+1
            enddo

            do isrc = 0, solver%diffside%dof-1
              if (solver%diffside%is_inward(i1+isrc)) then
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(src, k, i, j) - xediff(src, k, i, j+1))
              else
                xabso(i0,k,i,j) = xabso(i0,k,i,j) + (xediff(src, k, i, j+1) - xediff(src, k, i, j))
              endif
              src = src+1
            enddo

            if( ldebug ) then
              if( isnan(xabso(i0,k,i,j)) ) then
                print *,'nan in flxdiv',k,i,j,'::',xabso(i0,k,i,j)
              endif
            endif
          endif ! 1d/3D
        enddo
      enddo
    enddo

    if(solution%lsolar_rad) then
      call restoreVecPointer(ledir, xedir1d, xedir)
      call DMRestoreLocalVector(C_dir%da, ledir, ierr) ; call CHKERR(ierr)
    endif

    call restoreVecPointer(lediff, xediff1d, xediff)
    call DMRestoreLocalVector(C_diff%da, lediff, ierr) ; call CHKERR(ierr)

    call restoreVecPointer(solution%abso, xabso1d ,xabso)
    call VecPointwiseMult(solution%abso, solution%abso, solver%abso_scalevec, ierr); call CHKERR(ierr)

  end associate
end subroutine


!> @brief call PETSc Krylov Subspace Solver
!> @details solve with ksp and save residual history of solver
!> \n -- this may be handy later to decide next time if we have to calculate radiation again
!> \n if we did not get convergence, we try again with standard GMRES and a resetted(zero) initial guess -- if that doesnt help, we got a problem!
subroutine solve(solver, ksp, b, x, iter, ksp_residual_history)
  class(t_solver) :: solver
  type(tKSP) :: ksp
  type(tVec) :: b
  type(tVec) :: x
  integer(iintegers), intent(out) :: iter
  real(ireals), allocatable, intent(inout), optional :: ksp_residual_history(:)

  KSPConvergedReason :: reason

  KSPType :: old_ksp_type

  logical :: lskip_ksp_solve, lflg

  if(solver%myid.eq.0.and.ldebug) print *,'Solving Matrix'

  if(present(ksp_residual_history)) then
    if(.not.allocated( ksp_residual_history) ) allocate(ksp_residual_history(100) )
    ksp_residual_history = -1
    call KSPSetResidualHistory(ksp, ksp_residual_history, 100_iintegers, PETSC_TRUE, ierr); call CHKERR(ierr)
  endif

  lskip_ksp_solve=.False.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
    "-skip_ksp_solve" , lskip_ksp_solve, lflg , ierr) ;call CHKERR(ierr)
  if(lskip_ksp_solve) then
    call VecCopy(b, x, ierr); call CHKERR(ierr)
    return
  endif

  call hegedus_trick(ksp, b, x)
  call KSPSolve(ksp,b,x,ierr) ;call CHKERR(ierr)
  call KSPGetIterationNumber(ksp,iter,ierr) ;call CHKERR(ierr)
  call KSPGetConvergedReason(ksp,reason,ierr) ;call CHKERR(ierr)

  ! if(reason.eq.KSP_DIVERGED_ITS) then
  !   if(solver%myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
  !   return
  ! endif

  if(reason.le.0) then
    if(solver%myid.eq.0) print *,solver%myid,'Resetted initial guess to zero and try again with gmres:'
    call VecSet(x,zero,ierr) ;call CHKERR(ierr)
    call KSPGetType(ksp,old_ksp_type,ierr); call CHKERR(ierr)
    call KSPSetType(ksp,KSPGMRES,ierr) ;call CHKERR(ierr)
    call KSPSetUp(ksp,ierr) ;call CHKERR(ierr)
    call KSPSolve(ksp,b,x,ierr) ;call CHKERR(ierr)
    call KSPGetIterationNumber(ksp,iter,ierr) ;call CHKERR(ierr)
    call KSPGetConvergedReason(ksp,reason,ierr) ;call CHKERR(ierr)

    ! And return to normal solver...
    call KSPSetType(ksp,old_ksp_type,ierr) ;call CHKERR(ierr)
    call KSPSetFromOptions(ksp,ierr) ;call CHKERR(ierr)
    call KSPSetUp(ksp,ierr) ;call CHKERR(ierr)
    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Solver took',iter,'iterations and converged',reason.gt.0,'because',reason
  endif

  if(reason.le.0) then
    if(solver%myid.eq.0) print *,'***** SOLVER did NOT converge :( ********',reason
    call exit()
  endif
end subroutine

!> @brief initialize PETSc Krylov Subspace Solver
!> @details default KSP solver is a FGMRES with BJCAOBI // ILU(1)
!> \n -- the default does however not scale well -- and we configure petsc solvers per commandline anyway
!> \n -- see documentation for details on how to do so
subroutine setup_ksp(atm, ksp, C, A, prefix)
  type(t_atmosphere) :: atm
  type(tKSP), intent(inout), allocatable :: ksp
  type(t_coord) :: C
  type(tMat) :: A
  type(tPC)  :: prec
  logical :: linit

  type(tMatNullSpace) :: nullspace
  type(tVec) :: nullvecs(0)
  character(len=*),optional :: prefix

  real(ireals),parameter :: rtol=1e-5_ireals, rel_atol=1e-4_ireals
  integer(iintegers),parameter  :: maxiter=1000

  integer(mpiint) :: myid, numnodes

  real(ireals) :: atol

  logical,parameter :: lset_geometry=.True.  ! this may be necessary in order to use geometric multigrid
! logical,parameter :: lset_geometry=.False.  ! this may be necessary in order to use geometric multigrid
! logical,parameter :: lset_nullspace=.True. ! set constant nullspace?
  logical,parameter :: lset_nullspace=.False. ! set constant nullspace?

  linit = allocated(ksp)
  if(linit) return

  call mpi_comm_rank(C%comm, myid, ierr)     ; call CHKERR(ierr)
  call mpi_comm_size(C%comm, numnodes, ierr) ; call CHKERR(ierr)

  call imp_allreduce_min(C%comm, &
    rel_atol * real(C%dof*C%glob_xm*C%glob_ym*C%glob_zm, ireals) &
    * count(.not.atm%l1d)/real(size(atm%l1d), ireals), atol)
  atol = max(1e-8_ireals, atol)

  if(myid.eq.0.and.ldebug) &
    print *,'Setup KSP -- tolerances:',rtol,atol,&
      '::',rel_atol,(C%dof*C%glob_xm*C%glob_ym*C%glob_zm),count(.not.atm%l1d),one*size(atm%l1d)

  allocate(ksp)
  call KSPCreate(C%comm,ksp,ierr) ;call CHKERR(ierr)
  if(present(prefix) ) call KSPAppendOptionsPrefix(ksp,trim(prefix),ierr) ;call CHKERR(ierr)

  call KSPSetType(ksp,KSPGMRES,ierr)  ;call CHKERR(ierr)
  call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr) ;call CHKERR(ierr)
  call KSPGetPC  (ksp,prec,ierr)  ;call CHKERR(ierr)
  if(numnodes.eq.0) then
    call PCSetType (prec,PCILU,ierr);call CHKERR(ierr)
  else
    call PCSetType (prec,PCBJACOBI,ierr);call CHKERR(ierr)
  endif

  call KSPSetTolerances(ksp,rtol,atol,PETSC_DEFAULT_REAL,maxiter,ierr);call CHKERR(ierr)

  call KSPSetConvergenceTest(ksp, MyKSPConverged, 0, PETSC_NULL_FUNCTION, ierr)

  call KSPSetOperators(ksp,A,A,ierr) ;call CHKERR(ierr)
  call KSPSetDM(ksp,C%da,ierr) ;call CHKERR(ierr)
  call KSPSetDMActive(ksp,PETSC_FALSE,ierr) ;call CHKERR(ierr)

  call KSPSetUp(ksp,ierr) ;call CHKERR(ierr)

  !if(numnodes.eq.0) then
  !  call PCFactorSetLevels(prec,ilu_default_levels,ierr);call CHKERR(ierr)
  !else
  !  call PCBJacobiGetSubKSP(prec,pcbjac_n_local,pcbjac_iglob,PETSC_NULL_KSP,ierr);call CHKERR(ierr)
  !  if(.not.allocated(pcbjac_ksps)) allocate(pcbjac_ksps(pcbjac_n_local))
  !  call PCBJacobiGetSubKSP(prec,pcbjac_n_local,pcbjac_iglob,pcbjac_ksps,ierr);call CHKERR(ierr)

  !  do isub=1,pcbjac_n_local
  !    call KSPSetType(pcbjac_ksps(isub) ,KSPPREONLY,ierr)              ;call CHKERR(ierr)
  !    call KSPGetPC  (pcbjac_ksps(isub), pcbjac_sub_pc,ierr)        ;call CHKERR(ierr)
  !    call PCSetType (pcbjac_sub_pc, PCILU,ierr)                    ;call CHKERR(ierr)
  !    call PCFactorSetLevels(pcbjac_sub_pc,ilu_default_levels,ierr) ;call CHKERR(ierr)
  !  enddo

  !endif

  if(lset_geometry) call set_coordinates(atm, C,ierr);call CHKERR(ierr)

  if(lset_nullspace) then
    call MatNullSpaceCreate( C%comm, PETSC_TRUE, i0, nullvecs, nullspace, ierr) ; call CHKERR(ierr)
    call MatSetNearNullSpace(A, nullspace, ierr);call CHKERR(ierr)
  endif

  call KSPSetFromOptions(ksp,ierr) ;call CHKERR(ierr)

  linit = .True.
  if(myid.eq.0.and.ldebug) print *,'Setup KSP done'

  contains
    !> @brief define physical coordinates for DMDA to allow for geometric multigrid
    subroutine set_coordinates(atm, C,ierr)
      type(t_atmosphere) :: atm
      type(t_coord) :: C
      PetscErrorCode,intent(out) :: ierr

      Vec :: coordinates
      PetscReal,pointer,dimension(:,:,:,:) :: xv  =>null()
      PetscReal,pointer,dimension(:)       :: xv1d=>null()

      DM  :: coordDA
      integer(iintegers) :: mstart,nstart,pstart,m,n,p, i,j,k

      call DMDASetUniformCoordinates(C%da, zero, one, &
        zero, atm%dx * real(C%glob_xm, ireals), &
        zero, atm%dy * real(C%glob_ym, ireals), ierr );call CHKERR(ierr)
      call DMGetCoordinateDM(C%da,coordDA,ierr);call CHKERR(ierr)
      call DMGetCoordinates(C%da, coordinates,ierr);call CHKERR(ierr)

      call DMDAGetCorners(coordDA, mstart,nstart,pstart,m,n,p,ierr);call CHKERR(ierr)
      !print *,'coordinates',C%xs,C%xe,mstart,nstart,pstart,m,n,p

      call VecGetArrayF90(coordinates,xv1d,ierr) ;call CHKERR(ierr)
      xv(0:2 , mstart:mstart+m-1  , nstart:nstart+n-1   , pstart:pstart+p-1   ) => xv1d

      xv(0,mstart+m-1,:,:) = zero ! surface boundary condition
      do k=pstart,pstart+p-1
        do j=nstart,nstart+n-1
          do i=mstart+m-2, mstart, -1
            xv(0,i,j,k) = xv(0,i+1,j,k) + atm%dz(atmk(atm, i),j,k)
          enddo
        enddo
      enddo

      xv => null()
      call VecRestoreArrayF90(coordinates,xv1d,ierr) ;call CHKERR(ierr)
      xv1d => null()
    end subroutine
  end subroutine


  !> @brief override convergence tests -- the normal KSPConverged returns bad solution if no iterations are needed for convergence
  subroutine MyKSPConverged(ksp, n, rnorm, flag, dummy, ierr)
    ! Input Parameters:
    !    ksp   - iterative context
    !    n     - iteration number
    !    rnorm - 2-norm (preconditioned) residual value (may be estimated)
    !    dummy - optional user-defined monitor context (unused here)
    type(tKSP)        :: ksp
    integer(mpiint)   :: ierr
    integer(iintegers):: n, dummy
    KSPConvergedReason:: flag
    real(ireals)      :: rnorm

    real(ireals)       :: rtol,atol,dtol
    real(ireals),save  :: initial_rnorm
    integer(iintegers) :: maxits

    call KSPGetTolerances(ksp, rtol,atol,dtol,maxits,ierr)

    flag=0
    if(n.eq.0) then
      initial_rnorm=max(tiny(rnorm),rnorm)
      return
    endif

    if(rnorm / initial_rnorm .le.rtol) then
      flag=2
      return
    endif

    if(rnorm.le.atol) then
      flag=3
      return
    endif

    if(n.gt.maxits) then
      flag=-3
      return
    endif

    if(rnorm/initial_rnorm.ge.dtol) then
      flag=-4
      return
    endif

    if(isnan(rnorm)) then
      flag=-9
      return
    endif
    if(dummy.eq.1) print *,'DEBUG: stupid print to remove unused variable warning', dummy
  end subroutine


  !> @brief set solar incoming radiation at Top_of_Atmosphere
  !> @details todo: in case we do not have periodic boundaries, we should shine light in from the side of the domain...
  subroutine setup_incSolar(solver, incSolar,edirTOA)
    class(t_solver), intent(inout)  :: solver
    type(tVec), intent(inout)       :: incSolar
    real(ireals),intent(in)         :: edirTOA

    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    real(ireals) :: fac
    integer(iintegers) :: i,j,src

    fac = edirTOA * solver%atm%dx*solver%atm%dy / real(solver%dirtop%area_divider, ireals)

    call VecSet(incSolar,zero,ierr) ;call CHKERR(ierr)

    call getVecPointer(incSolar, solver%C_dir%da, x1d, x4d)

    do j=solver%C_dir%ys,solver%C_dir%ye
      do i=solver%C_dir%xs,solver%C_dir%xe
        do src=1, solver%dirtop%dof
          x4d(src-1,solver%C_dir%zs,i,j) = fac * solver%sun%costheta(solver%C_dir%zs,i,j)
        enddo
      enddo
    enddo

    call restoreVecPointer(incSolar, x1d, x4d)

    if(solver%myid.eq.0 .and. ldebug) print *,solver%myid,'Setup of IncSolar done', edirTOA

  end subroutine

  !> @brief build direct radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_dir_coeff(solver, sun, A,C)
    class(t_solver)     :: solver
    type(t_suninfo)     :: sun
    type(tMat)          :: A
    type(t_coord)       :: C

    integer(iintegers) :: i,j,k

    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'setup_direct_matrix ...'

    !      call MatZeroEntries(A, ierr) ;call CHKERR(ierr) !TODO necessary?
    !      call mat_set_diagonal(A,C)

    do j=C%ys,C%ye
      do i=C%xs,C%xe
        do k=C%zs,C%ze-1

          if( solver%atm%l1d(atmk(solver%atm,k),i,j) ) then
            call set_eddington_coeff(solver%atm, A,k, i,j)
          else
            call set_pprts_coeff(solver, C, A,k,i,j)
          endif

        enddo
      enddo
    enddo

    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'setup_direct_matrix done'

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)

  contains
    subroutine set_pprts_coeff(solver, C,A,k,i,j)
      class(t_solver)               :: solver
      type(t_coord),intent(in)      :: C
      type(tMat),intent(inout)      :: A
      integer(iintegers),intent(in) :: i,j,k

      MatStencil         :: row(4,C%dof)  ,col(4,C%dof)
      real(irealLUT)     :: v(C%dof**2), norm

      integer(iintegers) :: dst,src, xinc, yinc, isrc, idst

      xinc = sun%xinc(k,i,j)
      yinc = sun%yinc(k,i,j)

      do idst = 1, solver%dirtop%dof
        dst = idst
        row(MatStencil_j,dst) = i
        row(MatStencil_k,dst) = j
        row(MatStencil_i,dst) = k+1
        row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
      enddo

      do idst = 1, solver%dirside%dof
        dst = idst + solver%dirtop%dof
        row(MatStencil_j,dst) = i+xinc
        row(MatStencil_k,dst) = j
        row(MatStencil_i,dst) = k
        row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
      enddo

     do idst = 1, solver%dirside%dof
       dst = idst + solver%dirtop%dof + solver%dirside%dof
       row(MatStencil_j,dst) = i
       row(MatStencil_k,dst) = j+yinc
       row(MatStencil_i,dst) = k
       row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
     enddo

      do isrc = 1, solver%dirtop%dof
        src = isrc
        col(MatStencil_j,src) = i
        col(MatStencil_k,src) = j
        col(MatStencil_i,src) = k
        col(MatStencil_c,src) = src-i1 ! Define transmission towards the lower/upper lid
      enddo

      do isrc = 1, solver%dirside%dof
        src = isrc + solver%dirtop%dof
        col(MatStencil_j,src) = i+1-xinc
        col(MatStencil_k,src) = j
        col(MatStencil_i,src) = k
        col(MatStencil_c,src) = src-i1 ! Define transmission towards the left/right lid
      enddo

     do isrc = 1, solver%dirside%dof
       src = isrc + solver%dirtop%dof + solver%dirside%dof
       col(MatStencil_j,src) = i
       col(MatStencil_k,src) = j+1-yinc
       col(MatStencil_i,src) = k
       col(MatStencil_c,src) = src-i1 ! Define transmission towards the front/back lid
     enddo


      call PetscLogEventBegin(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)
      call get_coeff(solver, &
        solver%atm%kabs(atmk(solver%atm,k),i,j), &
        solver%atm%ksca(atmk(solver%atm,k),i,j), &
        solver%atm%g(atmk(solver%atm,k),i,j), &
        solver%atm%dz(atmk(solver%atm,k),i,j), .True., v, &
        solver%atm%l1d(atmk(solver%atm,k),i,j), &
        [real(sun%symmetry_phi(k,i,j), irealLUT), real(sun%theta(k,i,j), irealLUT)], &
        lswitch_east=xinc.eq.0, lswitch_north=yinc.eq.0)
      call PetscLogEventEnd(solver%logs%get_coeff_dir2dir, ierr); call CHKERR(ierr)

      call MatSetValuesStencil(A, C%dof, row, C%dof, col, real(-v, ireals), INSERT_VALUES, ierr) ;call CHKERR(ierr)

      if(ldebug) then
        do src=1,C%dof
          norm = sum( v(src:C%dof**2:C%dof) )
          if( norm.gt.one+10._ireals*sqrt(epsilon(one)) ) then
            print *,'direct sum(dst==',dst,') gt one',norm
            print *,'direct coeff',norm,'::',v
            call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
          endif
        enddo
      endif

    end subroutine

    subroutine set_eddington_coeff(atm,A,k,i,j)
      type(t_atmosphere), intent(inout) :: atm
      type(tMat),intent(inout)          :: A
      integer(iintegers),intent(in)     :: i,j,k

      MatStencil :: row(4,1), col(4,1)
      real(irealLUT) :: v(1)
      integer(iintegers) :: src

      if(luse_eddington) then
        v = real(atm%a33(atmk(atm,k),i,j), irealLUT)
      else
        call get_coeff(solver, &
          atm%kabs(atmk(atm,k),i,j), &
          atm%ksca(atmk(atm,k),i,j), &
          atm%g(atmk(atm,k),i,j), &
          atm%dz(atmk(atm,k),i,j),.True., v, &
          atm%l1d(atmk(atm,k),i,j), &
          [real(sun%symmetry_phi(k,i,j), irealLUT), real(sun%theta(k,i,j), irealLUT)] )
      endif

      col(MatStencil_j,i1) = i      ; col(MatStencil_k,i1) = j       ; col(MatStencil_i,i1) = k
      row(MatStencil_j,i1) = i      ; row(MatStencil_k,i1) = j       ; row(MatStencil_i,i1) = k+1

      do src=i1,solver%dirtop%dof
        col(MatStencil_c,i1) = src-i1 ! Source may be the upper/lower lid:
        row(MatStencil_c,i1) = src-i1 ! Define transmission towards the lower/upper lid
        call MatSetValuesStencil(A, i1, row, i1, col, real(-v, ireals), INSERT_VALUES, ierr); call CHKERR(ierr)
      enddo
    end subroutine

  end subroutine set_dir_coeff

  !> @brief setup source term for diffuse radiation
  !> @details this is either direct radiation scattered into one of the diffuse coeffs:
  !> \n direct source term is
  !> \n   direct radiation times the dir2diff coeffs
  !> \n or it may be that we have a source term due to thermal emission --
  !> \n   to determine emissivity of box, we use the forward transport coefficients backwards
  !> \n   a la: transmissivity $T = \sum(coeffs)$ and therefore emissivity $E = 1 - T$
  subroutine setup_b(solver, solution, b)
    class(t_solver)         :: solver
    type(t_state_container) :: solution
    type(tVec) :: local_b, b

    real(ireals),pointer,dimension(:,:,:,:) :: xsrc=>null()
    real(ireals),pointer,dimension(:) :: xsrc1d=>null()

    associate(  atm     => solver%atm, &
                C_dir   => solver%C_dir, &
                C_diff  => solver%C_diff)

    if(solver%myid.eq.0.and.ldebug) print *,'src Vector Assembly...'

    call DMGetLocalVector(C_diff%da,local_b,ierr) ;call CHKERR(ierr)
    call VecSet(local_b,zero,ierr) ;call CHKERR(ierr)

    call getVecPointer(local_b, C_diff%da, xsrc1d, xsrc)

    if(solution%lsolar_rad) &
      call set_solar_source(solver, solution%edir)

    if(allocated(atm%planck) ) &
      call set_thermal_source()

    if(solver%myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients ...done'

    call restoreVecPointer(local_b, xsrc1d, xsrc )

    call VecSet(b,zero,ierr) ;call CHKERR(ierr) ! reset global Vec

    call DMLocalToGlobalBegin(C_diff%da, local_b, ADD_VALUES, b,ierr) ;call CHKERR(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
    call DMLocalToGlobalEnd  (C_diff%da, local_b, ADD_VALUES, b,ierr) ;call CHKERR(ierr)

    call DMRestoreLocalVector(C_diff%da,local_b,ierr) ;call CHKERR(ierr)

    if(solver%myid.eq.0.and.ldebug) print *,'src Vector Assembly done'
  end associate
  contains
    subroutine set_thermal_source()

      real(ireals) :: Ax,Ay,Az,emis,b0,b1,btop,bbot,bfac
      real(irealLUT) :: diff2diff1d(4)
      real(irealLUT) :: diff2diff(solver%C_diff%dof**2)
      integer(iintegers) :: k,i,j,src,iside, ak

      associate(  atm     => solver%atm, &
                C_diff  => solver%C_diff)

      if(solver%myid.eq.0.and.ldebug) &
        print *,'Assembly of SRC-Vector ... setting thermal source terms min/max planck', &
        minval(atm%planck), maxval(atm%planck)

      Az = atm%dx * atm%dy / real(solver%difftop%area_divider, ireals)

      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          do k=C_diff%zs,C_diff%ze-1
            ak = atmk(atm,k)
            b0 = atm%planck(ak  ,i,j)
            b1 = atm%planck(ak+1,i,j)

            if( atm%l1d(ak,i,j) ) then

              bfac = pi * Az / real(solver%difftop%streams, ireals)

              if(atm%lcollapse.and.k.eq.i0) then
                btop = atm%Btop(i,j) * bfac
                bbot = atm%Bbot(i,j) * bfac

              else
                if(luse_eddington ) then
                  emis = min(one, max(zero, one-atm%a11(ak,i,j)-atm%a12(ak,i,j)))
                else
                  call get_coeff(solver, &
                    atm%kabs(ak,i,j), &
                    atm%ksca(ak,i,j), &
                    atm%g(ak,i,j), &
                    atm%dz(ak,i,j), &
                    .False., diff2diff1d, &
                    atm%l1d(ak,i,j))

                  emis = one-real(diff2diff1d(1)+diff2diff1d(2), ireals)
                endif

                call B_eff(b1, b0, atm%kabs(ak,i,j), btop)
                call B_eff(b0, b1, atm%kabs(ak,i,j), bbot)
                btop = btop * bfac * emis
                bbot = bbot * bfac * emis
              endif

              do src = 0, solver%difftop%dof-1
                if (solver%difftop%is_inward(i1+src)) then !Edn
                  xsrc(src, k+1, i, j) = xsrc(src, k+1, i, j) + bbot
                else !E_up
                  xsrc(src, k  , i, j) = xsrc(src, k  , i, j) + btop
                endif
              enddo

            else ! Tenstream source terms
              Ax = atm%dy*atm%dz(ak,i,j) / real(solver%diffside%area_divider, ireals)
              Ay = atm%dx*atm%dz(ak,i,j) / real(solver%diffside%area_divider, ireals)

              call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
              call get_coeff(solver, &
                atm%kabs(ak,i,j), &
                atm%ksca(ak,i,j), &
                atm%g(ak,i,j), &
                atm%dz(ak,i,j), &
                .False., diff2diff, &
                atm%l1d(ak,i,j) )
              call PetscLogEventEnd(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)

              call B_eff(b1, b0, atm%kabs(ak,i,j), btop)
              call B_eff(b0, b1, atm%kabs(ak,i,j), bbot)

              src = 0
              bfac = pi * Az / real(solver%difftop%streams, ireals)
              do iside=1,solver%difftop%dof
                emis = one-sum(diff2diff(src+1:C_diff%dof**2:C_diff%dof))
                if (solver%difftop%is_inward(iside) .eqv. .False.) then ! outgoing means Eup
                  xsrc(src, k  , i, j) = xsrc(src, k  , i, j) + btop * bfac * emis
                else
                  xsrc(src, k+1, i, j) = xsrc(src, k+1, i, j) + bbot * bfac * emis
                endif
                src = src+1
              enddo

              bfac = pi * Ax / real(solver%diffside%streams, ireals)
              do iside=1,solver%diffside%dof
                emis = one-sum(diff2diff(src+1:C_diff%dof**2:C_diff%dof))
                if(iside.gt.solver%diffside%dof/2) then ! upward streams
                  emis = btop * emis
                else
                  emis = bbot * emis
                endif
                if (solver%diffside%is_inward(iside) .eqv. .False.) then ! outgoing means towards +x
                  xsrc(src, k, i  , j) = xsrc(src, k, i  , j) + emis * bfac
                else
                  xsrc(src, k, i+1, j) = xsrc(src, k, i+1, j) + emis * bfac
                endif
                src = src+1
              enddo

              bfac = pi * Ay / real(solver%diffside%streams, ireals)
              do iside=1,solver%diffside%dof
                emis = one-sum(diff2diff(src+1:C_diff%dof**2:C_diff%dof))
                if(iside.gt.solver%diffside%dof/2) then ! upward streams
                  emis = btop * emis
                else
                  emis = bbot * emis
                endif
                if (solver%diffside%is_inward(iside) .eqv. .False.) then ! outgoing means towards +y
                  xsrc(src, k, i, j  ) = xsrc(src, k, i, j  ) + emis * bfac
                else
                  xsrc(src, k, i, j+1) = xsrc(src, k, i, j+1) + emis * bfac
                endif
                src = src+1
              enddo

            endif ! 1D or Tenstream?

          enddo ! k
        enddo
      enddo

      ! Thermal emission at surface
      k = C_diff%ze
      if(allocated(atm%Bsrfc)) then
        do j=C_diff%ys,C_diff%ye
          do i=C_diff%xs,C_diff%xe
            do src = 0, solver%difftop%dof-1
              if (.not.solver%difftop%is_inward(i1+src)) then !Eup
                xsrc(src,k,i,j) = xsrc(src,k,i,j) + atm%Bsrfc(i,j) &
                  * Az * (one-atm%albedo(i,j)) * pi / real(solver%difftop%streams, ireals)
              endif
            enddo
          enddo
        enddo
      else
        ak = atmk(atm,k)
        do j=C_diff%ys,C_diff%ye
          do i=C_diff%xs,C_diff%xe
            do src = 0, solver%difftop%dof-1
              if (.not.solver%difftop%is_inward(i1+src)) then !Eup
                xsrc(src,k,i,j) = xsrc(src,k,i,j) + atm%planck(ak,i,j) &
                  * Az * (one-atm%albedo(i,j)) * pi / real(solver%difftop%streams, ireals)
              endif
            enddo
          enddo
        enddo
      endif
    end associate
    end subroutine

    subroutine set_solar_source(solver, edir)
      class(t_solver)     :: solver
      type(tVec)          :: edir

      real(irealLUT)        :: dir2diff(solver%C_dir%dof*solver%C_diff%dof)
      real(ireals)        :: solrad
      integer(iintegers)  :: idof, idofdst, idiff
      integer(iintegers)  :: k, i, j, src, dst

      type(tVec) :: ledir
      real(ireals),pointer,dimension(:,:,:,:) :: xedir=>null()
      real(ireals),pointer,dimension(:)       :: xedir1d=>null()

      logical :: lsun_east,lsun_north

      associate(atm     => solver%atm, &
                C_dir   => solver%C_dir, &
                C_diff  => solver%C_diff, &
                sun     => solver%sun)

      ! Copy ghosted values for direct vec
      call DMGetLocalVector(C_dir%da ,ledir ,ierr)                      ; call CHKERR(ierr)
      call VecSet(ledir ,zero,ierr)                                     ; call CHKERR(ierr)
      call DMGlobalToLocalBegin(C_dir%da ,edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C_dir%da ,edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)

      call getVecPointer(ledir, C_dir%da, xedir1d, xedir)

      if(solver%myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector .. setting solar source', &
        sum(xedir(i0,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye)) / &
        size(xedir(i0,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye))

      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          do k=C_diff%zs,C_diff%ze-1

            if( any (xedir(:,k,i,j) .gt. epsilon(one)) ) then
              if( atm%l1d(atmk(atm,k),i,j) ) then
                dir2diff = zero

                if(luse_eddington ) then
                  ! Only transport the 4 tiles from dir0 to the Eup and Edn
                  do src=1,solver%dirtop%dof

                    do idiff=1,solver%difftop%dof
                      if (solver%difftop%is_inward(idiff)) then
                        ! fetch all diffuse downward fluxes at k+1
                        xsrc(idiff-1,k+1,i,j) = xsrc(idiff-1,k+1,i,j) + &
                          xedir(src-1,k,i,j) * atm%a23(atmk(atm,k),i,j) / real(solver%difftop%streams, ireals)
                      else
                        ! fetch all diffuse upward fluxes at k
                        xsrc(idiff-1,k,i,j) = xsrc(idiff-1,k,i,j) + &
                          xedir(src-1,k,i,j) * atm%a13(atmk(atm,k),i,j) / real(solver%difftop%streams, ireals)
                      endif
                    enddo
                  enddo

                else
                  !call get_coeff(atm%op(atmk(atm,k),i,j), atm%dz(atmk(atm,k),i,j),.False., twostr_coeff, atm%l1d(atmk(atm,k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta])
                  call CHKERR(1_mpiint, 'set solar source only implemented for use with eddington coeff')
                endif


              else ! Tenstream source terms
                lsun_east  = sun%xinc(k,i,j).eq.i0
                lsun_north = sun%yinc(k,i,j).eq.i0

                call PetscLogEventBegin(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)
                call get_coeff(solver, &
                  atm%kabs(atmk(atm,k),i,j), &
                  atm%ksca(atmk(atm,k),i,j), &
                  atm%g(atmk(atm,k),i,j), &
                  atm%dz(atmk(atm,k),i,j), &
                  .False., dir2diff, &
                  atm%l1d(atmk(atm,k),i,j), &
                  [real(sun%symmetry_phi(k,i,j), irealLUT), real(sun%theta(k,i,j),irealLUT)], &
                  lswitch_east=lsun_east, lswitch_north=lsun_north)
                call PetscLogEventEnd(solver%logs%get_coeff_dir2diff, ierr); call CHKERR(ierr)

                dst = 0
                do idofdst=1,solver%difftop%dof
                  src = 1
                  do idof = 1, solver%dirtop%dof
                    solrad = xedir( src-1, k , i , j )

                    if (solver%difftop%is_inward(idofdst)) then
                      xsrc(dst,k+1,i,j)= xsrc(dst,k+1,i,j) + solrad * dir2diff((dst)*C_dir%dof+src)
                    else
                      xsrc(dst,k  ,i,j)= xsrc(dst,k  ,i,j) + solrad * dir2diff((dst)*C_dir%dof+src)
                    endif
                    src = src+1
                  enddo

                  do idof = 1, solver%dirside%dof
                    solrad = xedir( src-1, k , i+i1-sun%xinc(k,i,j) , j )
                    if (solver%difftop%is_inward(idofdst)) then
                      xsrc(dst,k+1,i,j)= xsrc(dst,k+1,i,j) + solrad * dir2diff((dst)*C_dir%dof+src)
                    else
                      xsrc(dst,k  ,i,j)= xsrc(dst,k  ,i,j) + solrad * dir2diff((dst)*C_dir%dof+src)
                    endif
                    src = src+1
                  enddo

                  do idof = 1, solver%dirside%dof
                    solrad = xedir( src-1 , k , i , j+i1-sun%yinc(k,i,j) )
                    if (solver%difftop%is_inward(idofdst)) then
                      xsrc(dst,k+1,i,j)= xsrc(dst,k+1,i,j) + solrad * dir2diff((dst)*C_dir%dof+src)
                    else
                      xsrc(dst,k  ,i,j)= xsrc(dst,k  ,i,j) + solrad * dir2diff((dst)*C_dir%dof+src)
                    endif
                    src = src+1
                  enddo
                  dst = dst+1
                enddo

                do idofdst = 1, solver%diffside%dof
                  src = 1
                  do idof = 1, solver%dirtop%dof
                    solrad = xedir( src-1 , k , i , j )
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst, k, i+1, j) = xsrc(dst, k, i+1, j) + solrad * dir2diff((dst)*C_dir%dof + src)
                    else
                      xsrc(dst, k, i  , j) = xsrc(dst, k, i  , j) + solrad * dir2diff((dst)*C_dir%dof + src)
                    endif
                    src = src+1
                  enddo

                  do idof = 1, solver%dirside%dof
                    solrad = xedir(src-1, k, i+i1-sun%xinc(k,i,j) , j )
                    if (solver%diffside%is_inward(idofdst)) then
                      xsrc(dst, k, i+1, j) = xsrc(dst, k, i+1, j) + solrad * dir2diff((dst)*C_dir%dof + src)
                    else
                      xsrc(dst, k, i  , j) = xsrc(dst, k, i  , j) + solrad * dir2diff((dst)*C_dir%dof + src)
                    endif
                    src = src+1
                  enddo

                  do idof = 1, solver%dirside%dof
                    solrad = xedir( src-1 , k , i , j+i1-sun%yinc(k,i,j) )
                    if (solver%diffside%is_inward(idofdst)) then
                      xsrc(dst, k, i+1, j) = xsrc(dst, k, i+1, j) + solrad * dir2diff((dst)*C_dir%dof + src)
                    else
                      xsrc(dst, k, i  , j) = xsrc(dst, k, i  , j) + solrad * dir2diff((dst)*C_dir%dof + src)
                    endif
                    src = src+1
                  enddo
                  dst = dst+1
                enddo

                do idofdst = 1, solver%diffside%dof
                  src = 1
                  do idof = 1, solver%dirtop%dof
                    solrad = xedir(src-1, k, i, j)
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst, k, i, j+1) = xsrc(dst, k, i, j+1) + solrad * dir2diff((dst)*C_dir%dof + src)
                    else
                      xsrc(dst, k, i, j  ) = xsrc(dst, k, i, j  ) + solrad * dir2diff((dst)*C_dir%dof + src)
                    endif
                    src = src+1
                  enddo

                  do idof = 1, solver%dirside%dof
                    solrad = xedir(src-1, k, i+i1-sun%xinc(k,i,j) , j)
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst, k, i, j+1) = xsrc(dst, k, i, j+1) + solrad * dir2diff((dst)*C_dir%dof + src)
                    else
                      xsrc(dst, k, i, j  ) = xsrc(dst, k, i, j  ) + solrad * dir2diff((dst)*C_dir%dof + src)
                    endif
                    src = src+1
                  enddo

                  do idof = 1, solver%dirside%dof
                    solrad = xedir( src-1 , k , i , j+i1-sun%yinc(k,i,j) )
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst, k, i, j+1) = xsrc(dst, k, i, j+1) + solrad * dir2diff((dst)*C_dir%dof + src)
                    else
                      xsrc(dst, k, i, j  ) = xsrc(dst, k, i, j  ) + solrad * dir2diff((dst)*C_dir%dof + src)
                    endif
                    src = src+1
                  enddo
                  dst = dst+1
                enddo

                if(ldebug) then
                  do src=1,C_dir%dof
                  if(sum(dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof)) .gt. one .or. &
                    sum(dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof)) .lt. zero   ) &
                    print *,'DEBUG Found dir2diff gt one:',src,'::', &
                      sum(dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof)),&
                      ':',dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof) ,&
                      '   :::::::     ', dir2diff
                  enddo
                endif

              endif ! 1D or Tenstream?
            endif ! if solar

          enddo
        enddo
      enddo

      ! Ground Albedo reflecting direct radiation, the diffuse part is considered by the solver(Matrix)
      k = C_diff%ze
      do j = C_diff%ys, C_diff%ye
        do i = C_diff%xs, C_diff%xe
          do dst = 0, solver%difftop%dof-1
            if (.not. solver%difftop%is_inward(dst+1)) then
              do src = 0, solver%dirtop%dof-1
                xsrc(dst,k,i,j) = xsrc(dst, k, i, j) + xedir(src,k,i,j) * atm%albedo(i,j) &
                  / real(solver%difftop%streams, ireals)
              enddo
            endif
          enddo
        enddo
      enddo

      call restoreVecPointer(ledir, xedir1d, xedir )
      call DMRestoreLocalVector(C_dir%da, ledir, ierr); call CHKERR(ierr)
      end associate
    end subroutine
  end subroutine setup_b


  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(solver, kabs, ksca, g, dz, ldir, coeff, &
      lone_dimensional, angles, lswitch_east, lswitch_north)
    class(t_solver), intent(inout)    :: solver
    real(ireals),intent(in)           :: kabs, ksca, g, dz
    logical,intent(in)                :: ldir
    real(irealLUT),intent(out)        :: coeff(:)

    logical,intent(in)                :: lone_dimensional
    real(irealLUT),intent(in),optional:: angles(2)
    logical,intent(in),optional       :: lswitch_east, lswitch_north

    real(irealLUT) :: aspect_zx, tauz, w0
    integer(mpiint) :: ierr

    aspect_zx = real(dz / solver%atm%dx, irealLUT)
    aspect_zx = max(solver%OPP%OPP_LUT%diffconfig%dims(3)%vrange(1),aspect_zx) !DEBUG

    tauz = max(solver%OPP%OPP_LUT%diffconfig%dims(1)%vrange(1), &
      min(solver%OPP%OPP_LUT%diffconfig%dims(1)%vrange(2), real((kabs+ksca) * dz, irealLUT)))
    w0 = real(ksca / max(kabs+ksca, epsilon(kabs)), irealLUT)
    w0 = max(solver%OPP%OPP_LUT%diffconfig%dims(2)%vrange(1), &
      min(solver%OPP%OPP_LUT%diffconfig%dims(2)%vrange(2), w0))

    if(lone_dimensional) then
      call CHKERR(1_mpiint, 'currently, we dont support using LUT Twostream for l1d layers')
      !call OPP_1_2%get_coeff (aspect, tauz, w0, g,ldir,coeff,angles)
    else
      call solver%OPP%get_coeff(tauz, w0, real(g, irealLUT), aspect_zx, ldir, coeff, ierr, &
        angles, lswitch_east, lswitch_north); call CHKERR(ierr)
    endif
  end subroutine


  !> @brief nca wrapper to call NCA of Carolin Klinger
  !> @details This is supposed to work on a 1D solution which has to be calculated beforehand
  !> \n the wrapper copies fluxes and optical properties on one halo and then gives that to NCA
  !> \n the result is the 3D approximation of the absorption, considering neighbouring information
  subroutine nca_wrapper(solver, ediff, abso)
    use m_ts_nca, only : ts_nca
    class(t_solver) :: solver
    type(tVec) :: ediff, abso
    type(tVec) :: gnca ! global nca vector
    type(tVec) :: lnca ! local nca vector with ghost values -- in dimension 0 and 1 are fluxes followed by dz,planck,kabs

    real(ireals),pointer,dimension(:,:,:,:) :: xv  =>null()
    real(ireals),pointer,dimension(:)       :: xv1d=>null()
    real(ireals),pointer,dimension(:,:,:,:) :: xvlnca  =>null(), xvgnca  =>null()
    real(ireals),pointer,dimension(:)       :: xvlnca1d=>null(), xvgnca1d=>null()
    real(ireals),pointer,dimension(:,:,:,:) :: xhr  =>null()
    real(ireals),pointer,dimension(:)       :: xhr1d=>null()
    integer(iintegers) :: k

    integer(iintegers),parameter :: E_up=i0, E_dn=i1, idz=i2, iplanck=i3, ikabs=i4, ihr=i5

    associate(  atm     => solver%atm, &
                C_one   => solver%C_one, &
                C_diff  => solver%C_diff)

     !call CHKERR(1_mpiint, 'nca_wrapper not implemented')
     !print *, 'DEBUG: Stupid print statement to prevent unused compiler warnings', ediff, abso
     if(C_diff%dof.lt.i6) call CHKERR(1_mpiint, 'For NCA, need a solver with at least 6 diffuse streams to copy over some data')

     ! put additional values into a local ediff vec .. TODO: this is a rather dirty hack but is straightforward

     ! get ghost values for dz, planck, kabs and fluxes, ready to give it to NCA
     call DMGetGlobalVector(solver%C_diff%da ,gnca ,ierr) ; call CHKERR(ierr)

     call getVecPointer(gnca ,solver%C_diff%da ,xvgnca1d, xvgnca)
     xvgnca(  idz    , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%dz
     xvgnca(  iplanck, C_diff%zs:C_diff%ze  , C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%planck
     xvgnca(  ikabs  , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%kabs


     ! Copy Edn and Eup to local convenience vector
     call getVecPointer(ediff ,C_diff%da ,xv1d, xv)
     xvgnca( E_up,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_up,:,:,:)
     xvgnca( E_dn,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_dn,:,:,:)
     call restoreVecPointer(ediff, xv1d, xv)

     call restoreVecPointer(gnca, xvgnca1d, xvgnca )


     ! retrieve ghost values into l(ocal) nca vec
     call DMGetLocalVector (C_diff%da ,lnca ,ierr) ; call CHKERR(ierr)
     call VecSet(lnca, zero, ierr); call CHKERR(ierr)

     call DMGlobalToLocalBegin(C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)
     call DMGlobalToLocalEnd  (C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)

     call DMRestoreGlobalVector(C_diff%da, gnca, ierr); call CHKERR(ierr)

     ! call NCA
     call getVecPointer(lnca ,C_diff%da ,xvlnca1d, xvlnca)

     call ts_nca( atm%dx, atm%dy,                    &
       xvlnca(   idz        , : , : , :), &
       xvlnca(   iplanck    , : , : , :), &
       xvlnca(   ikabs      , : , : , :), &
       xvlnca(   E_dn       , : , : , :), &
       xvlnca(   E_up       , : , : , :), &
       xvlnca(   ihr        , : , : , :))


     ! return absorption
     call getVecPointer( abso, C_one%da ,xhr1d, xhr)

     do k=C_one%zs,C_one%ze
       xhr(i0,k,:,:) = xvlnca( ihr , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) / &
         xvlnca( idz , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye)
     enddo
     call restoreVecPointer(abso, xhr1d, xhr )

     !return convenience vector that holds optical properties
     call restoreVecPointer(lnca, xvlnca1d, xvlnca )
     call DMRestoreLocalVector(C_diff%da, lnca, ierr); call CHKERR(ierr)
   end associate
  end subroutine


  !> @brief build diffuse radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_diff_coeff(solver, A,C)
    class(t_solver) :: solver
    type(tMat)      :: A
    type(t_coord)   :: C

    integer(iintegers) :: i,j,k

    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Setting coefficients for diffuse Light'

    !      call MatZeroEntries(A, ierr) ;call CHKERR(ierr) !TODO necessary?
    !      call mat_set_diagonal(A,C)

    do j=C%ys,C%ye
      do i=C%xs,C%xe
        do k=C%zs,C%ze-1

          if( solver%atm%l1d(atmk(solver%atm, k),i,j) ) then
            call set_eddington_coeff(solver%atm, A, k,i,j)
          else
            call set_pprts_coeff(solver, C, A, k,i,j, ierr); call CHKERR(ierr)
          endif

        enddo
      enddo
    enddo

    call set_albedo_coeff(solver, C, A )

    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'setup_diffuse_matrix done'

    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Final diffuse Matrix Assembly:'
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)


    call PetscObjectViewFromOptions(solver%Mdiff, PETSC_NULL_MAT, "-show_Mdiff", ierr); call CHKERR(ierr)
  contains
    subroutine set_pprts_coeff(solver,C,A,k,i,j,ierr)
      class(t_solver)               :: solver
      type(t_coord),intent(in)      :: C
      type(tMat),intent(inout)      :: A
      integer(iintegers),intent(in) :: k,i,j
      integer(mpiint),intent(out)   :: ierr

      MatStencil :: row(4,0:C%dof-1)  ,col(4,0:C%dof-1)
      real(irealLUT) :: v(C%dof**2),norm

      integer(iintegers) :: dst,src,idof

      src = 0
      do idof=1, solver%difftop%dof
        if (solver%difftop%is_inward(idof)) then
          col(MatStencil_j,src) = i
          col(MatStencil_k,src) = j
          col(MatStencil_i,src) = k
          col(MatStencil_c,src) = src
        else
          col(MatStencil_j,src) = i
          col(MatStencil_k,src) = j
          col(MatStencil_i,src) = k+1
          col(MatStencil_c,src) = src
        endif
        src = src + 1
      enddo

      do idof=1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          col(MatStencil_j,src) = i
          col(MatStencil_k,src) = j
          col(MatStencil_i,src) = k
          col(MatStencil_c,src) = src
        else
          col(MatStencil_j,src) = i+1
          col(MatStencil_k,src) = j
          col(MatStencil_i,src) = k
          col(MatStencil_c,src) = src
        endif
        src = src + 1
      enddo

      do idof=1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          col(MatStencil_j,src) = i
          col(MatStencil_k,src) = j
          col(MatStencil_i,src) = k
          col(MatStencil_c,src) = src
        else
          col(MatStencil_j,src) = i
          col(MatStencil_k,src) = j+1
          col(MatStencil_i,src) = k
          col(MatStencil_c,src) = src
        endif
        src = src + 1
      enddo

      dst = 0
      do idof=1, solver%difftop%dof
        if (solver%difftop%is_inward(idof)) then
          row(MatStencil_j,dst) = i
          row(MatStencil_k,dst) = j
          row(MatStencil_i,dst) = k+1
          row(MatStencil_c,dst) = dst
        else
          row(MatStencil_j,dst) = i
          row(MatStencil_k,dst) = j
          row(MatStencil_i,dst) = k
          row(MatStencil_c,dst) = dst
        endif
        dst = dst + 1
      enddo

      do idof=1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          row(MatStencil_j,dst) = i+1
          row(MatStencil_k,dst) = j
          row(MatStencil_i,dst) = k
          row(MatStencil_c,dst) = dst
        else
          row(MatStencil_j,dst) = i
          row(MatStencil_k,dst) = j
          row(MatStencil_i,dst) = k
          row(MatStencil_c,dst) = dst
        endif
        dst = dst + 1
      enddo

      do idof=1, solver%diffside%dof
        if (solver%diffside%is_inward(idof)) then
          row(MatStencil_j,dst) = i
          row(MatStencil_k,dst) = j+1
          row(MatStencil_i,dst) = k
          row(MatStencil_c,dst) = dst
        else
          row(MatStencil_j,dst) = i
          row(MatStencil_k,dst) = j
          row(MatStencil_i,dst) = k
          row(MatStencil_c,dst) = dst
        endif
        dst = dst + 1
      enddo

      call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
      call get_coeff(solver, &
        solver%atm%kabs(atmk(solver%atm, k),i,j), &
        solver%atm%ksca(atmk(solver%atm, k),i,j), &
        solver%atm%g(atmk(solver%atm, k),i,j), &
        solver%atm%dz(atmk(solver%atm, k),i,j), &
        .False., v, solver%atm%l1d(atmk(solver%atm, k),i,j))
      call PetscLogEventEnd(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)

      call MatSetValuesStencil(A,C%dof, row,C%dof, col , real(-v, ireals) ,INSERT_VALUES,ierr) ;call CHKERR(ierr)

      if(ldebug) then
        do src=1,C%dof
          norm = sum( v(src:C%dof**2:C%dof))
          if(norm.gt.one+10._ireals*epsilon(one)) then
            if(norm.gt.one+10._ireals*sqrt(epsilon(one))) then
              print *,'diffuse sum(src==',src,') gt one',norm,'=>', v(src:C%dof**2:C%dof)
              print *,'get_coeff', solver%atm%kabs(atmk(solver%atm, k),i,j), &
                solver%atm%ksca(atmk(solver%atm, k),i,j), &
                solver%atm%g(atmk(solver%atm, k),i,j), &
                solver%atm%dz(atmk(solver%atm, k),i,j), &
                .False., solver%atm%l1d(atmk(solver%atm, k),i,j), '=>', v
              call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
            endif
          endif
        enddo
      endif

    end subroutine

    subroutine set_eddington_coeff(atm,A,k,i,j)
      type(t_atmosphere)            :: atm
      type(tMat),intent(inout)      :: A
      integer(iintegers),intent(in) :: k,i,j
      integer(iintegers) :: src, dst, idof

      MatStencil   :: row(4,0:solver%difftop%dof-1)  ,col(4,0:solver%difftop%dof-1)
      real(ireals) :: v(solver%difftop%dof**2)

      do idof=1, solver%difftop%dof
        src = idof-1
        if (solver%difftop%is_inward(idof)) then
          col(MatStencil_j,src) = i
          col(MatStencil_k,src) = j
          col(MatStencil_i,src) = k
          col(MatStencil_c,src) = src
        else
          col(MatStencil_j,src) = i
          col(MatStencil_k,src) = j
          col(MatStencil_i,src) = k+1
          col(MatStencil_c,src) = src
        endif
      enddo

      do idof=1, solver%difftop%dof
        dst = idof-1
        if (solver%difftop%is_inward(idof)) then
          row(MatStencil_j,dst) = i
          row(MatStencil_k,dst) = j
          row(MatStencil_i,dst) = k+1
          row(MatStencil_c,dst) = dst
        else
          row(MatStencil_j,dst) = i
          row(MatStencil_k,dst) = j
          row(MatStencil_i,dst) = k
          row(MatStencil_c,dst) = dst
        endif
      enddo

      ! for each destination, find all transmission coeffs
      v = zero
      do dst = 0, solver%difftop%dof-1
          do src = 0, solver%difftop%dof-1
              if(col(MatStencil_i,src).eq.row(MatStencil_i,dst)) then ! for reflection, has to be the same k layers

                if(src.ne.inv_dof(dst)) cycle ! in 1D has to be the inverse stream
                v(i1 + dst*solver%difftop%dof + src) = atm%a12(atmk(atm,k),i,j)
                !print *,i,j,k,'setting r ',itoa(src)//' (k='//itoa(col(MatStencil_i,src))//') ->', &
                !  itoa(dst)//' (k='//itoa(row(MatStencil_i,dst))//')'// &
                !  ':', i1 + dst*solver%difftop%dof + src, v(i1 + dst*solver%difftop%dof + src), &
                !  'invdof',src,dst,inv_dof(dst)
              else

                if(src.ne.dst) cycle ! in 1D has to be the same
                v(i1 + dst*solver%difftop%dof + src) = atm%a11(atmk(atm,k),i,j)
                !print *,i,j,k,'setting t ',itoa(src)//' (k='//itoa(col(MatStencil_i,src))//') ->', &
                !  itoa(dst)//' (k='//itoa(row(MatStencil_i,dst))//') :', i1 + dst*solver%difftop%dof + src, v(i1 + dst*solver%difftop%dof + src)
              endif ! which k-lev
          enddo
      enddo

      call MatSetValuesStencil(A, solver%difftop%dof, row, solver%difftop%dof, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)
    end subroutine
    pure function inv_dof(dof) ! returns the dof that is the same stream but the opposite direction
      integer(iintegers), intent(in) :: dof
      integer(iintegers) :: inv_dof, inc
      if(solver%difftop%is_inward(1)) then ! starting with downward streams
        inc = 1
      else
        inc = -1
      endif
      if(solver%difftop%is_inward(i1+dof)) then ! downward stream
        inv_dof = dof + inc
      else
        inv_dof = dof - inc
      endif
    end function

    !> @brief insert lower boundary condition, i.e. diffuse reflection of downward radiation
    subroutine set_albedo_coeff(solver,C,A)
      class(t_solver), intent(in) :: solver
      type(t_coord), intent(in)   :: C
      type(tMat), intent(inout)   :: A

      MatStencil :: row(4,1) ,col(4,1)
      integer(iintegers) :: i, j, src, dst

      ! Set surface albedo values
      col(MatStencil_i,i1) = C%ze
      row(MatStencil_i,i1) = C%ze

      do j=C%ys,C%ye
        row(MatStencil_k,i1) = j
        col(MatStencil_k,i1) = j

        do i=C%xs,C%xe
          row(MatStencil_j,i1) = i
          col(MatStencil_j,i1) = i

          do dst = 1, solver%difftop%dof
            if (.not.solver%difftop%is_inward(dst)) then
              do src = 1, solver%difftop%dof
                if (solver%difftop%is_inward(src)) then
                  col(MatStencil_c,i1) = src-1
                  row(MatStencil_c,i1) = dst-1
                  !print *,solver%myid, 'i '//itoa(i)//' j '//itoa(j), ' Setting albedo for dst '//itoa(dst)//' src '//itoa(src)
                  call MatSetValuesStencil(A, i1, row, i1, col, &
                    [-solver%atm%albedo(i,j) / real(solver%difftop%streams, ireals)], &
                    INSERT_VALUES, ierr) ;call CHKERR(ierr)
                endif
              enddo
            endif
          enddo

        enddo
      enddo
    end subroutine
  end subroutine

  subroutine pprts_get_result(solver, redn, reup, rabso, redir, opt_solution_uid )
    class(t_solver) :: solver
    real(ireals),dimension(:,:,:),intent(inout),allocatable          :: redn,reup,rabso
    real(ireals),dimension(:,:,:),intent(inout),allocatable,optional :: redir
    integer(iintegers),optional,intent(in) :: opt_solution_uid

    integer(iintegers)  :: uid, iside
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    uid = get_arg(0_iintegers, opt_solution_uid)

    if(ldebug .and. solver%myid.eq.0) print *,'calling pprts_get_result',present(redir),'for uid',uid

    if(solver%solutions(uid)%lchanged) &
      call CHKERR(1_mpiint, 'tried to get results from unrestored solution -- call restore_solution first')

    if(present(opt_solution_uid)) then
      if(.not.solver%solutions(uid)%lWm2_diff) &
        call CHKERR(1_mpiint, 'solution vecs for diffuse radiation are not in W/m2 ... this is not what I expected')
    endif

    if(allocated(redn)) then
      if(.not.all(shape(redn).eq.[solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym])) then
        print *,'Shape redn', shape(redn), 'vs', [solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym]
        call CHKERR(1_mpiint, 'the shape of edn result array which you provided does not conform to output size. '// &
          'Either call with unallocated object or make sure it has the correct size')
      endif
    else
      allocate(redn(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    endif

    if(allocated(reup)) then
      if(.not.all(shape(reup).eq.[solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym])) then
        print *,'Shape reup', shape(reup), 'vs', [solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym]
        call CHKERR(1_mpiint, 'the shape of eup result array which you provided does not conform to output size. '// &
          'Either call with unallocated object or make sure it has the correct size')
      endif
    else
      allocate(reup(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    endif

    if(allocated(rabso)) then
      if(.not.all(shape(rabso).eq.[solver%C_one%zm, solver%C_one%xm, solver%C_one%ym])) then
        print *,'Shape rabso', shape(rabso), 'vs', [solver%C_one%zm, solver%C_one%xm, solver%C_one%ym]
        call CHKERR(1_mpiint, 'the shape of absorption result array which you provided does not conform to output size. '//&
          'Either call with unallocated object or make sure it has the correct size')
      endif
    else
      allocate(rabso(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    endif

    if(present(redir)) then
      if(allocated(redir)) then
        if(.not.all(shape(redir).eq.[solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym])) then
          print *,'Shape redir', shape(redir), 'vs', [solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym]
          call CHKERR(1_mpiint, 'pprts_get_result :: you should not call it with an allocated redir array')
        endif
      else
        allocate(redir(solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym))
      endif

      if( .not. solver%solutions(uid)%lsolar_rad ) then
        print *,'Hey, You called pprts_get_result for uid '//itoa(uid)// &
          ' and provided an array for direct radiation.'// &
          ' However in this particular band we haven`t computed direct radiation.'// &
          ' I will return with edir=0 but are you sure this is what you intended?'
        redir = zero
      else
        if(.not.solver%solutions(uid)%lWm2_dir) &
          call CHKERR(1_mpiint, 'tried to get result from a result vector(dir) which is not in [W/m2]')

        call getVecPointer(solver%solutions(uid)%edir, solver%C_dir%da, x1d, x4d)
        ! average of direct radiation of all fluxes through top faces
        redir = sum(x4d(0:solver%dirtop%dof-1, :, :, :), dim=1) / real(solver%dirtop%area_divider, ireals)
        call restoreVecPointer(solver%solutions(uid)%edir, x1d, x4d)

        if(ldebug) then
          if(solver%myid.eq.0) print *,'Edir vertically first column',redir(:, lbound(redir,2), lbound(redir,3))
          if(any(redir.lt.-one)) then
            print *,'Found direct radiation smaller than 0 in dir result... that should not happen',minval(redir)
            call CHKERR(1_mpiint)
          endif
        endif
      endif
    endif

    if(.not.solver%solutions(uid)%lWm2_diff) &
      call CHKERR(1_mpiint, 'tried to get result from a result vector(diff) which is not in [W/m2]')


    redn = zero
    reup = zero
    call getVecPointer(solver%solutions(uid)%ediff, solver%C_diff%da, x1d, x4d)
    do iside=1,solver%difftop%dof
      if(solver%difftop%is_inward(iside)) then
        redn = redn + x4d(iside-1, :, :, :) / real(solver%difftop%area_divider, ireals)
      else
        reup = reup + x4d(iside-1, :, :, :) / real(solver%difftop%area_divider, ireals)
      endif
    enddo
    call restoreVecPointer(solver%solutions(uid)%ediff,x1d,x4d)

    if(solver%myid.eq.0 .and. ldebug .and. present(redir)) &
      print *,'mean surface Edir',meanval(redir(ubound(redir,1),:,:))
    if(solver%myid.eq.0 .and. ldebug) print *,'mean surface Edn',meanval(redn(ubound(redn,1), :,:))
    if(solver%myid.eq.0 .and. ldebug) print *,'mean surface Eup',meanval(reup(ubound(reup,1), :,:))

    if(ldebug .and. solver%solutions(uid)%lsolar_rad) then
      if(any(redn.lt.-one)) then
        print *,'Found radiation smaller than 0 in edn result... that should not happen',minval(redn)
        call exit(1)
      endif
      if(any(reup.lt.-one)) then
        print *,'Found radiation smaller than 0 in eup result... that should not happen',minval(reup)
        call exit(1)
      endif
    endif

    call getVecPointer(solver%solutions(uid)%abso, solver%C_one%da, x1d, x4d)
    rabso = x4d(i0,:,:,:)
    call restoreVecPointer(solver%solutions(uid)%abso,x1d,x4d)
    if(solver%myid.eq.0 .and. ldebug) print *,'get_result done'
  end subroutine

  subroutine pprts_get_result_toZero(solver, gedn, geup, gabso, gedir, opt_solution_uid)
    ! after solving equations -- retrieve the results for edir,edn,eup and absorption
    ! only zeroth node gets the results back.
    class(t_solver)   :: solver

    real(ireals),intent(inout),dimension(:,:,:),allocatable :: gedn
    real(ireals),intent(inout),dimension(:,:,:),allocatable :: geup
    real(ireals),intent(inout),dimension(:,:,:),allocatable :: gabso
    real(ireals),intent(inout),dimension(:,:,:),allocatable,optional :: gedir
    integer(iintegers),optional,intent(in) :: opt_solution_uid

    real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso

    call PetscLogEventBegin(solver%logs%scatter_to_Zero, ierr); call CHKERR(ierr)
    if(solver%myid.eq.0) then
      call check_arr_size(solver%C_one_atm1, gedn )
      call check_arr_size(solver%C_one_atm1, geup )
      call check_arr_size(solver%C_one_atm , gabso)
      if(present(gedir)) call check_arr_size(solver%C_one_atm1, gedir)
    endif

    if(present(gedir)) then
      call pprts_get_result(solver,redn,reup,rabso,redir=redir,opt_solution_uid=opt_solution_uid)
      call gather_all_toZero(solver%C_one_atm1, redir, gedir)
    else
      call pprts_get_result(solver,redn,reup,rabso,opt_solution_uid=opt_solution_uid)
    endif

    call gather_all_toZero(solver%C_one_atm1, redn , gedn )
    call gather_all_toZero(solver%C_one_atm1, reup , geup )
    call gather_all_toZero(solver%C_one_atm , rabso, gabso)

    if(solver%myid.eq.0 .and. ldebug) then
      print *,'Retrieving results:'
      if(present(gedir)) print *,sum(gedir)/size(gedir)
      print *,sum(gedn) /size(gedn)
      print *,sum(geup) /size(geup)
      print *,sum(gabso)/size(gabso)
    endif
    call PetscLogEventEnd(solver%logs%scatter_to_Zero, ierr); call CHKERR(ierr)

  contains
    subroutine check_arr_size(C,inp)
      type(t_coord),intent(in) :: C
      real(ireals),intent(in), allocatable :: inp(:,:,:)
      if(allocated(inp)) then
        if(size(inp).ne.C%dof*C%glob_zm*C%glob_xm*C%glob_ym) then
          print *,'pprts_get_result_toZero was called with an already allocated output array but it has the wrong size.'
          print *,'while I could just re-allocate it, this may not be just what you intended. Please do that yourself.'
          print *,'Size of global Dimensions of simulation:',C%dof*C%glob_zm*C%glob_xm*C%glob_ym, '.vs. your input:',size(inp)
          call CHKERR(1_mpiint, 'pprts_get_result_toZero :: should not be called with already allocated array with wrong size')
        endif
      endif
    end subroutine
  end subroutine

  subroutine gather_all_toZero(C, inp, outp)
    type(t_coord),intent(in) :: C
    real(ireals),intent(in), allocatable :: inp(:,:,:) ! local array from get_result
    real(ireals),intent(inout),allocatable :: outp(:,:,:) ! global sized array on rank 0

    type(tVec) :: vec, lvec_on_zero

    integer(mpiint) :: myid, ierr
    call mpi_comm_rank(C%comm, myid, ierr)

    if(ldebug) then
      print *,myid,'exchange_var',allocated(inp), allocated(outp)
      print *,myid,'exchange_var shape',shape(inp)
    endif

    call DMGetGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)
    call f90VecToPetsc(inp, C%da, vec)
    call petscGlobalVecToZero(vec, C%da, lvec_on_zero)
    call DMRestoreGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)

    if(myid.eq.0) then
      call petscVecToF90(lvec_on_zero, C%da, outp, opt_l_only_on_rank0=.True.)
      call VecDestroy(lvec_on_zero, ierr); call CHKERR(ierr)
    endif
  end subroutine

  subroutine destroy_pprts(solver, lfinalizepetsc)
    class(t_solver)   :: solver
    logical,optional :: lfinalizepetsc
    logical :: lfinalize
    integer(iintegers) :: uid
    lfinalize = get_arg(.False., lfinalizepetsc)

    if(solver%linitialized) then
      if(allocated(solver%ksp_solar_dir)) then
        call KSPDestroy(solver%ksp_solar_dir, ierr) ;call CHKERR(ierr)
        deallocate(solver%ksp_solar_dir)
      endif
      if(allocated(solver%ksp_solar_diff)) then
        call KSPDestroy(solver%ksp_solar_diff, ierr) ;call CHKERR(ierr)
        deallocate(solver%ksp_solar_diff)
      endif
      if(allocated(solver%ksp_thermal_diff)) then
        call KSPDestroy(solver%ksp_thermal_diff, ierr) ;call CHKERR(ierr)
        deallocate(solver%ksp_thermal_diff)
      endif

      if(.not. ltwostr_only) then
        call VecDestroy(solver%incSolar , ierr) ;call CHKERR(ierr)
        call VecDestroy(solver%b        , ierr) ;call CHKERR(ierr)
        deallocate(solver%incSolar)
        deallocate(solver%b)
      endif
      call destroy_matrices(solver)

      do uid=lbound(solver%solutions,1),ubound(solver%solutions,1)
        call destroy_solution(solver%solutions(uid))
      enddo

      if(allocated(solver%dir_scalevec_Wm2_to_W)) then
        call VecDestroy(solver%dir_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        deallocate(solver%dir_scalevec_Wm2_to_W)
      endif

      if(allocated(solver%diff_scalevec_Wm2_to_W)) then
        call VecDestroy(solver%diff_scalevec_Wm2_to_W, ierr); call CHKERR(ierr)
        deallocate(solver%diff_scalevec_Wm2_to_W)
      endif

      if(allocated(solver%dir_scalevec_W_to_Wm2)) then
        call VecDestroy(solver%dir_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        deallocate(solver%dir_scalevec_W_to_Wm2)
      endif

      if(allocated(solver%diff_scalevec_W_to_Wm2)) then
        call VecDestroy(solver%diff_scalevec_W_to_Wm2, ierr); call CHKERR(ierr)
        deallocate(solver%diff_scalevec_W_to_Wm2)
      endif

      if(allocated(solver%abso_scalevec)) then
        call VecDestroy(solver%abso_scalevec, ierr); call CHKERR(ierr)
        deallocate(solver%abso_scalevec)
      endif

      if(allocated(solver%atm%hgrad)) then
        call VecDestroy(solver%atm%hgrad, ierr); call CHKERR(ierr)
        deallocate(solver%atm%hgrad)
      endif
      if(allocated(solver%atm)) deallocate(solver%atm)

      if(allocated(solver%sun%symmetry_phi)) deallocate(solver%sun%symmetry_phi)
      if(allocated(solver%sun%theta       )) deallocate(solver%sun%theta       )
      if(allocated(solver%sun%phi         )) deallocate(solver%sun%phi         )
      if(allocated(solver%sun%costheta    )) deallocate(solver%sun%costheta    )
      if(allocated(solver%sun%sintheta    )) deallocate(solver%sun%sintheta    )
      if(allocated(solver%sun%xinc        )) deallocate(solver%sun%xinc        )
      if(allocated(solver%sun%yinc        )) deallocate(solver%sun%yinc        )

      if(allocated(solver%OPP)) call solver%OPP%destroy()
      if(allocated(solver%C_dir     )) call DMDestroy(solver%C_dir%da ,ierr); deallocate(solver%C_dir )
      if(allocated(solver%C_diff    )) call DMDestroy(solver%C_diff%da,ierr); deallocate(solver%C_diff)
      if(allocated(solver%C_one     )) call DMDestroy(solver%C_one%da ,ierr); deallocate(solver%C_one )
      if(allocated(solver%C_one1    )) call DMDestroy(solver%C_one1%da,ierr); deallocate(solver%C_one1)
      if(allocated(solver%C_two1    )) call DMDestroy(solver%C_two1%da,ierr); deallocate(solver%C_two1)
      if(allocated(solver%C_one_atm )) call DMDestroy(solver%C_one_atm%da ,ierr); deallocate(solver%C_one_atm)
      if(allocated(solver%C_one_atm1)) call DMDestroy(solver%C_one_atm1%da,ierr); deallocate(solver%C_one_atm1)

      if(allocated(solver%difftop%is_inward)) deallocate(solver%difftop%is_inward)
      if(allocated(solver%diffside%is_inward)) deallocate(solver%diffside%is_inward)
      if(allocated(solver%dirtop%is_inward)) deallocate(solver%dirtop%is_inward)
      if(allocated(solver%dirside%is_inward)) deallocate(solver%dirside%is_inward)

      solver%comm = -1
      solver%linitialized=.False.
      if(solver%myid.eq.0 .and. ldebug)print *,'Destroyed TenStream'

      if(lfinalize) then
        call PetscFinalize(ierr) ;call CHKERR(ierr)
        if(solver%myid.eq.0 .and. ldebug)print *,'Finalized Petsc'
      endif
    endif
  end subroutine

  subroutine destroy_matrices(solver)
    class(t_solver) :: solver

    if(solver%myid.eq.0 .and. ldebug) print *,'Trying to destroy matrices...', allocated(solver%Mdir), allocated(solver%Mdiff)
    if(allocated(solver%Mdir)) then
      call MatDestroy(solver%Mdir , ierr) ;call CHKERR(ierr)
      deallocate(solver%Mdir)
    endif
    if(allocated(solver%Mdiff)) then
      call MatDestroy(solver%Mdiff, ierr) ;call CHKERR(ierr)
      deallocate(solver%Mdiff)
    endif
  end subroutine
end module

