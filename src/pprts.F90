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

  use m_data_parameters, only : ireals, iintegers, irealLUT, &
    & init_mpi_data_parameters, mpiint,                      &
    & zero, one, nan, pi,                                    &
    & nil, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi,      &
    & default_str_len

  use m_helper_functions, only : CHKERR, CHKWARN, deg2rad, rad2deg, imp_allreduce_min, &
    & imp_bcast, imp_allreduce_max, delta_scale, mpi_logical_and, meanval, get_arg, approx, &
    & inc, cstr, toStr, imp_allreduce_mean, imp_min_mean_max, &
    & normalize_vec, vec_proj_on_plane, angle_between_two_normed_vec, cross_3d, &
    & rotation_matrix_world_to_local_basis, deallocate_allocatable, &
    & ind_1d_to_nd

  use m_schwarzschild, only: schwarzschild, B_eff
  use m_optprop, only: t_optprop, &
    & t_optprop_1_2, t_optprop_3_6, t_optprop_3_10, &
    & t_optprop_8_10, t_optprop_3_16, t_optprop_8_16, t_optprop_8_18, &
    & t_optprop_3_10_ann
  use m_eddington, only : eddington_coeff_zdun

  use m_tenstream_options, only : read_commandline_options, ltwostr, luse_eddington, twostr_ratio, &
    options_max_solution_err, options_max_solution_time, ltwostr_only, luse_twostr_guess,        &
    lcalc_nca, lskip_thermal, lschwarzschild, ltopography, &
    lmcrts

  use m_petsc_helpers, only : petscGlobalVecToZero, scatterZerotoPetscGlobal, &
    petscGlobalVecToAll, &
    petscVecToF90, f90VecToPetsc, getVecPointer, restoreVecPointer, hegedus_trick

  use m_mcrts_dmda, only : solve_mcrts

  use m_pprts_base, only : t_solver, t_solver_1_2, t_solver_3_6, t_solver_3_10, &
    t_solver_8_10, t_solver_3_16, t_solver_8_16, t_solver_8_18, &
    t_pprts_shell_ctx, &
    t_coord, t_suninfo, t_atmosphere, compute_gradient, atmk, &
    t_state_container, prepare_solution, destroy_solution, &
    t_dof, t_solver_log_events, setup_log_events, &
    set_dmda_cell_coordinates

  use m_buildings, only: t_pprts_buildings, &
    & PPRTS_TOP_FACE  , &
    & PPRTS_BOT_FACE  , &
    & PPRTS_LEFT_FACE , &
    & PPRTS_RIGHT_FACE, &
    & PPRTS_REAR_FACE , &
    & PPRTS_FRONT_FACE

  use m_pprts_external_solvers, only: twostream, schwarz, pprts_rayli_wrapper, disort

  implicit none
  private

  public :: &
    & init_pprts, &
    & set_optical_properties, set_global_optical_properties, &
    & solve_pprts, set_angles, pprts_get_result, &
    & pprts_get_result_toZero, &
    & gather_all_toZero, &
    & gather_all_to_all, &
    & scale_flx

  logical,parameter :: ldebug=.False.
  logical,parameter :: lcyclic_bc=.True.
  logical,parameter :: lprealloc=.True.

  integer(iintegers),parameter :: minimal_dimension=3 ! this is the minimum number of gridpoints in x or y direction

  interface matcreateshell
    subroutine matcreateshell(comm, mloc, nloc, m, n, ctx, mat, ierr)
      import iintegers, mpiint, t_pprts_shell_ctx, tMat
      integer(mpiint) :: comm, ierr
      integer(iintegers) :: mloc, nloc, m, n
      type(t_pprts_shell_ctx) :: ctx
      type(tMat) :: mat
    end subroutine
  end interface
  interface matshellgetcontext
    subroutine matshellgetcontext(mat, ctx_ptr, ierr)
      import mpiint, t_pprts_shell_ctx, tMat
      type(tMat) :: mat
      type(t_pprts_shell_ctx),pointer :: ctx_ptr
      integer(mpiint) :: ierr
    end subroutine
  end interface
  interface matshellsetcontext
    subroutine matshellsetcontext(mat, ctx, ierr)
      import mpiint, t_pprts_shell_ctx, tMat
      type(tMat) :: mat
      type(t_pprts_shell_ctx) :: ctx
      integer(mpiint) :: ierr
    end subroutine
  end interface
  interface
    subroutine mat_mult_sub(A, x, b, ierr)
      import tMat, tVec, mpiint
      type(tMat), intent(in) :: A
      type(tVec), intent(in) :: x
      type(tVec), intent(inout) :: b
      integer(mpiint), intent(out) :: ierr
    end subroutine
  end interface

  contains

  !> @brief Main routine to setup PPRTS solver
  !> @details This will setup the PETSc DMDA grid and set other grid information, needed for the pprts
  !> \n Nx, Ny Nz are either global domain size or have to be local sizes if present(nxproc,nyproc)
  !> \n where nxproc and nyproc then are the number of pixel per rank for all ranks -- i.e. sum(nxproc) != Nx_global
  subroutine init_pprts(icomm, Nz,Nx,Ny, dx,dy, sundir, solver, dz1d, dz3d, nxproc, nyproc, collapseindex, solvername)
    MPI_Comm, intent(in)          :: icomm         !< @param MPI_Communicator for this solver
    integer(iintegers),intent(in) :: Nz            !< @param[in] Nz     Nz is the number of layers and Nz+1 would be the number of levels
    integer(iintegers),intent(in) :: Nx            !< @param[in] Nx     number of boxes in x-direction
    integer(iintegers),intent(in) :: Ny            !< @param[in] Ny     number of boxes in y-direction
    real(ireals),intent(in)       :: dx            !< @param[in] dx     physical size of grid in [m]
    real(ireals),intent(in)       :: dy            !< @param[in] dy     physical size of grid in [m]
    real(ireals),intent(in)       :: sundir(:)     !< @param[in] cartesian sun direction (pointing away from the sun), dim(3)

    class(t_solver), intent(inout)         :: solver         !< @param[inout] solver
    real(ireals),optional,intent(in)       :: dz1d(:)        !< @param[in]    dz1d    if given, dz1d is used everywhere on the rank
    real(ireals),optional,intent(in)       :: dz3d(:,:,:)    !< @param[in]    dz3d    if given, dz3d has to be local domain size, cannot have global shape
    integer(iintegers),optional,intent(in) :: nxproc(:)      !< @param[in]    nxproc  if given, Nx has to be the local size, dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    integer(iintegers),optional,intent(in) :: nyproc(:)      !< @param[in]    nyproc  if given, Ny has to be the local size, dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    integer(iintegers),optional,intent(in) :: collapseindex  !< @param[in]    collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them
    character(len=*), optional, intent(in) :: solvername     !< @param[in] primarily for logging purposes, name will be prefix to logging stages

    integer(iintegers) :: k,i,j
    logical :: lpetsc_is_initialized

    integer(mpiint) :: ierr

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

      allocate(solver%solutions(-1:1000))
      solver%lenable_solutions_err_estimates = &
        options_max_solution_err.gt.zero .and. options_max_solution_time.gt.zero

      if(.not.approx(dx,dy)) &
        call CHKERR(1_mpiint, 'dx and dy currently have to be the same '//toStr(dx)//' vs '//toStr(dy))

      if(ldebug.and.solver%myid.eq.0) then
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

      call print_domain_geometry_summary(solver, opt_lview=ldebug)

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
    call set_angles(solver, sundir)

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

      if (.not.allocated(solver%atm%hhl)) then
        allocate(solver%atm%hhl)
        call compute_vertical_height_levels(dz=solver%atm%dz, C_hhl=solver%C_one_atm1_box, vhhl=solver%atm%hhl)
      endif

      if(.not.allocated(solver%atm%hgrad)) then
        call compute_gradient(            &
          & comm  = solver%comm,          &
          & atm   = solver%atm,           &
          & C_hhl = solver%C_one_atm1_box,&
          & vhhl  = solver%atm%hhl,       &
          & C_grad= solver%C_two1,        &
          & vgrad = solver%atm%hgrad      )
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
        if(solver%atm%lcollapse) then
          solver%atm%icollapse=collapseindex
          ierr = count(.not.solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs),:,:))
          call CHKWARN(ierr, 'Found non 1D cells in an area that will be collapsed.'// &
            & 'This will change the results! '//&
            & 'collapse index: '//toStr(collapseindex))
          solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs),:,:) = .True. ! if need to be collapsed, they have to be 1D.
          if(ldebug) print *,'Using icollapse:',collapseindex, solver%atm%lcollapse
        endif
      endif

    end subroutine
  end subroutine

  !> @brief Print a summary of the domain shape
  subroutine print_domain_geometry_summary(solver, opt_lview)
    class(t_solver), intent(in) :: solver
    logical, intent(in), optional :: opt_lview

    logical :: lview, lflg
    real(ireals), dimension(3) :: mdz, mhhl
    integer(mpiint) :: myid, numnodes, ierr
    integer(iintegers) :: k
    real(ireals), pointer :: hhl(:,:,:,:)=>null(), hhl1d(:)=>null()

    lview = get_arg(.False., opt_lview)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_view_geometry",&
      lview, lflg, ierr) ;call CHKERR(ierr)
    if(.not.lview) return

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    if(myid.eq.0) print *,''
    if(myid.eq.0) print *,' * Cell Domain Corners'//new_line('')//&
      & cstr(' rank  '      , 'blue')// &
      & cstr(' zs:ze       ', 'green')// &
      & cstr(' xs:xe       ', 'blue' )// &
      & cstr(' ys:ye'       ,'green')

    do k = 0, numnodes-1
      if(myid.eq.k) then
        associate(C => solver%C_one)
          print *,' '//cstr(toStr(myid), 'blue')//&
            & '     '//cstr(toStr(C%zs), 'green')//' : '//cstr(toStr(C%ze), 'green')// &
            & '     '//cstr(toStr(C%xs), 'blue' )//' : '//cstr(toStr(C%xe), 'blue' )// &
            & '     '//cstr(toStr(C%ys), 'green')//' : '//cstr(toStr(C%ye), 'green')
        end associate
      endif
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    enddo
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    if(myid.eq.0) print *,''

    if(myid.eq.0) print *,'*  k(layer)'// &
      & ' '//cstr('all(1d) / any(1d)', 'blue' )//' '// &
      & ' '//cstr('dz(min/mean/max)', 'red')//' [m]                   '// &
      & ' '//cstr('height [km]', 'blue')
    associate( &
        & atm => solver%atm, &
        & C_one_atm => solver%C_one_atm, &
        & C_one_atm1 => solver%C_one_atm1 )

      call getVecPointer(C_one_atm1%da, atm%hhl, hhl1d, hhl, readonly=.True.)
      do k=C_one_atm%zs,C_one_atm%ze

        call imp_min_mean_max(solver%comm, atm%dz(k,:,:), mdz)
        call imp_min_mean_max(solver%comm, (hhl(i0,k,:,:)+hhl(i0,k+1,:,:))/2, mhhl)

        if(myid.eq.0) &
          & print *, k, &
          & cstr(toStr(all(atm%l1d(k,:,:))),'blue'), '        ', &
          & cstr(toStr(any(atm%l1d(k,:,:))),'blue'), '        ', &
          & cstr(toStr(mdz),'red'), ' ', cstr(toStr(mhhl*1e-3_ireals),'blue')
      enddo
      call restoreVecPointer(C_one_atm1%da, atm%hhl, hhl1d, hhl, readonly=.True.)

      if(myid.eq.0) print *,' * atm dx/dy '//toStr(atm%dx)//' , '//toStr(atm%dy)
    end associate
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
  end subroutine


  !> @brief Construct PETSC grid information for regular DMDA
  !> @details setup DMDA grid containers for direct, diffuse and absorption grid
  !>  \n and fill user context containers(t_coord) which include local as well as global array sizes
  !>  \n every mpi rank has to call this
  subroutine setup_grid(solver, Nz_in, Nx, Ny, nxproc, nyproc, collapseindex)
      class(t_solver), intent(inout) :: solver
      integer(iintegers), intent(in) :: Nz_in, Nx, Ny                  !< @param[in] local number of grid boxes -- in the vertical we have Nz boxes and Nz+1 levels
      integer(iintegers), intent(in), optional :: nxproc(:), nyproc(:) !< @param[in] size of local domains on each node
      integer(iintegers), optional, intent(in) :: collapseindex        !< @param[in] collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them

      integer(iintegers) :: Nz
      DMBoundaryType, parameter :: &
        & bp=DM_BOUNDARY_PERIODIC, &
        & bn=DM_BOUNDARY_NONE,     &
        & bm=DM_BOUNDARY_MIRROR!,   &
        !& bg=DM_BOUNDARY_GHOSTED

      DMBoundaryType :: boundaries(3)
      integer(iintegers), allocatable :: nxprocp1(:), nyprocp1(:) !< last entry one larger for vertices
      integer(mpiint) :: ierr

      if(lcyclic_bc) then
        boundaries = [bn, bp, bp]
      else
        boundaries = [bn, bm, bm]
      endif

      Nz = Nz_in
      if(present(collapseindex)) then
        if(collapseindex.gt.1) then
          Nz = Nz_in-collapseindex+i1
        endif
      endif

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,&
        & 'Setting up the DMDA grid for',Nz,Nx,Ny,'using',solver%numnodes,'nodes'

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA C_diff', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%C_diff, Nz+1,Nx,Ny, boundaries, &
        & solver%difftop%dof + 2* solver%diffside%dof, nxproc,nyproc)

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA C_dir', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%C_dir, Nz+1,Nx,Ny, boundaries, &
        & solver%dirtop%dof + 2* solver%dirside%dof, nxproc,nyproc)

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA C_one', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%C_one , Nz  , Nx,Ny, boundaries, &
        & i1, nxproc,nyproc)
      call setup_dmda(solver%comm, solver%C_one1, Nz+1, Nx,Ny, boundaries, &
        & i1, nxproc,nyproc)

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA C_two', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%C_two1, Nz+1, Nx,Ny, boundaries, &
        & i2, nxproc,nyproc)

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA atm', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%C_one_atm , Nz_in  , Nx,Ny, boundaries, &
        & i1, nxproc,nyproc)

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA atm1', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%C_one_atm1, Nz_in+1, Nx,Ny, boundaries, &
        & i1, nxproc,nyproc)

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA atm1_box', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%C_one_atm1_box, Nz_in+1, Nx,Ny, boundaries, &
        & i1, nxproc,nyproc, DMDA_STENCIL_BOX)

      if(ldebug) then
        if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA srfc_one', 'green')
        call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      endif
      call setup_dmda(solver%comm, solver%Csrfc_one, i1, Nx,Ny, boundaries, &
        & i1, nxproc,nyproc)

      if(present(nxproc)) then
        allocate(nxprocp1(size(nxproc)), nyprocp1(size(nyproc)))
        nxprocp1(:) = nxproc
        nyprocp1(:) = nyproc
        nxprocp1(size(nxprocp1)) = nxprocp1(size(nxprocp1)) +1
        nyprocp1(size(nyprocp1)) = nyprocp1(size(nyprocp1)) +1

        if(ldebug) then
          if(solver%myid.eq.0) print *,solver%myid, cstr('Configuring DMDA vert_one_atm1', 'green')
          call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
        endif
        call setup_dmda(solver%comm, solver%Cvert_one_atm1, Nz_in+1, Nx+1,Ny+1, boundaries, &
          & i1, nxprocp1,nyprocp1)
      else
        call setup_dmda(solver%comm, solver%Cvert_one_atm1, Nz_in+1, Nx+1,Ny+1, boundaries, &
          & i1)
      endif

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'DMDA grid ready'
    end subroutine

    subroutine setup_dmda(icomm, C, Nz, Nx, Ny, boundary, dof, nxproc, nyproc, stencil_type)
      integer(mpiint), intent(in) :: icomm
      type(t_coord), allocatable :: C
      integer(iintegers), intent(in) :: Nz,Nx,Ny,dof
      DMBoundaryType, intent(in) :: boundary(:)
      integer(iintegers), intent(in), optional :: nxproc(:), nyproc(:) ! size of local domains on each node
      DMDAStencilType, intent(in), optional :: stencil_type
      DMDAStencilType :: opt_stencil_type

      integer(iintegers), parameter :: stencil_size=1
      integer(mpiint) :: myid, ierr
      call mpi_comm_rank(icomm, myid, ierr); call CHKERR(ierr)

      if(any([present(nxproc), present(nyproc)])&
        & .and..not.all([present(nxproc), present(nyproc)])) &
        & call CHKERR(1_mpiint, 'have to have both, nxproc AND nyproc present or none')

      if(ldebug.and.present(nxproc).and.present(nyproc)) &
        & print *, myid, 'setup_dmda nxproc', nxproc, 'nyproc', nyproc

      opt_stencil_type = get_arg(DMDA_STENCIL_STAR, stencil_type)

      allocate(C)

      C%comm = icomm

      C%dof = i1*dof
      if(present(nxproc) .and. present(nyproc) ) then
        call DMDACreate3d( C%comm ,                                                     &
          boundary(1)             , boundary(2)                , boundary(3)          , &
          opt_stencil_type        ,                                                     &
          Nz                      , sum(nxproc)                , sum(nyproc)          , &
          i1                      , size(nxproc,kind=iintegers), size(nyproc,kind=iintegers), &
          C%dof                   , stencil_size               ,                        &
          [Nz]                    , nxproc                     , nyproc               , &
          C%da                    , ierr)
        call CHKERR(ierr)
      else
        call DMDACreate3d( C%comm ,                                                   &
          boundary(1)             , boundary(2)              , boundary(3)          , &
          opt_stencil_type        ,                                                   &
          i1*Nz                   , Nx                       , Ny                   , &
          i1                      , PETSC_DECIDE             , PETSC_DECIDE         , &
          C%dof                   , stencil_size             ,                        &
          [Nz]                    , PETSC_NULL_INTEGER       , PETSC_NULL_INTEGER   , &
          C%da                    , ierr) ;call CHKERR(ierr)
      endif

      ! need this first setfromoptions call because of a bug which happens with intel compilers
      call DMSetFromOptions(C%da, ierr)                                    ; call CHKERR(ierr)

      call DMSetMatType(C%da, MATAIJ, ierr)                                ; call CHKERR(ierr)
      if(lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ; call CHKERR(ierr)
      call DMSetFromOptions(C%da, ierr)                                    ; call CHKERR(ierr)
      call DMSetup(C%da,ierr)                                              ; call CHKERR(ierr)

      if(ldebug) call DMView(C%da, PETSC_VIEWER_STDOUT_WORLD ,ierr)        ; call CHKERR(ierr)
      call setup_coords(C)
    contains
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

        allocate(C%neighbors(0:3**C%dim-1) )
        call DMDAGetNeighbors(C%da,C%neighbors,ierr) ;call CHKERR(ierr)
        call mpi_comm_size(icomm, numnodes, ierr); call CHKERR(ierr)
        if(numnodes.gt.i1) then
          if(ldebug.and.(C%dim.eq.3)) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors(10),  &
            C%neighbors(4) ,  &
            C%neighbors(16),  &
            C%neighbors(22),  &
            'while I am ',C%neighbors(13)

          if(ldebug.and.(C%dim.eq.2)) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors(1), &
            C%neighbors(3), &
            C%neighbors(7), &
            C%neighbors(5), &
            'while I am ',C%neighbors(4)
        endif
        if(C%glob_xm.lt.i3) call CHKERR(1_mpiint, 'Global domain is too small in x-direction (Nx='//toStr(C%glob_xm)// &
          '). However, need at least 3 because of horizontal ghost cells')
        if(C%glob_ym.lt.i3) call CHKERR(1_mpiint, 'Global domain is too small in y-direction (Ny='//toStr(C%glob_ym)// &
          '). However, need at least 3 because of horizontal ghost cells')

        if(C%xm.lt.i1) call CHKERR(1_mpiint, 'Local domain is too small in x-direction (Nx='//toStr(C%xm)// &
          '). However, need at least 1')
        if(C%ym.lt.i1) call CHKERR(1_mpiint, 'Local domain is too small in y-direction (Ny='//toStr(C%ym)// &
          '). However, need at least 1')
      end subroutine
    end subroutine


    !> @brief Determine height levels by summing up the atm%dz with the assumption that TOA is at a constant value
    !>        or a max_height is given in the option database
    subroutine compute_vertical_height_levels(dz, C_hhl, vhhl)
      type(t_coord), intent(in) :: C_hhl
      real(ireals), intent(in) :: dz(:,:,:)
      type(tVec), intent(inout) :: vhhl

      type(tVec) :: g_hhl
      real(ireals), pointer :: hhl(:,:,:,:)=>null(), hhl1d(:)=>null()
      real(ireals) :: max_height, global_max_height

      integer(mpiint) :: comm, ierr
      integer(iintegers) :: k,i,j
      integer(iintegers) :: kk,ii,jj
      logical :: lflg

      max_height = zero

      call DMGetGlobalVector(C_hhl%da, g_hhl, ierr) ;call CHKERR(ierr)
      call VecSet(g_hhl, zero, ierr); call CHKERR(ierr)

      call getVecPointer(C_hhl%da, g_hhl, hhl1d, hhl)

      do j = C_hhl%ys, C_hhl%ye
        jj = i1 + j - C_hhl%ys
        do i = C_hhl%xs, C_hhl%xe
          ii = i1 + i - C_hhl%xs

          hhl(i0, C_hhl%zs, i, j ) = zero
          do k = C_hhl%zs, C_hhl%ze-1
            kk = i1 + k - C_hhl%zs
            hhl(i0,k+1,i,j) = hhl(i0,k,i,j) - dz(kk,ii,jj)
            max_height = min(max_height, hhl(i0,k+1,i,j))
          enddo
        enddo
      enddo

      call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_global_max_height", &
        global_max_height, lflg , ierr) ;call CHKERR(ierr)
      if(.not.lflg) then
        call PetscObjectGetComm(C_hhl%da, comm, ierr); call CHKERR(ierr)
        call imp_allreduce_min(comm, max_height, global_max_height)
      endif
      hhl(i0, :, :, :) = hhl(i0, :, :, :) - global_max_height

      call restoreVecPointer(C_hhl%da, g_hhl, hhl1d, hhl)

      call DMCreateLocalVector(C_hhl%da, vhhl, ierr) ;call CHKERR(ierr)
      call DMGlobalToLocalBegin(C_hhl%da, g_hhl, INSERT_VALUES, vhhl,ierr) ;call CHKERR(ierr)
      call DMGlobalToLocalEnd(C_hhl%da, g_hhl, INSERT_VALUES, vhhl,ierr) ;call CHKERR(ierr)
      call DMRestoreGlobalVector(C_hhl%da, g_hhl, ierr) ;call CHKERR(ierr)

      call PetscObjectViewFromOptions(vhhl, PETSC_NULL_VEC, "-pprts_show_hhl", ierr); call CHKERR(ierr)
    end subroutine

    !> @brief initialize basic memory structs like PETSc vectors and matrices
    subroutine init_memory(C_dir, C_diff, incSolar,b)
      type(t_coord), intent(in) :: C_dir, C_diff
      type(tVec), intent(inout), allocatable :: b,incSolar
      integer(mpiint) :: ierr

      if(ltwostr_only) return

      if(.not.allocated(incSolar)) allocate(incSolar)
      if(.not.allocated(b)) allocate(b)

      call DMCreateGlobalVector(C_dir%da,incSolar,ierr) ; call CHKERR(ierr)
      call DMCreateGlobalVector(C_diff%da,b,ierr)       ; call CHKERR(ierr)

      call VecSet(incSolar,zero,ierr) ; call CHKERR(ierr)
      call VecSet(b,zero,ierr)        ; call CHKERR(ierr)
    end subroutine

    subroutine set_angles(solver, sundir)
      class(t_solver), intent(inout)   :: solver
      real(ireals),intent(in)          :: sundir(:)      !< @param[in] cartesian sun direction
      logical :: luse_ann, lflg
      integer(mpiint) :: ierr

      luse_ann = .False.
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_use_ANN", &
        luse_ann, lflg , ierr) ;call CHKERR(ierr)

      if(.not.solver%linitialized) then
        print *,solver%myid,'You tried to set angles in the PPRTS solver.  &
          & This should be called right after init_pprts'
        ierr=1; call CHKERR(ierr)
      endif

      call setup_suninfo(solver, sundir, solver%sun, solver%C_one)

      if(ltwostr_only .or. lmcrts) return ! dont need anything here, we just compute Twostream anyway

      ! init box montecarlo model
      select type(solver)
        class is (t_solver_1_2)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_1_2::solver%OPP)

        class is (t_solver_3_6)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_3_6::solver%OPP)

        class is (t_solver_3_10)
           if(.not.allocated(solver%OPP) ) then
             if(luse_ann) then
               allocate(t_optprop_3_10_ann::solver%OPP)
             else
               allocate(t_optprop_3_10::solver%OPP)
             endif
           endif

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

      ! reset Matrices to generate new preallocation
  end subroutine

  !> @brief set direction where sun stands
  !> @details save sun azimuth and zenith angle
  !>   \n sun azimuth is reduced to the range of [0,90] and the transmission of direct radiation is contributed for by a integer increment,
  !>   \n determining which neighbouring box is used in the horizontal direction
  subroutine setup_suninfo(solver, sundir, sun, C_one)
    class(t_solver), intent(in) :: solver
    real(ireals),intent(in) :: sundir(:)
    type(t_suninfo),intent(inout) :: sun
    type(t_coord), intent(in) :: C_one

    real(ireals),pointer :: grad   (:,:,:,:) =>null()
    real(ireals),pointer :: grad_1d(:)       =>null()
    real(ireals) :: loc_sundir(3), proj_sundir(3), e_x(3), e_y(3), e_z(3), Mrot(3,3)
    real(ireals) :: az, zenith
    integer(iintegers) :: k,i,j
    integer(mpiint) :: ierr

    if(.not.allocated(solver%atm%hgrad)) call CHKERR(1_mpiint, 'atm%hgrad not initialized!')

    call alloc_sun_rfield(sun%symmetry_phi)
    call alloc_sun_rfield(sun%theta)
    call alloc_sun_rfield(sun%phi)
    call alloc_sun_rfield(sun%costheta)
    call alloc_sun_rfield(sun%sintheta)
    call alloc_sun_ifield(sun%xinc)
    call alloc_sun_ifield(sun%yinc)

    sun%sundir = sundir
    call normalize_vec(sun%sundir, ierr)

    ! default unit vectors (if not ltopography)
    e_x = [one, zero, zero]
    e_y = [zero, one, zero]
    e_z = [zero, zero, one]
    loc_sundir = sun%sundir

    call getVecPointer(solver%C_two1%da, solver%atm%hgrad, grad_1d, grad)
    do j=C_one%ys,C_one%ye
      do i=C_one%xs,C_one%xe
        do k=C_one%zs,C_one%ze

          if(sun%luse_topography) then
            e_x = [one, zero, grad(i0, k, i, j)]
            e_y = [zero, one, grad(i1, k, i, j)]
            call normalize_vec(e_x, ierr)
            call normalize_vec(e_y, ierr)
            e_z = cross_3d(e_x, e_y) ! upward pointing normal on top face
            Mrot = rotation_matrix_world_to_local_basis(e_x, e_y, e_z)
            loc_sundir = matmul(Mrot, sun%sundir)
          endif

          proj_sundir = vec_proj_on_plane(loc_sundir, -e_z)
          call normalize_vec(proj_sundir, ierr)

          az = angle_between_two_normed_vec(proj_sundir, -e_y)
          az = az * sign(one, dot_product(proj_sundir, -e_x))
          sun%phi(k,i,j) = rad2deg(az)

          zenith = angle_between_two_normed_vec(loc_sundir, -e_z)
          sun%theta(k,i,j) = rad2deg(zenith)
          if(sun%theta(k,i,j).ge.90._ireals) sun%theta(k,i,j) = -one
        enddo
      enddo
    enddo
    call restoreVecPointer(solver%C_two1%da, solver%atm%hgrad, grad_1d, grad)

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

    call print_suninfo_summary(solver, ldebug)

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

  !> @brief Print a summary of the sun_info
  subroutine print_suninfo_summary(solver, opt_lview)
    class(t_solver), intent(in) :: solver
    logical, intent(in), optional :: opt_lview

    logical :: lview, lflg
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: k
    real(ireals), dimension(3) :: mtheta, mphi, msundir

    lview = get_arg(.False., opt_lview)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_view_suninfo",&
      lview, lflg, ierr) ;call CHKERR(ierr)
    if(.not.lview) return

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
    associate(sun => solver%sun, C_one => solver%C_one)

      do k = 1, 3
        call imp_min_mean_max(solver%comm, sun%sundir(k), msundir)
        if(myid.eq.0) then
          if(.not. all(approx(sun%sundir(k), msundir))) then
            print *, 'Sundir rank=0 :: ', sun%sundir
            ierr = int(k, mpiint)
            call CHKERR(ierr, 'sundir component '//toStr(k)//' does not match across ranks '//&
              & 'min mean max sundir('//toStr(k)//') = '//toStr(msundir))
          endif
        endif
      enddo

      if(solver%myid.eq.0) &
        & print *, ' * '//cstr('sundir ['//toStr(sun%sundir)//'] ', 'red'), &
        & ' count(xinc) ', count(sun%xinc.eq.0),count(sun%xinc.eq.i1), &
        & ' count(yinc) ', count(sun%yinc.eq.0),count(sun%yinc.eq.i1)

      if(myid.eq.0) print *,' * k(layer)  '// &
        & ' '//cstr('theta(min/mean/max)', 'red')//'                    '// &
        & ' '//cstr('phi', 'blue')
      do k = C_one%zs, C_one%ze
        call imp_min_mean_max(solver%comm, sun%phi(k,:,:), mphi)
        call imp_min_mean_max(solver%comm, sun%theta(k,:,:), mtheta)
        if(myid.eq.0) &
          & print *, k, &
          & cstr(toStr(mtheta),'red'), ' ', cstr(toStr(mphi),'blue')
      enddo

    end associate
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
  end subroutine

  !> @brief setup topography information
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information and
  !>   \n tweak the local sun angles to bend the rays.

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
    integer(mpiint) :: ierr

    if(allocated(A)) return

    allocate(A)
    call DMCreateMatrix(C%da, A, ierr) ;call CHKERR(ierr)

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
    call MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr) ;call CHKERR(ierr)

    ! pressure mesh  may wiggle a bit and change atm%l1d -- keep the nonzeros flexible
    call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr) ;call CHKERR(ierr)

    ! call MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr) ;call CHKERR(ierr) ! dont throw away the zero -- this completely destroys preallocation performance

    call MatSetUp(A,ierr) ;call CHKERR(ierr)
  end subroutine

  subroutine mat_set_diagonal(A,vdiag)
    type(tMat) :: A
    real(ireals),intent(in),optional :: vdiag

    integer(iintegers) :: is, ie, irow
    real(ireals) :: v
    integer(mpiint) :: ierr

    v = get_arg(one, vdiag)

    call MatGetOwnershipRange(A, is, ie, ierr); call CHKERR(ierr)
    do irow = is, ie-1
      call MatSetValue(A, irow, irow, v, INSERT_VALUES, ierr); call CHKERR(ierr)
    enddo
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

    call getVecPointer(C%da, v_o_nnz, xo1d, xo)
    call getVecPointer(C%da, v_d_nnz, xd1d, xd)

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

    call restoreVecPointer(C%da, v_o_nnz, xo1d, xo)
    call restoreVecPointer(C%da, v_d_nnz, xd1d, xd)

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

    call getVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.True.)
    call getVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.True.)

    call VecGetLocalSize(g_d_nnz,vsize,ierr) ;call CHKERR(ierr)
    allocate(o_nnz(0:vsize-1))
    allocate(d_nnz(0:vsize-1))

    o_nnz=int(xo1d, kind=iintegers)
    d_nnz=int(xd1d, kind=iintegers) + i1  ! +1 for diagonal entries

    call restoreVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.True.)
    call restoreVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.True.)

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

    call getVecPointer(C%da, v_o_nnz, xo1d, xo)
    call getVecPointer(C%da, v_d_nnz, xd1d, xd)


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

    call restoreVecPointer(C%da, v_o_nnz, xo1d, xo)
    call restoreVecPointer(C%da, v_d_nnz, xd1d, xd)

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

    call getVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.True.)
    call getVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.True.)

    call VecGetLocalSize(g_d_nnz, vsize, ierr) ;call CHKERR(ierr)
    allocate(o_nnz(0:vsize-1))
    allocate(d_nnz(0:vsize-1))

    o_nnz=int(xo1d, kind=iintegers)
    d_nnz=int(xd1d, kind=iintegers) + i1  ! +1 for diagonal entries

    call restoreVecPointer(C%da, g_o_nnz, xo1d, xo, readonly=.True.)
    call restoreVecPointer(C%da, g_d_nnz, xd1d, xd, readonly=.True.)

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

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call MatGetInfo(A,MAT_LOCAL,info,ierr); call CHKWARN(ierr)
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
    logical :: lpprts_no_absorption, lpprts_no_scatter
    integer(mpiint) :: ierr

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
                              ' '//toStr(solver%myid)//&
                              ' kabs min '//toStr(minval(kabs))//' max '//toStr(maxval(kabs))//&
                              ' ksca min '//toStr(minval(ksca))//' max '//toStr(maxval(ksca))//&
                              ' g    min '//toStr(minval(g   ))//' max '//toStr(maxval(g   )))
      endif
      if( any(isnan([kabs,ksca,g]))) then
        call CHKERR(1_mpiint, 'set_optical_properties :: found NaN value in optical properties!'//&
                              ' NaN in kabs? '//toStr(any(isnan(kabs)))// &
                              ' NaN in ksca? '//toStr(any(isnan(ksca)))// &
                              ' NaN in g   ? '//toStr(any(isnan(g   ))))
      endif
    endif
    if(ldebug) then
      if( (any([atm%kabs,atm%ksca,atm%g].lt.zero)) .or. (any(isnan([atm%kabs,atm%ksca,atm%g]))) ) then
        call CHKERR(1_mpiint, 'set_optical_properties :: found illegal value in optical properties! '//&
                              ' '//toStr(solver%myid)//&
                              ' kabs min '//toStr(minval(atm%kabs))//' max '//toStr(maxval(atm%kabs))//&
                              ' ksca min '//toStr(minval(atm%ksca))//' max '//toStr(maxval(atm%ksca))//&
                              ' g    min '//toStr(minval(atm%g   ))//' max '//toStr(maxval(atm%g   )))
      endif
    endif

    lpprts_no_absorption = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_no_absorption", &
      lpprts_no_absorption, lflg , ierr) ;call CHKERR(ierr)

    if(lpprts_no_absorption) then
      atm%kabs = 0
    endif

    lpprts_no_scatter = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_no_scatter", &
      lpprts_no_scatter, lflg , ierr) ;call CHKERR(ierr)

    if(lpprts_no_scatter) then
      atm%ksca = 0
    endif


    lpprts_delta_scale = get_arg(.True., ldelta_scaling)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_delta_scale", &
      lpprts_delta_scale, lflg , ierr) ;call CHKERR(ierr)

    pprts_delta_scale_max_g=.85_ireals-epsilon(pprts_delta_scale_max_g)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_delta_scale_max_g", &
      pprts_delta_scale_max_g, lflg , ierr) ;call CHKERR(ierr)

    if(lpprts_delta_scale) then
      call delta_scale(atm%kabs, atm%ksca, atm%g, max_g=pprts_delta_scale_max_g)
    else
      if(solver%myid.eq.0.and.lflg) print *,"Skipping Delta scaling of optprops"
      if(any(atm%g.ge.0.85_ireals)) &
        call CHKWARN(1_mpiint, 'Skipping delta scaling but now we have values of '// &
        'g > '//toStr(pprts_delta_scale_max_g)// &
        ' (max='//toStr(maxval(atm%g))//')')
    endif

    call print_optical_properties_summary(solver, opt_lview=ldebug)

    if(ltwostr_only) then
      return ! twostream should not depend on eddington coeffs... it will have to calculate it on its own.
    endif

    if( (any([atm%kabs,atm%ksca,atm%g].lt.zero)) .or. (any(isnan([atm%kabs,atm%ksca,atm%g]))) ) then
      call CHKERR(1_mpiint, 'set_optical_properties :: found illegal value in delta_scaled optical properties! abort!')
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
                sun%costheta(i0,i,j), &
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

    subroutine print_optical_properties_summary(solver, opt_lview)
      class(t_solver), intent(in) :: solver
      logical, intent(in), optional :: opt_lview

      logical :: lview, lflg
      real(ireals), dimension(3) :: mkabs, mksca, mg, malbedo, mplck
      integer(mpiint) :: myid, ierr
      integer(iintegers) :: k

      lview = get_arg(.False., opt_lview)
      call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-pprts_view_optprop",&
        lview, lflg, ierr) ;call CHKERR(ierr)
      if(.not.lview) return


      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(solver%comm, myid, ierr); call CHKERR(ierr)
      if(myid.eq.0) print *,'*      k(layer)  '// &
        '     '//cstr('kabs(min/mean/max)', 'blue')//'                  '// &
        '     '//cstr('ksca              ', 'red' )//'                  '// &
        '     '//cstr('g                 ', 'blue')//'                  '// &
        '     '//cstr('plck              ', 'red')
      associate( atm => solver%atm, C_one_atm => solver%C_one_atm )

        do k=C_one_atm%zs,C_one_atm%ze

          if(allocated(atm%kabs)) call imp_min_mean_max(solver%comm, atm%kabs(k,:,:), mkabs)
          if(allocated(atm%ksca)) call imp_min_mean_max(solver%comm, atm%ksca(k,:,:), mksca)
          if(allocated(atm%g   )) call imp_min_mean_max(solver%comm, atm%g   (k,:,:), mg   )
          if(allocated(atm%planck)) then
            call imp_min_mean_max(solver%comm, atm%planck(k,:,:), mplck)
          else
            mplck = nan
          endif

          if(myid.eq.0) &
            & print *,k,cstr(toStr(mkabs),'blue'), cstr(toStr(mksca),'red'), cstr(toStr(mg),'blue'), cstr(toStr(mplck),'red')

        enddo

        if(allocated(atm%albedo)) then
          call imp_min_mean_max(solver%comm, atm%albedo(:,:), malbedo)
          if(myid.eq.0) print *, ' * Albedo (min/mean,max)', malbedo
        endif


        if(myid.eq.0) then
          print *,' * Number of 1D layers: ', count(atm%l1d) , size(atm%l1d),'(',(100._ireals* count(atm%l1d) )/size(atm%l1d),'%)'
          print *,' * shape local optprop:', shape(atm%kabs)
        endif

      end associate
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
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
    integer(mpiint) :: ierr

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

  subroutine solve_pprts(solver, lthermal, lsolar, edirTOA, opt_solution_uid, opt_solution_time, opt_buildings)
    class(t_solver), intent(inout)          :: solver
    logical, intent(in)                     :: lthermal, lsolar
    real(ireals),intent(in)                 :: edirTOA
    integer(iintegers),optional, intent(in) :: opt_solution_uid
    real(ireals),      optional, intent(in) :: opt_solution_time
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    integer(iintegers) :: uid
    logical            :: lflg, derived_lsolar, luse_rayli, lrayli_snapshot, luse_disort
    integer(mpiint) :: ierr

    if(.not.allocated(solver%atm)) call CHKERR(1_mpiint, 'atmosphere is not allocated?!')
    if(.not.allocated(solver%atm%kabs)) &
      call CHKERR(1_mpiint, 'atmosphere%kabs is not allocated! - maybe you need to call set_optical_properties() first')
    if(.not.allocated(solver%atm%ksca)) &
      call CHKERR(1_mpiint, 'atmosphere%ksca is not allocated! - maybe you need to call set_optical_properties() first')

    uid = get_arg(0_iintegers, opt_solution_uid)

    associate( solution => solver%solutions(uid) )

    if(lthermal .and. (.not.allocated(solver%atm%planck))) &
      & call CHKERR(1_mpiint, 'you asked to compute thermal radiation but '// &
      & 'did not provide planck emissions in set_optical_properties')

    derived_lsolar = mpi_logical_and(solver%comm, &
      & lsolar.and. &
      & edirTOA.gt.zero.and. &
      & any(solver%sun%theta.ge.zero))

    if(ldebug) print *,'uid', uid, 'lsolar', derived_lsolar, &
      & 'edirTOA', edirTOA, 'any_theta in [0..90]?', any(solver%sun%theta.ge.zero)

    if(derived_lsolar.and.lthermal) then
      print *,'uid', uid, 'lthermal', lthermal, 'lsolar', lsolar, derived_lsolar, &
        & 'edirTOA', edirTOA, 'any_theta in [0..90]?', any(solver%sun%theta.ge.zero)
      call CHKERR(1_mpiint, 'Somehow ended up with a request to compute solar and thermal radiation in one call.'//new_line('')// &
        & '       This is currently probably not going to work or at least has to be tested further...'//new_line('')// &
        & '       I recommend you call it one after another')
    endif

    if(.not.solution%lset) then
      call prepare_solution(solver%C_dir%da, solver%C_diff%da, solver%C_one%da, &
        lsolar=derived_lsolar, lthermal=lthermal, solution=solution, uid=uid)
    else
      if(solution%lsolar_rad.neqv.derived_lsolar) then
        call destroy_solution(solution)
        call prepare_solution(solver%C_dir%da, solver%C_diff%da, solver%C_one%da, &
          lsolar=derived_lsolar, lthermal=lthermal, solution=solution, uid=uid)
      endif
    endif

    if((solution%lthermal_rad.eqv..False.).and. (solution%lsolar_rad.eqv..False.)) then ! nothing else to do
      if(allocated(solution%edir)) then
        call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
      endif
      if(allocated(solution%ediff)) then
        call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
      endif
      solution%lchanged=.True.
      goto 99 ! quick exit
    endif

    ! --------- Skip Thermal Computation (-lskip_thermal) --
    if(lskip_thermal .and. (solution%lsolar_rad.eqv..False.)) then ! nothing else to do
      if(ldebug .and. solver%myid.eq.0) print *,'skipping thermal calculation -- returning zero flux'
      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)
      solution%lchanged=.True.
      goto 99 ! quick exit
    endif

    ! --------- Calculate Radiative Transfer with RayLi ------------
    call PetscLogEventBegin(solver%logs%solve_mcrts, ierr)
    luse_rayli = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-pprts_use_rayli", luse_rayli, lflg, ierr) ; call CHKERR(ierr)
    lrayli_snapshot = .False.
    call PetscOptionsHasName(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-rayli_snapshot", lrayli_snapshot, ierr) ; call CHKERR(ierr)
    call pprts_rayli_wrapper(luse_rayli, lrayli_snapshot, solver, edirTOA, solution, opt_buildings)
    call PetscLogEventEnd(solver%logs%solve_mcrts, ierr)
    if(luse_rayli) goto 99

    ! --------- Calculate Radiative Transfer with Disort ------------
    luse_disort = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-pprts_use_disort", luse_disort, lflg, ierr) ; call CHKERR(ierr)
    if(luse_disort) then
      call PetscLogEventBegin(solver%logs%solve_disort, ierr)
      call disort(solver, edirTOA, solution)
      call PetscLogEventEnd(solver%logs%solve_disort, ierr)
      goto 99
    endif

    ! --------- Calculate 1D Radiative Transfer ------------
    if(  ltwostr &
      .or. (solution%lthermal_rad .and. lcalc_nca) &
      .or. (solution%lthermal_rad .and. lschwarzschild) ) then


      if( solution%lthermal_rad .and. lschwarzschild ) then
        call PetscLogEventBegin(solver%logs%solve_schwarzschild, ierr)
        call schwarz(solver, solution)
        call PetscLogEventEnd(solver%logs%solve_schwarzschild, ierr)
      else
        call PetscLogEventBegin(solver%logs%solve_twostream, ierr)
        call twostream(solver, edirTOA, solution, opt_buildings)
        call PetscLogEventEnd(solver%logs%solve_twostream, ierr)
      endif

      if(ldebug .and. solver%myid.eq.0) print *,'1D calculation done'


      if( ltwostr_only ) goto 99
      if( solution%lthermal_rad .and. lcalc_nca ) goto 99
      if( solution%lthermal_rad .and. lschwarzschild ) goto 99
    endif

    if( lmcrts ) then
      call PetscLogEventBegin(solver%logs%solve_mcrts, ierr)
      call solve_mcrts(solver, edirTOA, solution)
      call PetscLogEventEnd(solver%logs%solve_mcrts, ierr)
      goto 99
    endif

    call pprts(solver, edirTOA, solution, opt_buildings)

    99 continue ! this is the quick exit final call where we clean up before the end of the routine

    call restore_solution(solver, solution, opt_solution_time)

    end associate
  end subroutine

  !> @brief call the matrix assembly and petsc solve routines for pprts solvers
  subroutine setup_matshell(solver, C, A, mat_mult_subroutine)
    class(t_solver), target, intent(inout) :: solver
    type(t_coord), intent(in) :: C
    type(tMat), allocatable, intent(inout) :: A
    procedure(mat_mult_sub) :: mat_mult_subroutine

    integer(iintegers) :: Nlocal, Nglobal
    integer(mpiint) :: ierr

    solver%shell_ctx%solver => solver

    if(.not.allocated(A)) then
      allocate(A)
      Nlocal = C%zm * C%xm * C%ym * C%dof
      Nglobal = C%glob_zm * C%glob_xm * C%glob_ym * C%dof
      call MatCreateShell(solver%comm, Nlocal, Nlocal, Nglobal, Nglobal, solver%shell_ctx, A, ierr); call CHKERR(ierr)
      call MatShellSetContext(A, solver%shell_ctx, ierr); call CHKERR(ierr)
      call MatShellSetOperation(A, MATOP_MULT, mat_mult_subroutine, ierr); call CHKERR(ierr)
    endif
  end subroutine

  !> @brief define the matmult function for dir2dir computations (used for MatShell)
  subroutine op_mat_mult_edir(A, x, b, ierr)
    type(tMat), intent(in) :: A
    type(tVec), intent(in) :: x
    type(tVec), intent(inout) :: b
    integer(mpiint), intent(out) :: ierr

    type(t_pprts_shell_ctx), pointer :: ctx_ptr
    class(t_solver), pointer :: solver
    real(ireals),pointer,dimension(:,:,:,:) :: xx=>null(), xb=>null()
    real(ireals),pointer,dimension(:) :: xx1d=>null(), xb1d=>null()
    real(irealLUT), target, allocatable :: coeff(:)
    real(irealLUT), pointer :: v(:,:) ! dim(src, dst)
    integer(iintegers) :: k,i,j
    integer(iintegers) :: idst, isrc, src, dst, xinc, yinc
    type(tVec) :: lb, lx

    nullify(ctx_ptr)
    call matshellgetcontext(A, ctx_ptr, ierr); call CHKERR(ierr)

    nullify(solver)
    solver => ctx_ptr%solver
    if(.not.associated(solver)) call CHKERR(1_mpiint, 'mat shell context has not been set!')

    call PetscObjectViewFromOptions(x, PETSC_NULL_VEC, "-show_shell_x", ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(b, PETSC_NULL_VEC, "-show_shell_b", ierr); call CHKERR(ierr)

    associate( sun => solver%sun, &
               atm => solver%atm, &
               C   => solver%C_dir)

      allocate(coeff(C%dof**2))
      v(0:C%dof-1, 0:C%dof-1) => coeff(1:C%dof**2)

      call DMGetLocalVector(C%da, lx, ierr); call CHKERR(ierr)
      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call VecSet(lb, zero, ierr); call CHKERR(ierr)

      call DMGlobalToLocalBegin(C%da, x, INSERT_VALUES, lx, ierr) ;call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C%da, x, INSERT_VALUES, lx, ierr) ;call CHKERR(ierr)

      call getVecPointer(C%da, lx, xx1d, xx, readonly=.True.)
      call getVecPointer(C%da, lb, xb1d, xb)

      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-1

            if( atm%l1d(atmk(atm,k),i,j) ) then
              do idst = 0, solver%dirtop%dof-1
                xb(idst, k+i1, i, j) = xb(idst, k+i1, i, j) - xx(idst, k, i, j) * atm%a33(atmk(atm,k),i,j)
              enddo
            else

              xinc = sun%xinc(k,i,j)
              yinc = sun%yinc(k,i,j)

              call get_coeff(solver, &
                atm%kabs(atmk(solver%atm,k),i,j), &
                atm%ksca(atmk(solver%atm,k),i,j), &
                atm%g(atmk(solver%atm,k),i,j), &
                atm%dz(atmk(solver%atm,k),i,j), .True., coeff, &
                atm%l1d(atmk(solver%atm,k),i,j), &
                [real(sun%symmetry_phi(k,i,j), irealLUT), real(sun%theta(k,i,j), irealLUT)], &
                lswitch_east=xinc.eq.0, lswitch_north=yinc.eq.0)

              dst = 0
              do idst = 0, solver%dirtop%dof-1
                src = 0
                do isrc = 0, solver%dirtop%dof-1
                   xb(dst, k+i1, i, j) = xb(dst, k+i1, i, j) - xx(src, k, i, j) * v(src, dst)
                   src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                   xb(dst, k+i1, i, j) = xb(dst, k+i1, i, j) - xx(src, k, i+1-xinc, j) * v(src, dst)
                   src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                   xb(dst, k+i1, i, j) = xb(dst, k+i1, i, j) - xx(src, k, i, j+1-yinc) * v(src, dst)
                   src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%dirside%dof-1
                src = 0
                do isrc = 0, solver%dirtop%dof-1
                   xb(dst, k, i+xinc, j) = xb(dst, k, i+xinc, j) - xx(src, k, i, j) * v(src, dst)
                   src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                   xb(dst, k, i+xinc, j) = xb(dst, k, i+xinc, j) - xx(src, k, i+1-xinc, j) * v(src, dst)
                   src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                   xb(dst, k, i+xinc, j) = xb(dst, k, i+xinc, j) - xx(src, k, i, j+1-yinc) * v(src, dst)
                   src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%dirside%dof-1
                src = 0
                do isrc = 0, solver%dirtop%dof-1
                   xb(dst, k, i, j+yinc) = xb(dst, k, i, j+yinc) - xx(src, k, i, j) * v(src, dst)
                   src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                   xb(dst, k, i, j+yinc) = xb(dst, k, i, j+yinc) - xx(src, k, i+1-xinc, j) * v(src, dst)
                   src = src+1
                enddo
                do isrc = 0, solver%dirside%dof-1
                   xb(dst, k, i, j+yinc) = xb(dst, k, i, j+yinc) - xx(src, k, i, j+1-yinc) * v(src, dst)
                   src = src+1
                enddo
                dst = dst+1
              enddo
            endif

          enddo
        enddo
      enddo

      call restoreVecPointer(C%da, lb, xb1d, xb)
      call restoreVecPointer(C%da, lx, xx1d, xx, readonly=.True.)
      call DMRestoreLocalVector(C%da, lx, ierr); call CHKERR(ierr)

      call VecSet(b, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(C%da, lb, ADD_VALUES, b, ierr) ;call CHKERR(ierr)
      call DMLocalToGlobalEnd  (C%da, lb, ADD_VALUES, b, ierr) ;call CHKERR(ierr)

      ! vec copy simulates diagonal entries
      call VecAXPY(b, one, x, ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
    end associate
    call PetscObjectViewFromOptions(b, PETSC_NULL_VEC, "-show_shell_rb", ierr); call CHKERR(ierr)

    ierr = 0
  end subroutine

  subroutine op_mat_mult_ediff(A, x, b, ierr)
    type(tMat), intent(in) :: A
    type(tVec), intent(in) :: x
    type(tVec), intent(inout) :: b
    integer(mpiint), intent(out) :: ierr

    type(t_pprts_shell_ctx), pointer :: ctx_ptr
    class(t_solver), pointer :: solver
    real(ireals),pointer,dimension(:,:,:,:) :: xx=>null(), xb=>null()
    real(ireals),pointer,dimension(:) :: xx1d=>null(), xb1d=>null()
    real(irealLUT), target, allocatable :: coeff(:)
    real(irealLUT), pointer :: v(:,:) ! dim(src,dst)
    integer(iintegers) :: k,i,j
    integer(iintegers) :: idst, isrc, dst, src
    integer(iintegers) :: mdst, msrc
    integer(iintegers), allocatable :: row(:,:), col(:,:)
    type(tVec) :: lb, lx

    nullify(ctx_ptr)
    call matshellgetcontext(A, ctx_ptr, ierr); call CHKERR(ierr)

    nullify(solver)
    solver => ctx_ptr%solver
    if(.not.associated(solver)) call CHKERR(1_mpiint, 'mat shell context has not been set!')

    call PetscObjectViewFromOptions(x, PETSC_NULL_VEC, "-show_shell_x", ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(b, PETSC_NULL_VEC, "-show_shell_b", ierr); call CHKERR(ierr)

    associate( atm => solver%atm, &
               C   => solver%C_diff)

      allocate(coeff(C%dof**2))
      v(0:C%dof-1, 0:C%dof-1) => coeff(1:C%dof**2)

      allocate(row(4,0:C%dof-1), col(4,0:C%dof-1))

      call DMGetLocalVector(C%da, lx, ierr); call CHKERR(ierr)
      call DMGetLocalVector(C%da, lb, ierr); call CHKERR(ierr)
      call VecSet(lb, zero, ierr); call CHKERR(ierr)

      call DMGlobalToLocalBegin(C%da, x, INSERT_VALUES, lx, ierr) ;call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C%da, x, INSERT_VALUES, lx, ierr) ;call CHKERR(ierr)

      call getVecPointer(C%da, lx, xx1d, xx, readonly=.True.)
      call getVecPointer(C%da, lb, xb1d, xb)

      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-1

            if( atm%l1d(atmk(atm,k),i,j) ) then

              do idst = 0, solver%difftop%dof-1
                if (solver%difftop%is_inward(i1+idst)) then ! edn
                  xb(idst, k+i1, i, j) = xb(idst, k+i1, i, j) - xx(idst, k   , i, j) * atm%a11(atmk(atm,k),i,j)
                  xb(idst, k+i1, i, j) = xb(idst, k+i1, i, j) - xx(inv_dof(idst), k+i1, i, j) * atm%a12(atmk(atm,k),i,j)
                else ! eup
                  xb(idst, k   , i, j) = xb(idst, k   , i, j) - xx(idst, k+i1, i, j) * atm%a11(atmk(atm,k),i,j)
                  xb(idst, k   , i, j) = xb(idst, k   , i, j) - xx(inv_dof(idst), k   , i, j) * atm%a12(atmk(atm,k),i,j)
                endif
              enddo

            else

              call PetscLogEventBegin(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)
              call get_coeff(solver, &
                atm%kabs(atmk(solver%atm,k),i,j), &
                atm%ksca(atmk(solver%atm,k),i,j), &
                atm%g(atmk(solver%atm,k),i,j), &
                atm%dz(atmk(solver%atm,k),i,j), .False., coeff, &
                atm%l1d(atmk(solver%atm,k),i,j))
              call PetscLogEventEnd(solver%logs%get_coeff_diff2diff, ierr); call CHKERR(ierr)

              dst = 0
              do idst = 0, solver%difftop%dof-1
                mdst = merge(k+1, k, solver%difftop%is_inward(i1+idst))
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  xb(dst, mdst, i, j) = xb(dst, mdst, i, j) - xx(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  xb(dst, mdst, i, j) = xb(dst, mdst, i, j) - xx(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  xb(dst, mdst, i, j) = xb(dst, mdst, i, j) - xx(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%diffside%dof-1
                mdst = merge(i+1, i, solver%diffside%is_inward(i1+idst))
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  xb(dst, k, mdst, j) = xb(dst, k, mdst, j) - xx(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  xb(dst, k, mdst, j) = xb(dst, k, mdst, j) - xx(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  xb(dst, k, mdst, j) = xb(dst, k, mdst, j) - xx(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

              do idst = 0, solver%diffside%dof-1
                mdst = merge(j+1, j, solver%diffside%is_inward(i1+idst))
                src = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  xb(dst, k, i, mdst) = xb(dst, k, i, mdst) - xx(src, msrc, i, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  xb(dst, k, i, mdst) = xb(dst, k, i, mdst) - xx(src, k, msrc, j) * v(src, dst)
                  src = src+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  xb(dst, k, i, mdst) = xb(dst, k, i, mdst) - xx(src, k, i, msrc) * v(src, dst)
                  src = src+1
                enddo
                dst = dst+1
              enddo

            endif
          enddo

          ! Albedo:
          do idst = 0, solver%difftop%dof-1
            if (.not.solver%difftop%is_inward(i1+idst)) then ! Eup
              xb(idst, C%ze, i, j) = xb(idst, C%ze, i, j) - xx(inv_dof(idst), C%ze, i, j) * solver%atm%albedo(i,j)
            endif
          enddo

        enddo
      enddo

      call restoreVecPointer(C%da, lb, xb1d, xb)
      call restoreVecPointer(C%da, lx, xx1d, xx, readonly=.True.)
      call DMRestoreLocalVector(C%da, lx, ierr); call CHKERR(ierr)

      call VecSet(b, zero, ierr); call CHKERR(ierr)
      call DMLocalToGlobalBegin(C%da, lb, ADD_VALUES, b, ierr) ;call CHKERR(ierr)
      call DMLocalToGlobalEnd  (C%da, lb, ADD_VALUES, b, ierr) ;call CHKERR(ierr)

      ! vec copy simulates diagonal entries
      call VecAXPY(b, one, x, ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(C%da, lb, ierr); call CHKERR(ierr)
    end associate
    call PetscObjectViewFromOptions(b, PETSC_NULL_VEC, "-show_shell_rb", ierr); call CHKERR(ierr)

    ierr = 0
  contains
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
  end subroutine

  !> @brief call the matrix assembly and petsc solve routines for pprts solvers
  subroutine pprts(solver, edirTOA, solution, opt_buildings)
    class(t_solver), intent(inout) :: solver
    real(ireals), intent(in) :: edirTOA
    type(t_state_container), intent(inout) :: solution
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    logical :: lflg, lskip_diffuse_solve, lshell_pprts
    real(ireals) :: b_norm, rtol, atol
    integer(mpiint) :: ierr

    lshell_pprts = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-pprts_shell", lshell_pprts, lflg,ierr) ; call CHKERR(ierr)

    ! --------- scale from [W/m**2] to [W] -----------------
    call scale_flx(solver, solution, lWm2=.False. )

    ! ---------------------------- Edir  -------------------
    if( solution%lsolar_rad ) then
      call PetscLogEventBegin(solver%logs%compute_Edir, ierr)

      call setup_incSolar(solver, solver%incSolar, edirTOA)

      if(lshell_pprts) then
        call setup_matshell(solver, solver%C_dir, solver%Mdir, op_mat_mult_edir)
      else
        call init_Matrix(solver, solver%C_dir, solver%Mdir, setup_direct_preallocation)
        call PetscLogEventBegin(solver%logs%setup_Mdir, ierr)
        call set_dir_coeff(solver, solver%sun, solver%Mdir, solver%C_dir, opt_buildings)
        call PetscLogEventEnd(solver%logs%setup_Mdir, ierr)
        if(ldebug) call mat_info(solver%comm, solver%Mdir)
      endif

      call PetscLogEventBegin(solver%logs%solve_Mdir, ierr)
      call setup_ksp(solver, solver%ksp_solar_dir, solver%C_dir, solver%Mdir, "solar_dir_")
      call solve(solver, &
        solver%ksp_solar_dir, &
        solver%incSolar, &
        solution%edir, &
        solution%Niter_dir, &
        solution%dir_ksp_residual_history)

      call PetscLogEventEnd(solver%logs%solve_Mdir, ierr)

      solution%lchanged=.True.
      solution%lWm2_dir=.False.
      call PetscObjectSetName(solution%edir,'debug_edir',ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_debug_edir", ierr); call CHKERR(ierr)
      call PetscLogEventEnd(solver%logs%compute_Edir, ierr)
    endif


    ! ---------------------------- Source Term -------------
    call PetscLogEventBegin(solver%logs%compute_Ediff, ierr); call CHKERR(ierr)
    call setup_b(solver, solution, solver%b, opt_buildings)

    if( solution%lsolar_rad ) then
      call determine_ksp_tolerances(solver%C_diff, solver%atm%l1d, rtol, atol, solver%ksp_solar_diff)
    else
      call determine_ksp_tolerances(solver%C_diff, solver%atm%l1d, rtol, atol, solver%ksp_thermal_diff)
    endif
    call VecNorm(solver%b, NORM_1, b_norm, ierr); call CHKERR(ierr)
    lskip_diffuse_solve = b_norm.lt.atol
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-skip_diffuse_solve", lskip_diffuse_solve, lflg , ierr) ;call CHKERR(ierr)

    if(lskip_diffuse_solve) then
      if(ldebug) &
        & call CHKWARN(1_mpiint, 'Skipping diffuse solver :: b_norm='//toStr(b_norm))
      call VecCopy(solver%b, solution%ediff, ierr); call CHKERR(ierr)
      solution%Niter_diff = i0
      solution%diff_ksp_residual_history = atol
    else
      ! ---------------------------- Ediff -------------------
      if(lshell_pprts) then
        call setup_matshell(solver, solver%C_diff, solver%Mdiff, op_mat_mult_ediff)
      else
        call init_Matrix(solver, solver%C_diff, solver%Mdiff, setup_diffuse_preallocation)

        call PetscLogEventBegin(solver%logs%setup_Mdiff, ierr)
        call set_diff_coeff(solver, solver%Mdiff, solver%C_diff, opt_buildings)
        call PetscLogEventEnd(solver%logs%setup_Mdiff, ierr)
        if(ldebug) call mat_info(solver%comm, solver%Mdiff)
      endif

      call PetscLogEventBegin(solver%logs%solve_Mdiff, ierr)
      if( solution%lsolar_rad ) then
        call setup_ksp(solver, solver%ksp_solar_diff, solver%C_diff, solver%Mdiff, "solar_diff_")
        call solve(solver, &
          solver%ksp_solar_diff, &
          solver%b, &
          solution%ediff, &
          solution%Niter_diff, &
          solution%diff_ksp_residual_history)
      else
        call setup_ksp(solver, solver%ksp_thermal_diff, solver%C_diff, solver%Mdiff, "thermal_diff_")
        call solve(solver, &
          solver%ksp_thermal_diff, &
          solver%b, &
          solution%ediff, &
          solution%Niter_diff, &
          solution%diff_ksp_residual_history)
      endif
      call PetscLogEventEnd(solver%logs%solve_Mdiff, ierr)
    endif

    solution%lchanged=.True.
    solution%lWm2_diff=.False. !Tenstream solver returns fluxes as [W]

    call PetscLogEventEnd(solver%logs%compute_Ediff, ierr)
  end subroutine

  !> @brief renormalize fluxes with the size of a face(sides or lid)
  subroutine scale_flx(solver, solution, lWm2)
    class(t_solver), intent(inout)        :: solver
    type(t_state_container),intent(inout) :: solution   !< @param solution container with computed fluxes
    logical,intent(in)                    :: lWm2  !< @param determines direction of scaling, if true, scale to W/m**2
    integer(mpiint) :: ierr

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
      integer(iintegers)   :: i, j, k, d, iside, ak
      real(ireals)         :: Ax, Ay, Az, fac

      associate( atm => solver%atm )

      if(solver%myid.eq.0.and.ldebug) print *,'rescaling direct fluxes',C%zm,C%xm,C%ym
      call getVecPointer(C%da, v, xv1d, xv)

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
            ak = atmk(atm, k)

            ! First the faces in x-direction
            Ax = solver%atm%dy*solver%atm%dz(ak,i,j) / real(solver%dirside%area_divider, ireals)
            fac = Ax
            do iside=0,solver%dirside%dof-1
              d = solver%dirtop%dof + iside
              xv(d,k,i,j) = fac
            enddo

            ! Then the rest of the faces in y-direction
            Ay = atm%dy*atm%dz(ak,i,j) / real(solver%dirside%area_divider, ireals)
            fac = Ay
            do iside=0,solver%dirside%dof-1
              d = solver%dirtop%dof + solver%dirside%dof + iside
              xv(d,k,i,j) = fac
            enddo
          enddo
          ! the side faces underneath the surface are always scaled by unity
          xv(solver%dirtop%dof:ubound(xv,dim=1),k,i,j) = one
        enddo
      enddo
      call restoreVecPointer(C%da, v, xv1d, xv)
      end associate
    end subroutine
    subroutine gen_scale_diff_flx_vec(solver, v, C)
      class(t_solver)      :: solver
      type(tVec)           :: v
      type(t_coord)        :: C
      real(ireals),pointer :: xv(:,:,:,:)=>null()
      real(ireals),pointer :: xv1d(:)=>null()

      integer(iintegers)  :: iside, src, i, j, k, ak
      real(ireals)        :: Az, Ax, Ay, fac

      if(solver%myid.eq.0.and.ldebug) print *,'rescaling fluxes',C%zm,C%xm,C%ym
      call getVecPointer(C%da, v, xv1d, xv)

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
            ak = atmk(solver%atm, k)
            ! faces in x-direction
            Ax = solver%atm%dy*solver%atm%dz(ak,i,j) / real(solver%diffside%area_divider, ireals)
            fac = Ax

            do iside=1,solver%diffside%dof
              src = solver%difftop%dof + iside -1
              xv(src ,k,i,j) = fac
            enddo

            ! faces in y-direction
            Ay = solver%atm%dx*solver%atm%dz(ak,i,j) / real(solver%difftop%area_divider, ireals)
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
      call restoreVecPointer(C%da, v, xv1d, xv )

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
    integer(mpiint) :: ierr

    if( .not. solution%lset ) call CHKERR(1_mpiint, 'cant restore solution that was not initialized')
    if( .not. allocated(solution%abso) ) call CHKERR(1_mpiint, 'cant restore solution that was not initialized')

    if( .not. solution%lchanged ) return

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

    logical :: by_flx_divergence, lflg
    integer(mpiint) :: ierr

    if(allocated(solution%edir)) then
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_flxdiv_edir", ierr); call CHKERR(ierr)
    endif
    call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-show_flxdiv_ediff", ierr); call CHKERR(ierr)

    if(.not.allocated(solver%abso_scalevec)) &
      & call gen_abso_scalevec()

    if( (solution%lsolar_rad.eqv..False.) .and. lcalc_nca ) then ! if we should calculate NCA (Klinger), we can just return afterwards
      call scale_flx(solver, solution, lWm2=.True.)
      call nca_wrapper(solver, solution%ediff, solution%abso)
      return
    endif

    if(solver%myid.eq.0.and.ldebug) print *,'Calculating flux divergence solar?',solution%lsolar_rad,'NCA?',lcalc_nca
    call VecSet(solution%abso,zero,ierr) ;call CHKERR(ierr)

    ! make sure to bring the fluxes into [W] before the absorption calculation
    call scale_flx(solver, solution, lWm2=.False.)

    ! if there are no 3D layers globally, we should skip the ghost value copying....
    !lhave_no_3d_layer = mpi_logical_and(solver%comm, all(atm%l1d.eqv..True.))
    if(ltwostr_only) then
      call compute_1D_absorption()
      goto 99
    endif

    by_flx_divergence = .True.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-absorption_by_flx_divergence",&
      by_flx_divergence, lflg, ierr) ;call CHKERR(ierr)

    if(by_flx_divergence) then
      call compute_absorption_by_flx_divergence()
    else
      call compute_absorption_by_coeff_divergence()
    endif

    99 continue ! cleanup
    call VecPointwiseMult(solution%abso, solution%abso, solver%abso_scalevec, ierr); call CHKERR(ierr)

  contains

    subroutine compute_1D_absorption()
      associate(                    &
          C_dir   => solver%C_dir,  &
          C_diff  => solver%C_diff, &
          C_one   => solver%C_one   )
        if(ldebug.and.solver%myid.eq.0) print *,'lhave_no_3d_layer => will use 1D absorption computation'

        if(solution%lsolar_rad) call getVecPointer(C_dir%da, solution%edir, xedir1d, xedir)

        call getVecPointer(C_diff%da, solution%ediff, xediff1d, xediff, readonly=.True.)
        call getVecPointer(C_one%da , solution%abso , xabso1d, xabso)

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

        if(solution%lsolar_rad) call restoreVecPointer(C_dir%da, solution%edir, xedir1d, xedir )

        call restoreVecPointer(C_diff%da, solution%ediff, xediff1d, xediff, readonly=.True.)
        call restoreVecPointer(C_one%da , solution%abso, xabso1d ,xabso)
      end associate
    end subroutine

    subroutine compute_absorption_by_coeff_divergence()
      real(ireals) :: cdiv
      real(irealLUT), target, allocatable :: cdir2dir(:), cdir2diff(:), cdiff2diff(:)
      real(irealLUT), pointer :: dir2dir(:,:) ! dim(dst, src)
      real(irealLUT), pointer :: dir2diff(:,:) ! dim(dst, src)
      real(irealLUT), pointer :: diff2diff(:,:) ! dim(dst, src)
      integer(iintegers) :: xinc, yinc, idof, msrc

      call getVecPointer(solver%C_one%da, solution%abso, xabso1d, xabso)

      associate(                    &
          sun     => solver%sun,    &
          atm     => solver%atm,    &
          C_dir   => solver%C_dir,  &
          C_diff  => solver%C_diff, &
          C_one   => solver%C_one   )

        allocate(cdir2dir(C_dir%dof**2))
        dir2dir(0:C_dir%dof-1, 0:C_dir%dof-1) => cdir2dir(1:C_dir%dof**2)

        allocate(cdir2diff(C_dir%dof * C_diff%dof))
        dir2diff(0:C_dir%dof-1, 0:C_diff%dof-1) => cdir2diff(1:C_dir%dof*C_diff%dof)

        if(solution%lsolar_rad) then
          ! Copy ghosted values for direct vec
          call DMGetLocalVector(C_dir%da ,ledir ,ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C_dir%da, solution%edir , INSERT_VALUES, ledir ,ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd  (C_dir%da, solution%edir , INSERT_VALUES, ledir ,ierr); call CHKERR(ierr)
          call getVecPointer(C_dir%da, ledir, xedir1d ,xedir, readonly=.True.)

          do j=C_one%ys,C_one%ye
            do i=C_one%xs,C_one%xe
              do k=C_one%zs,C_one%ze

                ! Divergence = Incoming - Outgoing

                if(atm%l1d(atmk(atm, k),i,j)) then ! one dimensional i.e. twostream

                  cdiv = max(zero, one - atm%a33(atmk(atm,k),i,j) - atm%a13(atmk(atm,k),i,j) - atm%a23(atmk(atm,k),i,j))
                  do isrc = 0, solver%dirtop%dof-1
                    xabso(i0,k,i,j) = xabso(i0,k,i,j) + xedir(isrc, k, i, j) * cdiv
                  enddo

                else

                  xinc = sun%xinc(k,i,j)
                  yinc = sun%yinc(k,i,j)

                  call get_coeff(solver, &
                    atm%kabs(atmk(solver%atm,k),i,j), &
                    atm%ksca(atmk(solver%atm,k),i,j), &
                    atm%g(atmk(solver%atm,k),i,j), &
                    atm%dz(atmk(solver%atm,k),i,j), .True., cdir2dir, &
                    atm%l1d(atmk(solver%atm,k),i,j), &
                    [real(sun%symmetry_phi(k,i,j), irealLUT), real(sun%theta(k,i,j), irealLUT)], &
                    lswitch_east=xinc.eq.0, lswitch_north=yinc.eq.0)

                  call get_coeff(solver, &
                    atm%kabs(atmk(solver%atm,k),i,j), &
                    atm%ksca(atmk(solver%atm,k),i,j), &
                    atm%g(atmk(solver%atm,k),i,j), &
                    atm%dz(atmk(solver%atm,k),i,j), .False., cdir2diff, &
                    atm%l1d(atmk(solver%atm,k),i,j), &
                    [real(sun%symmetry_phi(k,i,j), irealLUT), real(sun%theta(k,i,j), irealLUT)], &
                    lswitch_east=xinc.eq.0, lswitch_north=yinc.eq.0)

                  idof = 0
                  do isrc = 0, solver%dirtop%dof-1
                    cdiv = one - sum(dir2dir(isrc,:)) - sum(dir2diff(isrc,:))
                    xabso(i0,k,i,j) = xabso(i0,k,i,j) + xedir(isrc, k, i, j) * cdiv
                    idof = idof+1
                  enddo
                  do isrc = 0, solver%dirside%dof-1
                    cdiv = one - sum(dir2dir(idof,:)) - sum(dir2diff(idof,:))
                    xabso(i0,k,i,j) = xabso(i0,k,i,j) + xedir(idof, k, i+1-xinc, j) * cdiv
                    idof = idof+1
                  enddo
                  do isrc = 0, solver%dirside%dof-1
                    cdiv = one - sum(dir2dir(idof,:)) - sum(dir2diff(idof,:))
                    xabso(i0,k,i,j) = xabso(i0,k,i,j) + xedir(idof, k, i, j+1-yinc) * cdiv
                    idof = idof+1
                  enddo

                endif

              enddo
            enddo
          enddo

          call restoreVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.True.)
          call DMRestoreLocalVector(C_dir%da, ledir, ierr) ; call CHKERR(ierr)
        endif
      end associate

      associate(                    &
          atm     => solver%atm,    &
          C_diff  => solver%C_diff, &
          C_one   => solver%C_one   )

        ! Copy ghosted values for diffuse vec
        call DMGetLocalVector(C_diff%da,lediff,ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin (C_diff%da, solution%ediff  , INSERT_VALUES, lediff,ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd   (C_diff%da, solution%ediff  , INSERT_VALUES, lediff,ierr); call CHKERR(ierr)
        call getVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.True.)

        allocate(cdiff2diff(C_diff%dof**2))
        diff2diff(0:C_diff%dof-1, 0:C_diff%dof-1) => cdiff2diff(1:C_diff%dof**2)

        do j=C_one%ys,C_one%ye
          do i=C_one%xs,C_one%xe
            do k=C_one%zs,C_one%ze

              ! Divergence = Incoming - Outgoing

              if(atm%l1d(atmk(atm, k),i,j)) then ! one dimensional i.e. twostream

                cdiv = max(zero, one - atm%a11(atmk(atm,k),i,j) - atm%a12(atmk(atm,k),i,j))
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  xabso(i0,k,i,j) = xabso(i0,k,i,j) + xediff(isrc, msrc, i, j) * cdiv
                enddo

              else

                call get_coeff(solver, &
                  atm%kabs(atmk(solver%atm,k),i,j), &
                  atm%ksca(atmk(solver%atm,k),i,j), &
                  atm%g(atmk(solver%atm,k),i,j), &
                  atm%dz(atmk(solver%atm,k),i,j), .False., cdiff2diff, &
                  atm%l1d(atmk(solver%atm,k),i,j))

                idof = 0
                do isrc = 0, solver%difftop%dof-1
                  msrc = merge(k, k+1, solver%difftop%is_inward(i1+isrc))
                  cdiv = one - sum(diff2diff(idof,:))
                  xabso(i0,k,i,j) = xabso(i0,k,i,j) + xediff(idof, msrc, i, j) * cdiv
                  idof = idof+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(i, i+1, solver%diffside%is_inward(i1+isrc))
                  cdiv = one - sum(diff2diff(idof,:))
                  xabso(i0,k,i,j) = xabso(i0,k,i,j) + xediff(idof, k, msrc, j) * cdiv
                  idof = idof+1
                enddo
                do isrc = 0, solver%diffside%dof-1
                  msrc = merge(j, j+1, solver%diffside%is_inward(i1+isrc))
                  cdiv = one - sum(diff2diff(idof,:))
                  xabso(i0,k,i,j) = xabso(i0,k,i,j) + xediff(idof, k, i, msrc) * cdiv
                  idof = idof+1
                enddo

              endif

            enddo
          enddo
        enddo

        call restoreVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.True.)
        call DMRestoreLocalVector(C_diff%da, lediff, ierr); call CHKERR(ierr)
      end associate

      call restoreVecPointer(solver%C_one%da, solution%abso, xabso1d ,xabso)
    end subroutine

    subroutine compute_absorption_by_flx_divergence()
      associate(                    &
          atm     => solver%atm,    &
          C_dir   => solver%C_dir,  &
          C_diff  => solver%C_diff, &
          C_one   => solver%C_one   )

        if(solution%lsolar_rad) then
          ! Copy ghosted values for direct vec
          call DMGetLocalVector(C_dir%da ,ledir ,ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C_dir%da, solution%edir , INSERT_VALUES, ledir ,ierr); call CHKERR(ierr)
          call DMGlobalToLocalEnd  (C_dir%da, solution%edir , INSERT_VALUES, ledir ,ierr); call CHKERR(ierr)
          call getVecPointer(C_dir%da, ledir, xedir1d ,xedir, readonly=.True.)
        endif

        ! Copy ghosted values for diffuse vec
        call DMGetLocalVector(C_diff%da,lediff,ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin (C_diff%da, solution%ediff  , INSERT_VALUES, lediff,ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd   (C_diff%da, solution%ediff  , INSERT_VALUES, lediff,ierr); call CHKERR(ierr)
        call getVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.True.)

        call getVecPointer(C_one%da, solution%abso, xabso1d, xabso)

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
          call restoreVecPointer(C_dir%da, ledir, xedir1d, xedir, readonly=.True.)
          call DMRestoreLocalVector(C_dir%da, ledir, ierr) ; call CHKERR(ierr)
        endif

        call restoreVecPointer(C_diff%da, lediff, xediff1d, xediff, readonly=.True.)
        call DMRestoreLocalVector(C_diff%da, lediff, ierr) ; call CHKERR(ierr)

        call restoreVecPointer(C_one%da, solution%abso, xabso1d ,xabso)
      end associate
    end subroutine

    subroutine gen_abso_scalevec()
      associate(atm     => solver%atm, C_one   => solver%C_one)

        allocate(solver%abso_scalevec)
        call VecDuplicate(solution%abso, solver%abso_scalevec, ierr); call CHKERR(ierr)
        call getVecPointer(C_one%da, solver%abso_scalevec, xabso1d, xabso)
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

        call restoreVecPointer(C_one%da, solver%abso_scalevec, xabso1d ,xabso)
      end associate
    end subroutine
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
    real(ireals), intent(inout), optional :: ksp_residual_history(:)

    KSPConvergedReason :: reason

    character(len=default_str_len) :: prefix
    KSPType :: old_ksp_type

    logical :: lskip_ksp_solve, laccept_incomplete_solve, lflg
    integer(mpiint) :: ierr

    if(solver%myid.eq.0.and.ldebug) print *,'Solving Matrix'

    if(present(ksp_residual_history)) then
      call KSPSetResidualHistory(ksp, ksp_residual_history, &
        & size(ksp_residual_history, kind=iintegers), PETSC_TRUE, ierr); call CHKERR(ierr)
    endif

    lskip_ksp_solve=.False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
      "-skip_ksp_solve" , lskip_ksp_solve, lflg , ierr) ;call CHKERR(ierr)
    if(lskip_ksp_solve) then
      call VecCopy(b, x, ierr); call CHKERR(ierr)
      return
    endif

    call hegedus_trick(ksp, b, x)
    call KSPSolve(ksp,b,x,ierr); call CHKERR(ierr)
    call KSPGetIterationNumber(ksp,iter,ierr); call CHKERR(ierr)
    call KSPGetConvergedReason(ksp,reason,ierr); call CHKERR(ierr)

    ! if(reason.eq.KSP_DIVERGED_ITS) then
    !   if(solver%myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
    !   return
    ! endif
    laccept_incomplete_solve = .False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-accept_incomplete_solve", laccept_incomplete_solve, lflg, ierr); call CHKERR(ierr)
    if(laccept_incomplete_solve) return

    if(reason.le.0) then
      call KSPGetOptionsPrefix(ksp, prefix, ierr); call CHKERR(ierr)
      if(solver%myid.eq.0) &
        & call CHKWARN(int(reason, mpiint), trim(prefix)//' :: Resetted initial guess to zero and try again with gmres')
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
      call CHKERR(int(reason, mpiint), &
        & '***** SOLVER did NOT converge :( -- '// &
        & 'if you know what you are doing, you can use the option -accept_incomplete_solve to continue')
    endif
  end subroutine

  !> @brief initialize PETSc Krylov Subspace Solver
  !> @details default KSP solver is a FGMRES with BJCAOBI // ILU(1)
  !> \n -- the default does however not scale well -- and we configure petsc solvers per commandline anyway
  !> \n -- see documentation for details on how to do so
  subroutine setup_ksp(solver, ksp, C, A, prefix)
    class(t_solver) :: solver
    type(tKSP), intent(inout), allocatable :: ksp
    type(t_coord) :: C
    type(tMat) :: A
    type(tPC)  :: prec
    logical :: linit

    character(len=*), intent(in), optional :: prefix

    integer(iintegers),parameter  :: maxiter=1000

    integer(mpiint) :: myid, numnodes

    real(ireals) :: rtol, atol

    logical :: prec_is_set
    integer(mpiint) :: ierr
    character(len=default_str_len) :: kspprefix

    linit = allocated(ksp)
    if(linit) return

    call set_dmda_cell_coordinates(solver, solver%atm, C%da, ierr); call CHKERR(ierr)

    call mpi_comm_rank(C%comm, myid, ierr)    ; call CHKERR(ierr)
    call mpi_comm_size(C%comm, numnodes, ierr); call CHKERR(ierr)


    allocate(ksp)
    call KSPCreate(C%comm, ksp, ierr); call CHKERR(ierr)
    if(present(prefix)) then
      call KSPAppendOptionsPrefix(ksp, trim(prefix), ierr); call CHKERR(ierr)
    endif
    call KSPGetOptionsPrefix(ksp, kspprefix, ierr); call CHKERR(ierr)

    call KSPSetType(ksp, KSPBCGS, ierr); call CHKERR(ierr)
    call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr); call CHKERR(ierr)

    prec_is_set = .False.
    call PetscOptionsHasName(PETSC_NULL_OPTIONS, trim(kspprefix), '-pc_type', prec_is_set, ierr); call CHKERR(ierr)

    if(.not.prec_is_set) then
      !call CHKWARN(1_mpiint, 'no preconditioner setting found, applying defaults')
      call KSPGetPC(ksp, prec, ierr); call CHKERR(ierr)
      if(numnodes.eq.0) then
        call PCSetType(prec, PCILU, ierr); call CHKERR(ierr)
      else
        call PCSetType(prec, PCASM, ierr); call CHKERR(ierr)
        call PCASMSetOverlap(prec, i1, ierr); call CHKERR(ierr)
      endif
    endif

    call determine_ksp_tolerances(C, solver%atm%l1d, rtol, atol)
    call KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT_REAL, maxiter, ierr); call CHKERR(ierr)

    call KSPSetConvergenceTest(ksp, MyKSPConverged, 0, PETSC_NULL_FUNCTION, ierr); call CHKERR(ierr)

    call KSPSetOperators(ksp, A, A, ierr); call CHKERR(ierr)
    call KSPSetDM(ksp, C%da, ierr); call CHKERR(ierr)
    call KSPSetDMActive(ksp, PETSC_FALSE, ierr); call CHKERR(ierr)

    call KSPSetFromOptions(ksp, ierr); call CHKERR(ierr)
    call KSPSetUp(ksp, ierr); call CHKERR(ierr)

    if(.not.prec_is_set) then
      default_preconditioner_settings: block
        integer(iintegers) :: i, asm_N, asm_iter
        integer(iintegers) :: first_local
        type(tKSP), allocatable :: asm_ksps(:)
        type(tPC) :: subpc
        logical :: lflg

        asm_iter = 2
        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
          "-ts_ksp_iter", asm_iter, lflg, ierr) ;call CHKERR(ierr)

        call PCASMGetSubKSP(prec, asm_N, first_local, PETSC_NULL_KSP, ierr); call CHKERR(ierr)
        allocate(asm_ksps(asm_N))
        call PCASMGetSubKSP(prec, asm_N, first_local, asm_ksps, ierr); call CHKERR(ierr)
        do i = 1, asm_N
          call KSPSetType(asm_ksps(i), KSPRICHARDSON, ierr); call CHKERR(ierr)
          call KSPRichardsonSetSelfScale(asm_ksps(i), PETSC_TRUE, ierr); call CHKERR(ierr)
          call KSPSetNormType(asm_ksps(i), KSP_NORM_PRECONDITIONED, ierr); call CHKERR(ierr)
          call KSPSetTolerances(asm_ksps(i), &
            & PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
            & asm_iter, ierr); call CHKERR(ierr)
          call KSPGetPC(asm_ksps(i), subpc, ierr); call CHKERR(ierr)
          call PCSetType(subpc, PCSOR, ierr); call CHKERR(ierr)
          call KSPSetFromOptions(asm_ksps(i), ierr); call CHKERR(ierr)
        enddo
      end block default_preconditioner_settings
    endif

    linit = .True.
    if(myid.eq.0.and.ldebug) print *,'Setup KSP done'
  end subroutine


  !> @brief: determine tolerances for solvers
  subroutine determine_ksp_tolerances(C, l1d, rtol, atol, ksp)
    type(t_coord), intent(in) :: C
    logical, intent(in) :: l1d(:,:,:)
    real(ireals), intent(out) :: rtol, atol
    type(tKSP), intent(in), allocatable, optional:: ksp
    real(ireals) :: rel_atol=1e-4_ireals
    real(ireals) :: unconstrained_fraction
    integer(mpiint) :: myid, ierr
    integer(iintegers) :: maxit
    real(ireals) :: dtol

    if(present(ksp)) then
      if(allocated(ksp)) then
        call KSPGetTolerances(ksp, rtol, atol, dtol, maxit, ierr); call CHKERR(ierr)
        if(ldebug) then
          call mpi_comm_rank(C%comm, myid, ierr); call CHKERR(ierr)
          if(myid.eq.0) print *,cstr('Read tolerances from ksp','red'), rtol, atol, dtol, maxit
        endif
        return
      endif
    endif

    rtol = 1e-5_ireals
    unconstrained_fraction = real(count(.not.l1d), ireals) / real(size(l1d), ireals)
    call imp_allreduce_min(C%comm, &
      & rel_atol &
      & * real(C%dof*C%glob_xm*C%glob_ym*C%glob_zm, ireals) &
      & * unconstrained_fraction, atol)
    atol = max(1e-8_ireals, atol)

    if(ldebug) then
      call mpi_comm_rank(C%comm, myid, ierr); call CHKERR(ierr)
      if(myid.eq.0) &
        print *,'KSP ', &
        & '-- tolerances:', rtol, atol, &
        & ':: rel_atol', rel_atol, &
        & ':: total dof', C%dof * C%glob_xm * C%glob_ym * C%glob_zm, &
        & ':: unconstrained fraction', unconstrained_fraction
    endif
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
    if(.False.) dummy = dummy + 1 ! stupid statement to remove unused variable warning
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
    integer(mpiint) :: ierr

    fac = edirTOA * solver%atm%dx*solver%atm%dy / real(solver%dirtop%area_divider, ireals)

    call VecSet(incSolar,zero,ierr) ;call CHKERR(ierr)

    call getVecPointer(solver%C_dir%da, incSolar, x1d, x4d)

    do j=solver%C_dir%ys,solver%C_dir%ye
      do i=solver%C_dir%xs,solver%C_dir%xe
        do src=1, solver%dirtop%dof
          x4d(src-1,solver%C_dir%zs,i,j) = fac * solver%sun%costheta(solver%C_dir%zs,i,j)
        enddo
      enddo
    enddo

    call restoreVecPointer(solver%C_dir%da, incSolar, x1d, x4d)

    if(solver%myid.eq.0 .and. ldebug) print *,solver%myid,'Setup of IncSolar done', edirTOA, &
      & '(', fac*solver%sun%costheta(solver%C_dir%zs,solver%C_dir%xs,solver%C_dir%ys), ')'

  end subroutine

  !> @brief build direct radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_dir_coeff(solver, sun, A, C, opt_buildings)
    class(t_solver), intent(in) :: solver
    type(t_suninfo), intent(in) :: sun
    type(tMat), intent(inout)   :: A
    type(t_coord), intent(in)   :: C
    type(t_pprts_buildings), optional, intent(in) :: opt_buildings

    integer(iintegers) :: i,j,k
    integer(mpiint) :: ierr

    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'setup_direct_matrix ...'

    call MatZeroEntries(A, ierr) ;call CHKERR(ierr)
    call mat_set_diagonal(A)

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

    if(present(opt_buildings)) then
      call MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)
      call MatAssemblyEnd  (A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)
      call set_buildings_coeff(solver, C, opt_buildings, A, ierr); call CHKERR(ierr)
    endif

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, "-show_Mdir", ierr); call CHKERR(ierr)
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
          if( norm.gt.one+10._ireals*epsilon(norm) ) then ! could renormalize
            if( norm.gt.one+10._ireals*sqrt(epsilon(norm)) ) then ! fatally off
              print *,'direct sum(dst==',dst,') gt one',norm
              print *,'direct coeff',norm,'::',v
              call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
            endif ! fatally off
          endif ! could renormalize
        enddo
      endif
    end subroutine

    subroutine set_eddington_coeff(atm,A,k,i,j)
      type(t_atmosphere), intent(in) :: atm
      type(tMat),intent(inout)       :: A
      integer(iintegers),intent(in)  :: i,j,k

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

    !> @brief   apply blocking of direct radiation from buildings
    !> @details Goal: set all src dof on a buildings face towards all dst dof to zero
    !> \n       albedo is not used in the dir2dir case, we only set blocking of radiation
    subroutine set_buildings_coeff(solver, C, opt_buildings, A, ierr)
      class(t_solver)                     :: solver
      type(t_coord),intent(in)            :: C
      type(t_pprts_buildings), intent(in) :: opt_buildings
      type(tMat),intent(inout)            :: A
      integer(mpiint), intent(out)        :: ierr

      MatStencil         :: row(4,0:C%dof-1)  ,col(4,1)
      integer(iintegers) :: m, isrc, idst, dst, idx(4)
      integer(iintegers) :: xinc, yinc, zinc
      real(ireals) :: v(C%dof)

      v(:) = -zero

      ierr = 0

      associate( B => opt_buildings )
        do m = 1, size(B%iface)
          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) -1 + [C%zs, C%xs, C%ys]

          associate(k => idx(2), i => idx(3), j => idx(4))

            !print *, m, 'face', B%iface(m), 'idx', idx, 'Ag', B%albedo(m), 'x/y inc', sun%xinc(k,i,j), sun%yinc(k,i,j)
            !lsun_east  = sun%xinc(idx(k,i,j).eq.i0
            !lsun_north = sun%yinc(idx(k,i,j).eq.i0

            xinc = sun%xinc(k,i,j)
            yinc = sun%yinc(k,i,j)
            zinc = 0
            if(idx(1).eq.PPRTS_BOT_FACE) zinc = 1

            dst = 0
            do idst = 0, solver%dirtop%dof-1
              row(MatStencil_j,dst) = i
              row(MatStencil_k,dst) = j
              row(MatStencil_i,dst) = k+1+zinc
              row(MatStencil_c,dst) = dst
              dst = dst + 1
            enddo

            do idst = 1, solver%dirside%dof
              row(MatStencil_j,dst) = i+xinc
              row(MatStencil_k,dst) = j
              row(MatStencil_i,dst) = k+zinc
              row(MatStencil_c,dst) = dst ! Define transmission towards the left/right lid
              dst = dst + 1
            enddo

            do idst = 1, solver%dirside%dof
              row(MatStencil_j,dst) = i
              row(MatStencil_k,dst) = j+yinc
              row(MatStencil_i,dst) = k+zinc
              row(MatStencil_c,dst) = dst ! Define transmission towards the front/back lid
              dst = dst + 1
            enddo

            select case(idx(1))
            case(PPRTS_TOP_FACE, PPRTS_BOT_FACE)
              do isrc = 0, solver%dirtop%dof-1
                col(MatStencil_j,1) = i
                col(MatStencil_k,1) = j
                col(MatStencil_i,1) = k+zinc
                col(MatStencil_c,1) = isrc
                call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
              enddo

            case(PPRTS_LEFT_FACE, PPRTS_RIGHT_FACE)
              do isrc = solver%dirtop%dof, solver%dirtop%dof + solver%dirside%dof - 1
                col(MatStencil_j,1) = i+1-xinc
                col(MatStencil_k,1) = j
                col(MatStencil_i,1) = k+zinc
                col(MatStencil_c,1) = isrc
                call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
              enddo

            case(PPRTS_REAR_FACE, PPRTS_FRONT_FACE)
              do isrc = solver%dirtop%dof + solver%dirside%dof, solver%dirtop%dof + solver%dirside%dof + solver%dirside%dof - 1
                col(MatStencil_j,1) = i
                col(MatStencil_k,1) = j+1-yinc
                col(MatStencil_i,1) = k+zinc
                col(MatStencil_c,1) = isrc
                call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
              enddo

            case default
              call CHKERR(1_mpiint, 'wrong building fidx '//toStr(idx(1)))
            end select
          end associate
        enddo
      end associate
    end subroutine

  end subroutine set_dir_coeff

  !> @brief setup source term for diffuse radiation
  !> @details this is either direct radiation scattered into one of the diffuse coeffs:
  !> \n direct source term is
  !> \n   direct radiation times the dir2diff coeffs
  !> \n or it may be that we have a source term due to thermal emission --
  !> \n   to determine emissivity of box, we use the forward transport coefficients backwards
  !> \n   a la: transmissivity $T = \sum(coeffs)$ and therefore emissivity $E = 1 - T$
  subroutine setup_b(solver, solution, b, opt_buildings)
    class(t_solver)        , intent(in)    :: solver
    type(t_state_container), intent(in)    :: solution
    type(tVec)             , intent(inout) :: b
    type(t_pprts_buildings), intent(in), optional :: opt_buildings

    type(tVec) :: local_b
    type(tVec) :: local_edir

    real(ireals),pointer,dimension(:,:,:,:) :: xsrc=>null()
    real(ireals),pointer,dimension(:) :: xsrc1d=>null()
    integer(mpiint) :: ierr

    associate(  atm     => solver%atm, &
                C_dir   => solver%C_dir, &
                C_diff  => solver%C_diff)

    if(solver%myid.eq.0.and.ldebug) print *,'src Vector Assembly...'

    call DMGetLocalVector(C_diff%da,local_b,ierr) ;call CHKERR(ierr)
    call VecSet(local_b,zero,ierr) ;call CHKERR(ierr)

    call getVecPointer(C_diff%da, local_b, xsrc1d, xsrc)

    if(solution%lsolar_rad) then
      ! Copy ghosted values for direct vec
      call DMGetLocalVector(C_dir%da, local_edir, ierr); call CHKERR(ierr)
      call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
      call DMGlobalToLocalBegin(C_dir%da, solution%edir, ADD_VALUES, local_edir ,ierr); call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C_dir%da, solution%edir, ADD_VALUES, local_edir ,ierr); call CHKERR(ierr)

      call set_solar_source(solver, local_edir)
    endif

    if(solution%lthermal_rad) &
      call set_thermal_source()

    if(present(opt_buildings)) then
      if(solution%lsolar_rad) then
        call set_buildings_reflection(solver, local_edir, opt_buildings)
      endif

      if(solution%lthermal_rad) then
        call set_buildings_emission(solver, opt_buildings, ierr); call CHKERR(ierr)
      endif
    endif

    if(solver%myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients ...done'

    if(solution%lsolar_rad) then
      call DMRestoreLocalVector(C_dir%da, local_edir, ierr); call CHKERR(ierr)
    endif

    call restoreVecPointer(C_diff%da, local_b, xsrc1d, xsrc )

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

      if(.not.allocated(atm%planck)) then
        ierr = 1
        call CHKERR(ierr, 'have thermal computation but planck data is not allocated')
      endif
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

    subroutine set_solar_source(solver, local_edir)
      class(t_solver), intent(in) :: solver
      type(tVec)     , intent(in) :: local_edir

      real(irealLUT)      :: dir2diff(solver%C_dir%dof*solver%C_diff%dof)
      real(ireals)        :: solrad
      integer(iintegers)  :: idof, idofdst, idiff
      integer(iintegers)  :: k, i, j, src, dst

      real(ireals),pointer,dimension(:,:,:,:) :: xedir=>null()
      real(ireals),pointer,dimension(:)       :: xedir1d=>null()

      logical :: lsun_east,lsun_north

      associate(atm     => solver%atm, &
                C_dir   => solver%C_dir, &
                C_diff  => solver%C_diff, &
                sun     => solver%sun)

      call getVecPointer(C_dir%da, local_edir, xedir1d, xedir)

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

      call restoreVecPointer(C_dir%da, local_edir, xedir1d, xedir )
      end associate
    end subroutine

    subroutine set_buildings_reflection(solver, local_edir, buildings)
      class(t_solver)        , intent(in)   :: solver
      type(tVec)             , intent(in)   :: local_edir
      type(t_pprts_buildings), intent(in)   :: buildings

      real(ireals),pointer,dimension(:,:,:,:) :: xedir=>null()
      real(ireals),pointer,dimension(:)       :: xedir1d=>null()

      integer(iintegers) :: m, idx(4), isrc, idiff, dof_offset

      associate(                           &
          & B => buildings,                &
          & C => solver%C_dir              )

        if(solution%lsolar_rad) then
          call getVecPointer(C%da, local_edir, xedir1d, xedir)
          do m = 1, size(B%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) -1 + [C%zs, C%xs, C%ys]

            associate(k => idx(2), i => idx(3), j => idx(4))

              select case(idx(1))
              case(PPRTS_TOP_FACE)
                do idiff = 0, solver%difftop%dof-1
                  if (.not.solver%difftop%is_inward(i1+idiff)) xsrc(idiff,k,i,j) = 0
                enddo
                do isrc = 0, solver%dirtop%dof-1
                  do idiff = 0, solver%difftop%dof-1
                    if (.not.solver%difftop%is_inward(i1+idiff)) then ! eup on upper face
                      xsrc(idiff,k,i,j) = xsrc(idiff,k,i,j) &
                        & + xedir(isrc,k,i,j) * B%albedo(m) / real(solver%difftop%streams, ireals)
                    endif
                  enddo
                enddo

              case(PPRTS_BOT_FACE)
                do idiff = 0, solver%difftop%dof-1
                  if (solver%difftop%is_inward(i1+idiff)) xsrc(idiff,k+1,i,j) = 0
                enddo
                do isrc = 0, solver%dirtop%dof-1
                  do idiff = 0, solver%difftop%dof-1
                    if (solver%difftop%is_inward(i1+idiff)) then ! edn on lower face
                      xsrc(idiff,k+1,i,j) = xsrc(idiff,k+1,i,j) &
                        & + xedir(isrc,k+1,i,j) * B%albedo(m) / real(solver%difftop%streams, ireals)
                    endif
                  enddo
                enddo

              case(PPRTS_LEFT_FACE)
                dof_offset = solver%difftop%dof
                do idiff = 0, solver%diffside%dof-1
                  if (.not.solver%diffside%is_inward(i1+idiff)) xsrc(dof_offset+idiff,k,i,j) = 0
                enddo

                do isrc = solver%dirtop%dof, solver%dirtop%dof + solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof-1
                    if (.not.solver%diffside%is_inward(i1+idiff)) then ! e_left on left face
                      xsrc(dof_offset+idiff,k,i,j) = xsrc(dof_offset+idiff,k,i,j) &
                        & + xedir(isrc,k,i,j) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    endif
                  enddo
                enddo

              case(PPRTS_RIGHT_FACE)
                dof_offset = solver%difftop%dof
                do idiff = 0, solver%diffside%dof-1
                  if (solver%diffside%is_inward(i1+idiff)) xsrc(dof_offset+idiff,k,i+1,j) = 0
                enddo

                do isrc = solver%dirtop%dof, solver%dirtop%dof + solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof-1
                    if (solver%diffside%is_inward(i1+idiff)) then ! e_right on right face
                      xsrc(dof_offset+idiff,k,i+1,j) = xsrc(dof_offset+idiff,k,i+1,j) &
                        & + xedir(isrc,k,i+1,j) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    endif
                  enddo
                enddo

              case(PPRTS_REAR_FACE )
                dof_offset = solver%difftop%dof + solver%diffside%dof
                do idiff = 0, solver%diffside%dof-1
                  if (.not.solver%diffside%is_inward(i1+idiff)) xsrc(dof_offset+idiff,k,i,j) = 0
                enddo

                do isrc = solver%dirtop%dof + solver%dirside%dof, solver%dirtop%dof + 2*solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof-1
                    if (.not.solver%diffside%is_inward(i1+idiff)) then ! e_backward
                      xsrc(dof_offset+idiff,k,i,j) = xsrc(dof_offset+idiff,k,i,j) &
                        & + xedir(isrc,k,i,j) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    endif
                  enddo
                enddo

              case(PPRTS_FRONT_FACE)
                dof_offset = solver%difftop%dof + solver%diffside%dof
                do idiff = 0, solver%diffside%dof-1
                  if (solver%diffside%is_inward(i1+idiff)) xsrc(dof_offset+idiff,k,i,j+1) = 0
                enddo

                do isrc = solver%dirtop%dof + solver%dirside%dof, solver%dirtop%dof + 2*solver%dirside%dof - 1
                  do idiff = 0, solver%diffside%dof-1
                    if (solver%diffside%is_inward(i1+idiff)) then ! e_forward
                      xsrc(dof_offset+idiff,k,i,j+1) = xsrc(dof_offset+idiff,k,i,j+1) &
                        & + xedir(isrc,k,i,j+1) * B%albedo(m) / real(solver%diffside%streams, ireals)
                    endif
                  enddo
                enddo

              case default
                call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1)+1))
              end select
            end associate
          enddo ! loop over building faces
          call restoreVecPointer(C%da, local_edir, xedir1d, xedir)
        endif ! is_solar
      end associate
    end subroutine

    subroutine set_buildings_emission(solver, buildings, ierr)
      class(t_solver)        , intent(in)   :: solver
      type(t_pprts_buildings), intent(in)   :: buildings

      integer(iintegers) :: m, idx(4), idiff, dof_offset, ak
      real(ireals) :: emis, Az, Ax, Ay
      integer(mpiint) :: ierr

      if(.not.allocated(buildings%planck)) then
        ierr = 1
        call CHKERR(ierr, 'tried to set building emissions but planck data is not allocated')
      endif

      associate( atm => solver%atm, B => buildings, C => solver%C_diff )

        Az = atm%dx * atm%dy / real(solver%difftop%area_divider, ireals)

        do m = 1, size(B%iface)

          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) -1 + [C%zs, C%xs, C%ys]

          associate(k => idx(2), i => idx(3), j => idx(4))

            ak = atmk(atm, k)
            Ax = atm%dy * atm%dz(ak,i,j) / real(solver%diffside%area_divider, ireals)
            Ay = atm%dx * atm%dz(ak,i,j) / real(solver%diffside%area_divider, ireals)

            emis = pi * B%planck(m) * (one-B%albedo(m))

            select case(idx(1))
            case(PPRTS_TOP_FACE)
              do idiff = 0, solver%difftop%dof-1
                if (.not.solver%difftop%is_inward(i1+idiff)) then ! eup on upper face
                  xsrc(idiff,k,i,j) = Az * emis / real(solver%difftop%streams, ireals)
                else                                              ! edn on upper face
                  xsrc(idiff,k,i,j) = 0
                endif
              enddo

            case(PPRTS_BOT_FACE)
              do idiff = 0, solver%difftop%dof-1
                if (solver%difftop%is_inward(i1+idiff)) then ! edn on lower face
                  xsrc(idiff,k+1,i,j) = Az * emis / real(solver%difftop%streams, ireals)
                else                                         ! eup on lower face
                  xsrc(idiff,k+1,i,j) = 0
                endif
              enddo

            case(PPRTS_LEFT_FACE)
              dof_offset = solver%difftop%dof
              do idiff = 0, solver%diffside%dof-1
                if (.not.solver%diffside%is_inward(i1+idiff)) then ! e_left on left face
                  xsrc(dof_offset+idiff,k,i,j) = Ax * emis / real(solver%diffside%streams, ireals)
                else
                  xsrc(dof_offset+idiff,k,i,j) = 0
                endif
              enddo

            case(PPRTS_RIGHT_FACE)
              dof_offset = solver%difftop%dof
              do idiff = 0, solver%diffside%dof-1
                if (solver%diffside%is_inward(i1+idiff)) then ! e_right on right face
                  xsrc(dof_offset+idiff,k,i+1,j) = Ax * emis / real(solver%diffside%streams, ireals)
                else
                  xsrc(dof_offset+idiff,k,i+1,j) = 0
                endif
              enddo

            case(PPRTS_REAR_FACE )
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idiff = 0, solver%diffside%dof-1
                if (.not.solver%diffside%is_inward(i1+idiff)) then ! e_backward
                  xsrc(dof_offset+idiff,k,i,j) = Ay * emis / real(solver%diffside%streams, ireals)
                else
                  xsrc(dof_offset+idiff,k,i,j) = 0
                endif
              enddo

            case(PPRTS_FRONT_FACE)
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idiff = 0, solver%diffside%dof-1
                if (solver%diffside%is_inward(i1+idiff)) then ! e_forward
                  xsrc(dof_offset+idiff,k,i,j+1) = Ay * emis / real(solver%diffside%streams, ireals)
                else
                  xsrc(dof_offset+idiff,k,i,j+1) = 0
                endif
              enddo

            case default
              call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1)+1))
            end select
          end associate
        enddo
      end associate

    end subroutine

  end subroutine setup_b


  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(solver, kabs, ksca, g, dz, ldir, coeff, &
      lone_dimensional, angles, lswitch_east, lswitch_north)
    class(t_solver), intent(in)       :: solver
    real(ireals),intent(in)           :: kabs, ksca, g, dz
    logical,intent(in)                :: ldir
    real(irealLUT),intent(out)        :: coeff(:)

    logical,intent(in)                :: lone_dimensional
    real(irealLUT),intent(in),optional:: angles(2)
    logical,intent(in),optional       :: lswitch_east, lswitch_north

    real(irealLUT) :: aspect_zx, tauz, w0
    integer(mpiint) :: ierr

    aspect_zx = real(dz / solver%atm%dx, irealLUT)
    w0        = real(ksca / max(kabs+ksca, epsilon(kabs)), irealLUT)
    tauz      = real((kabs+ksca) * dz, irealLUT)

    if(present(angles)) then
      aspect_zx = max(solver%OPP%dev%dirconfig%dims(3)%vrange(1), aspect_zx)
      tauz = max(solver%OPP%dev%dirconfig%dims(1)%vrange(1), &
           & min(solver%OPP%dev%dirconfig%dims(1)%vrange(2), tauz ))
      w0   = max(solver%OPP%dev%dirconfig%dims(2)%vrange(1), &
           & min(solver%OPP%dev%dirconfig%dims(2)%vrange(2), w0))
    else
      aspect_zx = max(solver%OPP%dev%diffconfig%dims(3)%vrange(1), aspect_zx)
      tauz = max(solver%OPP%dev%diffconfig%dims(1)%vrange(1), &
           & min(solver%OPP%dev%diffconfig%dims(1)%vrange(2), tauz ))
      w0   = max(solver%OPP%dev%diffconfig%dims(2)%vrange(1), &
           & min(solver%OPP%dev%diffconfig%dims(2)%vrange(2), w0))
    endif

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
    integer(mpiint) :: ierr

    associate(  atm     => solver%atm, &
                C_one   => solver%C_one, &
                C_diff  => solver%C_diff)

     !call CHKERR(1_mpiint, 'nca_wrapper not implemented')
     !print *, 'DEBUG: Stupid print statement to prevent unused compiler warnings', ediff, abso
     if(C_diff%dof.lt.i6) call CHKERR(1_mpiint, 'For NCA, need a solver with at least 6 diffuse streams to copy over some data')

     ! put additional values into a local ediff vec .. TODO: this is a rather dirty hack but is straightforward

     ! get ghost values for dz, planck, kabs and fluxes, ready to give it to NCA
     call DMGetGlobalVector(solver%C_diff%da ,gnca ,ierr) ; call CHKERR(ierr)

     call getVecPointer(solver%C_diff%da, gnca, xvgnca1d, xvgnca)
     xvgnca(  idz    , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%dz
     xvgnca(  iplanck, C_diff%zs:C_diff%ze  , C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%planck
     xvgnca(  ikabs  , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%kabs


     ! Copy Edn and Eup to local convenience vector
     call getVecPointer(C_diff%da, ediff, xv1d, xv)
     xvgnca( E_up,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_up,:,:,:)
     xvgnca( E_dn,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_dn,:,:,:)
     call restoreVecPointer(C_diff%da, ediff, xv1d, xv)

     call restoreVecPointer(C_diff%da, gnca, xvgnca1d, xvgnca )


     ! retrieve ghost values into l(ocal) nca vec
     call DMGetLocalVector (C_diff%da ,lnca ,ierr) ; call CHKERR(ierr)
     call VecSet(lnca, zero, ierr); call CHKERR(ierr)

     call DMGlobalToLocalBegin(C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)
     call DMGlobalToLocalEnd  (C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)

     call DMRestoreGlobalVector(C_diff%da, gnca, ierr); call CHKERR(ierr)

     ! call NCA
     call getVecPointer(C_diff%da, lnca, xvlnca1d, xvlnca)

     call ts_nca( atm%dx, atm%dy,                    &
       xvlnca(   idz        , : , : , :), &
       xvlnca(   iplanck    , : , : , :), &
       xvlnca(   ikabs      , : , : , :), &
       xvlnca(   E_dn       , : , : , :), &
       xvlnca(   E_up       , : , : , :), &
       xvlnca(   ihr        , : , : , :))


     ! return absorption
     call getVecPointer(C_one%da, abso, xhr1d, xhr)

     do k=C_one%zs,C_one%ze
       xhr(i0,k,:,:) = xvlnca( ihr , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) / &
         xvlnca( idz , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye)
     enddo
     call restoreVecPointer(C_one%da, abso, xhr1d, xhr )

     !return convenience vector that holds optical properties
     call restoreVecPointer(C_diff%da, lnca, xvlnca1d, xvlnca )
     call DMRestoreLocalVector(C_diff%da, lnca, ierr); call CHKERR(ierr)
   end associate
  end subroutine


  !> @brief build diffuse radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_diff_coeff(solver, A, C, opt_buildings)
    class(t_solver) :: solver
    type(tMat)      :: A
    type(t_coord)   :: C
    type(t_pprts_buildings), intent(in), optional :: opt_buildings

    integer(iintegers) :: i,j,k
    integer(mpiint) :: ierr

    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Setting coefficients for diffuse Light'

    call MatZeroEntries(A, ierr) ;call CHKERR(ierr)
    call mat_set_diagonal(A)

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

    if(present(opt_buildings)) then
      call MatAssemblyBegin(A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)
      call MatAssemblyEnd  (A, MAT_FLUSH_ASSEMBLY, ierr); call CHKERR(ierr)
      call set_buildings_coeff(solver, C, opt_buildings, A, ierr); call CHKERR(ierr)
    endif

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

      MatStencil :: row(4,0:C%dof-1), col(4,0:C%dof-1)
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
          if(norm.gt.one+10._ireals*epsilon(norm)) then ! could renormalize
            if(norm.gt.one+10._ireals*sqrt(epsilon(norm))) then ! fatally off
              print *,'diffuse sum(src==',src,') gt one',norm,'=>', v(src:C%dof**2:C%dof)
              print *,'get_coeff', solver%atm%kabs(atmk(solver%atm, k),i,j), &
                solver%atm%ksca(atmk(solver%atm, k),i,j), &
                solver%atm%g(atmk(solver%atm, k),i,j), &
                solver%atm%dz(atmk(solver%atm, k),i,j), &
                .False., solver%atm%l1d(atmk(solver%atm, k),i,j), '=>', v
              call CHKERR(1_mpiint, 'omg.. shouldnt be happening')
            endif ! fatal
          endif ! could renormalize
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
                !print *,i,j,k,'setting r ',toStr(src)//' (k='//toStr(col(MatStencil_i,src))//') ->', &
                !  toStr(dst)//' (k='//toStr(row(MatStencil_i,dst))//')'// &
                !  ':', i1 + dst*solver%difftop%dof + src, v(i1 + dst*solver%difftop%dof + src), &
                !  'invdof',src,dst,inv_dof(dst)
              else

                if(src.ne.dst) cycle ! in 1D has to be the same
                v(i1 + dst*solver%difftop%dof + src) = atm%a11(atmk(atm,k),i,j)
                !print *,i,j,k,'setting t ',toStr(src)//' (k='//toStr(col(MatStencil_i,src))//') ->', &
                !  toStr(dst)//' (k='//toStr(row(MatStencil_i,dst))//') :', i1 + dst*solver%difftop%dof + src, v(i1 + dst*solver%difftop%dof + src)
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
                  !print *,solver%myid, 'i '//toStr(i)//' j '//toStr(j), &
                  !  & ' Setting albedo for dst '//toStr(row)//' src '//toStr(col), &
                  !  & ':', -solver%atm%albedo(i,j) / real(solver%difftop%streams, ireals)
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

    !> @brief   apply blocking of diffuse radiation from buildings and do lambertian reflections
    !> @details Goal: first set all src dof on a buildings face towards all dst dof to zero
    !> \n       note, this assumes that the building is the full cell,
    !> \n       i.e. the albedo is applied on the outside of the cell
    subroutine set_buildings_coeff(solver, C, buildings, A, ierr)
    class(t_solver)                     :: solver
      type(t_coord),intent(in)            :: C
      type(t_pprts_buildings), intent(in) :: buildings
      type(tMat),intent(inout)            :: A
      integer(mpiint), intent(out)        :: ierr

      MatStencil         :: row(4,0:C%dof-1)  ,col(4,1)
      integer(iintegers) :: m, isrc, idst, dst, idx(4), dof_offset
      real(ireals) :: v(0:C%dof-1)

      ierr = 0

      associate( B => buildings )
        do m = 1, size(B%iface)
          v(:) = -zero

          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) -1 + [C%zs, C%xs, C%ys]

            associate(k => idx(2), i => idx(3), j => idx(4))
              dst = 0
              do idst = 0, solver%difftop%dof-1
                row(MatStencil_j,dst) = i
                row(MatStencil_k,dst) = j
                row(MatStencil_c,dst) = dst
                if(solver%diffside%is_inward(i1+idst)) then ! edn on bot face
                  row(MatStencil_i,dst) = k+1
                else ! eup on top face
                  row(MatStencil_i,dst) = k
                endif
                dst = dst + 1
              enddo

              do idst = 0, solver%diffside%dof-1
                row(MatStencil_i,dst) = k
                row(MatStencil_k,dst) = j
                row(MatStencil_c,dst) = dst
                if(solver%diffside%is_inward(i1+idst)) then ! e_right
                  row(MatStencil_j,dst) = i+1
                else ! e_left
                  row(MatStencil_j,dst) = i
                endif
                dst = dst + 1
              enddo

              do idst = 0, solver%diffside%dof-1
                row(MatStencil_i,dst) = k
                row(MatStencil_j,dst) = i
                row(MatStencil_c,dst) = dst
                if(solver%diffside%is_inward(i1+idst)) then ! e_forward
                  row(MatStencil_k,dst) = j+1
                else ! e_back
                  row(MatStencil_k,dst) = j
                endif
                dst = dst + 1
              enddo

              select case(idx(1))
              case(PPRTS_TOP_FACE)
                do idst = 0, solver%difftop%dof-1
                  if(.not.solver%difftop%is_inward(i1+idst)) v(idst) = -B%albedo(m) / real(solver%difftop%streams, ireals) ! eup
                enddo
                do isrc = 0, solver%difftop%dof-1
                  if(solver%difftop%is_inward(i1+isrc)) then ! edn
                    col(MatStencil_i,1) = k
                    col(MatStencil_j,1) = i
                    col(MatStencil_k,1) = j
                    col(MatStencil_c,1) = isrc
                    call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
                  endif
                enddo

              case(PPRTS_BOT_FACE)
                do idst = 0, solver%difftop%dof-1
                  if(solver%difftop%is_inward(i1+idst)) v(idst) = -B%albedo(m) / real(solver%difftop%streams, ireals) ! eup
                enddo
                do isrc = 0, solver%difftop%dof-1
                  if(.not.solver%difftop%is_inward(i1+isrc)) then ! eup
                    col(MatStencil_i,1) = k+1
                    col(MatStencil_j,1) = i
                    col(MatStencil_k,1) = j
                    col(MatStencil_c,1) = isrc
                    call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
                  endif
                enddo

              case(PPRTS_LEFT_FACE)
                dof_offset = solver%difftop%dof
                do idst = 0, solver%diffside%dof-1
                  if(.not.solver%diffside%is_inward(i1+idst)) &
                    & v(dof_offset+idst) = -B%albedo(m) / real(solver%diffside%streams, ireals)
                enddo
                do isrc = 0, solver%diffside%dof -1
                  if(solver%diffside%is_inward(i1+isrc)) then ! e_right
                    col(MatStencil_j,1) = i
                    col(MatStencil_k,1) = j
                    col(MatStencil_i,1) = k
                    col(MatStencil_c,1) = dof_offset + isrc
                    call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
                  endif
                enddo

              case(PPRTS_RIGHT_FACE)
                dof_offset = solver%difftop%dof
                do idst = 0, solver%diffside%dof-1
                  if(solver%diffside%is_inward(i1+idst)) &
                    & v(dof_offset+idst) = -B%albedo(m) / real(solver%diffside%streams, ireals)
                enddo
                do isrc = 0, solver%diffside%dof -1
                  if(.not.solver%diffside%is_inward(i1+isrc)) then ! e_left
                    col(MatStencil_j,1) = i+1
                    col(MatStencil_k,1) = j
                    col(MatStencil_i,1) = k
                    col(MatStencil_c,1) = dof_offset + isrc
                    call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
                  endif
                enddo

              case(PPRTS_REAR_FACE)
                dof_offset = solver%difftop%dof + solver%diffside%dof
                do idst = 0, solver%diffside%dof-1
                  if(.not.solver%diffside%is_inward(i1+idst)) &
                    & v(dof_offset+idst) = -B%albedo(m) / real(solver%diffside%streams, ireals)
                enddo
                do isrc = 0, solver%diffside%dof -1
                  if(solver%diffside%is_inward(i1+isrc)) then ! e_forward
                    col(MatStencil_j,1) = i
                    col(MatStencil_k,1) = j
                    col(MatStencil_i,1) = k
                    col(MatStencil_c,1) = dof_offset + isrc
                    call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
                  endif
                enddo

              case(PPRTS_FRONT_FACE)
                dof_offset = solver%difftop%dof + solver%diffside%dof
                do idst = 0, solver%diffside%dof-1
                  if(solver%diffside%is_inward(i1+idst)) &
                    & v(dof_offset+idst) = -B%albedo(m) / real(solver%diffside%streams, ireals)
                enddo
                do isrc = 0, solver%diffside%dof -1
                  if(.not.solver%diffside%is_inward(i1+isrc)) then ! e_backward
                    col(MatStencil_j,1) = i
                    col(MatStencil_k,1) = j+1
                    col(MatStencil_i,1) = k
                    col(MatStencil_c,1) = dof_offset + isrc
                    call MatSetValuesStencil(A, C%dof, row, i1, col, v, INSERT_VALUES, ierr); call CHKERR(ierr)
                  endif
                enddo

              case default
                call CHKERR(1_mpiint, 'wrong building fidx '//toStr(idx(1)))
              end select
            end associate

        enddo
      end associate
    end subroutine

  end subroutine

  subroutine pprts_get_result(solver, redn, reup, rabso, redir, opt_solution_uid, opt_buildings)
    class(t_solver) :: solver
    real(ireals),dimension(:,:,:),intent(inout),allocatable          :: redn,reup,rabso
    real(ireals),dimension(:,:,:),intent(inout),allocatable,optional :: redir
    integer(iintegers),optional,intent(in) :: opt_solution_uid
    type(t_pprts_buildings), optional, intent(inout) :: opt_buildings

    integer(iintegers)  :: uid, iside
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()
    integer(mpiint) :: ierr

    uid = get_arg(0_iintegers, opt_solution_uid)

    associate(solution => solver%solutions(uid))
      if(solution%lsolar_rad) then
        call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-pprts_show_edir", ierr); call CHKERR(ierr)
      endif
      call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-pprts_show_ediff", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, "-pprts_show_abso", ierr); call CHKERR(ierr)

      if(ldebug .and. solver%myid.eq.0) print *,'calling pprts_get_result',present(redir),'for uid',uid

      if(solution%lchanged) &
        & call CHKERR(1_mpiint, 'tried to get results from unrestored solution -- call restore_solution first')

      if(present(opt_solution_uid)) then
        if(.not.solution%lWm2_diff) &
          & call CHKERR(1_mpiint, 'solution vecs for diffuse radiation are not in W/m2 ... this is not what I expected')
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

        if( .not. solution%lsolar_rad ) then
          if(ldebug) then
            call CHKWARN(1_mpiint, 'Hey, You called pprts_get_result for uid '//toStr(uid)// &
              ' and provided an array for direct radiation.'// &
              ' However in this particular band we haven`t computed direct radiation.'// &
              ' I will return with edir=0 but are you sure this is what you intended?')
          endif
          redir = zero
        else
          if(.not.solution%lWm2_dir) &
            & call CHKERR(1_mpiint, 'tried to get result from a result vector(dir) which is not in [W/m2]')

          call getVecPointer(solver%C_dir%da, solution%edir, x1d, x4d)
          ! average of direct radiation of all fluxes through top faces
          redir = sum(x4d(0:solver%dirtop%dof-1, :, :, :), dim=1) / real(solver%dirtop%area_divider, ireals)
          call restoreVecPointer(solver%C_dir%da, solution%edir, x1d, x4d)

          if(ldebug) then
            if(solver%myid.eq.0) print *,'Edir vertically first column',redir(:, lbound(redir,2), lbound(redir,3))
            if(any(redir.lt.-one)) then
              print *,'Found direct radiation smaller than 0 in dir result... that should not happen',minval(redir)
              call CHKERR(1_mpiint)
            endif
          endif
        endif
      endif

      if(.not.solution%lWm2_diff) &
        & call CHKERR(1_mpiint, 'tried to get result from a result vector(diff) which is not in [W/m2]')


      redn = zero
      reup = zero
      call getVecPointer(solver%C_diff%da, solution%ediff, x1d, x4d)
      do iside=1,solver%difftop%dof
        if(solver%difftop%is_inward(iside)) then
          redn = redn + x4d(iside-1, :, :, :) / real(solver%difftop%area_divider, ireals)
        else
          reup = reup + x4d(iside-1, :, :, :) / real(solver%difftop%area_divider, ireals)
        endif
      enddo
      call restoreVecPointer(solver%C_diff%da, solution%ediff,x1d,x4d)

      if(solver%myid.eq.0 .and. ldebug .and. present(redir)) &
        print *,'mean surface Edir',meanval(redir(ubound(redir,1),:,:))
      if(solver%myid.eq.0 .and. ldebug) print *,'mean surface Edn',meanval(redn(ubound(redn,1), :,:))
      if(solver%myid.eq.0 .and. ldebug) print *,'mean surface Eup',meanval(reup(ubound(reup,1), :,:))

      if(ldebug .and. solution%lsolar_rad) then
        if(any(redn.lt.-one)) then
          print *,'Found radiation smaller than 0 in edn result... that should not happen',minval(redn)
          call exit(1)
        endif
        if(any(reup.lt.-one)) then
          print *,'Found radiation smaller than 0 in eup result... that should not happen',minval(reup)
          call exit(1)
        endif
      endif

      call getVecPointer(solver%C_one%da, solution%abso, x1d, x4d, readonly=.True.)
      rabso = x4d(i0,:,:,:)
      call restoreVecPointer(solver%C_one%da, solution%abso, x1d, x4d, readonly=.True.)

      if(present(opt_buildings)) then
        call fill_buildings()
      endif

      if(solver%myid.eq.0 .and. ldebug) print *,'get_result done'
    end associate

  contains

    subroutine fill_buildings
      integer(iintegers) :: m, idx(4), dof_offset, idof
      type(tVec) :: ledir, lediff

      associate(                               &
          & solution => solver%solutions(uid), &
          & B => opt_buildings,                &
          & C => solver%C_dir                  )

        if(solution%lsolar_rad) then
          if(.not.allocated(B%edir)) allocate(B%edir(size(B%iface)))

          call DMGetLocalVector(C%da, ledir, ierr); call CHKERR(ierr)
          call DMGlobalToLocalBegin(C%da, solution%edir, INSERT_VALUES, ledir, ierr) ;call CHKERR(ierr)
          call DMGlobalToLocalEnd  (C%da, solution%edir, INSERT_VALUES, ledir, ierr) ;call CHKERR(ierr)

          call getVecPointer(C%da, ledir, x1d, x4d, readonly=.True.)
          ! average of direct radiation of all fluxes through top faces
          !redir = sum(x4d(0:solver%dirtop%dof-1, :, :, :), dim=1) / real(solver%dirtop%area_divider, ireals)

          do m = 1, size(B%iface)
            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) -1 + [C%zs, C%xs, C%ys]

            associate(k => idx(2), i => idx(3), j => idx(4))

              select case(idx(1))
              case(PPRTS_TOP_FACE)
                dof_offset = 0
                B%edir(m) = sum(x4d(dof_offset:dof_offset+solver%dirtop%dof-1, k, i, j)) &
                  & / real(solver%dirtop%area_divider, ireals)
              case(PPRTS_BOT_FACE)
                dof_offset = 0
                B%edir(m) = sum(x4d(dof_offset:dof_offset+solver%dirtop%dof-1, k+1, i, j)) &
                  & / real(solver%dirtop%area_divider, ireals)
              case(PPRTS_LEFT_FACE)
                dof_offset = solver%dirtop%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset+solver%dirside%dof-1, k, i, j)) &
                  & / real(solver%dirside%area_divider, ireals)
              case(PPRTS_RIGHT_FACE)
                dof_offset = solver%dirtop%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset+solver%dirside%dof-1, k, i+1, j)) &
                  & / real(solver%dirside%area_divider, ireals)
              case(PPRTS_REAR_FACE )
                dof_offset = solver%dirtop%dof + solver%dirside%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset+solver%dirside%dof-1, k, i, j)) &
                  & / real(solver%dirside%area_divider, ireals)
              case(PPRTS_FRONT_FACE)
                dof_offset = solver%dirtop%dof + solver%dirside%dof
                B%edir(m) = sum(x4d(dof_offset:dof_offset+solver%dirside%dof-1, k, i, j+1)) &
                  & / real(solver%dirside%area_divider, ireals)
              case default
                call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1)+1))
              end select
            end associate
          enddo

          call restoreVecPointer(C%da, ledir, x1d, x4d, readonly=.True.)
          call DMRestoreLocalVector(C%da, ledir, ierr); call CHKERR(ierr)
        endif
      end associate


      associate(                               &
          & solution => solver%solutions(uid), &
          & B => opt_buildings,                &
          & C => solver%C_diff                  )
        if(.not.allocated(B%incoming)) allocate(B%incoming(size(B%iface)))
        if(.not.allocated(B%outgoing)) allocate(B%outgoing(size(B%iface)))

        call DMGetLocalVector(C%da, lediff, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(C%da, solution%ediff, INSERT_VALUES, lediff, ierr) ;call CHKERR(ierr)
        call DMGlobalToLocalEnd  (C%da, solution%ediff, INSERT_VALUES, lediff, ierr) ;call CHKERR(ierr)
        call getVecPointer(C%da, lediff, x1d, x4d, readonly=.True.)

        do m = 1, size(B%iface)
          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) -1 + [C%zs, C%xs, C%ys]

          associate(k => idx(2), i => idx(3), j => idx(4))

            B%incoming(m) = zero
            B%outgoing(m) = zero

            select case(idx(1))

            case(PPRTS_TOP_FACE)
              do idof = 0, solver%difftop%dof-1
                if(solver%difftop%is_inward(i1+idof)) then ! e_down
                  B%incoming(m) = B%incoming(m) + x4d(idof, k, i, j)
                else
                  B%outgoing(m) = B%outgoing(m) + x4d(idof, k, i, j)
                endif
              enddo
              B%incoming(m) = B%incoming(m) / real(solver%difftop%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%difftop%area_divider, ireals)

            case(PPRTS_BOT_FACE)
              do idof = 0, solver%difftop%dof-1
                if(solver%difftop%is_inward(i1+idof)) then ! e_down
                  B%outgoing(m) = B%outgoing(m) + x4d(idof, k+1, i, j)
                else
                  B%incoming(m) = B%incoming(m) + x4d(idof, k+1, i, j)
                endif
              enddo
              B%incoming(m) = B%incoming(m) / real(solver%difftop%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%difftop%area_divider, ireals)

            case(PPRTS_LEFT_FACE)
              dof_offset = solver%difftop%dof
              do idof = 0, solver%diffside%dof-1
                if(solver%diffside%is_inward(i1+idof)) then ! e_right
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i, j)
                else
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i, j)
                endif
              enddo
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case(PPRTS_RIGHT_FACE)
              dof_offset = solver%difftop%dof
              do idof = 0, solver%diffside%dof-1
                if(solver%diffside%is_inward(i1+idof)) then ! e_right
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i+1, j)
                else
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i+1, j)
                endif
              enddo
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case(PPRTS_REAR_FACE )
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idof = 0, solver%diffside%dof-1
                if(solver%diffside%is_inward(i1+idof)) then ! e_forward
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i, j)
                else
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i, j)
                endif
              enddo
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case(PPRTS_FRONT_FACE)
              dof_offset = solver%difftop%dof + solver%diffside%dof
              do idof = 0, solver%diffside%dof-1
                if(solver%diffside%is_inward(i1+idof)) then ! e_forward
                  B%outgoing(m) = B%outgoing(m) + x4d(dof_offset + idof, k, i, j+1)
                else
                  B%incoming(m) = B%incoming(m) + x4d(dof_offset + idof, k, i, j+1)
                endif
              enddo
              B%incoming(m) = B%incoming(m) / real(solver%diffside%area_divider, ireals)
              B%outgoing(m) = B%outgoing(m) / real(solver%diffside%area_divider, ireals)

            case default
              call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1)+1))
            end select
          end associate
        enddo

        call restoreVecPointer(C%da, solution%ediff, x1d, x4d, readonly=.True.)
        call DMRestoreLocalVector(C%da, lediff, ierr); call CHKERR(ierr)
      end associate
    end subroutine
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
    integer(mpiint) :: ierr

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
      call petscVecToF90(lvec_on_zero, C%da, outp, only_on_rank0=.True.)
      call VecDestroy(lvec_on_zero, ierr); call CHKERR(ierr)
    endif
  end subroutine

  subroutine gather_all_to_all(C, inp, outp)
    type(t_coord),intent(in) :: C
    real(ireals),intent(in), allocatable :: inp(:,:,:) ! local array from get_result
    real(ireals),intent(inout),allocatable :: outp(:,:,:) ! global sized array on rank 0

    type(tVec) :: vec, lvec

    integer(mpiint) :: myid, ierr
    call mpi_comm_rank(C%comm, myid, ierr)

    if(ldebug) then
      print *,myid,'exchange_var',allocated(inp), allocated(outp)
      print *,myid,'exchange_var shape',shape(inp)
    endif

    call DMGetGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)
    call f90VecToPetsc(inp, C%da, vec)
    call petscGlobalVecToAll(vec, C%da, lvec)
    call DMRestoreGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)

    call petscVecToF90(lvec, C%da, outp)
    call VecDestroy(lvec, ierr); call CHKERR(ierr)
  end subroutine
end module

