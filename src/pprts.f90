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

  use m_data_parameters, only : ireals, iintegers,               &
    init_mpi_data_parameters, mpiint,                            &
    zero, one, nil, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi, &
    default_str_len

  use m_helper_functions, only : CHKERR, deg2rad, rad2deg, norm, imp_allreduce_min, &
    imp_bcast, imp_allreduce_max, delta_scale, mpi_logical_and, mean, get_arg, approx, &
    inc

  use m_twostream, only: delta_eddington_twostream
  use m_schwarzschild, only: schwarzschild
  use m_optprop, only: t_optprop, t_optprop_1_2, t_optprop_3_6, t_optprop_8_10, t_optprop_3_10
  use m_eddington, only : eddington_coeff_zdun

  use m_tenstream_options, only : read_commandline_options, ltwostr, luse_eddington, twostr_ratio, &
    options_max_solution_err, options_max_solution_time, ltwostr_only, luse_twostr_guess,        &
    options_phi, lforce_phi, options_theta, lforce_theta, &
    lcalc_nca, lskip_thermal, lschwarzschild, ltopography

  use m_petsc_helpers, only : petscGlobalVecToZero, scatterZerotoPetscGlobal, &
    petscVecToF90, f90VecToPetsc, getVecPointer, restoreVecPointer

  implicit none
  private

  public :: t_solver, t_solver_1_2, t_solver_3_6, t_solver_8_10, t_solver_3_10, init_pprts, &
            set_optical_properties, set_global_optical_properties, &
            solve_pprts, set_angles, destroy_pprts, pprts_get_result, &
            pprts_get_result_toZero, t_coord

  integer(iintegers), parameter :: E_up=0, E_dn=1

  type t_coord
    integer(iintegers)      :: xs,xe                   ! local domain start and end indices
    integer(iintegers)      :: ys,ye                   ! local domain start and end indices
    integer(iintegers)      :: zs,ze                   ! local domain start and end indices
    integer(iintegers)      :: xm,ym,zm                ! size of local domain
    integer(iintegers)      :: gxs,gys,gzs             ! domain indices including ghost points
    integer(iintegers)      :: gxe,gye,gze             !
    integer(iintegers)      :: gxm,gym,gzm             ! size of local domain including ghosts
    integer(iintegers)      :: glob_xm,glob_ym,glob_zm ! global domain size
    integer(iintegers)      :: dof,dim                 ! degrees of freedom of Petsc Domain, dimension of dmda
    type(tDM)               :: da                      ! The Domain Decomposition Object
    PetscMPIInt,allocatable :: neighbors(:)            ! all 3d neighbours((x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
    integer(mpiint)         :: comm                    ! mpi communicatior for this DMDA
  end type

  type t_opticalprops
    real(ireals) :: kabs,ksca,g
  end type

  type t_atmosphere
    type(t_opticalprops) , allocatable , dimension(:,:,:) :: op
    real(ireals)         , allocatable , dimension(:,:,:) :: planck
    real(ireals)         , allocatable , dimension(:,:,:) :: a11, a12, a21, a22, a13, a23, a33
    real(ireals)         , allocatable , dimension(:,:,:) :: g1,g2
    real(ireals)         , allocatable , dimension(:,:,:) :: dz
    logical              , allocatable , dimension(:,:,:) :: l1d
    real(ireals)         , allocatable , dimension(:,:)   :: albedo
    real(ireals)                                          :: dx,dy
    integer(iintegers)                                    :: icollapse=1
    logical                                               :: lcollapse = .False.
  end type

  type t_sunangles
    real(ireals)        :: symmetry_phi
    integer(iintegers)  :: yinc,xinc
    real(ireals)        :: theta, phi, costheta, sintheta
  end type

  type t_suninfo
    type(t_sunangles),allocatable :: angles(:,:,:) ! defined on DMDA grid
    logical                       :: luse_topography=.False.
  end type

  KSP,save :: kspdir, kspdiff
  logical,save :: linit_kspdir=.False., linit_kspdiff=.False.

  type t_state_container
    integer(iintegers)  :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
    Vec                 :: edir,ediff,abso

    logical             :: lset        = .False. ! initialized?
    logical             :: lsolar_rad  = .False. ! direct radiation calculated?
    logical             :: lchanged    = .True.  ! did the flux change recently? -- call restore_solution to bring it in a coherent state

    ! save state of solution vectors... they are either in [W](true) or [W/m**2](false)
    logical             :: lintegrated_dir=.True. , lintegrated_diff=.True.

    !save error statistics
    real(ireals)        :: time   (30) = -one
    real(ireals)        :: maxnorm(30) = zero
    real(ireals)        :: twonorm(30) = zero
    real(ireals),allocatable :: ksp_residual_history(:)
  end type

  type t_dof
    integer(iintegers) :: dof
    logical, allocatable :: is_inward(:)
  end type

  type, abstract :: t_solver
    integer(mpiint)                 :: comm, myid, numnodes     ! mpi communicator, my rank and number of ranks in comm
    type(t_coord), allocatable      :: C_dir, C_diff, C_one, C_one1, C_one_atm, C_one_atm1
    type(t_atmosphere),allocatable  :: atm
    type(t_suninfo)                 :: sun
    type(tMat),allocatable          :: Mdir,Mdiff
    class(t_optprop), allocatable   :: OPP

    type(t_dof)                     :: difftop, diffside, dirtop, dirside

    logical                         :: lenable_solutions_err_estimates=.True.  ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
    type(tVec),allocatable          :: incSolar,b

    logical                         :: linitialized=.False.
    type(t_state_container)         :: solutions(-1000:1000)
  end type

  type, extends(t_solver) :: t_solver_1_2
  end type
  type, extends(t_solver) :: t_solver_8_10
  end type
  type, extends(t_solver) :: t_solver_3_6
  end type
  type, extends(t_solver) :: t_solver_3_10
  end type

  logical,parameter :: ldebug=.True.
  logical,parameter :: lcycle_dir=.True.
  logical,parameter :: lprealloc=.True.

  integer(iintegers),parameter :: minimal_dimension=3 ! this is the minimum number of gridpoints in x or y direction

  integer(mpiint) :: ierr

  contains

  !> @brief Main routine to setup TenStream solver
  !> @details This will setup the PETSc DMDA grid and set other grid information, needed for the TenStream
  !> \n Nx, Ny Nz are either global domain size or have to be local sizes if present(nxproc,nyproc)
  !> \n where nxproc and nyproc then are the number of pixel per rank for all ranks -- i.e. sum(nxproc) != Nx_global
  subroutine init_pprts(icomm, Nz,Nx,Ny, dx,dy, phi0, theta0, solver, dz1d, dz3d, nxproc, nyproc, collapseindex)
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

    integer(iintegers) :: k,i,j
    !    character(default_str_len),parameter :: tenstreamrc='./.tenstreamrc'

    if(.not.solver%linitialized) then

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

        class is (t_solver_8_10)

          allocate(solver%difftop%is_inward(2))
          solver%difftop%is_inward = [.False.,.True.]

          allocate(solver%diffside%is_inward(4))
          solver%diffside%is_inward = [.False.,.True.,.False.,.True.]

          allocate(solver%dirtop%is_inward(4))
          solver%dirtop%is_inward = .True.

          allocate(solver%dirside%is_inward(2))
          solver%dirside%is_inward = .True.

        class is (t_solver_3_10)

          allocate(solver%difftop%is_inward(2))
          solver%difftop%is_inward = [.False.,.True.]

          allocate(solver%diffside%is_inward(4))
          solver%diffside%is_inward = [.False.,.True.,.False.,.True.]

          allocate(solver%dirtop%is_inward(1))
          solver%dirtop%is_inward = [.True.]

          allocate(solver%dirside%is_inward(1))
          solver%dirside%is_inward = [.True.]
        class default
        stop 'init pprts: unexpected type for solver'
      end select

      solver%difftop%dof= size(solver%difftop%is_inward)
      solver%diffside%dof= size(solver%diffside%is_inward)
      solver%dirtop%dof= size(solver%dirtop%is_inward)
      solver%dirside%dof= size(solver%dirside%is_inward)

      call init_mpi_data_parameters(icomm)

      solver%comm = icomm
      call mpi_comm_rank(solver%comm, solver%myid, ierr)     ; call CHKERR(ierr)
      call mpi_comm_size(solver%comm, solver%numnodes, ierr) ; call CHKERR(ierr)

      if(ldebug.and.solver%myid.eq.0) then
        print *,'Solver dirtop:', solver%dirtop%is_inward, ':', solver%dirtop%dof
        print *,'Solver dirside:', solver%dirside%is_inward, ':', solver%dirside%dof
        print *,'Solver difftop:', solver%difftop%is_inward, ':', solver%difftop%dof
        print *,'Solver diffside:', solver%diffside%is_inward, ':', solver%diffside%dof
      endif

      !      call PetscInitialize(tenstreamrc ,ierr) ;call CHKERR(ierr)
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr) ;call CHKERR(ierr)
#ifdef _XLF
      call PetscPopSignalHandler(ierr); call CHKERR(ierr) ! in case of xlf ibm compilers, remove petsc signal handler -- otherwise we dont get fancy signal traps from boundschecking or FPE's
#endif

      call read_commandline_options()

      if(present(nxproc) .and. present(nyproc) ) then
        if(ldebug.and.solver%myid.eq.0) print *,'nxproc',shape(nxproc),'::',nxproc
        if(ldebug.and.solver%myid.eq.0) print *,'nyproc',shape(nyproc),'::',nyproc
        if(present(collapseindex)) then
          call setup_grid(solver, Nz, Nx, Ny, nxproc,nyproc, collapseindex=collapseindex)
        else
          call setup_grid(solver, Nz, Nx, Ny, nxproc,nyproc)
        endif
      else
        if(present(collapseindex)) then
          call setup_grid(solver, Nz, max(minimal_dimension, Nx), max(minimal_dimension, Ny), collapseindex=collapseindex)
        else
          call setup_grid(solver, Nz, max(minimal_dimension, Nx), max(minimal_dimension, Ny) )
        endif

        solver%linitialized = .True.
      endif

      call setup_atm()

      ! init work vectors
      call init_memory(solver%C_dir, solver%C_diff, solver%incSolar, solver%b)

      ! init petsc logging facilities
      !call setup_logging()

    else
      print *,solver%myid,'You tried to initialize already initialized Tenstream          &
        &solver. This should not be done. If you need to reinitialize the grids, &
        &call destroy_tenstream() first.'
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

      if(.not.allocated(solver%atm%dz) ) allocate(solver%atm%dz( solver%C_one_atm%zs:solver%C_one_atm%ze, solver%C_one_atm%xs:solver%C_one_atm%xe, solver%C_one_atm%ys:solver%C_one_atm%ye ))

      if(present(dz1d)) then
        do j=solver%C_one_atm%ys,solver%C_one_atm%ye
        do i=solver%C_one_atm%xs,solver%C_one_atm%xe
        solver%atm%dz(:,i,j) = dz1d
        enddo
        enddo
      endif
      if(present(dz3d)) then
        if( any( shape(dz3d).ne.shape(solver%atm%dz) ) ) then
          print *,'Whoops I got a 3D dz definition but this does not correspond to the grid definition :: shapes: ', shape(dz3d), ' vs. ',shape(solver%atm%dz)
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
        allocate(solver%atm%l1d( solver%C_one_atm%zs:solver%C_one_atm%ze, solver%C_one_atm%xs:solver%C_one_atm%xe, solver%C_one_atm%ys:solver%C_one_atm%ye ) )
      endif

      !TODO if we have a horiz. staggered grid, this may lead to the point where one 3d box has a outgoing sideward flux but the adjacent
      !1d box does not send anything back --> therefore huge absorption :( -- would need to introduce mirror boundary conditions for
      !sideward fluxes in 1d boxes
      do j=solver%C_one_atm%ys,solver%C_one_atm%ye
        do i=solver%C_one_atm%xs,solver%C_one_atm%xe
          solver%atm%l1d(solver%C_one_atm%ze,i,j) = twostr_ratio*solver%atm%dz(solver%C_one_atm%ze,i,j).gt.solver%atm%dx
          do k=solver%C_one_atm%ze-1,solver%C_one_atm%zs,-1
            solver%atm%l1d(k,i,j) = twostr_ratio*solver%atm%dz(k,i,j).gt.solver%atm%dx
          enddo
        enddo
      enddo
      if(ltwostr_only) solver%atm%l1d = .True.

      if(present(collapseindex)) then
        solver%atm%lcollapse=.True.
        solver%atm%icollapse=collapseindex
        solver%atm%l1d(solver%C_one_atm%zs:atmk(solver%atm, solver%C_one%zs),:,:) = .True. ! if need to be collapsed, they have to be 1D.
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

      if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Setting up the DMDA grid for ',Nz,Nx,Ny,'using ',solver%numnodes,' nodes'

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

        allocate(C)

        C%comm = icomm

        C%dof = i1*dof
        if(present(nxproc) .and. present(nyproc) ) then
          call DMDACreate3d( C%comm ,                                                 &
            bn                      , boundary                 , boundary           , &
            DMDA_STENCIL_STAR       ,                                                 &
            Nz                      , i1*sum(nxproc)           , i1*sum(nyproc)     , &
            i1                      , i1*size(nxproc)          , i1*size(nyproc)    , &
            C%dof                   , stencil_size             ,                      &
            Nz                      , nxproc                   , nyproc             , &
            C%da                    , ierr)
          call CHKERR(ierr)
        else
          call DMDACreate3d( C%comm ,                                                 &
            bn                      , boundary                 , boundary           , &
            DMDA_STENCIL_STAR       ,                                                 &
            i1*Nz                   , Nx                       , Ny                 , &
            i1                      , PETSC_DECIDE             , PETSC_DECIDE       , &
            C%dof                   , stencil_size             ,                      &
            Nz                      , PETSC_NULL_INTEGER       , PETSC_NULL_INTEGER , &
            C%da                    , ierr) ;call CHKERR(ierr)
        endif

        call DMSetup(C%da,ierr)                                              ; call CHKERR(ierr)
        call DMSetMatType(C%da, MATAIJ, ierr)                                ; call CHKERR(ierr)
        if(lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ; call CHKERR(ierr)

        call DMSetFromOptions(C%da, ierr)                                    ; call CHKERR(ierr)
        if(ldebug) call DMView(C%da, PETSC_VIEWER_STDOUT_WORLD ,ierr)        ; call CHKERR(ierr)
        call setup_coords(C)
      end subroutine
      subroutine setup_coords(C)
        type(t_coord) :: C
        DMBoundaryType :: bx, by, bz
        DMDAStencilType :: st
        integer(iintegers) :: stencil_width, nproc_x, nproc_y, nproc_z, Ndof

        call DMDAGetInfo(C%da, C%dim,                             &
          C%glob_zm, C%glob_xm, C%glob_ym,                        &
          nproc_x, nproc_y, nproc_z, Ndof, stencil_width,         &
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
        !if(numnodes.gt.i1) then
        !  if(ldebug.and.(C%dim.eq.3)) print *,'PETSC id',solver%myid,C%dim,'Neighbors are',C%neighbors([10,4,16,22]),'while I am ',C%neighbors(13)
        !  if(ldebug.and.(C%dim.eq.2)) print *,'PETSC id',solver%myid,C%dim,'Neighbors are',C%neighbors([1,3,7,5]),'while I am ',C%neighbors(4)
        !endif
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
      class(t_solver), intent(inout)    :: solver
      real(ireals),intent(in)          :: phi0           !< @param[in] phi0   solar azmiuth and zenith angle
      real(ireals),intent(in)          :: theta0         !< @param[in] theta0 solar azmiuth and zenith angle
      real(ireals),optional,intent(in) :: phi2d  (:,:)   !< @param[in] phi2d   if given, horizontally varying azimuth
      real(ireals),optional,intent(in) :: theta2d(:,:)   !< @param[in] theta2d if given, and zenith angle

      logical :: lchanged_theta, lchanged_phi

      if(.not.solver%linitialized) then
          print *,solver%myid,'You tried to set angles in the Tenstream solver.  &
              & This should be called right after init_tenstream'
          ierr=1; call CHKERR(ierr)
      endif

      if(allocated(solver%sun%angles)) then ! was initialized
        if(present(theta2d)) then
          lchanged_theta = any(.not.approx(theta2d, solver%sun%angles(1,:,:)%theta))
        else
          lchanged_theta = any(.not.approx(theta0, solver%sun%angles(1,:,:)%theta))
        endif
        if(present(phi2d)) then
          lchanged_phi = any(.not.approx(phi2d, solver%sun%angles(1,:,:)%phi))
        else
          lchanged_phi = any(.not.approx(phi0, solver%sun%angles(1,:,:)%phi))
        endif
        if(solver%myid.eq.0 .and. ldebug) print *,'tenstr set_angles -- changed angles?',lchanged_theta, lchanged_phi
        if(.not. lchanged_theta .and. .not. lchanged_phi) then
          return
        endif
      endif

      if ( present(phi2d) .and. present(theta2d) ) then
          call setup_suninfo(solver, phi0, theta0, solver%sun, solver%C_one, phi2d=phi2d, theta2d=theta2d)
      elseif ( present(phi2d) ) then
          call setup_suninfo(solver, phi0, theta0, solver%sun, solver%C_one, phi2d=phi2d)
      elseif ( present(theta2d) ) then
          call setup_suninfo(solver, phi0, theta0, solver%sun, solver%C_one, theta2d=theta2d)
      else
          call setup_suninfo(solver, phi0, theta0, solver%sun, solver%C_one)
      endif

      if(ltwostr_only) return ! dont need anything here, we just compute Twostream anyway

      ! init box montecarlo model
      select type(solver)
        class is (t_solver_1_2)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_1_2::solver%OPP)

        class is (t_solver_3_6)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_3_6::solver%OPP)

        class is (t_solver_8_10)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_8_10::solver%OPP)

        class is (t_solver_3_10)
           if(.not.allocated(solver%OPP) ) allocate(t_optprop_3_10::solver%OPP)

        class default
        stop 'init pprts: unexpected type for solver'
      end select

      call solver%OPP%init(pack(solver%sun%angles%symmetry_phi,.True.), &
        pack(solver%sun%angles%theta,.True.),solver%comm )

  !     allocate(solver%OPP)
  !     call solver%OPP%init(pack(solver%sun%angles%symmetry_phi,.True.), &
  !                          pack(solver%sun%angles%theta,.True.),solver%comm, &
  !                          Ndiff=sum(solver%dofdiff), Ndir=sum(solver%dofdir))

  !     if(any(solver%atm%l1d.eqv..False.)) call OPP_3_6%init(pack(solver%sun%angles%symmetry_phi,.True.),pack(solver%sun%angles%theta,.True.),solver%comm)
  !  !   if(any(solver%atm%l1d.eqv..False.)) call OPP_8_10%init(pack(solver%sun%angles%symmetry_phi,.True.),pack(solver%sun%angles%theta,.True.),solver%comm)
  !     if(.not.luse_eddington)      call OPP_1_2%init (pack(solver%sun%angles%symmetry_phi,.True.),pack(solver%sun%angles%theta,.True.),solver%comm)

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

    integer(iintegers) :: k

    if(.not.allocated(sun%angles)) &
        allocate(sun%angles(C_one%zs:C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye))

    if(lforce_phi) then
        sun%angles(:,:,:)%phi = options_phi
    else
        if(present(phi2d)) then
            do k = C_one%zs, C_one%ze
                sun%angles(k,:,:)%phi = phi2d
            enddo
        else
            sun%angles(:,:,:)%phi   = phi0
        endif
    endif

    if(lforce_theta) then
      sun%angles(:,:,:)%theta = options_theta
    else
        if(present(theta2d)) then
            do k = C_one%zs, C_one%ze
                sun%angles(k,:,:)%theta = theta2d
            enddo
        else
            sun%angles(:,:,:)%theta = theta0
        endif
    endif

    if(sun%luse_topography) then
      call setup_topography(solver%atm, solver%C_one, solver%C_one1, sun)
    endif

    sun%angles(:,:,:)%costheta = max( cos(deg2rad(sun%angles(:,:,:)%theta)), zero)
    sun%angles(:,:,:)%sintheta = max( sin(deg2rad(sun%angles(:,:,:)%theta)), zero)

    ! use symmetry for direct beam: always use azimuth [0,90] an just reverse the order where we insert the coeffs
    sun%angles(:,:,:)%symmetry_phi = sym_rot_phi(sun%angles(:,:,:)%phi)

    where(sin(deg2rad(sun%angles%phi)).gt.zero ) ! phi between 0 and 180 degreee
        sun%angles%xinc=i0
    else where
        sun%angles%xinc=i1
    end where

    where(cos(deg2rad(sun%angles%phi)).lt.zero) ! phi between 90 and 270 degree
        sun%angles%yinc=i1
    else where
        sun%angles%yinc=i0
    end where

    print *,solver%myid,'setup_dir_inc done'

    if(ldebug) print *,solver%myid,'setup_dir_inc done', &
      count(sun%angles%xinc.eq.0),count(sun%angles%xinc.eq.i1), &
      count(sun%angles%yinc.eq.0),count(sun%angles%yinc.eq.i1), &
      '::', minval(sun%angles%xinc), maxval(sun%angles%xinc), minval(sun%angles%yinc), maxval(sun%angles%yinc)

    !print *, 'phi', sun%angles(1,1,1)%phi, 'sym_phi', sun%angles(1,1,1)%symmetry_phi
    contains
        elemental function sym_rot_phi(phi)
            real(ireals) :: sym_rot_phi
            real(ireals),intent(in) :: phi
            ! ''swap'' phi axis down to the range of [0,180]
            sym_rot_phi = acos(cos(deg2rad(phi)))
            !print *,'1st phi swap',phi,' :: ',sym_rot_phi,'=',phi*pi/180,cos(phi*pi/180),acos(cos(phi*pi/180))
            ! and then mirror it onto range [0,90]
            sym_rot_phi = int( rad2deg( asin(sin(sym_rot_phi)) ) )
            !print *,'2nd phi swap',phi,' :: ',sym_rot_phi,'=',sin(sym_rot_phi),asin(sin(sym_rot_phi)),asin(sin(sym_rot_phi)) /pi * 180,int(asin(sin(sym_rot_phi)) /pi * 180)
        end function
  end subroutine

  !> @brief setup topography information
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information and
  !>   \n tweak the local sun angles to bend the rays.
  subroutine setup_topography(atm, C_one, C_one1, sun)
    type(t_atmosphere), intent(in) :: atm
    type(t_coord), intent(in) :: C_one, C_one1
    type(t_suninfo),intent(inout) :: sun

    type(tVec) :: vgrad_x, vgrad_y
    real(ireals),Pointer :: grad_x(:,:,:,:)=>null(), grad_x1d(:)=>null()
    real(ireals),Pointer :: grad_y(:,:,:,:)=>null(), grad_y1d(:)=>null()

    integer(iintegers) :: i,j,k
    real(ireals) :: newtheta, newphi, xsun(3)
    real(ireals) :: rotmat(3, 3), newxsun(3)

    if(.not.allocated(sun%angles)) stop 'You called  setup_topography() &
        &but the sun struct is not yet up, make sure setup_suninfo is called before'
    if(.not.allocated(atm%dz)) stop 'You called  setup_topography() &
        &but the atm struct is not yet up, make sure we have atm%dz before'

    call compute_gradient(atm, C_one1, vgrad_x, vgrad_y)

    call getVecPointer(vgrad_x , C_one1%da, grad_x1d, grad_x)
    call getVecPointer(vgrad_y , C_one1%da, grad_y1d, grad_y)


    rotmat = reshape((/ one , zero, zero,  &
                        zero, one , zero,  &
                        nil , nil , one/), &
                     (/3, 3/), order=(/2, 1/) )

    do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe

            ! if we are at the global boundary we have to take that the gradient does not get too
            ! steep. That wouldnt make sense for cyclic boundaries... ! todo: check what we should do?
            !if(i.eq.i0 .or. i.eq.C_one%glob_xm-i1 .or. j.eq.i0 .or. j.eq.C_one%glob_ym-i1) then
            !   !print *,solver%myid,'global edge at:',i,j
            !    cycle
            !endif

            do k=C_one%zs,C_one%ze


                ! Vector of sun direction
                xsun(1) = sin(deg2rad(sun%angles(k,i,j)%theta))*sin(deg2rad(sun%angles(k,i,j)%phi))
                xsun(2) = sin(deg2rad(sun%angles(k,i,j)%theta))*cos(deg2rad(sun%angles(k,i,j)%phi))
                xsun(3) = cos(deg2rad(sun%angles(k,i,j)%theta))

                xsun = xsun / norm(xsun)

                rotmat(3, 1) = grad_x(i0, k, i, j)
                rotmat(3, 2) = grad_y(i0, k, i, j)

                newxsun = matmul(rotmat, xsun)

                newtheta = rad2deg(atan2(sqrt(newxsun(1)**2 + newxsun(2)**2), newxsun(3)))

                !newphi in meteorologiecal definitions: clockwise from y-axis
                newphi = rad2deg(atan2(newxsun(1), newxsun(2)))

                ! if(i.eq.C_one1%xs) print *,solver%myid,i,j,k, '::',hhl(i0,k+1,i,j-1:j+1),'::', grad ,'::',sun%angles(k,i,j)%theta, newtheta, '::', sun%angles(k,i,j)%phi, newphi

                sun%angles(k,i,j)%theta = max(zero, min( 90._ireals, newtheta ))
                sun%angles(k,i,j)%phi = newphi
            enddo
        enddo
    enddo

    call restoreVecPointer(vgrad_x, grad_x1d, grad_x)
    call restoreVecPointer(vgrad_y, grad_y1d, grad_y)

    call DMRestoreLocalVector(C_one1%da, vgrad_x, ierr);  call CHKERR(ierr)
    call DMRestoreLocalVector(C_one1%da, vgrad_y, ierr);  call CHKERR(ierr)
  end subroutine

  !> @brief compute gradient from dz3d
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information
  subroutine compute_gradient(atm, C_one1, vgrad_x, vgrad_y)
    type(t_atmosphere),intent(in) :: atm
    type(t_coord), intent(in) :: C_one1
    type(tVec) :: vgrad_x, vgrad_y

    type(tVec) :: vhhl
    real(ireals),Pointer :: hhl(:,:,:,:)=>null(), hhl1d(:)=>null()
    real(ireals),Pointer :: grad_x(:,:,:,:)=>null(), grad_x1d(:)=>null()
    real(ireals),Pointer :: grad_y(:,:,:,:)=>null(), grad_y1d(:)=>null()

    integer(iintegers) :: i,j,k

    real(ireals) :: zm(4), maxheight, global_maxheight

    if(.not.allocated(atm%dz)) stop 'You called  compute_gradient()&
      &but the atm struct is not yet up, make sure we have atm%dz before'

    call DMGetLocalVector(C_one1%da, vhhl, ierr) ;call CHKERR(ierr)
    call getVecPointer(vhhl, C_one1%da, hhl1d, hhl)

    hhl(i0, C_one1%ze, :, :) = zero
    do j=C_one1%ys,C_one1%ye
      do i=C_one1%xs,C_one1%xe
        hhl(i0, C_one1%ze, i, j ) = zero

        do k=C_one1%ze-1,C_one1%zs,-1
          hhl(i0,k,i,j) = hhl(i0,k+1,i,j)+atm%dz(atmk(atm, k),i,j)
        enddo
      enddo
    enddo


    maxheight = maxval(hhl(i0,C_one1%zs,C_one1%xs:C_one1%xe,C_one1%ys:C_one1%ye))
    call imp_allreduce_max(C_one1%comm, maxheight, global_maxheight)

    do j=C_one1%ys,C_one1%ye
      do i=C_one1%xs,C_one1%xe
        hhl(i0, :, i, j) = hhl(i0, :, i, j) + global_maxheight - hhl(i0, C_one1%zs, i, j)
      enddo
    enddo

    call restoreVecPointer(vhhl, hhl1d, hhl)

    call DMLocalToLocalBegin(C_one1%da, vhhl, INSERT_VALUES, vhhl,ierr) ;call CHKERR(ierr)
    call DMLocalToLocalEnd(C_one1%da, vhhl, INSERT_VALUES, vhhl,ierr) ;call CHKERR(ierr)

    call getVecPointer(vhhl , C_one1%da, hhl1d, hhl)

    call DMGetLocalVector(C_one1%da, vgrad_x, ierr) ;call CHKERR(ierr)
    call DMGetLocalVector(C_one1%da, vgrad_y, ierr) ;call CHKERR(ierr)

    call getVecPointer(vgrad_x , C_one1%da, grad_x1d, grad_x)
    call getVecPointer(vgrad_y , C_one1%da, grad_y1d, grad_y)

    do j=C_one1%ys,C_one1%ye
      do i=C_one1%xs,C_one1%xe
        do k=C_one1%zs,C_one1%ze-1
          ! Mean heights of adjacent columns
          zm(1) = (hhl(i0,k,i-1,j) + hhl(i0,k+1,i-1,j)) / 2
          zm(2) = (hhl(i0,k,i+1,j) + hhl(i0,k+1,i+1,j)) / 2

          zm(3) = (hhl(i0,k,i,j-1) + hhl(i0,k+1,i,j-1)) / 2
          zm(4) = (hhl(i0,k,i,j+1) + hhl(i0,k+1,i,j+1)) / 2

          ! Gradient of height field
          grad_x(i0, k, i, j) = -(zm(2)-zm(1)) / (2._ireals*atm%dx)
          grad_y(i0, k, i, j) = -(zm(4)-zm(3)) / (2._ireals*atm%dy)
        enddo
        zm(1) = hhl(i0, C_one1%ze, i-1, j)
        zm(2) = hhl(i0, C_one1%ze, i+1, j)
        zm(3) = hhl(i0, C_one1%ze, i, j-1)
        zm(4) = hhl(i0, C_one1%ze, i, j+1)
        ! Gradient of height field
        grad_x(i0, C_one1%ze, i, j) = -(zm(2)-zm(1)) / (2._ireals*atm%dx)
        grad_y(i0, C_one1%ze, i, j) = -(zm(4)-zm(3)) / (2._ireals*atm%dy)
      enddo
    enddo

    call restoreVecPointer(vhhl, hhl1d, hhl)
    call DMRestoreLocalVector(C_one1%da ,vhhl ,ierr);  call CHKERR(ierr)

    call restoreVecPointer(vgrad_x, grad_x1d, grad_x)
    call restoreVecPointer(vgrad_y, grad_y1d, grad_y)
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
        stop 'init_Matrix::setPreallocation : poor mans preallocation should really not be used...'
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

    call MatSetFromOptions(A,ierr) ;call CHKERR(ierr)
    call MatSetUp(A,ierr) ;call CHKERR(ierr)

    call mat_set_diagonal(A,C)
  contains
    subroutine mat_set_diagonal(A,C,vdiag)
      type(tMat) :: A
      type(t_coord),intent(in) :: C
      real(ireals),intent(in),optional :: vdiag

      integer(iintegers) :: i,j,k,dof
      MatStencil :: row(4,1), col(4,1)
      real(ireals) :: v(1)

      !TODO -- we should use this form... however this does somehow corrupt preallocation? - maybe fix this
      !        Vec :: diag
      !        call DMCreateGlobalVector(C%da,diag,ierr) ;call CHKERR(ierr)
      !        call VecSet(diag,one,ierr) ;call CHKERR(ierr)
      !        call MatDiagonalSet( A, diag, INSERT_VALUES,ierr) ;call CHKERR(ierr)
      !        call VecDestroy(diag,ierr) ;call CHKERR(ierr)

      if(present(vdiag)) then
        v(1) = vdiag
      else
        v(1)=one
      endif

      do j=C%ys,C%ye
        row(MatStencil_k,1) = j
        col(MatStencil_k,1) = j

        do i=C%xs,C%xe
          row(MatStencil_j,1) = i
          col(MatStencil_j,1) = i

          do k=C%zs,C%ze
            row(MatStencil_i,1) = k
            col(MatStencil_i,1) = k

            do dof=0,C%dof-1
              row(MatStencil_c,1) = dof
              col(MatStencil_c,1) = dof

              call MatSetValuesStencil(A,i1, row,i1, col , v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)
            enddo
          enddo
        enddo
      enddo

      call MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY,ierr) ;call CHKERR(ierr)
      call MatAssemblyEnd  (A,MAT_FLUSH_ASSEMBLY,ierr) ;call CHKERR(ierr)
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
            xinc = solver%sun%angles(k,i,j)%xinc
            yinc = solver%sun%angles(k,i,j)%yinc

            do idst = 1, solver%dirtop%dof
              dst = idst
              row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
            enddo

            do idst = 1, solver%dirside%dof
              dst = idst + solver%dirtop%dof
              row(MatStencil_j,dst) = i+xinc   ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
            enddo

            do idst = 1, solver%dirside%dof
              dst = idst + solver%dirtop%dof + solver%dirside%dof
              row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j+yinc   ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
            enddo

            do isrc = 1, solver%dirtop%dof
              src = isrc
              col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k; col(MatStencil_c,src) = src-i1 ! Define transmission towards the lower/upper lid
            enddo

            do isrc = 1, solver%dirside%dof
              src = isrc + solver%dirtop%dof
              col(MatStencil_j,src) = i+1-xinc   ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Define transmission towards the left/right lid
            enddo

            do isrc = 1, solver%dirside%dof
              src = isrc + solver%dirtop%dof + solver%dirside%dof
              col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j+1-yinc   ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Define transmission towards the front/back lid
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

    if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
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
            call inc( xd(E_dn, k+1, i,j), 2*one )
            call inc( xd(E_up, k  , i,j), 2*one )
          else

            ! Diffuse Coefficients Code Begin
            do idof=1, solver%difftop%dof
              src = idof-1
              if (solver%difftop%is_inward(idof)) then
                col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
              else
                col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k+1   ; col(MatStencil_c,src) = src
              endif
            enddo

            do idof=1, solver%diffside%dof
              src = solver%difftop%dof + idof -1
              if (solver%diffside%is_inward(idof)) then
                col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
              else
                col(MatStencil_j,src) = i+1  ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
              endif
            enddo

            do idof=1, solver%diffside%dof
              src = solver%difftop%dof + solver%diffside%dof + idof -1
              if (solver%diffside%is_inward(idof)) then
                col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
              else
                col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j+1   ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
              endif
            enddo

            do idof=1, solver%difftop%dof
              dst = idof-1
              if (solver%difftop%is_inward(idof)) then
                row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k+1   ; row(MatStencil_c,dst) = dst
              else
                row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
              endif
            enddo

            do idof=1, solver%diffside%dof
              dst = solver%difftop%dof + idof-1
              if (solver%diffside%is_inward(idof)) then
                row(MatStencil_j,dst) = i+1  ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
              else
                row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
              endif
            enddo

            do idof=1, solver%diffside%dof
              dst = solver%difftop%dof + solver%diffside%dof + idof-1
              if (solver%diffside%is_inward(idof)) then
                row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j+1   ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
              else
                row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
              endif
            enddo
            ! Diffuse Coefficients Code End

            do idst = 0,C%dof-1
              icnt = icnt+1
              dst_id = myid
              if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 .and. row(MatStencil_j,idst).lt.C%xs ) dst_id = C%neighbors(10) ! have real neighbor west  and is not local entry
              if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 .and. row(MatStencil_j,idst).gt.C%xe ) dst_id = C%neighbors(16) ! have real neighbor east  and is not local entry
              if( C%neighbors( 4).ne.myid .and. C%neighbors( 4).ge.i0 .and. row(MatStencil_k,idst).lt.C%ys ) dst_id = C%neighbors( 4) ! have real neighbor south and is not local entry
              if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 .and. row(MatStencil_k,idst).gt.C%ye ) dst_id = C%neighbors(22) ! have real neighbor north and is not local entry

              do isrc = 0,C%dof-1
                src_id = myid

                if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 .and. col(MatStencil_j,isrc).lt.C%xs ) src_id = C%neighbors(10) ! have real neighbor west  and is not local entry
                if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 .and. col(MatStencil_j,isrc).gt.C%xe ) src_id = C%neighbors(16) ! have real neighbor east  and is not local entry
                if( C%neighbors( 4).ne.myid .and. C%neighbors( 4).ge.i0 .and. col(MatStencil_k,isrc).lt.C%ys ) src_id = C%neighbors( 4) ! have real neighbor south and is not local entry
                if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 .and. col(MatStencil_k,isrc).gt.C%ye ) src_id = C%neighbors(22) ! have real neighbor north and is not local entry

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
            call inc( xd(src, C%ze, i,j), one )
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

    if(myid.eq.0 .and. ldebug) print *,myid,'diffuse d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)

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
  pure function atmk(atm, k) ! return vertical index value for DMDA grid on atmosphere grid
    integer(iintegers) :: atmk
    integer(iintegers),intent(in) :: k
    type(t_atmosphere),intent(in) :: atm
    atmk = k+atm%icollapse-1
  end function


  subroutine set_optical_properties(solver, albedo, local_kabs, local_ksca, local_g, local_planck, local_albedo_2d)
    class(t_solver)                                   :: solver
    real(ireals), intent(in)                          :: albedo
    real(ireals),intent(in),dimension(:,:,:),optional :: local_kabs, local_ksca, local_g ! dimensions (Nz  , Nx, Ny)
    real(ireals),intent(in),dimension(:,:,:),optional :: local_planck                    ! dimensions (Nz+1, Nx, Ny) layer quantity plus surface layer
    real(ireals),intent(in),dimension(:,:),optional   :: local_albedo_2d                 ! dimensions (Nx, Ny)

    real(ireals)        :: tau, kext, w0, g
    integer(iintegers)  :: k, i, j

    associate( atm => solver%atm, &
        C_one_atm => solver%C_one_atm, &
        C_one_atm1 => solver%C_one_atm1, &
        sun => solver%sun, &
        C_one => solver%C_one)

    if(.not.allocated(atm%op) )  allocate( atm%op       (C_one_atm%zs :C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye) )

    if(.not.allocated(atm%albedo)) allocate(atm%albedo(C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
    atm%albedo = albedo
    if(present(local_albedo_2d)) atm%albedo = local_albedo_2d

    if(present(local_kabs) ) atm%op(:,:,:)%kabs = local_kabs
    if(present(local_ksca) ) atm%op(:,:,:)%ksca = local_ksca
    if(present(local_g   ) ) atm%op(:,:,:)%g    = local_g

    if(present(local_planck) ) then
      if (.not.allocated(atm%planck) ) allocate( atm%planck   (C_one_atm1%zs:C_one_atm1%ze, C_one_atm1%xs:C_one_atm1%xe, C_one_atm1%ys:C_one_atm1%ye ) )
      atm%planck = local_planck
      if(atm%lcollapse) then
          !TODO: this does not work at the moment
          print *,'You are trying to collapse the atmosphere in the thermal &
                    &spectral range... this is not possible at the moment or at least not &
                    &tested.'
          ierr = 1; call CHKERR(ierr)
      endif
    else
      if(allocated(atm%planck)) deallocate(atm%planck)
    endif

  !    if(ldebug) then
      if( (any([local_kabs,local_ksca,local_g].lt.zero)) .or. (any(isnan([local_kabs,local_ksca,local_g]))) ) then
        print *,solver%myid,'set_optical_properties :: found illegal value in local_optical properties! abort!'
        do k=C_one_atm%zs,C_one_atm%ze
          print *,solver%myid,k,'local_kabs',local_kabs(k,:,:)
          print *,solver%myid,k,'local_ksca',local_ksca(k,:,:)
        enddo
      endif
  !    endif
  !    if(ldebug) then
      if( (any([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g].lt.zero)) .or. (any(isnan([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g]))) ) then
        print *,solver%myid,'set_optical_properties :: found illegal value in optical properties! abort!'
      endif
  !    endif

    if(ldebug.and.solver%myid.eq.0) then
      if(present(local_kabs) ) then
        print *,'atm_kabs     ',maxval(atm%op%kabs  )  ,shape(atm%op%kabs  )
      endif
      if(present(local_ksca) ) then
        print *,'atm_ksca     ',maxval(atm%op%ksca  )  ,shape(atm%op%ksca  )
      endif
      if(present(local_g) ) then
        print *,'atm_g        ',maxval(atm%op%g     )  ,shape(atm%op%g     )
      endif
      if(present(local_planck) ) then
        print *,'atm_planck   ',maxval(atm%planck   )  ,shape(atm%planck   )
      endif

      print *,'Number of 1D layers: ', count(atm%l1d) , size(atm%l1d),'(',(100._ireals* count(atm%l1d) )/size(atm%l1d),'%)'
      if(present(local_kabs)) print *,'init local optprop:', shape(local_kabs), '::', shape(atm%op)
    endif

    call delta_scale(atm%op(:,:,:)%kabs, atm%op(:,:,:)%ksca, atm%op(:,:,:)%g )

    if(ltwostr_only) then
      if(ldebug .and. solver%myid.eq.0) then
        do k=C_one_atm%zs,C_one_atm%ze
          if(present(local_planck)) then
            print *,solver%myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
              minval(atm%op(k,:,:)%kabs), minval(atm%op(k,:,:)%ksca), minval(atm%op(k,:,:)%g),&
              maxval(atm%op(k,:,:)%kabs), maxval(atm%op(k,:,:)%ksca), maxval(atm%op(k,:,:)%g),&
              '::',minval(atm%planck (k,:,:)),maxval(atm%planck        (k, :,:))
          else
            print *,solver%myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
              minval(atm%op(k,:,:)%kabs), minval(atm%op(k,:,:)%ksca), minval(atm%op(k,:,:)%g),&
              maxval(atm%op(k,:,:)%kabs), maxval(atm%op(k,:,:)%ksca), maxval(atm%op(k,:,:)%g)
          endif
        enddo
      endif
      return ! twostream should not depend on eddington coeffs... it will have to calculate it on its own.
    endif

    ! Make space for deltascaled optical properties
    !if(.not.allocated(atm%delta_op) ) allocate( atm%delta_op (C_one%zs :C_one%ze,C_one%xs :C_one%xe , C_one%ys :C_one%ye ) )
    !atm%delta_op = atm%op
    !call delta_scale(atm%delta_op(:,:,:)%kabs, atm%delta_op(:,:,:)%ksca, atm%delta_op(:,:,:)%g ) !todo should we instead use strong deltascaling? -- what gives better results? or is it as good?

    if(ldebug) then
      if( (any([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g].lt.zero)) .or. (any(isnan([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g]))) ) then
        print *,solver%myid,'set_optical_properties :: found illegal value in delta_scaled optical properties! abort!'
      endif
    endif

    if(luse_eddington) then
      if(.not.allocated(atm%a11) ) allocate(atm%a11 (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))   ! allocate space for twostream coefficients
      if(.not.allocated(atm%a12) ) allocate(atm%a12 (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a21) ) allocate(atm%a21 (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a22) ) allocate(atm%a22 (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a13) ) allocate(atm%a13 (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a23) ) allocate(atm%a23 (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%a33) ) allocate(atm%a33 (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%g1 ) ) allocate(atm%g1  (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
      if(.not.allocated(atm%g2 ) ) allocate(atm%g2  (C_one_atm%zs:C_one_atm%ze ,C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye))
    endif

    if(luse_eddington) then
      do j=C_one_atm%ys,C_one_atm%ye
        do i=C_one_atm%xs,C_one_atm%xe
          do k=C_one_atm%zs,C_one_atm%ze
            if( atm%l1d(k,i,j) ) then
              kext = atm%op(k,i,j)%kabs + atm%op(k,i,j)%ksca
              w0   = atm%op(k,i,j)%ksca / kext
              tau  = atm%dz(k,i,j)* kext
              g    = atm%op(k,i,j)%g
              call eddington_coeff_zdun ( tau , w0, g, sun%angles(C_one_atm%zs,i,j)%costheta, &
                atm%a11(k,i,j),          &
                atm%a12(k,i,j),          &
                atm%a13(k,i,j),          &
                atm%a23(k,i,j),          &
                atm%a33(k,i,j),          &
                atm%g1(k,i,j),           &
                atm%g2(k,i,j) )
            else
              !TODO :: we should really not have this memeory accesible at all....
              !     :: the fix would be trivial at the moment, as long as all 1d layers start at same 'k',
              !     :: however if that is to change i.e. on staggered grid, we would need staggered array constructs....
              if(ldebug) then
                atm%a11(k,i,j) = nil
                atm%a12(k,i,j) = nil
                atm%a13(k,i,j) = nil
                atm%a23(k,i,j) = nil
                atm%a33(k,i,j) = nil
                atm%g1(k,i,j)  = nil
                atm%g2(k,i,j)  = nil
              endif !ldebug
            endif !l1d
          enddo !k

          ! Set symmetric up and down transport coefficients. If only for one
          ! layer, they would be indeed symmetric. If we collapse the atmosphere
          ! however, we have different transmission if we come from top or from
          ! bottom.
          atm%a21 = atm%a12
          atm%a22 = atm%a11
          if(atm%lcollapse) then
              call adding(atm%a11(C_one_atm%zs:atmk(atm, C_one%zs), i, j), &
                          atm%a12(C_one_atm%zs:atmk(atm, C_one%zs), i, j), &
                          atm%a21(C_one_atm%zs:atmk(atm, C_one%zs), i, j), &
                          atm%a22(C_one_atm%zs:atmk(atm, C_one%zs), i, j), &
                          atm%a13(C_one_atm%zs:atmk(atm, C_one%zs), i, j), &
                          atm%a23(C_one_atm%zs:atmk(atm, C_one%zs), i, j), &
                          atm%a33(C_one_atm%zs:atmk(atm, C_one%zs), i, j))
          endif !lcollapse
        enddo !i
      enddo !j
    endif

    if(ldebug .and. solver%myid.eq.0) then
      do k=C_one_atm%zs,C_one_atm%ze
        if(present(local_planck)) then
          print *,solver%myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
            minval(atm%op(k,:,:)%kabs), minval(atm%op(k,:,:)%ksca), minval(atm%op(k,:,:)%g),&
            maxval(atm%op(k,:,:)%kabs), maxval(atm%op(k,:,:)%ksca), maxval(atm%op(k,:,:)%g),&
            '::',minval(atm%planck (k,:,:)),maxval(atm%planck        (k, :,:))
        else
          print *,solver%myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
            minval(atm%op(k,:,:)%kabs), minval(atm%op(k,:,:)%ksca), minval(atm%op(k,:,:)%g),&
            maxval(atm%op(k,:,:)%kabs), maxval(atm%op(k,:,:)%ksca), maxval(atm%op(k,:,:)%g),&
            '::',minval(atm%a33 (k,:,:)),maxval(atm%a33(k,:,:))
        endif
      enddo
    endif
    end associate

    contains
        subroutine adding(a11,a12,a21,a22,a13,a23,a33)
            real(ireals),intent(inout),dimension(:) :: a11,a12,a21,a22,a13,a23,a33
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

    if(.not.solver%linitialized) then
      print *,solver%myid,'You tried to set global optical properties but tenstream environment seems not to be initialized.... please call init first!'
      call exit(1)
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

    if(lhave_planck) then
      call set_optical_properties(solver, local_albedo, local_kabs, local_ksca, local_g, local_planck)
    else
      call set_optical_properties(solver, local_albedo, local_kabs, local_ksca, local_g)
    endif


  contains
    subroutine local_optprop()
      type(tVec) :: local_vec

      if(solver%myid.eq.0.and.ldebug .and. lhave_kabs) &
        print *,solver%myid,'copying optprop: global to local :: shape kabs',shape(global_kabs),'xstart/end',solver%C_one_atm%xs,solver%C_one_atm%xe,'ys/e',solver%C_one_atm%ys,solver%C_one_atm%ye

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
      if(any(shape(arr).lt.minimal_dimension) ) stop 'set_optprop -> extend_arr :: dimension is smaller than we support... please think of something here'
    end subroutine

  end subroutine


  subroutine solve_pprts(solver, edirTOA, opt_solution_uid, opt_solution_time)
    class(t_solver), intent(inout)          :: solver
    real(ireals),intent(in)                 :: edirTOA
    integer(iintegers),optional,intent(in)  :: opt_solution_uid
    real(ireals),      optional,intent(in)  :: opt_solution_time

    integer(iintegers) :: uid
    logical            :: lsolar

    associate(  solutions => solver%solutions, &
                C_dir     => solver%C_dir,     &
                C_diff    => solver%C_diff,    &
                Mdir      => solver%Mdir,      &
                Mdiff     => solver%Mdiff     )

    if(solver%lenable_solutions_err_estimates .and. present(opt_solution_uid)) then
      uid = opt_solution_uid
    else
      uid = i0 ! default solution is uid==0
    endif

    lsolar = mpi_logical_and(solver%comm, edirTOA.gt.zero .and. any(solver%sun%angles%theta.ge.zero))

    call prepare_solution(solver,  solutions(uid), uid, lsolar=lsolar ) ! setup solution vectors

    ! --------- Skip Thermal Computation (-lskip_thermal) --
    if(lskip_thermal .and. (solutions(uid)%lsolar_rad.eqv..False.) ) then !
      if(ldebug .and. solver%myid.eq.0) print *,'skipping thermal calculation -- returning zero flux'
      call VecSet(solutions(uid)%ediff, zero, ierr); call CHKERR(ierr)
      solutions(uid)%lchanged=.True.
      call restore_solution(solver, solutions(uid))
      return
    endif

    ! --------- Calculate 1D Radiative Transfer ------------
    if(  ltwostr                                                         &
      .or. all(solver%atm%l1d.eqv..True.)                                       &
      .or. ((solutions(uid)%lsolar_rad.eqv..False.) .and. lcalc_nca)     &
      .or. ((solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild) ) then

      if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild ) then
        call schwarz(solver, solutions(uid))
      else
        call twostream(solver, edirTOA,  solutions(uid) )
      endif

      if(ldebug .and. solver%myid.eq.0) print *,'1D calculation done'

      if(present(opt_solution_time) ) then
        call restore_solution(solver, solutions(uid),opt_solution_time)
      else
        call restore_solution(solver, solutions(uid))
      endif

        if( ltwostr_only ) return
      if( all(solver%atm%l1d.eqv..True.) ) return
      if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lcalc_nca ) return
      if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild ) return
    endif

    ! --------- scale from [W/m**2] to [W] -----------------
    call scale_flx(solver, solutions(uid), lWm2_to_W=.True. )

    ! ---------------------------- Edir  -------------------
    if( solutions(uid)%lsolar_rad ) then

      call setup_incSolar(solver, solver%incSolar, edirTOA)
      call set_dir_coeff(solver, solver%sun, solver%Mdir, C_dir)
      if(ldebug) call mat_info(solver%comm, solver%Mdir)

      call setup_ksp(solver%atm, kspdir, C_dir, solver%Mdir, linit_kspdir, "dir_")

      call solve(solver, kspdir, solver%incSolar, solutions(uid)%edir)
      solutions(uid)%lchanged=.True.
      solutions(uid)%lintegrated_dir=.True.
      call PetscObjectSetName(solutions(uid)%edir,'debug_edir',ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solutions(uid)%edir, PETSC_NULL_VEC, "-show_debug_edir", ierr); call CHKERR(ierr)
    endif

    ! ---------------------------- Source Term -------------
    call setup_b(solver, solutions(uid),solver%b)

    ! ---------------------------- Ediff -------------------
    call set_diff_coeff(solver, Mdiff,C_diff)
    if(ldebug) call mat_info(solver%comm, solver%Mdiff)

    call setup_ksp(solver%atm, kspdiff,C_diff,Mdiff,linit_kspdiff, "diff_")

    call solve(solver, kspdiff, solver%b, solutions(uid)%ediff,uid)
    solutions(uid)%lchanged=.True.
    solutions(uid)%lintegrated_diff=.True. !Tenstream solver returns fluxes as [W]

    if(present(opt_solution_time) ) then
      call restore_solution(solver, solutions(uid),opt_solution_time)
    else
      call restore_solution(solver, solutions(uid))
    endif

    end associate

  end subroutine

  subroutine prepare_solution(solver, solution, uid, lsolar)
    class(t_solver)               :: solver
    type(t_state_container)       :: solution
    integer(iintegers),intent(in) :: uid
    logical,intent(in)            :: lsolar

    solver%solutions(uid)%uid = uid ! dirty hack to give the solution a unique hash for example to write it out to disk

    if( lsolar .and.  (solution%lsolar_rad.eqv..False.) ) then
      ! set up the direct vectors in any case. This may be necessary even if solution was
      ! already initialized once... e.g. in case we calculated thermal before with same uid
      call DMCreateGlobalVector(solver%C_dir%da ,solution%edir  ,ierr)  ; call CHKERR(ierr)
      call VecSet(solution%edir,zero,ierr); call CHKERR(ierr)
      solution%lsolar_rad = .True.
    endif

    if(ldebug .and. solver%myid.eq.0) print *,'Prepare solution',uid,'::',lsolar,solution%lsolar_rad,'set?',solution%lset

    if(solution%lset) return ! already set-up

    call DMCreateGlobalVector(solver%C_diff%da,solution%ediff ,ierr)  ; call CHKERR(ierr)
    call DMCreateGlobalVector(solver%C_one%da ,solution%abso  ,ierr)  ; call CHKERR(ierr)

    call VecSet(solution%ediff,zero,ierr)    ; call CHKERR(ierr)
    call VecSet(solution%abso,zero,ierr)     ; call CHKERR(ierr)

    solution%lset = .True.
  end subroutine

  !> @brief TODO: renormalize fluxes with the size of a face(sides or lid)
  subroutine scale_flx(solver, solution, lWm2_to_W)
    class(t_solver), intent(in)           :: solver
    type(t_state_container),intent(inout) :: solution   !< @param solution container with computed fluxes
    logical,intent(in)                    :: lWm2_to_W  !< @param determines direction of scaling, if true, scale from W/m**2 to W

    if(solution%lsolar_rad) then
      if(solution%lintegrated_dir .neqv. lWm2_to_W) then
        call scale_dir_flx_vec(solver, solution%edir, solver%C_dir, lWm2_to_W)
        solution%lintegrated_dir = lWm2_to_W
      endif
    endif

    if(solution%lintegrated_diff .neqv. lWm2_to_W) then
      call scale_flx_vec(solver, solution%ediff, solver%C_diff, lWm2_to_W)
      solution%lintegrated_diff = lWm2_to_W
    endif

  contains
    subroutine scale_dir_flx_vec(solver, v, C, lWm2_to_W)
      class(t_solver)      :: solver
      type(tVec)           :: v
      type(t_coord)        :: C
      real(ireals),pointer :: xv  (:,:,:,:) =>null()
      real(ireals),pointer :: xv1d(:)       =>null()
      integer(iintegers)   :: i, j, k, d, iside
      real(ireals)         :: Ax, Ay, Az, fac
      logical,intent(in)   :: lWm2_to_W ! determines direction of scaling, if true, scale from W/m**2 to W

      associate(  atm     => solver%atm,    &
                  C_one1  => solver%C_one1)

      if(solver%myid.eq.0.and.ldebug) print *,'rescaling direct fluxes',C%zm,C%xm,C%ym
      call getVecPointer(v ,C%da ,xv1d, xv)

      Az  = solver%atm%dx*solver%atm%dy / solver%dirtop%dof  ! size of a direct stream in m**2
      if(lWm2_to_W) then
        fac = Az
      else
        fac = one/Az
      endif

      ! Scaling top faces
      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze
            do iside=1,solver%dirtop%dof
              d = iside-1
              xv(d,k,i,j) = xv(d,k,i,j) * fac
            enddo
          enddo
        enddo
      enddo

      ! Scaling side faces
      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-1
            if(.not.atm%l1d(atmk(atm, k),i,j)) then
              ! First the faces in x-direction
              Ax = solver%atm%dy*solver%atm%dz(k,i,j) / (solver%dirside%dof/2)
              if(lWm2_to_W) then
                fac = Ax
              else
                fac = one/Ax
              endif
              do iside=1,solver%dirside%dof/2
                d = solver%dirtop%dof + iside-1
                xv(d,k,i,j) = xv(d,k,i,j) * fac
              enddo

              ! Then the rest of the faces in y-direction
              Ay = atm%dy*atm%dz(k,i,j) / (solver%dirside%dof/2)
              if(lWm2_to_W) then
                fac = Ay
              else
                fac = one/Ay
              endif
              do iside=1,solver%dirside%dof/2
                d = solver%dirtop%dof + solver%dirside%dof/2 + iside-1
                xv(d,k,i,j) = xv(d,k,i,j) * fac
              enddo
            endif
          enddo
        enddo
      enddo

      call restoreVecPointer(v, xv1d, xv )

      end associate
    end subroutine
    subroutine scale_flx_vec(solver, v, C, lWm2_to_W)
      class(t_solver)      :: solver
      type(tVec)           :: v
      type(t_coord)        :: C
      real(ireals),pointer :: xv(:,:,:,:)=>null()
      real(ireals),pointer :: xv1d(:)=>null()
      logical,intent(in)   :: lWm2_to_W ! determines direction of scaling, if true, scale from W/m**2 to W

      integer(iintegers)  :: iside, src, i, j, k
      real(ireals)        :: Az, Ax, Ay, fac
      !Vec                 :: vgrad_x, vgrad_y
      !PetscScalar,Pointer :: grad_x(:,:,:,:)=>null(), grad_x1d(:)=>null()
      !PetscScalar,Pointer :: grad_y(:,:,:,:)=>null(), grad_y1d(:)=>null()
      !real(ireals)        :: grad(3)  ! is the cos(zenith_angle) of the tilted box in case of topography

      if(solver%myid.eq.0.and.ldebug) print *,'rescaling fluxes',C%zm,C%xm,C%ym
      call getVecPointer(v ,C%da ,xv1d, xv)


      if(C%dof.eq.i3 .or. C%dof.eq.i8) then
        print *,'scale_flx_vec is just for diffuse radia'
      endif


      ! Scaling top faces
      Az = solver%atm%dx*solver%atm%dy/(solver%difftop%dof/2)
      if(lWm2_to_W) then
        fac = Az
      else
        fac = one/Az
      endif

      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze
            do iside=1,solver%difftop%dof
              src = iside -1
              xv(src ,k,i,j) = xv(src ,k,i,j) * fac                  ! diffuse radiation
            enddo
          enddo
        enddo
      enddo

      ! Scaling side faces
      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-i1
            if(.not.solver%atm%l1d(atmk(solver%atm, k),i,j)) then

              ! faces in x-direction
              Ax = solver%atm%dy*solver%atm%dz(k,i,j)/(solver%difftop%dof/2)
              if(lWm2_to_W) then
                fac = Ax
              else
                fac = one/Ax
              endif

              do iside=1,solver%diffside%dof/2
                src = solver%difftop%dof + iside -1
                xv(src ,k,i,j) = xv(src ,k,i,j) * fac
              enddo

              ! faces in y-direction
              Ay = solver%atm%dx*solver%atm%dz(k,i,j)/(solver%difftop%dof/2)
              if(lWm2_to_W) then
                fac = Ay
              else
                fac = one/Ay
              endif

              do iside=1,solver%diffside%dof/2
                src = solver%difftop%dof + solver%diffside%dof/2 + iside -1
                xv(src ,k,i,j) = xv(src ,k,i,j) * fac
              enddo

            endif
          enddo
        enddo
      enddo

      if(solver%sun%luse_topography) then ! This is direct rad and we use topography !todo do we need this?
        stop('Dont knwo how i should rescale topography! - exit...')
        !select case (C%dof)

        !case(i8)
        !  call compute_gradient(atm, C_one1, vgrad_x, vgrad_y)

        !  call getVecPointer(vgrad_x , C_one1%da, grad_x1d, grad_x)
        !  call getVecPointer(vgrad_y , C_one1%da, grad_y1d, grad_y)

        !  do j=C%ys,C%ye
        !    do i=C%xs,C%xe
        !      do k=C%zs,C%ze
        !        grad(1) = grad_x(i0,k,i,j)
        !        grad(2) = grad_y(i0,k,i,j)
        !        grad(3) = one

        !        xv(i0:i3,k,i,j) = xv(i0:i3,k,i,j) / norm(grad)
        !      enddo
        !    enddo
        !  enddo

        !  call restoreVecPointer(vgrad_x, grad_x1d, grad_x)
        !  call restoreVecPointer(vgrad_y, grad_y1d, grad_y)

        !  call DMRestoreLocalVector(C_one1%da, vgrad_x, ierr);  call CHKERR(ierr)
        !  call DMRestoreLocalVector(C_one1%da, vgrad_y, ierr);  call CHKERR(ierr)

        !case(i6)
        !  ! Dont rescale diffuse fluxes

        !case(i10)

        !  ! Dont rescale diffuse fluxes

        !  ! call compute_gradient(atm, C_one1, vgrad_x, vgrad_y)

        !  ! call getVecPointer(vgrad_x , C_one1%da, grad_x1d, grad_x)
        !  ! call getVecPointer(vgrad_y , C_one1%da, grad_y1d, grad_y)

        !  !do j=C%ys,C%ye
        !  !  do i=C%xs,C%xe
        !  !    do k=C%zs,C%ze
        !  !      grad(1) = grad_x(i0,k,i,j)
        !  !      grad(2) = grad_y(i0,k,i,j)
        !  !      grad(3) = one

        !  !      xv([E_up, E_dn],k,i,j) = xv([E_up, E_dn],k,i,j) / norm(grad)
        !  !    enddo
        !  !  enddo
        !  !enddo

        !case default
        !  stop('Dont know how I should topography rescale this! - exiting...')
        !end select
      endif

      call restoreVecPointer(v, xv1d, xv )

    end subroutine
  end subroutine

  subroutine restore_solution(solver, solution,time)
    ! restore_solution:: if flux have changed, we need to update absorption, save the residual history
    class(t_solver)         :: solver
    type(t_state_container) :: solution
    real(ireals),intent(in),optional :: time

    character(default_str_len) :: vecname
    real(ireals) :: norm1, norm2, norm3
    type(tVec) :: abso_old

    if( .not. solution%lset ) &
      stop 'cant restore solution that was not initialized'

    if( .not. solution%lchanged ) &
      stop 'cant restore solution which was not changed'

    if(present(time) .and. solver%lenable_solutions_err_estimates) then ! Create working vec to determine difference between old and new absorption vec
      call DMGetGlobalVector(solver%C_one%da, abso_old, ierr) ; call CHKERR(ierr)
      call VecCopy( solution%abso, abso_old, ierr)     ; call CHKERR(ierr)
    endif

    ! make sure to bring the fluxes into [W] before the absorption calculation
    call scale_flx(solver, solution, lWm2_to_W=.True. )

    ! update absorption
    call calc_flx_div(solver, solution)

    if(ldebug .and. solver%myid.eq.0) &
      print *,'Saving Solution ',solution%uid

    ! make sure to bring the fluxes into [W/m**2]
    call scale_flx(solver, solution, lWm2_to_W=.False. )

    if(ldebug .and. solver%myid.eq.0) &
      print *,'Saving Solution done'
    solution%lchanged=.False.

    if(present(time) .and. solver%lenable_solutions_err_estimates) then ! Compute norm between old absorption and new one
      call VecAXPY(abso_old , -one, solution%abso , ierr)    ; call CHKERR(ierr) ! overwrite abso_old with difference to new one
      call VecNorm(abso_old ,  NORM_1, norm1, ierr)          ; call CHKERR(ierr)
      call VecNorm(abso_old ,  NORM_2, norm2, ierr)          ; call CHKERR(ierr)
      call VecNorm(abso_old ,  NORM_INFINITY, norm3, ierr)   ; call CHKERR(ierr)

      call DMRestoreGlobalVector(solver%C_one%da, abso_old, ierr)   ; call CHKERR(ierr)

      ! Save norm for later analysis
      solution%maxnorm = eoshift ( solution%maxnorm, shift = -1) !shift all values by 1 to the right
      solution%twonorm = eoshift ( solution%twonorm, shift = -1) !shift all values by 1 to the right
      solution%time    = eoshift ( solution%time   , shift = -1) !shift all values by 1 to the right

      solution%maxnorm( 1 ) = norm3
      solution%twonorm( 1 ) = norm2
      solution%time( 1 )    = time

      if(ldebug .and. solver%myid.eq.0) &
        print *,'Updating error statistics for solution ',solution%uid,'at time ',time,'::',solution%time(1),':: norm',norm1,norm2,norm3,'[W] :: hr_norm approx:',norm3*86.1,'[K/d]'

    endif !present(time) .and. solver%lenable_solutions_err_estimates


    if(solution%lsolar_rad) then
      write(vecname,FMT='("edir",I0)') solution%uid
      call PetscObjectSetName(solution%edir,vecname,ierr) ; call CHKERR(ierr)
      call PetscObjectViewFromOptions(solution%edir, PETSC_NULL_VEC, "-show_edir", ierr); call CHKERR(ierr)
    endif

    write(vecname,FMT='("ediff",I0)') solution%uid
    call PetscObjectSetName(solution%ediff,vecname,ierr) ; call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%ediff, PETSC_NULL_VEC, "-show_ediff", ierr); call CHKERR(ierr)

    write(vecname,FMT='("abso",I0)') solution%uid
    call PetscObjectSetName(solution%abso,vecname,ierr) ; call CHKERR(ierr)
    call PetscObjectViewFromOptions(solution%abso, PETSC_NULL_VEC, "-show_abso", ierr); call CHKERR(ierr)

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
    real(ireals) :: mu0,incSolar

    associate(atm         => solver%atm, &
              C_diff      => solver%C_diff, &
              C_dir       => solver%C_dir, &
              C_one_atm   => solver%C_one_atm, &
              C_one_atm1  => solver%C_one_atm1)

    if(solution%lsolar_rad) &
      call VecSet(solution%edir ,zero,ierr); call CHKERR(ierr)

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

        mu0 = solver%sun%angles(C_one_atm1%zs,i,j)%costheta
        incSolar = edirTOA* mu0
        if(solver%myid.eq.0 .and. ldebug) then
          print *,' CALCULATING DELTA EDDINGTON TWOSTREAM ::',solver%sun%angles(C_one_atm1%zs,i,j)%theta,':',incSolar
        endif

        kext = atm%op(:,i,j)%kabs + atm%op(:,i,j)%ksca
        dtau = atm%dz(:,i,j)* kext
        w0   = atm%op(:,i,j)%ksca / kext
        g    = atm%op(:,i,j)%g

        if(allocated(atm%planck) ) then
          call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup, planck=atm%planck(:,i,j) )
        else
          call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup )
        endif

        if(solution%lsolar_rad) then
          do src=i0,solver%dirtop%dof-1
            xv_dir(src,C_dir%zs+1:C_dir%ze,i,j) = S(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze)
            xv_dir(src,C_dir%zs           ,i,j) = S(C_one_atm1%zs)
          enddo
        endif

        xv_diff(E_up,C_diff%zs+1:C_diff%ze,i,j) = Eup(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze)
        xv_diff(E_up,C_diff%zs            ,i,j) = Eup(C_one_atm1%zs)
        xv_diff(E_dn,C_diff%zs+1:C_diff%ze,i,j) = Edn(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze)
        xv_diff(E_dn,C_diff%zs            ,i,j) = Edn(C_one_atm1%zs)
      enddo
    enddo

    if(solution%lsolar_rad) &
      call restoreVecPointer(solution%edir, xv_dir1d, xv_dir  )
    call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

    !Twostream solver returns fluxes as [W]
    solution%lintegrated_dir  = .False.
    solution%lintegrated_diff = .False.
    ! and mark solution that it is not up to date
    solution%lchanged         = .True.

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
    integer(iintegers) :: i,j

    real(ireals),allocatable :: dtau(:),Edn(:),Eup(:)


    associate( C_diff   => solver%C_diff, &
                atm     => solver%atm, &
                C_one   => solver%C_one, &
                C_one1  => solver%C_one1)

    if(solution%lsolar_rad) stop 'Tried calling schwarschild solver for solar calculation -- stopping!'
    if( .not. allocated(atm%planck) ) stop 'Tried calling schwarschild solver but no planck was given -- stopping!'

    call VecSet(solution%ediff,zero,ierr); call CHKERR(ierr)

    allocate( dtau(C_diff%zm-1) )

    call getVecPointer(solution%ediff ,C_diff%da ,xv_diff1d, xv_diff)

    allocate( Eup(C_diff%zm) )
    allocate( Edn(C_diff%zm) )

    if(solver%myid.eq.0 .and. ldebug) print *,' CALCULATING schwarzschild ::'

    do j=C_diff%ys,C_diff%ye
      do i=C_diff%xs,C_diff%xe

        dtau = atm%dz(atmk(atm, C_one%zs):C_one%ze,i,j)* atm%op(atmk(atm, C_one%zs):C_one%ze,i,j)%kabs

        call schwarzschild( dtau ,atm%albedo(i,j), Edn,Eup, atm%planck(atmk(atm, C_one1%zs):C_one1%ze,i,j) )

        xv_diff(E_up,:,i,j) = Eup(:)
        xv_diff(E_dn,:,i,j) = Edn(:)
      enddo
    enddo

    call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

    !Schwarzschild solver returns fluxes as [W/m^2]
    solution%lintegrated_dir  = .False.
    solution%lintegrated_diff = .False.
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
    class(t_solver)                       :: solver
    type(t_state_container)               :: solution

    real(ireals),pointer,dimension(:,:,:,:)  :: xediff=>null(),xedir=>null(),xabso=>null()
    real(ireals),pointer,dimension(:)        :: xediff1d=>null(),xedir1d=>null(),xabso1d=>null()

    integer(iintegers)        :: offset, isrc, src
    integer(iintegers)        :: i, j, k, xinc, yinc
    type(tVec)                :: ledir,lediff ! local copies of vectors, including ghosts
    real(ireals)              :: div(3)
    real(ireals), allocatable :: div2(:)
    real(ireals)              :: Volume,Az
    logical                   :: lhave_no_3d_layer

    associate(  atm     => solver%atm, &
                C_dir   => solver%C_dir, &
                C_diff  => solver%C_diff, &
                C_one   => solver%C_one, &
                C_one1  => solver%C_one1)

    if(solution%lsolar_rad .and. (solution%lintegrated_dir .eqv..False.)) stop 'tried calculating absorption but dir  vector was in [W/m**2], not in [W], scale first!'
    if(                          (solution%lintegrated_diff.eqv..False.)) stop 'tried calculating absorption but diff vector was in [W/m**2], not in [W], scale first!'

    if( (solution%lsolar_rad.eqv..False.) .and. lcalc_nca ) then ! if we should calculate NCA (Klinger), we can just return afterwards
      call scale_flx(solver, solution, lWm2_to_W=.False.)
      call nca_wrapper(solution%ediff, solution%abso)
      call scale_flx(solver, solution, lWm2_to_W=.True.)
      return
    endif

    if(solver%myid.eq.0.and.ldebug) print *,'Calculating flux divergence',solution%lsolar_rad.eqv..False.,lcalc_nca
    call VecSet(solution%abso,zero,ierr) ;call CHKERR(ierr)

    ! if there are no 3D layers globally, we should skip the ghost value copying....
    lhave_no_3d_layer = mpi_logical_and(solver%comm, all(atm%l1d.eqv..True.))
    if(lhave_no_3d_layer) then
      call scale_flx(solver, solution, lWm2_to_W=.True.)

      if(solution%lsolar_rad) call getVecPointer(solution%edir, C_dir%da ,xedir1d ,xedir )

      call getVecPointer(solution%ediff, C_diff%da, xediff1d, xediff)
      call getVecPointer(solution%abso, C_one%da, xabso1d, xabso)

      ! calculate absorption by flux divergence
      Az = atm%dx * atm%dy

      do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe
          do k=C_one%zs,C_one%ze
            Volume = Az     * atm%dz(atmk(atm, k),i,j)
            ! Divergence    =                       Incoming                -       Outgoing
            div(1) = zero
            if(solution%lsolar_rad) then
              do src=i0,solver%dirtop%dof-1
                div(1) = div(1) + xedir(src, k, i, j )  - xedir(src , k+i1 , i, j )
              enddo
            endif

            div(2) = ( xediff(E_up  ,k+1,i  ,j  )  - xediff(E_up  ,k  ,i  ,j  )  )
            div(3) = ( xediff(E_dn  ,k  ,i  ,j  )  - xediff(E_dn  ,k+1,i  ,j  )  )

            xabso(i0,k,i,j) = sum(div) / Volume
          enddo
        enddo
      enddo

      if(solution%lsolar_rad) call restoreVecPointer(solution%edir, xedir1d, xedir )

      call restoreVecPointer(solution%ediff, xediff1d, xediff)
      call restoreVecPointer(solution%abso, xabso1d ,xabso)

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
    Az = atm%dx * atm%dy

    allocate(div2(solver%C_dir%dof + solver%C_diff%dof))

    do j=C_one%ys,C_one%ye
      do i=C_one%xs,C_one%xe
        do k=C_one%zs,C_one%ze

          Volume = Az     * atm%dz(atmk(atm, k),i,j)
          ! Divergence = Incoming - Outgoing

          if(atm%l1d(atmk(atm, k),i,j)) then ! one dimensional i.e. twostream
            div(1) = zero
            if(solution%lsolar_rad) then
              do src=1,solver%dirtop%dof
                div(1) = div(1) + xedir(src-1, k, i, j )  - xedir(src-1 , k+i1 , i, j )
              enddo
            endif

            div(2) = ( xediff(E_up  ,k+1,i  ,j  )  - xediff(E_up  ,k  ,i  ,j  )  )
            div(3) = ( xediff(E_dn  ,k  ,i  ,j  )  - xediff(E_dn  ,k+1,i  ,j  )  )

            xabso(i0,k,i,j) = sum(div) / Volume

          else ! 3D-radiation
            offset = solver%dirtop%dof + solver%dirside%dof*2

            if(solution%lsolar_rad) then

              xinc = solver%sun%angles(k,i,j)%xinc
              yinc = solver%sun%angles(k,i,j)%yinc

              do isrc = 1,solver%dirtop%dof
                src = isrc
                div2(src) = xedir(src-1, k, i        , j         ) - xedir(src-1, k+i1 , i      , j      )
              enddo

              do isrc = 1, solver%dirside%dof
                src = isrc + solver%dirtop%dof
                div2(src) = xedir(src-1, k, i+1-xinc , j         ) - xedir(src-1, k    , i+xinc , j      )
              enddo

              do isrc = 1, solver%dirside%dof
                src = isrc + solver%dirtop%dof + solver%dirside%dof
                div2(src) = xedir(src-1, k, i        , j+i1-yinc ) - xedir(src-1, k    , i      , j+yinc )
              enddo

            else
              div2(1:offset) = zero
            endif

            ! diffuse part of absorption
            do isrc = 1, solver%difftop%dof
              src = isrc
              if (solver%difftop%is_inward(isrc) .eqv. .True.) then
                div2(offset+src) = xediff(src-1, k   , i   , j   ) - xediff(src-1, k+1 , i   , j   )
              else
                div2(offset+src) = xediff(src-1, k+1 , i   , j   ) - xediff(src-1, k   , i   , j   )
              endif
            enddo

            do isrc = 1, solver%diffside%dof
              src = isrc + solver%difftop%dof
              if (solver%diffside%is_inward(isrc) .eqv. .True.) then
                div2(offset+src) = xediff(src-1, k   , i   , j   ) - xediff(src-1, k   , i+1 , j   )
              else
                div2(offset+src) = xediff(src-1, k   , i+1 , j   ) - xediff(src-1, k   , i   , j   )
              endif
            enddo

            do isrc = 1, solver%diffside%dof
              src = isrc + solver%difftop%dof + solver%diffside%dof
              if (solver%diffside%is_inward(isrc) .eqv. .True.) then
                div2(offset+src) = xediff(src-1, k   , i   , j   ) - xediff(src-1, k   , i   , j+1 )
              else
                div2(offset+src) = xediff(src-1, k   , i   , j+1 ) - xediff(src-1, k   , i   , j   )
              endif
            enddo

            xabso(i0,k,i,j) = sum(div2) / Volume
            if( isnan(xabso(i0,k,i,j)) ) then
              print *,'nan in flxdiv',k,i,j,'::',xabso(i0,k,i,j),Volume,'::',div2
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

  end associate
end subroutine


!> @brief call PETSc Krylov Subspace Solver
!> @details solve with ksp and save residual history of solver
!> \n -- this may be handy later to decide next time if we have to calculate radiation again
!> \n if we did not get convergence, we try again with standard GMRES and a resetted(zero) initial guess -- if that doesnt help, we got a problem!
subroutine solve(solver, ksp,b,x,solution_uid)
  class(t_solver) :: solver
  type(tKSP) :: ksp
  type(tVec) :: b
  type(tVec) :: x
  integer(iintegers),optional,intent(in) :: solution_uid

  KSPConvergedReason :: reason
  integer(iintegers) :: iter

  KSPType :: old_ksp_type

  if(solver%myid.eq.0.and.ldebug) print *,'Solving Matrix'

  if(present(solution_uid)) then
    if(.not.allocated( solver%solutions(solution_uid)%ksp_residual_history) ) allocate(solver%solutions(solution_uid)%ksp_residual_history(100) )
    solver%solutions(solution_uid)%ksp_residual_history = -1
    call KSPSetResidualHistory(ksp, solver%solutions(solution_uid)%ksp_residual_history, 100_iintegers, .True.,ierr)
  endif

  call KSPSolve(ksp,b,x,ierr) ;call CHKERR(ierr)
  call KSPGetIterationNumber(ksp,iter,ierr) ;call CHKERR(ierr)
  call KSPGetConvergedReason(ksp,reason,ierr) ;call CHKERR(ierr)

  ! if(reason.eq.KSP_DIVERGED_ITS) then
  !   if(solver%myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
  !   return
  ! endif

  if(reason.le.0) then
    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Resetted initial guess to zero and try again with gmres:'
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
    if(solver%myid.eq.0.and.ldebug) print *,solver%myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
  endif

  if(reason.le.0) then
    if(solver%myid.eq.0) print *,'***** SOLVER did NOT converge :( ********',reason
    call exit()
  endif
end subroutine

!> @brief initialize PETSc Krylov Subspace Solver
!> @details default KSP solver is a FBCGS with BJCAOBI // ILU(1)
!> \n -- the default does however not scale well -- and we configure petsc solvers per commandline anyway
!> \n -- see documentation for details on how to do so
subroutine setup_ksp(atm, ksp,C,A,linit, prefix)
  type(t_atmosphere) :: atm
  type(tKSP) :: ksp
  type(t_coord) :: C
  type(tMat) :: A
  type(tPC)  :: prec
  logical :: linit

  type(tMatNullSpace) :: nullspace
  type(tVec) :: nullvecs(0)
  character(len=*),optional :: prefix

  real(ireals),parameter :: rtol=sqrt(epsilon(rtol))*10, rel_atol=1e-4_ireals
  integer(iintegers),parameter  :: maxiter=1000

  integer(mpiint) :: myid, numnodes

  real(ireals) :: atol

  logical,parameter :: lset_geometry=.True.  ! this may be necessary in order to use geometric multigrid
! logical,parameter :: lset_geometry=.False.  ! this may be necessary in order to use geometric multigrid
! logical,parameter :: lset_nullspace=.True. ! set constant nullspace?
  logical,parameter :: lset_nullspace=.False. ! set constant nullspace?

  if(linit) return

  call mpi_comm_rank(C%comm, myid, ierr)     ; call CHKERR(ierr)
  call mpi_comm_size(C%comm, numnodes, ierr) ; call CHKERR(ierr)

  call imp_allreduce_min(C%comm, rel_atol*(C%dof*C%glob_xm*C%glob_ym*C%glob_zm) * count(.not.atm%l1d)/(one*size(atm%l1d)), atol)
  atol = max(1e-8_ireals, atol)

  if(myid.eq.0.and.ldebug) &
    print *,'Setup KSP -- tolerances:',rtol,atol,'::',rel_atol,(C%dof*C%glob_xm*C%glob_ym*C%glob_zm),count(.not.atm%l1d),one*size(atm%l1d)

  call KSPCreate(C%comm,ksp,ierr) ;call CHKERR(ierr)
  if(present(prefix) ) call KSPAppendOptionsPrefix(ksp,trim(prefix),ierr) ;call CHKERR(ierr)

  call KSPSetType(ksp,KSPBCGS,ierr)  ;call CHKERR(ierr)
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

      call DMDASetUniformCoordinates(C%da,zero,one, zero, atm%dx*C%glob_xm,zero, atm%dy*C%glob_ym, ierr );call CHKERR(ierr)
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
      initial_rnorm=max(epsilon(rnorm),rnorm)
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

    fac = edirTOA * (solver%atm%dx*solver%atm%dy) / solver%dirtop%dof

    call VecSet(incSolar,zero,ierr) ;call CHKERR(ierr)

    call getVecPointer(incSolar, solver%C_dir%da, x1d, x4d)

    do j=solver%C_dir%ys,solver%C_dir%ye
      do i=solver%C_dir%xs,solver%C_dir%xe
        do src=1, solver%dirtop%dof
          x4d(src-1,solver%C_dir%zs,i,j) = fac * solver%sun%angles(solver%C_dir%zs,i,j)%costheta
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
      real(ireals)       :: v(C%dof**2),norm

      integer(iintegers) :: dst,src, xinc, yinc, isrc, idst

      xinc = sun%angles(k,i,j)%xinc
      yinc = sun%angles(k,i,j)%yinc

      do idst = 1, solver%dirtop%dof
        dst = idst
        row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
      enddo

      do idst = 1, solver%dirside%dof
        dst = idst + solver%dirtop%dof
        row(MatStencil_j,dst) = i+xinc   ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
      enddo

     do idst = 1, solver%dirside%dof
       dst = idst + solver%dirtop%dof + solver%dirside%dof
       row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j+yinc   ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
     enddo

      do isrc = 1, solver%dirtop%dof
        src = isrc
        col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k; col(MatStencil_c,src) = src-i1 ! Define transmission towards the lower/upper lid
      enddo

      do isrc = 1, solver%dirside%dof
        src = isrc + solver%dirtop%dof
        col(MatStencil_j,src) = i+1-xinc   ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Define transmission towards the left/right lid
      enddo

     do isrc = 1, solver%dirside%dof
       src = isrc + solver%dirtop%dof + solver%dirside%dof
       col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j+1-yinc   ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Define transmission towards the front/back lid
     enddo


      call get_coeff(solver, solver%atm%op(atmk(solver%atm,k),i,j), solver%atm%dz(atmk(solver%atm,k),i,j), .True., v, &
                             solver%atm%l1d(atmk(solver%atm,k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta], &
                             lswitch_east=xinc.eq.0, lswitch_north=yinc.eq.0)

      call MatSetValuesStencil(A, C%dof, row, C%dof, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)

      if(ldebug) then
        do src=1,C%dof
          norm = sum( v(src:C%dof**2:C%dof) )
          if( norm.gt.one+10._ireals*epsilon(one) ) then
            print *,'direct sum(dst==',dst,') gt one',norm
            print *,'direct coeff',norm,'::',v
            stop 'omg.. shouldnt be happening'
            ierr = -5
            return
          endif
        enddo
      endif

    end subroutine

    subroutine set_eddington_coeff(atm,A,k,i,j)
      type(t_atmosphere), intent(inout) :: atm
      type(tMat),intent(inout)          :: A
      integer(iintegers),intent(in)     :: i,j,k

      MatStencil :: row(4,1), col(4,1)
      real(ireals) :: v(1)
      integer(iintegers) :: src

      if(luse_eddington) then
        v = atm%a33(atmk(atm,k),i,j)
      else
        call get_coeff(solver, atm%op(atmk(atm,k),i,j), atm%dz(atmk(atm,k),i,j),.True., v, atm%l1d(atmk(atm,k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta] )
      endif

      col(MatStencil_j,i1) = i      ; col(MatStencil_k,i1) = j       ; col(MatStencil_i,i1) = k
      row(MatStencil_j,i1) = i      ; row(MatStencil_k,i1) = j       ; row(MatStencil_i,i1) = k+1

      do src=i1,solver%dirtop%dof
        col(MatStencil_c,i1) = src-i1 ! Source may be the upper/lower lid:
        row(MatStencil_c,i1) = src-i1 ! Define transmission towards the lower/upper lid
        call MatSetValuesStencil(A,i1, row,i1, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)
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

      real(ireals) :: Ax,Ay,Az,b0 !,c1,c2,c3,b1,dtau
      real(ireals) :: diff2diff1d(4)
      real(ireals) :: diff2diff(solver%C_diff%dof**2),v(solver%C_diff%dof**2)
      integer(iintegers) :: k,i,j,src,iside

      associate(  atm     => solver%atm, &
                C_diff  => solver%C_diff)

      if(solver%myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector ... setting thermal source terms', minval(atm%planck), maxval(atm%planck)
      Az = atm%dx*atm%dy

      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          do k=C_diff%zs,C_diff%ze-1

            if( atm%l1d(atmk(atm,k),i,j) ) then

              if(luse_eddington ) then

                b0 = atm%planck(atmk(atm,k),i,j) * (one-atm%a11(atmk(atm,k),i,j)-atm%a12(atmk(atm,k),i,j))
                xsrc(E_up   ,k  ,i,j) = xsrc(E_up   ,k  ,i,j) + b0 *Az*pi
                xsrc(E_dn   ,k+1,i,j) = xsrc(E_dn   ,k+1,i,j) + b0 *Az*pi

              else
                call get_coeff(solver, atm%op(atmk(atm,k),i,j), atm%dz(atmk(atm,k),i,j),.False., diff2diff1d, atm%l1d(atmk(atm,k),i,j))

                b0 = atm%planck(atmk(atm,k),i,j) * pi
                xsrc(E_up   ,k  ,i,j) = xsrc(E_up   ,k  ,i,j) +  b0  *(one-diff2diff1d(1)-diff2diff1d(2) ) *Az
                xsrc(E_dn   ,k+1,i,j) = xsrc(E_dn   ,k+1,i,j) +  b0  *(one-diff2diff1d(1)-diff2diff1d(2) ) *Az
              endif

            else ! Tenstream source terms
              Ax = atm%dy*atm%dz(atmk(atm,k),i,j)
              Ay = atm%dx*atm%dz(atmk(atm,k),i,j)

              call get_coeff(solver, atm%op(atmk(atm,k),i,j), atm%dz(atmk(atm,k),i,j),.False., diff2diff, atm%l1d(atmk(atm,k),i,j) )
              ! reorder from destination ordering to src ordering
              do src=1,C_diff%dof
                v(src:C_diff%dof**2:C_diff%dof) = diff2diff( i1+(src-i1)*C_diff%dof : src*C_diff%dof )
              enddo

              b0 = atm%planck(atmk(atm,k),i,j) * pi
              do src=1,solver%difftop%dof
                if (solver%difftop%is_inward(src) .eqv. .False.) then
                  xsrc(src-1, k, i, j) = xsrc(src-1, k, i, j) +  b0 *(one-sum(v((src-1)*C_diff%dof+1:src*C_diff%dof )  )  ) *Az
                else
                  xsrc(src-1, k+1, i, j) = xsrc(src-1, k+1, i, j) +  b0 *(one-sum(v((src-1)*C_diff%dof+1:src*C_diff%dof )  )  ) *Az
                endif
              enddo

              do iside=1,solver%diffside%dof
                src = iside+solver%difftop%dof
                if (solver%diffside%is_inward(iside) .eqv. .False.) then
                  xsrc(src-1, k, i, j) = xsrc(src-1, k, i, j) +  b0 *(one-sum(v((src-1)*C_diff%dof+1:src*C_diff%dof )  )  ) *Ax
                else
                  xsrc(src-1, k, i+1, j) = xsrc(src-1, k, i+1, j) +  b0 *(one-sum(v((src-1)*C_diff%dof+1:src*C_diff%dof )  )  ) *Ax
                endif
              enddo

              do iside=1,solver%diffside%dof
                src = iside + solver%difftop%dof + solver%diffside%dof
                if (solver%diffside%is_inward(iside) .eqv. .False.) then
                  xsrc(src-1, k, i, j) = xsrc(src-1, k, i, j) +  b0 *(one-sum(v((src-1)*C_diff%dof+1:src*C_diff%dof )  )  ) *Ay
                else
                  xsrc(src-1, k, i, j+1) = xsrc(src-1, k, i, j+1) +  b0 *(one-sum(v((src-1)*C_diff%dof+1:src*C_diff%dof )  )  ) *Ay
                endif
              enddo

            endif ! 1D or Tenstream?

          enddo
        enddo

      enddo ! k

      ! Thermal emission at surface
      k = C_diff%ze
      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          xsrc(E_up   ,k,i,j) = xsrc(E_up   ,k,i,j) + atm%planck(atmk(atm,k),i,j)*Az *(one-atm%albedo(i,j))*pi
        enddo
      enddo
    end associate
    end subroutine

    subroutine set_solar_source(solver, edir)
      class(t_solver)     :: solver
      type(tVec)          :: edir

      real(ireals)        :: dir2diff(solver%C_dir%dof*solver%C_diff%dof), solrad
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

      if(solver%myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector .. setting solar source',sum(xedir(i0,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye))/size(xedir(i0,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye))

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
                        xsrc(idiff-1,k+1,i,j) = xsrc(idiff-1,k+1,i,j) +  xedir(src-1,k,i,j) * atm%a23(atmk(atm,k),i,j)/(solver%difftop%dof/2)
                      else
                        ! fetch all diffuse upward fluxes at k
                        xsrc(idiff-1,k,i,j) = xsrc(idiff-1,k,i,j) +  xedir(src-1,k,i,j) * atm%a13(atmk(atm,k),i,j)/(solver%difftop%dof/2)
                      endif
                    enddo
                  enddo

                else
                  !call get_coeff(atm%op(atmk(atm,k),i,j), atm%dz(atmk(atm,k),i,j),.False., twostr_coeff, atm%l1d(atmk(atm,k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta])
                  stop 'set solar source only implemented for use with eddington coeff'
                endif


              else ! Tenstream source terms
                lsun_east  = sun%angles(k,i,j)%xinc.eq.i0
                lsun_north = sun%angles(k,i,j)%yinc.eq.i0

                call get_coeff(solver, atm%op(atmk(atm,k),i,j), atm%dz(atmk(atm,k),i,j),.False., dir2diff, atm%l1d(atmk(atm,k),i,j), &
                               [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta], lswitch_east=lsun_east, lswitch_north=lsun_north)

                do idofdst=1,solver%difftop%dof
                  dst = idofdst
                  do src=1, solver%dirtop%dof
                    solrad = xedir( src-1 , k , i , j )

                    if (solver%difftop%is_inward(idofdst)) then
                      xsrc(dst-1,k+1,i,j)= xsrc(dst-1,k+1,i,j) + solrad * dir2diff((dst-1)*C_dir%dof+src)
                    else
                      xsrc(dst-1,k,i,j)= xsrc(dst-1,k,i,j) + solrad * dir2diff((dst-1)*C_dir%dof+src)
                    endif
                  enddo

                  do idof = 1, solver%dirside%dof
                    src = idof + solver%dirtop%dof
                    solrad = xedir( src-1, k , i+i1-sun%angles(k,i,j)%xinc , j )
                    if (solver%difftop%is_inward(dst)) then
                      xsrc(dst-1,k+1,i,j)= xsrc(dst-1,k+1,i,j) + solrad * dir2diff((dst-1)*C_dir%dof+src)
                    else
                      xsrc(dst-1,k,i,j)= xsrc(dst-1,k,i,j) + solrad * dir2diff((dst-1)*C_dir%dof+src)
                    endif
                  enddo

                  do idof = 1, solver%dirside%dof
                    src = idof + solver%dirtop%dof + solver%dirside%dof
                    solrad = xedir( src-1 , k , i , j+i1-sun%angles(k,i,j)%yinc )
                    if (solver%difftop%is_inward(dst)) then
                      xsrc(dst-1,k+1,i,j)= xsrc(dst-1,k+1,i,j) + solrad * dir2diff((dst-1)*C_dir%dof+src)
                    else
                      xsrc(dst-1,k,i,j)= xsrc(dst-1,k,i,j) + solrad * dir2diff((dst-1)*C_dir%dof+src)
                    endif
                  enddo
                enddo

                do idofdst = 1, solver%diffside%dof
                  dst = idofdst + solver%difftop%dof

                  do idof = 1, solver%dirtop%dof
                    src = idof
                    solrad = xedir( src-1 , k , i , j )
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst-1, k, i+1, j) = xsrc(dst-1, k, i+1, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    else
                      xsrc(dst-1, k, i, j) = xsrc(dst-1, k, i, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    endif
                  enddo

                  do idof = 1, solver%dirside%dof
                    src = idof + solver%dirtop%dof
                    solrad = xedir(src-1, k, i+i1-sun%angles(k,i,j)%xinc , j )
                    if (solver%diffside%is_inward(idofdst)) then
                      xsrc(dst-1, k, i+1, j) = xsrc(dst-1, k, i+1, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    else
                      xsrc(dst-1, k, i, j) = xsrc(dst-1, k, i, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    endif
                  enddo

                  do idof = 1, solver%dirside%dof
                    src = idof + solver%dirtop%dof + solver%dirside%dof
                    solrad = xedir( src-1 , k , i , j+i1-sun%angles(k,i,j)%yinc )
                    if (solver%diffside%is_inward(idofdst)) then
                      xsrc(dst-1, k, i+1, j) = xsrc(dst-1, k, i+1, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    else
                      xsrc(dst-1, k, i, j) = xsrc(dst-1, k, i, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    endif
                  enddo
                enddo

                do idofdst = 1, solver%diffside%dof
                  dst = idofdst + solver%difftop%dof + solver%diffside%dof
                  do src = 1, solver%dirtop%dof
                    solrad = xedir(src-1, k, i, j)
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst-1, k, i, j+1) = xsrc(dst-1, k, i, j+1) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    else
                      xsrc(dst-1, k, i, j) = xsrc(dst-1, k, i, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    endif
                  enddo

                  do idof = 1, solver%dirside%dof
                    src = idof + solver%dirtop%dof
                    solrad = xedir(src-1, k, i+i1-sun%angles(k,i,j)%xinc , j)
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst-1, k, i, j+1) = xsrc(dst-1, k, i, j+1) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    else
                      xsrc(dst-1, k, i, j) = xsrc(dst-1, k, i, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    endif
                  enddo

                  do idof = 1, solver%dirside%dof
                    src = idof + solver%dirtop%dof + solver%dirside%dof
                    solrad = xedir( src-1 , k , i , j+i1-sun%angles(k,i,j)%yinc )
                    if(solver%diffside%is_inward(idofdst)) then
                      xsrc(dst-1, k, i, j+1) = xsrc(dst-1, k, i, j+1) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    else
                      xsrc(dst-1, k, i, j) = xsrc(dst-1, k, i, j) + solrad * dir2diff((dst-1)*C_dir%dof + src)
                    endif
                  enddo
                enddo

                if(ldebug) then
                  do src=1,C_dir%dof
                  if(sum(dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof)) .gt. one .or. &
                    sum(dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof)) .lt. zero   ) &
                    print *,'DEBUG Found dir2diff gt one:',src,'::',sum(dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof)),':',dir2diff(src : C_dir%dof*C_diff%dof : C_dir%dof) ,'   :::::::     ', dir2diff
                  enddo
                endif

              endif ! 1D or Tenstream?
            endif ! if solar

          enddo
        enddo
      enddo

      ! Ground Albedo reflecting direct radiation, the diffuse part is considered by the solver(Matrix)
      k = C_diff%ze
      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          do dst=1, solver%difftop%dof
            if (.not. solver%difftop%is_inward(dst)) then
              do src=1, solver%dirtop%dof
                xsrc(dst-1,k,i,j) = xsrc(dst-1, k, i, j) + xedir(src-1,k,i,j) * atm%albedo(i,j) *(2/solver%difftop%dof)
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
  subroutine get_coeff(solver, op,dz,ldir,coeff,lone_dimensional,angles, lswitch_east, lswitch_north)
    class(t_solver), intent(inout)    :: solver
    type(t_opticalprops),intent(in)   :: op
    real(ireals),intent(in)           :: dz
    logical,intent(in)                :: ldir
    real(ireals),intent(out)          :: coeff(:)

    logical,intent(in)                :: lone_dimensional
    real(ireals),intent(in),optional  :: angles(2)
    logical,intent(in),optional       :: lswitch_east, lswitch_north

    real(ireals) :: aspect, tauz, w0

    aspect = dz / solver%atm%dx
    tauz = (op%kabs+op%ksca) * dz
    w0 = op%ksca / (op%kabs+op%ksca)

    if(lone_dimensional) then
      stop
      !call OPP_1_2%get_coeff (aspect, tauz, w0, op%g,ldir,coeff,angles)
    else
      call solver%OPP%get_coeff(aspect, tauz, w0, op%g,ldir,coeff,angles, lswitch_east, lswitch_north)
    endif

  end subroutine


  !> @brief nca wrapper to call NCA of Carolin Klinger
  !> @details This is supposed to work on a 1D solution which has to be calculated beforehand
  !> \n the wrapper copies fluxes and optical properties on one halo and then gives that to NCA
  !> \n the result is the 3D approximation of the absorption, considering neighbouring information
  subroutine nca_wrapper(ediff, abso)
    use m_ts_nca, only : ts_nca
    type(tVec) :: ediff, abso
!    type(tVec) :: gnca ! global nca vector
!    type(tVec) :: lnca ! local nca vector with ghost values -- in dimension 0 and 1 are fluxes followed by dz,planck,kabs
!
!    real(ireals),pointer,dimension(:,:,:,:) :: xv  =>null()
!    real(ireals),pointer,dimension(:)       :: xv1d=>null()
!    real(ireals),pointer,dimension(:,:,:,:) :: xvlnca  =>null(), xvgnca  =>null()
!    real(ireals),pointer,dimension(:)       :: xvlnca1d=>null(), xvgnca1d=>null()
!    real(ireals),pointer,dimension(:,:,:,:) :: xhr  =>null()
!    real(ireals),pointer,dimension(:)       :: xhr1d=>null()
!    integer(iintegers) :: k
!
!    integer(iintegers),parameter :: idz=i2, iplanck=i3, ikabs=i4, ihr=i5

    stop 'nca_wrapper not implemented'
    print *, 'DEBUG: Stupid print statement to prevent unused compiler warnings', ediff, abso

!   ! put additional values into a local ediff vec .. TODO: this is a rather dirty hack but is straightforward

!   ! get ghost values for dz, planck, kabs and fluxes, ready to give it to NCA
!   call DMGetGlobalVector(C_diff%da ,gnca ,ierr) ; call CHKERR(ierr)

!   call getVecPointer(gnca ,C_diff%da ,xvgnca1d, xvgnca)
!   xvgnca(  idz    , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%dz
!   xvgnca(  iplanck, C_diff%zs:C_diff%ze  , C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%planck
!   xvgnca(  ikabs  , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%op%kabs


!   ! Copy Edn and Eup to local convenience vector
!   call getVecPointer(ediff ,C_diff%da ,xv1d, xv)
!   xvgnca( E_up,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_up,:,:,:)
!   xvgnca( E_dn,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_dn,:,:,:)
!   call restoreVecPointer(ediff, xv1d, xv)

!   call restoreVecPointer(gnca, xvgnca1d, xvgnca )


!   ! retrieve ghost values into l(ocal) nca vec
!   call DMGetLocalVector (C_diff%da ,lnca ,ierr) ; call CHKERR(ierr)
!   call VecSet(lnca, zero, ierr); call CHKERR(ierr)

!   call DMGlobalToLocalBegin(C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)
!   call DMGlobalToLocalEnd  (C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)

!   call DMRestoreGlobalVector(C_diff%da, gnca, ierr); call CHKERR(ierr)

!   ! call NCA
!   call getVecPointer(lnca ,C_diff%da ,xvlnca1d, xvlnca)

!   call ts_nca( atm%dx, atm%dy,                    &
!     xvlnca(   idz        , : , : , :), &
!     xvlnca(   iplanck    , : , : , :), &
!     xvlnca(   ikabs      , : , : , :), &
!     xvlnca(   E_dn       , : , : , :), &
!     xvlnca(   E_up       , : , : , :), &
!     xvlnca(   ihr        , : , : , :))


!   ! return absorption
!   call getVecPointer( abso, C_one%da ,xhr1d, xhr)

!   do k=C_one%zs,C_one%ze
!     xhr(i0,k,:,:) = xvlnca( ihr , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) / xvlnca( idz , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye)
!   enddo
!   call restoreVecPointer(abso, xhr1d, xhr )

!   !return convenience vector that holds optical properties
!   call restoreVecPointer(lnca, xvlnca1d, xvlnca )
!   call DMRestoreLocalVector(C_diff%da, lnca, ierr); call CHKERR(ierr)

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
      real(ireals) :: v(C%dof**2),norm

      integer(iintegers) :: dst,src,idof

      do idof=1, solver%difftop%dof
        src = idof-1
        if (solver%difftop%is_inward(idof)) then
          col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
        else
          col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k+1   ; col(MatStencil_c,src) = src
        endif
      enddo

      do idof=1, solver%diffside%dof
        src = solver%difftop%dof + idof -1
        if (solver%diffside%is_inward(idof)) then
          col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
        else
          col(MatStencil_j,src) = i+1  ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
        endif
      enddo

      do idof=1, solver%diffside%dof
      src = solver%difftop%dof + solver%diffside%dof + idof -1
        if (solver%diffside%is_inward(idof)) then
          col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
        else
          col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j+1   ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = src
        endif
      enddo

      do idof=1, solver%difftop%dof
        dst = idof-1
        if (solver%difftop%is_inward(idof)) then
          row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k+1   ; row(MatStencil_c,dst) = dst
        else
          row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
        endif
      enddo

      do idof=1, solver%diffside%dof
        dst = solver%difftop%dof + idof-1
        if (solver%diffside%is_inward(idof)) then
          row(MatStencil_j,dst) = i+1  ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
        else
          row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
        endif
      enddo

      do idof=1, solver%diffside%dof
        dst = solver%difftop%dof + solver%diffside%dof + idof-1
        if (solver%diffside%is_inward(idof)) then
          row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j+1   ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
        else
          row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = dst
        endif
      enddo

      call get_coeff(solver, solver%atm%op(atmk(solver%atm, k),i,j), solver%atm%dz(atmk(solver%atm, k),i,j),.False., v, solver%atm%l1d(atmk(solver%atm, k),i,j))
      call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)

      if(ldebug) then
        do src=1,C%dof
          norm = sum( v(src:C%dof**2:C%dof) )
          if( norm.gt.one+10._ireals*epsilon(one) ) then
            print *,'diffuse sum(src==',src,') gt one',norm
            stop 'omg.. shouldnt be happening'
            ierr = -5
            return
          endif
        enddo
      endif

    end subroutine

    subroutine set_eddington_coeff(atm,A,k,i,j)
      type(t_atmosphere)            :: atm
      type(tMat),intent(inout)      :: A
      integer(iintegers),intent(in) :: k,i,j

      MatStencil         :: row(4,2)  ,col(4,2)
      real(ireals)       :: v(4),twostr_coeff(4) ! v ==> a12,a11,a22,a21
      integer(iintegers) :: src,dst

      if(luse_eddington ) then
        v = [ atm%a12(atmk(atm,k),i,j), atm%a11(atmk(atm,k),i,j), atm%a22(atmk(atm,k),i,j), atm%a21(atmk(atm,k),i,j)]
      else
          call get_coeff(solver, atm%op(atmk(atm,k),i,j), atm%dz(atmk(atm,k),i,j),.False., twostr_coeff, atm%l1d(atmk(atm,k),i,j)) !twostr_coeff ==> a12,a11,a12,a11 !todo: check if this is the wrong order
          v = [ twostr_coeff(1), twostr_coeff(2) , twostr_coeff(2) , twostr_coeff(1) ]
      endif

      src = 1; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_dn
      src = 2; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k+1   ; col(MatStencil_c,src) = E_up

      dst = 1; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_up
      dst = 2; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k+1   ; row(MatStencil_c,dst) = E_dn

      call MatSetValuesStencil(A,i2, row, i2, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)


    end subroutine

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
                  call MatSetValuesStencil(A,i1, row, i1, col , [-solver%atm%albedo(i,j)/(solver%difftop%dof/2)] ,INSERT_VALUES,ierr) ;call CHKERR(ierr)
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
    if(.not.solver%lenable_solutions_err_estimates) uid = 0

    if(ldebug .and. solver%myid.eq.0) print *,'calling pprts_get_result',present(redir),'for uid',uid

    if(solver%solutions(uid)%lchanged) stop 'tried to get results from unrestored solution -- call restore_solution first'

    if(allocated(redn )) stop 'pprts_get_result :: you should not call it with an allocated redn  array'
    if(allocated(reup )) stop 'pprts_get_result :: you should not call it with an allocated reup  array'
    if(allocated(rabso)) stop 'pprts_get_result :: you should not call it with an allocated rabso array'

    if(present(redir)) then
      if(allocated(redir)) stop 'pprts_get_result :: you should not call it with an allocated redir array'
      allocate(redir(solver%C_dir%zm, solver%C_dir%xm, solver%C_dir%ym))

      if( .not. solver%solutions(uid)%lsolar_rad ) then
        print *,'Hey, You called pprts_get_result for uid',uid,'but in this particular band we dont have direct radiation calculated... I will return with edir=0 but are you sure this is what you intended?'
        redir = zero
      else
        if( solver%solutions(uid)%lintegrated_dir ) stop 'tried to get result from integrated result vector(dir)'
        if(solver%atm%lcollapse) stop 'pprts_get_result :: lcollapse needs to be implemented'

        call getVecPointer(solver%solutions(uid)%edir, solver%C_dir%da, x1d, x4d)
        redir = sum(x4d(0:solver%dirtop%dof-1, :, :, :), dim=1) / real(solver%dirtop%dof, kind=ireals)  ! average of direct radiation of all fluxes through top faces
        call restoreVecPointer(solver%solutions(uid)%edir, x1d, x4d)

        if(ldebug) then
          if(solver%myid.eq.0) print *,'Edir vertically first column',redir(:, lbound(redir,2), lbound(redir,3))
          if(any(redir.lt.-one)) then
            print *,'Found direct radiation smaller than 0 in dir result... that should not happen',minval(redir)
            call exit(1)
          endif
        endif
      endif
    endif

    if(solver%solutions(uid)%lintegrated_diff) stop 'tried to get result from integrated result vector(diff)'

    allocate(redn(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))
    allocate(reup(solver%C_diff%zm, solver%C_diff%xm, solver%C_diff%ym))

    if(solver%atm%lcollapse) then
      stop 'pprts_get_result :: lcollapse needs to be implemented'
    else
      redn = zero
      reup = zero
      call getVecPointer(solver%solutions(uid)%ediff, solver%C_diff%da, x1d, x4d)
      do iside=1,solver%difftop%dof
        if(solver%difftop%is_inward(iside)) then
          redn = redn + x4d(iside-1, :, :, :)
        else
          reup = reup + x4d(iside-1, :, :, :)
        endif
      enddo
      call restoreVecPointer(solver%solutions(uid)%ediff,x1d,x4d)
      reup = reup / (solver%difftop%dof / 2)
      redn = redn / (solver%difftop%dof / 2)
    endif

    if(solver%myid.eq.0 .and. ldebug .and. present(redir)) &
      print *,'mean surface Edir',mean(redir(ubound(redir,1),:,:))
    if(solver%myid.eq.0 .and. ldebug) print *,'mean surface Edn',mean(redn(ubound(redn,1), :,:))
    if(solver%myid.eq.0 .and. ldebug) print *,'mean surface Eup',mean(reup(ubound(reup,1), :,:))

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

    allocate(rabso(solver%C_one%zm, solver%C_one%xm, solver%C_one%ym))
    call getVecPointer(solver%solutions(uid)%abso, solver%C_one%da, x1d, x4d)
    if(solver%atm%lcollapse) then
      rabso(1:atmk(solver%atm, solver%C_one%zs)+1, :, :) = zero
      rabso(atmk(solver%atm, solver%C_one%zs)+2 :solver%C_one_atm%ze+1, :, :) = x4d(i0,solver%C_one%zs+1:solver%C_one%ze,:,:)
    else
      rabso = x4d(i0,:,:,:)
    endif
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

    if(solver%myid.eq.0) then
      call check_arr_size(solver%C_one_atm1, gedn )
      call check_arr_size(solver%C_one_atm1, geup )
      call check_arr_size(solver%C_one_atm , gabso)
      if(present(gedir)) call check_arr_size(solver%C_one_atm1, gedir)
    endif

    if(present(gedir)) then
      call pprts_get_result(solver,redn,reup,rabso,redir=redir,opt_solution_uid=opt_solution_uid)
      call exchange_var(solver%C_one_atm1, redir, gedir)
    else
      call pprts_get_result(solver,redn,reup,rabso,opt_solution_uid=opt_solution_uid)
    endif

    call exchange_var(solver%C_one_atm1, redn , gedn )
    call exchange_var(solver%C_one_atm1, reup , geup )
    call exchange_var(solver%C_one_atm , rabso, gabso)

    if(solver%myid.eq.0 .and. ldebug) then
      print *,'Retrieving results:'
      if(present(gedir)) print *,sum(gedir)/size(gedir)
      print *,sum(gedn) /size(gedn)
      print *,sum(geup) /size(geup)
      print *,sum(gabso)/size(gabso)
    endif

  contains
    subroutine exchange_var(C, inp, outp)
      type(t_coord),intent(in) :: C
      real(ireals),intent(in), allocatable :: inp(:,:,:) ! local array from get_result
      real(ireals),intent(inout),allocatable :: outp(:,:,:) ! global sized array on rank 0

      type(tVec) :: vec, lvec_on_zero
      print *,solver%myid,'exchange_var',allocated(inp), allocated(outp)
      print *,solver%myid,'exchange_var shape',shape(inp)

      call DMGetGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)
      call f90VecToPetsc(inp, C%da, vec)
      call petscGlobalVecToZero(vec, C%da, lvec_on_zero)
      call DMRestoreGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)

      if(solver%myid.eq.0) then
        call petscVecToF90(lvec_on_zero, C%da, outp, opt_l_only_on_rank0=.True.)
        call VecDestroy(lvec_on_zero, ierr); call CHKERR(ierr)
      endif
    end subroutine
    subroutine check_arr_size(C,inp)
      type(t_coord),intent(in) :: C
      real(ireals),intent(in), allocatable :: inp(:,:,:)
      if(allocated(inp)) then
        if(size(inp).ne.C%dof*C%glob_zm*C%glob_xm*C%glob_ym) then
          print *,'pprts_get_result_toZero was called with an already allocated output array but it has the wrong size.'
          print *,'while I could just re-allocate it, this may not be just what you intended. Please do that yourself.'
          print *,'Size of global Dimensions of simulation:',C%dof*C%glob_zm*C%glob_xm*C%glob_ym, '.vs. your input:',size(inp)
          stop 'pprts_get_result_toZero :: should not be called with already allocated array with wrong size'
        endif
      endif
    end subroutine
  end subroutine

  subroutine destroy_pprts(solver, lfinalizepetsc)
    class(t_solver)   :: solver
    logical,optional :: lfinalizepetsc
    logical :: lfinalize
    integer(iintegers) :: uid
    lfinalize = get_arg(.False., lfinalizepetsc)

    if(solver%linitialized) then
      if(linit_kspdir) then
        call KSPDestroy(kspdir , ierr) ;call CHKERR(ierr); linit_kspdir =.False.
      endif
      if(linit_kspdiff) then
        call KSPDestroy(kspdiff, ierr) ;call CHKERR(ierr); linit_kspdiff=.False.
      endif

      if(.not. ltwostr_only) then
        call VecDestroy(solver%incSolar , ierr) ;call CHKERR(ierr)
        call VecDestroy(solver%b        , ierr) ;call CHKERR(ierr)
        deallocate(solver%incSolar)
        deallocate(solver%b)
      endif
      call destroy_matrices(solver)

      do uid=lbound(solver%solutions,1),ubound(solver%solutions,1)
        if( solver%solutions(uid)%lset ) then
          if(solver%solutions(uid)%lsolar_rad) then
            call VecDestroy(solver%solutions(uid)%edir , ierr) ;call CHKERR(ierr)
            solver%solutions(uid)%lsolar_rad = .False.
          endif

          call VecDestroy(solver%solutions(uid)%ediff    , ierr) ;call CHKERR(ierr)
          call VecDestroy(solver%solutions(uid)%abso     , ierr) ;call CHKERR(ierr)

          if(allocated(solver%solutions(uid)%ksp_residual_history)) &
            deallocate(solver%solutions(uid)%ksp_residual_history)

          solver%solutions(uid)%lset = .False.
        endif
      enddo

      if(allocated(solver%atm)) deallocate(solver%atm)

      if(allocated(solver%sun%angles)) deallocate(solver%sun%angles)

      call solver%OPP%destroy()
      call DMDestroy(solver%C_dir%da ,ierr); deallocate(solver%C_dir )
      call DMDestroy(solver%C_diff%da,ierr); deallocate(solver%C_diff)
      call DMDestroy(solver%C_one%da ,ierr); deallocate(solver%C_one )
      call DMDestroy(solver%C_one1%da,ierr); deallocate(solver%C_one1)
      call DMDestroy(solver%C_one_atm%da ,ierr); deallocate(solver%C_one_atm)
      call DMDestroy(solver%C_one_atm1%da,ierr); deallocate(solver%C_one_atm1)

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

