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
    imp_comm, myid, numnodes, init_mpi_data_parameters, mpiint,  &
    zero, one, nil, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi, &
    default_str_len
  
  use m_helper_functions, only : CHKERR, deg2rad, rad2deg, norm, imp_allreduce_max, &
    delta_scale
  
  use m_optprop, only: t_optprop, t_optprop_1_2, t_optprop_3_6, t_optprop_8_10
  use m_eddington, only : eddington_coeff_zdun

  use m_tenstream_options, only : read_commandline_options, ltwostr, luse_eddington, twostr_ratio, &
    options_max_solution_err, options_max_solution_time, ltwostr_only, luse_twostr_guess,        &
    options_phi, lforce_phi, options_theta, lforce_theta, &
    lcalc_nca, lskip_thermal, lschwarzschild, ltopography

  implicit none
  private

  public :: t_solver_1_2, t_solver_3_6, t_solver_8_10, init_pprts, &
            set_optical_properties
  
  type t_opticalprops
    real(ireals) :: kabs,ksca,g
  end type

  type t_atmosphere
    type(t_opticalprops) , allocatable , dimension(:,:,:) :: op
    real(ireals)    , allocatable , dimension(:,:,:) :: planck
    real(ireals)    , allocatable , dimension(:,:,:) :: a11, a12, a21, a22, a13, a23, a33
    real(ireals)    , allocatable , dimension(:,:,:) :: g1,g2
    real(ireals)    , allocatable , dimension(:,:,:) :: dz
    logical         , allocatable , dimension(:,:,:) :: l1d
    real(ireals)    , allocatable , dimension(:,:)   :: albedo
    real(ireals) :: dx,dy
    integer(iintegers) :: icollapse=1
    logical :: lcollapse = .False.
  end type


  type t_coord
    PetscInt :: xs,xe                 ! local domain start and end indices
    PetscInt :: ys,ye                 ! local domain start and end indices
    PetscInt :: zs,ze                 ! local domain start and end indices
    PetscInt :: xm,ym,zm              ! size of local domain
    PetscInt :: gxs,gys,gzs           ! domain indices including ghost points
    PetscInt :: gxe,gye,gze           !
    PetscInt :: gxm,gym,gzm           ! size of local domain including ghosts
    PetscInt :: glob_xm,glob_ym,glob_zm ! global domain size
    PetscInt :: dof,dim               ! degrees of freedom of Petsc Domain, dimension of dmda
    type(tDM) :: da                   ! The Domain Decomposition Object
    PetscMPIInt,allocatable :: neighbors(:) ! all 3d neighbours( (x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
  end type

  type t_sunangles
    real(ireals) :: symmetry_phi
    integer(iintegers) :: yinc,xinc
    real(ireals) :: theta, phi, costheta, sintheta
  end type
  type t_suninfo
    type(t_sunangles),allocatable :: angles(:,:,:) ! defined on DMDA grid
    logical :: luse_topography=.False.
  end type

  type, abstract :: t_solver
    type(t_coord), allocatable :: C_dir, C_diff, C_one, C_one1, C_one_atm, C_one_atm1
    type(t_atmosphere),allocatable :: atm
    type(t_suninfo) :: sun
    type(tMat),allocatable :: Mdir,Mdiff
    class(t_optprop), allocatable :: OPP

    integer(iintegers) :: dofdiff(2), dofdir(2)
    logical :: lenable_solutions_err_estimates=.True.  ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
    type(tVec),allocatable :: incSolar,b

    logical :: linitialized=.False.
  end type

  type, extends(t_solver) :: t_solver_1_2
  end type 
  type, extends(t_solver) :: t_solver_8_10
  end type 
  type, extends(t_solver) :: t_solver_3_6
  end type 


  logical,parameter :: ldebug=.True.
  logical,parameter :: lcycle_dir=.True.
  logical,parameter :: lprealloc=.False.

  integer(iintegers),parameter :: minimal_dimension=3 ! this is the minimum number of gridpoints in x or y direction

  PetscErrorCode :: ierr


contains

  !> @brief Main routine to setup TenStream solver
  !> @details This will setup the PETSc DMDA grid and set other grid information, needed for the TenStream
  !> \n Nx, Ny Nz are either global domain size or have to be local sizes if present(nxproc,nyproc)
  !> \n where nxproc and nyproc then are the number of pixel per rank for all ranks -- i.e. sum(nxproc) != Nx_global
  subroutine init_pprts(icomm, Nz,Nx,Ny, dx,dy, phi0, theta0, solver, dz1d, dz3d, nxproc, nyproc, collapseindex)
    MPI_Comm, intent(in)          :: icomm   !< @param MPI_Communicator which should be used -- this will be used for PETSC_COMM_WORLD
    integer(iintegers),intent(in) :: Nz      !< @param[in] Nz     Nz is the number of layers and Nz+1 would be the number of levels
    integer(iintegers),intent(in) :: Nx      !< @param[in] Nx     number of boxes in x-direction
    integer(iintegers),intent(in) :: Ny      !< @param[in] Ny     number of boxes in y-direction
    real(ireals),intent(in)       :: dx      !< @param[in] dx     physical size of grid in [m]
    real(ireals),intent(in)       :: dy      !< @param[in] dy     physical size of grid in [m]
    real(ireals),intent(in)       :: phi0    !< @param[in] phi0   solar azmiuth and zenith angle
    real(ireals),intent(in)       :: theta0  !< @param[in] theta0 solar azmiuth and zenith angle

    class(t_solver), intent(inout):: solver         !< @param[inout] solver  
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
        solver%dofdiff = [2,0]
        solver%dofdir = [1,0]

        class is (t_solver_3_6)
        solver%dofdiff = [2,2]
        solver%dofdir = [1,1]

        class is (t_solver_8_10)
        solver%dofdiff = [2,4]
        solver%dofdir = [4,2]

        class default
        stop 'init pprts: unexpected type for solver'
      end select

      PETSC_COMM_WORLD = icomm

      !      call PetscInitialize(tenstreamrc ,ierr) ;call CHKERR(ierr)
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr) ;call CHKERR(ierr)
#ifdef _XLF
      call PetscPopSignalHandler(ierr); call CHKERR(ierr) ! in case of xlf ibm compilers, remove petsc signal handler -- otherwise we dont get fancy signal traps from boundschecking or FPE's
#endif
      call init_mpi_data_parameters(icomm)

      call read_commandline_options()

      if(present(nxproc) .and. present(nyproc) ) then
        if(ldebug.and.myid.eq.0) print *,'nxproc',shape(nxproc),'::',nxproc
        if(ldebug.and.myid.eq.0) print *,'nyproc',shape(nyproc),'::',nyproc
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
      print *,myid,'You tried to initialize already initialized Tenstream          &
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
      PetscInt,intent(in) :: Nz_in,Nx,Ny !< @param[in] local number of grid boxes -- in the vertical we have Nz boxes and Nz+1 levels
      integer(iintegers),optional :: nxproc(:), nyproc(:) ! size of local domains on each node
      integer(iintegers),optional,intent(in) :: collapseindex  !< @param[in] collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them

      DMBoundaryType :: bp=DM_BOUNDARY_PERIODIC, bn=DM_BOUNDARY_NONE, bg=DM_BOUNDARY_GHOSTED
      integer(iintegers) :: Nz

      Nz = Nz_in
      if(present(collapseindex)) Nz = Nz_in-collapseindex+i1

      if(myid.eq.0.and.ldebug) print *,myid,'Setting up the DMDA grid for ',Nz,Nx,Ny,'using ',numnodes,' nodes'

      if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_diff'
      call setup_dmda(solver%C_diff, Nz+1,Nx,Ny, bp, sum(solver%dofdiff))

      if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_dir'
      if(lcycle_dir) then
        call setup_dmda(solver%C_dir, Nz+1,Nx,Ny, bp, sum(solver%dofdir))
      else
        call setup_dmda(solver%C_dir, Nz+1,Nx,Ny, bg, sum(solver%dofdir))
      endif

      if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_one'
      call setup_dmda(solver%C_one , Nz  , Nx,Ny,  bp, i1)
      call setup_dmda(solver%C_one1, Nz+1, Nx,Ny,  bp, i1)

      if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA atm'
      call setup_dmda(solver%C_one_atm , Nz_in  , Nx,Ny,  bp, i1)
      call setup_dmda(solver%C_one_atm1, Nz_in+1, Nx,Ny,  bp, i1)

      if(myid.eq.0.and.ldebug) print *,myid,'DMDA grid ready'
    contains
      subroutine setup_dmda(C, Nz, Nx, Ny, boundary, dof)
        type(t_coord),allocatable :: C
        PetscInt,intent(in) :: Nz,Nx,Ny,dof
        DMBoundaryType,intent(in) :: boundary

        PetscInt,parameter :: stencil_size=1

        allocate(C)

        C%dof = i1*dof
        if(present(nxproc) .and. present(nyproc) ) then
          call DMDACreate3d( imp_comm,                                      &
            bn                      , boundary             , boundary     , &
            DMDA_STENCIL_STAR       ,                                       &
            Nz                      , i1*sum(nxproc)          , i1*sum(nyproc)  , &
            i1                      , i1*size(nxproc)         , i1*size(nyproc) , &
            C%dof                   , stencil_size         ,                &
            Nz                      , nxproc               , nyproc       , &
            C%da                    , ierr)
          call CHKERR(ierr)
        else
          call DMDACreate3d( imp_comm,                                            &
            bn                      , boundary             , boundary           , &
            DMDA_STENCIL_STAR       ,                                             &
            i1*Nz                   , Nx                   , Ny                 , &
            i1                      , PETSC_DECIDE         , PETSC_DECIDE       , &
            C%dof                   , stencil_size         ,                      &
            Nz                      , PETSC_NULL_INTEGER   , PETSC_NULL_INTEGER , &
            C%da                    , ierr) ;call CHKERR(ierr)
        endif

        call DMSetup(C%da,ierr) ;call CHKERR(ierr)
        call DMSetMatType(C%da, MATAIJ, ierr); call CHKERR(ierr)
        if(lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ;call CHKERR(ierr)

        call DMSetFromOptions(C%da, ierr) ; call CHKERR(ierr)
        if(ldebug) call DMView(C%da, PETSC_VIEWER_STDOUT_WORLD ,ierr)
        call setup_coords(C)
      end subroutine
      subroutine setup_coords(C)
        type(t_coord) :: C

        call DMDAGetInfo(C%da,C%dim,                               &
          C%glob_zm,C%glob_xm,C%glob_ym,                           &
          PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
          PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
          PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
          ierr) ;call CHKERR(ierr)

        call DMDAGetCorners(C%da, C%zs, C%xs,C%ys, C%zm, C%xm,C%ym, ierr) ;call CHKERR(ierr)
        C%xe = C%xs+C%xm-1
        C%ye = C%ys+C%ym-1
        C%ze = C%zs+C%zm-1

        call DMDAGetGhostCorners(C%da,C%gzs,C%gxs,C%gys,C%gzm,C%gxm,C%gym,ierr) ;call CHKERR(ierr)
        C%gxe = C%gxs+C%gxm-1
        C%gye = C%gys+C%gym-1
        C%gze = C%gzs+C%gzm-1

        if(ldebug) then
          print *,myid,'Domain Corners z:: ',C%zs,':',C%ze,' (',C%zm,' entries)','global size',C%glob_zm
          print *,myid,'Domain Corners x:: ',C%xs,':',C%xe,' (',C%xm,' entries)','global size',C%glob_xm
          print *,myid,'Domain Corners y:: ',C%ys,':',C%ye,' (',C%ym,' entries)','global size',C%glob_ym
        endif

        allocate(C%neighbors(0:3**C%dim-1) )
        call DMDAGetNeighbors(C%da,C%neighbors,ierr) ;call CHKERR(ierr)
        !if(numnodes.gt.i1) then
        !  if(ldebug.and.(C%dim.eq.3)) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([10,4,16,22]),'while I am ',C%neighbors(13)
        !  if(ldebug.and.(C%dim.eq.2)) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([1,3,7,5]),'while I am ',C%neighbors(4)
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
          print *,myid,'You tried to set angles in the Tenstream solver.  &
              & This should be called right after init_tenstream'
          ierr=1; call CHKERR(ierr)
      endif

      if(allocated(solver%sun%angles)) then ! was initialized
        if(present(theta2d)) then
          lchanged_theta = .not. all(theta2d.eq.solver%sun%angles(1,:,:)%theta)
        else
          lchanged_theta = .not. all(theta0.eq.solver%sun%angles(1,:,:)%theta)
        endif
        if(present(phi2d)) then
          lchanged_phi = .not. all(phi2d.eq.solver%sun%angles(1,:,:)%phi)
        else
          lchanged_phi = .not. all(phi0.eq.solver%sun%angles(1,:,:)%phi)
        endif
        if(myid.eq.0 .and. ldebug) print *,'tenstr set_angles -- changed angles?',lchanged_theta, lchanged_phi
        if(.not. lchanged_theta .and. .not. lchanged_phi) then
          return
        endif
      endif


      if ( present(phi2d) .and. present(theta2d) ) then
          call setup_suninfo(phi0, theta0, solver%sun, solver%C_one, phi2d=phi2d, theta2d=theta2d)
      elseif ( present(phi2d) ) then
          call setup_suninfo(phi0, theta0, solver%sun, solver%C_one, phi2d=phi2d)
      elseif ( present(theta2d) ) then
          call setup_suninfo(phi0, theta0, solver%sun, solver%C_one, theta2d=theta2d)
      else
          call setup_suninfo(phi0, theta0, solver%sun, solver%C_one)
      endif


      ! init box montecarlo model
        select type(solver)
          class is (t_solver_1_2)
             if(.not.allocated(solver%OPP) ) allocate(t_optprop_1_2::solver%OPP)

          class is (t_solver_3_6)
             if(.not.allocated(solver%OPP) ) allocate(t_optprop_3_6::solver%OPP)

          class is (t_solver_8_10)
             if(.not.allocated(solver%OPP) ) allocate(t_optprop_8_10::solver%OPP)

          class default
          stop 'init pprts: unexpected type for solver'
        end select

        call solver%OPP%init(pack(solver%sun%angles%symmetry_phi,.True.), &
          pack(solver%sun%angles%theta,.True.),imp_comm )

!     allocate(solver%OPP)
!     call solver%OPP%init(pack(solver%sun%angles%symmetry_phi,.True.), &
!                          pack(solver%sun%angles%theta,.True.),imp_comm, &
!                          Ndiff=sum(solver%dofdiff), Ndir=sum(solver%dofdir))

!     if(any(solver%atm%l1d.eqv..False.)) call OPP_3_6%init(pack(solver%sun%angles%symmetry_phi,.True.),pack(solver%sun%angles%theta,.True.),imp_comm)
!  !   if(any(solver%atm%l1d.eqv..False.)) call OPP_8_10%init(pack(solver%sun%angles%symmetry_phi,.True.),pack(solver%sun%angles%theta,.True.),imp_comm)
!     if(.not.luse_eddington)      call OPP_1_2%init (pack(solver%sun%angles%symmetry_phi,.True.),pack(solver%sun%angles%theta,.True.),imp_comm)

      call init_Matrix(solver%Mdir , solver%C_dir )
      call init_Matrix(solver%Mdiff, solver%C_diff)
  end subroutine


  !> @brief set direction where sun stands
  !> @details save sun azimuth and zenith angle
  !>   \n sun azimuth is reduced to the range of [0,90] and the transmission of direct radiation is contributed for by a integer increment,
  !>   \n determining which neighbouring box is used in the horizontal direction
  subroutine setup_suninfo(phi0, theta0, sun, C_one, phi2d, theta2d, solver)
    real(ireals),intent(in) :: phi0, theta0
    type(t_suninfo),intent(inout) :: sun
    type(t_coord), intent(in) :: C_one
    real(ireals), optional, intent(in) :: phi2d(:,:), theta2d(:,:)
    class(t_solver), optional, intent(in) :: solver


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
      if (.not.present(solver)) then
        print *,'You want to use topography but didnt call setup_suninfo with solver argument'
        call exit
      endif
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

    !if(ldebug.and.myid.eq.0) print *,myid,'setup_dir_inc done', &
    if(ldebug) print *,myid,'setup_dir_inc done', &
      count(sun%angles%xinc.eq.0),count(sun%angles%xinc.eq.i1), &
      count(sun%angles%yinc.eq.0),count(sun%angles%yinc.eq.i1), &
      '::', minval(sun%angles%xinc), maxval(sun%angles%xinc), minval(sun%angles%yinc), maxval(sun%angles%yinc)

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

    Vec :: vgrad_x, vgrad_y
    PetscScalar,Pointer :: grad_x(:,:,:,:)=>null(), grad_x1d(:)=>null()
    PetscScalar,Pointer :: grad_y(:,:,:,:)=>null(), grad_y1d(:)=>null()

    integer(iintegers) :: i,j,k
    real(ireals) :: newtheta, newphi, xsun(3)
    real(ireals) :: rotmat(3, 3), newxsun(3)

    if(.not.allocated(sun%angles)) stop 'You called  setup_topography() &
        &but the sun struct is not yet up, make sure setup_suninfo is called before'
    if(.not.allocated(atm%dz)) stop 'You called  setup_topography() &
        &but the atm struct is not yet up, make sure we have atm%dz before'

    call compute_gradient(atm, C_one1, vgrad_x, vgrad_y)

    call getVecPointer(vgrad_x , C_one1, grad_x1d, grad_x)
    call getVecPointer(vgrad_y , C_one1, grad_y1d, grad_y)


    rotmat = reshape((/ one , zero, zero,  &
                        zero, one , zero,  &
                        nil , nil , one/), &
                     (/3, 3/), order=(/2, 1/) )

    do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe

            ! if we are at the global boundary we have to take that the gradient does not get too
            ! steep. That wouldnt make sense for cyclic boundaries... ! todo: check what we should do?
            !if(i.eq.i0 .or. i.eq.C_one%glob_xm-i1 .or. j.eq.i0 .or. j.eq.C_one%glob_ym-i1) then
            !   !print *,myid,'global edge at:',i,j
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

                ! if(i.eq.C_one1%xs) print *,myid,i,j,k, '::',hhl(i0,k+1,i,j-1:j+1),'::', grad ,'::',sun%angles(k,i,j)%theta, newtheta, '::', sun%angles(k,i,j)%phi, newphi

                sun%angles(k,i,j)%theta = max(zero, min( 90._ireals, newtheta ))
                sun%angles(k,i,j)%phi = newphi
            enddo
        enddo
    enddo

    call restoreVecPointer(vgrad_x , C_one1, grad_x1d, grad_x)
    call restoreVecPointer(vgrad_y , C_one1, grad_y1d, grad_y)

    call DMRestoreLocalVector(C_one1%da, vgrad_x, ierr);  call CHKERR(ierr)
    call DMRestoreLocalVector(C_one1%da, vgrad_y, ierr);  call CHKERR(ierr)

    if(myid.eq.0 .and. ldebug) print *,'min,max theta',minval(sun%angles%theta), maxval(sun%angles%theta)
    if(myid.eq.0 .and. ldebug) print *,'min,max phi',minval(sun%angles%phi), maxval(sun%angles%phi)
  end subroutine

subroutine getVecPointer(vec,C,x1d,x4d)
  type(tVec) :: vec
  type(t_coord),intent(in) :: C
  PetscScalar,intent(inout),pointer,dimension(:,:,:,:) :: x4d
  PetscScalar,intent(inout),pointer,dimension(:) :: x1d

  integer(iintegers) :: N
  logical :: lghosted

  if(associated(x1d).or.associated(x4d)) then
    print *,'ERROR : getVecPointer : input vector already associated!!',associated(x1d),associated(x4d)
    call sleep(30)
    call exit(1)
  endif

  call VecGetLocalSize(vec,N,ierr)

  if( N .eq. C%dof*C%xm*C%ym*C%zm ) then
    lghosted=.False.
  else if( N .eq. C%dof*C%gxm*C%gym*C%gzm ) then
    lghosted=.True.
  else
    stop 'Local Vector dimensions does not conform to DMDA size'
  endif

  call VecGetArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
  if(lghosted) then
    x4d(0:C%dof-1 , C%gzs:C%gze, C%gxs:C%gxe , C%gys:C%gye ) => x1d
  else                                                     
    x4d(0:C%dof-1 , C%zs:C%ze  , C%xs:C%xe   , C%ys:C%ye   ) => x1d
  endif

end subroutine
subroutine restoreVecPointer(vec,C,x1d,x4d)
  type(tVec) :: vec
  type(t_coord),intent(in) :: C
  PetscScalar,intent(inout),pointer,dimension(:,:,:,:) :: x4d
  PetscScalar,intent(inout),pointer,dimension(:) :: x1d

  if(.not.associated(x1d).or..not.associated(x4d)) then
    print *,'ERROR : restoreVecPointer : input vector not yet associated!!',associated(x1d),associated(x4d)
    call exit(1)
  endif

  x4d => null()
  call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
  x1d => null()
end subroutine

  !> @brief compute gradient from dz3d
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information
  subroutine compute_gradient(atm, C_one1, vgrad_x, vgrad_y)
    type(t_atmosphere),intent(in) :: atm
    type(t_coord), intent(in) :: C_one1
    type(tVec) :: vgrad_x, vgrad_y

    type(tVec) :: vhhl
    PetscScalar,Pointer :: hhl(:,:,:,:)=>null(), hhl1d(:)=>null()
    PetscScalar,Pointer :: grad_x(:,:,:,:)=>null(), grad_x1d(:)=>null()
    PetscScalar,Pointer :: grad_y(:,:,:,:)=>null(), grad_y1d(:)=>null()

    integer(iintegers) :: i,j,k

    real(ireals) :: zm(4), maxheight, global_maxheight

    if(.not.allocated(atm%dz)) stop 'You called  compute_gradient()&
      &but the atm struct is not yet up, make sure we have atm%dz before'

    call DMGetLocalVector(C_one1%da, vhhl, ierr) ;call CHKERR(ierr)
    call getVecPointer(vhhl , C_one1, hhl1d, hhl)

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
    call imp_allreduce_max(imp_comm, maxheight, global_maxheight)

    do j=C_one1%ys,C_one1%ye
      do i=C_one1%xs,C_one1%xe
        hhl(i0, :, i, j) = hhl(i0, :, i, j) + global_maxheight - hhl(i0, C_one1%zs, i, j)
      enddo
    enddo

    call restoreVecPointer(vhhl , C_one1, hhl1d, hhl)

    call DMLocalToLocalBegin(C_one1%da, vhhl, INSERT_VALUES, vhhl,ierr) ;call CHKERR(ierr)
    call DMLocalToLocalEnd(C_one1%da, vhhl, INSERT_VALUES, vhhl,ierr) ;call CHKERR(ierr)

    call getVecPointer(vhhl , C_one1, hhl1d, hhl)

    call DMGetLocalVector(C_one1%da, vgrad_x, ierr) ;call CHKERR(ierr)
    call DMGetLocalVector(C_one1%da, vgrad_y, ierr) ;call CHKERR(ierr)

    call getVecPointer(vgrad_x , C_one1, grad_x1d, grad_x)
    call getVecPointer(vgrad_y , C_one1, grad_y1d, grad_y)

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

    call restoreVecPointer(vhhl , C_one1, hhl1d, hhl)
    call DMRestoreLocalVector(C_one1%da ,vhhl ,ierr);  call CHKERR(ierr)

    call restoreVecPointer(vgrad_x , C_one1, grad_x1d, grad_x)
    call restoreVecPointer(vgrad_y , C_one1, grad_y1d, grad_y)
  end subroutine

  !> @brief create PETSc matrix and inserts diagonal elements
  !> @details one important step for performance is to set preallocation of matrix structure.
  !>  \n  i.e. determine the number of local vs. remote number of entries in each row.
  !>  \n  DMDA actually provides a preallocation but this assumes that degrees of freedom on neighbouring boxes are fully connected to local ones.
  !>  \n  this does of course drastically overestimate non-zeros as we need only the streams that actually send radiation in the respective direction.
  !>  \n  at the moment preallocation routines determine nonzeros by manually checking bounadries --
  !>  \n  !todo we should really use some form of iterating through the entries as it is done in the matrix assembly routines and just flag the rows
  subroutine init_Matrix(A,C)!,prefix)
    Mat, allocatable, intent(inout) :: A
    type(t_coord) :: C
    !        character(len=*),optional :: prefix

    PetscInt,dimension(:),allocatable :: o_nnz,d_nnz!,dnz

    if(.not.allocated(A)) then
      allocate(A)
      call DMCreateMatrix(C%da, A, ierr) ;call CHKERR(ierr)
    endif

    call MatSeqAIJSetPreallocation(A, C%dof+i1, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)
    call MatMPIAIJSetPreallocation(A, C%dof+1, PETSC_NULL_INTEGER, C%dof, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)

    call mat_info(A)

    ! If matrix is resetted, keep nonzero pattern and allow to non-zero allocations -- those should not be many
    ! call MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr) ;call CHKERR(ierr)

    ! pressure mesh  may wiggle a bit and change atm%l1d -- keep the nonzeros flexible
    !call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr) ;call CHKERR(ierr)
    

    ! call MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr) ;call CHKERR(ierr) ! dont throw away the zero -- this completely destroys preallocation performance

    call MatSetFromOptions(A,ierr) ;call CHKERR(ierr)
    call MatSetUp(A,ierr) ;call CHKERR(ierr)

    call mat_info(A)
    call mat_set_diagonal(A,C)
  contains
    subroutine mat_set_diagonal(A,C,vdiag)
      Mat :: A
      type(t_coord),intent(in) :: C
      real(ireals),intent(in),optional :: vdiag

      PetscInt :: i,j,k,dof
      MatStencil :: row(4,1), col(4,1)
      PetscScalar :: v(1)

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

      if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients diagonally'
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

      if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients diagonally ... done'
    end subroutine
   
  end subroutine

  subroutine mat_info(A)
    Mat :: A
    MatInfo :: info(MAT_INFO_SIZE) 
    real(ireals) :: mal, nz_allocated, nz_used, nz_unneeded
    PetscInt :: m,n

    return !TODO see doxy details...
    call MatGetInfo(A,MAT_LOCAL,info,ierr) ;call CHKERR(ierr)
    mal = info(MAT_INFO_MALLOCS)
    nz_allocated = info(MAT_INFO_NZ_ALLOCATED)
    nz_used   = info(MAT_INFO_NZ_USED)
    nz_unneeded = info(MAT_INFO_NZ_UNNEEDED)

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
    class(t_solver) :: solver
    real(ireals), intent(in) :: albedo
    real(ireals),intent(in),dimension(:,:,:),optional :: local_kabs, local_ksca, local_g ! dimensions (Nz  , Nx, Ny)
    real(ireals),intent(in),dimension(:,:,:),optional :: local_planck                    ! dimensions (Nz+1, Nx, Ny) layer quantity plus surface layer
    real(ireals),intent(in),dimension(:,:),optional   :: local_albedo_2d                 ! dimensions (Nx, Ny)
    real(ireals) :: tau,kext,w0,g
    integer(iintegers) :: k,i,j

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
        print *,myid,'set_optical_properties :: found illegal value in local_optical properties! abort!'
        do k=C_one_atm%zs,C_one_atm%ze
          print *,myid,k,'local_kabs',local_kabs(k,:,:)
          print *,myid,k,'local_ksca',local_ksca(k,:,:)
        enddo
      endif
!    endif
!    if(ldebug) then
      if( (any([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g].lt.zero)) .or. (any(isnan([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g]))) ) then
        print *,myid,'set_optical_properties :: found illegal value in optical properties! abort!'
      endif
!    endif

    if(ldebug.and.myid.eq.0) then
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
      if(ldebug .and. myid.eq.0) then
        do k=C_one_atm%zs,C_one_atm%ze
          if(present(local_planck)) then
            print *,myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
              minval(atm%op(k,:,:)%kabs), minval(atm%op(k,:,:)%ksca), minval(atm%op(k,:,:)%g),&
              maxval(atm%op(k,:,:)%kabs), maxval(atm%op(k,:,:)%ksca), maxval(atm%op(k,:,:)%g),&
              '::',minval(atm%planck (k,:,:)),maxval(atm%planck        (k, :,:))
          else    
            print *,myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
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
        print *,myid,'set_optical_properties :: found illegal value in delta_scaled optical properties! abort!'
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

    if(ldebug .and. myid.eq.0) then
      do k=C_one_atm%zs,C_one_atm%ze
        if(present(local_planck)) then
          print *,myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
            minval(atm%op(k,:,:)%kabs), minval(atm%op(k,:,:)%ksca), minval(atm%op(k,:,:)%g),&
            maxval(atm%op(k,:,:)%kabs), maxval(atm%op(k,:,:)%ksca), maxval(atm%op(k,:,:)%g),&
            '::',minval(atm%planck (k,:,:)),maxval(atm%planck        (k, :,:))
        else    
          print *,myid,'Optical Properties:',k,'dz',atm%dz(k,C_one_atm%xs,C_one_atm%ys),atm%l1d(k,C_one_atm%xs,C_one_atm%ys),'k',&
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
end module

