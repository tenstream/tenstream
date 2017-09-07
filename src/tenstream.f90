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

!> \mainpage The Tenstream radiative transfer solver
!!  - The driving routines are: \subpage driving_routines
!!  - Important files to look into:
!!    *  \subpage optprop_parameters
!!  - At the moment there is no spectral integration parametrization included.

!> \page driving_routines Driving Routines
!! A typical use of the tenstream library would be to call the following functions in that order:
!!
!!  * m_tenstream::init_tenstream (needs only be called if any of the below changes)
!!    - setup the grid information for parallelization
!!    - setup the matrix structures
!!    - setup the LUT tables or neural networks for the coefficients
!!    - setup derivatives for sun directions
!!
!!  * m_tenstream::set_optical_properties
!!    - setup the optical properties for next calculation
!!    - call set_global_optical_properties if only rank 0 has the global info
!!
!!  * m_tenstream::solve_tenstream
!!    - setup the matrices and solve the equation system
!!
!!  * m_tenstream::tenstream_get_result
!!    - retrieve the result from solution vectors


module m_tenstream

#ifdef _XLF
  use ieee_arithmetic
#define isnan ieee_is_nan
#endif

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers,               &
    imp_comm, myid, numnodes, init_mpi_data_parameters, mpiint,  &
    zero, one, nil, i0, i1, i2, i3, i4, i5, i6, i7, i8, i10, pi, &
    default_str_len

  use m_twostream, only: delta_eddington_twostream
  use m_schwarzschild, only: schwarzschild
  use m_helper_functions, only: norm, rad2deg, deg2rad, &
    approx, rmse, delta_scale,                          &
    imp_bcast, mpi_logical_and, imp_allreduce_min,      &
    imp_allreduce_max, cumsum, inc, CHKERR

  use m_eddington, only : eddington_coeff_zdun
  use m_optprop_parameters, only : ldelta_scale
  use m_optprop, only : t_optprop_1_2,t_optprop_8_10
  use m_tenstream_options, only : read_commandline_options, ltwostr, luse_eddington, twostr_ratio, &
    options_max_solution_err, options_max_solution_time, ltwostr_only, luse_twostr_guess,        &
    options_phi, lforce_phi, options_theta, lforce_theta, &
    lcalc_nca, lskip_thermal, lschwarzschild, ltopography

  implicit none

  private
  public :: init_tenstream, set_angles,                                   &
    set_global_optical_properties, set_optical_properties,                &
    solve_tenstream, destroy_tenstream,                                   &
    getVecPointer,restoreVecPointer, get_mem_footprint,                   &
    tenstream_get_result, tenstream_get_result_toZero, need_new_solution, &
    t_coord,C_dir,C_diff,C_one,C_one1,C_one_atm, C_one_atm1

  PetscInt,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9

  logical,parameter :: ldebug=.False.
!  logical,parameter :: ldebug=.True.
  logical,parameter :: lcycle_dir=.True.
  logical,parameter :: lprealloc=.True.

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
    DM :: da                          ! The Domain Decomposition Object
    PetscMPIInt,allocatable :: neighbors(:) ! all 3d neighbours( (x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
  end type

  type(t_coord), allocatable, save :: C_dir, C_diff, C_one, C_one1, C_one_atm, C_one_atm1

  PetscErrorCode :: ierr

  type t_optprop
    real(ireals) :: kabs,ksca,g
  end type

  type t_atmosphere
    type(t_optprop) , allocatable , dimension(:,:,:) :: op
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
  type(t_atmosphere),allocatable,save :: atm


  type(t_optprop_1_2),save  :: OPP_1_2
  type(t_optprop_8_10),save :: OPP_8_10

  type t_sunangles
    real(ireals) :: symmetry_phi
    integer(iintegers) :: yinc,xinc
    real(ireals) :: theta, phi, costheta, sintheta
  end type
  type t_suninfo
    type(t_sunangles),allocatable :: angles(:,:,:) ! defined on DMDA grid
    logical :: luse_topography=.False.
  end type
  type(t_suninfo) :: sun

  PetscLogStage,save,allocatable :: logstage(:)

  Mat,allocatable :: Mdir,Mdiff

  Vec,allocatable,save :: incSolar,b

  KSP,save :: kspdir, kspdiff
  logical,save :: linit_kspdir=.False., linit_kspdiff=.False.

  logical,save :: ltenstream_is_initialized=.False.

  integer(iintegers),parameter :: minimal_dimension=3 ! this is the minimum number of gridpoints in x or y direction

  logical,parameter :: lenable_solutions_err_estimates=.True.  ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
  real(ireals),parameter :: time_debug_solutions=zero ! if enabled, we calculate new solutions but do not return the update solutions to the host model.(set to zero to disable)

  type t_state_container
    integer(iintegers) :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
    Vec :: edir,ediff,abso

    logical :: lset        = .False. ! initialized?
    logical :: lsolar_rad  = .False. ! direct radiation calculated?
    logical :: lchanged    = .True.  ! did the flux change recently? -- call restore_solution to bring it in a coherent state

    ! save state of solution vectors... they are either in [W](true) or [W/m**2](false)
    logical :: lintegrated_dir=.True. , lintegrated_diff=.True.

    !save error statistics
    real(ireals) :: time   (30) = -one
    real(ireals) :: maxnorm(30) = zero
    real(ireals) :: twonorm(30) = zero
    real(ireals),allocatable :: ksp_residual_history(:)
  end type
  type(t_state_container),save :: solutions(-1000:1000)

contains

  !> @brief Construct PETSC grid information for regular DMDA
  !> @details setup DMDA grid containers for direct, diffuse and absorption grid
  !>  \n and fill user context containers(t_coord) which include local as well as global array sizes
  !>  \n every mpi rank has to call this
  subroutine setup_grid(Nz_in,Nx,Ny,nxproc,nyproc, collapseindex)
    PetscInt,intent(in) :: Nz_in,Nx,Ny !< @param[in] local number of grid boxes -- in the vertical we have Nz boxes and Nz+1 levels
    integer(iintegers),optional :: nxproc(:), nyproc(:) ! size of local domains on each node
    integer(iintegers),optional,intent(in) :: collapseindex  !< @param[in] collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them

    DMBoundaryType :: bp=DM_BOUNDARY_PERIODIC, bn=DM_BOUNDARY_NONE, bg=DM_BOUNDARY_GHOSTED
    integer(iintegers) :: Nz

    Nz = Nz_in
    if(present(collapseindex)) Nz = Nz_in-collapseindex+i1

    if(myid.eq.0.and.ldebug) print *,myid,'Setting up the DMDA grid for ',Nz,Nx,Ny,'using ',numnodes,' nodes'

    if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_diff'
    call setup_dmda(C_diff, Nz+1,Nx,Ny, bp, i10)

    if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_dir'
    if(lcycle_dir) then
      call setup_dmda(C_dir, Nz+1,Nx,Ny, bp, i8)
    else
      call setup_dmda(C_dir, Nz+1,Nx,Ny, bg, i8)
    endif

    if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_one'
    call setup_dmda(C_one , Nz  , Nx,Ny,  bp, i1)
    call setup_dmda(C_one1, Nz+1, Nx,Ny,  bp, i1)

    if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA atm'
    call setup_dmda(C_one_atm , Nz_in  , Nx,Ny,  bp, i1)
    call setup_dmda(C_one_atm1, Nz_in+1, Nx,Ny,  bp, i1)

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

  !> @brief print information on PETSc Mat :: size and allocated rows
  !> @details TODO: currently broken -- need to clear what real_kind the returned values should have
  !> \n -- petsc doc says those should be double precision irrespective of petsc real type?? check!
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

    if(lprealloc) then
      ! Determine perfect preallocation
      if(numnodes.gt.1) then
        select case (C%dof)
        case(i3)
          call setup_dir3_preallocation(d_nnz,o_nnz,C)
        case(i8)
          call setup_dir8_preallocation(d_nnz,o_nnz,C)
        case(i10)
          call setup_diff10_preallocation(d_nnz,o_nnz,C)
        case default
          stop('Dont know which preallocation routine I shall call! - exiting...')
        end select

        call MatMPIAIJSetPreallocation(A, C%dof+1, d_nnz, C%dof, o_nnz, ierr) ;call CHKERR(ierr)

        deallocate(o_nnz)
        deallocate(d_nnz)

      endif
      call MatSeqAIJSetPreallocation(A, C%dof+i1, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)
    endif

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

    subroutine setup_dir3_preallocation(d_nnz,o_nnz,C)
      PetscInt,allocatable :: d_nnz(:)
      PetscInt,allocatable :: o_nnz(:)
      type(t_coord) :: C
      Vec :: gv_o_nnz,gv_d_nnz
      Vec :: lv_nnz
      PetscScalar,Pointer :: xl(:,:,:,:)=>null(),xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
      PetscScalar,Pointer :: xl1d(:)=>null(),xo1d(:)=>null(),xd1d(:)=>null()

      PetscInt :: i,j,k,dst, xinc, yinc

      call DMGetGlobalVector(C%da, gv_d_nnz,ierr) ;call CHKERR(ierr)
      call DMGetGlobalVector(C%da, gv_o_nnz,ierr) ;call CHKERR(ierr)
      call VecSet(gv_d_nnz, zero, ierr) ;call CHKERR(ierr)
      call VecSet(gv_o_nnz, zero, ierr) ;call CHKERR(ierr)

      call DMGetLocalVector(C%da, lv_nnz,ierr) ;call CHKERR(ierr)

      call getVecPointer(lv_nnz,C,xl1d,xl)

      xl=zero

      do j=C%ys,C%ye
        do i=C%xs,C%xe        
          do k=C%zs,C%ze-1
            if( atm%l1d(atmk(k),i,j) ) then
              do dst=i0,i3
                call inc( xl(dst, k+1, i,j), one )
              enddo
            else
              xinc = sun%angles(k,i,j)%xinc
              yinc = sun%angles(k,i,j)%yinc
              dst = i0 ;call inc( xl(dst , k+1 , i      , j      ) , one*C%dof )
              dst = i1 ;call inc( xl(dst , k+1 , i      , j      ) , one*C%dof )
              dst = i2 ;call inc( xl(dst , k+1 , i      , j      ) , one*C%dof )
              dst = i3 ;call inc( xl(dst , k+1 , i      , j      ) , one*C%dof )
              dst = i4 ;call inc( xl(dst , k   , i+xinc , j      ) , one*C%dof )
              dst = i5 ;call inc( xl(dst , k   , i+xinc , j      ) , one*C%dof )
              dst = i6 ;call inc( xl(dst , k   , i      , j+yinc ) , one*C%dof )
              dst = i7 ;call inc( xl(dst , k   , i      , j+yinc ) , one*C%dof )

            endif

          enddo 
        enddo 
      enddo
      call restoreVecPointer(lv_nnz,C,xl1d,xl)

      ! Now we have a local vector with ghosts that counts all entries
      ! We still have to seperate the diagonal and the off-diagonal part:
      ! For that, transmit the ghost parts with ADD_VALUES to _o_nnz
      call DMLocalToGlobalBegin(C%da,lv_nnz,ADD_VALUES,gv_o_nnz,ierr) ;call CHKERR(ierr)
      call DMLocalToGlobalEnd  (C%da,lv_nnz,ADD_VALUES,gv_o_nnz,ierr) ;call CHKERR(ierr)

      call getVecPointer(lv_nnz,C,xl1d,xl)
      call getVecPointer(gv_o_nnz,C,xo1d,xo)
      call getVecPointer(gv_d_nnz,C,xd1d,xd)

      do j=C%ys,C%ye
        do i=C%xs,C%xe        
          do k=C%zs,C%ze
            xo(:,k,i,j) = xo(:,k,i,j)-xl(:,k,i,j)
            xd(:,k,i,j) = -xo(:,k,i,j) + C%dof+i1
!            print *,myid,k,i,j,'off',int(xo(:,k,i,j)),'on',int(xd(:,k,i,j))
          enddo
        enddo
      enddo

      allocate(o_nnz(0:C%dof*C%zm*C%xm*C%ym-1))
      allocate(d_nnz(0:C%dof*C%zm*C%xm*C%ym-1))
      o_nnz=int(xo1d)
      d_nnz=int(xd1d)

      call restoreVecPointer(lv_nnz,C,xl1d,xl)
      call restoreVecPointer(gv_o_nnz,C,xo1d,xo)
      call restoreVecPointer(gv_d_nnz,C,xd1d,xd)

      call DMRestoreGlobalVector(C%da, gv_d_nnz,ierr) ;call CHKERR(ierr)
      call DMRestoreGlobalVector(C%da, gv_o_nnz,ierr) ;call CHKERR(ierr)
      call DMRestoreLocalVector (C%da, lv_nnz  ,ierr) ;call CHKERR(ierr)
    end subroutine 


    subroutine setup_dir8_preallocation(d_nnz,o_nnz,C)
      PetscInt,allocatable :: d_nnz(:)
      PetscInt,allocatable :: o_nnz(:)
      type(t_coord) :: C
      Vec :: v_o_nnz,v_d_nnz
      Vec :: g_o_nnz,g_d_nnz
      PetscScalar,Pointer :: xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
      PetscScalar,Pointer :: xo1d(:)=>null(),xd1d(:)=>null()

      PetscInt :: vsize,i,j,k, isrc, idst, xinc, yinc, src, dst, icnt

      logical :: llocal_src, llocal_dst

      MatStencil :: row(4,C%dof)  ,col(4,C%dof)

      if(myid.eq.0.and.ldebug) print *,myid,'building direct o_nnz for mat with',C%dof,'dof'
      call DMGetLocalVector(C%da,v_o_nnz,ierr) ;call CHKERR(ierr)
      call DMGetLocalVector(C%da,v_d_nnz,ierr) ;call CHKERR(ierr)

      call getVecPointer(v_o_nnz,C,xo1d,xo)
      call getVecPointer(v_d_nnz,C,xd1d,xd)

      
      xd = i0 
      xo = i0

      icnt = -1
      do j=C%ys,C%ye
          do i=C%xs,C%xe        
              do k=C%zs,C%ze-1

                  if( atm%l1d(atmk(k),i,j) ) then
                      do idst=i0,i3
                          call inc( xd(idst, k+1, i,j), one )
                      enddo
                  else
                      xinc = sun%angles(k,i,j)%xinc
                      yinc = sun%angles(k,i,j)%yinc

                      dst = 1 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
                      dst = 2 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
                      dst = 3 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
                      dst = 4 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
                      dst = 5 ; row(MatStencil_j,dst) = i+xinc   ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
                      dst = 6 ; row(MatStencil_j,dst) = i+xinc   ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
                      dst = 7 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j+yinc   ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
                      dst = 8 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j+yinc   ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid

                      src = 1 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
                      src = 2 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
                      src = 3 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
                      src = 4 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
                      src = 5 ; col(MatStencil_j,src) = i+1-xinc ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
                      src = 6 ; col(MatStencil_j,src) = i+1-xinc ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
                      src = 7 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j+1-yinc ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:
                      src = 8 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j+1-yinc ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:

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

      call restoreVecPointer(v_o_nnz,C,xo1d,xo)
      call restoreVecPointer(v_d_nnz,C,xd1d,xd)

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

      call getVecPointer(g_o_nnz,C,xo1d,xo)
      call getVecPointer(g_d_nnz,C,xd1d,xd)

      !call mpi_barrier(imp_comm, ierr)
      !icnt = -1
      !do j=C%ys,C%ye
      !    do i=C%xs,C%xe        
      !        do k=C%zs,C%ze
      !            do idst = 1,C%dof
      !                icnt = icnt+1
      !            enddo
      !            print *,myid,icnt,k,i,j,'::',int(xd(:, k,i,j)),'off',int(xo(:, k,i,j))
      !        enddo
      !    enddo
      !enddo
      !call mpi_barrier(imp_comm, ierr)

      call VecGetLocalSize(g_d_nnz,vsize,ierr) ;call CHKERR(ierr)
      allocate(o_nnz(0:vsize-1))
      allocate(d_nnz(0:vsize-1))

      o_nnz=int(xo1d)
      d_nnz=int(xd1d) + i1  ! +1 for diagonal entries

      call restoreVecPointer(g_o_nnz,C,xo1d,xo)
      call restoreVecPointer(g_d_nnz,C,xd1d,xd)

      call DMRestoreGlobalVector(C%da,g_o_nnz,ierr) ;call CHKERR(ierr)
      call DMRestoreGlobalVector(C%da,g_d_nnz,ierr) ;call CHKERR(ierr)

      if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
    end subroutine 
    subroutine setup_diff10_preallocation(d_nnz,o_nnz,C)
      PetscInt,allocatable :: d_nnz(:)
      PetscInt,allocatable :: o_nnz(:)
      type(t_coord) :: C

      Vec :: v_o_nnz,v_d_nnz
      Vec :: g_o_nnz,g_d_nnz
      PetscScalar,Pointer :: xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
      PetscScalar,Pointer :: xo1d(:)=>null(),xd1d(:)=>null()

      PetscInt :: vsize,i,j,k, isrc, idst, src, dst, icnt, src_id, dst_id

      MatStencil :: row(4,C%dof)  ,col(4,C%dof)

      if(myid.eq.0.and.ldebug) print *,myid,'building diffuse o_nnz for mat with',C%dof,'dof'
      call DMGetLocalVector(C%da,v_o_nnz,ierr) ;call CHKERR(ierr)
      call DMGetLocalVector(C%da,v_d_nnz,ierr) ;call CHKERR(ierr)

      call getVecPointer(v_o_nnz,C,xo1d,xo)
      call getVecPointer(v_d_nnz,C,xd1d,xd)

      
      xd = i0 
      xo = i0

      icnt = -1
      do j=C%ys,C%ye
        do i=C%xs,C%xe        
          do k=C%zs,C%ze-1

            if( atm%l1d(atmk(k),i,j) ) then
              call inc( xd(E_dn, k+1, i,j), 2*one )
              call inc( xd(E_up, k  , i,j), 2*one )

            else

              !           Sources
              !         ________E_dn___
              !        |               |
              !    E_ri|               |E_le_p
              !        |               |  
              !        |               |
              !    E_ri|               |E_le_m
              !        |_______________|  
              !                 E_up

              src = 1; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_dn
              src = 2; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k+1   ; col(MatStencil_c,src) = E_up 
              src = 3; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ri_m
              src = 4; col(MatStencil_j,src) = i+1  ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_le_m
              src = 5; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ri_p
              src = 6; col(MatStencil_j,src) = i+1  ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_le_p
              src = 7; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_fw_m
              src = 8; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j+1   ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ba_m
              src = 9; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_fw_p
              src =10; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j+1   ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ba_p

              !           Destinations
              !         ________E_up___
              !        |               |
              !    E_le|               |E_ri_p
              !        |               |  
              !        |               |
              !    E_le|               |E_ri_m
              !        |_______________|  
              !                 E_dn

              dst =  1; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_up  
              dst =  2; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k+1   ; row(MatStencil_c,dst) = E_dn  
              dst =  3; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_le_m
              dst =  4; row(MatStencil_j,dst) = i+1  ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ri_m
              dst =  5; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_le_p
              dst =  6; row(MatStencil_j,dst) = i+1  ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ri_p
              dst =  7; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ba_m
              dst =  8; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j+1   ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_fw_m
              dst =  9; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ba_p
              dst = 10; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j+1   ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_fw_p

              do idst = 1,C%dof
                icnt = icnt+1
                dst_id = myid
                if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 .and. row(MatStencil_j,idst).lt.C%xs ) dst_id = C%neighbors(10) ! have real neighbor west  and is not local entry
                if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 .and. row(MatStencil_j,idst).gt.C%xe ) dst_id = C%neighbors(16) ! have real neighbor east  and is not local entry
                if( C%neighbors( 4).ne.myid .and. C%neighbors( 4).ge.i0 .and. row(MatStencil_k,idst).lt.C%ys ) dst_id = C%neighbors( 4) ! have real neighbor south and is not local entry
                if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 .and. row(MatStencil_k,idst).gt.C%ye ) dst_id = C%neighbors(22) ! have real neighbor north and is not local entry

                do isrc = 1,C%dof
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
                  !if(myid.eq.0 .and. k.eq.C%zs) print *,myid,icnt,k,i,j,'::',idst,isrc,'::',src_id, dst_id,':: local?', src_id .eq. dst_id, '::', xd(row(4,idst),row(3,idst),row(2,idst),row(1,idst)), xo(row(4,idst),row(3,idst),row(2,idst),row(1,idst))

                enddo
              enddo

            endif ! atm_1d / 3d?

          enddo ! k

          ! Surface entries E_up 
          call inc( xd(E_up, C%ze, i,j), one )

        enddo
      enddo

      call restoreVecPointer(v_o_nnz,C,xo1d,xo)
      call restoreVecPointer(v_d_nnz,C,xd1d,xd)

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

      call getVecPointer(g_o_nnz,C,xo1d,xo)
      call getVecPointer(g_d_nnz,C,xd1d,xd)

      !call mpi_barrier(imp_comm, ierr)
      !icnt = 0
      !do j=C%ys,C%ye
      !    do i=C%xs,C%xe        
      !        do k=C%zs,C%ze
      !            print *,myid,icnt,k,i,j,'::',int(xd(:, k,i,j)),'off',int(xo(:, k,i,j))
      !            do idst = 1,C%dof
      !                icnt = icnt+1
      !            enddo
      !        enddo
      !    enddo
      !enddo
      !call mpi_barrier(imp_comm, ierr)

      call VecGetLocalSize(g_d_nnz,vsize,ierr) ;call CHKERR(ierr)
      allocate(o_nnz(0:vsize-1))
      allocate(d_nnz(0:vsize-1))

      o_nnz=int(xo1d)
      d_nnz=int(xd1d) + i1  ! +1 for diagonal entries

      call restoreVecPointer(g_o_nnz,C,xo1d,xo)
      call restoreVecPointer(g_d_nnz,C,xd1d,xd)

      call DMRestoreGlobalVector(C%da,g_o_nnz,ierr) ;call CHKERR(ierr)
      call DMRestoreGlobalVector(C%da,g_d_nnz,ierr) ;call CHKERR(ierr)

      if(myid.eq.0 .and. ldebug) print *,myid,'diffuse d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)

    end subroutine 

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
  end subroutine init_matrix

  !> @brief retrieve transport coefficients from optprop module
  !> @detail this may get the coeffs from a LUT or ANN or whatever and return diff2diff or dir2diff or dir2dir coeffs
  subroutine get_coeff(op,dz,ldir,coeff,lone_dimensional,angles)
    type(t_optprop),intent(in) :: op
    real(ireals),intent(in) :: dz
    logical,intent(in) :: ldir
    real(ireals),intent(out) :: coeff(:)

    real(ireals) :: aspect, tauz, w0

    logical,intent(in) :: lone_dimensional
    real(ireals),intent(in),optional :: angles(2)

    call PetscLogStagePush(logstage(7),ierr) ;call CHKERR(ierr)

    aspect = dz / atm%dx
    tauz = (op%kabs+op%ksca) * dz
    w0 = op%ksca / (op%kabs+op%ksca)

    if(lone_dimensional) then
      call OPP_1_2%get_coeff (aspect, tauz, w0, op%g,ldir,coeff,angles)
    else
      call OPP_8_10%get_coeff(aspect, tauz, w0, op%g,ldir,coeff,angles)
    endif

    call PetscLogStagePop(ierr) ;call CHKERR(ierr)
  end subroutine


  !> @brief compute gradient from dz3d
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information
  subroutine compute_gradient(atm, vgrad_x, vgrad_y)
    type(t_atmosphere),intent(in) :: atm
    Vec :: vgrad_x, vgrad_y

    Vec :: vhhl
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
          hhl(i0,k,i,j) = hhl(i0,k+1,i,j)+atm%dz(atmk(k),i,j)
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


  !> @brief setup topography information
  !> @details integrate dz3d from to top of atmosphere to bottom.
  !>   \n then build horizontal gradient of height information and
  !>   \n tweak the local sun angles to bend the rays.
  subroutine setup_topography(sun)
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

    call compute_gradient(atm, vgrad_x, vgrad_y)

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



  !> @brief set direction where sun stands
  !> @details save sun azimuth and zenith angle
  !>   \n sun azimuth is reduced to the range of [0,90] and the transmission of direct radiation is contributed for by a integer increment,
  !>   \n determining which neighbouring box is used in the horizontal direction
  subroutine setup_suninfo(phi0, theta0, sun, phi2d, theta2d)
    real(ireals),intent(in) :: phi0, theta0
    real(ireals), optional, intent(in) :: phi2d(:,:), theta2d(:,:)

    type(t_suninfo),intent(inout) :: sun
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

    if(sun%luse_topography) call setup_topography(sun)

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




  !> @brief build direct radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_dir_coeff(A,C)
    Mat :: A
    type(t_coord) :: C

    PetscInt :: i,j,k

    call PetscLogStagePush(logstage(2),ierr) ;call CHKERR(ierr)
    if(myid.eq.0.and.ldebug) print *,myid,'setup_direct_matrix ...'

    !      call MatZeroEntries(A, ierr) ;call CHKERR(ierr) !TODO necessary?
    !      call mat_set_diagonal(A,C)

    do j=C%ys,C%ye
      do i=C%xs,C%xe        
        do k=C%zs,C%ze-1

          if( atm%l1d(atmk(k),i,j) ) then
            call set_eddington_coeff(A,k, i,j)
          else
            call set_tenstream_coeff(C, A,k,i,j)
          endif

        enddo 
      enddo 
    enddo

    if(myid.eq.0.and.ldebug) print *,myid,'setup_direct_matrix done'

    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)  

    call PetscLogStagePop(ierr) ;call CHKERR(ierr)

  contains 
    subroutine set_tenstream_coeff(C,A,k,i,j)
      type(t_coord),intent(in) :: C
      Mat,intent(inout) :: A
      integer(iintegers),intent(in) :: i,j,k

      MatStencil :: row(4,C%dof)  ,col(4,C%dof)
      PetscScalar :: v(C%dof**2),norm

      integer(iintegers) :: dst,src, xinc, yinc

      xinc = sun%angles(k,i,j)%xinc
      yinc = sun%angles(k,i,j)%yinc

      dst = 1 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
      dst = 2 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
      dst = 3 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
      dst = 4 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k+1; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
      dst = 5 ; row(MatStencil_j,dst) = i+xinc   ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
      dst = 6 ; row(MatStencil_j,dst) = i+xinc   ; row(MatStencil_k,dst) = j        ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
      dst = 7 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j+yinc   ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
      dst = 8 ; row(MatStencil_j,dst) = i        ; row(MatStencil_k,dst) = j+yinc   ; row(MatStencil_i,dst) = k  ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid

      src = 1 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
      src = 2 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
      src = 3 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
      src = 4 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
      src = 5 ; col(MatStencil_j,src) = i+1-xinc ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
      src = 6 ; col(MatStencil_j,src) = i+1-xinc ; col(MatStencil_k,src) = j        ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
      src = 7 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j+1-yinc ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:
      src = 8 ; col(MatStencil_j,src) = i        ; col(MatStencil_k,src) = j+1-yinc ; col(MatStencil_i,src) = k  ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:

      call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j), .True., v, atm%l1d(atmk(k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta])

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

    subroutine set_eddington_coeff(A,k,i,j)
      Mat,intent(inout) :: A
      integer(iintegers),intent(in) :: i,j,k

      MatStencil :: row(4,1)  ,col(4,1)
      PetscScalar :: v(1)
      integer(iintegers) :: src

      if(luse_eddington) then
        v = atm%a33(atmk(k),i,j)
      else
        call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j),.True., v, atm%l1d(atmk(k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta] )
      endif

      col(MatStencil_j,i1) = i      ; col(MatStencil_k,i1) = j       ; col(MatStencil_i,i1) = k    
      row(MatStencil_j,i1) = i      ; row(MatStencil_k,i1) = j       ; row(MatStencil_i,i1) = k+1  

      do src=1,4
        col(MatStencil_c,i1) = src-i1 ! Source may be the upper/lower lid:
        row(MatStencil_c,i1) = src-i1 ! Define transmission towards the lower/upper lid
        call MatSetValuesStencil(A,i1, row,i1, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)
      enddo
    end subroutine

  end subroutine set_dir_coeff

  pure function atmk(k) ! return vertical index value for DMDA grid on atmosphere grid
      integer(iintegers) :: atmk
      integer(iintegers),intent(in) :: k
      atmk = k+atm%icollapse-1
  end function

  !> @brief set solar incoming radiation at Top_of_Atmosphere
  !> @details todo: in case we do not have periodic boundaries, we should shine light in from the side of the domain...
  subroutine setup_incSolar(incSolar,edirTOA)
    Vec :: incSolar
    real(ireals),intent(in) :: edirTOA

    PetscScalar,pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    PetscReal :: Az
    integer(iintegers) :: i,j
    Az = atm%dx*atm%dy

    call VecSet(incSolar,zero,ierr) ;call CHKERR(ierr)

    call getVecPointer(incSolar, C_dir, x1d, x4d)

    do j=C_dir%ys,C_dir%ye
      do i=C_dir%xs,C_dir%xe
        x4d(i0:i3,C_dir%zs,i,j) = edirTOA* Az * .25_ireals * sun%angles(C_dir%zs,i,j)%costheta
      enddo
    enddo

    call restoreVecPointer(incSolar,C_dir,x1d,x4d)

    if(myid.eq.0 .and. ldebug) print *,myid,'Setup of IncSolar done',edirTOA

  end subroutine

  !> @brief build diffuse radiation matrix
  !> @details will get the transfer coefficients for 1D and 3D Tenstream layers and input those into the matrix
  !>   \n get_coeff should provide coefficients in dst_order so that we can set  coeffs for a full block(i.e. all coeffs of one box)
  subroutine set_diff_coeff(A,C)
    Mat :: A
    type(t_coord) :: C

    PetscInt :: i,j,k

    call PetscLogStagePush(logstage(4),ierr) ;call CHKERR(ierr)

    if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients for diffuse Light'

    !      call MatZeroEntries(A, ierr) ;call CHKERR(ierr) !TODO necessary?
    !      call mat_set_diagonal(A,C)

    do j=C%ys,C%ye
      do i=C%xs,C%xe
        do k=C%zs,C%ze-1

          if( atm%l1d(atmk(k),i,j) ) then
            call set_eddington_coeff(A, k,i,j)
          else
            call set_tenstream_coeff(C, A, k,i,j, ierr); call CHKERR(ierr)
          endif

        enddo 
      enddo 
    enddo

    call set_albedo_coeff(C, A )

    if(myid.eq.0.and.ldebug) print *,myid,'setup_diffuse_matrix done'

    if(myid.eq.0.and.ldebug) print *,myid,'Final diffuse Matrix Assembly:'
    call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)
    call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;call CHKERR(ierr)

    call PetscLogStagePop(ierr) ;call CHKERR(ierr)

  contains
    subroutine set_tenstream_coeff(C,A,k,i,j,ierr)
      type(t_coord),intent(in) :: C
      Mat,intent(inout) :: A
      integer(iintegers),intent(in) :: k,i,j
      PetscErrorCode,intent(out) :: ierr

      MatStencil :: row(4,0:C%dof-1)  ,col(4,0:C%dof-1)
      PetscReal :: v(C%dof**2),norm

      integer(iintegers) :: dst,src

      ierr=0
      !           Sources
      !         ________E_dn___
      !        |               |
      !    E_ri|               |E_le_p
      !        |               |  
      !        |               |
      !    E_ri|               |E_le_m
      !        |_______________|  
      !                 E_up

      src = 0; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_dn
      src = 1; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k+1   ; col(MatStencil_c,src) = E_up 
      src = 2; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ri_m
      src = 3; col(MatStencil_j,src) = i+1  ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_le_m
      src = 4; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ri_p
      src = 5; col(MatStencil_j,src) = i+1  ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_le_p
      src = 6; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_fw_m
      src = 7; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j+1   ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ba_m
      src = 8; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_fw_p
      src = 9; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j+1   ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_ba_p

      !           Destinations
      !         ________E_up___
      !        |               |
      !    E_le|               |E_ri_p
      !        |               |  
      !        |               |
      !    E_le|               |E_ri_m
      !        |_______________|  
      !                 E_dn

      dst = 0; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_up  
      dst = 1; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k+1   ; row(MatStencil_c,dst) = E_dn  
      dst = 2; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_le_m
      dst = 3; row(MatStencil_j,dst) = i+1  ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ri_m
      dst = 4; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_le_p
      dst = 5; row(MatStencil_j,dst) = i+1  ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ri_p
      dst = 6; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ba_m
      dst = 7; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j+1   ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_fw_m
      dst = 8; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_ba_p
      dst = 9; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j+1   ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_fw_p

      call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j),.False., v, atm%l1d(atmk(k),i,j))

      call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)

      if(ldebug) then
        do src=1,C%dof
          norm = sum( v(src:C%dof**2:C%dof) )
          if( norm.gt.one+10._ireals*epsilon(one) ) then
            print *,'diffuse sum(dst==',dst,') gt one',norm
            stop 'omg.. shouldnt be happening'
            ierr = -5
            return
          endif
        enddo
      endif

    end subroutine

    subroutine set_eddington_coeff(A,k,i,j)
      Mat,intent(inout) :: A
      integer(iintegers),intent(in) :: k,i,j

      MatStencil :: row(4,2)  ,col(4,2)
      PetscReal :: v(4),twostr_coeff(4) ! v ==> a12,a11,a22,a21
      integer(iintegers) :: src,dst

      if(luse_eddington ) then
        v = [ atm%a12(atmk(k),i,j), atm%a11(atmk(k),i,j), atm%a22(atmk(k),i,j), atm%a21(atmk(k),i,j)]
      else
          call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j),.False., twostr_coeff, atm%l1d(atmk(k),i,j)) !twostr_coeff ==> a12,a11,a12,a11 !todo: check if this is the wrong order
          v = [ twostr_coeff(1), twostr_coeff(2) , twostr_coeff(2) , twostr_coeff(1) ]
      endif

      src = 1; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_dn
      src = 2; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k+1   ; col(MatStencil_c,src) = E_up 

      dst = 1; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_up  
      dst = 2; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k+1   ; row(MatStencil_c,dst) = E_dn  

      call MatSetValuesStencil(A,i2, row, i2, col , -v ,INSERT_VALUES,ierr) ;call CHKERR(ierr)


    end subroutine

    !> @brief insert lower boundary condition, i.e. diffuse reflection of downward radiation
    subroutine set_albedo_coeff(C,A)
      type(t_coord) :: C
      Mat,intent(inout) :: A

      MatStencil :: row(4,1)  ,col(4,1)
      integer(iintegers) :: i,j

      ! Set surface albedo values
      col(MatStencil_i,i1) = C%ze
      row(MatStencil_i,i1) = C%ze

      col(MatStencil_c,i1) = E_dn
      row(MatStencil_c,i1) = E_up

      do j=C%ys,C%ye
        row(MatStencil_k,i1) = j 
        col(MatStencil_k,i1) = j 

        do i=C%xs,C%xe
          row(MatStencil_j,i1) = i          
          col(MatStencil_j,i1) = i        

          call MatSetValuesStencil(A,i1, row, i1, col , [-atm%albedo(i,j)] ,INSERT_VALUES,ierr) ;call CHKERR(ierr)

        enddo
      enddo

    end subroutine
  end subroutine

  !> @brief setup source term for diffuse radiation
  !> @details this is either direct radiation scattered into one of the diffuse coeffs:
  !> \n direct source term is
  !> \n   direct radiation times the dir2diff coeffs
  !> \n or it may be that we have a source term due to thermal emission --
  !> \n   to determine emissivity of box, we use the forward transport coefficients backwards
  !> \n   a la: transmissivity $T = \sum(coeffs)$ and therefore emissivity $E = 1 - T$
  subroutine setup_b(solution,b)
    type(t_state_container) :: solution
    Vec :: local_b,b

    PetscScalar,pointer,dimension(:,:,:,:) :: xsrc=>null()
    PetscScalar,pointer,dimension(:) :: xsrc1d=>null()
    PetscInt :: k,i,j,src,dst

    if(myid.eq.0.and.ldebug) print *,'src Vector Assembly...'
    call PetscLogStagePush(logstage(6),ierr) ;call CHKERR(ierr)

    call DMGetLocalVector(C_diff%da,local_b,ierr) ;call CHKERR(ierr)
    call VecSet(local_b,zero,ierr) ;call CHKERR(ierr)

    call getVecPointer(local_b, C_diff, xsrc1d, xsrc)

    if(solution%lsolar_rad) &
      call set_solar_source(solution%edir)

    if(allocated(atm%planck) ) &
      call set_thermal_source()

    if(myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients ...done'

    call restoreVecPointer(local_b, C_diff, xsrc1d , xsrc )

    call VecSet(b,zero,ierr) ;call CHKERR(ierr) ! reset global Vec

    call DMLocalToGlobalBegin(C_diff%da, local_b, ADD_VALUES, b,ierr) ;call CHKERR(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
    call DMLocalToGlobalEnd  (C_diff%da, local_b, ADD_VALUES, b,ierr) ;call CHKERR(ierr)

    call DMRestoreLocalVector(C_diff%da,local_b,ierr) ;call CHKERR(ierr)

    call PetscLogStagePop(ierr) ;call CHKERR(ierr)
    if(myid.eq.0.and.ldebug) print *,'src Vector Assembly done'

  contains
    subroutine set_thermal_source()
      real(ireals) :: Ax,Ay,Az,b0 !,c1,c2,c3,b1,dtau
      real(ireals) :: diff2diff1d(4)
      real(ireals) :: diff2diff(C_diff%dof**2),v(C_diff%dof**2)

      if(myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector ... setting thermal source terms', minval(atm%planck), maxval(atm%planck)
      Az = atm%dx*atm%dy

      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          do k=C_diff%zs,C_diff%ze-1

            if( atm%l1d(atmk(k),i,j) ) then

              if(luse_eddington ) then

                b0 = atm%planck(atmk(k),i,j) * (one-atm%a11(atmk(k),i,j)-atm%a12(atmk(k),i,j))
                xsrc(E_up   ,k  ,i,j) = xsrc(E_up   ,k  ,i,j) + b0 *Az*pi
                xsrc(E_dn   ,k+1,i,j) = xsrc(E_dn   ,k+1,i,j) + b0 *Az*pi

              else

                call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j),.False., diff2diff1d, atm%l1d(atmk(k),i,j))

                b0 = atm%planck(atmk(k),i,j) * pi
                xsrc(E_up   ,k  ,i,j) = xsrc(E_up   ,k  ,i,j) +  b0  *(one-diff2diff1d(1)-diff2diff1d(2) ) *Az
                xsrc(E_dn   ,k+1,i,j) = xsrc(E_dn   ,k+1,i,j) +  b0  *(one-diff2diff1d(1)-diff2diff1d(2) ) *Az

              endif

            else ! Tenstream source terms
              Ax = atm%dy*atm%dz(atmk(k),i,j)
              Ay = atm%dx*atm%dz(atmk(k),i,j)

              call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j),.False., diff2diff, atm%l1d(atmk(k),i,j) )
              ! reorder from destination ordering to src ordering
              do src=1,C_diff%dof
                v(src:C_diff%dof**2:C_diff%dof) = diff2diff( i1+(src-i1)*C_diff%dof : src*C_diff%dof )
              enddo
              b0 = atm%planck(atmk(k),i,j) * pi
              xsrc(E_up   , k   , i   , j   ) = xsrc(E_up   , k   , i   , j   ) +  b0  *(one-sum( v( E_up  *C_diff%dof+i1 : E_up  *C_diff%dof+C_diff%dof )  )  ) *Az
              xsrc(E_dn   , k+1 , i   , j   ) = xsrc(E_dn   , k+1 , i   , j   ) +  b0  *(one-sum( v( E_dn  *C_diff%dof+i1 : E_dn  *C_diff%dof+C_diff%dof )  )  ) *Az
              xsrc(E_le_m , k   , i   , j   ) = xsrc(E_le_m , k   , i   , j   ) +  b0  *(one-sum( v( E_le_m*C_diff%dof+i1 : E_le_m*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
              xsrc(E_le_p , k   , i   , j   ) = xsrc(E_le_p , k   , i   , j   ) +  b0  *(one-sum( v( E_le_p*C_diff%dof+i1 : E_le_p*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
              xsrc(E_ri_m , k   , i+1 , j   ) = xsrc(E_ri_m , k   , i+1 , j   ) +  b0  *(one-sum( v( E_ri_m*C_diff%dof+i1 : E_ri_m*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
              xsrc(E_ri_p , k   , i+1 , j   ) = xsrc(E_ri_p , k   , i+1 , j   ) +  b0  *(one-sum( v( E_ri_p*C_diff%dof+i1 : E_ri_p*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
              xsrc(E_ba_m , k   , i   , j   ) = xsrc(E_ba_m , k   , i   , j   ) +  b0  *(one-sum( v( E_ba_m*C_diff%dof+i1 : E_ba_m*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
              xsrc(E_ba_p , k   , i   , j   ) = xsrc(E_ba_p , k   , i   , j   ) +  b0  *(one-sum( v( E_ba_p*C_diff%dof+i1 : E_ba_p*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
              xsrc(E_fw_m , k   , i   , j+1 ) = xsrc(E_fw_m , k   , i   , j+1 ) +  b0  *(one-sum( v( E_fw_m*C_diff%dof+i1 : E_fw_m*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
              xsrc(E_fw_p , k   , i   , j+1 ) = xsrc(E_fw_p , k   , i   , j+1 ) +  b0  *(one-sum( v( E_fw_p*C_diff%dof+i1 : E_fw_p*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
            endif ! 1D or Tenstream?

          enddo
        enddo

      enddo ! k

      ! Thermal emission at surface
      k = C_diff%ze
      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          xsrc(E_up   ,k,i,j) = xsrc(E_up   ,k,i,j) + atm%planck(atmk(k),i,j)*Az *(one-atm%albedo(i,j))*pi
        enddo
      enddo
    end subroutine
    subroutine set_solar_source(edir)
      Vec :: edir
      real(ireals) :: twostr_coeff(2)
      real(ireals) :: dir2diff(C_dir%dof*C_diff%dof), solrad(C_dir%dof)

      Vec :: ledir
      PetscScalar,pointer,dimension(:,:,:,:) :: xedir=>null()
      PetscScalar,pointer,dimension(:)       :: xedir1d=>null()

      logical :: lsun_east,lsun_north


      ! Copy ghosted values for direct vec
      call DMGetLocalVector(C_dir%da ,ledir ,ierr)                      ; call CHKERR(ierr)
      call VecSet(ledir ,zero,ierr)                                     ; call CHKERR(ierr)
      call DMGlobalToLocalBegin(C_dir%da ,edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C_dir%da ,edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)

      call getVecPointer(ledir, C_dir, xedir1d, xedir)

      if(myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector .. setting solar source',sum(xedir(0:3,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye))/size(xedir(0:3,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye))
      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          do k=C_diff%zs,C_diff%ze-1

            if( any (xedir(:,k,i,j) .gt. epsilon(one)) ) then
              if( atm%l1d(atmk(k),i,j) ) then
                dir2diff = zero
                if(luse_eddington ) then
                  ! Only transport the 4 tiles from dir0 to the Eup and Edn
                  do src=1,4
                    dir2diff(E_up  +i1+(src-1)*C_diff%dof) = atm%a13(atmk(k),i,j)
                    dir2diff(E_dn  +i1+(src-1)*C_diff%dof) = atm%a23(atmk(k),i,j)
                  enddo

                else
                  call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j),.False., twostr_coeff, atm%l1d(atmk(k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta])
                  do src=1,4
                    dir2diff(E_up  +i1+(src-1)*C_diff%dof) = twostr_coeff(1)
                    dir2diff(E_dn  +i1+(src-1)*C_diff%dof) = twostr_coeff(2)
                  enddo
                endif

                do src=1,C_dir%dof-4
                  xsrc(E_up   ,k,i,j)   = xsrc(E_up   ,k,i,j)   +  xedir(src-1,k,i,j)*dir2diff(E_up  +i1+(src-1)*C_diff%dof)
                  xsrc(E_dn   ,k+1,i,j) = xsrc(E_dn   ,k+1,i,j) +  xedir(src-1,k,i,j)*dir2diff(E_dn  +i1+(src-1)*C_diff%dof)
                enddo

              else ! Tenstream source terms
                lsun_east  = sun%angles(k,i,j)%xinc.eq.i0
                lsun_north = sun%angles(k,i,j)%yinc.eq.i0

                call get_coeff(atm%op(atmk(k),i,j), atm%dz(atmk(k),i,j),.False., dir2diff,  atm%l1d(atmk(k),i,j), [sun%angles(k,i,j)%symmetry_phi, sun%angles(k,i,j)%theta] )

                do src=1,C_dir%dof
                  select case(src)
                  case (1:4)
                    solrad = xedir( : , k , i , j )
                  case (5:6)
                    solrad = xedir( : , k , i+i1-sun%angles(k,i,j)%xinc , j )
                  case (7:8)
                    solrad = xedir( : , k , i , j+i1-sun%angles(k,i,j)%yinc )
                  case default
                    stop 'invalid dof for solar source term'
                  end select
                  xsrc(E_up   , k   , i   , j   ) = xsrc(E_up   , k   , i   , j   ) +  solrad(src) *dir2diff(E_up*C_dir%dof + src )
                  xsrc(E_dn   , k+1 , i   , j   ) = xsrc(E_dn   , k+1 , i   , j   ) +  solrad(src) *dir2diff(E_dn*C_dir%dof + src )

                  if(lsun_east) then ! if sun shines from right to left, switch ''reflection'' and ''transmission'' coefficients.
                    xsrc(E_le_m , k   , i   , j   ) = xsrc(E_le_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_ri_m*C_dir%dof + src)
                    xsrc(E_le_p , k   , i   , j   ) = xsrc(E_le_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_ri_p*C_dir%dof + src)
                    xsrc(E_ri_m , k   , i+1 , j   ) = xsrc(E_ri_m , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_le_m*C_dir%dof + src)
                    xsrc(E_ri_p , k   , i+1 , j   ) = xsrc(E_ri_p , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_le_p*C_dir%dof + src)
                  else
                    xsrc(E_le_m , k   , i   , j   ) = xsrc(E_le_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_le_m*C_dir%dof + src)
                    xsrc(E_le_p , k   , i   , j   ) = xsrc(E_le_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_le_p*C_dir%dof + src)
                    xsrc(E_ri_m , k   , i+1 , j   ) = xsrc(E_ri_m , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_ri_m*C_dir%dof + src)
                    xsrc(E_ri_p , k   , i+1 , j   ) = xsrc(E_ri_p , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_ri_p*C_dir%dof + src)
                  endif
                  if(lsun_north) then ! likewise if sun shines from forward to backwards,
                    xsrc(E_ba_m , k   , i   , j   ) = xsrc(E_ba_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_fw_m*C_dir%dof + src)
                    xsrc(E_ba_p , k   , i   , j   ) = xsrc(E_ba_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_fw_p*C_dir%dof + src)
                    xsrc(E_fw_m , k   , i   , j+1 ) = xsrc(E_fw_m , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_ba_m*C_dir%dof + src)
                    xsrc(E_fw_p , k   , i   , j+1 ) = xsrc(E_fw_p , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_ba_p*C_dir%dof + src)
                  else
                    xsrc(E_ba_m , k   , i   , j   ) = xsrc(E_ba_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_ba_m*C_dir%dof + src)
                    xsrc(E_ba_p , k   , i   , j   ) = xsrc(E_ba_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_ba_p*C_dir%dof + src)
                    xsrc(E_fw_m , k   , i   , j+1 ) = xsrc(E_fw_m , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_fw_m*C_dir%dof + src)
                    xsrc(E_fw_p , k   , i   , j+1 ) = xsrc(E_fw_p , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_fw_p*C_dir%dof + src)
                  endif

                  if(ldebug) then
                    do dst=0,C_diff%dof-1
                      if(sum(dir2diff( dst*C_dir%dof : dst*(C_dir%dof+1)-1 )) .gt. one .or. &
                        sum(dir2diff( dst*C_dir%dof : dst*(C_dir%dof+1)-1 )) .lt. zero   ) &
                        print *,'DEBUG Found dir2diff gt one:',src,'::',sum(dir2diff( dst*C_dir%dof : dst*(C_dir%dof+1)-1  ) ),'::',dir2diff(dst*C_dir%dof : dst*(C_dir%dof+1)-1) ,'   :::::::     ', dir2diff
                    enddo
                  endif
                enddo

              endif ! 1D or Tenstream?
            endif ! if solar

          enddo
        enddo
      enddo

      ! Ground Albedo reflecting direct radiation, the diffuse part is considered by the solver(Matrix)
      k = C_diff%ze
      do j=C_diff%ys,C_diff%ye
        do i=C_diff%xs,C_diff%xe
          xsrc(E_up   ,k,i,j) = sum(xedir(i0:i3,k,i,j))*atm%albedo(i,j)
        enddo
      enddo

      call restoreVecPointer(ledir ,C_dir ,xedir1d ,xedir )
      call DMRestoreLocalVector(C_dir%da ,ledir ,ierr)       ; call CHKERR(ierr)
    end subroutine
  end subroutine setup_b

  !> @brief Compute flux divergence, i.e. absorption
  !> @details from gauss divergence theorem, the divergence in the volume is the integral of the flux through the surface
  !> \n we therefore sum up the incoming and outgoing fluxes to compute the divergence
  subroutine calc_flx_div(solution)
    type(t_state_container) :: solution
    PetscReal,pointer,dimension(:,:,:,:) :: xediff=>null(),xedir=>null(),xabso=>null()
    PetscReal,pointer,dimension(:) :: xediff1d=>null(),xedir1d=>null(),xabso1d=>null()
    PetscInt :: i,j,k, xinc,yinc
    Vec :: ledir,lediff ! local copies of vectors, including ghosts
    PetscReal :: div(3),div2(13)
    PetscReal :: Volume,Az
    logical :: lhave_no_3d_layer

    if(solution%lsolar_rad .and. (solution%lintegrated_dir .eqv..False.)) stop 'tried calculating absorption but dir  vector was in [W/m**2], not in [W], scale first!'
    if(                          (solution%lintegrated_diff.eqv..False.)) stop 'tried calculating absorption but diff vector was in [W/m**2], not in [W], scale first!'

    if( (solution%lsolar_rad.eqv..False.) .and. lcalc_nca ) then ! if we should calculate NCA (Klinger), we can just return afterwards
      call scale_flx(solution, lWm2_to_W=.False.)
      call nca_wrapper(solution%ediff, solution%abso)
      call scale_flx(solution, lWm2_to_W=.True.)
      return
    endif

    if(myid.eq.0.and.ldebug) print *,'Calculating flux divergence',solution%lsolar_rad.eqv..False.,lcalc_nca
    call VecSet(solution%abso,zero,ierr) ;call CHKERR(ierr)

    ! if there are no 3D layers globally, we should skip the ghost value copying....
    lhave_no_3d_layer = mpi_logical_and(imp_comm, all(atm%l1d.eqv..True.))
    if(lhave_no_3d_layer) then

      if(solution%lsolar_rad) call getVecPointer(solution%edir, C_dir ,xedir1d ,xedir )

      call getVecPointer(solution%ediff, C_diff, xediff1d, xediff)
      call getVecPointer(solution%abso, C_one, xabso1d, xabso)

      ! calculate absorption by flux divergence
      Az = atm%dx * atm%dy

      do j=C_one%ys,C_one%ye         
        do i=C_one%xs,C_one%xe      
          do k=C_one%zs,C_one%ze
            Volume = Az     * atm%dz(atmk(k),i,j)
            ! Divergence    =                       Incoming                -       Outgoing
            if(solution%lsolar_rad) then
              div(1) = sum( xedir(i0:i3, k, i, j )  - xedir(i0:i3 , k+i1 , i, j ) ) 
            else 
              div(1) = zero
            endif

            div(2) = ( xediff(E_up  ,k+1,i  ,j  )  - xediff(E_up  ,k  ,i  ,j  )  ) 
            div(3) = ( xediff(E_dn  ,k  ,i  ,j  )  - xediff(E_dn  ,k+1,i  ,j  )  ) 

            xabso(i0,k,i,j) = sum(div) / Volume
          enddo                             
        enddo                             
      enddo   

      if(solution%lsolar_rad) call restoreVecPointer(solution%edir, C_dir, xedir1d, xedir )

      call restoreVecPointer(solution%ediff, C_diff,xediff1d,xediff)
      call restoreVecPointer(solution%abso , C_one ,xabso1d ,xabso )

      return
    endif

    if(solution%lsolar_rad) then
      ! Copy ghosted values for direct vec
      call DMGetLocalVector(C_dir%da ,ledir ,ierr)                      ; call CHKERR(ierr)
      call VecSet(ledir ,zero,ierr)                                     ; call CHKERR(ierr)
      call DMGlobalToLocalBegin(C_dir%da ,solution%edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)
      call DMGlobalToLocalEnd  (C_dir%da ,solution%edir ,ADD_VALUES,ledir ,ierr) ; call CHKERR(ierr)
      call getVecPointer(ledir, C_dir ,xedir1d ,xedir )
    endif

    ! Copy ghosted values for diffuse vec
    call DMGetLocalVector(C_diff%da,lediff,ierr)                      ; call CHKERR(ierr)
    call VecSet(lediff,zero,ierr)                                     ; call CHKERR(ierr)
    call DMGlobalToLocalBegin(C_diff%da,solution%ediff,ADD_VALUES,lediff,ierr) ; call CHKERR(ierr)
    call DMGlobalToLocalEnd  (C_diff%da,solution%ediff,ADD_VALUES,lediff,ierr) ; call CHKERR(ierr)
    call getVecPointer(lediff, C_diff, xediff1d, xediff)

    call getVecPointer(solution%abso, C_one, xabso1d, xabso)

    ! calculate absorption by flux divergence
    Az = atm%dx * atm%dy

    do j=C_one%ys,C_one%ye         
      do i=C_one%xs,C_one%xe      
        do k=C_one%zs,C_one%ze

          Volume = Az     * atm%dz(atmk(k),i,j)

          if(atm%l1d(atmk(k),i,j)) then ! one dimensional i.e. twostream
            ! Divergence    =                       Incoming                -       Outgoing
            if(solution%lsolar_rad) then
              div(1) = sum( xedir(i0:i3, k, i, j )  - xedir(i0:i3 , k+i1 , i, j ) ) 
            else 
              div(1) = zero
            endif

            div(2) = ( xediff(E_up  ,k+1,i  ,j  )  - xediff(E_up  ,k  ,i  ,j  )  ) 
            div(3) = ( xediff(E_dn  ,k  ,i  ,j  )  - xediff(E_dn  ,k+1,i  ,j  )  ) 

            xabso(i0,k,i,j) = sum(div) / Volume
            !              if(xabso(i0,k,i,j).lt.zero) print *,'1D abso<0 :: ',i,j,k,'::',xabso(i0,k,i,j),'::',div
          else

            ! Divergence    =                 Incoming                        -                   Outgoing
            if(solution%lsolar_rad) then
              xinc = sun%angles(k,i,j)%xinc
              yinc = sun%angles(k,i,j)%yinc
              div2( 1) = sum( xedir(i0:i3 , k, i         , j         )  - xedir(i0:i3 , k+i1 , i      , j       ) ) 
              div2( 2) = sum( xedir(i4:i5 , k, i+i1-xinc , j         )  - xedir(i4:i5 , k    , i+xinc , j       ) ) 
              div2( 3) = sum( xedir(i6:i7 , k, i         , j+i1-yinc )  - xedir(i6:i7 , k    , i      , j+yinc  ) ) 
            else 
              div2(1:3) = zero
            endif

            div2( 4) = ( xediff(E_up  ,k+1,i  ,j  )  - xediff(E_up  ,k  ,i  ,j  )  ) 
            div2( 5) = ( xediff(E_dn  ,k  ,i  ,j  )  - xediff(E_dn  ,k+1,i  ,j  )  ) 
            div2( 6) = ( xediff(E_le_m,k  ,i+1,j  )  - xediff(E_le_m,k  ,i  ,j  )  ) 
            div2( 7) = ( xediff(E_le_p,k  ,i+1,j  )  - xediff(E_le_p,k  ,i  ,j  )  ) 
            div2( 8) = ( xediff(E_ri_m,k  ,i  ,j  )  - xediff(E_ri_m,k  ,i+1,j  )  ) 
            div2( 9) = ( xediff(E_ri_p,k  ,i  ,j  )  - xediff(E_ri_p,k  ,i+1,j  )  ) 
            div2(10) = ( xediff(E_ba_m,k  ,i  ,j+1)  - xediff(E_ba_m,k  ,i  ,j  )  ) 
            div2(11) = ( xediff(E_ba_p,k  ,i  ,j+1)  - xediff(E_ba_p,k  ,i  ,j  )  ) 
            div2(12) = ( xediff(E_fw_m,k  ,i  ,j  )  - xediff(E_fw_m,k  ,i  ,j+1)  ) 
            div2(13) = ( xediff(E_fw_p,k  ,i  ,j  )  - xediff(E_fw_p,k  ,i  ,j+1)  ) 

            xabso(i0,k,i,j) = sum(div2) / Volume
            if(ldebug) then
              if( isnan(xabso(i0,k,i,j)) ) print *,'nan in flxdiv',k,i,j,'::',xabso(i0,k,i,j),Volume,'::',div2
            endif
          endif
        enddo                             
      enddo                             
    enddo   

    if(solution%lsolar_rad) then
      call restoreVecPointer(ledir          ,C_dir ,xedir1d ,xedir )
      call DMRestoreLocalVector(C_dir%da ,ledir ,ierr) ; call CHKERR(ierr)
    endif

    call restoreVecPointer(lediff         ,C_diff,xediff1d,xediff)
    call DMRestoreLocalVector(C_diff%da,lediff,ierr) ; call CHKERR(ierr)

    call restoreVecPointer(solution%abso  ,C_one ,xabso1d ,xabso )
  end subroutine


  !> @brief call PETSc Krylov Subspace Solver
  !> @details solve with ksp and save residual history of solver
  !> \n -- this may be handy later to decide next time if we have to calculate radiation again
  !> \n if we did not get convergence, we try again with standard GMRES and a resetted(zero) initial guess -- if that doesnt help, we got a problem!
  subroutine solve(ksp,b,x,solution_uid)
    KSP :: ksp
    Vec:: b
    Vec:: x
    integer(iintegers),optional,intent(in) :: solution_uid

    KSPConvergedReason :: reason
    PetscInt :: iter

    KSPType :: old_ksp_type

    if(myid.eq.0.and.ldebug) print *,'Solving Matrix'

    if(present(solution_uid)) then
      if(.not.allocated( solutions(solution_uid)%ksp_residual_history) ) allocate(solutions(solution_uid)%ksp_residual_history(100) )
      solutions(solution_uid)%ksp_residual_history = -1
      call KSPSetResidualHistory(ksp, solutions(solution_uid)%ksp_residual_history, 100_iintegers, .True.,ierr)
    endif

    call KSPSolve(ksp,b,x,ierr) ;call CHKERR(ierr)
    call KSPGetIterationNumber(ksp,iter,ierr) ;call CHKERR(ierr)
    call KSPGetConvergedReason(ksp,reason,ierr) ;call CHKERR(ierr)

    ! if(reason.eq.KSP_DIVERGED_ITS) then
    !   if(myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
    !   return
    ! endif

    if(reason.le.0) then
      if(myid.eq.0.and.ldebug) print *,myid,'Resetted initial guess to zero and try again with gmres:'
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
      if(myid.eq.0.and.ldebug) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
    endif

    if(reason.le.0) then
      if(myid.eq.0) print *,'***** SOLVER did NOT converge :( ********',reason
      call exit()
    endif
  end subroutine

  !> @brief initialize PETSc Krylov Subspace Solver
  !> @details default KSP solver is a FBCGS with BJCAOBI // ILU(1)
  !> \n -- the default does however not scale well -- and we configure petsc solvers per commandline anyway
  !> \n -- see documentation for details on how to do so
  subroutine setup_ksp(ksp,C,A,linit, prefix)
    KSP :: ksp
    type(t_coord) :: C
    Mat :: A
    PC  :: prec
    logical :: linit

    MatNullSpace :: nullspace
    Vec :: nullvecs(0)
    character(len=*),optional :: prefix

    PetscReal,parameter :: rtol=sqrt(epsilon(rtol))*10, rel_atol=1e-4_ireals
    PetscInt,parameter  :: maxiter=1000

    PetscInt,parameter :: ilu_default_levels=1
    PetscInt :: pcbjac_n_local, pcbjac_iglob ! number of local ksp contexts and index in global ksp-table
    KSP :: pcbjac_ksps(100) ! we dont have a good metric to check how many jacobi blocks there should be. in earlier versions of petsc we could let it decide. now this seems a problem at the moment. Lets just use a big number here... is only a list of pointers anyway...
    PC  :: pcbjac_sub_pc
    integer(iintegers) :: isub

    PetscReal :: atol

    logical,parameter :: lset_geometry=.True.  ! this may be necessary in order to use geometric multigrid
!    logical,parameter :: lset_geometry=.False.  ! this may be necessary in order to use geometric multigrid
    logical,parameter :: lset_nullspace=.True. ! set constant nullspace?
    !logical,parameter :: lset_nullspace=.False. ! set constant nullspace?

    if(linit) return
    call PetscLogStagePush(logstage(9),ierr) ;call CHKERR(ierr)

    call imp_allreduce_min(imp_comm, rel_atol*(C%dof*C%glob_xm*C%glob_ym*C%glob_zm) * count(.not.atm%l1d)/(one*size(atm%l1d)), atol)
    atol = max(1e-8_ireals, atol)

    if(myid.eq.0.and.ldebug) &
      print *,'Setup KSP -- tolerances:',rtol,atol,'::',rel_atol,(C%dof*C%glob_xm*C%glob_ym*C%glob_zm),count(.not.atm%l1d),one*size(atm%l1d)

    call KSPCreate(imp_comm,ksp,ierr) ;call CHKERR(ierr)
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

    if(numnodes.eq.0) then
      call PCFactorSetLevels(prec,ilu_default_levels,ierr);call CHKERR(ierr)
    else
      call PCBJacobiGetSubKSP(prec,pcbjac_n_local,pcbjac_iglob,pcbjac_ksps,ierr);call CHKERR(ierr)
      !if(.not.allocated(pcbjac_ksps)) allocate(pcbjac_ksps(pcbjac_n_local))
      !call PCBJacobiGetSubKSP(prec,pcbjac_n_local,pcbjac_iglob,pcbjac_ksps,ierr);call CHKERR(ierr)

      do isub=1,pcbjac_n_local
        call KSPSetType(pcbjac_ksps(isub) ,KSPPREONLY,ierr)              ;call CHKERR(ierr)
        call KSPGetPC  (pcbjac_ksps(isub), pcbjac_sub_pc,ierr)        ;call CHKERR(ierr)
        call PCSetType (pcbjac_sub_pc, PCILU,ierr)                    ;call CHKERR(ierr)
        call PCFactorSetLevels(pcbjac_sub_pc,ilu_default_levels,ierr) ;call CHKERR(ierr)
      enddo

    endif

    if(lset_geometry) call set_coordinates(C,ierr);call CHKERR(ierr)

    if(lset_nullspace) then
      call MatNullSpaceCreate( imp_comm, PETSC_TRUE, i0, nullvecs, nullspace, ierr) ; call CHKERR(ierr)
      call MatSetNearNullSpace(A, nullspace, ierr);call CHKERR(ierr)
    endif

    call KSPSetFromOptions(ksp,ierr) ;call CHKERR(ierr)

    linit = .True.
    if(myid.eq.0.and.ldebug) print *,'Setup KSP done'
    call PetscLogStagePop(ierr) ;call CHKERR(ierr)


  contains
    !> @brief define physical coordinates for DMDA to allow for geometric multigrid
    subroutine set_coordinates(C,ierr)
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
            xv(0,i,j,k) = xv(0,i+1,j,k) + atm%dz(atmk(i),j,k)
          enddo
        enddo
      enddo

      xv => null()
      call VecRestoreArrayF90(coordinates,xv1d,ierr) ;call CHKERR(ierr)
      xv1d => null()
    end subroutine
  end subroutine

  !> @brief override convergence tests -- the normal KSPConverged returns bad solution if no iterations are needed for convergence
  subroutine MyKSPConverged(ksp,n,rnorm,flag,dummy,ierr)
    ! Input Parameters:
    !    ksp   - iterative context
    !    n     - iteration number
    !    rnorm - 2-norm (preconditioned) residual value (may be estimated)
    !    dummy - optional user-defined monitor context (unused here)
    KSP               :: ksp
    PetscErrorCode    :: ierr
    PetscInt          :: n,dummy
    KSPConvergedReason:: flag
    PetscReal         :: rnorm

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
  end subroutine

  !> @brief assign string names to logging levels
  subroutine setup_logging()
      if(.not. allocated(logstage)) then
          allocate(logstage(13))
          call PetscLogStageRegister('total_tenstream' , logstage(1)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('setup_edir'      , logstage(2)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('calc_edir'       , logstage(3)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('setup_ediff'     , logstage(4)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('calc_ediff'      , logstage(5)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('setup_b'         , logstage(6)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('get_coeff'       , logstage(7)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('twostream'       , logstage(8)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('setup_ksp'       , logstage(9)     , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('write_hdf5'      , logstage(10)    , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('load_save_sol'   , logstage(11)    , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('nca'             , logstage(12)    , ierr) ;call CHKERR(ierr)
          call PetscLogStageRegister('schwarzschild'   , logstage(13)    , ierr) ;call CHKERR(ierr)

          if(myid.eq.0 .and. ldebug) print *, 'Logging stages' , logstage
      endif
  end subroutine

  !> @brief nca wrapper to call NCA of Carolin Klinger
  !> @details This is supposed to work on a 1D solution which has to be calculated beforehand
  !> \n the wrapper copies fluxes and optical properties on one halo and then gives that to NCA
  !> \n the result is the 3D approximation of the absorption, considering neighbouring information
  subroutine nca_wrapper(ediff,abso)
    use m_ts_nca, only : ts_nca
    Vec :: ediff,abso
    Vec :: gnca ! global nca vector 
    Vec :: lnca ! local nca vector with ghost values -- in dimension 0 and 1 are fluxes followed by dz,planck,kabs

    PetscReal,pointer,dimension(:,:,:,:) :: xv  =>null()
    PetscReal,pointer,dimension(:)       :: xv1d=>null()
    PetscReal,pointer,dimension(:,:,:,:) :: xvlnca  =>null(), xvgnca  =>null()
    PetscReal,pointer,dimension(:)       :: xvlnca1d=>null(), xvgnca1d=>null()
    PetscReal,pointer,dimension(:,:,:,:) :: xhr  =>null()
    PetscReal,pointer,dimension(:)       :: xhr1d=>null()
    integer(iintegers) :: k

    integer(iintegers),parameter :: idz=i2, iplanck=i3, ikabs=i4, ihr=i5

    call PetscLogStagePush(logstage(12),ierr) ;call CHKERR(ierr)

    ! put additional values into a local ediff vec .. TODO: this is a rather dirty hack but is straightforward

    ! get ghost values for dz, planck, kabs and fluxes, ready to give it to NCA
    call DMGetGlobalVector(C_diff%da ,gnca ,ierr) ; call CHKERR(ierr)

    call getVecPointer(gnca ,C_diff ,xvgnca1d, xvgnca)
    xvgnca(  idz    , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%dz
    xvgnca(  iplanck, C_diff%zs:C_diff%ze  , C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%planck
    xvgnca(  ikabs  , C_diff%zs:C_diff%ze-1, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye ) = atm%op%kabs


    ! Copy Edn and Eup to local convenience vector
    call getVecPointer(ediff ,C_diff ,xv1d, xv)
    xvgnca( E_up,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_up,:,:,:)
    xvgnca( E_dn,:,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) = xv( E_dn,:,:,:)
    call restoreVecPointer(ediff ,C_diff ,xv1d, xv)

    call restoreVecPointer(gnca ,C_diff ,xvgnca1d, xvgnca )


    ! retrieve ghost values into l(ocal) nca vec
    call DMGetLocalVector (C_diff%da ,lnca ,ierr) ; call CHKERR(ierr)
    call VecSet(lnca, zero, ierr); call CHKERR(ierr)

    call DMGlobalToLocalBegin(C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)
    call DMGlobalToLocalEnd  (C_diff%da ,gnca ,ADD_VALUES,lnca ,ierr) ; call CHKERR(ierr)

    call DMRestoreGlobalVector(C_diff%da, gnca, ierr); call CHKERR(ierr)

    ! call NCA
    call getVecPointer(lnca ,C_diff ,xvlnca1d, xvlnca)

    call ts_nca( atm%dx, atm%dy,                    &
      xvlnca(   idz        , : , : , :), &
      xvlnca(   iplanck    , : , : , :), &
      xvlnca(   ikabs      , : , : , :), &
      xvlnca(   E_dn       , : , : , :), &
      xvlnca(   E_up       , : , : , :), &
      xvlnca(   ihr        , : , : , :))


    ! return absorption
    call getVecPointer( abso, C_one ,xhr1d, xhr)

    do k=C_one%zs,C_one%ze
      xhr(i0,k,:,:) = xvlnca( ihr , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye) / xvlnca( idz , k,C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye)
    enddo
    call restoreVecPointer( abso ,C_one ,xhr1d, xhr )

    !return convenience vector that holds optical properties
    call restoreVecPointer(lnca ,C_diff ,xvlnca1d, xvlnca )
    call DMRestoreLocalVector(C_diff%da, lnca, ierr); call CHKERR(ierr)

    call PetscLogStagePop(ierr) ;call CHKERR(ierr)
  end subroutine

  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine schwarz(solution)
    type(t_state_container) :: solution

    PetscReal,pointer,dimension(:,:,:,:) :: xv_diff=>null()
    PetscReal,pointer,dimension(:)       :: xv_diff1d=>null()
    integer(iintegers) :: i,j

    real(ireals),allocatable :: dtau(:),Edn(:),Eup(:)

    if(solution%lsolar_rad) stop 'Tried calling schwarschild solver for solar calculation -- stopping!'
    if( .not. allocated(atm%planck) ) stop 'Tried calling schwarschild solver but no planck was given -- stopping!' 

    call PetscLogStagePush(logstage(13),ierr) ;call CHKERR(ierr)

    call VecSet(solution%ediff,zero,ierr); call CHKERR(ierr)

    allocate( dtau(C_diff%zm-1) )

    call getVecPointer(solution%ediff ,C_diff ,xv_diff1d, xv_diff)

    allocate( Eup(C_diff%zm) )
    allocate( Edn(C_diff%zm) )

    if(myid.eq.0 .and. ldebug) print *,' CALCULATING schwarzschild ::'

    do j=C_diff%ys,C_diff%ye         
      do i=C_diff%xs,C_diff%xe

        dtau = atm%dz(atmk(C_one%zs):C_one%ze,i,j)* atm%op(atmk(C_one%zs):C_one%ze,i,j)%kabs

        call schwarzschild( dtau ,atm%albedo(i,j), Edn,Eup, atm%planck(atmk(C_one1%zs):C_one1%ze,i,j) )

        xv_diff(E_up,:,i,j) = Eup(:) 
        xv_diff(E_dn,:,i,j) = Edn(:) 
      enddo
    enddo

    call restoreVecPointer(solution%ediff ,C_diff ,xv_diff1d, xv_diff )

    !Schwarzschild solver returns fluxes as [W/m^2]
    solution%lintegrated_dir  = .False.
    solution%lintegrated_diff = .False.
    ! and mark solution that it is not up to date
    solution%lchanged         = .True. 

    deallocate(Edn)
    deallocate(Eup)

    call PetscLogStagePop(ierr) ;call CHKERR(ierr)
  end subroutine

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine twostream(edirTOA,solution)
    real(ireals),intent(in) :: edirTOA
    type(t_state_container) :: solution

    PetscReal,pointer,dimension(:,:,:,:) :: xv_dir=>null(),xv_diff=>null()
    PetscReal,pointer,dimension(:) :: xv_dir1d=>null(),xv_diff1d=>null()
    integer(iintegers) :: i,j,src

    real(ireals),allocatable :: dtau(:),kext(:),w0(:),g(:),S(:),Edn(:),Eup(:)
    real(ireals) :: mu0,incSolar

    if(solution%lsolar_rad) &
      call VecSet(solution%edir ,zero,ierr); call CHKERR(ierr)

    call VecSet(solution%ediff,zero,ierr); call CHKERR(ierr)

    call PetscLogStagePush(logstage(8),ierr) ;call CHKERR(ierr)

    allocate( dtau(C_one_atm%zm) )
    allocate( kext(C_one_atm%zm) )
    allocate(   w0(C_one_atm%zm) )
    allocate(    g(C_one_atm%zm) )


    if(solution%lsolar_rad) &
      call getVecPointer(solution%edir  ,C_dir  ,xv_dir1d , xv_dir)
    call getVecPointer(solution%ediff ,C_diff ,xv_diff1d, xv_diff)

    allocate( S  (C_one_atm1%zs:C_one_atm1%ze) )
    allocate( Eup(C_one_atm1%zs:C_one_atm1%ze) )
    allocate( Edn(C_one_atm1%zs:C_one_atm1%ze) )


    do j=C_one_atm%ys,C_one_atm%ye
      do i=C_one_atm%xs,C_one_atm%xe

        mu0 = sun%angles(C_one_atm1%zs,i,j)%costheta
        incSolar = edirTOA* mu0
        if(myid.eq.0 .and. ldebug) print *,' CALCULATING DELTA EDDINGTON TWOSTREAM ::',sun%angles(C_one_atm1%zs,i,j)%theta,':',incSolar

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
          do src=i0,i3
            xv_dir(src,C_dir%zs+1:C_dir%ze,i,j) = S(atmk(C_one_atm1%zs)+1:C_one_atm1%ze)
            xv_dir(src,C_dir%zs           ,i,j) = S(C_one_atm1%zs)
          enddo
        endif

        xv_diff(E_up,C_diff%zs+1:C_diff%ze,i,j) = Eup(atmk(C_one_atm1%zs)+1:C_one_atm1%ze) 
        xv_diff(E_up,C_diff%zs            ,i,j) = Eup(C_one_atm1%zs) 
        xv_diff(E_dn,C_diff%zs+1:C_diff%ze,i,j) = Edn(atmk(C_one_atm1%zs)+1:C_one_atm1%ze) 
        xv_diff(E_dn,C_diff%zs            ,i,j) = Edn(C_one_atm1%zs) 
      enddo
    enddo

    if(solution%lsolar_rad) &
      call restoreVecPointer(solution%edir  ,C_dir  ,xv_dir1d , xv_dir  )
    call restoreVecPointer(solution%ediff ,C_diff ,xv_diff1d, xv_diff )

    !Twostream solver returns fluxes as [W]
    solution%lintegrated_dir  = .False.
    solution%lintegrated_diff = .False.
    ! and mark solution that it is not up to date
    solution%lchanged         = .True. 

    deallocate(S)
    deallocate(Edn)
    deallocate(Eup)

    call PetscLogStagePop(ierr) ;call CHKERR(ierr)
  end subroutine

  !> @brief renormalize fluxes with the size of a face(sides or lid)
  subroutine scale_flx(solution, lWm2_to_W)
    type(t_state_container),intent(inout) :: solution   !< @param solution container with computed fluxes
    logical,intent(in)                    :: lWm2_to_W  !< @param determines direction of scaling, if true, scale from W/m**2 to W 

    if(solution%lsolar_rad) then
      if(solution%lintegrated_dir .neqv. lWm2_to_W) then
        call scale_flx_vec(solution%edir, C_dir, lWm2_to_W)
        solution%lintegrated_dir = lWm2_to_W
      endif
    endif
    if(solution%lintegrated_diff .neqv. lWm2_to_W) then
      call scale_flx_vec(solution%ediff, C_diff, lWm2_to_W)
      solution%lintegrated_diff = lWm2_to_W
    endif

  contains
    subroutine scale_flx_vec(v, C, lWm2_to_W)
      Vec :: v
      type(t_coord) :: C
      PetscReal,pointer,dimension(:,:,:,:) :: xv  =>null()
      PetscReal,pointer,dimension(:)       :: xv1d=>null()
      PetscInt :: i,j,k
      PetscReal :: Ax,Ax2,Ay,Ay2,Az,Az4
      logical,intent(in) :: lWm2_to_W ! determines direction of scaling, if true, scale from W/m**2 to W

      Vec :: vgrad_x, vgrad_y
      PetscScalar,Pointer :: grad_x(:,:,:,:)=>null(), grad_x1d(:)=>null()
      PetscScalar,Pointer :: grad_y(:,:,:,:)=>null(), grad_y1d(:)=>null()
      real(ireals) :: grad(3)  ! is the cos(zenith_angle) of the tilted box in case of topography

      if(myid.eq.0.and.ldebug) print *,'rescaling fluxes',C%zm,C%xm,C%ym
      call getVecPointer(v ,C ,xv1d, xv)

      if(lWm2_to_W) then
        Az  = atm%dx*atm%dy
        Az4 = atm%dx*atm%dy*.25_ireals
      else
        Az  = one/(atm%dx*atm%dy)
        Az4 = one/(atm%dx*atm%dy*.25_ireals)
      endif

      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-i1

            if(C%dof.eq.i8) then ! This is 8 stream direct radiation
              xv(i0:i3,k,i,j) = xv(i0:i3,k,i,j) * Az4
            endif

            if(C%dof.eq.i10) then ! This is 10 stream diffuse radiation
              xv(E_up  ,k,i,j) = xv(E_up  ,k,i,j) * Az
              xv(E_dn  ,k,i,j) = xv(E_dn  ,k,i,j) * Az
            endif

            if(.not.atm%l1d(atmk(k),i,j)) then

              if(C%dof.eq.i8) then ! This is 8 stream direct radiation

                if(lWm2_to_W) then
                  Ax2 = atm%dy*atm%dz(k,i,j)*.5_ireals
                  Ay2 = atm%dx*atm%dz(k,i,j)*.5_ireals
                else
                  Ax2 = one/(atm%dy*atm%dz(k,i,j)*.5_ireals )
                  Ay2 = one/(atm%dx*atm%dz(k,i,j)*.5_ireals )
                endif
                xv(i4:i5,k,i,j) = xv(i4:i5,k,i,j) * Ax2 
                xv(i6:i7,k,i,j) = xv(i6:i7,k,i,j) * Ay2 
              endif

              if(C%dof.eq.i10) then ! This is 10 stream diffuse radiation
                if(lWm2_to_W) then
                  Ax  = atm%dy*atm%dz(k,i,j)
                  Ay  = atm%dx*atm%dz(k,i,j)
                else
                  Ax  = one/(atm%dy*atm%dz(k,i,j) )
                  Ay  = one/(atm%dx*atm%dz(k,i,j) )
                endif
                xv(E_le_m,k,i,j) = xv(E_le_m,k,i,j) * Ax
                xv(E_le_p,k,i,j) = xv(E_le_p,k,i,j) * Ax
                xv(E_ri_m,k,i,j) = xv(E_ri_m,k,i,j) * Ax
                xv(E_ri_p,k,i,j) = xv(E_ri_p,k,i,j) * Ax
                xv(E_ba_m,k,i,j) = xv(E_ba_m,k,i,j) * Ay
                xv(E_ba_p,k,i,j) = xv(E_ba_p,k,i,j) * Ay
                xv(E_fw_m,k,i,j) = xv(E_fw_m,k,i,j) * Ay
                xv(E_fw_p,k,i,j) = xv(E_fw_p,k,i,j) * Ay
              endif
            endif
          enddo
        enddo
      enddo

      k=C%ze
      do j=C%ys,C%ye
        do i=C%xs,C%xe

          if(C%dof.eq.i8) then ! This is 8 stream direct radiation
            xv (i0:i3 ,k,i,j) = xv (i0:i3 ,k,i,j) * Az4
          endif
          if(C%dof.eq.i10) then ! This is 10 stream diffuse radiation
            xv(E_up  ,k,i,j) = xv(E_up  ,k,i,j) * Az
            xv(E_dn  ,k,i,j) = xv(E_dn  ,k,i,j) * Az
          endif
        enddo
      enddo

      if(sun%luse_topography) then ! This is direct rad and we use topography !todo do we need this?
        select case (C%dof)

        case(i8)
          call compute_gradient(atm, vgrad_x, vgrad_y)

          call getVecPointer(vgrad_x , C_one1, grad_x1d, grad_x)
          call getVecPointer(vgrad_y , C_one1, grad_y1d, grad_y)

          do j=C%ys,C%ye
            do i=C%xs,C%xe
              do k=C%zs,C%ze
                grad(1) = grad_x(i0,k,i,j)
                grad(2) = grad_y(i0,k,i,j)
                grad(3) = one

                xv(i0:i3,k,i,j) = xv(i0:i3,k,i,j) / norm(grad)
              enddo
            enddo
          enddo

          call restoreVecPointer(vgrad_x , C_one1, grad_x1d, grad_x)
          call restoreVecPointer(vgrad_y , C_one1, grad_y1d, grad_y)

          call DMRestoreLocalVector(C_one1%da, vgrad_x, ierr);  call CHKERR(ierr)
          call DMRestoreLocalVector(C_one1%da, vgrad_y, ierr);  call CHKERR(ierr)

        case(i6)
          ! Dont rescale diffuse fluxes

        case(i10)

          ! Dont rescale diffuse fluxes

          ! call compute_gradient(atm, vgrad_x, vgrad_y)

          ! call getVecPointer(vgrad_x , C_one1, grad_x1d, grad_x)
          ! call getVecPointer(vgrad_y , C_one1, grad_y1d, grad_y)

          !do j=C%ys,C%ye
          !  do i=C%xs,C%xe
          !    do k=C%zs,C%ze
          !      grad(1) = grad_x(i0,k,i,j)
          !      grad(2) = grad_y(i0,k,i,j)
          !      grad(3) = one

          !      xv([E_up, E_dn],k,i,j) = xv([E_up, E_dn],k,i,j) / norm(grad)
          !    enddo
          !  enddo
          !enddo

        case default
          stop('Dont know how I should topography rescale this! - exiting...')
        end select

      endif

      call restoreVecPointer(v ,C ,xv1d, xv )
    end subroutine
  end subroutine

! Deprecated --  ! Deprecated -- this is probably not helping convergence....
! Deprecated --  subroutine set_diff_initial_guess(inp,guess,C)
! Deprecated --    Vec :: inp,guess
! Deprecated --    type(t_coord) :: C
! Deprecated --
! Deprecated --    Vec :: local_guess
! Deprecated --    PetscScalar,pointer,dimension(:,:,:,:) :: xinp=>null(),xguess=>null()
! Deprecated --    PetscScalar,pointer,dimension(:) :: xinp1d=>null(),xguess1d=>null()
! Deprecated --    PetscReal :: diff2diff(C_diff%dof**2)
! Deprecated --    PetscInt :: k,i,j,src
! Deprecated --
! Deprecated --    if(myid.eq.0.and.ldebug) print *,'setting initial guess...'
! Deprecated --    call DMGetLocalVector(C%da,local_guess,ierr) ;call CHKERR(ierr)
! Deprecated --    call VecSet(local_guess,zero,ierr) ;call CHKERR(ierr)
! Deprecated --
! Deprecated --    call getVecPointer(inp ,C ,xinp1d, xinp ,.False.)
! Deprecated --    call getVecPointer(local_guess ,C ,xguess1d, xguess ,.True.)
! Deprecated --
! Deprecated --    do j=C%ys,C%ye
! Deprecated --      do i=C%xs,C%xe
! Deprecated --        do k=C%zs,C%ze-1
! Deprecated --          if( .not. atm%l1d(k,i,j) ) then
! Deprecated --            call get_coeff(atm%op(k,i,j), atm%dz(k,i,j),.False., diff2diff, atm%l1d(k,i,j) )
! Deprecated --            do src=1,C%dof
! Deprecated --              xguess(E_up   , k   , i   , j   ) = xguess(E_up   , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_up  +i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_dn   , k+1 , i   , j   ) = xguess(E_dn   , k+1 , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_dn  +i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_le_m , k   , i   , j   ) = xguess(E_le_m , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_le_m+i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_le_p , k   , i   , j   ) = xguess(E_le_p , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_le_p+i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_ri_m , k   , i+1 , j   ) = xguess(E_ri_m , k   , i+1 , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ri_m+i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_ri_p , k   , i+1 , j   ) = xguess(E_ri_p , k   , i+1 , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ri_p+i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_ba_m , k   , i   , j   ) = xguess(E_ba_m , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ba_m+i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_ba_p , k   , i   , j   ) = xguess(E_ba_p , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ba_p+i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_fw_m , k   , i   , j+1 ) = xguess(E_fw_m , k   , i   , j+1 ) +  xinp(src-1 , k , i , j ) *diff2diff(E_fw_m+i1+(src-1 ) *C%dof )
! Deprecated --              xguess(E_fw_p , k   , i   , j+1 ) = xguess(E_fw_p , k   , i   , j+1 ) +  xinp(src-1 , k , i , j ) *diff2diff(E_fw_p+i1+(src-1 ) *C%dof )
! Deprecated --            enddo
! Deprecated --          endif
! Deprecated --        enddo
! Deprecated --      enddo
! Deprecated --    enddo
! Deprecated --
! Deprecated --    call restoreVecPointer(inp ,C ,xinp1d, xinp )
! Deprecated --    call restoreVecPointer(local_guess ,C ,xguess1d, xguess )
! Deprecated --
! Deprecated --    call VecSet(guess,zero,ierr) ;call CHKERR(ierr) ! reset global Vec
! Deprecated --
! Deprecated --    call DMLocalToGlobalBegin(C%da,local_guess,ADD_VALUES, guess,ierr) ;call CHKERR(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
! Deprecated --    call DMLocalToGlobalEnd  (C%da,local_guess,ADD_VALUES, guess,ierr) ;call CHKERR(ierr)
! Deprecated --
! Deprecated --    call DMRestoreLocalVector(C_diff%da, local_guess, ierr); call CHKERR(ierr)
! Deprecated --  end subroutine

  !> @brief initialize basic memory structs like PETSc vectors and matrices
  subroutine init_memory(incSolar,b)
    Vec, allocatable :: b,incSolar

    if(ltwostr_only) return

    if(.not.allocated(incSolar)) allocate(incSolar)
    if(.not.allocated(b)) allocate(b)

    call DMCreateGlobalVector(C_dir%da,incSolar,ierr) ; call CHKERR(ierr)
    call DMCreateGlobalVector(C_diff%da,b,ierr)       ; call CHKERR(ierr)

    call VecSet(incSolar,zero,ierr) ; call CHKERR(ierr)
    call VecSet(b,zero,ierr)        ; call CHKERR(ierr)
  end subroutine

  subroutine init_matrices()
    ! if already initialized, tear em down -- this usually
    ! happens when we change solar angles and need a new prealloc

    ! Notice that you can not call destroy_matrices() and init_matrices() and get new matrices,
    ! rather just call init_mat again and again -- PETSc somehow buffers the old matrix.
    ! It only really gets rebuilt if we also delete the DMDA's as in destroy_tenstream()

    call init_Matrix(Mdir ,C_dir )
    call init_Matrix(Mdiff,C_diff)
  end subroutine

  subroutine destroy_matrices()
    if(myid.eq.0 .and. ldebug) print *,'Trying to destroy matrices...', allocated(Mdir), allocated(Mdiff)
    if(allocated(Mdir)) then
      call MatDestroy(Mdir , ierr) ;call CHKERR(ierr)
      deallocate(Mdir)
    endif
    if(allocated(Mdiff)) then
      call MatDestroy(Mdiff, ierr) ;call CHKERR(ierr)
      deallocate(Mdiff)
    endif
  end subroutine

  subroutine set_angles(phi0, theta0, phi2d, theta2d)
    real(ireals),intent(in)          :: phi0           !< @param[in] phi0   solar azmiuth and zenith angle
    real(ireals),intent(in)          :: theta0         !< @param[in] theta0 solar azmiuth and zenith angle
    real(ireals),optional,intent(in) :: phi2d  (:,:)   !< @param[in] phi2d   if given, horizontally varying azimuth
    real(ireals),optional,intent(in) :: theta2d(:,:)   !< @param[in] theta2d if given, and zenith angle

    logical :: lchanged_theta, lchanged_phi

    if(.not.ltenstream_is_initialized) then
        print *,myid,'You tried to set angles in the Tenstream solver.  &
            & This should be called right after init_tenstream'
        ierr=1; call CHKERR(ierr)
    endif

    if(allocated(sun%angles)) then ! was initialized
      if(present(theta2d)) then
        lchanged_theta = .not. all(theta2d.eq.sun%angles(1,:,:)%theta)
      else
        lchanged_theta = .not. all(theta0.eq.sun%angles(1,:,:)%theta)
      endif
      if(present(phi2d)) then
        lchanged_phi = .not. all(phi2d.eq.sun%angles(1,:,:)%phi)
      else
        lchanged_phi = .not. all(phi0.eq.sun%angles(1,:,:)%phi)
      endif
      if(myid.eq.0 .and. ldebug) print *,'tenstr set_angles -- changed angles?',lchanged_theta, lchanged_phi
      if(.not. lchanged_theta .and. .not. lchanged_phi) then
        return
      endif
    endif


    if ( present(phi2d) .and. present(theta2d) ) then
        call setup_suninfo(phi0, theta0, sun, phi2d=phi2d, theta2d=theta2d)
    elseif ( present(phi2d) ) then
        call setup_suninfo(phi0, theta0, sun, phi2d=phi2d)
    elseif ( present(theta2d) ) then
        call setup_suninfo(phi0, theta0, sun, theta2d=theta2d)
    else
        call setup_suninfo(phi0, theta0, sun)
    endif


    ! init box montecarlo model
    if(any(atm%l1d.eqv..False.)) call OPP_8_10%init(pack(sun%angles%symmetry_phi,.True.),pack(sun%angles%theta,.True.),imp_comm)
    if(.not.luse_eddington)      call OPP_1_2%init (pack(sun%angles%symmetry_phi,.True.),pack(sun%angles%theta,.True.),imp_comm)

    call init_matrices()
  end subroutine

  !> @brief Main routine to setup TenStream solver
  !> @details This will setup the PETSc DMDA grid and set other grid information, needed for the TenStream
  !> \n Nx, Ny Nz are either global domain size or have to be local sizes if present(nxproc,nyproc)
  !> \n where nxproc and nyproc then are the number of pixel per rank for all ranks -- i.e. sum(nxproc) != Nx_global
  subroutine init_tenstream(icomm, Nz,Nx,Ny, dx,dy, phi0, theta0, dz1d, dz3d, nxproc, nyproc, collapseindex)
    MPI_Comm, intent(in)          :: icomm   !< @param MPI_Communicator which should be used -- this will be used for PETSC_COMM_WORLD
    integer(iintegers),intent(in) :: Nz      !< @param[in] Nz     Nz is the number of layers and Nz+1 would be the number of levels
    integer(iintegers),intent(in) :: Nx      !< @param[in] Nx     number of boxes in x-direction
    integer(iintegers),intent(in) :: Ny      !< @param[in] Ny     number of boxes in y-direction
    real(ireals),intent(in)       :: dx      !< @param[in] dx     physical size of grid in [m]
    real(ireals),intent(in)       :: dy      !< @param[in] dy     physical size of grid in [m]
    real(ireals),intent(in)       :: phi0    !< @param[in] phi0   solar azmiuth and zenith angle
    real(ireals),intent(in)       :: theta0  !< @param[in] theta0 solar azmiuth and zenith angle

    real(ireals),optional,intent(in)       :: dz1d(:)        !< @param[in] dz1d    if given, dz1d is used everywhere on the rank
    real(ireals),optional,intent(in)       :: dz3d(:,:,:)    !< @param[in] dz3d    if given, dz3d has to be local domain size, cannot have global shape
    integer(iintegers),optional,intent(in) :: nxproc(:)      !< @param[in] nxproc  if given, Nx has to be the local size, dimension of nxproc is number of ranks along x-axis, and entries in nxproc are the size of local Nx
    integer(iintegers),optional,intent(in) :: nyproc(:)      !< @param[in] nyproc  if given, Ny has to be the local size, dimension of nyproc is number of ranks along y-axis, and entries in nyproc are the number of local Ny
    integer(iintegers),optional,intent(in) :: collapseindex  !< @param[in] collapseindex if given, the upper n layers will be reduce to 1d and no individual output will be given for them

    integer(iintegers) :: k,i,j
    !    character(default_str_len),parameter :: tenstreamrc='./.tenstreamrc'

    if(.not.ltenstream_is_initialized) then

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
                call setup_grid( Nz, Nx, Ny, nxproc,nyproc, collapseindex=collapseindex)
            else
                call setup_grid( Nz, Nx, Ny, nxproc,nyproc)
            endif
        else
            if(present(collapseindex)) then
                call setup_grid( Nz, max(minimal_dimension, Nx), max(minimal_dimension, Ny), collapseindex=collapseindex)
            else
                call setup_grid( Nz, max(minimal_dimension, Nx), max(minimal_dimension, Ny) )
            endif
        endif

        call setup_atm()

        ! init work vectors
        call init_memory(incSolar,b)

        ! init petsc logging facilities
        call setup_logging()

        ltenstream_is_initialized=.True.
    else
        print *,myid,'You tried to initialize already initialized Tenstream          &
            &solver. This should not be done. If you need to reinitialize the grids, &
            &call destroy_tenstream() first.'
        ierr=1; call CHKERR(ierr)
    endif

    !Todo: this is just here so that we do not break the API. User could also
    !       directly use the set_angle routines?
    call set_angles(phi0, theta0)

  contains
      subroutine setup_atm()
          if(.not.allocated(atm)) allocate(atm)

          atm%dx  = dx
          atm%dy  = dy

          if(.not.allocated(atm%dz) ) allocate(atm%dz( C_one_atm%zs:C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye ))

          if(present(dz1d)) then
              do j=C_one_atm%ys,C_one_atm%ye
                  do i=C_one_atm%xs,C_one_atm%xe
                      atm%dz(:,i,j) = dz1d
                  enddo
              enddo
          endif
          if(present(dz3d)) then
              if( any( shape(dz3d).ne.shape(atm%dz) ) ) then
                  print *,'Whoops I got a 3D dz definition but this does not correspond to the grid definition :: shapes: ', shape(dz3d), ' vs. ',shape(atm%dz)
                  print *,'please know that providing a 3D dz profile has to be of local shape, it can not have global size'
                  call MPI_Abort(icomm,ierr,ierr)
              endif
              atm%dz = dz3d

          endif
          if(.not.present(dz1d) .and. .not.present(dz3d)) then
              print *,'have to give either dz1d or dz3d in routine call....'
              ierr=1; call CHKERR(ierr)
          endif

          sun%luse_topography = present(dz3d) .and. ltopography  ! if the user supplies 3d height levels and has set the topography option

          if(.not.allocated(atm%l1d)) then
              allocate(atm%l1d( C_one_atm%zs:C_one_atm%ze, C_one_atm%xs:C_one_atm%xe, C_one_atm%ys:C_one_atm%ye ) )
          endif

          !TODO if we have a horiz. staggered grid, this may lead to the point where one 3d box has a outgoing sideward flux but the adjacent
          !1d box does not send anything back --> therefore huge absorption :( -- would need to introduce mirror boundary conditions for
          !sideward fluxes in 1d boxes
          do j=C_one_atm%ys,C_one_atm%ye
              do i=C_one_atm%xs,C_one_atm%xe
                  atm%l1d(C_one_atm%ze,i,j) = twostr_ratio*atm%dz(C_one_atm%ze,i,j).gt.atm%dx
                  do k=C_one_atm%ze-1,C_one_atm%zs,-1
                      atm%l1d(k,i,j) = twostr_ratio*atm%dz(k,i,j).gt.atm%dx
                  enddo
              enddo
          enddo
          if(ltwostr_only) atm%l1d = .True.

          if(present(collapseindex)) then
              atm%lcollapse=.True.
              atm%icollapse=collapseindex
              atm%l1d(C_one_atm%zs:atmk(C_one%zs),:,:) = .True. ! if need to be collapsed, they have to be 1D.
          endif
      end subroutine

  end subroutine

  subroutine set_global_optical_properties(albedo, global_kabs, global_ksca, global_g, global_planck)
    real(ireals),intent(inout) :: albedo
    real(ireals),intent(inout),dimension(:,:,:),allocatable,optional :: global_kabs, global_ksca, global_g
    real(ireals),intent(inout),dimension(:,:,:),allocatable,optional :: global_planck
    real(ireals),dimension(:,:,:),allocatable :: local_kabs, local_ksca, local_g
    real(ireals),dimension(:,:,:),allocatable :: local_planck
    logical :: lhave_planck,lhave_kabs,lhave_ksca,lhave_g

    if(.not.ltenstream_is_initialized) then
      print *,myid,'You tried to set global optical properties but tenstream environment seems not to be initialized.... please call init first!'
      call exit(1)
    endif

    call imp_bcast(imp_comm, albedo, 0_mpiint)

    lhave_kabs   = present(global_kabs  ); call imp_bcast(imp_comm, lhave_kabs  , 0_mpiint)
    lhave_ksca   = present(global_ksca  ); call imp_bcast(imp_comm, lhave_ksca  , 0_mpiint)
    lhave_g      = present(global_g     ); call imp_bcast(imp_comm, lhave_g     , 0_mpiint)
    lhave_planck = present(global_planck); call imp_bcast(imp_comm, lhave_planck, 0_mpiint)

    ! Make sure that our domain has at least 3 entries in each dimension.... otherwise violates boundary conditions
    if(myid.eq.0) then
      if( lhave_kabs   ) call extend_arr(global_kabs)
      if( lhave_ksca   ) call extend_arr(global_ksca)
      if( lhave_g      ) call extend_arr(global_g)
      if( lhave_planck ) call extend_arr(global_planck)
    endif

    if( lhave_kabs   ) allocate( local_kabs   (C_one%zs :C_one%ze ,C_one%xs :C_one%xe , C_one%ys :C_one%ye  ) )
    if( lhave_ksca   ) allocate( local_ksca   (C_one%zs :C_one%ze ,C_one%xs :C_one%xe , C_one%ys :C_one%ye  ) )
    if( lhave_g      ) allocate( local_g      (C_one%zs :C_one%ze ,C_one%xs :C_one%xe , C_one%ys :C_one%ye  ) )
    if( lhave_planck ) allocate( local_planck (C_one1%zs:C_one1%ze,C_one1%xs:C_one1%xe, C_one1%ys:C_one1%ye ) )

    ! Scatter global optical properties to MPI nodes
    call local_optprop()
    ! Now global_fields are local to mpi subdomain.

    if(lhave_planck) then
      call set_optical_properties(albedo, local_kabs, local_ksca, local_g, local_planck)
    else
      call set_optical_properties(albedo, local_kabs, local_ksca, local_g)
    endif


  contains
    subroutine local_optprop()
      Vec :: local_vec
      PetscReal,pointer,dimension(:,:,:,:) :: xlocal_vec  =>null()
      PetscReal,pointer,dimension(:)       :: xlocal_vec1d=>null()

      if(myid.eq.0.and.ldebug .and. lhave_kabs) &
        print *,myid,'copying optprop: global to local :: shape kabs',shape(global_kabs),'xstart/end',C_one_atm%xs,C_one_atm%xe,'ys/e',C_one_atm%ys,C_one_atm%ye

      call DMGetGlobalVector(C_one_atm%da, local_vec, ierr) ; call CHKERR(ierr)

      if(lhave_kabs) then
        call scatterZerotoDM(global_kabs, C_one_atm, local_vec)
        call getVecPointer(local_vec ,C_one_atm ,xlocal_vec1d, xlocal_vec)
        local_kabs = xlocal_vec(0, C_one_atm%zs :C_one_atm%ze,C_one_atm%xs :C_one_atm%xe , C_one_atm%ys :C_one_atm%ye  )
        call restoreVecPointer(local_vec ,C_one_atm ,xlocal_vec1d, xlocal_vec )
      endif

      if(lhave_ksca) then
        call scatterZerotoDM(global_ksca,C_one_atm,local_vec)
        call getVecPointer(local_vec ,C_one_atm ,xlocal_vec1d, xlocal_vec)
        local_ksca = xlocal_vec(0,:,:,:)
        call restoreVecPointer(local_vec ,C_one_atm ,xlocal_vec1d, xlocal_vec )
      endif

      if(lhave_g) then
        call scatterZerotoDM(global_g,C_one_atm,local_vec)
        call getVecPointer(local_vec ,C_one_atm ,xlocal_vec1d, xlocal_vec)
        local_g = xlocal_vec(0,:,:,:)
        call restoreVecPointer(local_vec ,C_one_atm ,xlocal_vec1d, xlocal_vec )
      endif

      call DMRestoreGlobalVector(C_one_atm%da, local_vec, ierr) ; call CHKERR(ierr)

      if(lhave_planck) then
        call DMGetGlobalVector(C_one_atm1%da, local_vec, ierr) ; call CHKERR(ierr)
        call scatterZerotoDM(global_planck,C_one_atm1,local_vec)
        call getVecPointer(local_vec ,C_one_atm1 ,xlocal_vec1d, xlocal_vec)
        local_planck = xlocal_vec(0,:,:,:)
        call restoreVecPointer(local_vec ,C_one_atm1 ,xlocal_vec1d, xlocal_vec )
        call DMRestoreGlobalVector(C_one_atm1%da, local_vec, ierr) ; call CHKERR(ierr)
      endif
    end subroutine
  end subroutine

  subroutine scatterZerotoDM(arr, C, vec)
    real(ireals),allocatable,dimension(:,:,:),intent(in) :: arr
    type(t_coord),intent(in) :: C
    Vec :: vec

    VecScatter :: scatter_context
    Vec :: natural,local
    PetscScalar,Pointer :: xloc(:)=>null()

    if(ldebug) print *,myid,'scatterZerotoDM :: Create Natural Vec'
    call DMDACreateNaturalVector(C%da, natural, ierr); call CHKERR(ierr)

    if(ldebug) print *,myid,'scatterZerotoDM :: Create scatter ctx'
    call VecScatterCreateToZero(natural, scatter_context, local, ierr); call CHKERR(ierr)

    if(myid.eq.0) then
      call VecGetArrayF90(local,xloc,ierr) ;call CHKERR(ierr)
      if(ldebug) &
        print *,myid,'scatterZerotoDM :: shape of local',shape(xloc), 'shape of arr',shape(arr)
      xloc = reshape( arr , [ size(arr) ] )
      call VecRestoreArrayF90(local,xloc,ierr) ;call CHKERR(ierr)
    endif

    if(ldebug) print *,myid,'scatterZerotoDM :: scatter reverse....'
    call VecScatterBegin(scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
    call VecScatterEnd  (scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

    if(ldebug) print *,myid,'scatterZerotoDM :: natural to global....'
    call DMDANaturalToGlobalBegin(C%da,natural, INSERT_VALUES, vec, ierr); call CHKERR(ierr)
    call DMDANaturalToGlobalEnd  (C%da,natural, INSERT_VALUES, vec, ierr); call CHKERR(ierr)

    if(ldebug) print *,myid,'scatterZerotoDM :: destroying contexts....'
    call VecScatterDestroy(scatter_context, ierr); call CHKERR(ierr)
    call VecDestroy(local,ierr); call CHKERR(ierr)
    call VecDestroy(natural,ierr); call CHKERR(ierr)
    if(ldebug) print *,myid,'scatterZerotoDM :: done....'
  end subroutine

  subroutine extend_arr(arr)
    real(ireals),intent(inout),allocatable :: arr(:,:,:)
    real(ireals),allocatable :: tmp(:,:)
    integer(iintegers) :: dims(3),i

    if(.not. allocated(arr) ) print *,myid,'ERROR in SUBROUTINE extend_arr :: Cannot extend non allocated array!'

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

  subroutine set_optical_properties(albedo, local_kabs, local_ksca, local_g, local_planck, local_albedo_2d)
    real(ireals), intent(in) :: albedo
    real(ireals),intent(in),dimension(:,:,:),optional :: local_kabs, local_ksca, local_g ! dimensions (Nz  , Nx, Ny)
    real(ireals),intent(in),dimension(:,:,:),optional :: local_planck                    ! dimensions (Nz+1, Nx, Ny) layer quantity plus surface layer
    real(ireals),intent(in),dimension(:,:),optional   :: local_albedo_2d                 ! dimensions (Nx, Ny)
    real(ireals) :: tau,kext,w0,g
    integer(iintegers) :: k,i,j

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
              call adding(atm%a11(C_one_atm%zs:atmk(C_one%zs), i, j), &
                          atm%a12(C_one_atm%zs:atmk(C_one%zs), i, j), &
                          atm%a21(C_one_atm%zs:atmk(C_one%zs), i, j), &
                          atm%a22(C_one_atm%zs:atmk(C_one%zs), i, j), &
                          atm%a13(C_one_atm%zs:atmk(C_one%zs), i, j), &
                          atm%a23(C_one_atm%zs:atmk(C_one%zs), i, j), &
                          atm%a33(C_one_atm%zs:atmk(C_one%zs), i, j))
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

  subroutine solve_tenstream(edirTOA,opt_solution_uid,opt_solution_time)
    real(ireals),intent(in) :: edirTOA
    integer(iintegers),optional,intent(in) :: opt_solution_uid
    real(ireals),      optional,intent(in) :: opt_solution_time
    integer(iintegers) :: uid
    logical :: lsolar

    if(lenable_solutions_err_estimates .and. present(opt_solution_uid)) then
      uid = opt_solution_uid
    else
      uid = i0 ! default solution is uid==0
    endif

    lsolar = mpi_logical_and(imp_comm, edirTOA.gt.zero .and. any(sun%angles%theta.ge.zero))

    call prepare_solution( solutions(uid), uid, lsolar=lsolar ) ! setup solution vectors

    ! --------- Skip Thermal Computation (-lskip_thermal) --
    if(lskip_thermal .and. (solutions(uid)%lsolar_rad.eqv..False.) ) then !
      if(ldebug .and. myid.eq.0) print *,'skipping thermal calculation -- returning zero flux'
      call VecSet(solutions(uid)%ediff, zero, ierr); call CHKERR(ierr)
      solutions(uid)%lchanged=.True.
      call restore_solution(solutions(uid))
      return
    endif

    ! --------- Calculate 1D Radiative Transfer ------------
    if(  ltwostr                                                         &
      .or. all(atm%l1d.eqv..True.)                                       &
      .or. ((solutions(uid)%lsolar_rad.eqv..False.) .and. lcalc_nca)     &
      .or. ((solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild) ) then

    if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild ) then
      call schwarz(solutions(uid))
    else
      call twostream(edirTOA,  solutions(uid) )
    endif

    if(ldebug .and. myid.eq.0) print *,'1D calculation done'

    if(present(opt_solution_time) ) then
      call restore_solution(solutions(uid),opt_solution_time)
    else
      call restore_solution(solutions(uid))
    endif

    if( ltwostr_only ) return
    if( all(atm%l1d.eqv..True.) ) return
    if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lcalc_nca ) return
    if( (solutions(uid)%lsolar_rad.eqv..False.) .and. lschwarzschild ) return
  endif

  ! --------- scale from [W/m**2] to [W] -----------------
  call scale_flx(solutions(uid), lWm2_to_W=.True. )

  ! ---------------------------- Edir  -------------------
  if( solutions(uid)%lsolar_rad ) then

    call PetscLogStagePush(logstage(1),ierr) ;call CHKERR(ierr)
    call setup_incSolar(incSolar,edirTOA)
    call set_dir_coeff(Mdir,C_dir)

    call setup_ksp(kspdir,C_dir,Mdir,linit_kspdir, "dir_")

    call PetscLogStagePush(logstage(3),ierr) ;call CHKERR(ierr)
    call solve(kspdir,incSolar,solutions(uid)%edir)
    solutions(uid)%lchanged=.True.
    solutions(uid)%lintegrated_dir=.True.
    call PetscLogStagePop(ierr) ;call CHKERR(ierr)
  endif

  ! ---------------------------- Source Term -------------
  call setup_b(solutions(uid),b)

  ! ---------------------------- Ediff -------------------
  call set_diff_coeff(Mdiff,C_diff)
  call setup_ksp(kspdiff,C_diff,Mdiff,linit_kspdiff, "diff_")
  call PetscLogStagePush(logstage(5),ierr) ;call CHKERR(ierr)

  call solve(kspdiff, b, solutions(uid)%ediff,uid)
  solutions(uid)%lchanged=.True.
  solutions(uid)%lintegrated_diff=.True. !Tenstream solver returns fluxes as [W]

  call PetscLogStagePop(ierr) ;call CHKERR(ierr)

  if(present(opt_solution_time) ) then
    call restore_solution(solutions(uid),opt_solution_time)
  else
    call restore_solution(solutions(uid))
  endif
end subroutine

subroutine destroy_tenstream(lfinalizepetsc)
  logical,optional :: lfinalizepetsc
  logical :: lfinalize = .True.
  integer(iintegers) :: uid
  if(present(lfinalizepetsc)) lfinalize = lfinalizepetsc

  if(ltenstream_is_initialized) then
    if(linit_kspdir) then
      call KSPDestroy(kspdir , ierr) ;call CHKERR(ierr); linit_kspdir =.False.
    endif
    if(linit_kspdiff) then
      call KSPDestroy(kspdiff, ierr) ;call CHKERR(ierr); linit_kspdiff=.False.
    endif

    if(.not. ltwostr_only) then
      call VecDestroy(incSolar , ierr) ;call CHKERR(ierr)
      call VecDestroy(b        , ierr) ;call CHKERR(ierr)
      deallocate(incSolar)
      deallocate(b)
    endif
    call destroy_matrices()

    do uid=lbound(solutions,1),ubound(solutions,1)
        if( solutions(uid)%lset ) then
            if(solutions(uid)%lsolar_rad) then
                call VecDestroy(solutions(uid)%edir , ierr) ;call CHKERR(ierr)
                solutions(uid)%lsolar_rad = .False.
            endif

            call VecDestroy(solutions(uid)%ediff    , ierr) ;call CHKERR(ierr)
            call VecDestroy(solutions(uid)%abso     , ierr) ;call CHKERR(ierr)

            if(allocated(solutions(uid)%ksp_residual_history)) &
                deallocate(solutions(uid)%ksp_residual_history)

            solutions(uid)%lset = .False.
        endif
    enddo

    if(allocated(atm)) deallocate(atm)

    if(allocated(sun%angles)) deallocate(sun%angles)

    call OPP_1_2%destroy()
    call OPP_8_10%destroy()

    call DMDestroy(C_dir%da ,ierr); deallocate(C_dir )
    call DMDestroy(C_diff%da,ierr); deallocate(C_diff)
    call DMDestroy(C_one%da ,ierr); deallocate(C_one )
    call DMDestroy(C_one1%da,ierr); deallocate(C_one1)
    call DMDestroy(C_one_atm%da ,ierr); deallocate(C_one_atm)
    call DMDestroy(C_one_atm1%da,ierr); deallocate(C_one_atm1)

    ltenstream_is_initialized=.False.
    if(myid.eq.0 .and. ldebug)print *,'Destroyed TenStream'

    if(lfinalize) then
        call PetscFinalize(ierr) ;call CHKERR(ierr)
        deallocate(logstage)
        if(myid.eq.0 .and. ldebug)print *,'Finalized Petsc'
    endif
  endif
end subroutine

subroutine tenstream_get_result(redir,redn,reup,rabso, opt_solution_uid )
  real(ireals),dimension(:,:,:),intent(inout),allocatable :: redir
  real(ireals),dimension(:,:,:),intent(out)               :: redn,reup,rabso
  integer(iintegers),optional,intent(in) :: opt_solution_uid

  integer(iintegers) :: uid, lb_redir
  PetscScalar,pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()


  if(lenable_solutions_err_estimates .and. present(opt_solution_uid)) then
    uid = opt_solution_uid

    ! for residual history testing: old, un-updated solution should lie in solutions(-uid), return the old solution vector
    if( time_debug_solutions.gt.zero ) then ! two possibilities:

      ! first, we may be at a point where we want to overwrite the solution but keep everything else untouched
      if( (solutions(-uid)%lchanged.eqv..False.) .and. (solutions(uid)%time(1) .lt. solutions(-uid)%time(1)+time_debug_solutions) ) then  ! need consistent state in (-uid) and less time since last active solution than debug_sol_time
        print *,uid,'Getting Results from Vector -uid insted of ',time_debug_solutions,solutions(uid)%time(1),solutions(-uid)%time(1),'::', &
          solutions(-uid)%lchanged.eqv..False.,solutions(uid)%time(1) .lt. solutions(-uid)%time(1)+time_debug_solutions

        if(solutions(uid)%lsolar_rad) &
          call VecCopy(solutions(-uid)%edir , solutions(uid)%edir , ierr)
        call VecCopy(solutions(-uid)%ediff, solutions(uid)%ediff, ierr)
        call VecCopy(solutions(-uid)%abso , solutions(uid)%abso , ierr)


      else
        ! secondly, we really want to refresh the starting point -- here overwrite (-uid) with correct set:
        print *,uid,'Getting Results from normal Vector insted of (-uid) -- REFRESHING initial solution',time_debug_solutions,solutions(uid)%time(1),solutions(-uid)%time(1),'::', &
          solutions(-uid)%lchanged.eqv..False.,solutions(uid)%time(1) .lt. solutions(-uid)%time(1)+time_debug_solutions
        call copy_solution(solutions(uid),solutions(-uid))
      endif

    endif !time_debug_solutions

  else
    uid = i0 ! default solution is uid==0
  endif

  if(ldebug .and. myid.eq.0) print *,'calling tenstream_get_result',allocated(redir),'for uid',uid

  if(solutions(uid)%lchanged) stop 'tried to get results from unrestored solution -- call restore_solution first'

  if(allocated(redir)) then
    if( .not. solutions(uid)%lsolar_rad ) then
      redir = zero
    else
      if( solutions(uid)%lintegrated_dir ) stop 'tried to get result from integrated result vector(dir)'
      call getVecPointer(solutions(uid)%edir, C_dir, x1d, x4d)

      if(atm%lcollapse) then
          if(myid.eq.0 .and. ldebug) print *,'shape redir',shape(redir),lbound(redir,1),ubound(redir,1)
          lb_redir = lbound(redir,1)
          redir(lb_redir, :, :) = sum(x4d(i0:i3,C_dir%zs, :, :),dim=1)/4 ! return average of the 4 vertical tiles
          redir(lb_redir+1:atmk(C_dir%zs)+lb_redir, :, :) = zero
          redir(atmk(C_dir%zs)+1+lb_redir :C_one_atm1%ze+lb_redir, :, :) = sum(x4d(i0:i3,C_dir%zs+1:C_dir%ze,:,:),dim=1)/4 ! return average of the 4 vertical tiles
      else
          redir = sum(x4d(i0:i3, :, :, :),dim=1)/4
      endif
      if(ldebug) then
        if(myid.eq.0) print *,'Edir',redir(1,1,:)
        if(any(redir.lt.-one)) then 
          print *,'Found direct radiation smaller than 0 in dir result... that should not happen',minval(redir)
          call exit(1)
        endif
      endif
      call restoreVecPointer(solutions(uid)%edir,C_dir,x1d,x4d)
    endif
  endif

  if(solutions(uid)%lintegrated_diff) stop 'tried to get result from integrated result vector(diff)'
  call getVecPointer(solutions(uid)%ediff, C_diff, x1d, x4d)

  if(atm%lcollapse) then
      redn(1, :, :) = x4d(E_dn, C_diff%zs, :, :)
      reup(1, :, :) = x4d(E_up, C_diff%zs, :, :)
      redn(2:atmk(C_diff%zs)+1, :, :) = zero
      reup(2:atmk(C_diff%zs)+1, :, :) = zero
      redn(atmk(C_dir%zs)+2 :C_one_atm1%ze+1, :, :) = x4d(E_dn, C_diff%zs+1:C_diff%ze, :, :)
      reup(atmk(C_dir%zs)+2 :C_one_atm1%ze+1, :, :) = x4d(E_up, C_diff%zs+1:C_diff%ze, :, :)
  else
      redn = x4d(E_dn,:,:,:)
      reup = x4d(E_up,:,:,:)
  endif

  if(myid.eq.0 .and. ldebug) print *,'surface Edn',redn(ubound(redn,1), 1,1)
  if(myid.eq.0 .and. ldebug) print *,'surface Eup',reup(ubound(redn,1), 1,1)

  if(ldebug .and. solutions(uid)%lsolar_rad) then
    if(myid.eq.0) print *,' Edn',redn(1,1,:)
    if(myid.eq.0) print *,' Eup',reup(1,1,:)
    if(any(redn.lt.-one)) then 
      print *,'Found direct radiation smaller than 0 in edn result... that should not happen',minval(redn)
      call exit(1)
    endif
    if(any(reup.lt.-one)) then 
      print *,'Found direct radiation smaller than 0 in eup result... that should not happen',minval(reup)
      call exit(1)
    endif
  endif!ldebug
  call restoreVecPointer(solutions(uid)%ediff,C_diff,x1d,x4d)

  call getVecPointer(solutions(uid)%abso, C_one, x1d, x4d)
  if(atm%lcollapse) then
      rabso(1:atmk(C_one%zs)+1, :, :) = zero
      rabso(atmk(C_one%zs)+2 :C_one_atm%ze+1, :, :) = x4d(i0,C_one%zs+1:C_one%ze,:,:)
  else
      rabso = x4d(i0,:,:,:)
  endif
  call restoreVecPointer(solutions(uid)%abso,C_one,x1d,x4d)
end subroutine

      subroutine tenstream_get_result_toZero(res_edir,res_edn,res_eup,res_abso)
        ! after solving equations -- retrieve the results for edir,edn,eup and absorption
        ! only zeroth node gets the results back.

        real(ireals),intent(out),dimension(:,:,:) :: res_edir
        real(ireals),intent(out),dimension(:,:,:) :: res_edn
        real(ireals),intent(out),dimension(:,:,:) :: res_eup
        real(ireals),intent(out),dimension(:,:,:) :: res_abso

        real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso


        allocate( redir(C_one_atm1%zs:C_one_atm1%ze, C_dir%xs :C_dir%xe , C_dir%ys :C_dir%ye   )); redir=0
        allocate( redn (C_one_atm1%zs:C_one_atm1%ze, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye  )); redn =0
        allocate( reup (C_one_atm1%zs:C_one_atm1%ze, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye  )); reup =0
        allocate( rabso(C_one_atm%zs :C_one_atm%ze , C_one%xs :C_one%xe , C_one%ys :C_one%ye   )); rabso=0

        call tenstream_get_result(redir,redn,reup,rabso)

        call exchange_var(C_one_atm1, redir, res_edir)
        call exchange_var(C_one_atm1, redn , res_edn )
        call exchange_var(C_one_atm1, reup , res_eup )
        call exchange_var(C_one_atm , rabso, res_abso)

        if(myid.eq.0 .and. ldebug) then
          print *,'Retrieving results:',shape(res_edir)
          print *,sum(res_edir)/size(res_edir)
          print *,sum(res_edn) /size(res_edn)
          print *,sum(res_eup) /size(res_eup)
          print *,sum(res_abso)/size(res_abso)
        endif

        contains
            subroutine exchange_var(C, inp, outp)
                type(t_coord),intent(in) :: C
                real(ireals),intent(in) :: inp(:,:,:) ! local array from get_result
                real(ireals),intent(out) :: outp(:,:,:) ! global sized array on rank 0

                real(ireals),allocatable :: tmp(:,:,:,:)


                Vec :: vec
                PetscScalar,pointer,dimension(:,:,:,:) :: xinp=>null()
                PetscScalar,pointer,dimension(:) :: xinp1d=>null()

                call DMGetGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)
                call getVecPointer(vec ,C ,xinp1d, xinp)
                xinp(i0,:,:,:) = inp
                call restoreVecPointer(vec ,C ,xinp1d, xinp )

                call globalVec2Local(vec,C,tmp)

                call DMRestoreGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)

                if(myid.eq.0) outp = tmp(lbound(tmp,1), &
                                         lbound(tmp,2):lbound(tmp,2)+size(outp,1)-1,&
                                         lbound(tmp,3):lbound(tmp,3)+size(outp,2)-1,&
                                         lbound(tmp,4):lbound(tmp,4)+size(outp,3)-1)
            end subroutine

            subroutine globalVec2Local(vec,C,res)
                Vec :: vec
                real(ireals),allocatable :: res(:,:,:,:)
                type(t_coord) :: C

                Vec :: natural, local
                VecScatter :: scatter_context

                PetscScalar,Pointer :: xloc(:)

                if(allocated(res)) deallocate(res)
                if(myid.eq.0) allocate( res(C%dof,C%glob_zm,C%glob_xm,C%glob_ym) )

                call DMDACreateNaturalVector(C%da, natural, ierr); call CHKERR(ierr)

                call DMDAGlobalToNaturalBegin(C%da,vec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)
                call DMDAGlobalToNaturalEnd  (C%da,vec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)

                call VecScatterCreateToZero(natural, scatter_context, local, ierr); call CHKERR(ierr)

                call VecScatterBegin(scatter_context, natural, local, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
                call VecScatterEnd  (scatter_context, natural, local, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

                call VecScatterDestroy(scatter_context, ierr); call CHKERR(ierr)

                if(myid.eq.0) then
                    call VecGetArrayF90(local,xloc,ierr) ;call CHKERR(ierr)

                    res = reshape( xloc, (/ C%dof,C%glob_zm,C%glob_xm,C%glob_ym /) )

                    call VecRestoreArrayF90(local,xloc,ierr) ;call CHKERR(ierr)
                endif

                call VecDestroy(local,ierr); call CHKERR(ierr)
                call VecDestroy(natural,ierr); call CHKERR(ierr)
            end subroutine
        end subroutine

subroutine prepare_solution(solution, uid, lsolar)
  type(t_state_container) :: solution
  integer(iintegers),intent(in) :: uid
  logical,intent(in) :: lsolar

  solutions(uid)%uid = uid ! dirty hack to give the solution a unique hash for example to write it out to disk

  if( lsolar .and.  (solution%lsolar_rad.eqv..False.) ) then 
    ! set up the direct vectors in any case. This may be necessary even if solution was
    ! already initialized once... e.g. in case we calculated thermal before with same uid
    call DMCreateGlobalVector(C_dir%da ,solution%edir  ,ierr)  ; call CHKERR(ierr)
    call VecSet(solution%edir,zero,ierr); call CHKERR(ierr)
    solution%lsolar_rad = .True.
  endif

  if(ldebug .and. myid.eq.0) print *,'Prepare solution',uid,'::',lsolar,solution%lsolar_rad,'set?',solution%lset

  if(solution%lset) return ! already set-up 

  call DMCreateGlobalVector(C_diff%da,solution%ediff ,ierr)  ; call CHKERR(ierr)
  call DMCreateGlobalVector(C_one%da ,solution%abso  ,ierr)  ; call CHKERR(ierr)

  call VecSet(solution%ediff,zero,ierr)    ; call CHKERR(ierr)
  call VecSet(solution%abso,zero,ierr)     ; call CHKERR(ierr)

  solution%lset = .True.
end subroutine

function need_new_solution(uid,time)
  integer(iintegers),intent(in) :: uid
  real(ireals),intent(in),optional :: time
  logical :: need_new_solution

  integer,parameter :: Nfit=3   ! Number of used residuals
  integer,parameter :: Nporder=2 ! max order of polynomial
  real(ireals) :: t(Nfit),tm(Nfit),dt(Nfit-1),err(2, 2*(Nfit-1)), error_estimate
  real(ireals) :: polyc(Nporder+1),estimate(Nporder)

  character(default_str_len) :: reason
  integer, parameter :: out_unit=20

  integer(iintegers) :: k,ipoly

  ! Make time an optional argument here for
  ! convenience of the interface --
  ! otherwise the user needs to check if he
  ! opt_time is present and so on...
  if(.not. present(time)) then
    need_new_solution = .True.
    return
  endif

  if( .not. solutions(uid)%lset ) then !if we did not store a solution, return immediately
    need_new_solution=.True.
    write(reason,*) 'no solution yet'
    if(ldebug .and. myid.eq.0) print *,'new calc',need_new_solution,' bc ',reason,' t',time,uid
    return
  endif

  if(.not. lenable_solutions_err_estimates) then
    need_new_solution=.True.
    write(reason,*) 'err.est.thresh.inf.small'
    if(ldebug .and. myid.eq.0) print *,'new calc',need_new_solution,' bc ',reason,' t',time,uid
    return
  endif


  call PetscLogStagePush(logstage(11),ierr) ;call CHKERR(ierr)
  do k=1,Nfit
    t(k) = solutions(uid)%time(Nfit-k+1)
  enddo

  !            call parabola(t2, e1, t3, e2, t4, e3, time, error_estimate)
  !            call exponential(t3, e2, t4, e3, time, error_estimate)

  ! t is pts where the solution got updated
  ! tm is in between those time
  ! dt is the weight for the integral
  do k=1,Nfit-1
    tm(k) = (t(k)+t(k+1) )*.5_ireals
  enddo
  tm(Nfit) = ( t(Nfit) + time )*.5_ireals

  do k=1,Nfit-1
    dt(k) = tm(k+1) - tm(k)
  enddo

  ! setup error timeseries on which to fit -- first index is time, second is error
  do k=1,Nfit-1
    err(1, 2*k-1) = tm(k)
    err(1, 2*k  ) = t(k+1)
    err(2, 2*k-1) = zero
    err(2, 2*k  ) = solutions(uid)%maxnorm(Nfit-k)
  enddo

  ! try several polynomials and find max error:
  do ipoly=1,Nporder
    polyc(1:ipoly+1) = polyfit(err(1,:),err(2,:),ipoly, ierr) ! e.g. second order polynomial has 3 coefficients
    if(ierr.ne.0) then
      need_new_solution=.True.
      write(reason,*) 'problem fitting error curve',ierr
      call PetscLogStagePop(ierr) ;call CHKERR(ierr)
      if(ldebug .and. myid.eq.0) print *,'new calc',need_new_solution,' bc ',reason,' t',time,uid
      return
    endif
    estimate(ipoly)=zero
    do k=1,ipoly+1
      estimate(ipoly) = estimate(ipoly) + polyc(k)*time**(k-1)
    enddo
  enddo


  error_estimate = maxval(abs(estimate))

  if(myid.eq.0 .and. ldebug) then ! .and. error_estimate.le.zero) then
    print *,'DEBUG t',t
    print *,'DEBUG tm',tm
    print *,'DEBUG dt',dt
    print *,'DEBUG err',err(1,:),'::',err(2,:)
    print *,'DEBUG err_est',error_estimate,'::',estimate
    print *,'DEBUG polyc',polyc
  endif

  if(error_estimate.le.options_max_solution_err) then
    need_new_solution=.False.
    write(reason,*) 'ERR_TOL_IN_BOUND'
  else
    need_new_solution=.True.
    write(reason,*) 'ERR_TOL_EXCEEDED'
  endif

  if(any(t.lt.zero) ) then
    need_new_solution=.True.
    write(reason,*) 'FEW_SOLUTIONS'
  endif

  if(time-solutions(uid)%time(1) .gt. options_max_solution_time) then
    need_new_solution=.True.
    write(reason,*) 'MIN_TIME_EXCEEDED'
  endif

  if(time_debug_solutions.gt.zero) then
    if(.not.need_new_solution) then
      need_new_solution=.True. ! overwrite it and calculate anyway
      write(reason,*) 'MANUAL OVERRIDE'
      ! Hack to monitor error growth...
      ! We tell the user that he has to calculate radiation again.
      ! We will calculate and update the solution vectors...
      ! But once he tries to load em, we give him the old, un-updated solution.
      ! This happens in restore_solution i.e. we overwrite solution with old one.
    endif
  endif

  !            if(ldebug .and. myid.eq.0 .and. .not. need_new_solution ) &
  if(ldebug.and. myid.eq.0) then
    print *,''
    print *,''
    print *,'new calc',need_new_solution,' bc ',reason,' t',time,uid,'    ::     est.',error_estimate,'[W]',error_estimate*86.1,'[K/d]'
    if(allocated(solutions(uid)%ksp_residual_history) ) &
      print *,' residuals _solver ::', solutions(uid)%ksp_residual_history(1:4)
    if(uid.eq.501) then
      open (unit=out_unit,file="residuals.log",action="write",status="replace")
    else
      open (unit=out_unit,file="residuals.log",action="readwrite",status="unknown",position = "append")
    endif

    write (out_unit,*) uid,solutions(uid)%maxnorm
    write (out_unit,*) uid,solutions(uid)%time
    write (out_unit,*) uid,solutions(uid)%twonorm
    close (out_unit)
  endif
  call PetscLogStagePop(ierr) ;call CHKERR(ierr)

contains
  subroutine exponential(x1, y1, x2, y2, x3, y3)
    ! fit to the function y = exp(alpha * x) +beta
    real(ireals), intent(in) :: x1, y1, x2, y2, x3
    real(ireals), intent(out) :: y3
    real(ireals) :: alpha,beta
    beta = y1 - one
    alpha = log( max(y1,y2) -beta)/(x2-x1)
    y3 = exp( alpha * (x3-x1) ) + beta
    !                  if(myid.eq.0) print *,'exponential error_estimate:',x1,y1,x2,y2,'::',alpha,beta,'::',x3,y3
    !                  if(myid.eq.0) print *,''
  end subroutine
  subroutine parabola(x1, y1, x2, y2, x3, y3, x4, y4)
    ! Solve for coefficient in equation A.x**2 + B.x + C and evaluate polynomial at x4
    real(ireals), intent(in) :: x1, y1, x2, y2, x3, y3, x4
    real(ireals), intent(out) :: y4
    real(ireals) :: denom,A,B,C
    denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
    A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
    B     = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
    C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + &
      x1 * x2 * (x1 - x2) * y3) / denom
    y4 = x4*(A*x4+B)+C
    !                  print *,'parabola:',denom,A,B,C,'::',y4
  end subroutine
  function polyfit(vx, vy, d, ierr) !Rosetta Code http://rosettacode.org/wiki/Polynomial_regression#Fortran
    implicit none
    integer(iintegers), intent(in)            :: d
    real(ireals), dimension(d+1)              :: polyfit
    real(ireals), dimension(:), intent(in)    :: vx, vy
    PetscErrorCode,intent(out) :: ierr

    real(ireals), dimension(size(vx),d+1) :: X
    real(ireals), dimension(d+1,size(vx)) :: XT
    real(ireals), dimension(d+1,d+1)      :: XTX

    integer(iintegers) :: i, j

    integer :: n, lda, lwork
    integer :: info
    integer, dimension(d+1) :: ipiv
    real(ireals)      , dimension(d+1) :: work

    ierr=0

    n = d+1
    lda = n
    lwork = n

    do i=1,size(vx)
      if(any (approx( vx(i), vx(i+1:size(vx)) ) ) ) then ! polyfit cannot cope with same x values --> matrix gets singular
        if(ldebug .and. myid.eq.0) print *,'polyfit cannot cope with same x values --> matrix gets singular',vx,'::',vy
        polyfit=0
        polyfit(1) = nil
        ierr=1
        return
      endif
    enddo

    ! prepare the matrix
    do i = 0, d
      do j = 1, size(vx)
        X(j, i+1) = vx(j)**i
      end do
    end do

    XT  = transpose(X)
    XTX = matmul(XT, X)

    ! calls to LAPACK subs DGETRF and DGETRI
    if(sizeof(one).eq.4) then !single precision
      call SGETRF(n, n, XTX, lda, ipiv, info)
    else if(sizeof(one).eq.8) then !double precision
      call DGETRF(n, n, XTX, lda, ipiv, info)
    else
      print *,'dont know which lapack routine to call for reals with sizeof==',sizeof(one)
    endif
    if ( info /= 0 ) then
      if(myid.eq.0 .and. ldebug) print *, "problem with lapack lsqr :: 1 :: info",info
      ierr=2
      return
    end if
    if(sizeof(one).eq.4) then !single precision
      call SGETRI(n, XTX, lda, ipiv, work, lwork, info)
    else if(sizeof(one).eq.8) then !double precision
      call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    else
      print *,'dont know which lapack routine to call for reals with sizeof==',sizeof(one)
    endif
    if ( info /= 0 ) then
      if(myid.eq.0) print *, "problem with lapack lsqr :: 2 :: info",info
      ierr=3
      return
    end if

    polyfit = matmul( matmul(XTX, XT), vy)
  end function
end function

subroutine copy_solution(inp,out)
  type(t_state_container),intent(in) :: inp
  type(t_state_container),intent(out) :: out


  if(ldebug.and.myid.eq.0) &
    print *,'DEBUG: Copy solution state for',inp%uid,'lchanged',inp%lchanged

  ! make sure target solution is defined.
  call prepare_solution(out, inp%uid, inp%lsolar_rad)

  out%lintegrated_dir = inp%lintegrated_dir
  out%lintegrated_diff= inp%lintegrated_diff
  out%lchanged = inp%lchanged

  out%time    = inp%time
  out%maxnorm = inp%maxnorm
  out%twonorm = inp%twonorm
  if(.not. allocated(out%ksp_residual_history)) allocate(out%ksp_residual_history( size(inp%ksp_residual_history) ) )
  out%ksp_residual_history = inp%ksp_residual_history

  if(out%lsolar_rad) call VecCopy(inp%edir , out%edir , ierr)
  call VecCopy(inp%ediff, out%ediff, ierr)
  call VecCopy(inp%abso , out%abso , ierr)

end subroutine

subroutine restore_solution(solution,time)
  ! restore_solution:: if flux have changed, we need to update absorption, save the residual history
  type(t_state_container) :: solution
  real(ireals),intent(in),optional :: time

  character(default_str_len) :: vecname
  real(ireals) :: norm1,norm2,norm3
  Vec :: abso_old

  call PetscLogStagePush(logstage(11),ierr) ;call CHKERR(ierr)

  if( .not. solution%lset ) &
    stop 'cant restore solution that was not initialized'

  if( .not. solution%lchanged ) &
    stop 'cant restore solution which was not changed'

  if(present(time) .and. lenable_solutions_err_estimates) then ! Create working vec to determine difference between old and new absorption vec
    call DMGetGlobalVector(C_one%da, abso_old, ierr) ; call CHKERR(ierr)
    call VecCopy( solution%abso, abso_old, ierr)     ; call CHKERR(ierr)
  endif

  ! make sure to bring the fluxes into [W] for absorption calculation
  call scale_flx(solution, lWm2_to_W=.True. )

  ! update absorption
  call calc_flx_div(solution)

  if(ldebug .and. myid.eq.0) &
    print *,'Saving Solution ',solution%uid

  ! make sure to bring the fluxes into [W/m**2]
  call scale_flx(solution, lWm2_to_W=.False. )

  if(ldebug .and. myid.eq.0) &
    print *,'Saving Solution done'
  solution%lchanged=.False.

  if(present(time) .and. lenable_solutions_err_estimates) then ! Compute norm between old absorption and new one
    call VecAXPY(abso_old , -one, solution%abso , ierr)    ; call CHKERR(ierr) ! overwrite abso_old with difference to new one
    call VecNorm(abso_old ,  NORM_1, norm1, ierr)          ; call CHKERR(ierr)
    call VecNorm(abso_old ,  NORM_2, norm2, ierr)          ; call CHKERR(ierr)
    call VecNorm(abso_old ,  NORM_INFINITY, norm3, ierr)   ; call CHKERR(ierr)

    call DMRestoreGlobalVector(C_one%da, abso_old, ierr)   ; call CHKERR(ierr)

    ! Save norm for later analysis
    solution%maxnorm = eoshift ( solution%maxnorm, shift = -1) !shift all values by 1 to the right
    solution%twonorm = eoshift ( solution%twonorm, shift = -1) !shift all values by 1 to the right
    solution%time    = eoshift ( solution%time   , shift = -1) !shift all values by 1 to the right

    solution%maxnorm( 1 ) = norm3
    solution%twonorm( 1 ) = norm2
    solution%time( 1 )    = time

    if(ldebug .and. myid.eq.0) &
      print *,'Updating error statistics for solution ',solution%uid,'at time ',time,'::',solution%time(1),':: norm',norm1,norm2,norm3,'[W] :: hr_norm approx:',norm3*86.1,'[K/d]'

  endif !present(time) .and. lenable_solutions_err_estimates


  call PetscLogStagePop(ierr) ;call CHKERR(ierr)

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

  if(allocated(b)) then
    write(vecname,FMT='("b",I0)') solution%uid
    call PetscObjectSetName(b,vecname,ierr) ; call CHKERR(ierr)
    call PetscObjectViewFromOptions(b, PETSC_NULL_VEC, "-show_b", ierr); call CHKERR(ierr)
  endif

  if(allocated(incSolar)) then
    write(vecname,FMT='("incSolar",I0)') solution%uid
    call PetscObjectSetName(incSolar,vecname,ierr) ; call CHKERR(ierr)
    call PetscObjectViewFromOptions(incSolar, PETSC_NULL_VEC, "-show_incSolar", ierr); call CHKERR(ierr)
  endif
end subroutine

subroutine getVecPointer(vec,C,x1d,x4d)
  Vec :: vec
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
  Vec :: vec
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

function get_mem_footprint()
  real(ireals) :: get_mem_footprint
  PetscLogDouble :: memory_footprint, petsc_current_mem
  get_mem_footprint = zero

  call mpi_barrier(imp_comm, ierr)
  call PetscMemoryGetCurrentUsage(memory_footprint, ierr); call CHKERR(ierr)

  get_mem_footprint = memory_footprint / 1024. / 1024. / 1024.

!  call PetscMallocGetCurrentUsage(petsc_current_mem, ierr); call CHKERR(ierr)

!  if(ldebug) print *,myid,'Memory Footprint',memory_footprint, 'B', get_mem_footprint, 'G'
  return
end function
end module
