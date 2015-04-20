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

#include "finclude/petscdef.h"
  use petsc
  use m_data_parameters, only : ireals,iintegers,       &
      imp_comm, myid, numnodes,init_mpi_data_parameters,mpiint, &
      zero,one,nil,i0,i1,i2,i3,i4,i5,i6,i7,i8,i10,pi

  use m_twostream, only: delta_eddington_twostream
  use m_helper_functions, only: deg2rad,approx,rmse,delta_scale,imp_bcast,cumsum,inc
  use m_eddington, only : eddington_coeff_zdun
  use m_optprop_parameters, only : ldelta_scale
  use m_optprop, only : t_optprop_1_2,t_optprop_8_10
  use m_tenstream_options, only : read_commandline_options, ltwostr, luse_eddington, twostr_ratio, &
      options_max_solution_err, options_max_solution_time, ltwostr_only, luse_twostr_guess,        &
      lwriteall,lcalc_nca

  implicit none

  private
  public :: init_tenstream, set_global_optical_properties, set_optical_properties, solve_tenstream, destroy_tenstream,&
      getVecPointer,restoreVecPointer, &
      tenstream_get_result, need_new_solution, &
      t_coord,C_dir,C_diff,C_one

  PetscInt,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9

  logical,parameter :: ldebug=.False.
!  logical,parameter :: ldebug=.True.
  logical,parameter :: lcycle_dir=.True.
  logical,parameter :: lprealloc=.True.

  type t_coord
    PetscInt :: xs,xe,ys              ! local domain start and end indices
    PetscInt :: ye,zs,ze              ! 
    PetscInt :: xm,ym,zm              ! size of local domain
    PetscInt :: gxs,gys,gzs           ! domain indices including ghost points
    PetscInt :: gxe,gye,gze           ! 
    PetscInt :: gxm,gym,gzm           ! size of local domain including ghosts
    PetscInt :: glob_xm,glob_ym,glob_zm ! global domain size
    PetscInt :: dof,dim               ! degrees of freedom of Petsc Domain, dimension of dmda
    DM :: da                          ! The Domain Decomposition Object
    PetscMPIInt,allocatable :: neighbors(:) ! all 3d neighbours( (x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
  end type

  type(t_coord),save :: C_dir,C_diff,C_one,C_one1

  PetscErrorCode :: ierr

  type t_optprop
    real(ireals) :: kabs,ksca,g
  end type

  type t_atmosphere
    type(t_optprop) , allocatable , dimension(:,:,:) :: op
    type(t_optprop) , allocatable , dimension(:,:,:) :: delta_op
    real(ireals)    , allocatable , dimension(:,:,:) :: planck
    real(ireals)    , allocatable , dimension(:,:,:) :: a11, a12, a13, a23, a33
    real(ireals)    , allocatable , dimension(:,:,:) :: g1,g2
    real(ireals)    , allocatable , dimension(:,:,:) :: dz
    logical         , allocatable , dimension(:,:,:) :: l1d
    real(ireals) :: albedo
    real(ireals) :: dx,dy
  end type
  type(t_atmosphere),save :: atm


  type(t_optprop_1_2),save  :: OPP_1_2
  type(t_optprop_8_10),save :: OPP_8_10

  type t_suninfo
    real(ireals) :: symmetry_phi
    integer(iintegers) :: yinc,xinc
    real(ireals) :: theta,phi,costheta,sintheta
  end type
  type(t_suninfo),save :: sun

  PetscLogStage,save :: logstage(12)

  Mat,save :: Mdir,Mdiff

  Vec,save :: incSolar,b

  KSP,save :: kspdir, kspdiff
  logical,save :: linit_kspdir=.False., linit_kspdiff=.False.

  logical,save :: linitialized=.False.

  integer(iintegers),parameter :: minimal_dimension=3 ! this is the minimum number of gridpoints in x or y direction

  logical,parameter :: lenable_solutions=.True. ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
  type t_state_container
    integer(iintegers) :: uid ! dirty hack to give the solution a unique hash for example to write it out to disk -- this should be the same as the index in global solutions array
    Vec :: edir,ediff,abso
    logical :: lset        = .False. ! initialized?
    logical :: lsolar_rad  = .False. ! direct radiation calculated?
    logical :: lchanged    = .True.  ! did the flux change recently? -- call restore_solution to bring it in a coherent state

    ! save state of solution vectors... they are either in [W](true) or [W/m**2](false)
    logical :: lintegrated_dir=.True. , lintegrated_diff=.True.  

    !save error statistics
    real(ireals) :: time   (300) = -one
    real(ireals) :: maxnorm(300) = zero
    real(ireals) :: twonorm(300) = zero
    real(ireals),allocatable :: ksp_residual_history(:)
  end type
  type(t_state_container),save :: solutions(-1:1000)

contains 
  subroutine setup_grid(Nz,Nx,Ny,nxproc,nyproc)
      PetscInt,intent(in) :: Nz,Nx,Ny
      integer(iintegers),optional :: nxproc(:), nyproc(:) ! size of local domains on each node

      DMBoundaryType :: bp=DM_BOUNDARY_PERIODIC, bn=DM_BOUNDARY_NONE, bg=DM_BOUNDARY_GHOSTED

      if(myid.eq.0.and.ldebug) print *,myid,'Setting up the DMDA grid for ',Nz,Nx,Ny,'using ',numnodes,' nodes'

      if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_diff'
      call setup_dmda(C_diff, Nz+1,Nx,Ny, bp, i10)

      if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_dir'
      if(lcycle_dir) then
        call setup_dmda(C_dir, Nz+1,Nx,Ny, bp, i8)
      else
        call setup_dmda(C_dir, Nz+1,Nx,Ny, bg, i8)
      endif

      if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C1'
      call setup_dmda(C_one , Nz  , Nx,Ny,  bp, i1)
      call setup_dmda(C_one1, Nz+1, Nx,Ny,  bp, i1)

      if(myid.eq.0.and.ldebug) print *,myid,'DMDA grid ready'
    contains
      subroutine setup_dmda(C, Nz, Nx, Ny, boundary, dof)
          type(t_coord) :: C
          PetscInt,intent(in) :: Nz,Nx,Ny,dof
          DMBoundaryType,intent(in) :: boundary

          PetscInt,parameter :: stencil_size=1

          C%dof = i1*dof
          if(present(nxproc) .and. present(nyproc) ) then
            call DMDACreate3d( imp_comm ,                                       &
                bn                      , boundary             , boundary     , &
                DMDA_STENCIL_STAR       ,                                       &
                i1*Nz                   , sum(nxproc)          , sum(nyproc)  , &
                i1                      , size(nxproc)         , size(nyproc) , &
                C%dof                   , stencil_size         ,                &
                Nz                      , nxproc               , nyproc       , &
                C%da                    , ierr) ;CHKERRQ(ierr)
          else
            call DMDACreate3d( imp_comm ,                                             &
                bn                      , boundary             , boundary           , &
                DMDA_STENCIL_STAR       ,                                             &
                i1*Nz                   , Nx                   , Ny                 , &
                i1                      , PETSC_DECIDE         , PETSC_DECIDE       , &
                C%dof                   , stencil_size         ,                      &
                Nz                      , PETSC_NULL_INTEGER   , PETSC_NULL_INTEGER , &
                C%da                    , ierr) ;CHKERRQ(ierr)
          endif

          call DMSetup(C%da,ierr) ;CHKERRQ(ierr)
          call DMSetMatType(C%da, MATAIJ, ierr); CHKERRQ(ierr)
          if(lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ;CHKERRQ(ierr)

          call DMSetFromOptions(C%da, ierr) ; CHKERRQ(ierr)
          if(ldebug) call DMView(C%da, PETSC_VIEWER_STDOUT_WORLD ,ierr)
          call setup_coords(C)
      end subroutine
      subroutine setup_coords(C)
          type(t_coord) :: C

          call DMDAGetInfo(C%da,C%dim,                             &
              C%glob_zm,C%glob_xm,C%glob_ym,                           &
              PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
              PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
              PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
              ierr) ;CHKERRQ(ierr)

          call DMDAGetCorners(C%da, C%zs, C%xs,C%ys, C%zm, C%xm,C%ym, ierr) ;CHKERRQ(ierr)
          C%xe = C%xs+C%xm-1
          C%ye = C%ys+C%ym-1
          C%ze = C%zs+C%zm-1

          call DMDAGetGhostCorners(C%da,C%gzs,C%gxs,C%gys,C%gzm,C%gxm,C%gym,ierr) ;CHKERRQ(ierr)
          C%gxe = C%gxs+C%gxm-1
          C%gye = C%gys+C%gym-1
          C%gze = C%gzs+C%gzm-1

          if(ldebug) then
            print *,myid,'Domain Corners z:: ',C%zs,':',C%ze,' (',C%zm,' entries)','global size',C%glob_zm
            print *,myid,'Domain Corners x:: ',C%xs,':',C%xe,' (',C%xm,' entries)','global size',C%glob_xm
            print *,myid,'Domain Corners y:: ',C%ys,':',C%ye,' (',C%ym,' entries)','global size',C%glob_ym
          endif

          allocate(C%neighbors(0:3**C%dim-1) )
          call DMDAGetNeighbors(C%da,C%neighbors,ierr) ;CHKERRQ(ierr)
          if(ldebug.and.C%dim.eq.3) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([10,4,16,22]),'while I am ',C%neighbors(13)
          if(ldebug.and.C%dim.eq.2) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([1,3,7,5]),'while I am ',C%neighbors(4)
      end subroutine
  end subroutine

  subroutine mat_info(A)
      Mat :: A
      MatInfo :: info(MAT_INFO_SIZE) 
      real(ireals) :: mal, nz_allocated, nz_used, nz_unneeded
      PetscInt :: m,n

      return !TODO decide if this is the correct type for matgetinfo call, in docs it states we should use double precision? for the time being it may be safe to just return....
      call MatGetInfo(A,MAT_LOCAL,info,ierr) ;CHKERRQ(ierr)
      mal = info(MAT_INFO_MALLOCS)
      nz_allocated = info(MAT_INFO_NZ_ALLOCATED)
      nz_used   = info(MAT_INFO_NZ_USED)
      nz_unneeded = info(MAT_INFO_NZ_UNNEEDED)

      call MatGetOwnershipRange(A, m,n, ierr)

      if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_MALLOCS',mal,'MAT_INFO_NZ_ALLOCATED',nz_allocated
      if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_USED',nz_used,'MAT_INFO_NZ_unneded',nz_unneeded
      if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: Ownership range',m,n

  end subroutine
  subroutine init_Matrix(A,C)!,prefix)
      Mat :: A
      type(t_coord) :: C
      !        character(len=*),optional :: prefix

      PetscInt,dimension(:),allocatable :: o_nnz,d_nnz!,dnz

      call DMCreateMatrix(C%da, A, ierr) ;CHKERRQ(ierr)

      if(lprealloc) then
        ! Determine perfect preallocation
        if(numnodes.gt.1) then
          select case (C%dof)
          case(i3)
            call setup_dir_preallocation(d_nnz,o_nnz,C)
          case(i8)
!            call setup_dir_preallocation(d_nnz,o_nnz,C)
            call setup_dir8_preallocation(d_nnz,o_nnz,C)
          case(i10)
            call setup_diff10_preallocation(d_nnz,o_nnz,C)
          case default
            stop('Dont know which preallocation routine I shall call! - exiting...')
          end select

          call MatMPIAIJSetPreallocation(A, PETSC_NULL_INTEGER,d_nnz, PETSC_NULL_INTEGER, o_nnz, ierr) ;CHKERRQ(ierr)

          deallocate(o_nnz)
          deallocate(d_nnz)

        endif
        call MatSeqAIJSetPreallocation(A, C%dof+i1, PETSC_NULL_INTEGER, ierr) ;CHKERRQ(ierr)
      endif

      call mat_info(A)

      !        call MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr) ;CHKERRQ(ierr)
      !        call MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr) ;CHKERRQ(ierr)
      !        call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr) ;CHKERRQ(ierr)

      call MatSetFromOptions(A,ierr) ;CHKERRQ(ierr)
      call MatSetUp(A,ierr) ;CHKERRQ(ierr)

      call mat_info(A)
      call mat_set_diagonal(A,C)
  end subroutine
  subroutine mat_set_diagonal(A,C)
      Mat :: A
      type(t_coord),intent(in) :: C
      PetscInt :: i,j,k,dof
      MatStencil :: row(4,1), col(4,1)
      PetscScalar :: v(1)

      !TODO -- we should use this form... however this does somehow corrupt preallocation? - maybe fix this
      !        Vec :: diag
      !        call DMCreateGlobalVector(C%da,diag,ierr) ;CHKERRQ(ierr)
      !        call VecSet(diag,one,ierr) ;CHKERRQ(ierr)
      !        call MatDiagonalSet( A, diag, INSERT_VALUES,ierr) ;CHKERRQ(ierr)
      !        call VecDestroy(diag,ierr) ;CHKERRQ(ierr)

      v(1)=one

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

              call MatSetValuesStencil(A,i1, row,i1, col , v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr) 
            enddo 
          enddo 
        enddo 
      enddo

      call MatAssemblyBegin(A,MAT_FLUSH_ASSEMBLY,ierr) ;CHKERRQ(ierr)
      call MatAssemblyEnd  (A,MAT_FLUSH_ASSEMBLY,ierr) ;CHKERRQ(ierr)  

      if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients diagonally ... done'
  end subroutine
  subroutine setup_diff10_preallocation(d_nnz,o_nnz,C)
      PetscInt,allocatable :: o_nnz(:)
      PetscInt,allocatable :: d_nnz(:)
      type(t_coord) :: C
      Vec :: v_o_nnz,v_d_nnz
      PetscScalar,Pointer :: xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
      PetscScalar,Pointer :: xo1d(:)=>null(),xd1d(:)=>null()

      PetscInt :: vsize

      PetscInt,parameter :: ind(9)=[E_up,E_le_m,E_le_p,E_ri_m,E_ri_p,E_ba_m,E_ba_p,E_fw_m,E_fw_p]
      PetscInt :: i,j,k!d,li,lj,lk

      call DMGetGlobalVector(C%da,v_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMGetGlobalVector(C%da,v_d_nnz,ierr) ;CHKERRQ(ierr)

      call VecGetLocalSize(v_o_nnz,vsize,ierr) ;CHKERRQ(ierr)
      allocate(o_nnz(0:vsize-1))
      allocate(d_nnz(0:vsize-1))

      call getVecPointer(v_o_nnz,C,xo1d,xo)
      call getVecPointer(v_d_nnz,C,xd1d,xd)

      xo = i0
      xd = i1
      xd(E_up,C%ze,:,:) = i2 ! E_up at ground is only function of E_dn at grnd and +1 for self

      xd(ind ,C%zs:C%ze-1,:,:) = C%dof+i1
      xd(E_dn,C%zs+1:C%ze,:,:) = C%dof+i1

      ! Determine prealloc for E_dn and E_up
      if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh east
        ! E_le is not local
        xd(E_dn,C%zs+1:C%ze,C%xe,:) = xd(E_dn,C%zs+1:C%ze,C%xe,:) - i2
        xo(E_dn,C%zs+1:C%ze,C%xe,:) = xo(E_dn,C%zs+1:C%ze,C%xe,:) + i2
        xd(E_up,C%zs:C%ze-1,C%xe,:) = xd(E_up,C%zs:C%ze-1,C%xe,:) - i2
        xo(E_up,C%zs:C%ze-1,C%xe,:) = xo(E_up,C%zs:C%ze-1,C%xe,:) + i2
      endif
      if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0) then ! neigh north
        ! E_ba is not local
        xd(E_dn,C%zs+1:C%ze,:,C%ye) = xd(E_dn,C%zs+1:C%ze,:,C%ye) - i2
        xo(E_dn,C%zs+1:C%ze,:,C%ye) = xo(E_dn,C%zs+1:C%ze,:,C%ye) + i2
        xd(E_up,C%zs:C%ze-1,:,C%ye) = xd(E_up,C%zs:C%ze-1,:,C%ye) - i2
        xo(E_up,C%zs:C%ze-1,:,C%ye) = xo(E_up,C%zs:C%ze-1,:,C%ye) + i2
      endif
      if( C%neighbors(4).ne.myid .and. C%neighbors(4).ge.i0) then ! neigh south
        ! no foreign stream dependencies
      endif
      if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0) then ! neigh west
        ! no foreign stream dependencies
      endif
      ! Determine prealloc for E_le
      if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh east
        ! E_le is not local
        xd([E_le_m,E_le_p],C%zs:C%ze-1,C%xe,:) = xd([E_le_m,E_le_p],C%zs:C%ze-1,C%xe,:) - i2
        xo([E_le_m,E_le_p],C%zs:C%ze-1,C%xe,:) = xo([E_le_m,E_le_p],C%zs:C%ze-1,C%xe,:) + i2
      endif
      if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 ) then ! neigh north
        ! E_ba is not local
        xd([E_le_m,E_le_p],C%zs:C%ze-1,:,C%ye) = xd([E_le_m,E_le_p],C%zs:C%ze-1,:,C%ye) - i2
        xo([E_le_m,E_le_p],C%zs:C%ze-1,:,C%ye) = xo([E_le_m,E_le_p],C%zs:C%ze-1,:,C%ye) + i2
      endif
      if( C%neighbors(4).ne.myid .and. C%neighbors(4).ge.i0 ) then ! neigh south
        ! no foreign stream dependencies
      endif
      if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh west
        ! no foreign stream dependencies
      endif
      ! Determine prealloc for E_ri
      if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh east
        ! no foreign stream dependencies
      endif
      if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 ) then ! neigh north
        ! E_ba is not local
        xd([E_ri_m,E_ri_p],C%zs:C%ze-1,:,C%ye) = xd([E_ri_m,E_ri_p],C%zs:C%ze-1,:,C%ye) - i2
        xo([E_ri_m,E_ri_p],C%zs:C%ze-1,:,C%ye) = xo([E_ri_m,E_ri_p],C%zs:C%ze-1,:,C%ye) + i2
      endif
      if( C%neighbors(4).ne.myid .and. C%neighbors(4).ge.i0 ) then ! neigh south
        ! no foreign stream dependencies
      endif
      if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh west
        ! E_ri local dependencies are only self, and 2*E_le
        xd([E_ri_m,E_ri_p],C%zs:C%ze-1,C%xs,:) = i3
        xo([E_ri_m,E_ri_p],C%zs:C%ze-1,C%xs,:) = i8
      endif
      ! Determine prealloc for E_ba
      if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh east
        ! E_le is not local
        xd([E_ba_m,E_ba_p],C%zs:C%ze-1,C%xe,:) = xd([E_ba_m,E_ba_p],C%zs:C%ze-1,C%xe,:) - i2
        xo([E_ba_m,E_ba_p],C%zs:C%ze-1,C%xe,:) = xo([E_ba_m,E_ba_p],C%zs:C%ze-1,C%xe,:) + i2
      endif
      if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 ) then ! neigh north
        ! E_ba is not local
        xd([E_ba_m,E_ba_p],C%zs:C%ze-1,:,C%ye) = xd([E_ba_m,E_ba_p],C%zs:C%ze-1,:,C%ye) - i2
        xo([E_ba_m,E_ba_p],C%zs:C%ze-1,:,C%ye) = xo([E_ba_m,E_ba_p],C%zs:C%ze-1,:,C%ye) + i2
      endif
      if( C%neighbors(4).ne.myid .and. C%neighbors(4).ge.i0 ) then ! neigh south
        ! no foreign stream dependencies
      endif
      if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh west
        ! no foreign stream dependencies
      endif
      ! Determine prealloc for E_fw
      if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh east
        ! E_le is not local
        xd([E_fw_m,E_fw_p],C%zs:C%ze-1,C%xe,:) = xd([E_fw_m,E_fw_p],C%zs:C%ze-1,C%xe,:) - i2
        xo([E_fw_m,E_fw_p],C%zs:C%ze-1,C%xe,:) = xo([E_fw_m,E_fw_p],C%zs:C%ze-1,C%xe,:) + i2
      endif
      if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 ) then ! neigh north
        ! no foreign stream dependencies
      endif
      if( C%neighbors(4).ne.myid .and. C%neighbors(4).ge.i0 ) then ! neigh south
        ! E_fw local dependencies are only self, and 2*E_ba
        xd([E_fw_m,E_fw_p],C%zs:C%ze-1,:,C%ys) = i3
        xo([E_fw_m,E_fw_p],C%zs:C%ze-1,:,C%ys) = i8
      endif
      if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh west
        ! no foreign stream dependencies
      endif

      do j=C%ys,C%ye
        do i=C%xs,C%xe      
          do k=C%zs,C%ze-1
            if( atm%l1d(k,i,j) ) then

              xo(:  ,k,i,j) = i0
              xd(:  ,k,i,j) = i1
              xd(0:1,k,i,j) = i3

            endif
          enddo
        enddo
      enddo

      o_nnz=int(xo1d)
      d_nnz=int(xd1d)

      call restoreVecPointer(v_o_nnz,C,xo1d,xo)
      call restoreVecPointer(v_d_nnz,C,xd1d,xd)

      call DMRestoreGlobalVector(C%da,v_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMRestoreGlobalVector(C%da,v_d_nnz,ierr) ;CHKERRQ(ierr)

      if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
  end subroutine 
  subroutine setup_dir8_preallocation(d_nnz,o_nnz,C)
      PetscInt,allocatable :: d_nnz(:)
      PetscInt,allocatable :: o_nnz(:)
      type(t_coord) :: C
      Vec :: v_o_nnz,v_d_nnz
      PetscScalar,Pointer :: xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
      PetscScalar,Pointer :: xo1d(:)=>null(),xd1d(:)=>null()

      PetscInt :: vsize,i,j,k,s

      logical :: lsun_east,lsun_north
      lsun_east  = (sun%xinc.eq.i0)
      lsun_north = (sun%yinc.eq.i0 )

      if(myid.eq.0.and.ldebug) print *,myid,'building direct o_nnz for mat with',C%dof,'dof'
      call DMGetGlobalVector(C%da,v_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMGetGlobalVector(C%da,v_d_nnz,ierr) ;CHKERRQ(ierr)

      call VecGetLocalSize(v_o_nnz,vsize,ierr) ;CHKERRQ(ierr)
      allocate(o_nnz(0:vsize-1))
      allocate(d_nnz(0:vsize-1))

      call getVecPointer(v_o_nnz,C,xo1d,xo)
      call getVecPointer(v_d_nnz,C,xd1d,xd)

      xo = i0
      xd = C%dof+i1

      xd( i0:i3 ,C%zs+1:C%ze  , C%xs:C%xe, C%ys:C%ye ) = C%dof+i1 ! Edir_vertical depends on 3 values Edir_vertical,xaxis,yaxis :: starting with second entries(seen from top)
      xd( i4:i7 ,C%zs  :C%ze-1, C%xs:C%xe, C%ys:C%ye ) = C%dof+i1 ! Edir_xaxis,yaxis depends on 3 values Edir_vertical,xaxis,yaxis :: starting with first entries(seen from top)

      do j=C%ys,C%ye
        if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! real neigh east
          if( lsun_east ) then 
            ! if the sun is in the east, the channels in the last box are influenced by the 2nd channel which is a ghost
            xo(i0:i3, C%zs+1:C%ze, C%xe, j) = xo(i0:i3, C%zs+1:C%ze, C%xe, j)+i2 ! channel 1 from zs+1 to ze
            xd(i0:i3, C%zs+1:C%ze, C%xe, j) = xd(i0:i3, C%zs+1:C%ze, C%xe, j)-i2
            !                                if(myid.eq.0.and.ldebug.and.j.eq.C%ys) print *,myid,'Dir prealloc 0:3: lsun_east :: xo',xo(i0:i3, C%xe, j, C%zs+1), 'xd',xd(i0:i3, C%xe, j, C%zs+1)

            xo(i4:i7, C%zs:C%ze-1, C%xe, j) = xo(i4:i7, C%zs:C%ze-1, C%xe, j)+i2 ! channel 2 and 3 from zs to ze-1
            xd(i4:i7, C%zs:C%ze-1, C%xe, j) = xd(i4:i7, C%zs:C%ze-1, C%xe, j)-i2
            !                                if(myid.eq.0.and.ldebug.and.j.eq.C%ys) print *,myid,'Dir prealloc 4:7: lsun_east :: xo',xo(i4:i7, C%xe, j, C%zs), 'xd',xd(i4:i7, C%xe, j, C%zs)
            if(ldebug) then
              if(any(xd.lt.i0)) print *,myid,'lsun_east :: something wrong happened, we can not have preallocation to be less than 0 in xd!'
              if(any(xo.lt.i0)) print *,myid,'lsun_east :: something wrong happened, we can not have preallocation to be less than 0 in xo!'
            endif
          endif
        endif
      enddo
      do i=C%xs,C%xe
        if( C%neighbors(22).ne.myid .and. C%neighbors(22).ge.i0 ) then ! real neigh north
          if( lsun_north ) then 
            ! if the sun is in the north, the 3rd channel is a ghost
            !                                if(myid.eq.0.and.ldebug.and.i.eq.C%xs) print *,myid,'before Dir prealloc 0:3: lsun_north :: xo',xo(i0:i3, i, C%ye, C%zs+1), 'xd',xd(i0:i3, i, C%ye, C%zs+1)
            xo(i0:i3, C%zs+1:C%ze, i, C%ye) = xo(i0:i3, C%zs+1:C%ze, i, C%ye)+i2 ! channel 1 from zs+1 to ze
            xd(i0:i3, C%zs+1:C%ze, i, C%ye) = xd(i0:i3, C%zs+1:C%ze, i, C%ye)-i2
            !                                if(myid.eq.0.and.ldebug.and.i.eq.C%xs) print *,myid,'Dir prealloc 0:3: lsun_north :: xo',xo(i0:i3, i, C%ye, C%zs+1), 'xd',xd(i0:i3, i, C%ye, C%zs+1)

            xo(i4:i7, C%zs:C%ze-1,  i, C%ye) = xo(i4:i7, C%zs:C%ze-1, i, C%ye)+i2 ! channel 2 and 3 from zs
            xd(i4:i7, C%zs:C%ze-1,  i, C%ye) = xd(i4:i7, C%zs:C%ze-1, i, C%ye)-i2
            !                                if(myid.eq.0.and.ldebug.and.i.eq.C%xs) print *,myid,'Dir prealloc 4:7: lsun_norht :: xo',xo(i4:i7, i, C%ye, C%zs), 'xd',xd(i4:i7, i, C%ye, C%zs)
            if(ldebug) then
              if(any(xd.lt.i0)) print *,myid,'lsun_north :: something wrong happened, we can not have preallocation to be less than 0 in xd!'
              if(any(xo.lt.i0)) print *,myid,'lsun_north :: something wrong happened, we can not have preallocation to be less than 0 in xo!'
            endif
          endif
        endif
      enddo
      do j=C%ys,C%ye
        lsun_east  = (sun%xinc.eq.i0)

        if( C%neighbors(10).ne.myid.and. C%neighbors(10).ge.i0 ) then ! real neigh west
          if( .not. lsun_east ) then 
            ! if the sun is in the west, the 2nd channel is solemnly dependant on ghost values
            xo(i4:i5, C%zs:C%ze-1, C%xs, j) = C%dof
            xd(i4:i5, C%zs:C%ze-1, C%xs, j) = i1
          endif
        endif
      enddo
      do i=C%xs,C%xe
        lsun_north = (sun%yinc.eq.i0 )

        if( C%neighbors(4).ne.myid .and. C%neighbors(4).ge.i0 ) then ! real neigh south
          if( .not. lsun_north ) then 
            ! if the sun is in the south, the 3rd channel is solemnly dependant on ghost values
            xo(i6:i7, C%zs:C%ze-1, i, C%ys) = C%dof
            xd(i6:i7, C%zs:C%ze-1, i, C%ys) = i1
          endif
        endif
      enddo

      do j=C%ys,C%ye
        do i=C%xs,C%xe      
          do k=C%zs,C%ze-1
            if( atm%l1d(k,i,j) ) then
              xo(:  ,k,i,j) = i0
              xd(:  ,k,i,j) = i1
              xd(0:3,k,i,j) = i5
            endif
          enddo
        enddo
      enddo

      o_nnz=int(xo1d)
      d_nnz=int(xd1d)

      call restoreVecPointer(v_o_nnz,C,xo1d,xo)
      call restoreVecPointer(v_d_nnz,C,xd1d,xd)

      call DMRestoreGlobalVector(C%da,v_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMRestoreGlobalVector(C%da,v_d_nnz,ierr) ;CHKERRQ(ierr)

      if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
  end subroutine 
  subroutine setup_dir_preallocation(d_nnz,o_nnz,C)
      PetscInt,allocatable :: d_nnz(:)
      PetscInt,allocatable :: o_nnz(:)
      type(t_coord) :: C
      Vec :: gv_o_nnz,gv_d_nnz
      Vec :: lv_nnz
      PetscScalar,Pointer :: xl(:,:,:,:)=>null(),xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
      PetscScalar,Pointer :: xl1d(:)=>null(),xo1d(:)=>null(),xd1d(:)=>null()

      PetscInt :: vsize,i,j,k,src,dst

      logical :: lsun_east,lsun_north

      call DMGetGlobalVector(C%da, gv_d_nnz,ierr) ;CHKERRQ(ierr)
      call DMGetGlobalVector(C%da, gv_o_nnz,ierr) ;CHKERRQ(ierr)
      call VecSet(gv_d_nnz, zero, ierr) ;CHKERRQ(ierr)
      call VecSet(gv_o_nnz, zero, ierr) ;CHKERRQ(ierr)

      call DMGetLocalVector(C%da, lv_nnz,ierr) ;CHKERRQ(ierr)

      call getVecPointer(lv_nnz,C,xl1d,xl)

      xl=zero

      do j=C%ys,C%ye
        do i=C%xs,C%xe        
          do k=C%zs,C%ze-1
            if( atm%l1d(k,i,j) ) then
              do dst=i0,i3
                call inc( xl(dst, k+1, i,j), one )
              enddo
            else

               dst = i0 ;call inc( xl(dst , k+1 , i          , j          ) , one*C%dof )
               dst = i1 ;call inc( xl(dst , k+1 , i          , j          ) , one*C%dof )
               dst = i2 ;call inc( xl(dst , k+1 , i          , j          ) , one*C%dof )
               dst = i3 ;call inc( xl(dst , k+1 , i          , j          ) , one*C%dof )
               dst = i4 ;call inc( xl(dst , k   , i+sun%xinc , j          ) , one*C%dof )
               dst = i5 ;call inc( xl(dst , k   , i+sun%xinc , j          ) , one*C%dof )
               dst = i6 ;call inc( xl(dst , k   , i          , j+sun%yinc ) , one*C%dof )
               dst = i7 ;call inc( xl(dst , k   , i          , j+sun%yinc ) , one*C%dof )
                     
            endif

          enddo 
        enddo 
      enddo
      call restoreVecPointer(lv_nnz,C,xl1d,xl)

      ! Now we have a local vector with ghosts that counts all entries
      ! We still have to seperate the diagonal and the off-diagonal part:
      ! For that, transmit the ghost parts with ADD_VALUES to _o_nnz
      call DMLocalToGlobalBegin(C%da,lv_nnz,ADD_VALUES,gv_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMLocalToGlobalEnd  (C%da,lv_nnz,ADD_VALUES,gv_o_nnz,ierr) ;CHKERRQ(ierr)

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

      call DMRestoreGlobalVector(C%da, gv_d_nnz,ierr) ;CHKERRQ(ierr)
      call DMRestoreGlobalVector(C%da, gv_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMRestoreLocalVector (C%da, lv_nnz  ,ierr) ;CHKERRQ(ierr)
  end subroutine 
  subroutine setup_diff_preallocation(d_nnz,o_nnz,C)
      PetscInt,allocatable :: d_nnz(:)
      PetscInt,allocatable :: o_nnz(:)
      type(t_coord) :: C
      Vec :: gv_o_nnz,gv_d_nnz
      Vec :: lv_nnz
      PetscScalar,Pointer :: xl(:,:,:,:)=>null(),xo(:,:,:,:)=>null(),xd(:,:,:,:)=>null()
      PetscScalar,Pointer :: xl1d(:)=>null(),xo1d(:)=>null(),xd1d(:)=>null()

      PetscInt :: vsize,i,j,k,src,dst

      logical :: lsun_east,lsun_north

      call DMGetGlobalVector(C%da, gv_d_nnz,ierr) ;CHKERRQ(ierr)
      call DMGetGlobalVector(C%da, gv_o_nnz,ierr) ;CHKERRQ(ierr)
      call VecSet(gv_d_nnz, zero, ierr) ;CHKERRQ(ierr)
      call VecSet(gv_o_nnz, zero, ierr) ;CHKERRQ(ierr)

      call DMGetLocalVector(C%da, lv_nnz,ierr) ;CHKERRQ(ierr)

      call getVecPointer(lv_nnz,C,xl1d,xl)

      xl=zero

      do j=C%ys,C%ye
        do i=C%xs,C%xe        
          do k=C%zs,C%ze-1
            if( atm%l1d(k,i,j) ) then
              do dst=i0,i1
                call inc( xl(dst, k+1, i,j), one )
                call inc( xl(dst, k  , i,j), one )
              enddo
            else

               dst = E_up  ;call inc( xl(dst , k   , i   , j   ) , one*C%dof )
               dst = E_dn  ;call inc( xl(dst , k+1 , i   , j   ) , one*C%dof )
               dst = E_le_m;call inc( xl(dst , k   , i   , j   ) , one*C%dof )
               dst = E_ri_m;call inc( xl(dst , k   , i+1 , j   ) , one*C%dof-1 )
               dst = E_le_p;call inc( xl(dst , k   , i   , j   ) , one*C%dof )
               dst = E_ri_p;call inc( xl(dst , k   , i+1 , j   ) , one*C%dof-1 )
               dst = E_ba_m;call inc( xl(dst , k   , i   , j   ) , one*C%dof )
               dst = E_fw_m;call inc( xl(dst , k   , i   , j+1 ) , one*C%dof-1 )
               dst = E_ba_p;call inc( xl(dst , k   , i   , j   ) , one*C%dof )
               dst = E_fw_p;call inc( xl(dst , k   , i   , j+1 ) , one*C%dof-1 )
                     
              print *,myid,'diff',k,i,j,'xl',int(xl(:,k,i,j))
            endif

          enddo 
        enddo 
      enddo
      call restoreVecPointer(lv_nnz,C,xl1d,xl)

      ! Now we have a local vector with ghosts that counts all entries
      ! We still have to seperate the diagonal and the off-diagonal part:
      ! For that, transmit the ghost parts with ADD_VALUES to _o_nnz
      call DMLocalToGlobalBegin(C%da,lv_nnz,ADD_VALUES,gv_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMLocalToGlobalEnd  (C%da,lv_nnz,ADD_VALUES,gv_o_nnz,ierr) ;CHKERRQ(ierr)

      call getVecPointer(lv_nnz,C,xl1d,xl)
      call getVecPointer(gv_o_nnz,C,xo1d,xo)
      call getVecPointer(gv_d_nnz,C,xd1d,xd)

      do j=C%ys,C%ye
        do i=C%xs,C%xe        
          do k=C%zs,C%ze
            xo(:,k,i,j) = xo(:,k,i,j)-xl(:,k,i,j)
            xd(:,k,i,j) = + C%dof+i1  -xo(:,k,i,j) 
            print *,myid,'diff',k,i,j,'on',int(xd(:,k,i,j)),'off',int(xo(:,k,i,j))
          enddo
        enddo
      enddo

      if( C%neighbors(10) .eq. C%neighbors(13) ) then ! neighbour in x-dir is actually me again....
        xd(:,:,C%xs,:) = C%dof+1
        xd(:,:,C%xe,:) = C%dof+1
      endif
      if( C%neighbors(4) .eq. C%neighbors(13) ) then ! neighbour in y-dir is actually me again....
        xd(:,:,:,C%ys) = C%dof+1
        xd(:,:,:,C%ye) = C%dof+1
      endif

      allocate(o_nnz(0:C%dof*C%zm*C%xm*C%ym-1))
      allocate(d_nnz(0:C%dof*C%zm*C%xm*C%ym-1))
      o_nnz=int(xo1d)
      d_nnz=int(xd1d)

      call restoreVecPointer(lv_nnz,C,xl1d,xl)
      call restoreVecPointer(gv_o_nnz,C,xo1d,xo)
      call restoreVecPointer(gv_d_nnz,C,xd1d,xd)

      call DMRestoreGlobalVector(C%da, gv_d_nnz,ierr) ;CHKERRQ(ierr)
      call DMRestoreGlobalVector(C%da, gv_o_nnz,ierr) ;CHKERRQ(ierr)
      call DMRestoreLocalVector (C%da, lv_nnz,  ierr) ;CHKERRQ(ierr)
  end subroutine 

  subroutine get_coeff(op,dz,dir,coeff,lone_dimensional,angles)
      type(t_optprop),intent(in) :: op
      real(ireals),intent(in) :: dz
      logical,intent(in) :: dir
      real(ireals),intent(out) :: coeff(:)

      logical,intent(in) :: lone_dimensional
      real(ireals),intent(in),optional :: angles(2)

      call PetscLogStagePush(logstage(7),ierr) ;CHKERRQ(ierr)

      if(lone_dimensional) then
        call OPP_1_2%get_coeff(dz,op%kabs,op%ksca,op%g,dir,coeff,angles)
      else
        call OPP_8_10%get_coeff(dz,op%kabs,op%ksca,op%g,dir,coeff,angles)
      endif

      if(ldebug) then
        if( any(isnan(coeff)) .or. any(coeff.lt.zero) .or. any(coeff.gt.one) ) print *,'Wrong coeff',coeff,'op',op
      endif
      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
  end subroutine
  function sym_rot_phi(phi0)
      real(ireals) :: sym_rot_phi
      real(ireals),intent(in) :: phi0
      ! ''swap'' phi axis down to the range of [0,180] 
      sym_rot_phi = acos(cos(phi0*pi/180))
      !  print *,'1st phi0 swap',phi0,' :: ',sym_rot_phi,'=',phi0*pi/180,cos(phi0*pi/180),acos(cos(phi0*pi/180))
      ! and then mirror it onto range [0,90]
      sym_rot_phi = int( asin(sin(sym_rot_phi)) /pi * 180 )
      !  print *,'2nd phi0 swap',phi0,' :: ',sym_rot_phi,'=',sin(sym_rot_phi),asin(sin(sym_rot_phi)),asin(sin(sym_rot_phi)) /pi * 180,int(asin(sin(sym_rot_phi)) /pi * 180)
  end function
  subroutine setup_suninfo(phi0,theta0,sun)
      real(ireals),intent(in) :: phi0,theta0
      type(t_suninfo),intent(out) :: sun

      sun%phi   = phi0
      sun%theta = theta0
      sun%costheta = max( cos(deg2rad(theta0)), zero)
      sun%sintheta = max( sin(deg2rad(theta0)), zero)

      ! use symmetry for direct beam: always use azimuth [0,90] an just reverse the order where we insert the coeffs
      sun%symmetry_phi = sym_rot_phi(phi0)
      sun%xinc=i0 ; sun%yinc=i0
      if(phi0.gt.180) sun%xinc=i1
      if(phi0.gt.90.and.phi0.lt.270) sun%yinc=i1
      if(ldebug.and.myid.eq.0) print *,'setup_dir_inc done', sun
  end subroutine

  subroutine set_dir_coeff(A,C)
      Mat :: A
      type(t_coord) :: C

      PetscInt :: i,j,k

      call PetscLogStagePush(logstage(2),ierr) ;CHKERRQ(ierr)
      if(myid.eq.0.and.ldebug) print *,myid,'setup_direct_matrix ...'

!      call MatZeroEntries(A, ierr) ;CHKERRQ(ierr) !TODO necessary?
!      call mat_set_diagonal(A,C)

      do j=C%ys,C%ye
        do i=C%xs,C%xe        
          do k=C%zs,C%ze-1

            if( atm%l1d(k,i,j) ) then
              call set_eddington_coeff(A,k, i,j)
            else
              call set_tenstream_coeff(C, A,k,i,j)
            endif

          enddo 
        enddo 
      enddo

      if(myid.eq.0.and.ldebug) print *,myid,'setup_direct_matrix done'

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)  

      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

    contains 
      subroutine set_tenstream_coeff(C,A,k,i,j)
          type(t_coord),intent(in) :: C
          Mat,intent(inout) :: A
          integer(iintegers),intent(in) :: i,j,k

          MatStencil :: row(4,C%dof)  ,col(4,C%dof)
          PetscScalar :: v(C%dof**2),coeffs(C%dof**2),norm

          PetscInt,parameter :: entries(8)=[0,8,16,24,32,40,48,56]

          integer(iintegers) :: dst,src

          dst = 1 ; row(MatStencil_j,dst) = i            ; row(MatStencil_k,dst) = j            ; row(MatStencil_i,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 2 ; row(MatStencil_j,dst) = i            ; row(MatStencil_k,dst) = j            ; row(MatStencil_i,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 3 ; row(MatStencil_j,dst) = i            ; row(MatStencil_k,dst) = j            ; row(MatStencil_i,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 4 ; row(MatStencil_j,dst) = i            ; row(MatStencil_k,dst) = j            ; row(MatStencil_i,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 5 ; row(MatStencil_j,dst) = i+sun%xinc   ; row(MatStencil_k,dst) = j            ; row(MatStencil_i,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
          dst = 6 ; row(MatStencil_j,dst) = i+sun%xinc   ; row(MatStencil_k,dst) = j            ; row(MatStencil_i,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
          dst = 7 ; row(MatStencil_j,dst) = i            ; row(MatStencil_k,dst) = j+sun%yinc   ; row(MatStencil_i,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
          dst = 8 ; row(MatStencil_j,dst) = i            ; row(MatStencil_k,dst) = j+sun%yinc   ; row(MatStencil_i,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid

          src = 1 ; col(MatStencil_j,src) = i            ; col(MatStencil_k,src) = j            ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 2 ; col(MatStencil_j,src) = i            ; col(MatStencil_k,src) = j            ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 3 ; col(MatStencil_j,src) = i            ; col(MatStencil_k,src) = j            ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 4 ; col(MatStencil_j,src) = i            ; col(MatStencil_k,src) = j            ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 5 ; col(MatStencil_j,src) = i+1-sun%xinc ; col(MatStencil_k,src) = j            ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
          src = 6 ; col(MatStencil_j,src) = i+1-sun%xinc ; col(MatStencil_k,src) = j            ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
          src = 7 ; col(MatStencil_j,src) = i            ; col(MatStencil_k,src) = j+1-sun%yinc ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:
          src = 8 ; col(MatStencil_j,src) = i            ; col(MatStencil_k,src) = j+1-sun%yinc ; col(MatStencil_i,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:

          call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.True., coeffs, atm%l1d(k,i,j), [sun%symmetry_phi, sun%theta])

          if(ldebug) then
            do src=1,C%dof
              norm = sum( coeffs((src-1)*C%dof+1:src*C%dof) )
              if( ldebug .and. real(norm).gt.real(one) ) then
                print *,'sum(src==',src,') gt one',norm
                stop 'omg.. shouldnt be happening'
                ierr=-5
                return
              endif
            enddo
          endif

          !TODO we should really save the coefficients in another order... does not make sense that we have to revert them from
          !source odered to destination order every time....
          do src=1,C%dof
            v(entries+src) = coeffs( i1+(src-i1)*C%dof : src*C%dof )
          enddo
          call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)

      end subroutine

      subroutine set_eddington_coeff(A,k,i,j)
          Mat,intent(inout) :: A
          integer(iintegers),intent(in) :: i,j,k

          MatStencil :: row(4,1)  ,col(4,1)
          PetscScalar :: v(1)
          integer(iintegers) :: src

          if(luse_eddington) then
            v = atm%a33(k,i,j)
          else
            call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.True., v, atm%l1d(k,i,j), [sun%symmetry_phi, sun%theta] )
          endif

          col(MatStencil_j,i1) = i      ; col(MatStencil_k,i1) = j       ; col(MatStencil_i,i1) = k    
          row(MatStencil_j,i1) = i      ; row(MatStencil_k,i1) = j       ; row(MatStencil_i,i1) = k+1  

          do src=1,4
            col(MatStencil_c,i1) = src-i1 ! Source may be the upper/lower lid:
            row(MatStencil_c,i1) = src-i1 ! Define transmission towards the lower/upper lid
            call MatSetValuesStencil(A,i1, row,i1, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)
          enddo
      end subroutine

  end subroutine

  subroutine setup_incSolar(incSolar,edirTOA)
      Vec :: incSolar
      real(ireals),intent(in) :: edirTOA

      PetscScalar,pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

      PetscReal :: Az
      Az = atm%dx*atm%dy

      call VecSet(incSolar,zero,ierr) ;CHKERRQ(ierr)

      call getVecPointer(incSolar,C_dir,x1d,x4d,.False.)

      x4d(i0:i3,C_dir%zs,:,:) = edirTOA* Az * .25_ireals * sun%costheta

      call restoreVecPointer(incSolar,C_dir,x1d,x4d)

      if(myid.eq.0 .and. ldebug) print *,myid,'Setup of IncSolar done',edirTOA* sun%costheta

  end subroutine

  subroutine set_diff_coeff(A,C)
      Mat :: A
      type(t_coord) :: C

      PetscInt :: i,j,k

      call PetscLogStagePush(logstage(4),ierr) ;CHKERRQ(ierr)

      if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients for diffuse Light'

!      call MatZeroEntries(A, ierr) ;CHKERRQ(ierr) !TODO necessary?
!      call mat_set_diagonal(A,C)

      do j=C%ys,C%ye
        do i=C%xs,C%xe
          do k=C%zs,C%ze-1

            if( atm%l1d(k,i,j) ) then
              call set_eddington_coeff(A, k,i,j)
            else
              call set_tenstream_coeff(C, A, k,i,j, ierr); CHKERRQ(ierr)
            endif

          enddo 
        enddo 
      enddo

      call set_albedo_coeff(C, A )

      if(myid.eq.0.and.ldebug) print *,myid,'setup_diffuse_matrix done'

      if(myid.eq.0.and.ldebug) print *,myid,'Final diffuse Matrix Assembly:'
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)

      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

    contains
      subroutine set_tenstream_coeff(C,A,k,i,j,ierr)
          type(t_coord),intent(in) :: C
          Mat,intent(inout) :: A
          integer(iintegers),intent(in) :: k,i,j
          PetscErrorCode,intent(out) :: ierr

          MatStencil :: row(4,0:C%dof-1)  ,col(4,0:C%dof-1)
          PetscReal :: v(C%dof**2),coeffs(C%dof**2),norm
          PetscInt,parameter :: entries(10)=[0,10,20,30,40,50,60,70,80,90]

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

          call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.False., coeffs, atm%l1d(k,i,j))

          if(ldebug) then
            do dst=1,C%dof
              norm = sum( coeffs((dst-1)*C%dof+1:dst*C%dof) )
              if( ldebug ) then
                if( real(norm).gt.real(one) ) then
                  print *,'diffuse sum(dst==',dst,') gt one',norm
                  !                coeffs((dst-1)*C%dof+1:dst*C%dof)  = coeffs((dst-1)*C%dof+1:dst*C%dof) / (norm+1e-8_ireals)
                  stop 'omg.. shouldnt be happening'
                  ierr = -5
                  return
                endif
              endif
            enddo
          endif

          do src=1,C%dof
            v(entries+src) = coeffs( i1+(src-i1)*C%dof : src*C%dof )
          enddo
          call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)

      end subroutine

      subroutine set_eddington_coeff(A,k,i,j)
          Mat,intent(inout) :: A
          integer(iintegers),intent(in) :: k,i,j

          MatStencil :: row(4,2)  ,col(4,2)
          PetscReal :: v(4),twostr_coeff(4) ! v ==> a12,a11,a11,a12
          integer(iintegers) :: src,dst

          if(luse_eddington ) then
            v = [ atm%a12(k,i,j), atm%a11(k,i,j), atm%a11(k,i,j), atm%a12(k,i,j)]
          else
            call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.False., twostr_coeff, atm%l1d(k,i,j)) !twostr_coeff ==> a12,a11,a12,a11
            v = [ twostr_coeff(1), twostr_coeff(2) , twostr_coeff(2) , twostr_coeff(1) ]
          endif

          src = 1; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k     ; col(MatStencil_c,src) = E_dn
          src = 2; col(MatStencil_j,src) = i    ; col(MatStencil_k,src) = j     ; col(MatStencil_i,src) = k+1   ; col(MatStencil_c,src) = E_up 

          dst = 1; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k     ; row(MatStencil_c,dst) = E_up  
          dst = 2; row(MatStencil_j,dst) = i    ; row(MatStencil_k,dst) = j     ; row(MatStencil_i,dst) = k+1   ; row(MatStencil_c,dst) = E_dn  

          call MatSetValuesStencil(A,i2, row, i2, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)


      end subroutine
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

              call MatSetValuesStencil(A,i1, row, i1, col , [-atm%albedo] ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)

            enddo
          enddo

      end subroutine
  end subroutine

  subroutine setup_b(solution,b)
      type(t_state_container) :: solution
      Vec :: local_b,b

      PetscScalar,pointer,dimension(:,:,:,:) :: xsrc=>null()
      PetscScalar,pointer,dimension(:) :: xsrc1d=>null()
      PetscInt :: k,i,j,src

      call PetscLogStagePush(logstage(6),ierr) ;CHKERRQ(ierr)

      if(myid.eq.0.and.ldebug) print *,'src Vector Assembly...'

      call DMGetLocalVector(C_diff%da,local_b,ierr) ;CHKERRQ(ierr)
      call VecSet(local_b,zero,ierr) ;CHKERRQ(ierr)

      call getVecPointer(local_b,C_diff,xsrc1d,xsrc)

      if(solution%lsolar_rad) &
          call set_solar_source(solution%edir)

      if(allocated(atm%planck) ) &
          call set_thermal_source()

      if(myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients ...done'

      call restoreVecPointer(local_b, C_diff, xsrc1d , xsrc )

      call VecSet(b,zero,ierr) ;CHKERRQ(ierr) ! reset global Vec

      call DMLocalToGlobalBegin(C_diff%da, local_b, ADD_VALUES, b,ierr) ;CHKERRQ(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
      call DMLocalToGlobalEnd  (C_diff%da, local_b, ADD_VALUES, b,ierr) ;CHKERRQ(ierr)

      call DMRestoreLocalVector(C_diff%da,local_b,ierr) ;CHKERRQ(ierr)

      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
      if(myid.eq.0.and.ldebug) print *,'src Vector Assembly done'

    contains
      subroutine set_thermal_source()
          PetscReal :: Ax,Ay,Az,c1,c2,c3,b0,b1,dtau
          real(ireals) :: diff2diff1d(4)
          PetscReal :: diff2diff(C_diff%dof**2)

          if(myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector ... setting thermal source terms'
          Az = atm%dx*atm%dy

          do j=C_diff%ys,C_diff%ye         
            do i=C_diff%xs,C_diff%xe    
              do k=C_diff%zs,C_diff%ze-1 

                if( atm%l1d(k,i,j) ) then

                  if(luse_eddington ) then
                    !see libradtran rodents
                    dtau = atm%dz(k,i,j) * ( atm%op(k,i,j)%kabs + atm%op(k,i,j)%ksca)
                    if( dtau.gt.0.01_ireals ) then
                      b0 = atm%planck(k,i,j)
                      b1 = ( atm%planck(k,i,j)-atm%planck(k+1,i,j) ) / dtau
                    else
                      b0 = .5_ireals*(atm%planck(k,i,j)+atm%planck(k+1,i,j))
                      b1 = zero
                    endif
                    c1 = atm%g1(k,i,j) * (b0 + b1*dtau)
                    c2 = atm%g2(k,i,j) * b1
                    c3 = atm%g1(k,i,j) * b0

                    xsrc(E_up   ,k  ,i,j) = xsrc(E_up   ,k  ,i,j) + ( - atm%a11(k,i,j)*(c1+c2) - atm%a12(k,i,j)*(c3-c2) + c2 + c3 )*Az*pi
                    xsrc(E_dn   ,k+1,i,j) = xsrc(E_dn   ,k+1,i,j) + ( - atm%a12(k,i,j)*(c1+c2) - atm%a11(k,i,j)*(c3-c2) + c1 - c2 )*Az*pi

                  else

                    call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.False., diff2diff1d, atm%l1d(k,i,j))

                    b0 = .5_ireals*(atm%planck(k,i,j)+atm%planck(k+1,i,j)) *pi
                    xsrc(E_up   ,k  ,i,j) = xsrc(E_up   ,k  ,i,j) +  b0  *(one-diff2diff1d(1)-diff2diff1d(2) ) *Az
                    xsrc(E_dn   ,k+1,i,j) = xsrc(E_dn   ,k+1,i,j) +  b0  *(one-diff2diff1d(1)-diff2diff1d(2) ) *Az

                  endif

                else ! Tenstream source terms
                  Ax = atm%dy*atm%dz(k,i,j)
                  Ay = atm%dx*atm%dz(k,i,j)

                  call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.False., diff2diff, atm%l1d(k,i,j) )
                  b0 = .5_ireals*(atm%planck(k,i,j)+atm%planck(k+i1,i,j)) *pi
                  xsrc(E_up   , k   , i   , j   ) = xsrc(E_up   , k   , i   , j   ) +  b0  *(one-sum( diff2diff( E_up  *C_diff%dof+i1 : E_up  *C_diff%dof+C_diff%dof )  )  ) *Az
                  xsrc(E_dn   , k+1 , i   , j   ) = xsrc(E_dn   , k+1 , i   , j   ) +  b0  *(one-sum( diff2diff( E_dn  *C_diff%dof+i1 : E_dn  *C_diff%dof+C_diff%dof )  )  ) *Az
                  xsrc(E_le_m , k   , i   , j   ) = xsrc(E_le_m , k   , i   , j   ) +  b0  *(one-sum( diff2diff( E_le_m*C_diff%dof+i1 : E_le_m*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
                  xsrc(E_le_p , k   , i   , j   ) = xsrc(E_le_p , k   , i   , j   ) +  b0  *(one-sum( diff2diff( E_le_p*C_diff%dof+i1 : E_le_p*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
                  xsrc(E_ri_m , k   , i+1 , j   ) = xsrc(E_ri_m , k   , i+1 , j   ) +  b0  *(one-sum( diff2diff( E_ri_m*C_diff%dof+i1 : E_ri_m*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
                  xsrc(E_ri_p , k   , i+1 , j   ) = xsrc(E_ri_p , k   , i+1 , j   ) +  b0  *(one-sum( diff2diff( E_ri_p*C_diff%dof+i1 : E_ri_p*C_diff%dof+C_diff%dof )  )  ) *Ax*.5_ireals
                  xsrc(E_ba_m , k   , i   , j   ) = xsrc(E_ba_m , k   , i   , j   ) +  b0  *(one-sum( diff2diff( E_ba_m*C_diff%dof+i1 : E_ba_m*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
                  xsrc(E_ba_p , k   , i   , j   ) = xsrc(E_ba_p , k   , i   , j   ) +  b0  *(one-sum( diff2diff( E_ba_p*C_diff%dof+i1 : E_ba_p*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
                  xsrc(E_fw_m , k   , i   , j+1 ) = xsrc(E_fw_m , k   , i   , j+1 ) +  b0  *(one-sum( diff2diff( E_fw_m*C_diff%dof+i1 : E_fw_m*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
                  xsrc(E_fw_p , k   , i   , j+1 ) = xsrc(E_fw_p , k   , i   , j+1 ) +  b0  *(one-sum( diff2diff( E_fw_p*C_diff%dof+i1 : E_fw_p*C_diff%dof+C_diff%dof )  )  ) *Ay*.5_ireals
                endif ! 1D or Tenstream?

              enddo
            enddo

          enddo ! k

          ! Thermal emission at surface
          k = C_diff%ze
          do j=C_diff%ys,C_diff%ye         
            do i=C_diff%xs,C_diff%xe    
              xsrc(E_up   ,k,i,j) = xsrc(E_up   ,k,i,j) + atm%planck(k,i,j)*Az *(one-atm%albedo)*pi
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
          lsun_east  = (sun%xinc.eq.i0)
          lsun_north = (sun%yinc.eq.i0 )


          ! Copy ghosted values for direct vec
          call DMGetLocalVector(C_dir%da ,ledir ,ierr)                      ; CHKERRQ(ierr)
          call VecSet(ledir ,zero,ierr)                                     ; CHKERRQ(ierr)
          call DMGlobalToLocalBegin(C_dir%da ,edir ,ADD_VALUES,ledir ,ierr) ; CHKERRQ(ierr)
          call DMGlobalToLocalEnd  (C_dir%da ,edir ,ADD_VALUES,ledir ,ierr) ; CHKERRQ(ierr)

          call getVecPointer(ledir,C_dir,xedir1d,xedir)

          if(myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector .. setting solar source',sum(xedir(0:3,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye))/size(xedir(0:3,C_dir%zs:C_dir%ze,C_dir%xs:C_dir%xe,C_dir%ys:C_dir%ye))
          do j=C_diff%ys,C_diff%ye         
            do i=C_diff%xs,C_diff%xe    
              do k=C_diff%zs,C_diff%ze-1 

                if( any (xedir(:,k,i,j) .gt. epsilon(one)) ) then
                  if( atm%l1d(k,i,j) ) then
                    dir2diff = zero
                    if(luse_eddington ) then
                      ! Only transport the 4 tiles from dir0 to the Eup and Edn
                      do src=1,4
                        dir2diff(E_up  +i1+(src-1)*C_diff%dof) = atm%a13(k,i,j)
                        dir2diff(E_dn  +i1+(src-1)*C_diff%dof) = atm%a23(k,i,j)
                      enddo

                    else
                      call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.False., twostr_coeff, atm%l1d(k,i,j), [sun%symmetry_phi, sun%theta])
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

                    call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.False., dir2diff,  atm%l1d(k,i,j), [sun%symmetry_phi, sun%theta] )

                    do src=1,C_dir%dof
                      select case(src)
                      case (1:4)
                        solrad = xedir( : , k , i , j )
                      case (5:6)
                        solrad = xedir( : , k , i+i1-sun%xinc , j )
                      case (7:8)
                        solrad = xedir( : , k , i , j+i1-sun%yinc )
                      case default
                        stop 'invalid dof for solar source term'
                      end select
                      xsrc(E_up   , k   , i   , j   ) = xsrc(E_up   , k   , i   , j   ) +  solrad(src) *dir2diff(E_up  +i1 + (src-1)*C_diff%dof )
                      xsrc(E_dn   , k+1 , i   , j   ) = xsrc(E_dn   , k+1 , i   , j   ) +  solrad(src) *dir2diff(E_dn  +i1 + (src-1)*C_diff%dof )

                      if(lsun_east) then ! if sun shines from right to left, switch ''reflection'' and ''transmission'' coefficients.
                        xsrc(E_le_m , k   , i   , j   ) = xsrc(E_le_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_ri_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_le_p , k   , i   , j   ) = xsrc(E_le_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_ri_p+i1 + (src-1)*C_diff%dof )
                        xsrc(E_ri_m , k   , i+1 , j   ) = xsrc(E_ri_m , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_le_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_ri_p , k   , i+1 , j   ) = xsrc(E_ri_p , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_le_p+i1 + (src-1)*C_diff%dof )
                      else
                        xsrc(E_le_m , k   , i   , j   ) = xsrc(E_le_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_le_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_le_p , k   , i   , j   ) = xsrc(E_le_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_le_p+i1 + (src-1)*C_diff%dof )
                        xsrc(E_ri_m , k   , i+1 , j   ) = xsrc(E_ri_m , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_ri_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_ri_p , k   , i+1 , j   ) = xsrc(E_ri_p , k   , i+1 , j   ) +  solrad(src) *dir2diff(E_ri_p+i1 + (src-1)*C_diff%dof )
                      endif
                      if(lsun_north) then ! likewise if sun shines from forward to backwards, 
                        xsrc(E_ba_m , k   , i   , j   ) = xsrc(E_ba_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_fw_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_ba_p , k   , i   , j   ) = xsrc(E_ba_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_fw_p+i1 + (src-1)*C_diff%dof )
                        xsrc(E_fw_m , k   , i   , j+1 ) = xsrc(E_fw_m , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_ba_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_fw_p , k   , i   , j+1 ) = xsrc(E_fw_p , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_ba_p+i1 + (src-1)*C_diff%dof )
                      else
                        xsrc(E_ba_m , k   , i   , j   ) = xsrc(E_ba_m , k   , i   , j   ) +  solrad(src) *dir2diff(E_ba_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_ba_p , k   , i   , j   ) = xsrc(E_ba_p , k   , i   , j   ) +  solrad(src) *dir2diff(E_ba_p+i1 + (src-1)*C_diff%dof )
                        xsrc(E_fw_m , k   , i   , j+1 ) = xsrc(E_fw_m , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_fw_m+i1 + (src-1)*C_diff%dof )
                        xsrc(E_fw_p , k   , i   , j+1 ) = xsrc(E_fw_p , k   , i   , j+1 ) +  solrad(src) *dir2diff(E_fw_p+i1 + (src-1)*C_diff%dof )
                      endif

                      if(ldebug) then
                        if(sum(dir2diff(i1+(src-1 ) *C_diff%dof : C_dir%dof+(src-1 ) *C_diff%dof ) ) .gt. one .or. &
                           sum(dir2diff(i1+(src-1 ) *C_diff%dof : C_dir%dof+(src-1 ) *C_diff%dof ) ) .lt. zero   ) &
                            print *,'DEBUG Found dir2diff gt one:',src,'::',sum(dir2diff(i1+(src-1 ) *C_diff%dof : C_dir%dof+(src-1 ) *C_diff%dof ) ),'::',dir2diff(i1+(src-1 ) *C_diff%dof : C_dir%dof+(src-1 ) *C_diff%dof ) 
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
              xsrc(E_up   ,k,i,j) = sum(xedir(i0:i3,k,i,j))*atm%albedo
            enddo
          enddo

          call restoreVecPointer(ledir ,C_dir ,xedir1d ,xedir )
          call DMRestoreLocalVector(C_dir%da ,ledir ,ierr)       ; CHKERRQ(ierr)
      end subroutine
  end subroutine

  subroutine calc_flx_div(solution)
      type(t_state_container) :: solution
      PetscReal,pointer,dimension(:,:,:,:) :: xediff=>null(),xedir=>null(),xabso=>null()
      PetscReal,pointer,dimension(:) :: xediff1d=>null(),xedir1d=>null(),xabso1d=>null()
      PetscInt :: i,j,k
      Vec :: ledir,lediff ! local copies of vectors, including ghosts
      PetscReal :: div(3),div2(13)
      PetscReal :: Volume,Ax,Ay,Az

      if(.not. solution%lintegrated_dir) stop 'tried calculating absorption but dir  vector was in [W/m**2], not in [W], scale first!'
      if(.not. solution%lintegrated_diff)stop 'tried calculating absorption but diff vector was in [W/m**2], not in [W], scale first!'

      if( lcalc_nca ) then ! if we should calculate NCA (Klinger), we can just return afterwards
        stop 'TODO: this is not supported atm -- i guess nca does not expect [W] but [W/m**2]'
        call nca_wrapper(solution%ediff,solution%abso)
        return
      endif

      if(myid.eq.0.and.ldebug) print *,'Calculating flux divergence'
      call VecSet(solution%abso,zero,ierr) ;CHKERRQ(ierr)

      ! TODO: if there are no 3D layers, we should skip the ghost value copying....

      if(solution%lsolar_rad) then
        ! Copy ghosted values for direct vec
        call DMGetLocalVector(C_dir%da ,ledir ,ierr)                      ; CHKERRQ(ierr)
        call VecSet(ledir ,zero,ierr)                                     ; CHKERRQ(ierr)
        call DMGlobalToLocalBegin(C_dir%da ,solution%edir ,ADD_VALUES,ledir ,ierr) ; CHKERRQ(ierr)
        call DMGlobalToLocalEnd  (C_dir%da ,solution%edir ,ADD_VALUES,ledir ,ierr) ; CHKERRQ(ierr)
        call getVecPointer(ledir, C_dir ,xedir1d ,xedir )
      endif

      ! Copy ghosted values for diffuse vec
      call DMGetLocalVector(C_diff%da,lediff,ierr)                      ; CHKERRQ(ierr)
      call VecSet(lediff,zero,ierr)                                     ; CHKERRQ(ierr)
      call DMGlobalToLocalBegin(C_diff%da,solution%ediff,ADD_VALUES,lediff,ierr) ; CHKERRQ(ierr)
      call DMGlobalToLocalEnd  (C_diff%da,solution%ediff,ADD_VALUES,lediff,ierr) ; CHKERRQ(ierr)
      call getVecPointer(lediff, C_diff, xediff1d, xediff)

      call getVecPointer(solution%abso, C_one, xabso1d, xabso)

      ! calculate absorption by flux divergence
      Az = atm%dx * atm%dy

      do j=C_one%ys,C_one%ye         
        do i=C_one%xs,C_one%xe      
          do k=C_one%zs,C_one%ze

            Volume = Az     * atm%dz(k,i,j)

            if(atm%l1d(k,i,j)) then ! one dimensional i.e. twostream
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
                div2( 1) = sum( xedir(i0:i3 , k, i             , j             )  - xedir(i0:i3 , k+i1 , i          , j           ) ) 
                div2( 2) = sum( xedir(i4:i5 , k, i+i1-sun%xinc , j             )  - xedir(i4:i5 , k    , i+sun%xinc , j           ) ) 
                div2( 3) = sum( xedir(i6:i7 , k, i             , j+i1-sun%yinc )  - xedir(i6:i7 , k    , i          , j+sun%yinc  ) ) 
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
        call DMRestoreLocalVector(C_dir%da ,ledir ,ierr) ; CHKERRQ(ierr)
      endif

      call restoreVecPointer(lediff         ,C_diff,xediff1d,xediff)
      call DMRestoreLocalVector(C_diff%da,lediff,ierr) ; CHKERRQ(ierr)

      call restoreVecPointer(solution%abso  ,C_one ,xabso1d ,xabso )
  end subroutine

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

      call KSPSolve(ksp,b,x,ierr) ;CHKERRQ(ierr)
      call KSPGetIterationNumber(ksp,iter,ierr) ;CHKERRQ(ierr)
      call KSPGetConvergedReason(ksp,reason,ierr) ;CHKERRQ(ierr)

      if(myid.eq.0.and.ldebug) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      !      if(reason.eq.KSP_DIVERGED_ITS) then
      !        if(myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
      !        return
      !      endif

      if(reason.le.0) then
        if(myid.eq.0.and.ldebug) print *,myid,'Resetted initial guess to zero and try again with gmres:'
        call VecSet(x,zero,ierr) ;CHKERRQ(ierr)
        call KSPGetType(ksp,old_ksp_type,ierr); CHKERRQ(ierr)
        call KSPSetType(ksp,KSPGMRES,ierr) ;CHKERRQ(ierr)
        call KSPSetUp(ksp,ierr) ;CHKERRQ(ierr)
        call KSPSolve(ksp,b,x,ierr) ;CHKERRQ(ierr)
        call KSPGetIterationNumber(ksp,iter,ierr) ;CHKERRQ(ierr)
        call KSPGetConvergedReason(ksp,reason,ierr) ;CHKERRQ(ierr)

        ! And return to normal solver...
        call KSPSetType(ksp,old_ksp_type,ierr) ;CHKERRQ(ierr)
        call KSPSetFromOptions(ksp,ierr) ;CHKERRQ(ierr)
        call KSPSetUp(ksp,ierr) ;CHKERRQ(ierr)
        if(myid.eq.0.and.ldebug) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      endif

      if(reason.le.0) then
        if(myid.eq.0) print *,'***** SOLVER did NOT converge :( ********',reason
        call exit()
      endif
  end subroutine

  subroutine setup_ksp(ksp,C,A,linit, prefix)
      KSP :: ksp
      type(t_coord) :: C
      Mat :: A
      PC  :: prec
      logical :: linit

      !      MatNullSpace :: nullspace
      !      Vec :: nullvecs(0)
      character(len=*),optional :: prefix

      PetscReal,parameter :: rtol=1e-5_ireals, atol=1e-8_ireals
      PetscInt,parameter  :: maxiter=10000

      PetscInt,parameter :: ilu_default_levels=1
      PetscInt :: pcbjac_n_local, pcbjac_iglob ! number of local ksp contexts and index in global ksp-table
      KSP,allocatable :: pcbjac_ksps(:)
      PC  :: pcbjac_sub_pc
      integer(iintegers) :: isub

      if(linit) return
      call PetscLogStagePush(logstage(8),ierr) ;CHKERRQ(ierr)
      if(myid.eq.0.and.ldebug) &
          print *,'Setup KSP -- tolerances:',rtol,atol*(C%dof*C%glob_xm*C%glob_ym*C%glob_zm) * count(.not.atm%l1d)/(one*size(atm%l1d))

      call KSPCreate(imp_comm,ksp,ierr) ;CHKERRQ(ierr)
      if(present(prefix) ) call KSPAppendOptionsPrefix(ksp,trim(prefix),ierr) ;CHKERRQ(ierr)

      call KSPSetType(ksp,KSPBCGS,ierr)  ;CHKERRQ(ierr)
      call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr) ;CHKERRQ(ierr) 
      call KSPGetPC  (ksp,prec,ierr)  ;CHKERRQ(ierr)
      if(numnodes.eq.0) then
        call PCSetType (prec,PCILU,ierr);CHKERRQ(ierr)
      else
        call PCSetType (prec,PCBJACOBI,ierr);CHKERRQ(ierr)
      endif

      call KSPSetTolerances(ksp,rtol,max(1e-30_ireals, atol*(C%dof*C%glob_xm*C%glob_ym*C%glob_zm) * count(.not.atm%l1d)/(one*size(atm%l1d)) ),PETSC_DEFAULT_REAL,maxiter,ierr);CHKERRQ(ierr)

      call KSPSetConvergenceTest(ksp,MyKSPConverged, PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)

      call KSPSetOperators(ksp,A,A,ierr) ;CHKERRQ(ierr)
      call KSPSetDM(ksp,C%da,ierr) ;CHKERRQ(ierr)
      call KSPSetDMActive(ksp,PETSC_FALSE,ierr) ;CHKERRQ(ierr)

      call KSPSetUp(ksp,ierr) ;CHKERRQ(ierr)

      if(numnodes.eq.0) then
        call PCFactorSetLevels(prec,ilu_default_levels,ierr);CHKERRQ(ierr)
      else
        call PCBJacobiGetSubKSP(prec,pcbjac_n_local,pcbjac_iglob,PETSC_NULL_OBJECT,ierr);CHKERRQ(ierr)
        if(.not.allocated(pcbjac_ksps)) allocate(pcbjac_ksps(pcbjac_n_local))
        call PCBJacobiGetSubKSP(prec,pcbjac_n_local,pcbjac_iglob,pcbjac_ksps,ierr);CHKERRQ(ierr)

        do isub=1,pcbjac_n_local
          call KSPSetType(pcbjac_ksps(isub) ,KSPPREONLY,ierr)              ;CHKERRQ(ierr)
          call KSPGetPC  (pcbjac_ksps(isub), pcbjac_sub_pc,ierr)        ;CHKERRQ(ierr)
          call PCSetType (pcbjac_sub_pc, PCILU,ierr)                    ;CHKERRQ(ierr)
          call PCFactorSetLevels(pcbjac_sub_pc,ilu_default_levels,ierr) ;CHKERRQ(ierr)
        enddo

      endif

      call KSPSetFromOptions(ksp,ierr) ;CHKERRQ(ierr)

      !      call MatNullSpaceCreate( imp_comm, PETSC_TRUE, PETSC_NULL_INTEGER, nullvecs, nullspace, ierr) ; CHKERRQ(ierr)
      !      call KSPSetNullspace(ksp, nullspace, ierr) ; CHKERRQ(ierr)

      linit = .True.
      if(myid.eq.0.and.ldebug) print *,'Setup KSP done'
      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
  end subroutine
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

  subroutine setup_logging()
      logical,save :: logstage_init=.False.
      if(logstage_init) return
      call PetscLogStageRegister('total_tenstream' , logstage(1)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('setup_edir'      , logstage(2)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('calc_edir'       , logstage(3)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('setup_ediff'     , logstage(4)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('calc_ediff'      , logstage(5)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('setup_b'         , logstage(6)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('get_coeff'       , logstage(7)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('twostream'       , logstage(8)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('setup_ksp'       , logstage(9)     , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('write_hdf5'      , logstage(10)    , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('load_save_sol'   , logstage(11)    , ierr) ;CHKERRQ(ierr)
      call PetscLogStageRegister('nca'             , logstage(12)    , ierr) ;CHKERRQ(ierr)

      if(myid.eq.0 .and. ldebug) print *, 'Logging stages' , logstage
      logstage_init=.True.
  end subroutine

  subroutine nca_wrapper(ediff,abso)
      use m_ts_nca, only : ts_nca
      Vec :: ediff,abso
      Vec :: lediff ! local ediff vector with ghost values -- in dimension 0 and 1 are fluxes followed by dz,planck,kabs

      PetscReal,pointer,dimension(:,:,:,:) :: xv  =>null()
      PetscReal,pointer,dimension(:)       :: xv1d=>null()
      integer(iintegers) :: k

      real(ireals),allocatable,dimension(:,:,:) :: dz_g, planck_g, kabs_g, hr, Edn_g, Eup_g

      call PetscLogStagePush(logstage(12),ierr) ;CHKERRQ(ierr)

      ! put additional values into a new ediff vec .. TODO: this is a rather dirty hack but is straightforward
      call getVecPointer(ediff ,C_diff ,xv1d, xv)
      xv(  i2, C_diff%zs:C_diff%ze-1, :,: ) = atm%dz
      xv(  i2, C_diff%ze            , :,: ) = zero
                                         
      xv(  i3, C_diff%zs:C_diff%ze  , :,: ) = atm%planck
                                         
      xv(  i4, C_diff%zs:C_diff%ze-1, :,: ) = atm%delta_op%kabs
      xv(  i4, C_diff%ze            , :,: ) = zero
      call restoreVecPointer(ediff ,C_diff ,xv1d, xv )


      ! get ghost values for dz, planck, kabs and fluxes and pack em in new arrays, ready to give it to NCA
      call DMGetLocalVector(C_diff%da ,lediff ,ierr)                   ; CHKERRQ(ierr)
      call VecSet(lediff,zero,ierr); CHKERRQ(ierr)

      call DMGlobalToLocalBegin(C_diff%da ,ediff ,ADD_VALUES,lediff ,ierr) ; CHKERRQ(ierr)
      call DMGlobalToLocalEnd(C_diff%da ,ediff ,ADD_VALUES,lediff ,ierr)   ; CHKERRQ(ierr)

      allocate ( dz_g     ( C_one%zs :C_one%ze ,  C_one%gxs : C_one%gxe ,  C_one%gys : C_one%gye  ) )
      allocate ( kabs_g   ( C_one%zs :C_one%ze ,  C_one%gxs : C_one%gxe ,  C_one%gys : C_one%gye  ) )
      allocate ( planck_g ( C_diff%zs:C_diff%ze,  C_diff%gxs: C_diff%gxe,  C_diff%gys: C_diff%gye ) )
      allocate ( Edn_g    ( C_diff%zs:C_diff%ze,  C_diff%gxs: C_diff%gxe,  C_diff%gys: C_diff%gye ) )
      allocate ( Eup_g    ( C_diff%zs:C_diff%ze,  C_diff%gxs: C_diff%gxe,  C_diff%gys: C_diff%gye ) )
      allocate ( hr       ( C_one%zs: C_one%ze ,  C_one%xs  : C_one%xe  ,  C_one%ys  : C_one%ye   ) )

      call getVecPointer(lediff ,C_diff ,xv1d, xv)
      do k= C_one%zs, C_one%ze 
        dz_g    (k,:,:) = xv(   i2, k , :,:)
        kabs_g  (k,:,:) = xv(   i4, k , :,:)
      enddo
      do k= C_diff%zs, C_diff%ze 
        Edn_g   (k,:,:) = xv( E_dn, k , :,:)
        Eup_g   (k,:,:) = xv( E_up, k , :,:)
        planck_g(k,:,:) = xv(   i3, k , :,:)
      enddo

      call restoreVecPointer(lediff ,C_diff ,xv1d, xv )

      call DMRestoreLocalVector(C_diff%da, lediff, ierr); CHKERRQ(ierr)

      hr=0
      call ts_nca(C_one%zm, C_one%xm, C_one%ym, dz_g, planck_g, kabs_g, hr, atm%dx, atm%dy, Edn_g, Eup_g)

      call getVecPointer( abso, C_one ,xv1d, xv)
      do k=C_one%zs,C_one%ze 
        xv(i0,k,:,:) = hr(k,:,:) / dz_g(k,:,:)
      enddo
      call restoreVecPointer( abso ,C_one ,xv1d, xv )

      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
  end subroutine

  subroutine twostream(edirTOA,solution)
      real(ireals),intent(in) :: edirTOA
      type(t_state_container) :: solution

      PetscReal,pointer,dimension(:,:,:,:) :: xv_dir=>null(),xv_diff=>null()
      PetscReal,pointer,dimension(:) :: xv_dir1d=>null(),xv_diff1d=>null()
      integer(iintegers) :: i,j,src

      real(ireals),allocatable :: dtau(:),kext(:),w0(:),g(:),S(:),Edn(:),Eup(:)
      real(ireals) :: mu0,incSolar

      call VecSet(solution%edir ,zero,ierr); CHKERRQ(ierr)
      call VecSet(solution%ediff,zero,ierr); CHKERRQ(ierr)

      call PetscLogStagePush(logstage(8),ierr) ;CHKERRQ(ierr)

      allocate( dtau(C_dir%zm-1) )
      allocate( kext(C_dir%zm-1) )
      allocate(   w0(C_dir%zm-1) )
      allocate(    g(C_dir%zm-1) )

      mu0 = sun%costheta
      incSolar = edirTOA* sun%costheta

      call getVecPointer(solution%edir  ,C_dir  ,xv_dir1d , xv_dir  ,.False.)
      call getVecPointer(solution%ediff ,C_diff ,xv_diff1d, xv_diff ,.False.)

      allocate( S(C_dir%zm ) )
      allocate( Eup(C_diff%zm) )
      allocate( Edn(C_diff%zm) )

      if(myid.eq.0) print *,' CALCULATING DELTA EDDINGTON TWOSTREAM ::',sun%theta,':',incSolar

      do j=C_dir%ys,C_dir%ye         
        do i=C_dir%xs,C_dir%xe

          kext = atm%op(:,i,j)%kabs + atm%op(:,i,j)%ksca
          dtau = atm%dz(:,i,j)* kext
          w0   = atm%op(:,i,j)%ksca / kext
          g    = atm%op(:,i,j)%g

          if(allocated(atm%planck) ) then
            call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo, S,Edn,Eup, planck=atm%planck(:,i,j) )
          else
            call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo, S,Edn,Eup )
          endif

          do src=i0,i3
            xv_dir(src,:,i,j) = S(:)
          enddo

          xv_diff(E_up,:,i,j) = Eup(:) 
          xv_diff(E_dn,:,i,j) = Edn(:) 
        enddo
      enddo

      call restoreVecPointer(solution%edir  ,C_dir  ,xv_dir1d , xv_dir  )
      call restoreVecPointer(solution%ediff ,C_diff ,xv_diff1d, xv_diff )

      !Twostream solver returns fluxes as [W]
      solution%lintegrated_dir  = .False.
      solution%lintegrated_diff = .False.
      ! and mark solution that it is not up to date
      solution%lchanged         = .False.

      deallocate(S)
      deallocate(Edn)
      deallocate(Eup)

      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
  end subroutine

  subroutine scale_flx(solution, lWm2_to_W)
    type(t_state_container) :: solution
    logical,intent(in) :: lWm2_to_W ! determines direction of scaling, if true, scale from W/m**2 to W 

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
    subroutine scale_flx_vec(v,C,lWm2_to_W)
        Vec :: v
        type(t_coord) :: C
        PetscReal,pointer,dimension(:,:,:,:) :: xv  =>null()
        PetscReal,pointer,dimension(:)       :: xv1d=>null()
        PetscInt :: i,j,k
        PetscReal :: Ax,Ax2,Ay,Ay2,Az,Az4
        logical,intent(in) :: lWm2_to_W ! determines direction of scaling, if true, scale from W/m**2 to W

        if(myid.eq.0.and.ldebug) print *,'rescaling fluxes'
        call getVecPointer(v ,C ,xv1d, xv ,.False. )

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

              if(.not.atm%l1d(k,i,j)) then

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

        call restoreVecPointer(v ,C ,xv1d, xv )
    end subroutine
end subroutine

  subroutine set_diff_initial_guess(inp,guess,C)
      ! Deprecated -- this is probably not helping convergence....
      Vec :: inp,guess
      type(t_coord) :: C

      Vec :: local_guess
      PetscScalar,pointer,dimension(:,:,:,:) :: xinp=>null(),xguess=>null()
      PetscScalar,pointer,dimension(:) :: xinp1d=>null(),xguess1d=>null()
      PetscReal :: diff2diff(C_diff%dof**2)
      PetscInt :: k,i,j,src

      if(myid.eq.0.and.ldebug) print *,'setting initial guess...'
      call DMGetLocalVector(C%da,local_guess,ierr) ;CHKERRQ(ierr)
      call VecSet(local_guess,zero,ierr) ;CHKERRQ(ierr)

      call getVecPointer(inp ,C ,xinp1d, xinp ,.False.)
      call getVecPointer(local_guess ,C ,xguess1d, xguess ,.True.)

      do j=C%ys,C%ye         
        do i=C%xs,C%xe    
          do k=C%zs,C%ze-1 
            if( .not. atm%l1d(k,i,j) ) then
              call get_coeff(atm%delta_op(k,i,j), atm%dz(k,i,j),.False., diff2diff, atm%l1d(k,i,j) )
              do src=1,C%dof
                xguess(E_up   , k   , i   , j   ) = xguess(E_up   , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_up  +i1+(src-1 ) *C%dof )
                xguess(E_dn   , k+1 , i   , j   ) = xguess(E_dn   , k+1 , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_dn  +i1+(src-1 ) *C%dof )
                xguess(E_le_m , k   , i   , j   ) = xguess(E_le_m , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_le_m+i1+(src-1 ) *C%dof )
                xguess(E_le_p , k   , i   , j   ) = xguess(E_le_p , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_le_p+i1+(src-1 ) *C%dof )
                xguess(E_ri_m , k   , i+1 , j   ) = xguess(E_ri_m , k   , i+1 , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ri_m+i1+(src-1 ) *C%dof )
                xguess(E_ri_p , k   , i+1 , j   ) = xguess(E_ri_p , k   , i+1 , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ri_p+i1+(src-1 ) *C%dof )
                xguess(E_ba_m , k   , i   , j   ) = xguess(E_ba_m , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ba_m+i1+(src-1 ) *C%dof )
                xguess(E_ba_p , k   , i   , j   ) = xguess(E_ba_p , k   , i   , j   ) +  xinp(src-1 , k , i , j ) *diff2diff(E_ba_p+i1+(src-1 ) *C%dof )
                xguess(E_fw_m , k   , i   , j+1 ) = xguess(E_fw_m , k   , i   , j+1 ) +  xinp(src-1 , k , i , j ) *diff2diff(E_fw_m+i1+(src-1 ) *C%dof )
                xguess(E_fw_p , k   , i   , j+1 ) = xguess(E_fw_p , k   , i   , j+1 ) +  xinp(src-1 , k , i , j ) *diff2diff(E_fw_p+i1+(src-1 ) *C%dof )
              enddo
            endif
          enddo
        enddo
      enddo

      call restoreVecPointer(inp ,C ,xinp1d, xinp )
      call restoreVecPointer(local_guess ,C ,xguess1d, xguess )

      call VecSet(guess,zero,ierr) ;CHKERRQ(ierr) ! reset global Vec

      call DMLocalToGlobalBegin(C%da,local_guess,ADD_VALUES, guess,ierr) ;CHKERRQ(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
      call DMLocalToGlobalEnd  (C%da,local_guess,ADD_VALUES, guess,ierr) ;CHKERRQ(ierr)

      call DMRestoreLocalVector(C_diff%da, local_guess, ierr); CHKERRQ(ierr)
  end subroutine

  subroutine init_memory(incSolar,b,Mdir,Mdiff)
      Vec :: b,incSolar
      Mat :: Mdiff,Mdir

      call DMCreateGlobalVector(C_dir%da,incSolar,ierr) ; CHKERRQ(ierr)
      call DMCreateGlobalVector(C_diff%da,b,ierr)       ; CHKERRQ(ierr)

      call VecSet(incSolar,zero,ierr) ; CHKERRQ(ierr)
      call VecSet(b,zero,ierr)        ; CHKERRQ(ierr)

      call init_Matrix(Mdir ,C_dir )!,"dir_")
      call init_Matrix(Mdiff,C_diff)!,"diff_")
  end subroutine

  subroutine init_tenstream(icomm, Nz,Nx,Ny, dx,dy, phi0,theta0,albedo, dz1d, dz3d, nxproc, nyproc)
      ! vector sizes Nx, Ny Nz are either global domain size or if present(nxproc,nyproc) local domain size
      integer,intent(in) :: icomm
      integer(iintegers),intent(in) :: Nz,Nx,Ny
      real(ireals),intent(in) :: dx,dy, phi0,theta0,albedo

      real(ireals),optional :: dz1d(Nz), dz3d(Nz,Nx,Ny)   ! dz3d has to be local domain size
      integer(iintegers),optional :: nxproc(:), nyproc(:) ! size of local domains on each node

      integer(iintegers) :: k,i,j
      !    character(len=30),parameter :: tenstreamrc='./.tenstreamrc'

      if(.not.linitialized) then

        call setup_petsc_comm
        !      call PetscInitialize(tenstreamrc ,ierr) ;CHKERRQ(ierr)
        call PetscInitialize(PETSC_NULL_CHARACTER ,ierr) ;CHKERRQ(ierr)
#ifdef _XLF
        call PetscPopSignalHandler(ierr); CHKERRQ(ierr) ! in case of xlf ibm compilers, remove petsc signal handler -- otherwise we dont get fancy signal traps from boundschecking or FPE's
#endif

        call init_mpi_data_parameters(PETSC_COMM_WORLD)

        call read_commandline_options()

        if(present(nxproc) .and. present(nyproc) ) then
          if(ldebug.and.myid.eq.0) print *,'nxproc',shape(nxproc),'::',nxproc
          if(ldebug.and.myid.eq.0) print *,'nyproc',shape(nyproc),'::',nyproc
          call setup_grid( Nz, Nx, Ny, nxproc,nyproc)
        else
          call setup_grid( Nz, max(minimal_dimension, Nx), max(minimal_dimension, Ny) )
        endif

      endif

      atm%dx  = dx
      atm%dy  = dy
      atm%albedo = albedo

      if(.not.allocated(atm%dz) ) allocate(atm%dz( C_one%zs:C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye ))

      if(present(dz1d)) then
        do j=C_one%ys,C_one%ye
          do i=C_one%xs,C_one%xe
            atm%dz(:,i,j) = dz1d
          enddo
        enddo
      else if(present(dz3d)) then
        atm%dz = dz3d
      else
        print *,'have to give either dz1d or dz3d in routine call....'
        call exit(1)
      endif

      if(.not.allocated(atm%l1d)) then
        allocate(atm%l1d( C_one%zs:C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye ) )
        atm%l1d = .False.
      endif

      do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe
          do k=C_one%ze-1,C_one%zs,-1
            if( atm%l1d(k,i,j) ) cycle ! if it was already marked 1D before, we must not change it to 3D -- otherwise have to recreate matrix routines. possible but not implemented at the moment !TODO

            atm%l1d(C_one%ze,i,j) = twostr_ratio*atm%dz(C_one%ze,i,j).gt.atm%dx

            if( atm%l1d(k+1,i,j) ) then !can only be 3D RT if below is a 3D layer
              atm%l1d(k,i,j)=.True.
            else
              !TODO this does actually not really make sense. I am not sure why this might not work. i.e. if it is necessary.
              atm%l1d(k,i,j) = twostr_ratio*atm%dz(k,i,j).gt.atm%dx
              if(atm%dz(k,i,j).lt.atm%dx/10._ireals) atm%l1d(k,i,j)=.True.
            endif
          enddo
        enddo
      enddo

      call setup_suninfo(phi0,theta0,sun)

      ! init box montecarlo model
      if(any(atm%l1d.eqv..False.)) call OPP_8_10%init(atm%dx,atm%dy,[sun%symmetry_phi],[sun%theta],imp_comm)
      if(.not.luse_eddington)      call OPP_1_2%init (atm%dx,atm%dy,[sun%symmetry_phi],[sun%theta],imp_comm) 

      if(.not.linitialized) then
        ! init matrices & vectors
        call init_memory(incSolar,b,Mdir,Mdiff)

        ! init petsc logging facilities
        call setup_logging()
      endif

      linitialized=.True.
    contains
      subroutine setup_petsc_comm()
          !This is the code snippet from Petsc FAQ to change from PETSC (C) domain splitting to MPI(Fortran) domain splitting 
          ! the numbers of processors per direction are (int) x_procs, y_procs, z_procs respectively
          ! (no parallelization in direction 'dir' means dir_procs = 1)

          !        MPI_Comm :: NewComm
          !        PetscInt :: x,y
          !
          !        integer(mpiint) :: orig_id,new_id,petsc_id,ierr ! id according to fortran decomposition

          PETSC_COMM_WORLD = icomm

          !        if(present(nxproc) .and. present(nyproc) ) then
          !          call MPI_COMM_RANK( icomm, orig_id, ierr )
          !
          !          ! calculate coordinates of cpus in MPI ordering:
          !          x = int(orig_id) / size(nyproc)
          !          y = modulo(orig_id ,size(nyproc))
          !
          !          ! set new rank according to PETSc ordering:
          !          petsc_id = y*size(nxproc) + x
          !
          !          ! create communicator with new ranks according to PETSc ordering:
          !          call MPI_Comm_split(MPI_COMM_WORLD, i1, petsc_id, NewComm, ierr)
          !
          !          ! override the default communicator (was MPI_COMM_WORLD as default)
          !          PETSC_COMM_WORLD = NewComm
          !        endif
          !        print *,'setup_petsc_comm: MPI_COMM_WORLD',orig_id,'calc_id',petsc_id,'PETSC_COMM_WORLD',new_id

      end subroutine

  end subroutine

  subroutine set_global_optical_properties(global_kabs, global_ksca, global_g, global_planck)
      real(ireals),intent(inout),dimension(:,:,:),allocatable,optional :: global_kabs, global_ksca, global_g
      real(ireals),intent(inout),dimension(:,:,:),allocatable,optional :: global_planck
      real(ireals),dimension(:,:,:),allocatable :: local_kabs, local_ksca, local_g
      real(ireals),dimension(:,:,:),allocatable :: local_planck
      logical :: lhave_planck,lhave_kabs,lhave_ksca,lhave_g

      if(.not.linitialized) then
        print *,myid,'You tried to set global optical properties but tenstream environment seems not to be initialized.... please call init first!'
        call exit(1)
      endif

      lhave_kabs   = present(global_kabs  ); call imp_bcast( lhave_kabs  , 0_mpiint, myid )
      lhave_ksca   = present(global_ksca  ); call imp_bcast( lhave_ksca  , 0_mpiint, myid )
      lhave_g      = present(global_g     ); call imp_bcast( lhave_g     , 0_mpiint, myid )
      lhave_planck = present(global_planck); call imp_bcast( lhave_planck, 0_mpiint, myid )

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
        call set_optical_properties(local_kabs, local_ksca, local_g, local_planck)
      else
        call set_optical_properties(local_kabs, local_ksca, local_g)
      endif


    contains
      subroutine local_optprop()
          Vec :: local_vec
          PetscReal,pointer,dimension(:,:,:,:) :: xlocal_vec  =>null()
          PetscReal,pointer,dimension(:)       :: xlocal_vec1d=>null()

          if(myid.eq.0.and.ldebug .and. lhave_kabs) &
              print *,myid,'copying optprop: global to local :: shape kabs',shape(global_kabs),'xstart/end',C_one%xs,C_one%xe,'ys/e',C_one%ys,C_one%ye

          call DMGetGlobalVector(C_one%da, local_vec, ierr) ; CHKERRQ(ierr)

          if(lhave_kabs) then
            call scatterZerotoDM(global_kabs,C_one,local_vec)
            call getVecPointer(local_vec ,C_one ,xlocal_vec1d, xlocal_vec,.False. )
            local_kabs = xlocal_vec(0, C_one%zs :C_one%ze,C_one%xs :C_one%xe , C_one%ys :C_one%ye  )
            call restoreVecPointer(local_vec ,C_one ,xlocal_vec1d, xlocal_vec )
          endif

          if(lhave_ksca) then
            call scatterZerotoDM(global_ksca,C_one,local_vec)
            call getVecPointer(local_vec ,C_one ,xlocal_vec1d, xlocal_vec,.False. )
            local_ksca = xlocal_vec(0,:,:,:)
            call restoreVecPointer(local_vec ,C_one ,xlocal_vec1d, xlocal_vec )
          endif

          if(lhave_g) then
            call scatterZerotoDM(global_g,C_one,local_vec)
            call getVecPointer(local_vec ,C_one ,xlocal_vec1d, xlocal_vec,.False. )
            local_g = xlocal_vec(0,:,:,:)
            call restoreVecPointer(local_vec ,C_one ,xlocal_vec1d, xlocal_vec )
          endif

          call DMRestoreGlobalVector(C_one%da, local_vec, ierr) ; CHKERRQ(ierr)

          if(lhave_planck) then
            call DMGetGlobalVector(C_one1%da, local_vec, ierr) ; CHKERRQ(ierr)
            call scatterZerotoDM(global_planck,C_one1,local_vec)
            call getVecPointer(local_vec ,C_one1 ,xlocal_vec1d, xlocal_vec,.False. )
            local_planck = xlocal_vec(0,:,:,:)
            call restoreVecPointer(local_vec ,C_one1 ,xlocal_vec1d, xlocal_vec )
            call DMRestoreGlobalVector(C_one1%da, local_vec, ierr) ; CHKERRQ(ierr)
          endif
      end subroutine
  end subroutine

  subroutine scatterZerotoDM(arr,C,vec)
      real(ireals),allocatable,dimension(:,:,:),intent(in) :: arr
      type(t_coord),intent(in) :: C
      Vec :: vec

      VecScatter :: scatter_context
      Vec :: natural,local
      PetscScalar,Pointer :: xloc(:)=>null()

      if(ldebug) print *,myid,'scatterZerotoDM :: start....'
      call VecSet(vec,zero,ierr)

      call DMDACreateNaturalVector(C%da, natural, ierr); CHKERRQ(ierr)
      call VecScatterCreateToZero(natural, scatter_context, local, ierr); CHKERRQ(ierr)

      if(myid.eq.0) then
        call VecGetArrayF90(local,xloc,ierr) ;CHKERRQ(ierr)
        if(ldebug) &
            print *,myid,'scatterZerotoDM :: shape of local',shape(xloc), 'shape of arr',shape(arr)
        xloc = reshape( arr , [ size(arr) ] )
        call VecRestoreArrayF90(local,xloc,ierr) ;CHKERRQ(ierr)
      endif

      if(ldebug) print *,myid,'scatterZerotoDM :: scatter reverse....'
      call VecScatterBegin(scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); CHKERRQ(ierr)
      call VecScatterEnd  (scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); CHKERRQ(ierr)

      if(ldebug) print *,myid,'scatterZerotoDM :: natural to global....'
      call DMDANaturalToGlobalBegin(C%da,natural, INSERT_VALUES, vec, ierr); CHKERRQ(ierr)
      call DMDANaturalToGlobalEnd  (C%da,natural, INSERT_VALUES, vec, ierr); CHKERRQ(ierr)

      if(ldebug) print *,myid,'scatterZerotoDM :: destroying contexts....'
      call VecScatterDestroy(scatter_context, ierr); CHKERRQ(ierr)
      call VecDestroy(local,ierr); CHKERRQ(ierr)
      call VecDestroy(natural,ierr); CHKERRQ(ierr)
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

  subroutine set_optical_properties(local_kabs, local_ksca, local_g, local_planck)
      real(ireals),intent(in),dimension(:,:,:),optional :: local_kabs, local_ksca, local_g
      real(ireals),intent(in),dimension(:,:,:),optional :: local_planck
      real(ireals) :: tau,kext,w0,g
      integer(iintegers) :: k,i,j

      if(.not.allocated(atm%op) )  allocate( atm%op       (C_one%zs :C_one%ze, C_one%xs:C_one%xe, C_one%ys:C_one%ye) )
      if(present(local_kabs) ) atm%op(:,:,:)%kabs = local_kabs
      if(present(local_ksca) ) atm%op(:,:,:)%ksca = local_ksca
      if(present(local_g   ) ) atm%op(:,:,:)%g    = local_g   

      if(present(local_planck) ) then 
        if (.not.allocated(atm%planck) ) allocate( atm%planck   (C_one1%zs:C_one1%ze, C_one1%xs:C_one1%xe, C_one1%ys:C_one1%ye ) )
        atm%planck = local_planck
      else
        if(allocated(atm%planck)) deallocate(atm%planck)
      endif

      if(ldebug) then
        if( (any([local_kabs,local_ksca,local_g].lt.zero)) .or. (any(isnan([local_kabs,local_ksca,local_g]))) ) then
          print *,myid,'set_optical_properties :: found illegal value in local_optical properties! abort!'
          do k=C_one%zs,C_one%ze
            print *,myid,k,'local_kabs',local_kabs(k,:,:)
            print *,myid,k,'local_ksca',local_ksca(k,:,:)
          enddo
        endif
      endif
      if(ldebug) then
        if( (any([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g].lt.zero)) .or. (any(isnan([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g]))) ) then
          print *,myid,'set_optical_properties :: found illegal value in optical properties! abort!'
        endif
      endif

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

      ! Make space for deltascaled optical properties
      if(.not.allocated(atm%delta_op) ) allocate( atm%delta_op (C_one%zs :C_one%ze,C_one%xs :C_one%xe , C_one%ys :C_one%ye ) )
      atm%delta_op = atm%op
      call delta_scale(atm%delta_op(:,:,:)%kabs, atm%delta_op(:,:,:)%ksca, atm%delta_op(:,:,:)%g ) !todo should we instead use strong deltascaling? -- what gives better results? or is it as good?

      if(ldebug) then
        if( (any([atm%delta_op(:,:,:)%kabs,atm%delta_op(:,:,:)%ksca,atm%delta_op(:,:,:)%g].lt.zero)) .or. (any(isnan([atm%delta_op(:,:,:)%kabs,atm%delta_op(:,:,:)%ksca,atm%delta_op(:,:,:)%g]))) ) then
          print *,myid,'set_optical_properties :: found illegal value in delta_scaled optical properties! abort!'
        endif
      endif

      if(luse_eddington) then
        if(.not.allocated(atm%a11) ) allocate(atm%a11 (C_one%zs:C_one%ze ,C_one%xs:C_one%xe, C_one%ys:C_one%ye))   ! allocate space for twostream coefficients
        if(.not.allocated(atm%a12) ) allocate(atm%a12 (C_one%zs:C_one%ze ,C_one%xs:C_one%xe, C_one%ys:C_one%ye)) 
        if(.not.allocated(atm%a13) ) allocate(atm%a13 (C_one%zs:C_one%ze ,C_one%xs:C_one%xe, C_one%ys:C_one%ye)) 
        if(.not.allocated(atm%a23) ) allocate(atm%a23 (C_one%zs:C_one%ze ,C_one%xs:C_one%xe, C_one%ys:C_one%ye)) 
        if(.not.allocated(atm%a33) ) allocate(atm%a33 (C_one%zs:C_one%ze ,C_one%xs:C_one%xe, C_one%ys:C_one%ye)) 
        if(.not.allocated(atm%g1 ) ) allocate(atm%g1  (C_one%zs:C_one%ze ,C_one%xs:C_one%xe, C_one%ys:C_one%ye)) 
        if(.not.allocated(atm%g2 ) ) allocate(atm%g2  (C_one%zs:C_one%ze ,C_one%xs:C_one%xe, C_one%ys:C_one%ye)) 
      endif

      if(luse_eddington) then
        do j=C_one%ys,C_one%ye
          do i=C_one%xs,C_one%xe
            do k=C_one%zs,C_one%ze
              if( atm%l1d(k,i,j) ) then
                kext = atm%delta_op(k,i,j)%kabs + atm%delta_op(k,i,j)%ksca
                w0   = atm%delta_op(k,i,j)%ksca / kext
                tau  = atm%dz(k,i,j)* kext
                g    = atm%delta_op(k,i,j)%g 
                call eddington_coeff_zdun ( tau , w0, g, sun%costheta, & 
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
                endif
              endif
            enddo
          enddo
        enddo
      endif

      if(ldebug .and. myid.eq.0) then
        do k=C_one%zs,C_one%ze
          if(present(local_planck)) then
            print *,myid,'Optical Properties:',k,'dz',atm%dz(k,0,0),atm%l1d(k,0,0),'k',&
                minval(atm%delta_op    (k,:,:)%kabs),minval(atm%delta_op (k, :,:)%ksca),minval(atm%delta_op (k, :,:)%g),&
                maxval(atm%delta_op    (k,:,:)%kabs),maxval(atm%delta_op (k, :,:)%ksca),maxval(atm%delta_op (k, :,:)%g),&
                '::',minval(atm%planck (k,:,:)),maxval(atm%planck        (k, :,:))
          else    
            print *,myid,'Optical Properties:',k,'dz',atm%dz(k,0,0),atm%l1d(k,0,0),'k',&
                minval(atm%delta_op (k, :,:)%kabs),minval(atm%delta_op (k,:,:)%ksca),minval(atm%delta_op (k,:,:)%g),&
                maxval(atm%delta_op (k, :,:)%kabs),maxval(atm%delta_op (k,:,:)%ksca),maxval(atm%delta_op (k,:,:)%g),&
                '::',minval(atm%a33 (k, :,:)),maxval(atm%a33           (k,:,:))
          endif
        enddo
      endif

  end subroutine

  subroutine solve_tenstream(edirTOA,opt_solution_uid,opt_solution_time)
      real(ireals),intent(in) :: edirTOA
      integer(iintegers),optional,intent(in) :: opt_solution_uid
      real(ireals),      optional,intent(in) :: opt_solution_time
      logical :: loaded
      integer(iintegers) :: uid,twostr_uid

      if(present(opt_solution_uid)) then
        uid = opt_solution_uid
      else
        uid = i0 ! default solution is uid==0
      endif

      call prepare_solution( solutions(uid), uid, lsolar=(edirTOA.gt.zero .and. sun%theta.ge.zero) ) ! setup solution vectors

      ! --------- Can we get an initial guess? ---------------
      if(ltwostr) then
        if(ltwostr_only) then
          twostr_uid = uid
        else
          twostr_uid = -1
        endif

        call prepare_solution(solutions(twostr_uid),twostr_uid, lsolar=.True. ) ! make space for twostream solution 
        call twostream(edirTOA,  solutions(twostr_uid) )
        solutions(twostr_uid)%lchanged=.True.
        if(ldebug .and. myid.eq.0) print *,'twostream calculation done'

        if(present(opt_solution_time) ) then 
          call restore_solution(solutions(twostr_uid),opt_solution_time)
        else
          call restore_solution(solutions(twostr_uid))
        endif

        if(ltwostr_only) return
      endif

      ! --------- scale from [W/m**2] to [W] -----------------
      call scale_flx(solutions(uid), lWm2_to_W=.True. )

      ! ---------------------------- Edir  -------------------
      if( solutions(uid)%lsolar_rad ) then

        call PetscLogStagePush(logstage(1),ierr) ;CHKERRQ(ierr)
        call setup_incSolar(incSolar,edirTOA)
        call set_dir_coeff(Mdir,C_dir)

        call setup_ksp(kspdir,C_dir,Mdir,linit_kspdir, "dir_")

        call PetscLogStagePush(logstage(3),ierr) ;CHKERRQ(ierr)
        call solve(kspdir,incSolar,solutions(uid)%edir)
        solutions(uid)%lchanged=.True.
        solutions(uid)%lintegrated_dir=.True.
        call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
      endif

      ! ---------------------------- Source Term -------------
      call setup_b(solutions(uid),b)

      ! ---------------------------- Ediff -------------------
      call set_diff_coeff(Mdiff,C_diff)

      call setup_ksp(kspdiff,C_diff,Mdiff,linit_kspdiff, "diff_")

      call PetscLogStagePush(logstage(5),ierr) ;CHKERRQ(ierr)

      call solve(kspdiff, b, solutions(uid)%ediff,uid)
      solutions(uid)%lchanged=.True.
      solutions(uid)%lintegrated_diff=.True. !Tenstream solver returns fluxes as [W]

      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

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

      if(linitialized) then 

        call KSPDestroy(kspdir , ierr) ;CHKERRQ(ierr); linit_kspdir =.False.
        call KSPDestroy(kspdiff, ierr) ;CHKERRQ(ierr); linit_kspdiff=.False.

        call VecDestroy(incSolar , ierr) ;CHKERRQ(ierr)
        call VecDestroy(b        , ierr) ;CHKERRQ(ierr)
        call MatDestroy(Mdir     , ierr) ;CHKERRQ(ierr)
        call MatDestroy(Mdiff    , ierr) ;CHKERRQ(ierr)

        do uid=lbound(solutions,1),ubound(solutions,1)
          if( solutions(uid)%lset ) then
            if(solutions(uid)%lsolar_rad) &
                call VecDestroy(solutions(uid)%edir     , ierr) ;CHKERRQ(ierr)
            call VecDestroy(solutions(uid)%ediff    , ierr) ;CHKERRQ(ierr)
            call VecDestroy(solutions(uid)%abso     , ierr) ;CHKERRQ(ierr)
          endif
        enddo

        deallocate(atm%op)
        deallocate(atm%delta_op)
        if(allocated(atm%planck)) deallocate(atm%planck)
        deallocate(atm%a11)
        deallocate(atm%a12)
        deallocate(atm%a13)
        deallocate(atm%a23)
        deallocate(atm%a33)
        deallocate(atm%dz)
        deallocate(atm%l1d)
        call OPP_1_2%destroy()
        call OPP_8_10%destroy()

        deallocate(C_dir%neighbors)  ; call DMDestroy(C_dir%da ,ierr)
        deallocate(C_diff%neighbors) ; call DMDestroy(C_diff%da,ierr)
        deallocate(C_one%neighbors)  ; call DMDestroy(C_one%da ,ierr)
        deallocate(C_one1%neighbors) ; call DMDestroy(C_one1%da,ierr)

        linitialized=.False.

        if(lfinalize) call PetscFinalize(ierr) ;CHKERRQ(ierr)
      endif
  end subroutine

  subroutine tenstream_get_result(redir,redn,reup,rabso, opt_solution_uid )
      real(ireals),dimension(:,:,:),intent(inout),allocatable :: redir,redn,reup,rabso
      integer(iintegers),optional,intent(in) :: opt_solution_uid

      integer(iintegers) :: uid
      PetscScalar,pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()


      if(present(opt_solution_uid)) then
        uid = opt_solution_uid
      else
        uid = i0 ! default solution is uid==0
      endif

      if(ldebug .and. myid.eq.0) print *,'calling tenstream_get_result',allocated(redir),allocated(redn),allocated(reup),allocated(rabso),'for uid',uid

      if(solutions(uid)%lchanged) stop 'tried to get results from unrestored solution -- call restore_solution first'

      if(allocated(redir)) then
        if( .not. solutions(uid)%lsolar_rad ) then
          redir = zero
        else
          if( solutions(uid)%lintegrated_dir ) stop 'tried to get result from integrated result vector(dir)'
          call getVecPointer(solutions(uid)%edir,C_dir,x1d,x4d)
          redir = sum(x4d(i0:i3,:,:,:),dim=1)/4 ! return average of the 4 vertical tiles
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

      if(allocated(redn).or.allocated(reup)) then
        if(solutions(uid)%lintegrated_diff) stop 'tried to get result from integrated result vector(diff)'
        call getVecPointer(solutions(uid)%ediff,C_diff,x1d,x4d,.False.)
        if(allocated(redn) )redn = x4d(i1,:,:,:)
        if(allocated(reup) )reup = x4d(i0,:,:,:)
        if(ldebug) then
          if(myid.eq.0) print *,' Edn',redn(1,1,:)
          if(myid.eq.0) print *,' Eup',reup(1,1,:)
          if(allocated(redir).and.allocated(redn)) then
            if(any(redn.lt.-one)) then 
              print *,'Found direct radiation smaller than 0 in edn result... that should not happen',minval(redn)
              call exit(1)
            endif
          endif
          if(allocated(redir).and.allocated(reup)) then
            if(any(reup.lt.-one)) then 
              print *,'Found direct radiation smaller than 0 in eup result... that should not happen',minval(reup)
              call exit(1)
            endif
          endif
        endif!ldebug
        call restoreVecPointer(solutions(uid)%ediff,C_diff,x1d,x4d)
      endif

      if(allocated(rabso)) then
        call getVecPointer(solutions(uid)%abso,C_one,x1d,x4d)
        rabso = x4d(i0,:,:,:)
        call restoreVecPointer(solutions(uid)%abso,C_one,x1d,x4d)
      endif
  end subroutine

  subroutine prepare_solution(solution, uid, lsolar)
    type(t_state_container) :: solution
    integer(iintegers),intent(in) :: uid
    logical,intent(in) :: lsolar

    solutions(uid)%uid = uid ! dirty hack to give the solution a unique hash for example to write it out to disk

    if( lsolar .and.  (solution%lsolar_rad.eqv..False.) ) then 
      ! set up the direct vectors in any case. This may be necessary even if solution was
      ! already initialized once... e.g. in case we calculated thermal before with same uid
      call DMCreateGlobalVector(C_dir%da ,solution%edir  ,ierr)  ; CHKERRQ(ierr)
      call VecSet(solution%edir,zero,ierr); CHKERRQ(ierr)
      solution%lsolar_rad = .True.
    endif

    if(ldebug .and. myid.eq.0) print *,'Prepare solution',uid,'::',lsolar,solution%lsolar_rad,'set?',solution%lset

    if(solution%lset) return ! already set-up 

    call DMCreateGlobalVector(C_diff%da,solution%ediff ,ierr)  ; CHKERRQ(ierr)
    call DMCreateGlobalVector(C_one%da ,solution%abso  ,ierr)  ; CHKERRQ(ierr)

    call VecSet(solution%ediff,zero,ierr)    ; CHKERRQ(ierr)
    call VecSet(solution%abso,zero,ierr)     ; CHKERRQ(ierr)

    solution%lset = .True.
  end subroutine

  function need_new_solution(uid,time)
      integer(iintegers),intent(in) :: uid
      real(ireals),intent(in) :: time
      logical :: need_new_solution

      integer,parameter :: Nfit=4 ! Number of used residuals
      real(ireals) :: t(Nfit),tm(Nfit),dt(Nfit-1),e(Nfit-1),integ_err(Nfit), error_estimate, polyc(3)

      character(len=30) :: reason
      integer, parameter :: out_unit=20

      integer(iintegers) :: k

      if( .not. solutions(uid)%lset ) then !if we did not store a solution, return immediately
        need_new_solution=.True.
        write(reason,*) 'no solution yet' 
        return 
      endif

      call PetscLogStagePush(logstage(11),ierr) ;CHKERRQ(ierr)
      do k=1,Nfit
        t(k) = solutions(uid)%time(Nfit-k+1)
      enddo
      do k=1,Nfit-1
        e(k) = solutions(uid)%maxnorm(Nfit-k)
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

      ! cumsum(e*dt) is the integral over error
      integ_err(1:Nfit-1) = cumsum(e*dt)
      polyc = polyfit(t(2:Nfit),integ_err(1:Nfit-1),size(polyc)-i1, ierr)
      if(ierr.ne.0) then 
        need_new_solution=.True.
        write(reason,*) 'problem fitting error curve'
        call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
        return
      endif

      ! Use fit coefficients to calculate estimate of error integral at t=time
      integ_err(Nfit)=0
      do k=1,size(polyc)
        integ_err(Nfit) = integ_err(Nfit) + polyc(k)*time**(k-1)
      enddo

      error_estimate = abs( integ_err(Nfit)-integ_err(Nfit-1) )/ max(epsilon(time), time - t(Nfit))

      if(myid.eq.0 .and. error_estimate.le.zero) then
        print *,'DEBUG t',t
        print *,'DEBUG e',e
        print *,'DEBUG tm',tm
        print *,'DEBUG dt',dt
        print *,'DEBUG integ_err',integ_err
        print *,'DEBUG err_est',error_estimate
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

      !            if(.not.need_new_solution) then
      !              need_new_solution=.True. ! overwrite it and calculate anyway
      !              write(reason,*) 'MANUAL OVERRIDE' 
      !            endif

      !            if(ldebug .and. myid.eq.0 .and. .not. need_new_solution ) &
      if(ldebug.and. myid.eq.0) then
        print *,''
        print *,''
        print *,'new calc',need_new_solution,' bc ',reason,' t',time,uid,' residuals _solver ::', solutions(uid)%ksp_residual_history(1:4),'    ::     est.',error_estimate,'[W]',error_estimate*86.1,'[K/d]'
        if(uid.eq.2) then
          open (unit=out_unit,file="residuals.log",action="write",status="replace")
        else
          open (unit=out_unit,file="residuals.log",action="readwrite",status="unknown",position = "append")
        endif

        write (out_unit,*) uid,solutions(uid)%maxnorm
        write (out_unit,*) uid,solutions(uid)%time
        write (out_unit,*) uid,solutions(uid)%twonorm
        close (out_unit)
      endif
      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

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
            if(myid.eq.0) print *, "problem with lapack lsqr :: 1 :: info",info
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

  subroutine restore_solution(solution,time)
  ! restore_solution:: if flux have changed, we need to update absorption, save the residual history
      type(t_state_container) :: solution
      real(ireals),intent(in),optional :: time

      character(100) :: vecname
      real(ireals) :: norm1,norm2,norm3
      Vec :: abso_old

      call PetscLogStagePush(logstage(11),ierr) ;CHKERRQ(ierr)

      if( .not. solution%lset ) &
          stop 'cant restore solution that was not initialized'

      if( .not. solution%lchanged ) &
          stop 'cant restore solution which was not changed'

      if(present(time) ) then ! Create working vec to determine difference between old and new absorption vec
        call DMGetGlobalVector(C_one%da, abso_old, ierr) ; CHKERRQ(ierr) 
        call VecCopy( solution%abso, abso_old, ierr)     ; CHKERRQ(ierr)
      endif

      ! make sure to bring the fluxes into [W] for absorption calculation
      call scale_flx(solution, lWm2_to_W=.True. ) 

      ! update absorption
      call calc_flx_div(solution)

      if(present(time) ) then ! Compute norm between old absorption and new one
        call VecAXPY(abso_old , -one, solution%abso , ierr)    ; CHKERRQ(ierr) ! overwrite abso_old with difference to new one
        call VecNorm(abso_old ,  NORM_1, norm1, ierr)          ; CHKERRQ(ierr)
        call VecNorm(abso_old ,  NORM_2, norm2, ierr)          ; CHKERRQ(ierr)
        call VecNorm(abso_old ,  NORM_INFINITY, norm3, ierr)   ; CHKERRQ(ierr)

        call DMRestoreGlobalVector(C_one%da, abso_old, ierr)   ; CHKERRQ(ierr) 

        ! Save norm for later analysis
        solution%maxnorm = eoshift ( solution%maxnorm, shift = -1) !shift all values by 1 to the right
        solution%twonorm = eoshift ( solution%twonorm, shift = -1) !shift all values by 1 to the right
        solution%time    = eoshift ( solution%time   , shift = -1) !shift all values by 1 to the right

        solution%maxnorm( 1 ) = norm3 
        solution%twonorm( 1 ) = norm2
        solution%time( 1 )    = time

        !if(myid.eq.0) &
        if(ldebug .and. myid.eq.0) &
            print *,'Updating error statistics for solution ',solution%uid,'at time ',time,'::',solution%time(1),':: norm',norm1,norm2,norm3,'[W] :: hr_norm approx:',norm3*86.1,'[K/d]'

      endif !present(time)

      if(ldebug .and. myid.eq.0) &
          print *,'Saving Solution ',solution%uid

      ! make sure to bring the fluxes into [W/m**2]
      call scale_flx(solution, lWm2_to_W=.False. ) 

      if(ldebug .and. myid.eq.0) &
          print *,'Saving Solution done'
      solution%lchanged=.False.

      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

      if(lwriteall) then
        if(solution%lsolar_rad) then
          write(vecname,FMT='("edir",I0)') solution%uid
          call PetscObjectSetName(solution%edir,vecname,ierr) ; CHKERRQ(ierr)
          call vec_to_hdf5(solution%edir)
        endif

        write(vecname,FMT='("ediff",I0)') solution%uid
        call PetscObjectSetName(solution%ediff,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(solution%ediff)

        write(vecname,FMT='("abso",I0)') solution%uid
        call PetscObjectSetName(solution%abso,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(solution%abso)

        write(vecname,FMT='("b",I0)') solution%uid
        call PetscObjectSetName(b,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(b)

        write(vecname,FMT='("incSolar",I0)') solution%uid
        call PetscObjectSetName(incSolar,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(incSolar)
      endif
  end subroutine

  subroutine getVecPointer(vec,C,x1d,x4d,local)
      Vec :: vec
      type(t_coord),intent(in) :: C
      PetscScalar,intent(inout),pointer,dimension(:,:,:,:) :: x4d
      PetscScalar,intent(inout),pointer,dimension(:) :: x1d
      logical,optional,intent(in) :: local

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

      call VecGetArrayF90(vec,x1d,ierr) ;CHKERRQ(ierr)
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
      call VecRestoreArrayF90(vec,x1d,ierr) ;CHKERRQ(ierr)
      x1d => null()
  end subroutine

subroutine vec_to_hdf5(v)
      Vec,intent(in) :: v
      character(10),parameter :: suffix='.h5'
      character(110) :: fname
      logical fexists
      PetscFileMode :: fmode
      character(100) :: vecname
      
      PetscViewer :: view

      call PetscObjectGetName(v,vecname,ierr) ;CHKERRQ(ierr)

      fname = 'vecdump' // trim(suffix)
      inquire(file=trim(fname), exist=fexists)
      
      if(fexists) then
        if(myid.eq.0 .and. ldebug)  print *,myid,'appending vector to hdf5 file ',trim(fname),' vecname: ',vecname
        fmode = FILE_MODE_APPEND
      else 
        if(myid.eq.0 .and. ldebug)  print *,myid,'writing vector to hdf5 file ',trim(fname),' vecname: ',vecname
        fmode = FILE_MODE_WRITE
      endif

      call PetscViewerHDF5Open(imp_comm,trim(fname),fmode, view, ierr) ;CHKERRQ(ierr)
      call VecView(v, view, ierr) ;CHKERRQ(ierr)
      call PetscViewerDestroy(view,ierr) ;CHKERRQ(ierr)

      if(myid.eq.0 .and. ldebug ) print *,myid,'writing to hdf5 file done'
end subroutine
end module
