module m_tenstream

      use m_data_parameters, only : ireals,iintegers,       &
        imp_comm, myid, numnodes,init_mpi_data_parameters,mpiint, &
        zero,one,nil,i0,i1,i2,i3,i4,i5,i6,i7,i8,i10,pi

      use m_twostream, only: delta_eddington_twostream
      use m_helper_functions, only: deg2rad,approx,rmse,delta_scale,imp_bcast
      use m_eddington, only : eddington_coeff_fab
      use m_optprop_parameters, only : ldelta_scale
      use m_optprop, only : t_optprop_1_2,t_optprop_8_10
      use m_tenstream_options, only : read_commandline_options, ltwostr, luse_twostr_guess, luse_eddington, twostr_ratio

      implicit none
#include "finclude/petsc.h90"

      private
      public :: init_tenstream, set_global_optical_properties, set_optical_properties, solve_tenstream, destroy_tenstream,&
                tenstream_get_result, &
                b,edir,ediff,abso,&
                edir_twostr,ediff_twostr,abso_twostr, &
                t_coord,C_dir,C_diff,C_one

      PetscInt,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9
      PetscInt,parameter :: istartpar=i1, jstartpar=i1

      logical,parameter :: ldebug=.False.,lcycle_dir=.True.
      logical,parameter :: lprealloc=.True.

      type t_coord
        PetscInt :: xs,xe,ys,ye,zs,ze     ! local domain start and end indices
        PetscInt :: xm,ym,zm              ! size of local domain
        PetscInt :: gxs,gys,gzs           ! domain indices including ghost points
        PetscInt :: gxm,gym,gzm           ! size of local domain including ghosts
        PetscInt :: glob_xm,glob_ym,glob_zm ! global domain size
        PetscInt :: dof,dim               ! degrees of freedom of Petsc Domain, dimension of dmda
        DM :: da                          ! The Domain Decomposition Object
        PetscMPIInt,allocatable :: neighbors(:) ! all 3d neighbours( (x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
      end type

      type(t_coord) :: C_dir,C_diff,C_one,C_one1

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
      type(t_atmosphere) :: atm


      type(t_optprop_1_2) OPP_1_2
      type(t_optprop_8_10) OPP_8_10

      type t_suninfo
        real(ireals) :: symmetry_phi
        integer(iintegers) :: yinc,xinc
        real(ireals) :: theta,phi,costheta,sintheta
      end type
      type(t_suninfo) :: sun

      PetscLogStage :: logstage(10)

      Mat :: Mdir,Mdiff

      Vec :: incSolar,b,edir,ediff,abso
      Vec :: edir_twostr,ediff_twostr,abso_twostr

      KSP :: kspdir, kspdiff
      logical :: linit_kspdir=.False., linit_kspdiff=.False.

      logical :: linitialized=.False.

      integer(iintegers),parameter :: minimal_dimension=3 ! this is the minimum number of gridpoints in x or y direction

      logical,parameter :: lenable_solutions=.True. ! if enabled, we can save and load solutions.... just pass an unique identifer to solve()... beware, this may use lots of memory
      type t_state_container
        Vec :: edir,ediff
        logical :: lset=.False.
      end type
      type(t_state_container),save :: solutions(1000)

      contains 
      subroutine setup_grid(Nx,Ny,Nz,nxproc,nyproc)
        PetscInt,intent(in) :: Nx,Ny,Nz
        integer(iintegers),optional :: nxproc(:), nyproc(:) ! size of local domains on each node

        DMBoundaryType :: bp=DM_BOUNDARY_PERIODIC, bn=DM_BOUNDARY_NONE, bg=DM_BOUNDARY_GHOSTED

        if(myid.eq.0.and.ldebug) print *,myid,'Setting up the DMDA grid for ',Nx,Ny,Nz,'using ',numnodes,' nodes'
        
        if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_diff'
        call setup_dmda(C_diff,Nx,Ny, Nz+1, bp, i10)

        if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_dir'
        if(lcycle_dir) then
          call setup_dmda(C_dir,Nx,Ny, Nz+1, bp, i8)
        else
          call setup_dmda(C_dir,Nx,Ny, Nz+1, bg, i8)
        endif

        if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C1'
        call setup_dmda(C_one,Nx,Ny, Nz, bp, i1)
        call setup_dmda(C_one1,Nx,Ny, Nz+1, bp, i1)

        if(myid.eq.0.and.ldebug) print *,myid,'DMDA grid ready'
        contains
          subroutine setup_dmda(C, Nx, Ny, Nz, boundary, dof)
              type(t_coord) :: C
              PetscInt,intent(in) :: Nx,Ny,Nz,dof

              DMBoundaryType :: boundary
              PetscInt,parameter :: stencil_size=1

              C%dof = i1*dof

              if(present(nxproc) .and. present(nyproc) ) then
                call DMDACreate3d( imp_comm  ,                                           &
                    boundary           , boundary           , bn                 , &
                    DMDA_STENCIL_BOX  ,                                            &
                    sum(nxproc)     , sum(nyproc)     , i1*Nz                , &
                    size(nxproc)          , size(nyproc)          , i1                 , &
                    C%dof              , stencil_size       ,                      &
                    nxproc , nyproc , Nz , &
                    C%da               , ierr) ;CHKERRQ(ierr)

              else
                call DMDACreate3d( imp_comm  ,                                           &
                    boundary           , boundary           , bn                 , &
                    DMDA_STENCIL_BOX  ,                                            &
                    Nx     , Ny     , i1*Nz                , &
                    PETSC_DECIDE          , PETSC_DECIDE          , i1                 , &
                    C%dof              , stencil_size       ,                      &
                    PETSC_NULL_INTEGER , PETSC_NULL_INTEGER , PETSC_NULL_INTEGER , &
                    C%da               , ierr) ;CHKERRQ(ierr)
              endif
              call setup_coords(C)
              call DMSetup(C%da,ierr) ;CHKERRQ(ierr)
!              call DMSetMatType(C%da, MATBAIJ, ierr); CHKERRQ(ierr)
              if(lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ;CHKERRQ(ierr)

              call DMSetFromOptions(C%da, ierr) ; CHKERRQ(ierr)
          end subroutine
        subroutine setup_coords(C)
          type(t_coord) :: C

          call DMDAGetInfo(C%da,C%dim,                             &
            C%glob_xm,C%glob_ym,C%glob_zm,                           &
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
            ierr) ;CHKERRQ(ierr)

          call DMDAGetCorners(C%da,C%xs,C%ys,C%zs, C%xm,C%ym,C%zm, ierr) ;CHKERRQ(ierr)
          call DMDAGetGhostCorners(C%da,C%gxs,C%gys,C%gzs,C%gxm,C%gym,C%gzm,ierr) ;CHKERRQ(ierr)
          C%xe = C%xs+C%xm-1
          C%ye = C%ys+C%ym-1
          C%ze = C%zs+C%zm-1
          if(ldebug) print *,myid,'Domain Corners x:: ',C%xs,':',C%xe,' (',C%xm,' entries)','global size',C%glob_xm
          if(ldebug) print *,myid,'Domain Corners y:: ',C%ys,':',C%ye,' (',C%ym,' entries)','global size',C%glob_ym
          if(ldebug) print *,myid,'Domain Corners z:: ',C%zs,':',C%ze,' (',C%zm,' entries)','global size',C%glob_zm

          allocate(C%neighbors(0:3**C%dim-1) )
          call DMDAGetNeighbors(C%da,C%neighbors,ierr) ;CHKERRQ(ierr)
          call DMSetUp(C%da,ierr) ;CHKERRQ(ierr)
          if(ldebug.and.C%dim.eq.3) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([10,12,16,14]),'while I am ',C%neighbors(13)
          if(ldebug.and.C%dim.eq.2) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([1,3,7,5]),'while I am ',C%neighbors(4)
        end subroutine
      end subroutine

      subroutine mat_info(A)
        Mat :: A
        double precision :: info(MAT_INFO_SIZE)
        double precision :: mal, nz_allocated, nz_used, nz_unneeded

        call MatGetInfo(A,MAT_LOCAL,info,ierr) ;CHKERRQ(ierr)
        mal = info(MAT_INFO_MALLOCS)
        nz_allocated = info(MAT_INFO_NZ_ALLOCATED)
        nz_used   = info(MAT_INFO_NZ_USED)
        nz_unneeded = info(MAT_INFO_NZ_UNNEEDED)

        if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_MALLOCS',mal,'MAT_INFO_NZ_ALLOCATED',nz_allocated
        if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_USED',nz_used,'MAT_INFO_NZ_unneded',nz_unneeded

      end subroutine
      subroutine init_Matrix(A,C,prefix)
        Mat :: A
        type(t_coord) :: C
        character(len=*),optional :: prefix

        PetscInt,dimension(:),allocatable :: o_nnz,d_nnz!,dnz

        call DMCreateMatrix(C%da, A, ierr) ;CHKERRQ(ierr)
        if(present(prefix) ) then
            !call MatAppendOptionsPrefix(A,trim(prefix),ierr) !!! Not in the Fortran API?
        endif

!        call MatCreate(PETSC_COMM_WORLD, A,ierr); CHKERRQ(ierr)
!        call MatSetSizes(A, C%xm*C%ym*C%zm*C%dof, C%xm*C%ym*C%zm*C%dof, C%glob_xm*C%glob_ym*C%glob_zm*C%dof, C%glob_xm*C%glob_ym*C%glob_zm*C%dof,ierr);CHKERRQ(ierr)
!        call MatSetType(A, MATAIJ, ierr); CHKERRQ(ierr)

        if(lprealloc) then
          ! Determine perfect preallocation
          if(numnodes.gt.1) then
            select case (C%dof)
            case(i3)
              call setup_dir_preallocation(d_nnz,o_nnz,C)
              !             call MatMPIAIJSetPreallocation(A, C%dof,PETSC_NULL_INTEGER, i2, PETSC_NULL_INTEGER, ierr) ;CHKERRQ(ierr) !TODO
            case(i8)
              call setup_dir8_preallocation(d_nnz,o_nnz,C)
              !             call MatMPIAIJSetPreallocation(A, C%dof,PETSC_NULL_INTEGER, i4, PETSC_NULL_INTEGER, ierr) ;CHKERRQ(ierr) !TODO
            case(i10)
              call setup_diff_preallocation(d_nnz,o_nnz,C)
              !             call MatMPIAIJSetPreallocation(A, C%dof,PETSC_NULL_INTEGER, i4, PETSC_NULL_INTEGER, ierr) ;CHKERRQ(ierr) !TODO
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

        call MatSetFromOptions(A,ierr) ;CHKERRQ(ierr)
        call MatSetUp(A,ierr) ;CHKERRQ(ierr)

!        call MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr) ;CHKERRQ(ierr)
!        call MatSetOption(A,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr) ;CHKERRQ(ierr)
!        call MatSetOption(A,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE,ierr) ;CHKERRQ(ierr)

        call mat_info(A)
        call mat_set_diagonal(A,C)
      end subroutine
      subroutine mat_set_diagonal(A,C)
        Mat :: A
        type(t_coord),intent(in) :: C
        PetscInt :: i,j,k,dof
        MatStencil :: row(4,1), col(4,1)
        PetscScalar :: v(1)

        !TODO -- we should use this form... however get an allocation ERROR - fix this
!        Vec :: diag

!        call DMCreateGlobalVector(C%da,diag,ierr) ;CHKERRQ(ierr)
!        call VecSet(diag,one,ierr) ;CHKERRQ(ierr)
!        call MatDiagonalSet( A, diag, INSERT_VALUES,ierr) ;CHKERRQ(ierr)

        v(1)=one
        
        if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients diagonally'
        do k=C%zs,C%ze
          row(MatStencil_k,1) = k
          col(MatStencil_k,1) = k

          do j=C%ys,C%ye 
            row(MatStencil_j,1) = j
            col(MatStencil_j,1) = j

            do i=C%xs,C%xe
              row(MatStencil_i,1) = i
              col(MatStencil_i,1) = i

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
      subroutine setup_diff_preallocation(d_nnz,o_nnz,C)
        PetscInt,allocatable :: o_nnz(:)
        PetscInt,allocatable :: d_nnz(:)
        type(t_coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize

        PetscScalar, pointer :: xx_v_o(:),xx_v_d(:)

        PetscInt,parameter :: ind(9)=[E_up,E_le_m,E_le_p,E_ri_m,E_ri_p,E_ba_m,E_ba_p,E_fw_m,E_fw_p]
        PetscInt :: i,j,k!d,li,lj,lk

        call DMCreateGlobalVector(C%da,v_o_nnz,ierr) ;CHKERRQ(ierr)
        call DMCreateGlobalVector(C%da,v_d_nnz,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C%da,v_o_nnz,xo,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C%da,v_d_nnz,xd,ierr) ;CHKERRQ(ierr)

        xo = i0
        xd = i1
        xd(E_up,:,:,C%ze) = i2 ! E_up at ground is only function of E_dn at grnd and +1 for self

        xd(ind,:,:,C%zs:C%ze-1) = C%dof+i1
        xd(E_dn,:,:,C%zs+1:C%ze) = C%dof+i1

! Determine prealloc for E_dn and E_up
        if( C%neighbors(14).ne.myid .and. C%neighbors(14).ge.i0 ) then ! neigh east
                ! E_le is not local
                xd(E_dn,C%xe,:,C%zs+1:C%ze) = xd(E_dn,C%xe,:,C%zs+1:C%ze) - i2
                xo(E_dn,C%xe,:,C%zs+1:C%ze) = xo(E_dn,C%xe,:,C%zs+1:C%ze) + i2
                xd(E_up,C%xe,:,C%zs:C%ze-1) = xd(E_up,C%xe,:,C%zs:C%ze-1) - i2
                xo(E_up,C%xe,:,C%zs:C%ze-1) = xo(E_up,C%xe,:,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0) then ! neigh north
                ! E_ba is not local
                xd(E_dn,:,C%ye,C%zs+1:C%ze) = xd(E_dn,:,C%ye,C%zs+1:C%ze) - i2
                xo(E_dn,:,C%ye,C%zs+1:C%ze) = xo(E_dn,:,C%ye,C%zs+1:C%ze) + i2
                xd(E_up,:,C%ye,C%zs:C%ze-1) = xd(E_up,:,C%ye,C%zs:C%ze-1) - i2
                xo(E_up,:,C%ye,C%zs:C%ze-1) = xo(E_up,:,C%ye,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0) then ! neigh south
                ! no foreign stream dependencies
        endif
        if( C%neighbors(12).ne.myid .and. C%neighbors(12).ge.i0) then ! neigh west
                ! no foreign stream dependencies
        endif
! Determine prealloc for E_le
        if( C%neighbors(14).ne.myid .and. C%neighbors(14).ge.i0 ) then ! neigh east
                ! E_le is not local
                xd([E_le_m,E_le_p],C%xe,:,C%zs:C%ze-1) = xd([E_le_m,E_le_p],C%xe,:,C%zs:C%ze-1) - i2
                xo([E_le_m,E_le_p],C%xe,:,C%zs:C%ze-1) = xo([E_le_m,E_le_p],C%xe,:,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh north
                ! E_ba is not local
                xd([E_le_m,E_le_p],:,C%ye,C%zs:C%ze-1) = xd([E_le_m,E_le_p],:,C%ye,C%zs:C%ze-1) - i2
                xo([E_le_m,E_le_p],:,C%ye,C%zs:C%ze-1) = xo([E_le_m,E_le_p],:,C%ye,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                ! no foreign stream dependencies
        endif
        if( C%neighbors(12).ne.myid .and. C%neighbors(12).ge.i0 ) then ! neigh west
                ! no foreign stream dependencies
        endif
! Determine prealloc for E_ri
        if( C%neighbors(14).ne.myid .and. C%neighbors(14).ge.i0 ) then ! neigh east
                ! no foreign stream dependencies
        endif
        if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh north
                ! E_ba is not local
                xd([E_ri_m,E_ri_p],:,C%ye,C%zs:C%ze-1) = xd([E_ri_m,E_ri_p],:,C%ye,C%zs:C%ze-1) - i2
                xo([E_ri_m,E_ri_p],:,C%ye,C%zs:C%ze-1) = xo([E_ri_m,E_ri_p],:,C%ye,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                ! no foreign stream dependencies
        endif
        if( C%neighbors(12).ne.myid .and. C%neighbors(12).ge.i0 ) then ! neigh west
                ! E_ri local dependencies are only self, and 2*E_le
                xd([E_ri_m,E_ri_p],C%xs,:,C%zs:C%ze-1) = i3
                xo([E_ri_m,E_ri_p],C%xs,:,C%zs:C%ze-1) = i8
        endif
! Determine prealloc for E_ba
        if( C%neighbors(14).ne.myid .and. C%neighbors(14).ge.i0 ) then ! neigh east
                ! E_le is not local
                xd([E_ba_m,E_ba_p],C%xe,:,C%zs:C%ze-1) = xd([E_ba_m,E_ba_p],C%xe,:,C%zs:C%ze-1) - i2
                xo([E_ba_m,E_ba_p],C%xe,:,C%zs:C%ze-1) = xo([E_ba_m,E_ba_p],C%xe,:,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh north
                ! E_ba is not local
                xd([E_ba_m,E_ba_p],:,C%ye,C%zs:C%ze-1) = xd([E_ba_m,E_ba_p],:,C%ye,C%zs:C%ze-1) - i2
                xo([E_ba_m,E_ba_p],:,C%ye,C%zs:C%ze-1) = xo([E_ba_m,E_ba_p],:,C%ye,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                ! no foreign stream dependencies
        endif
        if( C%neighbors(12).ne.myid .and. C%neighbors(12).ge.i0 ) then ! neigh west
                ! no foreign stream dependencies
        endif
! Determine prealloc for E_fw
        if( C%neighbors(14).ne.myid .and. C%neighbors(14).ge.i0 ) then ! neigh east
                ! E_le is not local
                xd([E_fw_m,E_fw_p],C%xe,:,C%zs:C%ze-1) = xd([E_fw_m,E_fw_p],C%xe,:,C%zs:C%ze-1) - i2
                xo([E_fw_m,E_fw_p],C%xe,:,C%zs:C%ze-1) = xo([E_fw_m,E_fw_p],C%xe,:,C%zs:C%ze-1) + i2
        endif
        if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh north
                ! no foreign stream dependencies
        endif
        if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                ! E_fw local dependencies are only self, and 2*E_ba
                xd([E_fw_m,E_fw_p],:,C%ys,C%zs:C%ze-1) = i3
                xo([E_fw_m,E_fw_p],:,C%ys,C%zs:C%ze-1) = i8
        endif
        if( C%neighbors(12).ne.myid .and. C%neighbors(12).ge.i0 ) then ! neigh west
                ! no foreign stream dependencies
        endif

        do k=C%zs,C%ze-1
          do j=C%ys,C%ye
            do i=C%xs,C%xe      
              if( atm%l1d(i,j,k) ) then

                xo(:,i,j,k) = i0
                xd(:,i,j,k) = i1
                xd(0:1,i,j,k) = i3

              endif
            enddo
          enddo
        enddo

        call DMDAVecRestoreArrayF90(C%da,v_o_nnz,xo,ierr) ;CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C%da,v_d_nnz,xd,ierr) ;CHKERRQ(ierr)
        
        call VecGetLocalSize(v_o_nnz,vsize,ierr) ;CHKERRQ(ierr)
        allocate(o_nnz(0:vsize-1))
        allocate(d_nnz(0:vsize-1))

        call VecGetArrayF90(v_o_nnz,xx_v_o,ierr) ;CHKERRQ(ierr)
        call VecGetArrayF90(v_d_nnz,xx_v_d,ierr) ;CHKERRQ(ierr)
        o_nnz=int(xx_v_o)
        d_nnz=int(xx_v_d)
        call VecRestoreArrayF90(v_o_nnz,xx_v_o,ierr) ;CHKERRQ(ierr)
        call VecRestoreArrayF90(v_d_nnz,xx_v_d,ierr) ;CHKERRQ(ierr)

        call VecDestroy(v_o_nnz,ierr) ;CHKERRQ(ierr)
        call VecDestroy(v_d_nnz,ierr) ;CHKERRQ(ierr)

        if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
      end subroutine 
      subroutine setup_dir8_preallocation(d_nnz,o_nnz,C)
        PetscInt,allocatable :: d_nnz(:)
        PetscInt,allocatable :: o_nnz(:)
        type(t_coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize,i,j,k

        PetscScalar, pointer :: xx_v_o(:),xx_v_d(:)

        logical :: lsun_east,lsun_north

        if(myid.eq.0.and.ldebug) print *,myid,'building direct o_nnz for mat with',C%dof,'dof'
        call DMCreateGlobalVector(C%da,v_o_nnz,ierr) ;CHKERRQ(ierr)
        call DMCreateGlobalVector(C%da,v_d_nnz,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C%da,v_o_nnz,xo,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C%da,v_d_nnz,xd,ierr) ;CHKERRQ(ierr)

        xo = i0
        xd = i1
        xd( i0:i3 ,:,:,C%zs+1 :C%ze   ) = C%dof+i1 ! Edir_vertical depends on 3 values Edir_vertical,xaxis,yaxis :: starting with second entries(seen from top)
        xd( i4:i7 ,:,:,C%zs   :C%ze-i1) = C%dof+i1 ! Edir_xaxis,yaxis depends on 3 values Edir_vertical,xaxis,yaxis :: starting with first entries(seen from top)

        do j=C%ys,C%ye
                 lsun_east  = (sun%xinc.eq.i0)

                 if( C%neighbors(14).ne.myid .and. C%neighbors(14).ge.i0 ) then ! neigh east
                         if( lsun_east ) then 
                                ! if the sun is in the east, the channels in the last box are influenced by the 2nd channel which is a ghost
                                xo(i0:i3, C%xe, j, C%zs+1:C%ze) = xo(i0:i3, C%xe, j, C%zs+1:C%ze)+i2 ! channel 1 from zs+1 to ze
                                xd(i0:i3, C%xe, j, C%zs+1:C%ze) = xd(i0:i3, C%xe, j, C%zs+1:C%ze)-i2

                                xo(i4:i7, C%xe, j, C%zs:C%ze-1) = xo(i4:i7, C%xe, j, C%zs:C%ze-1)+i2 ! channel 2 and 3 from zs
                                xd(i4:i7, C%xe, j, C%zs:C%ze-1) = xd(i4:i7, C%xe, j, C%zs:C%ze-1)-i2
                        endif
                endif
        enddo
        do i=C%xs,C%xe
                 lsun_north = (sun%yinc.eq.i0 )

                 if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh north
                         if( lsun_north ) then 
                                ! if the sun is in the north, the 3rd channel is a ghost
                                xo(i0:i3, i, C%ye, C%zs+1:C%ze) = xo(i0:i3, i, C%ye, C%zs+1:C%ze)+i2 ! channel 1 from zs+1 to ze
                                xd(i0:i3, i, C%ye, C%zs+1:C%ze) = xd(i0:i3, i, C%ye, C%zs+1:C%ze)-i2

                                xo(i4:i7,  i, C%ye, C%zs:C%ze-1) = xo(i4:i7, i, C%ye, C%zs:C%ze-1)+i2 ! channel 2 and 3 from zs
                                xd(i4:i7,  i, C%ye, C%zs:C%ze-1) = xd(i4:i7, i, C%ye, C%zs:C%ze-1)-i2
                        endif
                endif
        enddo
        do j=C%ys,C%ye
                 lsun_east  = (sun%xinc.eq.i0)

                 if( C%neighbors(12).ne.myid.and. C%neighbors(12).ge.i0 ) then ! neigh west
                         if( .not. lsun_east ) then 
                                ! if the sun is in the west, the 2nd channel is solemnly dependant on ghost values
                                xo(i4:i5, C%xs, j, C%zs:C%ze-1) = C%dof
                                xd(i4:i5, C%xs, j, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo
        do i=C%xs,C%xe
                 lsun_north = (sun%yinc.eq.i0 )

                 if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                         if( .not. lsun_north ) then 
                                ! if the sun is in the south, the 3rd channel is solemnly dependant on ghost values
                                xo(i6:i7, i, C%ys, C%zs:C%ze-1) = C%dof
                                xd(i6:i7, i, C%ys, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo

         do k=C%zs,C%ze-1
           do j=C%ys,C%ye
             do i=C%xs,C%xe      
               if( atm%l1d(i,j,k) ) then
                 xo(:,i,j,k) = i0
                 xd(:,i,j,k) = i1
                 xd(0:3,i,j,k) = i5
               endif
             enddo
           enddo
         enddo

        call DMDAVecRestoreArrayF90(C%da,v_o_nnz,xo,ierr) ;CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C%da,v_d_nnz,xd,ierr) ;CHKERRQ(ierr)
        
        call VecGetLocalSize(v_o_nnz,vsize,ierr) ;CHKERRQ(ierr)
        allocate(o_nnz(0:vsize-1))
        allocate(d_nnz(0:vsize-1))

        call VecGetArrayF90(v_o_nnz,xx_v_o,ierr) ;CHKERRQ(ierr)
        call VecGetArrayF90(v_d_nnz,xx_v_d,ierr) ;CHKERRQ(ierr)
        o_nnz=int(xx_v_o)
        d_nnz=int(xx_v_d)
        call VecRestoreArrayF90(v_o_nnz,xx_v_o,ierr) ;CHKERRQ(ierr)
        call VecRestoreArrayF90(v_d_nnz,xx_v_d,ierr) ;CHKERRQ(ierr)

        call VecDestroy(v_o_nnz,ierr) ;CHKERRQ(ierr)
        call VecDestroy(v_d_nnz,ierr) ;CHKERRQ(ierr)

        if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
      end subroutine 
      subroutine setup_dir_preallocation(d_nnz,o_nnz,C)
        PetscInt,allocatable :: d_nnz(:)
        PetscInt,allocatable :: o_nnz(:)
        type(t_coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize,i,j,k

        PetscScalar, pointer :: xx_v_o(:),xx_v_d(:)

        logical :: lsun_east,lsun_north

!        if(myid.eq.0.and.ldebug) print *,myid,'building direct o_nnz for mat with',C%dof,'dof'
        call DMCreateGlobalVector(C%da,v_o_nnz,ierr) ;CHKERRQ(ierr)
        call DMCreateGlobalVector(C%da,v_d_nnz,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C%da,v_o_nnz,xo,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C%da,v_d_nnz,xd,ierr) ;CHKERRQ(ierr)

        xo = i0
        xd = i1
        xd( i0    ,:,:,C%zs+1:C%ze   ) = C%dof+i1 ! Edir_vertical depends on 3 values Edir_vertical,xaxis,yaxis :: starting with second entries(seen from top)
        xd([i1,i2],:,:,C%zs   :C%ze-i1) = C%dof+i1 ! Edir_xaxis,yaxis depends on 3 values Edir_vertical,xaxis,yaxis :: starting with first entries(seen from top)

        do j=C%ys,C%ye
                 lsun_east  = (sun%xinc.eq.i0)

                 if( C%neighbors(14).ne.myid .and. C%neighbors(14).ge.i0 ) then ! neigh east
                         if( lsun_east ) then 
                                ! if the sun is in the east, the channels in the last box are influenced by the 2nd channel which is a ghost
                                xo(i0, C%xe, j, C%zs+1:C%ze) = xo(i0, C%xe, j, C%zs+1:C%ze)+i1 ! channel 1 from zs+1 to ze
                                xd(i0, C%xe, j, C%zs+1:C%ze) = xd(i0, C%xe, j, C%zs+1:C%ze)-i1

                                xo([i1,i2], C%xe, j, C%zs:C%ze-1) = xo([i1,i2], C%xe, j, C%zs:C%ze-1)+i1 ! channel 2 and 3 from zs
                                xd([i1,i2], C%xe, j, C%zs:C%ze-1) = xd([i1,i2], C%xe, j, C%zs:C%ze-1)-i1
                        endif
                endif
        enddo
        do i=C%xs,C%xe
                 lsun_north = (sun%yinc.eq.i0 )

                 if( C%neighbors(16).ne.myid .and. C%neighbors(16).ge.i0 ) then ! neigh north
                         if( lsun_north ) then 
                                ! if the sun is in the north, the 3rd channel is a ghost
                                xo(i0, i, C%ye, C%zs+1:C%ze) = xo(i0, i, C%ye, C%zs+1:C%ze)+i1 ! channel 1 from zs+1 to ze
                                xd(i0, i, C%ye, C%zs+1:C%ze) = xd(i0, i, C%ye, C%zs+1:C%ze)-i1

                                xo([i1,i2],  i, C%ye, C%zs:C%ze-1) = xo([i1,i2], i, C%ye, C%zs:C%ze-1)+i1 ! channel 2 and 3 from zs
                                xd([i1,i2],  i, C%ye, C%zs:C%ze-1) = xd([i1,i2], i, C%ye, C%zs:C%ze-1)-i1
                        endif
                endif
        enddo
        do j=C%ys,C%ye
                 lsun_east  = (sun%xinc.eq.i0)

                 if( C%neighbors(12).ne.myid.and. C%neighbors(12).ge.i0 ) then ! neigh west
                         if( .not. lsun_east ) then 
                                ! if the sun is in the west, the 2nd channel is solemnly dependant on ghost values
                                xo(i1, C%xs, j, C%zs:C%ze-1) = i3
                                xd(i1, C%xs, j, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo
        do i=C%xs,C%xe
                 lsun_north = (sun%yinc.eq.i0 )

                 if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                         if( .not. lsun_north ) then 
                                ! if the sun is in the south, the 3rd channel is solemnly dependant on ghost values
                                xo(i2, i, C%ys, C%zs:C%ze-1) = i3
                                xd(i2, i, C%ys, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo

         do k=C%zs,C%ze-1
           do j=C%ys,C%ye
             do i=C%xs,C%xe      ! i,j,k indices are defined on petsc global grid
               if( atm%l1d(i,j,k) ) then
                 xo(:,i,j,k) = i0
                 xd(:,i,j,k) = i1
                 xd(0,i,j,k) = i2
               endif
             enddo
           enddo
         enddo

        call DMDAVecRestoreArrayF90(C%da,v_o_nnz,xo,ierr) ;CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C%da,v_d_nnz,xd,ierr) ;CHKERRQ(ierr)
        
        call VecGetLocalSize(v_o_nnz,vsize,ierr) ;CHKERRQ(ierr)
        allocate(o_nnz(0:vsize-1))
        allocate(d_nnz(0:vsize-1))

        call VecGetArrayF90(v_o_nnz,xx_v_o,ierr) ;CHKERRQ(ierr)
        call VecGetArrayF90(v_d_nnz,xx_v_d,ierr) ;CHKERRQ(ierr)
        o_nnz=int(xx_v_o)
        d_nnz=int(xx_v_d)
        call VecRestoreArrayF90(v_o_nnz,xx_v_o,ierr) ;CHKERRQ(ierr)
        call VecRestoreArrayF90(v_d_nnz,xx_v_d,ierr) ;CHKERRQ(ierr)

        call VecDestroy(v_o_nnz,ierr) ;CHKERRQ(ierr)
        call VecDestroy(v_d_nnz,ierr) ;CHKERRQ(ierr)

        if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
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

        call MatZeroEntries(A, ierr) ;CHKERRQ(ierr)
        call mat_set_diagonal(A,C)

        do k=C%zs,C%ze-1
          do j=C%ys,C%ye
            do i=C%xs,C%xe        

              if( atm%l1d(i,j,k) ) then
                call set_eddington_coeff(A, i,j,k)
              else
                call set_tenstream_coeff(C, A, i,j,k,ierr); CHKERRQ(ierr)
              endif

            enddo 
          enddo 
        enddo

        if(myid.eq.0.and.ldebug) print *,myid,'setup_direct_matrix done'

        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)  

        call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

        contains 
          subroutine set_tenstream_coeff(C,A,i,j,k,ierr)
              type(t_coord),intent(in) :: C
              Mat,intent(inout) :: A
              integer(iintegers),intent(in) :: i,j,k
              PetscErrorCode, intent(out) :: ierr

              MatStencil :: row(4,C%dof)  ,col(4,C%dof)
              PetscReal :: v(C%dof**2),coeffs(C%dof**2),norm

              PetscInt,parameter :: entries(8)=[0,8,16,24,32,40,48,56]

              integer(iintegers) :: dst,src

              ierr=0

              dst = 1 ; row(MatStencil_i,dst) = i            ; row(MatStencil_j,dst) = j            ; row(MatStencil_k,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
              dst = 2 ; row(MatStencil_i,dst) = i            ; row(MatStencil_j,dst) = j            ; row(MatStencil_k,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
              dst = 3 ; row(MatStencil_i,dst) = i            ; row(MatStencil_j,dst) = j            ; row(MatStencil_k,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
              dst = 4 ; row(MatStencil_i,dst) = i            ; row(MatStencil_j,dst) = j            ; row(MatStencil_k,dst) = k+1 ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
              dst = 5 ; row(MatStencil_i,dst) = i+sun%xinc   ; row(MatStencil_j,dst) = j            ; row(MatStencil_k,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
              dst = 6 ; row(MatStencil_i,dst) = i+sun%xinc   ; row(MatStencil_j,dst) = j            ; row(MatStencil_k,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
              dst = 7 ; row(MatStencil_i,dst) = i            ; row(MatStencil_j,dst) = j+sun%yinc   ; row(MatStencil_k,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
              dst = 8 ; row(MatStencil_i,dst) = i            ; row(MatStencil_j,dst) = j+sun%yinc   ; row(MatStencil_k,dst) = k   ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid

              src = 1 ; col(MatStencil_i,src) = i            ; col(MatStencil_j,src) = j            ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
              src = 2 ; col(MatStencil_i,src) = i            ; col(MatStencil_j,src) = j            ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
              src = 3 ; col(MatStencil_i,src) = i            ; col(MatStencil_j,src) = j            ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
              src = 4 ; col(MatStencil_i,src) = i            ; col(MatStencil_j,src) = j            ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
              src = 5 ; col(MatStencil_i,src) = i+1-sun%xinc ; col(MatStencil_j,src) = j            ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
              src = 6 ; col(MatStencil_i,src) = i+1-sun%xinc ; col(MatStencil_j,src) = j            ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
              src = 7 ; col(MatStencil_i,src) = i            ; col(MatStencil_j,src) = j+1-sun%yinc ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:
              src = 8 ; col(MatStencil_i,src) = i            ; col(MatStencil_j,src) = j+1-sun%yinc ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:

              call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.True., coeffs, atm%l1d(i,j,k), [sun%symmetry_phi, sun%theta])

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

              do src=1,C%dof
                v(entries+src) = coeffs( i1+(src-i1)*C%dof : src*C%dof )
              enddo
              call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)

          end subroutine

          subroutine set_eddington_coeff(A,i,j,k)
              Mat,intent(inout) :: A
              integer(iintegers),intent(in) :: i,j,k

              MatStencil :: row(4,1)  ,col(4,1)
              PetscReal :: v(1)
              integer(iintegers) :: src

              if(luse_eddington) then
                v = atm%a33(i,j,k)
              else
                call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.True., v, atm%l1d(i,j,k), [sun%symmetry_phi, sun%theta] )
              endif

              col(MatStencil_i,i1) = i      ; col(MatStencil_j,i1) = j       ; col(MatStencil_k,i1) = k    
              row(MatStencil_i,i1) = i      ; row(MatStencil_j,i1) = j       ; row(MatStencil_k,i1) = k+1  

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

          PetscScalar,pointer,dimension(:,:,:,:) :: xincSolar

          PetscReal :: Az
          Az = atm%dx*atm%dy

          call VecSet(incSolar,zero,ierr) ;CHKERRQ(ierr)

          call DMDAVecGetArrayF90(C_dir%da,incSolar,xincSolar,ierr) ;CHKERRQ(ierr)

          xincSolar(i0:i3,:,:,C_dir%zs) = edirTOA* Az * .25_ireals * sun%costheta

          call DMDAVecRestoreArrayF90(C_dir%da,incSolar,xincSolar,ierr) ;CHKERRQ(ierr)

          if(myid.eq.0 .and. ldebug) print *,myid,'Setup of IncSolar done',edirTOA* sun%costheta

        end subroutine

        subroutine set_diff_coeff(A,C)
            Mat :: A
            type(t_coord) :: C

            PetscInt :: i,j,k

            call PetscLogStagePush(logstage(4),ierr) ;CHKERRQ(ierr)

            if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients for diffuse Light'

!            call MatZeroEntries(A, ierr) ;CHKERRQ(ierr)
!            call mat_set_diagonal(A,C)

            do k=C%zs,C%ze-1
              do j=C%ys,C%ye
                do i=C%xs,C%xe

                  if( atm%l1d(i,j,k) ) then
                    call set_eddington_coeff(A, i,j,k)
                  else
                    call set_tenstream_coeff(C, A, i,j,k, ierr); CHKERRQ(ierr)
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
            subroutine set_tenstream_coeff(C,A,i,j,k,ierr)
                type(t_coord),intent(in) :: C
                Mat,intent(inout) :: A
                integer(iintegers),intent(in) :: i,j,k
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

                src = 0; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_dn
                src = 1; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k+1   ; col(MatStencil_c,src) = E_up 
                src = 2; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_ri_m
                src = 3; col(MatStencil_i,src) = i+1  ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_le_m
                src = 4; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_ri_p
                src = 5; col(MatStencil_i,src) = i+1  ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_le_p
                src = 6; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_fw_m
                src = 7; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j+1   ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_ba_m
                src = 8; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_fw_p
                src = 9; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j+1   ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_ba_p

                !           Destinations
                !         ________E_up___
                !        |               |
                !    E_le|               |E_ri_p
                !        |               |  
                !        |               |
                !    E_le|               |E_ri_m
                !        |_______________|  
                !                 E_dn

                dst = 0; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_up  
                dst = 1; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k+1   ; row(MatStencil_c,dst) = E_dn  
                dst = 2; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_le_m
                dst = 3; row(MatStencil_i,dst) = i+1  ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_ri_m
                dst = 4; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_le_p
                dst = 5; row(MatStencil_i,dst) = i+1  ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_ri_p
                dst = 6; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_ba_m
                dst = 7; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j+1   ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_fw_m
                dst = 8; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_ba_p
                dst = 9; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j+1   ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_fw_p

                call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.False., coeffs, atm%l1d(i,j,k))

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

            subroutine set_eddington_coeff(A,i,j,k)
                Mat,intent(inout) :: A
                integer(iintegers),intent(in) :: i,j,k

                MatStencil :: row(4,2)  ,col(4,2)
                PetscReal :: v(4),twostr_coeff(4) ! v ==> a12,a11,a11,a12
                integer(iintegers) :: src,dst

                if(luse_eddington ) then
                  v = [ atm%a12(i,j,k), atm%a11(i,j,k), atm%a11(i,j,k), atm%a12(i,j,k)]
                else
                  call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.False., twostr_coeff, atm%l1d(i,j,k)) !twostr_coeff ==> a12,a11,a12,a11
                  v = [ twostr_coeff(1), twostr_coeff(2) , twostr_coeff(2) , twostr_coeff(1) ]
                endif

                src = 1; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k     ; col(MatStencil_c,src) = E_dn
                src = 2; col(MatStencil_i,src) = i    ; col(MatStencil_j,src) = j     ; col(MatStencil_k,src) = k+1   ; col(MatStencil_c,src) = E_up 

                dst = 1; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k     ; row(MatStencil_c,dst) = E_up  
                dst = 2; row(MatStencil_i,dst) = i    ; row(MatStencil_j,dst) = j     ; row(MatStencil_k,dst) = k+1   ; row(MatStencil_c,dst) = E_dn  

                call MatSetValuesStencil(A,i2, row, i2, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)


            end subroutine
            subroutine set_albedo_coeff(C,A)
                type(t_coord) :: C
                Mat,intent(inout) :: A

                MatStencil :: row(4,1)  ,col(4,1)
                integer(iintegers) :: i,j

                ! Set surface albedo values
                col(MatStencil_k,i1) = C%ze
                col(MatStencil_c,i1) = E_dn

                row(MatStencil_k,i1) = C%ze
                row(MatStencil_c,i1) = E_up

                do j=C%ys,C%ye
                  row(MatStencil_j,i1) = j 
                  col(MatStencil_j,i1) = j 

                  do i=C%xs,C%xe
                    row(MatStencil_i,i1) = i          
                    col(MatStencil_i,i1) = i        

                    call MatSetValuesStencil(A,i1, row, i1, col , -atm%albedo ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)

                  enddo
                enddo

            end subroutine
        end subroutine
 
subroutine setup_b(edir,b)
        Vec :: edir
        Vec :: local_b,b

        PetscScalar,pointer,dimension(:,:,:,:) :: xsrc,xedir
        PetscReal :: diff2diff(C_diff%dof**2), dir2diff(C_dir%dof*C_diff%dof),twostr_coeff(2)
        PetscInt :: i,j,k,src

        call PetscLogStagePush(logstage(6),ierr) ;CHKERRQ(ierr)
        
        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly...'

        call DMCreateLocalVector(C_diff%da,local_b,ierr) ;CHKERRQ(ierr)
        call VecSet(local_b,zero,ierr) ;CHKERRQ(ierr)

        call DMDAVecGetArrayF90(C_diff%da,local_b,xsrc,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C_dir%da,edir,xedir,ierr) ;CHKERRQ(ierr)

        call set_solar_source()
        if(allocated(atm%planck) ) call set_thermal_source()

        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients ...done'

        call DMDAVecRestoreArrayF90(C_dir%da,edir,xedir,ierr) ;CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C_diff%da,local_b,xsrc,ierr) ;CHKERRQ(ierr)

        call VecSet(b,zero,ierr) ;CHKERRQ(ierr) ! reset global Vec

        call DMLocalToGlobalBegin(C_diff%da,local_b,ADD_VALUES, b,ierr) ;CHKERRQ(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
        call DMLocalToGlobalEnd  (C_diff%da,local_b,ADD_VALUES, b,ierr) ;CHKERRQ(ierr)

        call VecDestroy(local_b,ierr) ;CHKERRQ(ierr)

        call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly done'
      contains
        subroutine set_thermal_source()
            PetscReal :: Ax,Ay,Az,c1,c2,c3,b0,b1,dtau

            if(myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector ... setting thermal source terms'
            Az = atm%dx*atm%dy

            do k=C_diff%zs,C_diff%ze-1 
              do j=C_diff%ys,C_diff%ye         
                do i=C_diff%xs,C_diff%xe    

                  if( atm%l1d(i,j,k) ) then

                    if(luse_eddington ) then
                      !see libradtran rodents
                      dtau = atm%dz(i,j,k) * ( atm%op(i,j,k)%kabs + atm%op(i,j,k)%ksca)
                      if( dtau.gt.0.01_ireals ) then
                        b0 = atm%planck(i,j,k)
                        b1 = ( atm%planck(i,j,k)-atm%planck(i,j,k+1) ) / dtau
                      else
                        b0 = .5_ireals*(atm%planck(i,j,k)+atm%planck(i,j,k+1))
                        b1 = zero
                      endif
                      c1 = atm%g1(i,j,k) * (b0 + b1*dtau)
                      c2 = atm%g2(i,j,k) * b1
                      c3 = atm%g1(i,j,k) * b0
                    else
                      stop 'boxmc coeff for thermal emission not supported at the moment'
                    endif
                    xsrc(E_up   ,i,j,k)   = xsrc(E_up   ,i,j,k)   + ( - atm%a11(i,j,k)*(c1+c2) - atm%a12(i,j,k)*(c3-c2) + c2 + c3 )*Az*pi
                    xsrc(E_dn   ,i,j,k+1) = xsrc(E_dn   ,i,j,k+1) + ( - atm%a12(i,j,k)*(c1+c2) - atm%a11(i,j,k)*(c3-c2) + c1 - c2 )*Az*pi

                  else ! Tenstream source terms
                    Ax = atm%dy*atm%dz(i,j,k)
                    Ay = atm%dx*atm%dz(i,j,k)

                    call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.False., diff2diff, atm%l1d(i,j,k) )
                    b0 = .5_ireals*(atm%planck(i,j,k)+atm%planck(i,j,k+i1)) *pi
                    xsrc(E_up   ,i,j,k)   = xsrc(E_up   ,i,j,k)   +  b0  *(one-sum( diff2diff( E_up  *C_diff%dof+i1 : E_up  *C_diff%dof+C_diff%dof ))) *Az
                    xsrc(E_dn   ,i,j,k+1) = xsrc(E_dn   ,i,j,k+1) +  b0  *(one-sum( diff2diff( E_dn  *C_diff%dof+i1 : E_dn  *C_diff%dof+C_diff%dof ))) *Az
                    xsrc(E_le_m ,i,j,k)   = xsrc(E_le_m ,i,j,k)   +  b0  *(one-sum( diff2diff( E_le_m*C_diff%dof+i1 : E_le_m*C_diff%dof+C_diff%dof ))) *Ax*.5_ireals
                    xsrc(E_le_p ,i,j,k)   = xsrc(E_le_p ,i,j,k)   +  b0  *(one-sum( diff2diff( E_le_p*C_diff%dof+i1 : E_le_p*C_diff%dof+C_diff%dof ))) *Ax*.5_ireals
                    xsrc(E_ri_m ,i+1,j,k) = xsrc(E_ri_m ,i+1,j,k) +  b0  *(one-sum( diff2diff( E_ri_m*C_diff%dof+i1 : E_ri_m*C_diff%dof+C_diff%dof ))) *Ax*.5_ireals
                    xsrc(E_ri_p ,i+1,j,k) = xsrc(E_ri_p ,i+1,j,k) +  b0  *(one-sum( diff2diff( E_ri_p*C_diff%dof+i1 : E_ri_p*C_diff%dof+C_diff%dof ))) *Ax*.5_ireals
                    xsrc(E_ba_m ,i,j,k)   = xsrc(E_ba_m ,i,j,k)   +  b0  *(one-sum( diff2diff( E_ba_m*C_diff%dof+i1 : E_ba_m*C_diff%dof+C_diff%dof ))) *Ay*.5_ireals
                    xsrc(E_ba_p ,i,j,k)   = xsrc(E_ba_p ,i,j,k)   +  b0  *(one-sum( diff2diff( E_ba_p*C_diff%dof+i1 : E_ba_p*C_diff%dof+C_diff%dof ))) *Ay*.5_ireals
                    xsrc(E_fw_m ,i,j+1,k) = xsrc(E_fw_m ,i,j+1,k) +  b0  *(one-sum( diff2diff( E_fw_m*C_diff%dof+i1 : E_fw_m*C_diff%dof+C_diff%dof ))) *Ay*.5_ireals
                    xsrc(E_fw_p ,i,j+1,k) = xsrc(E_fw_p ,i,j+1,k) +  b0  *(one-sum( diff2diff( E_fw_p*C_diff%dof+i1 : E_fw_p*C_diff%dof+C_diff%dof ))) *Ay*.5_ireals
                  endif ! 1D or Tenstream?

                enddo
              enddo

            enddo ! k

            ! Thermal emission at surface
            k = C_diff%ze
            do j=C_diff%ys,C_diff%ye         
              do i=C_diff%xs,C_diff%xe    
                xsrc(E_up   ,i,j,k) = xsrc(E_up   ,i,j,k) + atm%planck(i,j,k)*Az *(one-atm%albedo)*pi
              enddo
            enddo
        end subroutine
        subroutine set_solar_source()
            if(myid.eq.0.and.ldebug) print *,'Assembly of SRC-Vector .. setting solar source',sum(xedir(0:3,:,:,:))/size(xedir(0:3,:,:,:))
            do k=C_diff%zs,C_diff%ze-1 
              do j=C_diff%ys,C_diff%ye         
                do i=C_diff%xs,C_diff%xe    

                  if( any (xedir(:,i,j,k) .gt. epsilon(one)) ) then
                    if( atm%l1d(i,j,k) ) then
                      dir2diff = zero
                      if(luse_eddington ) then
                        ! Only transport the 4 tiles from dir0 to the Eup and Edn
                        do src=1,4
                          dir2diff(E_up  +i1+(src-1)*C_diff%dof) = atm%a13(i,j,k)
                          dir2diff(E_dn  +i1+(src-1)*C_diff%dof) = atm%a23(i,j,k)
                        enddo

                      else
                        call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.False., twostr_coeff, atm%l1d(i,j,k), [sun%symmetry_phi, sun%theta])
                        do src=1,4
                          dir2diff(E_up  +i1+(src-1)*C_diff%dof) = twostr_coeff(1)
                          dir2diff(E_dn  +i1+(src-1)*C_diff%dof) = twostr_coeff(2)
                        enddo
                      endif

                      do src=1,C_dir%dof-4
                        xsrc(E_up   ,i,j,k)   = xsrc(E_up   ,i,j,k)   +  xedir(src-1,i,j,k)*dir2diff(E_up  +i1+(src-1)*C_diff%dof)
                        xsrc(E_dn   ,i,j,k+1) = xsrc(E_dn   ,i,j,k+1) +  xedir(src-1,i,j,k)*dir2diff(E_dn  +i1+(src-1)*C_diff%dof)
                      enddo

                    else ! Tenstream source terms

                      call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.False., dir2diff,  atm%l1d(i,j,k), [sun%symmetry_phi, sun%theta] )

                      do src=1,C_dir%dof
                        xsrc(E_up   ,i,j,k)   = xsrc(E_up   ,i,j,k)   +  xedir(src-1,i,j,k)*dir2diff(E_up  +i1+(src-1)*C_diff%dof) 
                        xsrc(E_dn   ,i,j,k+1) = xsrc(E_dn   ,i,j,k+1) +  xedir(src-1,i,j,k)*dir2diff(E_dn  +i1+(src-1)*C_diff%dof) 
                        xsrc(E_le_m ,i,j,k)   = xsrc(E_le_m ,i,j,k)   +  xedir(src-1,i,j,k)*dir2diff(E_le_m+i1+(src-1)*C_diff%dof) 
                        xsrc(E_le_p ,i,j,k)   = xsrc(E_le_p ,i,j,k)   +  xedir(src-1,i,j,k)*dir2diff(E_le_p+i1+(src-1)*C_diff%dof) 
                        xsrc(E_ri_m ,i+1,j,k) = xsrc(E_ri_m ,i+1,j,k) +  xedir(src-1,i,j,k)*dir2diff(E_ri_m+i1+(src-1)*C_diff%dof) 
                        xsrc(E_ri_p ,i+1,j,k) = xsrc(E_ri_p ,i+1,j,k) +  xedir(src-1,i,j,k)*dir2diff(E_ri_p+i1+(src-1)*C_diff%dof) 
                        xsrc(E_ba_m ,i,j,k)   = xsrc(E_ba_m ,i,j,k)   +  xedir(src-1,i,j,k)*dir2diff(E_ba_m+i1+(src-1)*C_diff%dof) 
                        xsrc(E_ba_p ,i,j,k)   = xsrc(E_ba_p ,i,j,k)   +  xedir(src-1,i,j,k)*dir2diff(E_ba_p+i1+(src-1)*C_diff%dof) 
                        xsrc(E_fw_m ,i,j+1,k) = xsrc(E_fw_m ,i,j+1,k) +  xedir(src-1,i,j,k)*dir2diff(E_fw_m+i1+(src-1)*C_diff%dof) 
                        xsrc(E_fw_p ,i,j+1,k) = xsrc(E_fw_p ,i,j+1,k) +  xedir(src-1,i,j,k)*dir2diff(E_fw_p+i1+(src-1)*C_diff%dof) 
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
                xsrc(E_up   ,i,j,k) = sum(xedir(i0:i3,i,j,k))*atm%albedo
              enddo
            enddo
        end subroutine
    end subroutine

subroutine calc_flx_div(edir,ediff,abso)
        Vec :: edir,ediff,abso
        PetscReal,pointer,dimension(:,:,:,:) :: xediff,xedir,xabso
        PetscInt :: i,j,k!d,li,lj,lk
        Vec :: ledir,lediff ! local copies of vectors, including ghosts
        PetscReal :: div2(13)
        PetscReal :: Volume

        !        real(ireals) :: c_dir2dir(C_dir%dof**2)
        !        real(ireals) :: c_dir2diff(C_dir%dof*C_diff%dof)
        !        real(ireals) :: c_diff2diff(C_diff%dof**2)
        !        integer(iintegers) :: isrc

        if(myid.eq.0.and.ldebug) print *,'Calculating flux divergence'
        call VecSet(abso,zero,ierr) ;CHKERRQ(ierr)

        ! Copy ghosted values for direct vec
        call DMCreateLocalVector(C_dir%da ,ledir ,ierr)                   ; CHKERRQ(ierr)
        call VecSet(ledir ,zero,ierr)                                     ; CHKERRQ(ierr)
        call DMGlobalToLocalBegin(C_dir%da ,edir ,ADD_VALUES,ledir ,ierr) ; CHKERRQ(ierr)
        call DMGlobalToLocalEnd(C_dir%da ,edir ,ADD_VALUES,ledir ,ierr)   ; CHKERRQ(ierr)

        ! Copy ghosted values for diffuse vec
        call DMCreateLocalVector(C_diff%da,lediff,ierr)                   ; CHKERRQ(ierr)
        call VecSet(lediff,zero,ierr)                                     ; CHKERRQ(ierr)
        call DMGlobalToLocalBegin(C_diff%da,ediff,ADD_VALUES,lediff,ierr) ; CHKERRQ(ierr)
        call DMGlobalToLocalEnd(C_diff%da,ediff,ADD_VALUES,lediff,ierr)   ; CHKERRQ(ierr)

        ! calculate absorption by flux divergence
        call DMDAVecGetArrayF90(C_diff%da,lediff,xediff,ierr) ; CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C_dir%da ,ledir ,xedir ,ierr) ; CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C_one%da ,abso ,xabso ,ierr)  ; CHKERRQ(ierr)

        do k=C_one%zs,C_one%ze
          do j=C_one%ys,C_one%ye         
            do i=C_one%xs,C_one%xe      

              div2 = zero
              Volume = atm%dx * atm%dy * atm%dz(i,j,k)
              if(atm%l1d(i,j,k)) then ! one dimensional i.e. twostream
                ! Divergence    =                       Incoming                -       Outgoing
                div2( 1) = sum( xedir(i0:i3, i, j , k)  - xedir(i0:i3 , i, j, k+i1  ) )

                div2( 4) = ( xediff(E_up  ,i  ,j  ,k+1)  - xediff(E_up  ,i  ,j  ,k  )  )
                div2( 5) = ( xediff(E_dn  ,i  ,j  ,k  )  - xediff(E_dn  ,i  ,j  ,k+1)  )

                !                div2(1) = sum( xedir(i0:i3, i, j,k ) ) * (one-atm%a33(li,lj,lk)-atm%a13(li,lj,lk)-atm%a23(li,lj,lk))
                !                if((one-atm%a33(li,lj,lk)-atm%a13(li,lj,lk)-atm%a23(li,lj,lk)).lt.zero) then 
                !                  ! negative values are of course physically not meaningful, they do however happen because of inaccuracies in delta eddington twostream coefficients
                !                  ! we can try to estimate absorption by ...
                !                  div2(1) = sum( xedir(i0:i3, i, j,k ) ) * exp(- atm%op(li,lj,lk)%kabs*atm%dz(lk)/sun%costheta ) ! absorption of direct radiation
                !                  ! this neglects absorption of diffuse radiation! this is
                !                  ! probably bad!...
                !                endif

                !                if(xabso(i0,i,j,k).lt.zero) then
                !                  print *,'OhOh :: 1-D Attention: neg. abso',xabso(i0,i,j,k),' at ',i,j,k
                !                  print *,'div2',div2
                !                endif

              else

                !          Divergence     =                        Incoming                        -                   Outgoing

                div2( 1) = sum( xedir(i0:i3 , i             , j             , k)  - xedir(i0:i3 , i          , j          , k+i1  ) )
                div2( 2) = sum( xedir(i4:i5 , i+i1-sun%xinc , j             , k)  - xedir(i4:i5 , i+sun%xinc , j          , k) )
                div2( 3) = sum( xedir(i6:i7 , i             , j+i1-sun%yinc , k)  - xedir(i6:i7 , i          , j+sun%yinc , k) )

                div2( 4) = ( xediff(E_up  ,i  ,j  ,k+1)  - xediff(E_up  ,i  ,j  ,k  )  )
                div2( 5) = ( xediff(E_dn  ,i  ,j  ,k  )  - xediff(E_dn  ,i  ,j  ,k+1)  )
                div2( 6) = ( xediff(E_le_m,i+1,j  ,k  )  - xediff(E_le_m,i  ,j  ,k  )  )
                div2( 7) = ( xediff(E_le_p,i+1,j  ,k  )  - xediff(E_le_p,i  ,j  ,k  )  )
                div2( 8) = ( xediff(E_ri_m,i  ,j  ,k  )  - xediff(E_ri_m,i+1,j  ,k  )  )
                div2( 9) = ( xediff(E_ri_p,i  ,j  ,k  )  - xediff(E_ri_p,i+1,j  ,k  )  )
                div2(10) = ( xediff(E_ba_m,i  ,j+1,k  )  - xediff(E_ba_m,i  ,j  ,k  )  )
                div2(11) = ( xediff(E_ba_p,i  ,j+1,k  )  - xediff(E_ba_p,i  ,j  ,k  )  )
                div2(12) = ( xediff(E_fw_m,i  ,j  ,k  )  - xediff(E_fw_m,i  ,j+1,k  )  )
                div2(13) = ( xediff(E_fw_p,i  ,j  ,k  )  - xediff(E_fw_p,i  ,j+1,k  )  )

                !         Divergence can also be expressed as the sum of divergence for each stream -- probably more stable... TODO: not tested at all
                !                call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.True., c_dir2dir,   atm%l1d(lk), [sun%symmetry_phi, sun%theta])
                !                call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.False.,c_dir2diff,  atm%l1d(lk), [sun%symmetry_phi, sun%theta])
                !                call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.False.,c_diff2diff, atm%l1d(lk) )
                !
                !                do isrc=1,4
                !                  div2(1) = div2(1) + xedir(isrc-i1,i,j,k) * max(zero, one &                          ! Absorption of direct streams is One minus
                !                                    - sum(c_dir2dir ( (isrc-1) *C_dir%dof  +1 : isrc*C_dir%dof  ) ) & ! - transmission
                !                                    - sum(c_dir2diff( (isrc-1) *C_diff%dof +1 : isrc*C_diff%dof ) ) ) ! - diffuse sources
                !                enddo     
                !                do isrc=5,6
                !                  div2(2) = div2(2) + xedir(isrc-i1,i+1-sun%xinc,j,k) * max(zero, one &
                !                                    - sum(c_dir2dir ( (isrc-1) *C_dir%dof  +1 : isrc*C_dir%dof  ) ) &
                !                                    - sum(c_dir2diff( (isrc-1) *C_diff%dof +1 : isrc*C_diff%dof ) ) )
                !                enddo
                !                do isrc=7,8
                !                  div2(3) = div2(3) + xedir(isrc-i1,i,j+1-sun%yinc,k) * max(zero, one &
                !                                    - sum(c_dir2dir ( (isrc-1) *C_dir%dof  +1 : isrc*C_dir%dof  ) ) &
                !                                    - sum(c_dir2diff( (isrc-1) *C_diff%dof +1 : isrc*C_diff%dof ) ) )
                !                enddo
                !
                !                div2( 4) = xediff(E_up  ,i  ,j  ,k+1)* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2( 5) = xediff(E_dn  ,i  ,j  ,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2( 6) = xediff(E_le_m,i+1,j  ,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2( 7) = xediff(E_le_p,i+1,j  ,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2( 8) = xediff(E_ri_m,i  ,j  ,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2( 9) = xediff(E_ri_p,i  ,j  ,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2(10) = xediff(E_ba_m,i  ,j+1,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2(11) = xediff(E_ba_p,i  ,j+1,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2(12) = xediff(E_fw_m,i  ,j  ,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )
                !                div2(13) = xediff(E_fw_p,i  ,j  ,k  )* max(zero, one - sum( c_diff2diff( E_up*C_diff%dof +1 : E_up*C_diff%dof +C_diff%dof  )  ) )

                !                if(any(div2.lt.epsilon(xabso)*100._ireals) ) then
                !                  print *,'Attention: neg. abso at ',i,j,k
                !                  print *,'div2',div2
                !                endif
              endif
              xabso(i0,i,j,k) = sum(div2) / Volume
!              if(myid.eq.0) print *,'flxdiv',i,j,k,'::',xabso(i0,i,j,k),Volume,'::',div2 !TODO
            enddo                             
          enddo                             
        enddo   
        
        call DMDAVecRestoreArrayF90(C_one%da ,abso ,xabso ,ierr)  ; CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C_diff%da,lediff,xediff,ierr) ; CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C_dir%da ,ledir ,xedir ,ierr) ; CHKERRQ(ierr)
        
        call VecDestroy(lediff,ierr) ; CHKERRQ(ierr)
        call VecDestroy(ledir ,ierr) ; CHKERRQ(ierr)
end subroutine

subroutine solve(ksp,b,x)
      KSP :: ksp
      Vec:: b
      Vec:: x

      KSPConvergedReason :: reason
      PetscInt :: iter

      if(myid.eq.0.and.ldebug) print *,'Solving Matrix'

      call KSPSolve(ksp,b,x,ierr) ;CHKERRQ(ierr)
      call KSPGetIterationNumber(ksp,iter,ierr) ;CHKERRQ(ierr)
      call KSPGetConvergedReason(ksp,reason,ierr) ;CHKERRQ(ierr)

!      print *,'Source Vector:'
!      call VecView(b,PETSC_VIEWER_STDOUT_WORLD ,ierr) ;CHKERRQ(ierr)
!
!      print *,'Solution Vector:'
!      call VecView(x,PETSC_VIEWER_STDOUT_WORLD ,ierr) ;CHKERRQ(ierr)

      if(myid.eq.0.and.ldebug) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      if(reason.eq.KSP_DIVERGED_ITS) then
        if(myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
        return
      endif

      if(reason.le.0) then
              if(myid.eq.0.and.ldebug) print *,myid,'Resetted initial guess to zero and try again with gmres:'
              call VecSet(x,zero,ierr) ;CHKERRQ(ierr)
              call KSPSetType(ksp,KSPGMRES,ierr) ;CHKERRQ(ierr)
              call KSPSetUp(ksp,ierr) ;CHKERRQ(ierr)
              call KSPSolve(ksp,b,x,ierr) ;CHKERRQ(ierr)
              call KSPGetIterationNumber(ksp,iter,ierr) ;CHKERRQ(ierr)
              call KSPGetConvergedReason(ksp,reason,ierr) ;CHKERRQ(ierr)

              call KSPSetFromOptions(ksp,ierr) ;CHKERRQ(ierr)
              call KSPSetUp(ksp,ierr) ;CHKERRQ(ierr)
              if(myid.eq.0.and.ldebug) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      endif

      if(reason.le.0) then
              print *,'*********************************************************** SOLVER did NOT converge :( ***************************************************88'
              call exit()
      endif
end subroutine

subroutine setup_ksp(ksp,C,A,linit, prefix)
      KSP :: ksp
      type(t_coord) :: C
      Mat :: A
      PC  :: prec
      logical :: linit

      MatNullSpace :: nullspace
!      PetscReal,parameter :: rtol=epsilon(one)*10, atol=epsilon(one)*10
      character(len=*),optional :: prefix

      PetscReal,parameter :: rtol=1e-5, atol=1e-5
      PetscInt,parameter  :: maxiter=100

      if(linit) return
      call PetscLogStagePush(logstage(8),ierr) ;CHKERRQ(ierr)
      if(myid.eq.0.and.ldebug) print *,'Setup KSP'

      call KSPCreate(imp_comm,ksp,ierr) ;CHKERRQ(ierr)
      if(present(prefix) ) call KSPAppendOptionsPrefix(ksp,trim(prefix),ierr) ;CHKERRQ(ierr)

      call KSPSetType(ksp,KSPFBCGS,ierr)  ;CHKERRQ(ierr)
      call KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr) ;CHKERRQ(ierr) 
      call KSPGetPC  (ksp,prec,ierr)  ;CHKERRQ(ierr)
      if(numnodes.eq.0) then
        call PCSetType (prec,PCILU,ierr);CHKERRQ(ierr)
      else
        call PCSetType (prec,PCBJACOBI,ierr);CHKERRQ(ierr)
      endif

      call KSPSetTolerances(ksp,rtol,atol*(C%dof*C%glob_xm*C%glob_ym*C%glob_zm) * count(.not.atm%l1d)/(one*size(atm%l1d)) ,PETSC_DEFAULT_REAL,maxiter,ierr);CHKERRQ(ierr)

      call KSPSetConvergenceTest(ksp,MyKSPConverged, PETSC_NULL_OBJECT,PETSC_NULL_FUNCTION,ierr)

      call KSPSetFromOptions(ksp,ierr) ;CHKERRQ(ierr)

      call KSPSetOperators(ksp,A,A,ierr) ;CHKERRQ(ierr)
      call KSPSetDM(ksp,C%da,ierr) ;CHKERRQ(ierr)
      call KSPSetDMActive(ksp,PETSC_FALSE,ierr) ;CHKERRQ(ierr)

      call KSPSetUp(ksp,ierr) ;CHKERRQ(ierr)

      call MatNullSpaceCreate( imp_comm, PETSC_TRUE, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, nullspace, ierr) ; CHKERRQ(ierr)
      call KSPSetNullspace(ksp, nullspace, ierr) ; CHKERRQ(ierr)

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
      initial_rnorm=rnorm
      return
    endif

    if(rnorm/initial_rnorm.le.rtol) then
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
    call PetscLogStageRegister('write_hdf5'      , logstage(10)     , ierr) ;CHKERRQ(ierr)

    if(myid.eq.0) print *, 'Logging stages' , logstage
    logstage_init=.True.
end subroutine

subroutine twostream(edirTOA)
    real(ireals),intent(in) :: edirTOA

    PetscReal,pointer,dimension(:,:,:,:) :: xv_dir,xv_diff
    integer(iintegers) :: i,j,src

    real(ireals),allocatable :: dtau(:),kext(:),w0(:),g(:),S(:),Edn(:),Eup(:)
    real(ireals) :: mu0,Az,incSolar

    call PetscLogStagePush(logstage(8),ierr) ;CHKERRQ(ierr)
    Az = atm%dx*atm%dy

    allocate( dtau(C_dir%zm-1) )
    allocate( kext(C_dir%zm-1) )
    allocate(   w0(C_dir%zm-1) )
    allocate(    g(C_dir%zm-1) )

    if(sun%theta.le.zero) then
      mu0=zero
      incSolar=zero
    else
      mu0 = sun%costheta
      incSolar = edirTOA* Az * sun%costheta
    endif

    call VecSet(edir_twostr ,zero,ierr); CHKERRQ(ierr)
    call VecSet(ediff_twostr,zero,ierr); CHKERRQ(ierr)

    call DMDAVecGetArrayF90(C_dir%da ,edir_twostr ,xv_dir ,ierr) ;CHKERRQ(ierr)
    call DMDAVecGetArrayF90(C_diff%da,ediff_twostr,xv_diff,ierr) ;CHKERRQ(ierr)

    allocate( S(C_dir%zm ) )
    allocate( Eup(C_diff%zm) )
    allocate( Edn(C_diff%zm) )

    if(myid.eq.0) print *,' CALCULATING DELTA EDDINGTON TWOSTREAM',incSolar/Az

    do j=C_dir%ys,C_dir%ye         
        do i=C_dir%xs,C_dir%xe

          kext = atm%op(i,j,:)%kabs + atm%op(i,j,:)%ksca
          dtau = atm%dz(i,j,:)* kext
          w0   = atm%op(i,j,:)%ksca / kext
          g    = atm%op(i,j,:)%g

          if(allocated(atm%planck) ) then
            call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo, S,Edn,Eup, planck=atm%planck(i,j,:)*Az )
          else
            call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo, S,Edn,Eup )
          endif

          do src=i0,i3
            xv_dir(src,i,j,:) = .25_ireals * S(:)
          enddo

          xv_diff(E_up,i,j,:) = Eup(:) 
          xv_diff(E_dn,i,j,:) = Edn(:) 
        enddo
    enddo

    call DMDAVecRestoreArrayF90(C_dir%da  ,edir_twostr ,xv_dir  ,ierr) ;CHKERRQ(ierr)
    call DMDAVecRestoreArrayF90(C_diff%da ,ediff_twostr,xv_diff ,ierr) ;CHKERRQ(ierr)

    deallocate(S)
    deallocate(Edn)
    deallocate(Eup)

    call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
end subroutine

subroutine scale_flx(v,C)
        Vec :: v
        type(t_coord) :: C
        PetscReal,pointer,dimension(:,:,:,:) :: xv
        PetscInt :: i,j,k!d,li,lj,lk
        PetscReal :: Ax,Ay,Az

        if(myid.eq.0.and.ldebug) print *,'rescaling fluxes'
        call DMDAVecGetArrayF90(C%da,v,xv,ierr) ;CHKERRQ(ierr)
        ! rescale total energy fluxes to average quantities i.e. W/m**2 or W/m**3

        Az = atm%dx*atm%dy

        do k=C%zs,C%ze-i1
          do j=C%ys,C%ye         
            do i=C%xs,C%xe      

              if(C%dof.eq.i8) then ! This is 8 stream direct radiation
                xv(i0:i3,i,j,k) = xv(i0:i3,i,j,k) / ( Az*.25 ) !*costheta0
              endif

              if(C%dof.eq.i10) then ! This is 10 stream diffuse radiation
                xv(E_up  ,i,j,k) = xv(E_up  ,i,j,k) / Az
                xv(E_dn  ,i,j,k) = xv(E_dn  ,i,j,k) / Az
              endif

              if(.not.atm%l1d(i,j,k)) then

                Ax = atm%dy*atm%dz(i,j,k)
                Ay = atm%dx*atm%dz(i,j,k)

                if(C%dof.eq.i8) then ! This is 8 stream direct radiation
                  xv(i4:i5,i,j,k) = xv(i4:i5,i,j,k) / ( Ax*.5  ) !*sintheta0
                  xv(i6:i7,i,j,k) = xv(i6:i7,i,j,k) / ( Ay*.5  ) !*sintheta0
                endif

                if(C%dof.eq.i10) then ! This is 10 stream diffuse radiation
                  xv(E_le_m,i,j,k) = xv(E_le_m,i,j,k) / Ax
                  xv(E_le_p,i,j,k) = xv(E_le_p,i,j,k) / Ax
                  xv(E_ri_m,i,j,k) = xv(E_ri_m,i,j,k) / Ax
                  xv(E_ri_p,i,j,k) = xv(E_ri_p,i,j,k) / Ax
                  xv(E_ba_m,i,j,k) = xv(E_ba_m,i,j,k) / Ay
                  xv(E_ba_p,i,j,k) = xv(E_ba_p,i,j,k) / Ay
                  xv(E_fw_m,i,j,k) = xv(E_fw_m,i,j,k) / Ay
                  xv(E_fw_p,i,j,k) = xv(E_fw_p,i,j,k) / Ay
                endif
              endif
            enddo
          enddo
        enddo

        k=C%ze
        do j=C%ys,C%ye         
          do i=C%xs,C%xe      

            if(C%dof.eq.i8) then ! This is 8 stream direct radiation
              xv (i0:i3 ,i,j,k) = xv (i0:i3 ,i,j,k) / ( Az*.25 ) !*costheta0
            endif
            if(C%dof.eq.i10) then ! This is 10 stream diffuse radiation
              xv(E_up  ,i,j,k) = xv(E_up  ,i,j,k) / Az
              xv(E_dn  ,i,j,k) = xv(E_dn  ,i,j,k) / Az
            endif
          enddo
        enddo

        call DMDAVecRestoreArrayF90(C%da ,v ,xv ,ierr) ;CHKERRQ(ierr)
end subroutine

    subroutine set_diff_initial_guess(inp,guess,C)
        ! Deprecated -- this is probably not helping convergence....
        Vec :: inp,guess
        type(t_coord) :: C

        Vec :: local_guess
        PetscScalar,pointer,dimension(:,:,:,:) :: xinp,xguess
        PetscReal :: diff2diff(C_diff%dof**2)!, dir2diff(C_dir%dof*C_diff%dof)
        PetscInt :: i,j,k,src

        if(myid.eq.0.and.ldebug) print *,'setting initial guess...'
        call DMCreateLocalVector(C%da,local_guess,ierr) ;CHKERRQ(ierr)
        call VecSet(local_guess,zero,ierr) ;CHKERRQ(ierr)

        call DMDAVecGetArrayF90(C%da,inp,  xinp,  ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C%da,local_guess,xguess,ierr) ;CHKERRQ(ierr)

        do k=C%zs,C%ze-1 
          do j=C%ys,C%ye         
            do i=C%xs,C%xe    
              if( .not. atm%l1d(i,j,k) ) then
                call get_coeff(atm%delta_op(i,j,k), atm%dz(i,j,k),.False., diff2diff, atm%l1d(i,j,k) )
                !                  print *,'xinp before',xinp(:,i,j,k)
                !                  print *,'xguess before',xguess(:,i,j,k)
                do src=1,C%dof
                  xguess(E_up   ,i,j,k)   = xguess(E_up   ,i,j,k)   +  xinp(src-1,i,j,k)*diff2diff(E_up  +i1+(src-1)*C%dof) 
                  xguess(E_dn   ,i,j,k+1) = xguess(E_dn   ,i,j,k+1) +  xinp(src-1,i,j,k)*diff2diff(E_dn  +i1+(src-1)*C%dof) 
                  xguess(E_le_m ,i,j,k)   = xguess(E_le_m ,i,j,k)   +  xinp(src-1,i,j,k)*diff2diff(E_le_m+i1+(src-1)*C%dof) 
                  xguess(E_le_p ,i,j,k)   = xguess(E_le_p ,i,j,k)   +  xinp(src-1,i,j,k)*diff2diff(E_le_p+i1+(src-1)*C%dof) 
                  xguess(E_ri_m ,i+1,j,k) = xguess(E_ri_m ,i+1,j,k) +  xinp(src-1,i,j,k)*diff2diff(E_ri_m+i1+(src-1)*C%dof) 
                  xguess(E_ri_p ,i+1,j,k) = xguess(E_ri_p ,i+1,j,k) +  xinp(src-1,i,j,k)*diff2diff(E_ri_p+i1+(src-1)*C%dof) 
                  xguess(E_ba_m ,i,j,k)   = xguess(E_ba_m ,i,j,k)   +  xinp(src-1,i,j,k)*diff2diff(E_ba_m+i1+(src-1)*C%dof) 
                  xguess(E_ba_p ,i,j,k)   = xguess(E_ba_p ,i,j,k)   +  xinp(src-1,i,j,k)*diff2diff(E_ba_p+i1+(src-1)*C%dof) 
                  xguess(E_fw_m ,i,j+1,k) = xguess(E_fw_m ,i,j+1,k) +  xinp(src-1,i,j,k)*diff2diff(E_fw_m+i1+(src-1)*C%dof) 
                  xguess(E_fw_p ,i,j+1,k) = xguess(E_fw_p ,i,j+1,k) +  xinp(src-1,i,j,k)*diff2diff(E_fw_p+i1+(src-1)*C%dof) 
                enddo
                !                  print *,'xguess after',xguess(:,i,j,k)
              endif
            enddo
          enddo
        enddo

        call DMDAVecRestoreArrayF90(C%da,local_guess,xguess,ierr) ;CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C%da,inp,xinp,ierr) ;CHKERRQ(ierr)

        call VecSet(guess,zero,ierr) ;CHKERRQ(ierr) ! reset global Vec

        call DMLocalToGlobalBegin(C%da,local_guess,ADD_VALUES, guess,ierr) ;CHKERRQ(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
        call DMLocalToGlobalEnd  (C%da,local_guess,ADD_VALUES, guess,ierr) ;CHKERRQ(ierr)

        call VecDestroy(local_guess,ierr) ;CHKERRQ(ierr)
    end subroutine

subroutine init_memory(incSolar,b,edir,ediff,abso,Mdir,Mdiff,edir_twostr,ediff_twostr,abso_twostr)
        Vec :: b,ediff,edir,abso,incSolar,edir_twostr,ediff_twostr,abso_twostr
        Mat :: Mdiff,Mdir

        call DMCreateGlobalVector(C_dir%da,incSolar,ierr) ; CHKERRQ(ierr)
        call DMCreateGlobalVector(C_diff%da,b,ierr)       ; CHKERRQ(ierr)
        call DMCreateGlobalVector(C_dir%da,edir,ierr)     ; CHKERRQ(ierr)
        call DMCreateGlobalVector(C_diff%da,ediff,ierr)       ; CHKERRQ(ierr)
        call DMCreateGlobalVector(C_one%da,abso,ierr)       ; CHKERRQ(ierr)

        call VecSet(incSolar,zero,ierr) ; CHKERRQ(ierr)
        call VecSet(b,zero,ierr)        ; CHKERRQ(ierr)
        call VecSet(edir,zero,ierr)     ; CHKERRQ(ierr)
        call VecSet(ediff,zero,ierr)    ; CHKERRQ(ierr)
        call VecSet(abso,zero,ierr)     ; CHKERRQ(ierr)

        call init_Matrix(Mdir,C_dir,"dir_")
        call init_Matrix(Mdiff,C_diff,"diff_")

        if(ltwostr) then
          call DMCreateGlobalVector(C_dir%da,edir_twostr,ierr)     ; CHKERRQ(ierr)
          call DMCreateGlobalVector(C_diff%da,ediff_twostr,ierr)       ; CHKERRQ(ierr)
          call DMCreateGlobalVector(C_one%da,abso_twostr,ierr)       ; CHKERRQ(ierr)
          call VecSet(edir_twostr,zero,ierr)     ; CHKERRQ(ierr)
          call VecSet(ediff_twostr,zero,ierr)    ; CHKERRQ(ierr)
          call VecSet(abso_twostr,zero,ierr)     ; CHKERRQ(ierr)
        endif
end subroutine

subroutine init_tenstream(icomm, Nx,Ny,Nz, dx,dy, phi0,theta0,albedo, dz1d, dz3d, nxproc, nyproc)
    ! vector sizes Nx, Ny Nz are either global domain size or if present(nxproc,nyproc) local domain size
    integer,intent(in) :: icomm
    integer(iintegers),intent(in) :: Nx,Ny,Nz
    real(ireals),intent(in) :: dx,dy, phi0,theta0,albedo

    real(ireals),optional :: dz1d(Nz), dz3d(Nx,Ny,Nz)   ! dz3d has to be local domain size
    integer(iintegers),optional :: nxproc(:), nyproc(:) ! size of local domains on each node

    integer(iintegers) :: i,j,k

    if(.not.linitialized) then

      call setup_petsc_comm
      call PetscInitialize(PETSC_NULL_CHARACTER,ierr) ;CHKERRQ(ierr)
      call init_mpi_data_parameters(PETSC_COMM_WORLD)

      call read_commandline_options()

      if(present(nxproc) .and. present(nyproc) ) then
        if(ldebug.and.myid.eq.0) print *,'nxproc',shape(nxproc),'::',nxproc
        if(ldebug.and.myid.eq.0) print *,'nyproc',shape(nyproc),'::',nyproc
        call setup_grid( Nx, Ny, Nz, nxproc,nyproc)
      else
        call setup_grid( max(minimal_dimension, Nx), max(minimal_dimension, Ny), Nz)
      endif

    endif

    atm%dx  = dx
    atm%dy  = dy
    atm%albedo = albedo

    if(.not.allocated(atm%dz) ) allocate(atm%dz( C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze ))
    if(.not.allocated(atm%l1d)) allocate(atm%l1d( C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze ) )

    if(present(dz1d)) then
      do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe
          atm%dz(i,j,:) = dz1d
        enddo
      enddo
    else if(present(dz3d)) then
      atm%dz = dz3d
    else
      print *,'have to give either dz1d or dz3d in routine call....'
      call exit(1)
    endif


    atm%l1d(:,:,C_one%ze) = twostr_ratio*atm%dz(:,:,C_one%ze).gt.atm%dx
    do k=C_one%ze-1,C_one%zs,-1
      do j=C_one%ys,C_one%ye
        do i=C_one%xs,C_one%xe
          if( atm%l1d(i,j,k+1) ) then !can only be 3D RT if below is a 3D layer
            atm%l1d(i,j,k)=.True.
          else
            !TODO this does actually not really make sense. I am not sure why this might not work. i.e. if it is necessary.
            atm%l1d(i,j,k) = twostr_ratio*atm%dz(i,j,k).gt.atm%dx
            if(atm%dz(i,j,k).lt.atm%dx/10._ireals) atm%l1d(i,j,k)=.True.
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
      call init_memory(incSolar,b,edir,ediff,abso,Mdir,Mdiff,edir_twostr,ediff_twostr,abso_twostr)

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
      logical :: lhave_planck

      lhave_planck = present(global_planck); call imp_bcast( lhave_planck, 0_mpiint, myid )

      ! Make sure that our domain has at least 3 entries in each dimension.... otherwise violates boundary conditions
      if(myid.eq.0) then
        call extend_arr(global_kabs)
        call extend_arr(global_ksca)
        call extend_arr(global_g)
        if(present(global_planck)) call extend_arr(global_planck)
      endif

      if(present(global_kabs) ) allocate( local_kabs (C_one%xs :C_one%xe , C_one%ys :C_one%ye  , C_one%zs :C_one%ze) )
      if(present(global_ksca) ) allocate( local_ksca (C_one%xs :C_one%xe , C_one%ys :C_one%ye  , C_one%zs :C_one%ze) )
      if(present(global_g   ) ) allocate( local_g    (C_one%xs :C_one%xe , C_one%ys :C_one%ye  , C_one%zs :C_one%ze) )
      if(lhave_planck .and. present(global_planck) ) allocate( local_planck   (C_one1%xs:C_one1%xe, C_one1%ys:C_one1%ye , C_one1%zs:C_one1%ze) )

      ! Scatter global optical properties to MPI nodes
      call local_optprop()
      ! Now global_fields are local to mpi subdomain.

      call set_optical_properties(local_kabs, local_ksca, local_g, local_planck)

      contains
          subroutine local_optprop()
                Vec :: local_vec
                PetscReal,pointer,dimension(:,:,:,:) :: xlocal_vec

                if(myid.eq.0.and.ldebug) print *,myid,'copying optprop: global to local :: shape kabs',shape(global_kabs),'xstart/end',C_one%xs,C_one%xe,'ys/e',C_one%ys,C_one%ye

                call DMCreateGlobalVector(C_one%da, local_vec, ierr) ; CHKERRQ(ierr)
                call scatterZerotoDM(global_kabs,C_one,local_vec)
                call DMDAVecGetArrayF90(C_one%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)
                local_kabs = xlocal_vec(0,:,:,:)
                call DMDAVecRestoreArrayF90(C_one%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)

                call scatterZerotoDM(global_ksca,C_one,local_vec)
                call DMDAVecGetArrayF90(C_one%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)
                local_ksca = xlocal_vec(0,:,:,:)
                call DMDAVecRestoreArrayF90(C_one%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)

                call scatterZerotoDM(global_g,C_one,local_vec)
                call DMDAVecGetArrayF90(C_one%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)
                local_g = xlocal_vec(0,:,:,:)
                call DMDAVecRestoreArrayF90(C_one%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)

                call VecDestroy(local_vec,ierr) ; CHKERRQ(ierr)

                if(lhave_planck) then
                  call DMCreateGlobalVector(C_one1%da, local_vec, ierr) ; CHKERRQ(ierr)
                  call scatterZerotoDM(global_planck,C_one1,local_vec)
                  call DMDAVecGetArrayF90(C_one1%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)
                  local_planck = xlocal_vec(0,:,:,:)
                  call DMDAVecRestoreArrayF90(C_one1%da ,local_vec ,xlocal_vec ,ierr)  ; CHKERRQ(ierr)
                  call VecDestroy(local_vec,ierr) ; CHKERRQ(ierr)
                endif
          end subroutine
          subroutine scatterZerotoDM(arr,C,vec)
              real(ireals),allocatable,dimension(:,:,:) :: arr
              type(t_coord) :: C
              Vec :: vec

              VecScatter :: scatter_context
              Vec :: natural,local
              PetscScalar,Pointer :: xloc(:)

              call DMDACreateNaturalVector(C%da, natural, ierr); CHKERRQ(ierr)
              call VecScatterCreateToZero(natural, scatter_context, local, ierr); CHKERRQ(ierr)

              if(myid.eq.0) then
                call VecGetArrayF90(local,xloc,ierr) ;CHKERRQ(ierr)
!                print *,myid,'shape of local',shape(xloc), 'shape of arr',shape(arr)
                xloc = reshape( arr , [ size(arr) ] )
                call VecRestoreArrayF90(local,xloc,ierr) ;CHKERRQ(ierr)
              endif

              call VecScatterBegin(scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); CHKERRQ(ierr)
              call VecScatterEnd  (scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); CHKERRQ(ierr)

              call DMDANaturalToGlobalBegin(C%da,natural, INSERT_VALUES, vec, ierr); CHKERRQ(ierr)
              call DMDANaturalToGlobalEnd  (C%da,natural, INSERT_VALUES, vec, ierr); CHKERRQ(ierr)
              
              call VecScatterDestroy(scatter_context, ierr); CHKERRQ(ierr)
              call VecDestroy(local,ierr); CHKERRQ(ierr)
              call VecDestroy(natural,ierr); CHKERRQ(ierr)

          end subroutine
    end subroutine

    subroutine extend_arr(arr)
        real(ireals),intent(inout),allocatable :: arr(:,:,:)
        real(ireals),allocatable :: tmp(:,:)
        integer(iintegers) :: dims(3),i

        dims = shape(arr)
        if( dims(1) .eq. 1 ) then
          allocate( tmp(dims(2),dims(3) ), SOURCE=arr(1,:,:) )
          deallocate(arr) 
          allocate( arr(minimal_dimension, dims(2), dims(3) ) )
          do i=1,minimal_dimension
            arr(i, :, :) = tmp
          enddo
          deallocate(tmp)
        endif

        dims = shape(arr)
        if( dims(2) .eq. 1 ) then
          allocate( tmp(dims(1),dims(3) ), SOURCE=arr(:,1,:) )
          deallocate(arr) 
          allocate( arr(dims(1), minimal_dimension, dims(3) ) )
          do i=1,minimal_dimension
            arr(:,i , :) = tmp
          enddo
          deallocate(tmp)
        endif
        if(any(shape(arr).lt.minimal_dimension) ) stop 'set_optprop -> extend_arr :: dimension is smaller than we support... please think of something here'
    end subroutine

    subroutine set_optical_properties(local_kabs, local_ksca, local_g, local_planck)
      real(ireals),intent(in),dimension(:,:,:),optional :: local_kabs, local_ksca, local_g
      real(ireals),intent(in),dimension(:,:,:),optional :: local_planck
      real(ireals) :: tau,kext,w0,g
      integer(iintegers) :: i,j,k

      if(.not.allocated(atm%op) )  allocate( atm%op       (C_one%xs :C_one%xe , C_one%ys :C_one%ye  , C_one%zs :C_one%ze) )
      if(present(local_kabs) ) atm%op(:,:,:)%kabs = local_kabs
      if(present(local_ksca) ) atm%op(:,:,:)%ksca = local_ksca
      if(present(local_g   ) ) atm%op(:,:,:)%g    = local_g   

      if(present(local_planck) ) then 
        if (.not.allocated(atm%planck) ) allocate( atm%planck   (C_one1%xs:C_one1%xe, C_one1%ys:C_one1%ye , C_one1%zs:C_one1%ze) )
        atm%planck = local_planck
      else
        if(allocated(atm%planck)) deallocate(atm%planck)
      endif

      if(ldebug) then
        if( (any([local_kabs,local_ksca,local_g].lt.zero)) .or. (any(isnan([local_kabs,local_ksca,local_g]))) ) then
          print *,myid,'set_optical_properties :: found illegal value in local_optical properties! abort!'
          do k=1,ubound(local_kabs,3)
            print *,myid,k,'local_kabs',local_kabs(:,:,k)
            print *,myid,k,'local_ksca',local_ksca(:,:,k)
          enddo
        endif
      endif
      if(ldebug) then
        if( (any([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g].lt.zero)) .or. (any(isnan([atm%op(:,:,:)%kabs,atm%op(:,:,:)%ksca,atm%op(:,:,:)%g]))) ) then
          print *,myid,'set_optical_properties :: found illegal value in optical properties! abort!'
        endif
      endif
      !      print *,'local_kabs   ',maxval(local_kabs   )  ,shape(local_kabs   )
      !      print *,'local_ksca   ',maxval(local_ksca   )  ,shape(local_ksca   )
      !      print *,'local_g      ',maxval(local_g      )  ,shape(local_g      )
      !      print *,'local_planck ',maxval(local_planck )  ,shape(local_planck ), shape(atm%planck)
      !      print *,'atm_kabs   ',maxval(atm%op%kabs  )  ,shape(atm%op%kabs  )
      !      print *,'atm_ksca   ',maxval(atm%op%ksca  )  ,shape(atm%op%ksca  )
      !      print *,'atm_g      ',maxval(atm%op%g     )  ,shape(atm%op%g     )
      !      print *,'atm_planck ',maxval(atm%planck   )  ,shape(atm%planck   )
      !      print *,'atm_dz     ',maxval(atm%dz       )  ,shape(atm%dz       )
      !
      !      print *, count(atm%l1d) , size(atm%l1d)
      !      if(ldebug.and.myid.eq.0) print *,'init local optprop:', shape(local_kabs), '::', shape(atm%op)

      ! Make space for deltascaled optical properties
      if(.not.allocated(atm%delta_op) ) allocate( atm%delta_op (C_one%xs :C_one%xe , C_one%ys :C_one%ye  , C_one%zs :C_one%ze) )
      atm%delta_op = atm%op
      call delta_scale(atm%delta_op(:,:,:)%kabs, atm%delta_op(:,:,:)%ksca, atm%delta_op(:,:,:)%g ) !todo should we instead use strong deltascaling? -- what gives better results? or is it as good?

      if(ldebug) then
        if( (any([atm%delta_op(:,:,:)%kabs,atm%delta_op(:,:,:)%ksca,atm%delta_op(:,:,:)%g].lt.zero)) .or. (any(isnan([atm%delta_op(:,:,:)%kabs,atm%delta_op(:,:,:)%ksca,atm%delta_op(:,:,:)%g]))) ) then
          print *,myid,'set_optical_properties :: found illegal value in delta_scaled optical properties! abort!'
        endif
      endif
      if(luse_eddington) then
        if(.not.allocated(atm%a11) ) allocate(atm%a11 (C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze )) ! allocate space for twostream coefficients
        if(.not.allocated(atm%a12) ) allocate(atm%a12 (C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze )) 
        if(.not.allocated(atm%a13) ) allocate(atm%a13 (C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze )) 
        if(.not.allocated(atm%a23) ) allocate(atm%a23 (C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze )) 
        if(.not.allocated(atm%a33) ) allocate(atm%a33 (C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze )) 
        if(.not.allocated(atm%g1 ) ) allocate(atm%g1  (C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze )) 
        if(.not.allocated(atm%g2 ) ) allocate(atm%g2  (C_one%xs:C_one%xe, C_one%ys:C_one%ye, C_one%zs:C_one%ze )) 
      endif
        
      if(luse_eddington) then
        do k=lbound(atm%op,3),ubound(atm%op,3)
          do j=lbound(atm%op,2),ubound(atm%op,2)
            do i=lbound(atm%op,1),ubound(atm%op,1)
              if( atm%l1d(i,j,k) ) then
                kext = atm%delta_op(i,j,k)%kabs + atm%delta_op(i,j,k)%ksca
                w0   = atm%delta_op(i,j,k)%ksca / kext
                tau  = atm%dz(i,j,k)* kext
                g    = atm%delta_op(i,j,k)%g 
                call eddington_coeff_fab ( tau , w0, g, sun%costheta, & 
                    atm%a11(i,j,k),          &
                    atm%a12(i,j,k),          &
                    atm%a13(i,j,k),          &
                    atm%a23(i,j,k),          &
                    atm%a33(i,j,k),          &
                    atm%g1(i,j,k),           &
                    atm%g2(i,j,k) )
              endif
            enddo
          enddo
        enddo
      endif

      if(ldebug .and. myid.eq.0) then
        do k=lbound(atm%op,3),ubound(atm%op,3)
          if(present(local_planck)) then
            print *,myid,'Optical Properties:',k,'dz',atm%dz(0,0,k),atm%l1d(0,0,k),'k',&
                minval(atm%delta_op(:,:,k)%kabs),minval(atm%delta_op(:,:,k)%ksca),minval(atm%delta_op(:,:,k)%g),&
                maxval(atm%delta_op(:,:,k)%kabs),maxval(atm%delta_op(:,:,k)%ksca),maxval(atm%delta_op(:,:,k)%g),&
                '::',minval(atm%planck(:,:,k)),maxval(atm%planck(:,:,k))
          else    
            print *,myid,'Optical Properties:',k,'dz',atm%dz(0,0,k),atm%l1d(0,0,k),'k',&
                minval(atm%delta_op(:,:,k)%kabs),minval(atm%delta_op(:,:,k)%ksca),minval(atm%delta_op(:,:,k)%g),&
                maxval(atm%delta_op(:,:,k)%kabs),maxval(atm%delta_op(:,:,k)%ksca),maxval(atm%delta_op(:,:,k)%g),&
                '::',minval(atm%a33(:,:,k)),maxval(atm%a33(:,:,k))
          endif
        enddo
      endif

    end subroutine

    subroutine solve_tenstream(edirTOA,solution_uid)
        real(ireals),intent(in) :: edirTOA
        integer(iintegers),optional,intent(in) :: solution_uid
        logical :: loaded

            if(ltwostr) then
              call twostream(edirTOA)
              call calc_flx_div(edir_twostr,ediff_twostr,abso_twostr)

              call scale_flx(edir_twostr  ,C_dir)
              call scale_flx(ediff_twostr ,C_diff)

              if(luse_twostr_guess) then
                call VecCopy(edir_twostr  ,edir ,ierr) ;CHKERRQ(ierr)
                call VecCopy(ediff_twostr ,ediff,ierr) ;CHKERRQ(ierr)
                !                call set_diff_initial_guess(ediff_twostr, ediff, C_diff)
              endif

              if(myid.eq.0) print *,'twostream calculation done'
            endif

            ! ---------------------------- Edir  -------------------
            if(edirTOA.gt.zero .and. sun%theta.ge.zero) then

              if( present(solution_uid) ) loaded = load_solution(solution_uid)

              call PetscLogStagePush(logstage(1),ierr) ;CHKERRQ(ierr)
              call setup_incSolar(incSolar,edirTOA)
              call set_dir_coeff(Mdir,C_dir)

              call setup_ksp(kspdir,C_dir,Mdir,linit_kspdir, "dir_")

              call PetscLogStagePush(logstage(3),ierr) ;CHKERRQ(ierr)
              call solve(kspdir,incSolar,edir)
              call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
            else
              call VecSet(edir,zero,ierr)
            endif
            ! ---------------------------- Source Term -------------
            call setup_b(edir,b)

            ! ---------------------------- Ediff -------------------
            call set_diff_coeff(Mdiff,C_diff)

            call setup_ksp(kspdiff,C_diff,Mdiff,linit_kspdiff, "diff_")

            call PetscLogStagePush(logstage(5),ierr) ;CHKERRQ(ierr)
            call solve(kspdiff, b, ediff)
            call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

            if(present(solution_uid) ) call save_solution(solution_uid)

            ! ---------------------------- Absorption and Rescaling-
            call calc_flx_div(edir,ediff,abso)

            call scale_flx(edir,C_dir)
            call scale_flx(ediff,C_diff)
        end subroutine

        subroutine destroy_tenstream()
            if(ltwostr) then
              call VecDestroy(edir_twostr    , ierr) ;CHKERRQ(ierr)
              call VecDestroy(ediff_twostr   , ierr) ;CHKERRQ(ierr)
              call VecDestroy(abso_twostr    , ierr) ;CHKERRQ(ierr)
            endif

            call KSPDestroy(kspdir , ierr) ;CHKERRQ(ierr); linit_kspdir =.False.
            call KSPDestroy(kspdiff, ierr) ;CHKERRQ(ierr); linit_kspdiff=.False.

            call VecDestroy(incSolar , ierr) ;CHKERRQ(ierr)
            call VecDestroy(b        , ierr) ;CHKERRQ(ierr)
            call VecDestroy(edir     , ierr) ;CHKERRQ(ierr)
            call VecDestroy(ediff    , ierr) ;CHKERRQ(ierr)
            call VecDestroy(abso     , ierr) ;CHKERRQ(ierr)
            call MatDestroy(Mdir     , ierr) ;CHKERRQ(ierr)
            call MatDestroy(Mdiff    , ierr) ;CHKERRQ(ierr)

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
        end subroutine

        subroutine tenstream_get_result(redir,redn,reup,rabso)
            real(ireals),dimension(:,:,:),intent(out),optional :: redir,redn,reup,rabso
            PetscScalar,pointer,dimension(:,:,:,:) :: x

            if(present(redir)) then
              call DMDAVecGetArrayF90(C_dir%da,edir,x,ierr) ;CHKERRQ(ierr)
              redir = sum(x(i0:i3,:,:,:),dim=1)/4
              call DMDAVecRestoreArrayF90(C_dir%da,edir,x,ierr) ;CHKERRQ(ierr)
            endif

            if(present(redn).or.present(reup)) then
              call DMDAVecGetArrayF90(C_diff%da,ediff,x,ierr) ;CHKERRQ(ierr)
              if(present(redn) )redn = x(i1,:,:,:)
              if(present(reup) )reup = x(i0,:,:,:)
              call DMDAVecRestoreArrayF90(C_diff%da,ediff,x,ierr) ;CHKERRQ(ierr)
            endif

            if(present(rabso)) then
              call DMDAVecGetArrayF90(C_one%da,abso,x,ierr) ;CHKERRQ(ierr)
              rabso = x(i0,:,:,:)
              call DMDAVecRestoreArrayF90(C_one%da,abso,x,ierr) ;CHKERRQ(ierr)
            endif
        end subroutine

        function load_solution(uid) result(loaded)
            integer(iintegers),intent(in) :: uid
            logical :: loaded
            real(ireals) :: norm1,norm2
            if(.not.lenable_solutions) return

            if(uid.gt.size(solutions)) then
              print *,'unique identifier exceeds container size.... you might want to grow it...',uid,size(solutions)
              call exit(1)
            endif

            if( .not. solutions(uid)%lset ) then
              loaded = .False.
              return
            else
              if(ldebug.and.myid.eq.0) print *,'Loading Solution for uid',uid
              call VecCopy(solutions(uid)%edir , edir , ierr) ;CHKERRQ(ierr)
              call VecCopy(solutions(uid)%ediff, ediff, ierr) ;CHKERRQ(ierr)
            endif
            loaded = .True.

!            call VecNorm(edir,NORM_2,norm1,ierr)
!            call VecNorm(solutions(uid)%edir,NORM_2,norm2,ierr)
!            print *,'loading vectors edir norms',norm1,norm2
!            call VecNorm(ediff,NORM_2,norm1,ierr)
!            call VecNorm(solutions(uid)%ediff,NORM_2,norm2,ierr)
!            print *,'loading vectors ediff norms',norm1,norm2
        end function
        subroutine save_solution(uid)
            integer(iintegers),intent(in) :: uid
            character(100) :: vecname
            real(ireals) :: norm1,norm2
            if(.not.lenable_solutions) return

            if(uid.gt.size(solutions)) then
              print *,'unique identifier exceeds container size.... you might want to grow it...',uid,size(solutions)
              call exit(1)
            endif

            if( .not. solutions(uid)%lset ) then
!              if(myid.eq.0 .and. ldebug) print *,'duplicating vectors to store solution',uid
              call VecDuplicate(edir , solutions(uid)%edir , ierr) ;CHKERRQ(ierr)
              call VecDuplicate(ediff, solutions(uid)%ediff, ierr) ;CHKERRQ(ierr)
              solutions(uid)%lset=.True.
            endif

!            call VecNorm(edir,NORM_2,norm1,ierr)
!            call VecNorm(solutions(uid)%edir,NORM_2,norm2,ierr)
!            print *,'before saving vectors edir norms',norm1,norm2
!            call VecNorm(ediff,NORM_2,norm1,ierr)
!            call VecNorm(solutions(uid)%ediff,NORM_2,norm2,ierr)
!            print *,'before saving vectors ediff norms',norm1,norm2


            if(ldebug.and.myid.eq.0) print *,'Saving Solution for uid',uid
            call VecCopy(edir , solutions(uid)%edir , ierr) ;CHKERRQ(ierr)
            call VecCopy(ediff, solutions(uid)%ediff, ierr) ;CHKERRQ(ierr)

!            call VecNorm(edir,NORM_2,norm1,ierr)
!            call VecNorm(solutions(uid)%edir,NORM_2,norm2,ierr)
!            print *,'saving vectors edir norms',norm1,norm2
!            call VecNorm(ediff,NORM_2,norm1,ierr)
!            call VecNorm(solutions(uid)%ediff,NORM_2,norm2,ierr)
!            print *,'saving vectors ediff norms',norm1,norm2


!            write(vecname,FMT='("edir",I0)') uid
!            call PetscObjectSetName(edir,vecname,ierr) ; CHKERRQ(ierr)
!            call vec_to_hdf5(edir)
!
!            write(vecname,FMT='("ediff",I0)') uid
!            call PetscObjectSetName(ediff,vecname,ierr) ; CHKERRQ(ierr)
!            call vec_to_hdf5(ediff)
!
!            write(vecname,FMT='("abso",I0)') uid
!            call PetscObjectSetName(abso,vecname,ierr) ; CHKERRQ(ierr)
!            call vec_to_hdf5(abso)
!
!            write(vecname,FMT='("b",I0)') uid
!            call PetscObjectSetName(b,vecname,ierr) ; CHKERRQ(ierr)
!            call vec_to_hdf5(b)
!
!            write(vecname,FMT='("incSolar",I0)') uid
!            call PetscObjectSetName(incSolar,vecname,ierr) ; CHKERRQ(ierr)
!            call vec_to_hdf5(incSolar)
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
