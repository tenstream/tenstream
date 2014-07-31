module m_tenstream

      use m_data_parameters, only : ireals,iintegers,       &
        imp_comm, myid, numnodes,init_mpi_data_parameters,mpiint, &
        zero,one,nil,i0,i1,i2,i3,i4,i5,i6,i7,i8,i10,pi

      use m_twostream, only: delta_eddington_twostream
      use m_helper_functions, only: deg2rad,approx,imp_bcast,rmse,delta_scale
      use m_eddington, only : eddington_coeff_rb
      use m_optprop_parameters, only : ldelta_scale
      use m_optprop, only : t_optprop_1_2,t_optprop_8_10
      use m_tenstream_options, only : ltwostr, luse_twostr_guess, luse_eddington, twostr_ratio

      implicit none
#include "finclude/petsc.h90"

      private
      public :: init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
                edir,ediff,abso,&
                edir_twostr,ediff_twostr,abso_twostr

      PetscInt,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9
      PetscInt,parameter :: istartpar=i1, jstartpar=i1

      logical,parameter :: ldebug=.True.,lcycle_dir=.True.
      logical,parameter :: lprealloc=.True.

      type coord
        PetscInt :: xs,xe,ys,ye,zs,ze     ! local domain start and end indices
        PetscInt :: xm,ym,zm              ! size of local domain
        PetscInt :: gxs,gys,gzs           ! domain indices including ghost points
        PetscInt :: gxm,gym,gzm           ! size of local domain including ghosts
        PetscInt :: glob_xm,glob_ym,glob_zm ! global domain size
        PetscInt :: dof,dim               ! degrees of freedom of Petsc Domain, dimension of dmda
        DM :: da                          ! The Domain Decomposition Object
        PetscMPIInt,allocatable :: neighbors(:) ! all 3d neighbours( (x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
      end type

      type(coord) :: C_dir,C_diff,C_one

      PetscErrorCode :: ierr

      type t_optprop
        real(ireals) :: kabs,ksca,g
      end type

      type t_atmosphere
        type(t_optprop) , allocatable , dimension(:,:,:) :: op
        type(t_optprop) , allocatable , dimension(:,:,:) :: delta_op
        real(ireals)    , allocatable , dimension(:,:,:) :: a11, a12, a13, a23, a33
        real(ireals)    , allocatable , dimension(:) :: dz
        logical         , allocatable , dimension(:) :: l1d
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

      contains 
      subroutine setup_grid(Nx,Ny,Nz)
        PetscInt,intent(in) :: Nx,Ny,Nz
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

        if(myid.eq.0.and.ldebug) print *,myid,'DMDA grid ready'
        contains
          subroutine setup_dmda(C, Nx, Ny, Nz, boundary, dof)
              type(coord) :: C
              PetscInt,intent(in) :: Nx,Ny,Nz,dof

              DMBoundaryType :: boundary
              PetscInt,parameter :: stencil_size=1

              PetscInt :: decomp_x=PETSC_DECIDE, decomp_y=PETSC_DECIDE

              C%dof = i1*dof

              call DMDACreate3d( imp_comm  ,                                           &
                                boundary           , boundary           , bn                 , &
                                DMDA_STENCIL_BOX  ,                                            &
                                Nx     , Ny     , i1*Nz                , &
                                decomp_x          , decomp_y          , i1                 , &
                                C%dof              , stencil_size       ,                      &
                                PETSC_NULL_INTEGER , PETSC_NULL_INTEGER , PETSC_NULL_INTEGER , &
                                C%da               , ierr) ;CHKERRQ(ierr)
              call setup_coords(C)
              call DMSetup(C%da,ierr) ;CHKERRQ(ierr)
!              call DMSetMatType(C%da, MATAIJ, ierr); CHKERRQ(ierr)
              if(lprealloc) call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ;CHKERRQ(ierr)

              call DMSetFromOptions(C%da, ierr) ; CHKERRQ(ierr)
          end subroutine
        subroutine setup_coords(C)
          type(coord) :: C

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
        type(coord) :: C
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
        type(coord),intent(in) :: C
        PetscInt :: i,j,k,dof
        MatStencil :: row(4,1), col(4,1)
        PetscScalar :: v(1)

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
        if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients diagonally ... done'
      end subroutine
      subroutine setup_diff_preallocation(d_nnz,o_nnz,C)
        PetscInt,allocatable :: o_nnz(:)
        PetscInt,allocatable :: d_nnz(:)
        type(coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize

        PetscScalar, pointer :: xx_v_o(:),xx_v_d(:)

        PetscInt,parameter :: ind(9)=[E_up,E_le_m,E_le_p,E_ri_m,E_ri_p,E_ba_m,E_ba_p,E_fw_m,E_fw_p]
        PetscInt :: i,j,k,li,lj,lk

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
            do i=C%xs,C%xe      ! i,j,k indices are defined on petsc global grid
              li = istartpar+i-C%xs ! l-indices are used for local arrays, 
              lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
              lk = 1+k-C%zs         !
              if( atm%l1d(lk) ) then
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
        type(coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize,i,j,k,li,lj,lk

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
             do i=C%xs,C%xe      ! i,j,k indices are defined on petsc global grid
               li = istartpar+i-C%xs ! l-indices are used for local arrays, 
               lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
               lk = 1+k-C%zs         !
               if( atm%l1d(lk) ) then
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
        type(coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize,i,j,k,li,lj,lk

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
               li = istartpar+i-C%xs ! l-indices are used for local arrays, 
               lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
               lk = 1+k-C%zs         !
               if( atm%l1d(lk) ) then
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
        if(ldebug) print *,'setup_dir_inc done'
end subroutine

subroutine set_dir_coeff(A,C)
        Mat :: A
        type(coord) :: C

        PetscInt :: i,j,k,src,dst, li,lj,lk

        MatStencil :: row(4,C%dof)  ,col(4,C%dof)
        PetscReal :: v(C%dof**2),coeffs(C%dof**2),norm
        PetscInt,parameter :: entries(8)=[0,8,16,24,32,40,48,56]
        PetscInt,parameter :: entries1d(4)=[0,4,8,12]
        PetscReal :: twostr_coeff(1)

!        if(ldebug) print *,myid,'DEBUG(set_dir_coeff) DMDA',C%xs,C%xe,C%ys,C%ye,C%zs,C%ze-1
        call PetscLogStagePush(logstage(2),ierr) ;CHKERRQ(ierr)

        row=PETSC_NULL_INTEGER
        col=PETSC_NULL_INTEGER
        v=PETSC_NULL_REAL
        coeffs=PETSC_NULL_REAL

        do k=C%zs,C%ze-1
        do j=C%ys,C%ye
        do i=C%xs,C%xe        ! i,j,k indices are defined on petsc global grid
          li = istartpar+i-C%xs ! l-indices are used for local arrays, 
          lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
          lk = 1+k-C%zs         !
          !          print *,'setting direct coeff',li,lj,lk,'PETSC_grid',i,j,k,'optprop',lbound(optprop,1),ubound(optprop,1),':',lbound(optprop,2),ubound(optprop,2),':',lbound(optprop,3),ubound(optprop,3)

          dst = 1 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 2 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 3 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 4 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 5 ; row(MatStencil_i,dst) = i+sun%xinc ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
          dst = 6 ; row(MatStencil_i,dst) = i+sun%xinc ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
          dst = 7 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j+sun%yinc  ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
          dst = 8 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j+sun%yinc  ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid

          src = 1 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 2 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 3 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 4 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 5 ; col(MatStencil_i,src) = i+1-sun%xinc ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
          src = 6 ; col(MatStencil_i,src) = i+1-sun%xinc ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
          src = 7 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j+1-sun%yinc ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:
          src = 8 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j+1-sun%yinc ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:

          coeffs=PETSC_NULL_REAL
          if( atm%l1d(lk) ) then
            if(luse_eddington) then
              coeffs([0, 8+1, 16+2, 24+3 ]+i1) = atm%a33(li,lj,lk) ! only use the four vertical tiles with a33
            else
              call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.True., twostr_coeff, atm%l1d(lk), [sun%symmetry_phi, sun%theta] )
              coeffs([0, 8+1, 16+2, 24+3 ]+i1) = twostr_coeff(1) ! only use the four vertical tiles with direct transmission
            endif
          else
            call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.True., coeffs, atm%l1d(lk), [sun%symmetry_phi, sun%theta])
          endif

          if(ldebug) then
            do src=1,C%dof
              norm = sum( coeffs((src-1)*C%dof+1:src*C%dof) )
              if( ldebug .and. real(norm).gt.real(one) ) then
                  print *,'sum(src==',src,') gt one',norm
                  stop 'omg.. shouldnt be happening'
              endif
            enddo
          endif

          where(approx(coeffs,zero))
            coeffs = PETSC_NULL_REAL
          end where

          if( atm%l1d(lk) ) then
            ! reorder coeffs from src-ordered to dst-ordered
            do src=1,4
              v(entries1d+src) = coeffs( i1+(src-i1)*C%dof : i5+(src-i1)*C%dof )
            enddo

            call MatSetValuesStencil(A,i4, row,i4, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)
          else

            ! reorder coeffs from src-ordered to dst-ordered
            do src=1,C%dof
              v(entries+src) = coeffs( i1+(src-i1)*C%dof : src*C%dof )
            enddo
            call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)
          endif

        enddo ; enddo ; enddo

        if(myid.eq.0.and.ldebug) print *,myid,'setup_direct_matrix done'

        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)  

        call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
        end subroutine

        subroutine setup_incSolar(incSolar,edirTOA)
          Vec :: incSolar
          real(ireals),intent(in) :: edirTOA

          PetscInt :: i,j,li,lj
          PetscScalar,pointer,dimension(:,:,:,:) :: xincSolar

          PetscReal :: Az
          Az = atm%dx*atm%dy

          call VecSet(incSolar,zero,ierr) ;CHKERRQ(ierr)

          call DMDAVecGetArrayF90(C_dir%da,incSolar,xincSolar,ierr) ;CHKERRQ(ierr)
          do j=C_dir%ys,C_dir%ye
            do i=C_dir%xs,C_dir%xe        ! i,j,k indices are defined on petsc global grid
              li = istartpar+i-C_dir%xs ! l-indices are used for local arrays, 
              lj = jstartpar+j-C_dir%ys ! which are defined on the cosmo grid

              xincSolar(i0:i3,i,j,C_dir%zs) = edirTOA* Az * .25_ireals * sun%costheta
            enddo
          enddo

          call DMDAVecRestoreArrayF90(C_dir%da,incSolar,xincSolar,ierr) ;CHKERRQ(ierr)

          if(myid.eq.0 .and. ldebug) print *,myid,'Setup of IncSolar done',edirTOA

        end subroutine

subroutine set_diff_coeff(A,C)
  Mat :: A
  type(coord) :: C

  PetscInt :: i,j,k,src,dst, li,lj,lk

  MatStencil :: row(4,0:C%dof-1)  ,col(4,0:C%dof-1)
  PetscReal :: v(C%dof**2),coeffs(C%dof**2),norm,twostr_coeff(4)
  PetscInt,parameter :: entries(10)=[0,10,20,30,40,50,60,70,80,90]
!  PetscInt,parameter :: entries1d(2)=[0,2]

  ! if(ldebug) print *,myid,'DEBUG(set_dir_coeff) jspec',jspec,'igas',igas,'isub',isub
  call PetscLogStagePush(logstage(4),ierr) ;CHKERRQ(ierr)

  row=PETSC_NULL_INTEGER
  col=PETSC_NULL_INTEGER
  v  =PETSC_NULL_REAL

  if(myid.eq.0.and.ldebug) print *,myid,'Setting coefficients for diffuse Light'!,row,col,v

  ! i,j,k indices are defined on petsc global grid
  do k=C%zs,C%ze-1
    do j=C%ys,C%ye
      do i=C%xs,C%xe
        li = istartpar+i-C%xs ! l-indices are used for local arrays, 
        lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
        lk = 1+k-C%zs         !
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

        coeffs = PETSC_NULL_REAL
        if( atm%l1d(lk) ) then
          if(luse_eddington ) then
!            call rodents(optprop(li,lj,lk)%kext1*newgrid%dz(lk), optprop(li,lj,lk)%w1, optprop(li,lj,lk)%g1, costheta0, edd_coeff )
            coeffs([ 0,1, 10,11 ]+i1) = [ atm%a12(li,lj,lk), atm%a11(li,lj,lk), atm%a11(li,lj,lk), atm%a12(li,lj,lk) ]
          else
            call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.False., twostr_coeff, atm%l1d(lk))
            coeffs([ 0,1, 10,11 ]+i1) = twostr_coeff
          endif
        else
          call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.False., coeffs, atm%l1d(lk))
        endif

        where(approx(coeffs,zero))
          coeffs = PETSC_NULL_REAL
        end where

        if(ldebug) then
          do dst=1,C%dof
            norm = sum( coeffs((dst-1)*C%dof+1:dst*C%dof) )
            if( ldebug ) then
              if( real(norm).gt.real(one) ) then
!                coeffs((dst-1)*C%dof+1:dst*C%dof)  = coeffs((dst-1)*C%dof+1:dst*C%dof) / (norm+1e-8_ireals)
                print *,'diffuse sum(dst==',dst,') gt one',norm
                stop 'omg.. shouldnt be happening'
              endif
            endif
          enddo
        endif


        if( atm%l1d(lk) ) then
          do src=1,2
            v(entries+src) = coeffs( i1+(src-i1)*C%dof : i3+(src-i1)*C%dof )
          enddo
          call MatSetValuesStencil(A,i2, row,i2, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)
        else
          ! reorder coeffs from src-ordered to dst-ordered
          do src=1,C%dof
            v(entries+src) = coeffs( i1+(src-i1)*C%dof : src*C%dof )
          enddo
          call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr) ;CHKERRQ(ierr)
        endif


      enddo ; enddo ; enddo
      if(myid.eq.0.and.ldebug) print *,myid,'setup_diffuse_matrix done'

      if(myid.eq.0.and.ldebug) print *,myid,'Final diffuse Matrix Assembly:'
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ;CHKERRQ(ierr)
!      call MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE, ierr) ;CHKERRQ(ierr)
      call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
end subroutine
 
subroutine setup_b(edir,b)
        Vec :: edir
        Vec :: local_b,b

        PetscScalar,pointer,dimension(:,:,:,:) :: &
        xsrc,xedir
        PetscInt :: li,lj,lk
        PetscReal :: coeffs(C_dir%dof*C_diff%dof),twostr_coeff(2)
        
        PetscInt :: i,j,k,src
        call PetscLogStagePush(logstage(6),ierr) ;CHKERRQ(ierr)
        
        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly...'

        call DMCreateLocalVector(C_diff%da,local_b,ierr) ;CHKERRQ(ierr)
        call VecSet(local_b,zero,ierr) ;CHKERRQ(ierr)

        call DMDAVecGetArrayF90(C_diff%da,local_b,xsrc,ierr) ;CHKERRQ(ierr)
        call DMDAVecGetArrayF90(C_dir%da,edir,xedir,ierr) ;CHKERRQ(ierr)

        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients'
        do k=C_diff%zs,C_diff%ze-1 
          do j=C_diff%ys,C_diff%ye         
            do i=C_diff%xs,C_diff%xe    
              li = istartpar+i-C_diff%xs
              lj = jstartpar+j-C_diff%ys
              lk = i1+k-C_diff%zs        

              coeffs=PETSC_NULL_REAL
              if( atm%l1d(lk) ) then
                if(luse_eddington ) then
!                  call rodents(optprop(li,lj,lk)%kext1*newgrid%dz(lk), optprop(li,lj,lk)%w1, optprop(li,lj,lk)%g1, costheta0, edd_coeff )
                  ! Only transport the 4 tiles from dir0 to the Eup and Edn
                  do src=1,4
                    coeffs(E_up  +i1+(src-1)*C_diff%dof) = atm%a13(li,lj,lk)
                    coeffs(E_dn  +i1+(src-1)*C_diff%dof) = atm%a23(li,lj,lk)
                  enddo

                else
                  call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.False., twostr_coeff, atm%l1d(lk), [sun%symmetry_phi, sun%theta])
                  do src=1,4
                    coeffs(E_up  +i1+(src-1)*C_diff%dof) = twostr_coeff(1)
                    coeffs(E_dn  +i1+(src-1)*C_diff%dof) = twostr_coeff(2)
                  enddo
                endif

              else
                call get_coeff(atm%delta_op(li,lj,lk), atm%dz(lk),.False., coeffs, atm%l1d(lk), [sun%symmetry_phi, sun%theta] )
              endif

              if( atm%l1d(lk) ) then
                do src=1,C_dir%dof-4
                  xsrc(E_up   ,i,j,k)   = xsrc(E_up   ,i,j,k)   +  xedir(src-1,i,j,k)*coeffs(E_up  +i1+(src-1)*C_diff%dof)
                  xsrc(E_dn   ,i,j,k+1) = xsrc(E_dn   ,i,j,k+1) +  xedir(src-1,i,j,k)*coeffs(E_dn  +i1+(src-1)*C_diff%dof)
                enddo
              else
                do src=1,C_dir%dof
                  xsrc(E_up   ,i,j,k)   = xsrc(E_up   ,i,j,k)   +  xedir(src-1,i,j,k)*coeffs(E_up  +i1+(src-1)*C_diff%dof)
                  xsrc(E_dn   ,i,j,k+1) = xsrc(E_dn   ,i,j,k+1) +  xedir(src-1,i,j,k)*coeffs(E_dn  +i1+(src-1)*C_diff%dof)
                  xsrc(E_le_m ,i,j,k)   = xsrc(E_le_m ,i,j,k)   +  xedir(src-1,i,j,k)*coeffs(E_le_m+i1+(src-1)*C_diff%dof)
                  xsrc(E_le_p ,i,j,k)   = xsrc(E_le_p ,i,j,k)   +  xedir(src-1,i,j,k)*coeffs(E_le_p+i1+(src-1)*C_diff%dof)
                  xsrc(E_ri_m ,i+1,j,k) = xsrc(E_ri_m ,i+1,j,k) +  xedir(src-1,i,j,k)*coeffs(E_ri_m+i1+(src-1)*C_diff%dof)  
                  xsrc(E_ri_p ,i+1,j,k) = xsrc(E_ri_p ,i+1,j,k) +  xedir(src-1,i,j,k)*coeffs(E_ri_p+i1+(src-1)*C_diff%dof)  
                  xsrc(E_ba_m ,i,j,k)   = xsrc(E_ba_m ,i,j,k)   +  xedir(src-1,i,j,k)*coeffs(E_ba_m+i1+(src-1)*C_diff%dof)  
                  xsrc(E_ba_p ,i,j,k)   = xsrc(E_ba_p ,i,j,k)   +  xedir(src-1,i,j,k)*coeffs(E_ba_p+i1+(src-1)*C_diff%dof)  
                  xsrc(E_fw_m ,i,j+1,k) = xsrc(E_fw_m ,i,j+1,k) +  xedir(src-1,i,j,k)*coeffs(E_fw_m+i1+(src-1)*C_diff%dof)  
                  xsrc(E_fw_p ,i,j+1,k) = xsrc(E_fw_p ,i,j+1,k) +  xedir(src-1,i,j,k)*coeffs(E_fw_p+i1+(src-1)*C_diff%dof) 
                enddo
              endif

            enddo
          enddo
        enddo

        ! Ground Albedo reflecting direct radiation, the diffuse part is considered by the solver(Matrix)
        do i=C_diff%xs,C_diff%xe     
          do j=C_diff%ys,C_diff%ye     
            li = istartpar+i-C_diff%xs 
            lj = jstartpar+j-C_diff%ys 
            k = C_diff%ze
            xsrc(E_up   ,i,j,k) = sum(xedir(i0:i3,i,j,k))*atm%albedo
          enddo
        enddo

        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients ...done'

        call DMDAVecRestoreArrayF90(C_dir%da,edir,xedir,ierr) ;CHKERRQ(ierr)
        call DMDAVecRestoreArrayF90(C_diff%da,local_b,xsrc,ierr) ;CHKERRQ(ierr)

        call VecSet(b,zero,ierr) ;CHKERRQ(ierr) ! reset global Vec

        call DMLocalToGlobalBegin(C_diff%da,local_b,ADD_VALUES, b,ierr) ;CHKERRQ(ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
        call DMLocalToGlobalEnd  (C_diff%da,local_b,ADD_VALUES, b,ierr) ;CHKERRQ(ierr)

        call VecDestroy(local_b,ierr) ;CHKERRQ(ierr)

        call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly done'
end subroutine

subroutine calc_flx_div(edir,ediff,abso)
        Vec :: edir,ediff,abso
        PetscReal,pointer,dimension(:,:,:,:) :: xediff,xedir,xabso
        PetscInt :: i,j,k,li,lj,lk
        Vec :: ledir,lediff ! local copies of vectors, including ghosts
        PetscReal :: div2(13)
        PetscReal :: Volume

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
          lk = i1+k-C_one%zs
          Volume = atm%dx * atm%dy * atm%dz(lk)

          do j=C_one%ys,C_one%ye         
            lj = jstartpar+j-C_one%ys 

            do i=C_one%xs,C_one%xe      
              li = istartpar+i-C_one%xs 

              ! Divergence    =                       Incoming                -       Outgoing
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

              xabso(i0,i,j,k) = sum(div2) / Volume
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

      if(myid.eq.0) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      if(reason.eq.KSP_DIVERGED_ITS) then
        if(myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
        return
      endif

      if(reason.le.0) then
              if(myid.eq.0) print *,myid,'Resetted initial guess to zero and try again:'
              call VecSet(x,zero,ierr) ;CHKERRQ(ierr)
              call KSPSolve(ksp,b,x,ierr) ;CHKERRQ(ierr)
              call KSPGetIterationNumber(ksp,iter,ierr) ;CHKERRQ(ierr)
              call KSPGetConvergedReason(ksp,reason,ierr) ;CHKERRQ(ierr)
              if(myid.eq.0) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      endif

      if(reason.le.0) then
              print *,'*********************************************************** SOLVER did NOT converge :( ***************************************************88'
              call exit()
      endif
end subroutine
subroutine setup_ksp(ksp,C,A,linit, prefix)
      KSP :: ksp
      type(coord) :: C
      Mat :: A
      logical :: linit

      MatNullSpace :: nullspace
      PetscReal,parameter :: rtol=1e-7, atol=1e-25
      PetscInt,parameter :: maxiter=1000
      character(len=*),optional :: prefix

      if(linit) return
      if(myid.eq.0.and.ldebug) print *,'Setup KSP'

      call KSPCreate(imp_comm,ksp,ierr) ;CHKERRQ(ierr)
      if(present(prefix) ) call KSPAppendOptionsPrefix(ksp,trim(prefix),ierr) ;CHKERRQ(ierr)

      call KSPSetOperators(ksp,A,A,ierr) ;CHKERRQ(ierr)
      call KSPSetDM(ksp,C%da,ierr) ;CHKERRQ(ierr)
      call KSPSetDMActive(ksp,PETSC_FALSE,ierr) ;CHKERRQ(ierr)

      call KSPSetTolerances(ksp,rtol,atol*(C%dof*C%glob_xm*C%glob_ym*C%glob_zm),PETSC_DEFAULT_REAL,maxiter,ierr);CHKERRQ(ierr)

      call KSPSetFromOptions(ksp,ierr) ;CHKERRQ(ierr)
      call KSPSetUp(ksp,ierr) ;CHKERRQ(ierr)

      call MatNullSpaceCreate( imp_comm, PETSC_TRUE, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, nullspace, ierr) ; CHKERRQ(ierr)
      call KSPSetNullspace(ksp, nullspace, ierr) ; CHKERRQ(ierr)

      linit = .True.
end subroutine

subroutine setup_logging()
        call PetscLogStageRegister('total_tenstream' , logstage(1)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('setup_edir'      , logstage(2)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('calc_edir'       , logstage(3)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('setup_ediff'     , logstage(4)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('calc_ediff'      , logstage(5)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('setup_b'         , logstage(6)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('get_coeff'       , logstage(7)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('twostream_dir'   , logstage(8)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('twostream_diff'  , logstage(9)     , ierr) ;CHKERRQ(ierr)
        call PetscLogStageRegister('write_hdf5'      , logstage(10)     , ierr) ;CHKERRQ(ierr)

        if(myid.eq.0) print *, 'Logging stages' , logstage
end subroutine

subroutine twostream(edirTOA)
    real(ireals),intent(in) :: edirTOA

    PetscReal,pointer,dimension(:,:,:,:) :: xv_dir,xv_diff,xv_abso
    integer(iintegers) :: i,j,k,src

    real(ireals),allocatable :: dtau(:),w0(:),g(:),S(:),Edn(:),Eup(:)
    real(ireals) :: mu0,Az

    integer(iintegers) :: li,lj

    call PetscLogStagePush(logstage(8),ierr) ;CHKERRQ(ierr)
    if(myid.eq.0) print *,' CALCULATING DELTA EDDINGTON TWOSTREAM'

    allocate( dtau(C_dir%zm-1) )
    allocate(   w0(C_dir%zm-1) )
    allocate(    g(C_dir%zm-1) )

    mu0 = sun%costheta

    call VecSet(edir_twostr ,zero,ierr); CHKERRQ(ierr)
    call VecSet(ediff_twostr,zero,ierr); CHKERRQ(ierr)

    call DMDAVecGetArrayF90(C_dir%da ,edir_twostr ,xv_dir ,ierr) ;CHKERRQ(ierr)
    call DMDAVecGetArrayF90(C_diff%da,ediff_twostr,xv_diff,ierr) ;CHKERRQ(ierr)

    allocate( S(C_dir%zm ) )
    allocate( Eup(C_diff%zm) )
    allocate( Edn(C_diff%zm) )

    Az = atm%dx*atm%dy
    do j=C_dir%ys,C_dir%ye         
        lj = jstartpar+j-C_dir%ys 

        do i=C_dir%xs,C_dir%xe
          li = istartpar+i-C_dir%xs 

          dtau = atm%dz * (atm%op(li,lj,:)%kabs + atm%op(li,lj,:)%ksca) 
          g    = atm%op(li,lj,:)%g
          w0   = atm%op(li,lj,:)%ksca / dtau
          call delta_eddington_twostream(dtau,w0,g,mu0,edirTOA,atm%albedo, S,Edn,Eup)

          do src=i0,i3
            xv_dir(src,i,j,:) = S(:) *.25_ireals *Az
          enddo

          xv_diff(E_up,i,j,:) = Eup(:) * Az
          xv_diff(E_dn,i,j,:) = Edn(:) * Az
        enddo
    enddo

    deallocate(S)
    deallocate(Edn)
    deallocate(Eup)


!    call DMDAVecGetArrayF90(C_one%da ,abso_twostr ,xv_abso ,ierr) ;CHKERRQ(ierr)
!
!    do k=C_one%zs,C_one%ze         
!      do j=C_one%ys,C_one%ye         
!        do i=C_one%xs,C_one%xe
!          xv_abso(i0,i,j,k) = xv_abso(i0,i,j,k) + sum( xv_dir (0:3 ,i,j,k   ) - xv_dir(0:3 ,i,j,k+i1) )*.25_ireals
!          xv_abso(i0,i,j,k) = xv_abso(i0,i,j,k) +      xv_diff(E_up,i,j,k+i1) + xv_diff(E_dn,i,j,k   )
!          xv_abso(i0,i,j,k) = xv_abso(i0,i,j,k) -      xv_diff(E_up,i,j,k   ) - xv_diff(E_dn,i,j,k+i1)
!        enddo
!      enddo
!    enddo
!
!    call DMDAVecRestoreArrayF90(C_one%da  ,abso_twostr ,xv_abso  ,ierr) ;CHKERRQ(ierr)

    call DMDAVecRestoreArrayF90(C_dir%da  ,edir_twostr ,xv_dir  ,ierr) ;CHKERRQ(ierr)
    call DMDAVecRestoreArrayF90(C_diff%da ,ediff_twostr,xv_diff ,ierr) ;CHKERRQ(ierr)

    call PetscLogStagePop(ierr) ;CHKERRQ(ierr)
end subroutine

subroutine scale_flx(v,C)
        Vec :: v
        type(coord) :: C
        PetscReal,pointer,dimension(:,:,:,:) :: xv
        PetscInt :: i,j,k,li,lj,lk
        PetscReal :: Ax,Ay,Az

        if(myid.eq.0.and.ldebug) print *,'rescaling fluxes'
        call DMDAVecGetArrayF90(C%da,v,xv,ierr) ;CHKERRQ(ierr)
        ! rescale total energy fluxes to average quantities i.e. W/m**2 or W/m**3

        Az = atm%dx*atm%dy

        do k=C%zs,C%ze-i1
          lk = i1+k-C%zs
          Ax = atm%dy*atm%dz(lk)
          Ay = atm%dx*atm%dz(lk)

          do j=C%ys,C%ye         
            lj = jstartpar+j-C%ys 

            do i=C%xs,C%xe      
              li = istartpar+i-C%xs 

              if(C%dof.eq.i8) then ! This is 8 stream direct radiation
                xv(i0:i3,i,j,k) = xv(i0:i3,i,j,k) / ( Az*.25 ) !*costheta0
                xv(i4:i5,i,j,k) = xv(i4:i5,i,j,k) / ( Ax*.5  ) !*sintheta0
                xv(i6:i7,i,j,k) = xv(i6:i7,i,j,k) / ( Ay*.5  ) !*sintheta0
              endif

              if(C%dof.eq.i10) then ! This is 10 stream diffuse radiation
                xv(E_up  ,i,j,k) = xv(E_up  ,i,j,k) / Az
                xv(E_dn  ,i,j,k) = xv(E_dn  ,i,j,k) / Az
                xv(E_le_m,i,j,k) = xv(E_le_m,i,j,k) / Ax
                xv(E_le_p,i,j,k) = xv(E_le_p,i,j,k) / Ax
                xv(E_ri_m,i,j,k) = xv(E_ri_m,i,j,k) / Ax
                xv(E_ri_p,i,j,k) = xv(E_ri_p,i,j,k) / Ax
                xv(E_ba_m,i,j,k) = xv(E_ba_m,i,j,k) / Ay
                xv(E_ba_p,i,j,k) = xv(E_ba_p,i,j,k) / Ay
                xv(E_fw_m,i,j,k) = xv(E_fw_m,i,j,k) / Ay
                xv(E_fw_p,i,j,k) = xv(E_fw_p,i,j,k) / Ay
              endif
            enddo
          enddo
        enddo

        k=C%ze
        do j=C%ys,C%ye         
          lj = jstartpar+j-C_diff%ys 

          do i=C%xs,C%xe      
            li = istartpar+i-C_diff%xs 

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

    subroutine init_tenstream(icomm, Nx,Ny,Nz, dx,dy,hhl ,phi0,theta0)
        integer,intent(in) :: icomm
        integer(iintegers),intent(in) :: Nx,Ny,Nz
        real(ireals),intent(in) :: dx,dy,hhl(:), phi0,theta0

        call init_mpi_data_parameters(icomm)

        atm%dx  = dx
        atm%dy  = dy

        allocate(atm%dz(size(hhl)-1))
        atm%dz = hhl( 1:size(hhl)-1 ) - hhl( 2:size(hhl) ) 

        allocate(atm%l1d( size(atm%dz) ) )
        where( twostr_ratio*atm%dz.gt.atm%dx ) 
          atm%l1d = .True.
        else where
          atm%l1d = .False.
        end where

        call setup_grid(Nx,Ny,Nz)
        call setup_suninfo(phi0,theta0,sun)

      ! init boxmc
        call OPP_8_10%init(atm%dx,atm%dy,[sun%symmetry_phi],[sun%theta],imp_comm)
        if(.not.luse_eddington) call OPP_1_2%init(atm%dx,atm%dy,[sun%symmetry_phi],[sun%theta],imp_comm) 

      ! init matrices & vectors
        call init_memory(incSolar,b,edir,ediff,abso,Mdir,Mdiff,edir_twostr,ediff_twostr,abso_twostr)

        call setup_logging()
    end subroutine

    subroutine set_optical_properties(global_kabs, global_ksca, global_g)
      real(ireals),intent(inout),dimension(:,:,:),allocatable :: global_kabs, global_ksca, global_g
      real(ireals) :: tau,w0
      integer(iintegers) :: i,j,k

      ! Make sure that our domain has at least 3 entries in each dimension.... otherwise violates boundary conditions
      call extend_arr(global_kabs)
      call extend_arr(global_ksca)
      call extend_arr(global_g)

      if(.not.allocated(atm%op) )       allocate( atm%op       (C_one%xm, C_one%ym, C_one%zm) )
      if(.not.allocated(atm%delta_op) ) allocate( atm%delta_op (C_one%xm, C_one%ym, C_one%zm) ) ! allocate space and copy optical properties for delta scaling

      if(luse_eddington) then
        if(.not.allocated(atm%a11) ) allocate(atm%a11 (C_one%xm, C_one%ym, C_one%zm )) ! allocate space for twostream coefficients
        if(.not.allocated(atm%a12) ) allocate(atm%a12 (C_one%xm, C_one%ym, C_one%zm )) 
        if(.not.allocated(atm%a13) ) allocate(atm%a13 (C_one%xm, C_one%ym, C_one%zm )) 
        if(.not.allocated(atm%a23) ) allocate(atm%a23 (C_one%xm, C_one%ym, C_one%zm )) 
        if(.not.allocated(atm%a33) ) allocate(atm%a33 (C_one%xm, C_one%ym, C_one%zm )) 
      endif
        
      ! Scatter global optical properties to MPI nodes
      call local_optprop(global_kabs, global_ksca, global_g)
      ! Now atm%op(:,:,:)%k... is populated.

      atm%delta_op = atm%op
      call delta_scale(atm%delta_op(:,:,:)%kabs, atm%delta_op(:,:,:)%ksca, atm%delta_op(:,:,:)%g )

      if(luse_eddington) then
        do k=1,size(atm%dz)
          if( atm%l1d(k) ) then
            do j=1,ubound(atm%op,2)
              do i=1,ubound(atm%op,1)
                tau = atm%dz(k)* (atm%op(i,j,k)%kabs + atm%op(i,j,k)%ksca)
                w0  = atm%op(i,j,k)%ksca / tau
                call eddington_coeff_rb ( tau , atm%op(i,j,k)%g, w0, sun%costheta, & 
                  atm%a11(i,j,k),          &
                  atm%a12(i,j,k),          &
                  atm%a13(i,j,k),          &
                  atm%a23(i,j,k),          &
                  atm%a33(i,j,k) )
              enddo
            enddo

          endif
        enddo
      endif

      if(myid.eq.0) then
        do k=1,ubound(atm%op,3)
          print *,myid,'Optical Properties:',k,'dz',atm%dz(k),atm%l1d(k),'k',minval(atm%op(:,:,k)%kabs),minval(atm%op(:,:,k)%ksca),minval(atm%op(:,:,k)%g),maxval(atm%op(:,:,k)%kabs),maxval(atm%op(:,:,k)%ksca),maxval(atm%op(:,:,k)%g)
        enddo
      endif

      contains
          subroutine local_optprop(global_kabs, global_ksca, global_g)
                real(ireals),allocatable,dimension(:,:,:) :: global_kabs, global_ksca, global_g

                !TODO : this is very poorly done... we should not scatter the global optical properties to all nodes and then pick what is local, rather only copy local parts....

                if(myid.eq.0.and.ldebug) print *,myid,'copying optprop: global to local :: shape kabs',shape(global_kabs),'xstart/end',C_one%xs,C_one%xe,'ys/e',C_one%ys,C_one%ye

                atm%op(:,:,:)%kabs = global_kabs(C_one%xs+1:C_one%xe+1, C_one%ys+1:C_one%ye+1, :)
                atm%op(:,:,:)%ksca = global_ksca(C_one%xs+1:C_one%xe+1, C_one%ys+1:C_one%ye+1, :)
                atm%op(:,:,:)%g    = global_g   (C_one%xs+1:C_one%xe+1, C_one%ys+1:C_one%ye+1, :)

                deallocate(global_kabs)
                deallocate(global_ksca)
                deallocate(global_g   )
          end subroutine

          subroutine extend_arr(arr)
                real(ireals),intent(inout),allocatable :: arr(:,:,:)
                real(ireals),allocatable :: tmp(:,:)
                integer(iintegers) :: dims(3),i
                integer(iintegers),parameter :: mindim=2

                dims = shape(arr)
                if( dims(1) .eq. 1 ) then
                  allocate( tmp(dims(2),dims(3) ), SOURCE=arr(1,:,:) )
                  deallocate(arr) 
                  allocate( arr(mindim, dims(2), dims(3) ) )
                  do i=1,mindim
                    arr(i, :, :) = tmp
                  enddo
                  deallocate(tmp)
                endif

                if( dims(2) .eq. 1 ) then
                  allocate( tmp(dims(1),dims(3) ), SOURCE=arr(:,1,:) )
                  deallocate(arr) 
                  allocate( arr(dims(1), mindim, dims(3) ) )
                  do i=1,mindim
                    arr(:,i , :) = tmp
                  enddo
                  deallocate(tmp)
                endif
          end subroutine
    end subroutine

    subroutine solve_tenstream(edirTOA)
      real(ireals),intent(in) :: edirTOA

            if(ltwostr) then
              call twostream(edirTOA)
              call calc_flx_div(edir_twostr,ediff_twostr,abso_twostr)

              call scale_flx(edir_twostr  ,C_dir)
              call scale_flx(ediff_twostr ,C_diff)

              if(luse_twostr_guess) then
                call VecCopy(edir_twostr  ,edir ,ierr) ;CHKERRQ(ierr)
                call VecCopy(ediff_twostr ,ediff,ierr) ;CHKERRQ(ierr)
              endif
              if(myid.eq.0) print *,'twostream calculation done'
            endif

            ! ---------------------------- Edir  -------------------
            call PetscLogStagePush(logstage(1),ierr) ;CHKERRQ(ierr)
            call setup_incSolar(incSolar,edirTOA)
            call set_dir_coeff(Mdir,C_dir)

            call setup_ksp(kspdir,C_dir,Mdir,linit_kspdir, "dir_")

            call PetscLogStagePush(logstage(3),ierr) ;CHKERRQ(ierr)
            call solve(kspdir,incSolar,edir)
            call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

            ! ---------------------------- Source Term -------------
            call setup_b(edir,b)

            ! ---------------------------- Ediff -------------------
            call set_diff_coeff(Mdiff,C_diff)

            call setup_ksp(kspdiff,C_diff,Mdiff,linit_kspdiff, "diff_")

            call PetscLogStagePush(logstage(5),ierr) ;CHKERRQ(ierr)
            call solve(kspdiff, b, ediff)
            call PetscLogStagePop(ierr) ;CHKERRQ(ierr)

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
    end subroutine
end module
