module petsc_ts
#include "finclude/petscdef.h"
      use petsc
      use gridtransform
      use arrayio
      use eddington, only : rodents
      USE data_parameters, ONLY :   &
          ireals,    &
          iintegers
      use kato_data
      use tenstream_optprop_8_10

      implicit none

      PetscInt,parameter :: E_up=0, E_dn=1, E_le_m=2, E_le_p=4, E_ri_m=3, E_ri_p=5, E_ba_m=6, E_ba_p=8, E_fw_m=7, E_fw_p=9
      PetscInt,parameter :: i0=0,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i10=10
      PetscInt,parameter :: istartpar=i1, jstartpar=i1
      PetscReal,parameter:: zero=0,one=1.0, pi=3.141592653589793, nil=-9999999

      logical,parameter :: l_writeall=.True.

      logical,parameter :: ldebug=.True.,lcycle_dir=.True.,luse_twostream=.False.
      character(len=*),parameter :: basepath='/home/users/jakub/scratch/tenstream/'

      logical,parameter ::ltest=.False.
!      character(len=*),parameter :: ident='run_test'   ; double precision,parameter :: ident_dx=67, ident_dy=67, ident_dz=100,phi0=270 ; logical,parameter ::ltest=.True.
!                                                                                                                              
!      character(len=*),parameter :: ident='run_box2'   ; double precision,parameter :: ident_dx=500, ident_dy=500, ident_dz=100,phi0=270
!      character(len=*),parameter :: ident='run_box4'   ; double precision,parameter :: ident_dx=70, ident_dy=70, ident_dz=100,phi0=270
!      character(len=*),parameter :: ident='run_box4'   ; double precision,parameter :: ident_dx=2800, ident_dy=2800, ident_dz=100,phi0=270
!      character(len=*),parameter :: ident='run_cb'     ; double precision,parameter :: ident_dx=250, ident_dy=250, ident_dz=100,phi0=180
      character(len=*),parameter :: ident='run_cosmo1' ; double precision,parameter :: ident_dx=2800,ident_dy=2800,ident_dz=250,phi0=180
!      character(len=*),parameter :: ident='run_cosmo2' ; double precision,parameter :: ident_dx=2800,ident_dy=2800,ident_dz=250,phi0=180
!      character(len=*),parameter :: ident='run_cosmo3' ; double precision,parameter :: ident_dx=2800,ident_dy=2800,ident_dz=250,phi0=180
!      character(len=*),parameter :: ident='run_i3rc1'  ; double precision,parameter :: ident_dx=66.7,ident_dy=66.7,ident_dz=100,phi0=180

      double precision,parameter :: albedo=0.05, theta0=0, twostr_ratio=.5_ireals !phi=azimuth ; theta=sza


      PetscInt :: myid,numnodes

      type coord
        PetscInt :: xs,xe,ys,ye,zs,ze     ! local domain start and end indices
        PetscInt :: xm,ym,zm              ! size of local domain
        PetscInt :: gxs,gys,gzs           ! domain indices including ghost points
        PetscInt :: gxm,gym,gzm           ! size of local domain including ghosts
        PetscInt :: glob_xm,glob_ym,glob_zm ! global domain size
        PetscInt :: dof,dim               ! degrees of freedom of Petsc Domain, dimension of dmda
        DM :: da                          ! The Domain Decomposition Object
        PetscInt,allocatable :: neighbors(:)         ! all 3d neighbours( (x=-1,y=-1,z=-1), (x=0,y=-1,z=-1) ...), i.e. 14 is one self.
      end type

      type(coord) :: C_dir,C_diff,C_one

      PetscErrorCode :: ierr
      integer(iintegers) :: iierr

      type t_optprop
        real(ireals) :: bg(3)=nil,fg(4)=nil
        logical :: updated=.True.
        double precision :: kext1=nil ,kext2=nil ,ksca1=nil ,ksca2=nil, w1=nil, w2=nil, g1=nil, g2=nil
      end type
      type(t_optprop),allocatable :: pcc_optprop(:,:,:)

      PetscReal :: phi
      PetscInt :: yinc,xinc

      PetscLogStage :: logstage(7)

      contains 
     subroutine chkerr(errcode,str)
      PetscErrorCode,intent(in) :: errcode
      character(len=*),intent(in),optional :: str
      if(errcode.ne.0) then
              if(present(str)) then
                      print *,'******** ERROR :: ',errcode,':: ',trim(str)
                      call MPI_Abort(PETSC_COMM_WORLD,errcode,ierr)
              else
                      print *,'******** ERROR :: ',errcode
                      call MPI_Abort(PETSC_COMM_WORLD,errcode,ierr)
              endif
              call petscfinalize(ierr)
              call exit()
       endif
       end subroutine

elemental double precision function deg2rad(a)
                double precision,intent(in) :: a
                deg2rad = a*pi/180.
        end function

      subroutine setup_grid()
        PetscInt,parameter :: stencil_size=1
        integer :: mpierr
        DMBoundaryType :: bp=DM_BOUNDARY_PERIODIC, bn=DM_BOUNDARY_NONE, bg=DM_BOUNDARY_GHOSTED
        DMDAStencilType,parameter  :: stype=DMDA_STENCIL_BOX

        call MPI_Comm_size( PETSC_COMM_WORLD, numnodes, mpierr)
        if(myid.eq.0.and.ldebug) print *,myid,'Setting up the DMDA grid for ',newgrid%Nx,newgrid%Ny,newgrid%Nz,'using ',numnodes,' nodes'
        
        if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_diff'
        C_diff%dof=10
        call DMDACreate3d(PETSC_COMM_WORLD,bp,bp,bn,stype,       &
                              -newgrid%Nx, -newgrid%Ny, -(newgrid%Nz+1),  &
                              PETSC_DECIDE,PETSC_DECIDE,i1,    &
                              i1*C_diff%dof,stencil_size,              &
                              PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
                              C_diff%da,ierr)
        call setup_coords(C_diff)

        if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C_dir'
        C_dir%dof=8
        if(lcycle_dir) then
          call DMDACreate3d(PETSC_COMM_WORLD,bp,bp,bn,stype,       &
              -newgrid%Nx, -newgrid%Ny, -(newgrid%Nz+1),  &
              PETSC_DECIDE,PETSC_DECIDE,i1,    &
              i1*C_dir%dof,stencil_size,              &
              PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
              C_dir%da,ierr)
        else
          call DMDACreate3d(PETSC_COMM_WORLD,bg,bg,bn,stype,       &
              -newgrid%Nx, -newgrid%Ny, -(newgrid%Nz+1),  &
              PETSC_DECIDE,PETSC_DECIDE,i1,    &
              i1*C_dir%dof,stencil_size,              &
              PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
              C_dir%da,ierr)
        endif
        call setup_coords(C_dir)

        if(myid.eq.0.and.ldebug) print *,myid,'Configuring DMDA C1'
        C_one%dof=1
        call DMDACreate3d(PETSC_COMM_WORLD,bp,bp,bn,stype,       &
                              -newgrid%Nx, -newgrid%Ny, -(newgrid%Nz),  &
                              PETSC_DECIDE,PETSC_DECIDE,i1,    &
                              i1*C_one%dof,stencil_size,              &
                              PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
                              C_one%da,ierr)
        call setup_coords(C_one)

        if(myid.eq.0.and.ldebug) print *,myid,'DMDA grid ready'
        contains
        subroutine setup_coords(C)
          type(coord) :: C

          call DMDAGetInfo(C%da,C%dim,                             &
            C%glob_xm,C%glob_ym,C%glob_zm,                           &
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
            PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
            ierr)

          call DMDAGetCorners(C%da,C%xs,C%ys,C%zs, C%xm,C%ym,C%zm, ierr)
          call DMDAGetGhostCorners(C%da,C%gxs,C%gys,C%gzs,C%gxm,C%gym,C%gzm,ierr)
          C%xe = C%xs+C%xm-1
          C%ye = C%ys+C%ym-1
          C%ze = C%zs+C%zm-1
          if(ldebug) print *,myid,'Domain Corners x:: ',C%xs,':',C%xe,' (',C%xm,' entries)','global size',C%glob_xm
          if(ldebug) print *,myid,'Domain Corners y:: ',C%ys,':',C%ye,' (',C%ym,' entries)','global size',C%glob_ym
          if(ldebug) print *,myid,'Domain Corners z:: ',C%zs,':',C%ze,' (',C%zm,' entries)','global size',C%glob_zm

          allocate(C%neighbors(0:3**C%dim-1) )
          call DMDAGetNeighbors(C%da,C%neighbors,ierr)
          call DMSetUp(C%da,ierr)
          if(ldebug.and.C%dim.eq.3) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([10,12,16,14]),'while I am ',C%neighbors(13)
          if(ldebug.and.C%dim.eq.2) print *,'PETSC id',myid,C%dim,'Neighbors are',C%neighbors([1,3,7,5]),'while I am ',C%neighbors(4)
        end subroutine
      end subroutine

      subroutine mat_info(A)
        Mat :: A
        double  precision info(MAT_INFO_SIZE)
        double  precision mal, nz_allocated, nz_used, nz_unneeded

        call MatGetInfo(A,MAT_LOCAL,info,ierr)
        mal = info(MAT_INFO_MALLOCS)
        nz_allocated = info(MAT_INFO_NZ_ALLOCATED)
        nz_used   = info(MAT_INFO_NZ_USED)
        nz_unneeded = info(MAT_INFO_NZ_UNNEEDED)

        if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_MALLOCS',mal,'MAT_INFO_NZ_ALLOCATED',nz_allocated
        if(myid.eq.0.and.ldebug) print *,myid,'mat_info :: MAT_INFO_USED',nz_used,'MAT_INFO_NZ_unneded',nz_unneeded

      end subroutine
      subroutine init_matrix(A,C)
        Mat :: A
        type(coord) :: C

!        PetscInt :: Msize,Msize_glob
        PetscInt,dimension(:),allocatable :: o_nnz,d_nnz!,dnz

        call DMSetMatrixPreallocateOnly(C%da, PETSC_TRUE,ierr) ; call chkerr(ierr,'Init Matrix :: DMSetMatrixPreallocateOnly')
!        call MatCreate(PETSC_COMM_WORLD, A, ierr) ; call chkerr(ierr,'Init Matrix :: MatCreate')
!
!        Msize =  C%xm*C%ym*C%zm*C%dof 
!        Msize_glob = C%glob_xm*C%glob_ym*C%glob_zm*C%dof
!
!        call MatSetSizes(A, Msize, Msize, Msize_glob, Msize_glob, ierr) ; call chkerr(ierr,'Init Matrix :: MatSetSizes')

        call DMCreateMatrix(C%da, A, ierr) ; call chkerr(ierr,'Init Matrix :: CreateMatrix')

        call MatSetFromOptions(A,ierr) ; call chkerr(ierr,'Init Matrix :: Set Options')

        ! Determine perfect preallocation
        if(myid.eq.0.and.ldebug) print *,'determine preallocation for ',C%dof
        select case (C%dof)
                case(i3)
                        call setup_dir_preallocation(d_nnz,o_nnz,C)
                        call MatMPIAIJSetPreallocation(A, PETSC_NULL_INTEGER,d_nnz, PETSC_NULL_INTEGER, o_nnz, ierr) ; call chkerr(ierr,'Init Matrix :: MatMPIAIJSetPreallocation dir')
                case(i8)
                        call setup_dir8_preallocation(d_nnz,o_nnz,C)
                        call MatMPIAIJSetPreallocation(A, PETSC_NULL_INTEGER,d_nnz, PETSC_NULL_INTEGER, o_nnz, ierr) ; call chkerr(ierr,'Init Matrix :: MatMPIAIJSetPreallocation dir 8streams')
                case(i10)
                        call setup_diff_preallocation(d_nnz,o_nnz,C)
                        call MatMPIAIJSetPreallocation(A, PETSC_NULL_INTEGER,d_nnz, PETSC_NULL_INTEGER, o_nnz, ierr) ; call chkerr(ierr,'Init Matrix :: MatMPIAIJSetPreallocation diff')
                case default
                        ierr = C%dof
                        call chkerr(ierr,'Dont know which preallocation routine I shall call! - exiting...')
        end select

        call MatSeqAIJSetPreallocation(A, C%dof+i1, PETSC_NULL_INTEGER, ierr) ; call chkerr(ierr,'Init Matrix :: MatSeqAIJSetPreallocation')

        deallocate(o_nnz)
        deallocate(d_nnz)

        !call MatMPIAIJSetPreallocation(A, C%dof+1, PETSC_NULL_INTEGER, C%dof, PETSC_NULL_INTEGER, ierr)
        !call MatSeqAIJSetPreallocation(A, C%dof+1, PETSC_NULL_INTEGER, ierr) ; call chkerr(ierr,'Init Matrix :: MatSeqAIJSetPreallocation')

        call MatSetOption(A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr)

        call mat_set_diagonal(A,C)
      end subroutine
      subroutine mat_set_diagonal(A,C)
        Mat :: A
        type(coord),intent(in) :: C
        PetscInt :: i,j,k,dof
        MatStencil :: row(4,1), col(4,1)
        PetscReal :: v(1)

        character(100) :: errstr

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
                                        
                                        call MatSetValuesStencil(A,i1, row,i1, col , v ,INSERT_VALUES,ierr) 
                                        if(ldebug) then
                                                write(errstr,FMT='("mat_set_diagonal :: MatSetValuesStencil",I0,".",I0,".",I0,".",I0)' ) dof,i,j,k
                                                call chkerr(ierr,errstr)
                                        endif
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

        PetscInt :: vsize !,mrows,mcols, i,j,k,dof,li,lj,lk

        PetscScalar, pointer :: xx_v_o(:),xx_v_d(:)

        PetscInt,parameter :: ind(9)=[E_up,E_le_m,E_le_p,E_ri_m,E_ri_p,E_ba_m,E_ba_p,E_fw_m,E_fw_p]

        call DMCreateGlobalVector(C%da,v_o_nnz,ierr)
        call DMCreateGlobalVector(C%da,v_d_nnz,ierr)
        call DMDAVecGetArrayF90(C%da,v_o_nnz,xo,ierr)
        call DMDAVecGetArrayF90(C%da,v_d_nnz,xd,ierr)

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

!        if(luse_twostream) then
!          do k=C%zs,C%ze-1
!            do j=C%ys,C%ye
!              do i=C%xs,C%xe        ! i,j,k indices are defined on petsc global grid
!                li = istartpar+i-C%xs ! l-indices are used for local arrays, 
!                lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
!                lk = 1+k-C%zs         !
!                if( newgrid%dz(lk).gt.newgrid%dx(li) ) then
!                  xo(:,i,j,k) = i0
!                  xd(:,i,j,k) = i1
!                  xd([E_up,E_dn],i,j,k) = i3
!                endif
!              enddo
!            enddo
!          enddo
!        endif

        call DMDAVecRestoreArrayF90(C%da,v_o_nnz,xo,ierr)
        call DMDAVecRestoreArrayF90(C%da,v_d_nnz,xd,ierr)
        
        call VecGetLocalSize(v_o_nnz,vsize,ierr)
        allocate(o_nnz(0:vsize-1))
        allocate(d_nnz(0:vsize-1))

        call VecGetArrayF90(v_o_nnz,xx_v_o,ierr)
        call VecGetArrayF90(v_d_nnz,xx_v_d,ierr)
        o_nnz=int(xx_v_o)
        d_nnz=int(xx_v_d)
        call VecRestoreArrayF90(v_o_nnz,xx_v_o,ierr)
        call VecRestoreArrayF90(v_d_nnz,xx_v_d,ierr)

        call VecDestroy(v_o_nnz,ierr)
        call VecDestroy(v_d_nnz,ierr)

        if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
      end subroutine 
      subroutine setup_dir8_preallocation(d_nnz,o_nnz,C)
        PetscInt,allocatable :: d_nnz(:)
        PetscInt,allocatable :: o_nnz(:)
        type(coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize,i,j !mrows,mcols, i,j,k,li,lj,lk

        PetscScalar, pointer :: xx_v_o(:),xx_v_d(:)

        logical :: lsun_east,lsun_north

!        if(myid.eq.0.and.ldebug) print *,myid,'building direct o_nnz for mat with',C%dof,'dof'
        call DMCreateGlobalVector(C%da,v_o_nnz,ierr)
        call DMCreateGlobalVector(C%da,v_d_nnz,ierr)
        call DMDAVecGetArrayF90(C%da,v_o_nnz,xo,ierr)
        call DMDAVecGetArrayF90(C%da,v_d_nnz,xd,ierr)

        xo = i0
        xd = i1
        xd( i0:i3 ,:,:,C%zs+1 :C%ze   ) = C%dof+i1 ! Edir_vertical depends on 3 values Edir_vertical,xaxis,yaxis :: starting with second entries(seen from top)
        xd( i4:i7 ,:,:,C%zs   :C%ze-i1) = C%dof+i1 ! Edir_xaxis,yaxis depends on 3 values Edir_vertical,xaxis,yaxis :: starting with first entries(seen from top)

        do j=C%ys,C%ye
                 lsun_east  = (xinc.eq.i0)

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
                 lsun_north = (yinc.eq.i0 )

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
                 lsun_east  = (xinc.eq.i0)

                 if( C%neighbors(12).ne.myid.and. C%neighbors(12).ge.i0 ) then ! neigh west
                         if( .not. lsun_east ) then 
                                ! if the sun is in the west, the 2nd channel is solemnly dependant on ghost values
                                xo(i4:i5, C%xs, j, C%zs:C%ze-1) = C%dof
                                xd(i4:i5, C%xs, j, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo
        do i=C%xs,C%xe
                 lsun_north = (yinc.eq.i0 )

                 if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                         if( .not. lsun_north ) then 
                                ! if the sun is in the south, the 3rd channel is solemnly dependant on ghost values
                                xo(i6:i7, i, C%ys, C%zs:C%ze-1) = C%dof
                                xd(i6:i7, i, C%ys, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo

!        if(luse_twostream) then
!          do k=C%zs,C%ze-1
!            do j=C%ys,C%ye
!              do i=C%xs,C%xe        ! i,j,k indices are defined on petsc global grid
!                li = istartpar+i-C%xs ! l-indices are used for local arrays, 
!                lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
!                lk = 1+k-C%zs         !
!                if( newgrid%dz(lk).gt.newgrid%dx(li) ) then
!                  xo(:,i,j,k) = i0
!                  xd(:,i,j,k) = i1
!                  xd(0:3,i,j,k) = i2
!                endif
!              enddo
!            enddo
!          enddo
!        endif

        call DMDAVecRestoreArrayF90(C%da,v_o_nnz,xo,ierr)
        call DMDAVecRestoreArrayF90(C%da,v_d_nnz,xd,ierr)
        
        call VecGetLocalSize(v_o_nnz,vsize,ierr)
        allocate(o_nnz(0:vsize-1))
        allocate(d_nnz(0:vsize-1))

        call VecGetArrayF90(v_o_nnz,xx_v_o,ierr)
        call VecGetArrayF90(v_d_nnz,xx_v_d,ierr)
        o_nnz=int(xx_v_o)
        d_nnz=int(xx_v_d)
        call VecRestoreArrayF90(v_o_nnz,xx_v_o,ierr)
        call VecRestoreArrayF90(v_d_nnz,xx_v_d,ierr)

        call VecDestroy(v_o_nnz,ierr)
        call VecDestroy(v_d_nnz,ierr)

        if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
      end subroutine 
      subroutine setup_dir_preallocation(d_nnz,o_nnz,C)
        PetscInt,allocatable :: d_nnz(:)
        PetscInt,allocatable :: o_nnz(:)
        type(coord) :: C
        Vec :: v_o_nnz,v_d_nnz
        PetscScalar,Pointer :: xo(:,:,:,:),xd(:,:,:,:)

        PetscInt :: vsize,i,j !,mrows,mcols, i,j,k,li,lj,lk

        PetscScalar, pointer :: xx_v_o(:),xx_v_d(:)

        logical :: lsun_east,lsun_north

!        if(myid.eq.0.and.ldebug) print *,myid,'building direct o_nnz for mat with',C%dof,'dof'
        call DMCreateGlobalVector(C%da,v_o_nnz,ierr)
        call DMCreateGlobalVector(C%da,v_d_nnz,ierr)
        call DMDAVecGetArrayF90(C%da,v_o_nnz,xo,ierr)
        call DMDAVecGetArrayF90(C%da,v_d_nnz,xd,ierr)

        xo = i0
        xd = i1
        xd( i0    ,:,:,C%zs+1:C%ze   ) = C%dof+i1 ! Edir_vertical depends on 3 values Edir_vertical,xaxis,yaxis :: starting with second entries(seen from top)
        xd([i1,i2],:,:,C%zs   :C%ze-i1) = C%dof+i1 ! Edir_xaxis,yaxis depends on 3 values Edir_vertical,xaxis,yaxis :: starting with first entries(seen from top)

        do j=C%ys,C%ye
                 lsun_east  = (xinc.eq.i0)

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
                 lsun_north = (yinc.eq.i0 )

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
                 lsun_east  = (xinc.eq.i0)

                 if( C%neighbors(12).ne.myid.and. C%neighbors(12).ge.i0 ) then ! neigh west
                         if( .not. lsun_east ) then 
                                ! if the sun is in the west, the 2nd channel is solemnly dependant on ghost values
                                xo(i1, C%xs, j, C%zs:C%ze-1) = i3
                                xd(i1, C%xs, j, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo
        do i=C%xs,C%xe
                 lsun_north = (yinc.eq.i0 )

                 if( C%neighbors(10).ne.myid .and. C%neighbors(10).ge.i0 ) then ! neigh south
                         if( .not. lsun_north ) then 
                                ! if the sun is in the south, the 3rd channel is solemnly dependant on ghost values
                                xo(i2, i, C%ys, C%zs:C%ze-1) = i3
                                xd(i2, i, C%ys, C%zs:C%ze-1) = i1
                        endif
                endif
        enddo

!        if(luse_twostream) then
!          do k=C%zs,C%ze-1
!            do j=C%ys,C%ye
!              do i=C%xs,C%xe        ! i,j,k indices are defined on petsc global grid
!                li = istartpar+i-C%xs ! l-indices are used for local arrays, 
!                lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
!                lk = 1+k-C%zs         !
!                if( newgrid%dz(lk).gt.newgrid%dx(li) ) then
!                  xo(:,i,j,k) = i0
!                  xd(:,i,j,k) = i1
!                  xd(0,i,j,k) = i2
!                endif
!              enddo
!            enddo
!          enddo
!        endif

        call DMDAVecRestoreArrayF90(C%da,v_o_nnz,xo,ierr)
        call DMDAVecRestoreArrayF90(C%da,v_d_nnz,xd,ierr)
        
        call VecGetLocalSize(v_o_nnz,vsize,ierr)
        allocate(o_nnz(0:vsize-1))
        allocate(d_nnz(0:vsize-1))

        call VecGetArrayF90(v_o_nnz,xx_v_o,ierr)
        call VecGetArrayF90(v_d_nnz,xx_v_d,ierr)
        o_nnz=int(xx_v_o)
        d_nnz=int(xx_v_d)
        call VecRestoreArrayF90(v_o_nnz,xx_v_o,ierr)
        call VecRestoreArrayF90(v_d_nnz,xx_v_d,ierr)

        call VecDestroy(v_o_nnz,ierr)
        call VecDestroy(v_d_nnz,ierr)

        if(myid.eq.0 .and. ldebug) print *,myid,'direct d_nnz, ',sum(d_nnz),'o_nnz',sum(o_nnz),'together:',sum(d_nnz)+sum(o_nnz),'expected less than',vsize*(C%dof+1)
      end subroutine 

        subroutine get_coeff(op,dz,dir,coeff,angles)
               double precision,intent(out) :: coeff(:)
               type(t_optprop),intent(in) :: op
               double precision,intent(in) :: dz
               double precision,intent(in),optional :: angles(2)
               logical,intent(in) :: dir
               logical,parameter :: lround=.False.
               double precision,parameter :: fround=1e-4

               double precision :: kabs,ksca,g !,tmp_coeff(size(coeff))

               ! coeff mode determines the coeff by:
               ! 1: plane parallel approx
               ! 2: ipa
               ! 3: mix of ipa and plane-parallel
               ! 4: neural network with cloud fraction

!               integer(iintegers),parameter :: coeff_mode=2
               call PetscLogStagePush(logstage(7),ierr)

               kabs = op%bg(1) 
               ksca = op%bg(2) 
               g    = op%bg(3)
               call optprop_lookup_coeff(dz,kabs,ksca,g,dir,coeff,angles)

!               select case(coeff_mode)
!               case(1) ! plane parallel
!                       kext = op%kext1*(one-op%fg(4)) + op%kext2*op%fg(4)
!                       ksca = op%ksca1*(one-op%fg(4)) + op%ksca2*op%fg(4)
!                       w   = ksca/kext
!                       g   = (op%g1*op%ksca1 + op%g2*op%ksca2)/(op%ksca1+op%ksca2)
!
!                       call optprop_lookup_coeff(dz,kext*dz,w,g,dir,coeff,angles)
!
!               case(2) ! IPA
!                       call optprop_lookup_coeff(dz,op%kext1*dz,op%w1,op%g1,dir,coeff,angles)
!                       if(op%fg(4).ge.1e-3) then
!                         call optprop_lookup_coeff(dz,op%kext2*dz,op%w2,op%g2,dir,tmp_coeff    ,angles)
!                         coeff = coeff*(one-op%fg(4)) + tmp_coeff*op%fg(4)
!                       endif
!
!               case(3) ! mean( pp.,  IPA)
!                       call optprop_lookup_coeff(dz,op%kext1*dz,op%w1,op%g1,dir,tmp_coeff,angles)
!                       call optprop_lookup_coeff(dz,op%kext2*dz,op%w2,op%g2,dir,coeff    ,angles)
!
!                       coeff = coeff*(one-op%fg(4)) + tmp_coeff*op%fg(4)
!
!                       kext = op%kext1*(one-op%fg(4)) + op%kext2*op%fg(4)
!                       ksca = op%ksca1*(one-op%fg(4)) + op%ksca2*op%fg(4)
!                       w   = ksca/kext
!                       g   = (op%g1*op%ksca1 + op%g2*op%ksca2)/(op%ksca1+op%ksca2)
!                       call optprop_lookup_coeff(dz,kext*dz,w,g,dir,tmp_coeff,angles)
!
!                       coeff = (coeff+tmp_coeff)/2
!
!               case(4)
!                       print *,'ANN currently not implemented! Need to determine coefficients by other means, please change coeff_mode!'
!                       call exit()
!               case default
!                       print *,'coeff mode:',coeff_mode,' not supported??'
!                       call exit()
!               end select
               if(lround) coeff = max(zero, min(one, anint(coeff/fround) * fround ))
!               if(op%kext1.ge.1e-2.and.present(angles).and.dir) print *,'tod',op%kext1*dz,'coeff',coeff
               if(ldebug) then
                 if( any(isnan(coeff)) .or. any(coeff.lt.zero) .or. any(coeff.gt.one) ) print *,'Wrong coeff',coeff,'op',op
               endif
            call PetscLogStagePop(ierr)
        end subroutine

subroutine setup_dir_inc(phi0)
  double precision,intent(in) :: phi0
  ! use symmetry for direct beam: always use azimuth [0,90] an just reverse the order where we insert the coeffs
  ! ''swap'' phi axis down to the range of [0,180] 
  phi = acos(cos(phi0*pi/180))
  !        print *,'1st phi0 swap',phi0,phi,'=',phi0*pi/180,cos(phi0*pi/180),acos(cos(phi0*pi/180))
  ! and then mirror it onto range [0,90]
  phi = asin(sin(phi)) /pi * 180
  !        print *,'2nd phi0 swap',phi0,phi,'=',sin(phi),asin(sin(phi)),asin(sin(phi)) /pi * 180,int(asin(sin(phi)) /pi * 180)
  xinc=i0 ; yinc=i0
  if(phi0.gt.180) xinc=i1
  if(phi0.gt.90.and.phi0.lt.270) yinc=i1
end subroutine

subroutine set_dir_coeff(A,C)
        Mat :: A
        type(coord) :: C

        PetscInt :: i,j,k,src,dst, li,lj,lk

        MatStencil :: row(4,C%dof)  ,col(4,C%dof)
        PetscReal :: v(C%dof**2),coeffs(C%dof**2),norm
        PetscInt,parameter :: entries(8)=[0,8,16,24,32,40,48,56]
        PetscReal :: edd_coeff(5)

!        if(ldebug) print *,myid,'DEBUG(set_dir_coeff) DMDA',C%xs,C%xe,C%ys,C%ye,C%zs,C%ze-1

        row=PETSC_NULL_INTEGER
        col=PETSC_NULL_INTEGER
        v=PETSC_NULL_REAL

        do k=C%zs,C%ze-1
        do j=C%ys,C%ye
        do i=C%xs,C%xe        ! i,j,k indices are defined on petsc global grid
          li = istartpar+i-C%xs ! l-indices are used for local arrays, 
          lj = jstartpar+j-C%ys ! which are defined on the cosmo grid
          lk = 1+k-C%zs         !
          !          print *,'setting direct coeff',li,lj,lk,'PETSC_grid',i,j,k,'pcc_optprop',lbound(pcc_optprop,1),ubound(pcc_optprop,1),':',lbound(pcc_optprop,2),ubound(pcc_optprop,2),':',lbound(pcc_optprop,3),ubound(pcc_optprop,3)

          !         _______________
          !        |               |
          !        |           3   |  
          !      2 |               |  
          !        |               |
          !        |               |
          !        |_______________|  
          !                 1

          dst = 1 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 2 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 3 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 4 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k+1     ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the lower/upper lid
          dst = 5 ; row(MatStencil_i,dst) = i+xinc ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
          dst = 6 ; row(MatStencil_i,dst) = i+xinc ; row(MatStencil_j,dst) = j       ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the left/right lid
          dst = 7 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j+yinc  ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid
          dst = 8 ; row(MatStencil_i,dst) = i      ; row(MatStencil_j,dst) = j+yinc  ;  row(MatStencil_k,dst) = k       ; row(MatStencil_c,dst) = dst-i1 ! Define transmission towards the front/back lid

          src = 1 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 2 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 3 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 4 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the upper/lower lid:
          src = 5 ; col(MatStencil_i,src) = i+1-xinc ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
          src = 6 ; col(MatStencil_i,src) = i+1-xinc ; col(MatStencil_j,src) = j        ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the left/right lid
          src = 7 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j+1-yinc ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:
          src = 8 ; col(MatStencil_i,src) = i        ; col(MatStencil_j,src) = j+1-yinc ; col(MatStencil_k,src) = k   ; col(MatStencil_c,src) = src-i1 ! Source may be the front/back lid:


          !       if(myid.eq.0.and.li.eq.istartpar.and.lj.eq.jstartpar.and.lk.eq.1) optprop_debug=.True.
          if(luse_twostream .and. twostr_ratio*newgrid%dz(lk).gt.newgrid%dx(li) ) then
            call rodents(pcc_optprop(li,lj,lk)%kext1*newgrid%dz(lk), pcc_optprop(li,lj,lk)%w1, pcc_optprop(li,lj,lk)%g1, cos(deg2rad(theta0)), edd_coeff )
            coeffs=zero
            coeffs([0, 8+1, 16+2, 24+3 ]+i1) = edd_coeff(5) ! only use the four vertical tiles with a33
          else ! use only one dimensional direct radiation
            call get_coeff(pcc_optprop(li,lj,lk), newgrid%dz(lk),.True., coeffs, [phi, theta0])
          endif

          ! make sure that energy is conserved... better to be on the absorbing side
          if(ldebug) then
            do src=1,C%dof
              norm = sum( coeffs((src-1)*C%dof+1:src*C%dof) )
              !          if( norm.ge.one )  coeffs((src-1)*C%dof+1:src*C%dof)  = coeffs((src-1)*C%dof+1:src*C%dof) / (norm+1e-8)
              if( ldebug .and. norm.gt.one ) print *,'sum(src==',src,') ge one',norm
            enddo
          endif

          ! reorder coeffs from src-ordered to dst-ordered
          do src=1,C%dof
            v(entries+src) = coeffs( i1+(src-i1)*C%dof : src*C%dof )
          enddo

          call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr)

        enddo ; enddo ; enddo

        if(myid.eq.0.and.ldebug) print *,myid,'setup_direct_matrix done'

        call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr) ; call chkerr(ierr,'Set direct coeff :: MatAssemblyBegin')
        call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr) ; call chkerr(ierr,'Set direct coeff :: MatAssemblyEnd')

        end subroutine

        subroutine setup_incSolar(incSolar,kato,iq)
          Vec :: incSolar
          integer(iintegers) :: kato,iq

          PetscInt :: i,j,li,lj
          PetscScalar,pointer,dimension(:,:,:,:) :: xincSolar

          PetscReal :: edirTOA

          call VecSet(incSolar,zero,ierr)

          edirTOA = get_edirTOA(kato,iq,newgrid%z(1))

          call DMDAVecGetArrayF90(C_dir%da,incSolar,xincSolar,ierr)
          do j=C_dir%ys,C_dir%ye
            do i=C_dir%xs,C_dir%xe        ! i,j,k indices are defined on petsc global grid
              li = istartpar+i-C_dir%xs ! l-indices are used for local arrays, 
              lj = jstartpar+j-C_dir%ys ! which are defined on the cosmo grid

              xincSolar(i0:i3,i,j,C_dir%zs) = edirTOA* newgrid%dx(li)*newgrid%dy(lj) * max(zero,cos(deg2rad(theta0))) *.25
            enddo
          enddo

          call DMDAVecRestoreArrayF90(C_dir%da,incSolar,xincSolar,ierr)
          if(myid.eq.0.and.ldebug) print *,myid,'Setup of IncSolar done'

        end subroutine

subroutine set_diff_coeff(A,C)
  Mat :: A
  type(coord) :: C

  PetscInt :: i,j,k,src,dst, li,lj,lk

  MatStencil :: row(4,0:C%dof-1)  ,col(4,0:C%dof-1)
  PetscReal :: v(C%dof**2),coeffs(C%dof**2),norm,edd_coeff(5)
  PetscInt,parameter :: entries(10)=[0,10,20,30,40,50,60,70,80,90]

  ! if(ldebug) print *,myid,'DEBUG(set_dir_coeff) jspec',jspec,'igas',igas,'isub',isub

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

        if(luse_twostream .and. twostr_ratio*newgrid%dz(lk).gt.newgrid%dx(li) ) then
          coeffs = zero
          call rodents(pcc_optprop(li,lj,lk)%kext1*newgrid%dz(lk), pcc_optprop(li,lj,lk)%w1, pcc_optprop(li,lj,lk)%g1, cos(deg2rad(theta0)), edd_coeff )
          coeffs([ 0,1, 10,11 ]+i1) = [ edd_coeff(2), edd_coeff(1), edd_coeff(1), edd_coeff(2) ]
        else
          call get_coeff(pcc_optprop(li,lj,lk), newgrid%dz(lk),.False., coeffs )
        endif

        if(ldebug) then
          do dst=1,C%dof
            norm = sum( coeffs((dst-1)*C%dof+1:dst*C%dof) )
            if( norm.ge.one )  coeffs((dst-1)*C%dof+1:dst*C%dof)  = coeffs((dst-1)*C%dof+1:dst*C%dof) / (norm+1e-8_ireals)
            if( ldebug .and. norm.gt.one ) print *,'diffuse sum(dst==',dst,') ge one',norm
          enddo
        endif

        ! reorder coeffs from src-ordered to dst-ordered
        do src=1,C%dof
          v(entries+src) = coeffs( i1+(src-i1)*C%dof : src*C%dof )
        enddo

        call MatSetValuesStencil(A,C%dof, row,C%dof, col , -v ,INSERT_VALUES,ierr)

      enddo ; enddo ; enddo
      if(myid.eq.0.and.ldebug) print *,myid,'setup_diffuse_matrix done'

      if(myid.eq.0.and.ldebug) print *,myid,'Final diffuse Matrix Assembly:'
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)
      call MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE, ierr)
end subroutine
 
subroutine setup_b(edir,kato,iq,b)
        Vec :: edir
        Vec :: local_b,b
        integer(iintegers),intent(in) :: kato,iq

        PetscScalar,pointer,dimension(:,:,:,:) :: &
        xsrc,xedir
        PetscInt :: li,lj,lk
        PetscReal :: coeffs(C_dir%dof*C_diff%dof),edd_coeff(5)
        
        PetscInt :: i,j,k,src
        
        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly...'

        call VecSet(b,zero,ierr)

        call DMCreateLocalVector(C_diff%da,local_b,ierr)
        call VecSet(local_b,zero,ierr)

        call DMDAVecGetArrayF90(C_diff%da,local_b,xsrc,ierr)
        call DMDAVecGetArrayF90(C_dir%da,edir,xedir,ierr)

        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients'
        do k=C_diff%zs,C_diff%ze-1 
          do j=C_diff%ys,C_diff%ye         
            do i=C_diff%xs,C_diff%xe    
              li = istartpar+i-C_diff%xs
              lj = jstartpar+j-C_diff%ys
              lk = i1+k-C_diff%zs        

              if(luse_twostream .and. twostr_ratio*newgrid%dz(lk).gt.newgrid%dx(li) ) then
                call rodents(pcc_optprop(li,lj,lk)%kext1*newgrid%dz(lk), pcc_optprop(li,lj,lk)%w1, pcc_optprop(li,lj,lk)%g1, cos(deg2rad(theta0)), edd_coeff )
                xsrc(E_up   ,i,j,k)   = xsrc(E_up   ,i,j,k)   +  xedir(i0,i,j,k)*edd_coeff(3)
                xsrc(E_dn   ,i,j,k+1) = xsrc(E_dn   ,i,j,k+1) +  xedir(i0,i,j,k)*edd_coeff(4)

              else
                call get_coeff(pcc_optprop(li,lj,lk), newgrid%dz(lk),.False., coeffs, [phi, theta0] )

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

        do i=C_diff%xs,C_diff%xe     
          do j=C_diff%ys,C_diff%ye     
            li = istartpar+i-C_diff%xs 
            lj = jstartpar+j-C_diff%ys 
            k = C_diff%ze
            !          xsrc(   :   ,i,j,k) = nil
            xsrc(E_up   ,i,j,k) = sum(xedir(i0:i3,i,j,k))*albedo
          enddo
        enddo
        if(newgrid%z(1).le.15e3) then
          if(myid.eq.0) print *,'Updating diffuse downward stream with ',get_ednTOA(kato,iq,newgrid%z(1)),'W/m2'
          do i=C_diff%xs,C_diff%xe     
            do j=C_diff%ys,C_diff%ye     
              li = istartpar+i-C_diff%xs 
              lj = jstartpar+j-C_diff%ys 
              k = C_diff%zs
              !              xsrc(   :   ,i,j,k ) = nil
              xsrc(E_dn   ,i,j,k ) = get_ednTOA(kato,iq,newgrid%z(1))*newgrid%dx(li)*newgrid%dy(lj)
            enddo
          enddo
        endif

        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly... setting coefficients ...done'

        call DMDAVecRestoreArrayF90(C_dir%da,edir,xedir,ierr)
        call DMDAVecRestoreArrayF90(C_diff%da,local_b,xsrc,ierr)

        call DMLocalToGlobalBegin(C_diff%da,local_b,ADD_VALUES, b,ierr) ! USE ADD_VALUES, so that also ghosted entries get updated
        call DMLocalToGlobalEnd  (C_diff%da,local_b,ADD_VALUES, b,ierr)

        call VecDestroy(local_b,ierr)

        if(myid.eq.0.and.ldebug) print *,'src Vector Assembly done'
end subroutine

subroutine load_test_optprop(kato,iq)
      integer(iintegers),intent(in) :: kato,iq
      double precision,allocatable,dimension(:) :: hhl1d,dx,dy,dz
      PetscInt :: i,j,k

      type op
        double precision,allocatable,dimension(:,:,:) :: kabs,ksca,g
      end type
      type(op) :: optP
      integer(iintegers),parameter :: glob_Nx=5,glob_Ny=5,glob_Nz=5, zTOA=150000

      if(myid.eq.0.and.ldebug) print *,myid,'Creating Optical Properties here instead of taking them from kato',kato,iq

!      test_dz(1:5)=10
!      test_dz(6:10)=20
!      test_dz(11:15)=100
!      test_dz(15:20)=200

      allocate(optP%kabs(glob_Nx,glob_Ny,glob_Nz))
      allocate(optP%ksca(glob_Nx,glob_Ny,glob_Nz))
      allocate(optP%g   (glob_Nx,glob_Ny,glob_Nz))
      allocate(hhl1d (glob_Nz+1))
      allocate(dz (glob_Nz))

      hhl1d(ubound(hhl1d,1))=0
      do k=ubound(hhl1d,1),2,-1
        hhl1d(k-1)=hhl1d(k)+zTOA/glob_Nz
      enddo
      dz = hhl1d(1:glob_Nz)-hhl1d(2:glob_Nz+1)

      optP%g=0
      optP%kabs=1e-8
      optP%ksca=1e-8

      if(myid.eq.0) optP%ksca(1,1,1) = 1e-3
!      optP%ksca(1:5,1:5,5:15) = .1

      allocate(dx(ubound(optP%kabs,1))) ; dx=ident_dx
      allocate(dy(ubound(optP%kabs,2))) ; dy=ident_dy

      call init_grid_transform(dx,dy,dz,.True.)
      call grid_old_to_new(optP%kabs)
      call grid_old_to_new(optP%ksca)
      call grid_old_to_new(optP%g)

      if(allocated(pcc_optprop)) deallocate(pcc_optprop)
      allocate(pcc_optprop(ubound(optP%kabs,1),ubound(optP%kabs,2),ubound(optP%kabs,3) ))
      do k=1,ubound(optP%kabs,3)
        do j=1,ubound(optP%kabs,2)
          do i=1,ubound(optP%kabs,1)
            pcc_optprop(i,j,k)%bg(:)   = [ optP%kabs(i,j,k), optP%ksca(i,j,k), optP%g(i,j,k) ]
            pcc_optprop(i,j,k)%fg(1:3) = pcc_optprop(i,j,k)%bg(:)
            pcc_optprop(i,j,k)%fg(4)   = zero

            pcc_optprop(i,j,k)%kext1 = pcc_optprop(i,j,k)%bg(1) + pcc_optprop(i,j,k)%bg(2)
            pcc_optprop(i,j,k)%kext2 = pcc_optprop(i,j,k)%kext1 + pcc_optprop(i,j,k)%fg(1) + pcc_optprop(i,j,k)%fg(2)
            pcc_optprop(i,j,k)%ksca1 = pcc_optprop(i,j,k)%bg(2)
            pcc_optprop(i,j,k)%ksca2 = pcc_optprop(i,j,k)%bg(2)+pcc_optprop(i,j,k)%fg(2)

            pcc_optprop(i,j,k)%w1  = pcc_optprop(i,j,k)%ksca1/pcc_optprop(i,j,k)%kext1
            pcc_optprop(i,j,k)%w2  = pcc_optprop(i,j,k)%ksca2/pcc_optprop(i,j,k)%kext2
            pcc_optprop(i,j,k)%g1  = pcc_optprop(i,j,k)%bg(3)
            pcc_optprop(i,j,k)%g2  = (pcc_optprop(i,j,k)%bg(3)*pcc_optprop(i,j,k)%bg(2) + pcc_optprop(i,j,k)%fg(3)*pcc_optprop(i,j,k)%fg(2))/ (pcc_optprop(i,j,k)%bg(2)+pcc_optprop(i,j,k)%fg(2) )
          enddo
        enddo
      enddo
  end subroutine
  subroutine load_optprop(kato,iq)
      integer(iintegers),intent(in) :: kato,iq
      double precision,allocatable,dimension(:,:,:) :: tmp
      double precision,allocatable,dimension(:) :: hhl1d,dx,dy,dz
      PetscInt :: i,j,k
      double precision,parameter :: min_coeff=1e-30

      character(100) :: skato,siq
      character(100),parameter :: uvspec_file='/usr/users/jakub/data/tenstream/scenes/scenes.h5'

      type op
        double precision,allocatable,dimension(:,:,:) :: kabs,ksca,g
      end type
      type(op) :: optP

      call h5load([uvspec_file,ident,'hhl'],hhl1d,iierr)
      if(myid.eq.0) print *,'hhl',hhl1d

      write(skato,FMT='("kato",I0)') kato
      write(siq  ,FMT='("iq",I0)') iq
      call h5load([uvspec_file,ident,skato,siq,'kabs'],optP%kabs,iierr) 
      call h5load([uvspec_file,ident,skato,siq,'ksca'],optP%ksca,iierr) 
      call h5load([uvspec_file,ident,skato,siq,'g   '],optP%g,iierr) 

      if(size(optP%ksca).ne.size(optP%kabs).or.size(optP%kabs).ne.size(optP%g).or.size(hhl1d).ne.ubound(optP%kabs,3)+1) then !.or.(ubound(hhl1d,1)-1.ne.ubound(kabs,3))) then
        print *,'ERROR : shapes of optical properties do not match!!!'
        print *,'shape(kabs)',shape(optP%kabs)
        print *,'shape(ksca)',shape(optP%ksca)
        print *,'shape(g)',shape(optP%g)
        print *,'shape(hhl1d)',shape(hhl1d)

        call exit()
      endif

      optP%kabs = max(min_coeff, optP%kabs )
      optP%ksca = max(min_coeff, optP%ksca )
      optP%g    = max(min_coeff, optP%g    )


      allocate(dz(size(hhl1d)-1))
      dz = hhl1d(1:ubound(hhl1d,1)-1) - hhl1d(2:ubound(hhl1d,1)) 

      if(myid.eq.0) then
        do k=1,ubound(optP%kabs,3)
          print *,myid,'Optical Properties:',k,'hhl',hhl1d(k),'dz',dz(k),'k',minval(optP%kabs(:,:,k)),minval(optP%ksca(:,:,k)),minval(optP%g(:,:,k)),maxval(optP%kabs(:,:,k)),maxval(optP%ksca(:,:,k)),maxval(optP%g(:,:,k))
        enddo
      endif

      ! Make sure that our domain has at least 3 entries in each dimension.... otherwise violates boundary conditions
      if(ubound(optP%kabs,2).eq.1) then
        allocate(tmp(ubound(optP%kabs,1),3,ubound(optP%kabs,3) ) )
        tmp(:,1,:) = optP%kabs(:,1,:) ; tmp(:,2,:) = optP%kabs(:,1,:) ; tmp(:,3,:) = optP%kabs(:,1,:)
        deallocate(optP%kabs) ; allocate(optP%kabs(ubound(tmp,1),ubound(tmp,2),ubound(tmp,3) ) )
        optP%kabs = tmp ; deallocate(tmp)

        allocate(tmp(ubound(optP%ksca,1),3,ubound(optP%ksca,3) ) )
        tmp(:,1,:) = optP%ksca(:,1,:) ; tmp(:,2,:) = optP%ksca(:,1,:) ; tmp(:,3,:) = optP%ksca(:,1,:)
        deallocate(optP%ksca) ; allocate(optP%ksca(ubound(tmp,1),ubound(tmp,2),ubound(tmp,3) ) )
        optP%ksca = tmp ; deallocate(tmp)

        allocate(tmp(ubound(optP%g,1),3,ubound(optP%g,3) ) )
        tmp(:,1,:) = optP%g(:,1,:) ; tmp(:,2,:) = optP%g(:,1,:) ; tmp(:,3,:) = optP%g(:,1,:)
        deallocate(optP%g) ; allocate(optP%g(ubound(tmp,1),ubound(tmp,2),ubound(tmp,3) ) )
        optP%g = tmp ; deallocate(tmp)
      endif

      allocate(dx(ubound(optP%kabs,1))) ; dx=ident_dx
      allocate(dy(ubound(optP%kabs,2))) ; dy=ident_dy

      solution_dx = ident_dx
      solution_dy = ident_dy
      solution_dz = ident_dz

      call init_grid_transform(dx,dy,dz,.True.)
      call grid_old_to_new(optP%kabs)
      call grid_old_to_new(optP%ksca)
      call grid_old_to_new(optP%g)

      if(allocated(pcc_optprop)) deallocate(pcc_optprop)
      allocate(pcc_optprop(ubound(optP%kabs,1),ubound(optP%kabs,2),ubound(optP%kabs,3) ))
      do k=1,ubound(optP%kabs,3)
        do j=1,ubound(optP%kabs,2)
          do i=1,ubound(optP%kabs,1)
            pcc_optprop(i,j,k)%bg(:)   = [ optP%kabs(i,j,k), optP%ksca(i,j,k), optP%g(i,j,k) ]
            pcc_optprop(i,j,k)%fg(1:3) = nil
            pcc_optprop(i,j,k)%fg(4)   = zero

            pcc_optprop(i,j,k)%kext1 = pcc_optprop(i,j,k)%bg(1) + pcc_optprop(i,j,k)%bg(2)
            pcc_optprop(i,j,k)%kext2 = pcc_optprop(i,j,k)%kext1 + pcc_optprop(i,j,k)%fg(1) + pcc_optprop(i,j,k)%fg(2)
            pcc_optprop(i,j,k)%ksca1 = pcc_optprop(i,j,k)%bg(2)
            pcc_optprop(i,j,k)%ksca2 = pcc_optprop(i,j,k)%bg(2)+pcc_optprop(i,j,k)%fg(2)

            pcc_optprop(i,j,k)%w1  = pcc_optprop(i,j,k)%ksca1/pcc_optprop(i,j,k)%kext1
            pcc_optprop(i,j,k)%w2  = pcc_optprop(i,j,k)%ksca2/pcc_optprop(i,j,k)%kext2
            pcc_optprop(i,j,k)%g1  = pcc_optprop(i,j,k)%bg(3)
            pcc_optprop(i,j,k)%g2  = (pcc_optprop(i,j,k)%bg(3)*pcc_optprop(i,j,k)%bg(2) + pcc_optprop(i,j,k)%fg(3)*pcc_optprop(i,j,k)%fg(2))/ (pcc_optprop(i,j,k)%bg(2)+pcc_optprop(i,j,k)%fg(2) )
          enddo
        enddo
      enddo
end subroutine

subroutine local_optprop(glob_optP)
      type(t_optprop),allocatable :: glob_optP(:,:,:)
      type(t_optprop),allocatable :: loc_optP(:,:,:)
      integer :: k

      if(myid.eq.0.and.ldebug) print *,myid,'copying optprop: global to local :: shape kabs',shape(glob_optP),'xstart/end',C_one%xs,C_one%xe,'ys/e',C_one%ys,C_one%ye

      do k=1,3
        if( any(glob_optP%bg(k).lt.zero) .or. any(isnan(glob_optP%bg(k))) ) then
          print *,'local_optprop :: corrupt global optical properties: glob_optP:: '
          call exit
        endif
      enddo

      allocate(loc_optP(C_one%xm, C_one%ym, C_one%zm))

      loc_optP(:,:,:) = glob_optP(C_one%xs+1:C_one%xe+1, C_one%ys+1:C_one%ye+1, :)
      deallocate(glob_optP)
      allocate(glob_optP(C_one%xm, C_one%ym, C_one%zm))
      glob_optP=loc_optP
      deallocate(loc_optP)

      do k=1,3
        if( any(glob_optP%bg(k).lt.zero) .or. any(isnan(glob_optP%bg(k))) ) then
          print *,'local_optprop :: corrupt local optical properties: glob_optP%bg(k) :: '
          call exit
        endif
      enddo
  end subroutine

subroutine calc_flx_div(edir,ediff,abso)
        Vec :: edir,ediff,abso
        PetscReal,pointer,dimension(:,:,:,:) :: xediff,xedir,xabso
        PetscInt :: i,j,k,li,lj,lk
        Vec :: ledir,lediff ! local copies of vectors, including ghosts
        PetscReal :: div(13),Ax,Ay,Az

        if(myid.eq.0.and.ldebug) print *,'Calculating flux divergence'
        call VecSet(abso,zero,ierr)
        
        call DMCreateLocalVector(C_diff%da,lediff,ierr) ; call VecSet(lediff,zero,ierr)
        call DMCreateLocalVector(C_dir%da ,ledir ,ierr) ; call VecSet(ledir ,zero,ierr)
        call DMGlobalToLocalBegin(C_dir%da ,edir ,ADD_VALUES,ledir ,ierr)
        call DMGlobalToLocalEnd(C_dir%da ,edir ,ADD_VALUES,ledir ,ierr)

        call DMGlobalToLocalBegin(C_diff%da,ediff,ADD_VALUES,lediff,ierr)
        call DMGlobalToLocalEnd(C_diff%da,ediff,ADD_VALUES,lediff,ierr)

        ! calculate absorption by flux divergence
        call DMDAVecGetArrayF90(C_diff%da,lediff,xediff,ierr)
        call DMDAVecGetArrayF90(C_dir%da ,ledir ,xedir ,ierr)
        call DMDAVecGetArrayF90(C_one%da ,abso ,xabso ,ierr)

        do k=C_one%zs,C_one%ze
        do j=C_one%ys,C_one%ye         
        do i=C_one%xs,C_one%xe      
          li = istartpar+i-C_one%xs 
          lj = jstartpar+j-C_one%ys 
          lk = 1+k-C_one%zs
                Ax = newgrid%dy(lj)*newgrid%dz(lk)
                Ay = newgrid%dx(li)*newgrid%dz(lk)
                Az = newgrid%dx(li)*newgrid%dy(lj)

                ! Divergence    =                       Incoming                -       Outgoing
                div( 1) = sum( xedir(i0:i3, i,j,k)          - xedir(i0:i3, i,j,k+i1  ) ) *Az*.25
                div( 2) = sum( xedir(i4:i5, i+i1-xinc,j,k)  - xedir(i4:i5, i+xinc,j,k) ) *Ax*.5
                div( 3) = sum( xedir(i6:i7, i,j+i1-yinc,k)  - xedir(i6:i7, i,j+yinc,k) ) *Ay*.5
                                                                     
                div( 4) = Az* ( xediff(E_up  ,i  ,j  ,k+1)  - xediff(E_up  ,i  ,j  ,k  )  )
                div( 5) = Az* ( xediff(E_dn  ,i  ,j  ,k  )  - xediff(E_dn  ,i  ,j  ,k+1)  )
                div( 6) = Ax* ( xediff(E_le_m,i+1,j  ,k  )  - xediff(E_le_m,i  ,j  ,k  )  )
                div( 7) = Ax* ( xediff(E_le_p,i+1,j  ,k  )  - xediff(E_le_p,i  ,j  ,k  )  )
                div( 8) = Ax* ( xediff(E_ri_m,i  ,j  ,k  )  - xediff(E_ri_m,i+1,j  ,k  )  )
                div( 9) = Ax* ( xediff(E_ri_p,i  ,j  ,k  )  - xediff(E_ri_p,i+1,j  ,k  )  )
                div(10) = Ay* ( xediff(E_ba_m,i  ,j+1,k  )  - xediff(E_ba_m,i  ,j  ,k  )  )
                div(11) = Ay* ( xediff(E_ba_p,i  ,j+1,k  )  - xediff(E_ba_p,i  ,j  ,k  )  )
                div(12) = Ay* ( xediff(E_fw_m,i  ,j  ,k  )  - xediff(E_fw_m,i  ,j+1,k  )  )
                div(13) = Ay* ( xediff(E_fw_p,i  ,j  ,k  )  - xediff(E_fw_p,i  ,j+1,k  )  )

                xabso(i0,i,j,k) = sum(div)/ ( newgrid%dx(li)*newgrid%dy(lj)*newgrid%dz(lk) )
        enddo                             
        enddo                             
        enddo   
        
        call DMDAVecRestoreArrayF90(C_one%da ,abso ,xabso ,ierr)
        call DMDAVecRestoreArrayF90(C_diff%da,lediff,xediff,ierr)
        call DMDAVecRestoreArrayF90(C_dir%da ,ledir ,xedir ,ierr)
        
        call VecDestroy(lediff,ierr)
        call VecDestroy(ledir ,ierr)
end subroutine

subroutine scale_flx(edir,ediff)
        Vec :: edir,ediff
        PetscReal,pointer,dimension(:,:,:,:) :: xediff,xedir
        PetscInt :: i,j,k,li,lj,lk
        PetscReal :: Ax,Ay,Az

        if(myid.eq.0.and.ldebug) print *,'rescaling fluxes'
        call DMDAVecGetArrayF90(C_diff%da,ediff,xediff,ierr)
        call DMDAVecGetArrayF90(C_dir%da ,edir ,xedir ,ierr)
        ! rescale total energy fluxes to average quantities i.e. W/m**2 or W/m**3
        do k=C_one%zs,C_one%ze
        do j=C_one%ys,C_one%ye         
        do i=C_one%xs,C_one%xe      
          li = istartpar+i-C_one%xs 
          lj = jstartpar+j-C_one%ys 
          lk = i1+k-C_one%zs
                Ax = newgrid%dy(lj)*newgrid%dz(lk)
                Ay = newgrid%dx(li)*newgrid%dz(lk)
                Az = newgrid%dx(li)*newgrid%dy(lj)

                xedir(i0:i3,i,j,k) = xedir(i0:i3,i,j,k) / ( Az*.25 )
                xedir(i4:i5,i,j,k) = xedir(i4:i5,i,j,k) / ( Ax*.5  )
                xedir(i6:i7,i,j,k) = xedir(i6:i7,i,j,k) / ( Ay*.5  )

                xediff(E_up  ,i,j,k) = xediff(E_up  ,i,j,k) / Az
                xediff(E_dn  ,i,j,k) = xediff(E_dn  ,i,j,k) / Az
                xediff(E_le_m,i,j,k) = xediff(E_le_m,i,j,k) / Ax
                xediff(E_le_p,i,j,k) = xediff(E_le_p,i,j,k) / Ax
                xediff(E_ri_m,i,j,k) = xediff(E_ri_m,i,j,k) / Ax
                xediff(E_ri_p,i,j,k) = xediff(E_ri_p,i,j,k) / Ax
                xediff(E_ba_m,i,j,k) = xediff(E_ba_m,i,j,k) / Ay
                xediff(E_ba_p,i,j,k) = xediff(E_ba_p,i,j,k) / Ay
                xediff(E_fw_m,i,j,k) = xediff(E_fw_m,i,j,k) / Ay
                xediff(E_fw_p,i,j,k) = xediff(E_fw_p,i,j,k) / Ay
        enddo
        enddo
        enddo

        do j=C_diff%ys,C_diff%ye         
        do i=C_diff%xs,C_diff%xe      
                k=C_diff%ze
                li = istartpar+i-C_diff%xs 
                lj = jstartpar+j-C_diff%ys 
                lk = 1+k-C_diff%zs         
                Az = newgrid%dx(li)*newgrid%dy(lj)
                xedir (i0:i3 ,i,j,k) = xedir (i0:i3 ,i,j,k) / ( Az*.25 )
                xediff(E_up  ,i,j,k) = xediff(E_up  ,i,j,k) / Az
                xediff(E_dn  ,i,j,k) = xediff(E_dn  ,i,j,k) / Az
        enddo
        enddo

        call DMDAVecRestoreArrayF90(C_diff%da,ediff,xediff,ierr)
        call DMDAVecRestoreArrayF90(C_dir%da ,edir ,xedir ,ierr)
end subroutine

subroutine hdf5_to_vec(v,ie)
      Vec :: v
      PetscErrorCode,intent(out) :: ie

      character(10),parameter :: suffix='.h5'
      character(110) :: fname
      logical fexists
      PetscFileMode :: fmode
      PetscViewer :: view

      fname = trim(basepath) // 'ts.' // trim(ident)  // trim(suffix)
      inquire(file=trim(fname), exist=fexists)
      
      if(fexists) then
        if(myid.eq.0)  print *,myid,'loading vector-data from hdf5 file',trim(fname)
        fmode = FILE_MODE_READ
      else 
          ie=-1
          return
      endif

      call PetscViewerHDF5Open(PETSC_COMM_WORLD,trim(fname),fmode, view, ierr)
      call VecLoad(v, view, ie)
      call PetscViewerDestroy(view,ierr)
end subroutine
subroutine vec_from_hdf5(v,err_code)
      Vec :: v
      character(10),parameter :: suffix='.h5'
      character(110) :: fname
      logical fexists
      PetscFileMode :: fmode
      character(100) :: vecname,s_theta0

      PetscInt :: err_code
      
      PetscViewer view
      PetscErrorCode ierr

      call PetscObjectGetName(v,vecname,ierr)

      write(s_theta0,FMT='(".",I0)' ) int(theta0)
      
      fname = trim(basepath) // 'ts.' // trim(ident) // trim(s_theta0) // trim(suffix)
      inquire(file=trim(fname), exist=fexists)
      
      if(myid.eq.0)  print *,myid,'reading vector from hdf5 file ',trim(fname),' vecname: ',vecname
      fmode = FILE_MODE_READ

      call PetscViewerHDF5Open(PETSC_COMM_WORLD,trim(fname),fmode, view, ierr)
      call VecLoad(v, view, err_code)
      call PetscViewerDestroy(view,ierr)
      if(myid.eq.0.and.ldebug) print *,myid,'reading from hdf5 file done :: error:',err_code
end subroutine
subroutine vec_to_hdf5(v)
      Vec,intent(in) :: v
      character(10),parameter :: suffix='.h5'
      character(110) :: fname
      logical fexists
      PetscFileMode :: fmode
      character(100) :: vecname,s_theta0
      
      PetscViewer view
      PetscErrorCode ierr

      call PetscObjectGetName(v,vecname,ierr)

      write(s_theta0,FMT='(".",I0)' ) int(theta0)
      
      fname = trim(basepath) // 'ts.' // trim(ident) // trim(s_theta0) // trim(suffix)
      inquire(file=trim(fname), exist=fexists)
      
      if(fexists) then
        if(myid.eq.0)  print *,myid,'appending vector to hdf5 file ',trim(fname),' vecname: ',vecname
        fmode = FILE_MODE_APPEND
      else 
        if(myid.eq.0)  print *,myid,'writing vector to hdf5 file ',trim(fname),' vecname: ',vecname
        fmode = FILE_MODE_WRITE
      endif

      call PetscViewerHDF5Open(PETSC_COMM_WORLD,trim(fname),fmode, view, ierr)
      call VecView(v, view, ierr)
      call PetscViewerDestroy(view,ierr)
      if(myid.eq.0.and.ldebug) print *,myid,'writing to hdf5 file done'
end subroutine

subroutine solve(ksp,b,A,x,C)
      KSP :: ksp
      Mat:: A
      Vec:: b
      Vec:: x
      type(coord) :: C

      KSPConvergedReason :: reason
      PetscInt :: iter

      if(myid.eq.0.and.ldebug) print *,'Solving Matrix'

      call KSPSetOperators(ksp,A,A,ierr)
!      call KSPSetOperators(ksp,A,A,SAME_NONZERO_PATTERN,ierr)
!      call KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)

      call KSPSetFromOptions(ksp,ierr)
      call KSPSetDM(ksp,C%da,ierr)
      call KSPSetDMActive(ksp,PETSC_FALSE,ierr)
      call KSPSetUp(ksp,ierr)

      call KSPSolve(ksp,b,x,ierr)
      call KSPGetIterationNumber(ksp,iter,ierr)
      call KSPGetConvergedReason(ksp,reason,ierr)

!      print *,'Source Vector:'
!      call VecView(b,PETSC_VIEWER_STDOUT_WORLD ,ierr)
!
!      print *,'Solution Vector:'
!      call VecView(x,PETSC_VIEWER_STDOUT_WORLD ,ierr)

      if(myid.eq.0) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      if(reason.eq.KSP_DIVERGED_ITS) then
        if(myid.eq.0) print *,'We take our chances, that this is a meaningful output.... and just go on'
        return
      endif

      if(reason.le.0) then
              if(myid.eq.0) print *,myid,'Resetted initial guess to zero and try again:'
              call VecSet(x,zero,ierr)
              call KSPSolve(ksp,b,x,ierr)
              call KSPGetIterationNumber(ksp,iter,ierr)
              call KSPGetConvergedReason(ksp,reason,ierr)
              if(myid.eq.0) print *,myid,'Solver took ',iter,' iterations and converged',reason.gt.0,'because',reason
      endif

      if(reason.le.0) then
              print *,'*********************************************************** SOLVER did NOT converge :( ***************************************************88'
              call exit()
      endif
end subroutine
subroutine setup_ksp(ksp)
      KSP :: ksp
      call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)

end subroutine

end module

program main
        use petsc_ts
        implicit none

        Mat :: Mdiff,Mdir
        Vec :: b,x,edir,intedir,intx,incSolar,abso
!        PetscInt :: iter

        KSP :: kspdir,kspdiff

        integer(iintegers) :: iq
        integer(iintegers) :: kato

        character(100) :: vecname

        call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
        call MPI_COMM_RANK( PETSC_COMM_WORLD, myid, ierr )

        call PetscLogStageRegister('total_tenstream',logstage(1), ierr)
        call PetscLogStageRegister('setup_edir',logstage(2), ierr)
        call PetscLogStageRegister('calc_edir',logstage(3), ierr)
        call PetscLogStageRegister('setup_ediff',logstage(4), ierr)
        call PetscLogStageRegister('calc_ediff',logstage(5), ierr)
        call PetscLogStageRegister('setup_b',logstage(6), ierr)
        call PetscLogStageRegister('get_coeff',logstage(7), ierr)
        print *,'Loggin stages',logstage

        if(ltest) then
          call load_test_optprop(1_iintegers,0_iintegers)
        else
          call load_optprop(1_iintegers,0_iintegers)
        endif

        call setup_grid()

        ! Write the result vectors once, to ensure that we are able to write the results
        call DMCreateGlobalVector(C_diff%da,intx,ierr)
        write(vecname,FMT='("ediff.",I0,".",I0)') int(phi0),int(theta0); call PetscObjectSetName(intx,vecname,ierr)
        call VecSet(intx,zero,ierr)
!        call vec_to_hdf5(intx)

        call DMCreateGlobalVector(C_dir%da,intedir,ierr)
        write(vecname,FMT='("edir.",I0,".",I0)') int(phi0),int(theta0); call PetscObjectSetName(intedir,vecname,ierr)
        call VecSet(intedir,zero,ierr)
!        call vec_to_hdf5(intedir)

        call DMCreateGlobalVector(C_one%da,abso,ierr)
        write(vecname,FMT='("abso.",I0,".",I0)') int(phi0),int(theta0); call PetscObjectSetName(abso,vecname,ierr)
        call VecSet(abso,zero,ierr)
!        call vec_to_hdf5(abso)

        call init_optprop(newgrid%dx(1),newgrid%dy(1),PETSC_COMM_WORLD) ! i0 is LUT, i1 is ANN

        call setup_dir_inc(phi0)

        ! Create Objects to work with
        call init_Matrix(Mdir,C_dir)
        call init_Matrix(Mdiff,C_diff)
        call setup_ksp(kspdir)
        call KSPAppendOptionsPrefix(kspdir,"dir_",ierr)
        call setup_ksp(kspdiff)
        call KSPAppendOptionsPrefix(kspdiff,"diff_",ierr)

        call DMCreateGlobalVector(C_dir%da,edir,ierr)
        call DMCreateGlobalVector(C_dir%da,incSolar,ierr)
        call DMCreateGlobalVector(C_diff%da,b,ierr)
        call DMCreateGlobalVector(C_diff%da,x,ierr)

        call VecSet(edir,zero,ierr)
        call VecSet(x,zero,ierr)

        do kato=1,32
!                  do kato=11,11
          do iq=0,kato_bands(kato)
            if(myid.eq.0) print *,'                            ----------------------------------------------------------------------------------------'
            if(myid.eq.0) print *,'                            -------------------------- Calculate kato',kato,'iq',iq,'-------------------------------'
            if(myid.eq.0) print *,'                            ----------------------------------------------------------------------------------------'


            !First try reading them from file and we can resume from there:
            if(l_writeall) then
              write(vecname,FMT='("edir.",I0,".",I0)') kato,iq; call PetscObjectSetName(edir,vecname,ierr)
              call vec_from_hdf5(edir,ierr)
              write(vecname,FMT='("x.",I0,".",I0)') kato,iq; call PetscObjectSetName(x,vecname,ierr)
              call vec_from_hdf5(x,ierr)
            endif

            if(ltest) then
              call load_test_optprop(kato,iq)
            else
              call load_optprop(kato,iq)
            endif
            call local_optprop(pcc_optprop)

            call PetscLogStagePush(logstage(1),ierr)

            if(myid.eq.0) print *,'Calculate Edir'
            call PetscLogStagePush(logstage(2),ierr)
            call set_dir_coeff(Mdir,C_dir)
            call setup_incSolar(incSolar,kato,iq)
            if(l_writeall) then
              write(vecname,FMT='("incSolar.",I0,".",I0)') kato,iq; call PetscObjectSetName(incSolar,vecname,ierr)
              call vec_to_hdf5(incSolar)
            endif
            call PetscLogStagePop(ierr)

            call PetscLogStagePush(logstage(3),ierr)
            call solve(kspdir,incSolar,Mdir,edir,C_dir)
            if(l_writeall) then
              write(vecname,FMT='("edir.",I0,".",I0)') kato,iq; call PetscObjectSetName(edir,vecname,ierr)
              call vec_to_hdf5(edir)
            endif
            call PetscLogStagePop(ierr)

            if(myid.eq.0) print *,'Calculate diffuse src term'
            call PetscLogStagePush(logstage(6),ierr)
            call setup_b(edir,kato,iq,b)
            if(l_writeall) then
              write(vecname,FMT='("b.",I0,".",I0)') kato,iq; call PetscObjectSetName(b,vecname,ierr)
              call vec_to_hdf5(b)
            endif
            call PetscLogStagePop(ierr)

            if(myid.eq.0) print *,'Calculate diffuse Radition'
            call PetscLogStagePush(logstage(4),ierr)
            call set_diff_coeff(Mdiff,C_diff)
            call PetscLogStagePop(ierr)

            call PetscLogStagePush(logstage(5),ierr)
            call solve(kspdiff,b,Mdiff,x,C_diff)
            call PetscLogStagePop(ierr)

            if(l_writeall) then
              write(vecname,FMT='("x.",I0,".",I0)') kato,iq; call PetscObjectSetName(x,vecname,ierr)
              call vec_to_hdf5(x)
            endif
            call VecAXPY(intedir,one,edir,ierr)
            call VecAXPY(intx,one,x,ierr)

            call PetscLogStagePop(ierr)
          enddo
        enddo

        call VecDestroy(incSolar,ierr)
        call VecDestroy(edir,ierr)
        call VecDestroy(b,ierr)
        call VecDestroy(x,ierr)
        call MatDestroy(Mdir,ierr)
        call MatDestroy(Mdiff,ierr)

        call KSPDestroy(kspdir,ierr)
        call KSPDestroy(kspdiff,ierr)

        print *,'scale Energy from W to W/m^2 or W/m^3'
        call scale_flx(intedir,intx)

        print *,'Calculating divergence'
        call calc_flx_div(intedir,intx,abso)

        call vec_to_hdf5(abso)
        call vec_to_hdf5(intedir)
        call vec_to_hdf5(intx)

        print *,'Cleanup Result vectors'
        call VecDestroy(intedir,ierr)
        call VecDestroy(intx,ierr)
        call VecDestroy(abso,ierr)

        !!        if(myid.eq.0) read(*,*) iter
        call PetscFinalize(ierr)
end program
