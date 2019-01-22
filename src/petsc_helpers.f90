module m_petsc_helpers
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    zero, i0, i1, i2, i3

  use m_helper_functions, only : get_arg, CHKERR

  implicit none

  private
  public :: petscGlobalVecToZero, scatterZerotoPetscGlobal, &
    petscVecToF90, &
    f90VecToPetsc, &
    getVecPointer, restoreVecPointer, &
    dmda_convolve_ediff_srfc

  interface f90VecToPetsc
    module procedure f90VecToPetsc_3d, f90VecToPetsc_4d
  end interface
  interface petscVecToF90
    module procedure petscVecToF90_3d, petscVecToF90_4d
  end interface
  interface getVecPointer
    module procedure getVecPointer_2d, getVecPointer_3d
  end interface
  interface restoreVecPointer
    module procedure restoreVecPointer_2d, restoreVecPointer_3d
  end interface

  logical, parameter :: ldebug=.True.
contains

  !> @brief Scatter a petsc global vector into a local vector on Rank 0
  !> @details The local Vec will be in natural ordering.
  !>  \n you may use this routine e.g. to gather the results from mpi parallel vectors into a sequential calling program.
  !>  \n lVec will be created and you should VecDestroy it when you are done.
  subroutine petscGlobalVecToZero(gVec, dm, lVec)
    type(tVec), intent(in)    :: gVec
    type(tDM), intent(in)     :: dm
    type(tVec), intent(out)   :: lVec

    type(tVec) :: natural
    type(tVecScatter) :: scatter_context

    integer(mpiint) :: ierr

    call DMDACreateNaturalVector(dm, natural, ierr); call CHKERR(ierr)

    call DMDAGlobalToNaturalBegin(dm, gVec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)
    call DMDAGlobalToNaturalEnd  (dm, gVec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)

    call VecScatterCreateToZero(natural, scatter_context, lVec, ierr); call CHKERR(ierr)

    call VecScatterBegin(scatter_context, natural, lVec, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
    call VecScatterEnd  (scatter_context, natural, lVec, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

    call VecScatterDestroy(scatter_context, ierr); call CHKERR(ierr)

    call VecDestroy(natural,ierr); call CHKERR(ierr)
  end subroutine

  !> @brief Scatter a local array on rank0 vector into a petsc global vector
  !> @details you may use this routine e.g. to scatter the optical properties from a sequential calling program.
  subroutine scatterZerotoPetscGlobal(arr, dm, vec)
    real(ireals),allocatable,intent(in) :: arr(:,:,:)
    type(tDM), intent(in)     :: dm
    type(tVec) :: vec

    type(tVecScatter) :: scatter_context
    type(tVec) :: natural,local
    real(ireals), pointer :: xloc(:)=>null()
    integer(mpiint) :: comm, myid, ierr

    call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: Create Natural Vec'
    call DMDACreateNaturalVector(dm, natural, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: Create scatter ctx'
    call VecScatterCreateToZero(natural, scatter_context, local, ierr); call CHKERR(ierr)

    if(myid.eq.0) then
      if(.not. allocated(arr)) stop 'Cannot call scatterZerotoPetscGlobal with unallocated input array'
      if(ldebug) print *,myid,'scatterZerotoDM :: Copy data from Fortran array to Local Petsc Vec'
      call VecGetArrayF90(local,xloc,ierr) ;call CHKERR(ierr)
      xloc = reshape( arr , [ size(arr) ] )
      call VecRestoreArrayF90(local,xloc,ierr) ;call CHKERR(ierr)
    endif

    if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: scatter reverse....'
    call VecScatterBegin(scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)
    call VecScatterEnd  (scatter_context, local, natural, INSERT_VALUES, SCATTER_REVERSE, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: natural to global....'
    call DMDANaturalToGlobalBegin(dm,natural, INSERT_VALUES, vec, ierr); call CHKERR(ierr)
    call DMDANaturalToGlobalEnd  (dm,natural, INSERT_VALUES, vec, ierr); call CHKERR(ierr)

    if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: destroying contexts....'
    call VecScatterDestroy(scatter_context, ierr); call CHKERR(ierr)
    call VecDestroy(local,ierr); call CHKERR(ierr)
    call VecDestroy(natural,ierr); call CHKERR(ierr)
    if(ldebug.and.myid.eq.0) print *,myid,'scatterZerotoDM :: done....'
  end subroutine

  !> @brief Copies the data from a petsc vector into an allocatable array
  !> @details if flag opt_l_only_on_rank0 is True,
  !>     \n we assume this is just a local vector on rank 0, i.e. coming petscGlobalVecToZero()
  subroutine petscVecToF90_4d(vec, dm, arr, opt_l_only_on_rank0)
    type(tVec), intent(in)    :: vec
    type(tDM), intent(in)     :: dm
    real(ireals), intent(inout), allocatable :: arr(:,:,:,:)
    logical, intent(in), optional :: opt_l_only_on_rank0
    logical :: l_only_on_rank0

    integer(iintegers) :: vecsize
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: comm, myid, ierr
    integer(iintegers) :: dmdim, dof, zs, xs, ys, zm, xm, ym, glob_zm, glob_xm, glob_ym

    if(allocated(arr)) stop 'You shall not call petscVecToF90 with an already allocated array!'

    call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    l_only_on_rank0 = get_arg(.False., opt_l_only_on_rank0)

    if(l_only_on_rank0 .and. myid.ne.0) stop 'Only rank 0 should call the routine petscVecToF90 with opt_l_only_on_rank0=.T.'

    call DMDAGetInfo(dm, dmdim, glob_zm, glob_xm, glob_ym,        &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      dof               , PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      ierr) ;call CHKERR(ierr)

    call DMDAGetCorners(dm, zs, xs, ys, zm, xm, ym, ierr) ;call CHKERR(ierr)

    if(.not.l_only_on_rank0) then
      call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
      if(.not.allocated(arr)) allocate(arr(dof, zm, xm, ym))
    else
      call VecGetSize(vec, vecsize, ierr); call CHKERR(ierr)
      if(.not.allocated(arr)) allocate(arr(dof, glob_zm, glob_xm, glob_ym))
    endif

    if(vecsize.ne.size(arr)) then
      print *,'petscVecToF90 Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      stop 'petscVecToF90 Vecsizes dont match!'
    endif

    if(.not.l_only_on_rank0) then
      call getVecPointer(vec, dm, x1d, x4d)
      arr = x4d
      call restoreVecPointer(vec, x1d, x4d)
    else
      call VecGetArrayF90(vec,x1d,ierr); call CHKERR(ierr)
      arr = reshape( x1d, (/ dof, glob_zm, glob_xm, glob_ym /) )
      call VecRestoreArrayF90(vec,x1d,ierr); call CHKERR(ierr)
    endif
  end subroutine
  subroutine petscVecToF90_3d(vec, dm, arr, opt_l_only_on_rank0)
    type(tVec), intent(in)    :: vec
    type(tDM), intent(in)     :: dm
    real(ireals), intent(inout), allocatable :: arr(:,:,:)
    logical, intent(in), optional :: opt_l_only_on_rank0
    logical :: l_only_on_rank0

    integer(iintegers) :: vecsize
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: comm, myid, ierr
    integer(iintegers) :: dmdim, dof, zs, xs, ys, zm, xm, ym, glob_zm, glob_xm, glob_ym

    l_only_on_rank0 = get_arg(.False., opt_l_only_on_rank0)

    call DMDAGetInfo(dm, dmdim, glob_zm, glob_xm, glob_ym,        &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      dof               , PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      ierr) ;call CHKERR(ierr)

    call DMDAGetCorners(dm, zs, xs, ys, zm, xm, ym, ierr) ;call CHKERR(ierr)

    if(dof.ne.1) stop 'petscVecToF90_3d should only be called with anything else than DM%dof of 1'

    call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)


    if(l_only_on_rank0 .and. myid.ne.0) stop 'Only rank 0 should call the routine petscVecToF90 with opt_l_only_on_rank0=.T.'

    if(.not.l_only_on_rank0) then
      call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
      if(.not.allocated(arr)) allocate(arr(zm,xm,ym))
    else
      call VecGetSize(vec, vecsize, ierr); call CHKERR(ierr)
      if(.not.allocated(arr)) allocate(arr(glob_zm, glob_xm, glob_ym))
    endif

    if(vecsize.ne.size(arr)) then
      print *,'petscVecToF90 Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      stop 'petscVecToF90 Vecsizes dont match!'
    endif

    if(.not.l_only_on_rank0) then
      call getVecPointer(vec, dm, x1d, x4d)
      arr = x4d(i0,:,:,:)
      call restoreVecPointer(vec, x1d, x4d)
    else
      call VecGetArrayF90(vec,x1d,ierr); call CHKERR(ierr)
      arr = reshape( x1d, (/ glob_zm, glob_xm, glob_ym /) )
      call VecRestoreArrayF90(vec,x1d,ierr); call CHKERR(ierr)
    endif

  end subroutine

  !> @brief Copies the data from a fortran vector into an petsc global vec
  subroutine f90VecToPetsc_4d(arr, dm, vec)
    real(ireals), intent(in)  :: arr(:,:,:,:)
    type(tDM), intent(in)     :: dm
    type(tVec), intent(inout) :: vec

    integer(iintegers) :: vecsize
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: ierr

    call getVecPointer(vec, dm, x1d, x4d)
    x4d = arr
    call restoreVecPointer(vec, x1d, x4d)

    call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
    if(vecsize.ne.size(arr)) then
      print *,'f90VecToPetsc Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      call CHKERR(1_mpiint, 'f90VecToPetsc Vecsizes dont match!')
    endif
  end subroutine
  subroutine f90VecToPetsc_3d(arr, dm, vec)
    real(ireals), intent(in)  :: arr(:,:,:)
    type(tDM), intent(in)     :: dm
    type(tVec), intent(inout) :: vec

    integer(iintegers) :: vecsize
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: ierr

    call getVecPointer(vec, dm, x1d, x4d)
    x4d(i0,:,:,:) = arr
    call restoreVecPointer(vec, x1d, x4d)

    call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
    if(vecsize.ne.size(arr)) then
      print *,'f90VecToPetsc Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      stop 'f90VecToPetsc Vecsizes dont match!'
    endif
  end subroutine

  subroutine getVecPointer_3d(vec,dm,x1d,x4d)
    type(tVec) :: vec
    type(tDM), intent(in) :: dm
    real(ireals),intent(inout),pointer,dimension(:,:,:,:) :: x4d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d

    integer(iintegers) :: N, dmdim, dof, glob_zm, glob_xm, glob_ym
    integer(iintegers) :: zs, ze, xs, xe, ys, ye, zm, xm, ym
    integer(iintegers) :: gzs, gze, gxs, gxe, gys, gye, gzm, gxm, gym
    integer(mpiint) :: comm, myid, ierr
    logical :: lghosted

    if(associated(x1d).or.associated(x4d)) then
      print *,'ERROR : getVecPointer : input vector already associated!!',associated(x1d),associated(x4d)
      call CHKERR(1_mpiint, 'getVecPointer : input vector already associated')
    endif

    call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, '-show_getvecpointerdm', ierr); call CHKERR(ierr)
    call DMDAGetInfo(dm, dmdim, glob_zm, glob_xm, glob_ym,        &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      dof               , PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      ierr) ;call CHKERR(ierr)


    call DMDAGetCorners(dm, zs, xs, ys, zm, xm, ym, ierr) ;call CHKERR(ierr)
    call DMDAGetGhostCorners(dm,gzs,gxs,gys,gzm,gxm,gym,ierr) ;call CHKERR(ierr)
    xe = xs+xm-1
    ye = ys+ym-1
    ze = zs+zm-1
    gxe = gxs+gxm-1
    gye = gys+gym-1
    gze = gzs+gzm-1

    call VecGetLocalSize(vec,N,ierr)

    if( N .eq. dof*xm*ym*zm .or. N .eq. dof*glob_xm*glob_ym*glob_zm) then
      lghosted=.False.
    else if( N .eq. dof*gxm*gym*gzm ) then
      lghosted=.True.
    else
      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      print *,myid,'Size N:', N, dof*xm*ym*zm, dof*glob_xm*glob_ym*glob_zm, dof*gxm*gym*gzm
      stop 'Local Vector dimensions do not conform to DMDA size'
    endif

    call VecGetArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
    if(lghosted) then
      x4d(0:dof-1 , gzs:gze, gxs:gxe , gys:gye ) => x1d
    else
      x4d(0:dof-1 , zs:ze  , xs:xe   , ys:ye   ) => x1d
    endif

  end subroutine
  subroutine restoreVecPointer_3d(vec,x1d,x4d)
    type(tVec) :: vec
    real(ireals),intent(inout),pointer,dimension(:,:,:,:) :: x4d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d
    integer(mpiint) :: ierr

    if(.not.associated(x1d).or..not.associated(x4d)) then
      print *,'ERROR : restoreVecPointer : input vector not yet associated!!',associated(x1d),associated(x4d)
      call exit(1)
    endif

    x4d => null()
    call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
    x1d => null()
  end subroutine
  subroutine getVecPointer_2d(vec,dm,x1d,x3d)
    type(tVec) :: vec
    type(tDM), intent(in) :: dm
    real(ireals),intent(inout),pointer,dimension(:,:,:) :: x3d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d

    integer(iintegers) :: N, dmdim, dof, glob_xm, glob_ym
    integer(iintegers) :: xs, xe, ys, ye, xm, ym
    integer(iintegers) :: gxs, gxe, gys, gye, gxm, gym
    integer(mpiint) :: comm, myid, ierr
    logical :: lghosted

    if(associated(x1d).or.associated(x3d)) then
      print *,'ERROR : getVecPointer : input vector already associated!!',associated(x1d),associated(x3d)
      call CHKERR(1_mpiint, 'getVecPointer : input vector already associated')
    endif

    call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, '-show_getvecpointerdm', ierr); call CHKERR(ierr)
    call DMDAGetInfo(dm, dmdim, glob_xm, glob_ym, PETSC_NULL_INTEGER, &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      dof               , PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      ierr) ;call CHKERR(ierr)


    call DMDAGetCorners(dm, xs, ys, PETSC_NULL_INTEGER, xm, ym, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)
    call DMDAGetGhostCorners(dm,gxs,gys,PETSC_NULL_INTEGER,gxm,gym,PETSC_NULL_INTEGER,ierr) ;call CHKERR(ierr)
    xe = xs+xm-1
    ye = ys+ym-1
    gxe = gxs+gxm-1
    gye = gys+gym-1

    call VecGetLocalSize(vec,N,ierr)

    if( N .eq. dof*xm*ym .or. N .eq. dof*glob_xm*glob_ym) then
      lghosted=.False.
    else if( N .eq. dof*gxm*gym ) then
      lghosted=.True.
    else
      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      print *,myid,'Size N:', N, dof*xm*ym, dof*glob_xm*glob_ym, dof*gxm*gym
      stop 'Local Vector dimensions do not conform to DMDA size'
    endif

    call VecGetArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
    if(lghosted) then
      x3d(0:dof-1 , gxs:gxe , gys:gye ) => x1d
    else
      x3d(0:dof-1 , xs:xe   , ys:ye   ) => x1d
    endif

  end subroutine
  subroutine restoreVecPointer_2d(vec,x1d,x3d)
    type(tVec) :: vec
    real(ireals),intent(inout),pointer,dimension(:,:,:) :: x3d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d
    integer(mpiint) :: ierr

    if(.not.associated(x1d).or..not.associated(x3d)) then
      print *,'ERROR : restoreVecPointer : input vector not yet associated!!',associated(x1d),associated(x3d)
      call exit(1)
    endif

    x3d => null()
    call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
    x1d => null()
  end subroutine

  subroutine dmda_convolve_ediff_srfc(dm3d, kernel_width, arr)
    type(tDM), intent(in) :: dm3d
    integer(iintegers), intent(in) :: kernel_width
    real(ireals), intent(inout) :: arr(:,:,:) ! shape (dof,x,y)


    type(tVec) :: gvec, lvec
    type(tDM) :: dm2d

    real(ireals),pointer,dimension(:,:,:) :: x3d=>null()
    real(ireals),pointer,dimension(:) :: x1d=>null()
    real(ireals),pointer,dimension(:,:,:) :: g3d=>null()
    real(ireals),pointer,dimension(:) :: g1d=>null()

    integer(iintegers) :: Nx, Ny, Ndof, nprocz, nprocx, nprocy
    integer(iintegers), allocatable, dimension(:) :: Nxperproc, Nyperproc
    integer(iintegers) :: idof
    integer(mpiint) :: comm, myid, numnodes, ierr

    Ndof=size(arr, dim=1, kind=iintegers)

    call PetscObjectGetComm(dm3d, comm, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call DMDAGetInfo(dm3d, PETSC_NULL_INTEGER, &
      PETSC_NULL_INTEGER, Nx, Ny, &
      nprocz, nprocx, nprocy, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, &
      PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)

    allocate(Nxperproc(nprocx), Nyperproc(nprocy))

    call DMDAGetOwnershipRanges(dm3d, PETSC_NULL_INTEGER, Nxperproc, Nyperproc, ierr); call CHKERR(ierr)

    if(kernel_width.gt.int(min(size(arr,dim=2),size(arr,dim=3)), iintegers)) then
      call CHKERR(int(kernel_width, mpiint), 'smoothing kernel size is bigger than local domains...'// &
        'this would need more work to be done... go for it...'// &
        '(cascaded uniform filters with updates inbetween or just call this function more often)')
    endif

    call DMDACreate2d(comm, &
      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
      DMDA_STENCIL_BOX, Nx, Ny, &
      nprocx, nprocy, &
      Ndof, kernel_width, &
      Nxperproc, Nyperproc, &
      dm2d, ierr) ;call CHKERR(ierr)
    call DMSetup(dm2d, ierr); call CHKERR(ierr)

    call DMGetGlobalVector(dm2d, gvec, ierr); call CHKERR(ierr)
    call getVecPointer(gvec, dm2d, g1d, g3d)
    g3d(:,:,:) = arr(:,:,:)
    call restoreVecPointer(gvec, g1d, g3d)

    call DMGetLocalVector(dm2d, lvec, ierr); call CHKERR(ierr)
    call VecSet(lvec, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(dm2d, gvec, ADD_VALUES, lvec, ierr) ;call CHKERR(ierr)
    call DMGlobalToLocalEnd  (dm2d, gvec, ADD_VALUES, lvec, ierr) ;call CHKERR(ierr)
    call DMRestoreGlobalVector(dm2d, gvec, ierr); call CHKERR(ierr)

    call getVecPointer(lvec, dm2d, x1d, x3d)

    do idof = i0, Ndof-i1
      call box_filter(x3d(idof,:,:), kernel_width)
    enddo

    arr(:,:,:) = x3d(:, &
      lbound(x3d,2)+kernel_width:ubound(x3d,2)-kernel_width, &
      lbound(x3d,3)+kernel_width:ubound(x3d,3)-kernel_width)

    call restoreVecPointer(lvec, x1d, x3d)
    call DMRestoreLocalVector(dm2d, lvec, ierr) ;call CHKERR(ierr)
  end subroutine

  subroutine box_filter(im, kernel_width)
    real(kind=ireals),intent(inout) :: im(:,:)
    integer(kind=iintegers) :: kernel_width

    real(kind=ireals),dimension(lbound(im,1):ubound(im,1),lbound(im,2):ubound(im,2) ) :: tmp
    real(kind=ireals),dimension(lbound(im,2):ubound(im,2),lbound(im,1):ubound(im,1) ) :: tmp_T,im_T

    integer(iintegers), parameter :: iter=1

    im_T = transpose(im)
    call box_blur_1d(im_T,kernel_width,iter,tmp_T)
    tmp = transpose(tmp_T)
    call box_blur_1d(tmp,kernel_width,iter,im)
  end subroutine

  subroutine box_blur_1d(im,k,iter,img)
    real(kind=ireals),intent(in) :: im(:,:)
    real(kind=ireals),dimension(lbound(im,1):ubound(im,1),lbound(im,2):ubound(im,2)),intent(out) :: img
    integer(kind=iintegers),intent(in) ::  k
    integer(kind=iintegers),intent(in) :: iter

    real(kind=ireals),dimension(lbound(im,1):ubound(im,1),lbound(im,2):ubound(im,2)) :: tmp
    real(kind=ireals),parameter :: c=1
    integer(kind=iintegers) :: inc,l,r,ki,x,i,lb,ub

    logical,parameter :: conv_mirror=.True.


    img(:,:) = im(:,:)

    do x=1,iter

      ! first do initialize average with full kernel
      lb=lbound(tmp,2)
      ub=ubound(tmp,2)
      tmp(:,lb) = 0
      tmp(:,lb) = tmp(:,lb) + img(:,lb) ! first the '0'

      ! then walk left
      inc = -1
      l = lb

      do ki=1,k
        l = l+inc
        if(conv_mirror) then
          if(l.eq.lb-1) then
            l   = lb
            inc = 1
          endif
          if(l.eq.ub+1) then
            l   = ub
            inc = -1
          endif
        else !cyclic
          if(l.eq.lb-1) l = ub
        endif
        tmp(:,lb) = tmp(:,lb) + img(:,l)
      enddo

      ! and walk right
      r = lb
      inc = 1
      do ki=1,k
        r = r+inc
        if(conv_mirror) then
          if(r.eq.ub+1) then
            r = ub
            inc = -1
          endif
          if(r.eq.lb-1) then
            r = lb
            inc = 1
          endif
        else !cyclic
          if(r.eq.ub+1) r = lb
        endif
        tmp(:,lb) = tmp(:,lb) + img(:,r)
      enddo

      ! now we have initialized first entry, can now go on using just the two new
      ! entries

      do i=lb+1,ub

        l = i-1
        inc = -1
        do ki=1,k
          l = l+inc
          if(conv_mirror) then
            if(l.eq.lb-1) then
              l = lb
              inc = 1
            endif
            if(l.eq.ub+1) then
              l = ub
              inc = -1
            endif
          else !cyclic
            if(l.eq.lb-1) l = ub
          endif
        enddo !ki

        r = i
        inc = 1
        do ki=1,k
          r = r+inc
          if(conv_mirror) then
            if(r.eq.ub+1) then
              r = ub
              inc = -1
            endif
            if(r.eq.lb-1) then
              r = lb
              inc = 1
            endif
          else !cyclic
            if(r.eq.ub+1) r = lb
          endif
        enddo !ki

        tmp(:,i) = tmp(:,i-1) + img(:,r) - img(:,l)

      enddo !i

      img = tmp
    enddo !iter

    img = img * (c/(2*k+1))**iter

  end subroutine
end module
