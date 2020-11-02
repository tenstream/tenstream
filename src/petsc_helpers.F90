module m_petsc_helpers
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    zero, i0, i1, i2, i3, default_str_len, init_mpi_data_parameters

  use m_helper_functions, only : get_arg, CHKERR, char_to_upper

  implicit none

  private
  public :: petscGlobalVecToZero, scatterZerotoPetscGlobal, &
    petscGlobalVecToAll, &
    petscVecToF90, &
    f90VecToPetsc, &
    getVecPointer, restoreVecPointer, &
    dmda_convolve_ediff_srfc, &
    hegedus_trick, &
    gen_shared_subcomm, gen_shared_scatter_ctx, &
    is_local_vec

  interface f90VecToPetsc
    module procedure f90VecToPetsc_2d, f90VecToPetsc_3d, f90VecToPetsc_4d
  end interface
  interface petscVecToF90
    module procedure petscVecToF90_2d, petscVecToF90_3d, petscVecToF90_4d
  end interface
  interface getVecPointer
    module procedure getVecPointer_2d, getVecPointer_3d
  end interface
  interface restoreVecPointer
    module procedure restoreVecPointer_2d, restoreVecPointer_3d
  end interface

  logical, parameter :: ldebug=.False.

  ! in earlier versions, we had problems with compilers
  ! and therefore have a rewrite here, set to false to use it
  logical, parameter :: lgetVecPointer_use_petsc_func=.True.

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

  !> @brief Scatter a petsc global vector into a local vector with global size
  subroutine petscGlobalVecToAll(gVec, dm, lVec)
    type(tVec), intent(in)    :: gVec
    type(tDM), intent(in)     :: dm
    type(tVec), intent(out)   :: lVec

    type(tVec) :: natural
    type(tVecScatter) :: scatter_context

    integer(mpiint) :: ierr

    call DMDACreateNaturalVector(dm, natural, ierr); call CHKERR(ierr)

    call DMDAGlobalToNaturalBegin(dm, gVec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)
    call DMDAGlobalToNaturalEnd  (dm, gVec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)

    call VecScatterCreateToAll(natural, scatter_context, lVec, ierr); call CHKERR(ierr)

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
  subroutine petscVecToF90_4d(vec, dm, arr, only_on_rank0)
    type(tVec), intent(in)    :: vec
    type(tDM), intent(in)     :: dm
    real(ireals), intent(inout), allocatable :: arr(:,:,:,:)
    logical, intent(in), optional :: only_on_rank0
    logical :: l_only_on_rank0, l_has_global_dimensions

    integer(iintegers) :: vecsize
    VecType :: vtype
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: comm, myid, ierr
    integer(iintegers) :: zs, xs, ys, zm, xm, ym
    integer(iintegers) :: dims(4)

    integer(iintegers) :: dmdim, dof, glob_xm, glob_ym, glob_zm
    integer(iintegers) :: nprocz, nprocx, nprocy
    integer(iintegers) :: stencil_width, stencil_type
    integer(iintegers) :: boundary_z, boundary_x, boundary_y


    if(allocated(arr)) stop 'You shall not call petscVecToF90 with an already allocated array!'

    call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    l_only_on_rank0 = get_arg(.False., only_on_rank0)

    if(l_only_on_rank0 .and. myid.ne.0) &
      call CHKERR(myid, 'Only rank 0 should call the routine petscVecToF90 with only_on_rank0=.T.')

    call DMDAGetInfo(dm, dmdim,             &
      & glob_zm, glob_xm, glob_ym,          &
      & nprocz, nprocx, nprocy,             &
      & dof, stencil_width,                 &
      & boundary_z, boundary_x, boundary_y, &
      & stencil_type, ierr) ;call CHKERR(ierr)

    call DMDAGetCorners(dm, zs, xs, ys, zm, xm, ym, ierr) ;call CHKERR(ierr)

    l_has_global_dimensions = l_only_on_rank0
    call VecGetType(vec, vtype, ierr); call CHKERR(ierr)
    if(vtype.eq.VECSEQ) l_has_global_dimensions = .True.

    dims(:) = [dof, zm, xm, ym]
    if(l_has_global_dimensions) dims(:) = [dof, glob_zm, glob_xm, glob_ym]

    if(.not.allocated(arr)) allocate(arr(dims(1), dims(2), dims(3), dims(4)))

    call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
    if(vecsize.ne.size(arr)) then
      print *,'petscVecToF90 Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      stop 'petscVecToF90 Vecsizes dont match!'
    endif

    if(.not.l_has_global_dimensions) then
      call getVecPointer(dm, vec, x1d, x4d)
      arr = x4d
      call restoreVecPointer(dm, vec, x1d, x4d)
    else
      call VecGetArrayF90(vec,x1d,ierr); call CHKERR(ierr)
      arr = reshape( x1d, dims )
      call VecRestoreArrayF90(vec,x1d,ierr); call CHKERR(ierr)
    endif
  end subroutine
  subroutine petscVecToF90_3d(vec, dm, arr, only_on_rank0)
    type(tVec), intent(in)    :: vec
    type(tDM), intent(in)     :: dm
    real(ireals), intent(inout), allocatable :: arr(:,:,:)
    logical, intent(in), optional :: only_on_rank0
    logical :: l_only_on_rank0, l_has_global_dimensions

    VecType :: vtype
    integer(iintegers) :: vecsize
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: comm, myid, ierr
    integer(iintegers) :: zs, xs, ys, zm, xm, ym
    integer(iintegers) :: dims(3)

    integer(iintegers) :: dmdim, dof, glob_xm, glob_ym, glob_zm
    integer(iintegers) :: nprocz, nprocx, nprocy
    integer(iintegers) :: stencil_width, stencil_type
    integer(iintegers) :: boundary_z, boundary_x, boundary_y

    l_only_on_rank0 = get_arg(.False., only_on_rank0)

    call DMDAGetInfo(dm, dmdim,             &
      & glob_zm, glob_xm, glob_ym,          &
      & nprocz, nprocx, nprocy,             &
      & dof, stencil_width,                 &
      & boundary_z, boundary_x, boundary_y, &
      & stencil_type, ierr) ;call CHKERR(ierr)

    call DMDAGetCorners(dm, zs, xs, ys, zm, xm, ym, ierr) ;call CHKERR(ierr)

    if(dof.ne.1) &
      call CHKERR(1_mpiint, 'petscVecToF90_3d should only be called with anything else than DM%dof of 1')

    call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)


    if(l_only_on_rank0 .and. myid.ne.0) &
      call CHKERR(myid, 'Only rank 0 should call the routine petscVecToF90 with opt_l_only_on_rank0=.T.')

    l_has_global_dimensions = l_only_on_rank0
    call VecGetType(vec, vtype, ierr); call CHKERR(ierr)
    if(vtype.eq.VECSEQ) l_has_global_dimensions = .True.

    dims(:) = [zm, xm, ym]
    if(l_has_global_dimensions) dims(:) = [glob_zm, glob_xm, glob_ym]

    if(.not.allocated(arr)) allocate(arr(dims(1), dims(2), dims(3)))

    call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
    if(vecsize.ne.size(arr)) then
      print *,'petscVecToF90 Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      call CHKERR(1_mpiint, 'petscVecToF90 Vecsizes dont match!')
    endif

    if(.not.l_has_global_dimensions) then
      call getVecPointer(dm, vec, x1d, x4d)
      arr = x4d(i0,:,:,:)
      call restoreVecPointer(dm, vec, x1d, x4d)
    else
      call VecGetArrayF90(vec,x1d,ierr); call CHKERR(ierr)
      arr = reshape( x1d, dims )
      call VecRestoreArrayF90(vec,x1d,ierr); call CHKERR(ierr)
    endif
  end subroutine
  subroutine petscVecToF90_2d(vec, dm, arr, only_on_rank0)
    type(tVec), intent(in)    :: vec
    type(tDM), intent(in)     :: dm
    real(ireals), intent(inout), allocatable :: arr(:,:)
    logical, intent(in), optional :: only_on_rank0
    logical :: l_only_on_rank0, l_has_global_dimensions

    VecType :: vtype
    integer(iintegers) :: vecsize
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: comm, myid, ierr
    integer(iintegers) :: zs, xs, ys, zm, xm, ym
    integer(iintegers) :: dims(2)

    integer(iintegers) :: dmdim, dof, glob_xm, glob_ym, glob_zm
    integer(iintegers) :: nprocz, nprocx, nprocy
    integer(iintegers) :: stencil_width, stencil_type
    integer(iintegers) :: boundary_z, boundary_x, boundary_y

    l_only_on_rank0 = get_arg(.False., only_on_rank0)

    call DMDAGetInfo(dm, dmdim,             &
      & glob_zm, glob_xm, glob_ym,          &
      & nprocz, nprocx, nprocy,             &
      & dof, stencil_width,                 &
      & boundary_z, boundary_x, boundary_y, &
      & stencil_type, ierr) ;call CHKERR(ierr)

    call DMDAGetCorners(dm, zs, xs, ys, zm, xm, ym, ierr) ;call CHKERR(ierr)

    if(dof.ne.1) &
      call CHKERR(1_mpiint, 'petscVecToF90_3d should only be called with anything else than DM%dof of 1')

    call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)


    if(l_only_on_rank0 .and. myid.ne.0) &
      call CHKERR(myid, 'Only rank 0 should call the routine petscVecToF90 with opt_l_only_on_rank0=.T.')

    l_has_global_dimensions = l_only_on_rank0
    call VecGetType(vec, vtype, ierr); call CHKERR(ierr)
    if(vtype.eq.VECSEQ) l_has_global_dimensions = .True.

    dims(:) = [xm, ym]
    if(l_has_global_dimensions) dims(:) = [glob_xm, glob_ym]

    if(.not.allocated(arr)) allocate(arr(dims(1), dims(2)))

    call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
    if(vecsize.ne.size(arr)) then
      print *,'petscVecToF90 Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      call CHKERR(1_mpiint, 'petscVecToF90 Vecsizes dont match!')
    endif

    if(.not.l_has_global_dimensions) then
      call getVecPointer(dm, vec, x1d, x4d)
      arr = x4d(i0,i0,:,:)
      call restoreVecPointer(dm, vec, x1d, x4d)
    else
      call VecGetArrayF90(vec,x1d,ierr); call CHKERR(ierr)
      arr = reshape( x1d, dims )
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

    call getVecPointer(dm, vec, x1d, x4d)
    x4d = arr
    call restoreVecPointer(dm, vec, x1d, x4d)

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

    call getVecPointer(dm, vec, x1d, x4d)
    x4d(i0,:,:,:) = arr
    call restoreVecPointer(dm, vec, x1d, x4d)

    call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
    if(vecsize.ne.size(arr)) then
      print *,'f90VecToPetsc Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      call CHKERR(1_mpiint, 'f90VecToPetsc Vecsizes dont match!')
    endif
  end subroutine
  subroutine f90VecToPetsc_2d(arr, dm, vec)
    real(ireals), intent(in)  :: arr(:,:)
    type(tDM), intent(in)     :: dm
    type(tVec), intent(inout) :: vec

    integer(iintegers) :: vecsize
    real(ireals),pointer :: x1d(:)=>null(),x4d(:,:,:,:)=>null()

    integer(mpiint) :: ierr

    call getVecPointer(dm, vec, x1d, x4d)
    x4d(i0,i0,:,:) = arr
    call restoreVecPointer(dm, vec, x1d, x4d)

    call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
    if(vecsize.ne.size(arr)) then
      print *,'f90VecToPetsc Vecsizes dont match! petsc:', vecsize, 'f90 arr', size(arr)
      call CHKERR(1_mpiint, 'f90VecToPetsc Vecsizes dont match!')
    endif
  end subroutine

  subroutine getVecPointer_3d(dm,vec,x1d,x4d,readonly)
    type(tVec) :: vec
    type(tDM), intent(in) :: dm
    real(ireals),intent(inout),pointer,dimension(:,:,:,:) :: x4d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d
    logical, optional, intent(in) :: readonly

    integer(iintegers) :: dmdim, dof, glob_xm, glob_ym, glob_zm
    integer(iintegers) :: nprocz, nprocx, nprocy
    integer(iintegers) :: stencil_width, stencil_type
    integer(iintegers) :: boundary_z, boundary_x, boundary_y

    integer(iintegers) :: N
    integer(iintegers) :: zs, ze, xs, xe, ys, ye, zm, xm, ym
    integer(iintegers) :: gzs, gze, gxs, gxe, gys, gye, gzm, gxm, gym
    integer(mpiint) :: comm, myid, ierr
    logical :: lghosted

    if(associated(x1d).or.associated(x4d)) then
      print *,'ERROR : getVecPointer : input vector already associated!!',associated(x1d),associated(x4d)
      call CHKERR(1_mpiint, 'getVecPointer : input vector already associated')
    endif

    if(lgetVecPointer_use_petsc_func) then
      if(get_arg(.False.,readonly)) then
        call DMDAVecGetArrayReadF90(dm, vec, x4d, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call DMDAVecGetArrayF90(dm, vec, x4d, ierr); call CHKERR(ierr)
        call VecGetArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif
    else
      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, '-show_getvecpointerdm', ierr); call CHKERR(ierr)

      call DMDAGetInfo(dm, dmdim,             &
        & glob_zm, glob_xm, glob_ym,          &
        & nprocz, nprocx, nprocy,             &
        & dof, stencil_width,                 &
        & boundary_z, boundary_x, boundary_y, &
        & stencil_type, ierr) ;call CHKERR(ierr)

      call DMDAGetCorners(dm, zs, xs, ys, zm, xm, ym, ierr) ;call CHKERR(ierr)
      call DMDAGetGhostCorners(dm,gzs,gxs,gys,gzm,gxm,gym,ierr) ;call CHKERR(ierr)
      xe = xs+xm-1
      ye = ys+ym-1
      ze = zs+zm-1
      gxe = gxs+gxm-1
      gye = gys+gym-1
      gze = gzs+gzm-1

      call VecGetLocalSize(vec,N,ierr)

      if(N .eq. dof*xm*ym*zm) then
        lghosted=.False.
      else if( N .eq. dof*gxm*gym*gzm ) then
        lghosted=.True.
      else
        call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
        print *,myid,'Size N:', N, dof*xm*ym*zm, dof*glob_xm*glob_ym*glob_zm, dof*gxm*gym*gzm
        call CHKERR(1_mpiint, 'Local Vector dimensions do not conform to DMDA size')
        stop 'Local Vector dimensions do not conform to DMDA size'
      endif

      if(get_arg(.False.,readonly)) then
        call VecGetArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call VecGetArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif

      if(lghosted) then
        x4d(0:dof-1 , gzs:gze, gxs:gxe , gys:gye ) => x1d
      else
        x4d(0:dof-1 , zs:ze  , xs:xe   , ys:ye   ) => x1d
      endif
    endif

  end subroutine

  subroutine restoreVecPointer_3d(dm,vec,x1d,x4d,readonly)
    type(tDM), intent(in) :: dm
    type(tVec) :: vec
    real(ireals),intent(inout),pointer,dimension(:,:,:,:) :: x4d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d
    logical, optional, intent(in) :: readonly
    integer(mpiint) :: ierr

    if(.not.associated(x1d).or..not.associated(x4d)) then
      print *,'ERROR : restoreVecPointer : input vector not yet associated!!',associated(x1d),associated(x4d)
      call CHKERR(1_mpiint, 'input vector not yet associated')
    endif

    if(lgetVecPointer_use_petsc_func) then
      if(get_arg(.False.,readonly)) then
        call DMDAVecRestoreArrayReadF90(dm, vec, x4d, ierr); call CHKERR(ierr)
        call VecRestoreArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call DMDAVecRestoreArrayF90(dm, vec, x4d, ierr); call CHKERR(ierr)
        call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif
    else
      x4d => null()
      if(get_arg(.False.,readonly)) then
        call VecRestoreArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif
    endif
  end subroutine

  subroutine getVecPointer_2d(dm,vec,x1d,x3d, readonly)
    type(tDM), intent(in) :: dm
    type(tVec), intent(in) :: vec
    real(ireals),intent(inout),pointer,dimension(:,:,:) :: x3d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d
    logical, optional, intent(in) :: readonly

    integer(iintegers) :: N, xs, xe, ys, ye, xm, ym
    integer(iintegers) :: gxs, gxe, gys, gye, gxm, gym

    integer(iintegers) :: dmdim, dof, glob_xm, glob_ym, glob_zm
    integer(iintegers) :: nprocz, nprocx, nprocy
    integer(iintegers) :: stencil_width, stencil_type
    integer(iintegers) :: boundary_z, boundary_x, boundary_y

    integer(mpiint) :: comm, myid, ierr
    logical :: lghosted

    if(associated(x1d).or.associated(x3d)) then
      print *,'ERROR : getVecPointer : input vector already associated!!',associated(x1d),associated(x3d)
      call CHKERR(1_mpiint, 'getVecPointer : input vector already associated')
    endif

    if(lgetVecPointer_use_petsc_func) then
      if(get_arg(.False.,readonly)) then
        call DMDAVecGetArrayReadF90(dm, vec, x3d, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call DMDAVecGetArrayF90(dm, vec, x3d, ierr); call CHKERR(ierr)
        call VecGetArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif
    else

      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, '-show_getvecpointerdm', ierr); call CHKERR(ierr)
      call DMDAGetInfo(dm, dmdim,             &
        & glob_xm, glob_ym, glob_zm,          &
        & nprocx, nprocy, nprocz,             &
        & dof, stencil_width,                 &
        & boundary_x, boundary_y, boundary_z, &
        & stencil_type, ierr) ;call CHKERR(ierr)


      call DMDAGetCorners(dm, xs, ys, PETSC_NULL_INTEGER, xm, ym, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)
      call DMDAGetGhostCorners(dm,gxs,gys,PETSC_NULL_INTEGER,gxm,gym,PETSC_NULL_INTEGER,ierr) ;call CHKERR(ierr)
      xe = xs+xm-1
      ye = ys+ym-1
      gxe = gxs+gxm-1
      gye = gys+gym-1

      call VecGetLocalSize(vec,N,ierr)

      if(N .eq. dof*xm*ym) then
        lghosted=.False.
      else if( N .eq. dof*gxm*gym ) then
        lghosted=.True.
      else
        call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
        print *,myid,'Size N:', N, dof*xm*ym, dof*glob_xm*glob_ym, dof*gxm*gym
        call CHKERR(1_mpiint, 'Local Vector dimensions do not conform to DMDA size')
        stop 'Local Vector dimensions do not conform to DMDA size'
      endif

      if(get_arg(.False.,readonly)) then
        call VecGetArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call VecGetArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif
      if(lghosted) then
        x3d(0:dof-1 , gxs:gxe , gys:gye ) => x1d
      else
        x3d(0:dof-1 , xs:xe   , ys:ye   ) => x1d
      endif

    endif
  end subroutine
  subroutine restoreVecPointer_2d(dm, vec, x1d, x3d, readonly)
    type(tDM), intent(in) :: dm
    type(tVec), intent(in) :: vec
    real(ireals),intent(inout),pointer,dimension(:,:,:) :: x3d
    real(ireals),intent(inout),pointer,dimension(:) :: x1d
    logical, optional, intent(in) :: readonly
    integer(mpiint) :: ierr

    if(.not.associated(x1d).or..not.associated(x3d)) then
      print *,'ERROR : restoreVecPointer : input vector not yet associated!!',associated(x1d),associated(x3d)
      call exit(1)
    endif

    if(lgetVecPointer_use_petsc_func) then

      if(get_arg(.False.,readonly)) then
        call DMDAVecRestoreArrayReadF90(dm, vec, x3d, ierr); call CHKERR(ierr)
        call VecRestoreArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call DMDAVecRestoreArrayF90(dm, vec, x3d, ierr); call CHKERR(ierr)
        call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif

    else

      x3d => null()
      if(get_arg(.False.,readonly)) then
        call VecRestoreArrayReadF90(vec,x1d,ierr) ;call CHKERR(ierr)
      else
        call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
      endif
    endif
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

    integer(iintegers) :: dmdim, dof, glob_xm, glob_ym, glob_zm
    integer(iintegers) :: nprocz, nprocx, nprocy
    integer(iintegers) :: stencil_width, stencil_type
    integer(iintegers) :: boundary_z, boundary_x, boundary_y

    integer(iintegers), allocatable, dimension(:) :: Nxperproc, Nyperproc
    integer(iintegers) :: Ndof, idof
    integer(mpiint) :: comm, myid, numnodes, ierr

    Ndof=size(arr, dim=1, kind=iintegers)

    call PetscObjectGetComm(dm3d, comm, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call DMDAGetInfo(dm3d, dmdim,           &
      & glob_zm, glob_xm, glob_ym,          &
      & nprocz, nprocx, nprocy,             &
      & dof, stencil_width,                 &
      & boundary_z, boundary_x, boundary_y, &
      & stencil_type, ierr) ;call CHKERR(ierr)

    allocate(Nxperproc(nprocx), Nyperproc(nprocy))

    call DMDAGetOwnershipRanges(dm3d, PETSC_NULL_INTEGER, Nxperproc, Nyperproc, ierr); call CHKERR(ierr)

    if(kernel_width.gt.int(min(size(arr,dim=2),size(arr,dim=3)), iintegers)) then
      call CHKERR(int(kernel_width, mpiint), 'smoothing kernel size is bigger than local domains...'// &
        'this would need more work to be done... go for it...'// &
        '(cascaded uniform filters with updates inbetween or just call this function more often)')
    endif

    call DMDACreate2d(comm, &
      DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, &
      DMDA_STENCIL_BOX, glob_xm, glob_ym, &
      nprocx, nprocy, &
      Ndof, kernel_width, &
      Nxperproc, Nyperproc, &
      dm2d, ierr) ;call CHKERR(ierr)
    call DMSetup(dm2d, ierr); call CHKERR(ierr)

    call DMGetGlobalVector(dm2d, gvec, ierr); call CHKERR(ierr)
    call getVecPointer(dm2d, gvec, g1d, g3d)
    g3d(:,:,:) = arr(:,:,:)
    call restoreVecPointer(dm2d, gvec, g1d, g3d)

    call DMGetLocalVector(dm2d, lvec, ierr); call CHKERR(ierr)
    call VecSet(lvec, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(dm2d, gvec, ADD_VALUES, lvec, ierr) ;call CHKERR(ierr)
    call DMGlobalToLocalEnd  (dm2d, gvec, ADD_VALUES, lvec, ierr) ;call CHKERR(ierr)
    call DMRestoreGlobalVector(dm2d, gvec, ierr); call CHKERR(ierr)

    call getVecPointer(dm2d, lvec, x1d, x3d)

    do idof = i0, Ndof-i1
      call box_filter(x3d(idof,:,:), kernel_width)
    enddo

    arr(:,:,:) = x3d(:, &
      lbound(x3d,2)+kernel_width:ubound(x3d,2)-kernel_width, &
      lbound(x3d,3)+kernel_width:ubound(x3d,3)-kernel_width)

    call restoreVecPointer(dm2d, lvec, x1d, x3d)
    call DMRestoreLocalVector(dm2d, lvec, ierr) ;call CHKERR(ierr)
    call DMDestroy(dm2d, ierr); call CHKERR(ierr)
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

    img = img * (c/real(2*k+1, ireals))**iter

  end subroutine

  subroutine hegedus_trick(ksp, b, x)
    type(tksp), intent(in) :: ksp
    type(tVec), intent(in) :: b
    type(tVec), intent(inout) :: x

    type(tDM)  :: dm
    character(len=default_str_len) :: prefix
    type(tPC)  :: prec
    type(tMat) :: A, P
    type(tVec) :: Ax0, z
    real(ireals) :: znorm, norm
    logical :: lhegedus, lflg
    integer(mpiint) :: comm, myid, ierr

    call KSPGetOptionsPrefix(ksp, prefix, ierr); call CHKERR(ierr)
    if(prefix.eq.'solar_diff_') then
      lhegedus = .True.
    else
      lhegedus = .False.
    endif
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, prefix, &
      "-use_hegedus" , lhegedus , lflg , ierr) ;call CHKERR(ierr)

    if(lhegedus) then
      call KSPGetDM(ksp, dm, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(dm, Ax0, ierr); call CHKERR(ierr)
      call DMGetGlobalVector(dm, z  , ierr); call CHKERR(ierr)

      call KSPGetPC(ksp,prec,ierr); call CHKERR(ierr)
      call PCGetOperators(prec, A, P, ierr); call CHKERR(ierr)

      call MatMult(A, x, Ax0, ierr); call CHKERR(ierr)
      call MatMult(P, Ax0, z, ierr); call CHKERR(ierr)

      call VecDot(z, Ax0, znorm, ierr); call CHKERR(ierr)
      call VecDot(z, b, norm, ierr); call CHKERR(ierr)

      if(znorm.gt.epsilon(znorm)) then
        call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
        if(myid.eq.0.and.ldebug) print *,'hegedus_trick', norm, znorm, norm/znorm
        call VecScale(x, norm / znorm, ierr); call CHKERR(ierr)
      endif

      call DMRestoreglobalVector(dm, Ax0, ierr); call CHKERR(ierr)
      call DMRestoreglobalVector(dm, z  , ierr); call CHKERR(ierr)
    endif
  end subroutine

  subroutine gen_shared_subcomm(comm, subcomm, ierr)
    integer(mpiint), intent(in) :: comm
    integer(mpiint), intent(out):: subcomm, ierr
    integer(mpiint) :: mpinfo
    integer(mpiint) :: stype
    call init_mpi_data_parameters(comm)
    call select_split_type(stype, ierr); call CHKERR(ierr)
    call MPI_Info_create(mpinfo, ierr); call CHKERR(ierr)
    call MPI_Comm_split_type(comm, stype, 0_mpiint, mpinfo, subcomm, ierr); call CHKERR(ierr)
    call MPI_Info_free(mpinfo, ierr); call CHKERR(ierr)
  end subroutine

  subroutine select_split_type(stype, ierr)
    integer(mpiint), intent(out) :: stype
    integer(mpiint), intent(out) :: ierr
    character(len=default_str_len) :: csplit_type
    logical :: lflg
    ierr = 0
    csplit_type = 'MPI_COMM_TYPE_SHARED'
    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      & '-mpi_split_type', csplit_type, lflg, ierr) ; call CHKERR(ierr)
    call char_to_upper(csplit_type)
    select case(trim(csplit_type))
    case('MPI_COMM_TYPE_SHARED')
      stype = MPI_COMM_TYPE_SHARED
#ifdef HAVE_OMPI
    case('OMPI_COMM_TYPE_HWTHREAD')
      stype = OMPI_COMM_TYPE_HWTHREAD
    case('OMPI_COMM_TYPE_CORE')
      stype = OMPI_COMM_TYPE_CORE
    case('OMPI_COMM_TYPE_L1CACHE')
      stype = OMPI_COMM_TYPE_L1CACHE
    case('OMPI_COMM_TYPE_L2CACHE')
      stype = OMPI_COMM_TYPE_L2CACHE
    case('OMPI_COMM_TYPE_L3CACHE')
      stype = OMPI_COMM_TYPE_L3CACHE
    case('OMPI_COMM_TYPE_SOCKET')
      stype = OMPI_COMM_TYPE_SOCKET
    case('OMPI_COMM_TYPE_NUMA')
      stype = OMPI_COMM_TYPE_NUMA
    case('OMPI_COMM_TYPE_BOARD')
      stype = OMPI_COMM_TYPE_BOARD
    case('OMPI_COMM_TYPE_HOST')
      stype = OMPI_COMM_TYPE_HOST
    case('OMPI_COMM_TYPE_CU')
      stype = OMPI_COMM_TYPE_CU
#endif
    case default
      ierr = 1
      call CHKERR(1_mpiint, 'Dont know mpi_split_type `'//trim(csplit_type)//'`'//&
        & " Available options:"//new_line('')// &
        & " MPI_COMM_TYPE_SHARED"//new_line('')// &
#ifdef HAVE_OMPI
        & " OMPI_COMM_TYPE_HWTHREAD"//new_line('')// &
        & " OMPI_COMM_TYPE_CORE"//new_line('')// &
        & " OMPI_COMM_TYPE_L1CACHE"//new_line('')// &
        & " OMPI_COMM_TYPE_L2CACHE"//new_line('')// &
        & " OMPI_COMM_TYPE_L3CACHE"//new_line('')// &
        & " OMPI_COMM_TYPE_SOCKET"//new_line('')// &
        & " OMPI_COMM_TYPE_NUMA"//new_line('')// &
        & " OMPI_COMM_TYPE_BOARD"//new_line('')// &
        & " OMPI_COMM_TYPE_HOST"//new_line('')// &
        & " OMPI_COMM_TYPE_CU"//new_line('')// &
#endif
        & "")
    end select
  end subroutine

  subroutine gen_shared_scatter_ctx(gvec, svec, ctx, ierr, opt_is_in, opt_is_out)
    type(tVec), intent(in) :: gvec, svec
    type(tVecScatter), intent(out) :: ctx
    integer(mpiint), intent(out) :: ierr
    type(tIS), intent(in), optional :: opt_is_in, opt_is_out
    type(tIS) :: is_in, is_out

    integer(iintegers) :: N_out
    logical :: created_is_out

    created_is_out = .False.

    if(present(opt_is_out)) then
      is_out = opt_is_out
    else
        call VecGetSize(svec, N_out, ierr); call CHKERR(ierr)
        call ISCreateStride(PETSC_COMM_SELF, N_out, 0_iintegers, 1_iintegers, is_out, ierr); call CHKERR(ierr)
        created_is_out = .True.
    endif

    if(present(opt_is_in)) then
      is_in = opt_is_in
    else
      is_in = is_out
    endif

    call VecScatterCreate(gvec, is_in, svec, is_out, ctx, ierr); call CHKERR(ierr)
    if(created_is_out) then
      call ISDestroy(is_out, ierr); call CHKERR(ierr)
    endif
  end subroutine

  subroutine is_local_vec(da, vec, is_local, ierr)
    type(tDM), intent(in) :: da
    type(tVec), intent(in) :: vec
    logical, intent(out) :: is_local
    integer(mpiint), intent(out) :: ierr

    integer(iintegers) :: vsize, ncx, ncy, ncz, numCells

    call VecGetLocalSize(vec, vsize, ierr); call CHKERR(ierr)
    call DMDAGetNumCells(da, ncx, ncy, ncz, numCells, ierr); call CHKERR(ierr)

    is_local = vsize.eq.numCells
  end subroutine
end module
