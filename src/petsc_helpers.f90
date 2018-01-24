module m_petsc_helpers
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i3

  use m_helper_functions, only : get_arg, CHKERR

  implicit none

  private
  public :: petscGlobalVecToZero, scatterZerotoPetscGlobal, &
    petscVecToF90, &
    f90VecToPetsc, &
    getVecPointer, restoreVecPointer

  interface f90VecToPetsc
    module procedure f90VecToPetsc_3d, f90VecToPetsc_4d
  end interface
  interface petscVecToF90
    module procedure petscVecToF90_3d, petscVecToF90_4d
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
        stop 'f90VecToPetsc Vecsizes dont match!'
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

  subroutine getVecPointer(vec,dm,x1d,x4d)
  type(tVec) :: vec
  type(tDM), intent(in) :: dm
  PetscScalar,intent(inout),pointer,dimension(:,:,:,:) :: x4d
  PetscScalar,intent(inout),pointer,dimension(:) :: x1d

  integer(iintegers) :: N, dmdim, dof, glob_zm, glob_xm, glob_ym
  integer(iintegers) :: zs, ze, xs, xe, ys, ye, zm, xm, ym
  integer(iintegers) :: gzs, gze, gxs, gxe, gys, gye, gzm, gxm, gym
  integer(mpiint) :: comm, myid, ierr
  logical :: lghosted

  if(associated(x1d).or.associated(x4d)) then
    print *,'ERROR : getVecPointer : input vector already associated!!',associated(x1d),associated(x4d)
    call sleep(30)
    call exit(1)
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

  print *,'Domain Corners z:: ',zs,':',ze,' (',zm,' entries)','global size',glob_zm
  print *,'Domain Corners x:: ',xs,':',xe,' (',xm,' entries)','global size',glob_xm
  print *,'Domain Corners y:: ',ys,':',ye,' (',ym,' entries)','global size',glob_ym
  print *,'DOF',dof

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
  subroutine restoreVecPointer(vec,x1d,x4d)
  type(tVec) :: vec
  PetscScalar,intent(inout),pointer,dimension(:,:,:,:) :: x4d
  PetscScalar,intent(inout),pointer,dimension(:) :: x1d
  integer(mpiint) :: ierr

  if(.not.associated(x1d).or..not.associated(x4d)) then
    print *,'ERROR : restoreVecPointer : input vector not yet associated!!',associated(x1d),associated(x4d)
    call exit(1)
  endif

  x4d => null()
  call VecRestoreArrayF90(vec,x1d,ierr) ;call CHKERR(ierr)
  x1d => null()
  end subroutine

end module
