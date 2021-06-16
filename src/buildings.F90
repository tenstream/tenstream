module m_buildings
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & iintegers, &
    & ireals, &
    & mpiint, &
    & default_str_len

  use m_helper_functions, only: &
    & CHKERR, &
    & CHKWARN, &
    & toStr, &
    & ind_1d_to_nd, ind_nd_to_1d, ndarray_offsets, &
    & get_arg, &
    & deallocate_allocatable

  use m_pprts_base, only: &
    & t_solver

  implicit none

  private
  public :: &
    & buildingid2str,                     &
    & check_buildings_consistency,        &
    & clone_buildings,                    &
    & destroy_buildings,                  &
    & dump_pprts_buildings,               &
    & faceidx_by_cell_plus_offset,        &
    & init_buildings,                     &
    & t_plex_buildings,                   &
    & t_pprts_buildings,                  &
    & PPRTS_TOP_FACE,                     &
    & PPRTS_BOT_FACE,                     &
    & PPRTS_LEFT_FACE,                    &
    & PPRTS_RIGHT_FACE,                   &
    & PPRTS_REAR_FACE,                    &
    & PPRTS_FRONT_FACE

  type :: t_buildings_face_data
    ! iface_data is, depending on the implementing class:
    !  * local face indices in pprts mesh
    !  * face indices dmplex range [fStart, fEnd-1]
    ! subclasses will have a pointer to this memory
    ! so that it can be shared across multiple albedo types
    integer(iintegers), pointer :: data(:)
    integer(iintegers) :: ref_count = 0 ! data will be freed if reference count is 0
  end type

  type, abstract :: t_buildings
    type(t_buildings_face_data), pointer :: iface_data
    integer(iintegers), pointer :: iface(:) => null()

    ! albedo on a face, dim(Nbuilding_faces)
    real(ireals), allocatable :: albedo(:)

    ! face temperature to compute monochromatic planck emissivity, dim(Nbuilding_faces)
    ! should only be allocated if thermal computations are required
    real(ireals), allocatable :: temp(:)

    ! monochromatic planck emissivity on a face, dim(Nbuilding_faces)
    ! should only be allocated if thermal computations are required
    real(ireals), allocatable :: planck(:)

    ! irradiances on surfaces, if allocated -> dim(Nbuilding_faces)
    real(ireals), allocatable :: edir(:), incoming(:), outgoing(:)
  end type

  type, extends(t_buildings) :: t_pprts_buildings
    ! local face indices in pprts mesh,
    ! use faceidx_by_cell_plus_offset to translate from cell idx + face_ids (from PPRTS_XXX_FACE)

    ! offsets from iface idx sets to da indices, da_size e.g. [6, zm, xm, ym]
    integer(iintegers), allocatable :: da_offsets(:)
  end type

  type, extends(t_buildings) :: t_plex_buildings
    !integer(iintegers), pointer :: iface(:) => null() ! face indices dim(Nbuilding_faces) valid in range [fStart, fEnd-1]
  end type

  integer(iintegers), parameter :: &
    & PPRTS_TOP_FACE  =1, & ! z+0
    & PPRTS_BOT_FACE  =2, & ! z+1
    & PPRTS_LEFT_FACE =3, & ! x+0
    & PPRTS_RIGHT_FACE=4, & ! x+1
    & PPRTS_REAR_FACE =5, & ! y+0
    & PPRTS_FRONT_FACE=6    ! y+1

  interface destroy_buildings
    module procedure destroy_pprts_buildings
    module procedure destroy_plex_buildings
  end interface

contains

  subroutine init_buildings(buildings, da_sizes, Nfaces, ierr)
    type(t_pprts_buildings), intent(inout), allocatable :: buildings
    integer(iintegers), intent(in) :: da_sizes(:)
    integer(iintegers), intent(in) :: Nfaces
    integer(mpiint), intent(out) :: ierr
    ierr = 0
    if(allocated(buildings)) then
      ierr = int(size(buildings%iface)-Nfaces, mpiint)
      call CHKERR(ierr, 'buildings struct already allocated with different size ('//toStr(size(buildings%iface)))
      ierr = int(size(buildings%da_offsets)-size(da_sizes), mpiint)
      call CHKERR(ierr, 'buildings struct already allocated with different size of da_offsets '//&
        & ' '//toStr(size(buildings%da_offsets))//' vs '//toStr(size(da_sizes)))
    else
      allocate(buildings)
      allocate(buildings%iface_data)
      allocate(buildings%iface_data%data(Nfaces))
      allocate(buildings%da_offsets(size(da_sizes)))
      allocate(buildings%albedo(Nfaces))
    endif
    buildings%iface_data%data = -1_iintegers
    buildings%albedo = -1._ireals

    buildings%iface => buildings%iface_data%data
    buildings%iface_data%ref_count = 1

    call ndarray_offsets(da_sizes, buildings%da_offsets)
  end subroutine

  !> @brief clone a buildings struct
  !> @details iface_data is shared with reference counter on original buildings
  !> \n       data arrays(e.g. albedo etc) will be allocated on its own
  !> \n       and the user is responsible to set the values
  !> \n       if l_copy_data is True, additional data arrays are also copied
  subroutine clone_buildings(buildings, buildings_clone, l_copy_data, ierr)
    type(t_pprts_buildings), intent(in) :: buildings
    type(t_pprts_buildings), intent(inout), allocatable :: buildings_clone
    logical, intent(in) :: l_copy_data
    integer(mpiint), intent(out) :: ierr
    ierr = 0

    if(allocated(buildings_clone)) then
      ierr = int(size(buildings%iface)-size(buildings_clone%iface), mpiint)
      call CHKERR(ierr, 'clone already allocated but with a non-matching size faces'//&
        & '(orig:'//toStr(size(buildings%iface))//'.vs.clone:'//toStr(size(buildings_clone%iface)) )
      ierr = int(size(buildings%da_offsets)-size(buildings_clone%da_offsets), mpiint)
      call CHKERR(ierr, 'clone already allocated but with a non-matching size offsets'//&
        & ' '//toStr(size(buildings%da_offsets))//' vs '//toStr(size(buildings_clone%da_offsets)))
    else
      allocate(buildings_clone)
      if(allocated(buildings%da_offsets)) then
        allocate(buildings_clone%da_offsets(size(buildings%da_offsets)), source=buildings%da_offsets)
      endif
    endif

    associate( B => buildings, C => buildings_clone )
      C%iface      => B%iface
      C%iface_data => B%iface_data
      C%iface_data%ref_count = C%iface_data%ref_count+1

      if(l_copy_data) then
        !call copy_var(B%edir    , C%edir)
        !call copy_var(B%incoming, C%incoming)
        !call copy_var(B%outgoing, C%outgoing)
        call copy_var(B%albedo  , C%albedo)
        call copy_var(B%temp    , C%temp)
        call copy_var(B%planck  , C%planck)
      endif
    end associate
    contains
      subroutine copy_var(var_in, var_clone) ! allocate with same size of if already allocated, check that size matches
        real(ireals), intent(in), allocatable :: var_in(:)
        real(ireals), intent(inout), allocatable :: var_clone(:)
        if(allocated(var_in)) then
          if(.not.allocated(var_clone)) then
            allocate(var_clone(size(var_in)))
          else
            call CHKERR(int(size(var_in)-size(var_clone), mpiint), &
              & 'wrong size in clone ('//toStr(size(var_in))//'/'//toStr(size(var_clone))//')')
          endif
          var_clone = var_in
        else
          if(allocated(var_clone)) deallocate(var_clone)
        endif
      end subroutine
  end subroutine


  !> @brief destroy a buildings struct
  !> @details deallocates memory in this buildings struct
  !> \n       and frees shared iface data if reference count is 0
  subroutine destroy_pprts_buildings(buildings, ierr)
    type(t_pprts_buildings), allocatable, intent(inout) :: buildings
    integer(mpiint), intent(out) :: ierr
    ierr = 0

    if(.not.allocated(buildings)) ierr = 1
    call CHKERR(ierr, 'buildings not yet allocated')

    associate( B => buildings )

      nullify(B%iface)
      B%iface_data%ref_count = B%iface_data%ref_count -1

      if(B%iface_data%ref_count.eq.0) then
        deallocate(B%iface_data%data)
        nullify(B%iface_data%data)
      endif

      call deallocate_allocatable(B%da_offsets)
      call deallocate_allocatable(B%edir)
      call deallocate_allocatable(B%incoming)
      call deallocate_allocatable(B%outgoing)

      call deallocate_allocatable(B%albedo)
      call deallocate_allocatable(B%temp)
      call deallocate_allocatable(B%planck)
    end associate
    deallocate(buildings)
  end subroutine
  subroutine destroy_plex_buildings(buildings, ierr)
    type(t_plex_buildings), allocatable, intent(inout) :: buildings
    integer(mpiint), intent(out) :: ierr
    ierr = 0

    if(.not.allocated(buildings)) ierr = 1
    call CHKERR(ierr, 'buildings not yet allocated')

    associate( B => buildings )

      nullify(B%iface)
      B%iface_data%ref_count = B%iface_data%ref_count -1

      if(B%iface_data%ref_count.eq.0) then
        deallocate(B%iface_data%data)
        nullify(B%iface_data%data)
      endif

      call deallocate_allocatable(B%albedo)
      call deallocate_allocatable(B%temp)
      call deallocate_allocatable(B%planck)
    end associate
    deallocate(buildings)
  end subroutine

  !> @brief: find face index in a pprts dmda
  !> @details: this a helper routine to determine a flattened index in a 3D DMDA domain as in pprts
  !> \n        k,i,j are cell indices in 3D starting to count at 1
  !> \n        face_id has to be one of [PPRTS_TOP_FACE, ... ,PPRTS_FRONT_FACE]
  function faceidx_by_cell_plus_offset(da_offsets, k, i, j, face_id) result(fidx)
    integer(iintegers), intent(in) :: da_offsets(:)
    integer(iintegers), intent(in) :: k, i, j, face_id
    integer(iintegers) :: fidx
    fidx = ind_nd_to_1d(da_offsets, [face_id, k, i, j])
  end function

  subroutine check_buildings_consistency(buildings, Nz, Nx, Ny, ierr)
    type(t_pprts_buildings), intent(in) :: buildings
    integer(iintegers) :: Nx, Ny, Nz ! domain size, i.e. cell count
    integer(mpiint), intent(out) :: ierr

    logical, allocatable :: sides(:,:,:,:)
    integer(iintegers) :: m, idx(4), i,j,k

    ierr = 0

    call CHKERR(int(buildings%da_offsets(2)-6_iintegers, mpiint), &
      & 'da_offsets(2) should be 6 as in 6 sides but found '//toStr(buildings%da_offsets(2)))
    call CHKERR(int(buildings%da_offsets(3)-6*Nz, mpiint), &
      & 'da_offsets(3) should be 6*Nz but found '//toStr(buildings%da_offsets(3))// &
      & ' instead of '//toStr(Nz*6))

    call CHKERR(int(buildings%da_offsets(4)-6*Nz*Nx, mpiint), &
      & 'da_offsets(3) should be 6*Nz*Nx but found '//toStr(buildings%da_offsets(4))// &
      & ' instead of '//toStr(6*Nz*Nx))

    allocate(sides(6, Nz, Nx, Ny))
    sides = .False.

    do m = 1, size(buildings%iface)
      call ind_1d_to_nd(buildings%da_offsets, buildings%iface(m), idx)
      associate(d => idx(1), k => idx(2), i => idx(3), j => idx(4))
        sides(d,k,i,j) = .True.
      end associate
    enddo

    do j = 1, Ny
      do i = 1, Nx
        do k = 1, Nz
          if(any(sides(:,k,i,j))) then
            if(.not.all(sides(:,k,i,j))) then
              ierr = ierr + 1
              print *,'in cell k,i,j:', k,i,j, &
                & ' we found an semi open box'// new_line('')// &
                & ' - this may work for direct radiation or rayli results'// &
                & ' but diffuse radiation in the TenStream is probably wrong!', &
                & ' sides:', sides(:,k,i,j)
            endif
          endif
        enddo
      enddo
    enddo
  end subroutine


  !> @brief: dump the buildings information as xdmf
  !> @details: basename of the file will be expanded by .xmf postfix
  !> \n        if you are wondering, we use duplicates of vertices because it is easy
  subroutine dump_pprts_buildings(solver, buildings, fbasename, ierr, verbose)
    class(t_solver), intent(in) :: solver
    type(t_pprts_buildings), intent(inout) :: buildings
    character(len=*), intent(in) :: fbasename
    integer(mpiint), intent(out) :: ierr
    logical, optional, intent(in) :: verbose

    type(tDM) :: coordDA
    type(tVec) :: coordinates
    real(ireals), pointer, dimension(:,:,:,:) :: xv  =>null()
    real(ireals), pointer, dimension(:)       :: xv1d=>null()
    integer(iintegers) :: zs,zm, xs,xm, ys,ym

    integer(iintegers) :: m, idx(4), l, n
    real(ireals) :: verts(3,4)

    character(len=default_str_len) :: fname
    integer :: funit
    logical :: file_exists
    integer(mpiint) :: irank, numnodes

    ierr = 0

    fname = trim(fbasename)//'.xmf'

    if(solver%myid.eq.0 .and. get_arg(.False., verbose)) print *,'Writing xmf for buildings to ', trim(fname)

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    inquire(file=fname, exist=file_exists)
    if(file_exists) then
      call CHKWARN(1_mpiint, "skipping output because file already exists: "//trim(fname))
      return
    endif
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    if (solver%myid.eq.0) then
      open(newunit=funit, file=fname, status='new', action='write')
      call write_header()
      close(funit)
    endif
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    do irank = 0, numnodes-1
      if(irank.eq.solver%myid .and. size(buildings%iface).gt.0) then
        open(newunit=funit, file=fname, status='old', action='write', position='append')
        call write_grid()
        close(funit)
      endif
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    enddo

    if (solver%myid.eq.0) then
      open(newunit=funit, file=fname, status='old', action='write', position='append')
      call write_footer()
    endif
    close(funit)

    contains
      subroutine write_header()
        write (funit,FMT="(A)")'<?xml version="1.0" encoding="utf-8"?>'
        write (funit,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write (funit,*)'<Xdmf Version="3.0">'
        write (funit,*)'<Domain>'

        ! surface ground mesh
        write (funit,*)'<Grid Name="GroundMesh" GridType="Uniform">'
        write (funit,*)'<Topology TopologyType="3DCORECTMesh" NumberOfElements="1 ',&
          & solver%C_one%glob_xm+1, solver%C_one%glob_ym+1,'"/>'
        write (funit,*)'<Geometry GeometryType="ORIGIN_DXDYDZ">'
        write (funit,*)'<DataStructure Name="Origin" Dimensions="3" Format="XML"> 0 0 0'
        write (funit,*)'</DataStructure>'
        write (funit,*)'<DataStructure Name="Origin" Dimensions="3" Format="XML">'
        write (funit,*) solver%atm%dx, solver%atm%dy, 0
        write (funit,*)'</DataStructure>'
        write (funit,*)'</Geometry>'
        write (funit,*)'</Grid>'

        write (funit,*)'<Grid Name="BuildingsQuads" GridType="Collection">'
      end subroutine
      subroutine write_footer()
        write (funit,*) '</Grid>'
        write (funit,*) '</Domain>'
        write (funit,FMT="(A)") '</Xdmf>'
      end subroutine

      subroutine write_grid ()
        associate(               &
            & B => buildings,    &
            & C => solver%C_diff )

          call DMGetCoordinateDM(C%da, coordDA, ierr); call CHKERR(ierr)
          call DMDAGetGhostCorners(coordDA, zs, xs, ys, zm, xm, ym, ierr); call CHKERR(ierr)

          call DMGetCoordinatesLocal(C%da, coordinates, ierr); call CHKERR(ierr)
          call VecGetArrayF90(coordinates, xv1d, ierr); call CHKERR(ierr)
          xv(0:2, zs:zs+zm-1 ,xs:xs+xm-1 ,ys:ys+ym-1) => xv1d

          write (funit,*) '<Grid Name="Quads'//toStr(solver%myid)//'">'
          write (funit,*) '<Topology TopologyType="Quadrilateral" NumberOfElements="'//toStr(size(B%iface))//'">'
          write (funit,*) '<DataItem Format="XML" DataType="Int" '//&
            & 'Dimensions="'//toStr(size(B%iface))//' 4">'

          n = 0
          do m = 1, size(B%iface)
            write (funit,*) ( l, l = (m-1)*4, m*4-1 )
          enddo

          write (funit,*) '</DataItem>'
          write (funit,*) '</Topology>'

          write (funit,*) '<Geometry GeometryType="XYZ">'
          write (funit,*) '<DataItem Format="XML" '//&
            & 'Dimensions="', size(B%iface), 4, 3,'">'

          do m = 1, size(B%iface)

            call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
            idx(2:4) = idx(2:4) -1 + [C%zs, C%xs, C%ys]

            associate(k => idx(2), i => idx(3), j => idx(4))

              select case(idx(1))
              case(PPRTS_TOP_FACE)
                verts(:,1) = xv(:,k,i,j)
                verts(:,2) = xv(:,k,i+1,j)
                verts(:,3) = xv(:,k,i+1,j+1)
                verts(:,4) = xv(:,k,i,j+1)
              case(PPRTS_BOT_FACE)
                verts(:,1) = xv(:,k+1,i,j)
                verts(:,2) = xv(:,k+1,i+1,j)
                verts(:,3) = xv(:,k+1,i+1,j+1)
                verts(:,4) = xv(:,k+1,i,j+1)
              case(PPRTS_LEFT_FACE)
                verts(:,1) = xv(:,k  ,i,j)
                verts(:,2) = xv(:,k  ,i,j+1)
                verts(:,3) = xv(:,k+1,i,j+1)
                verts(:,4) = xv(:,k+1,i,j)
              case(PPRTS_RIGHT_FACE)
                verts(:,1) = xv(:,k  ,i+1,j)
                verts(:,2) = xv(:,k  ,i+1,j+1)
                verts(:,3) = xv(:,k+1,i+1,j+1)
                verts(:,4) = xv(:,k+1,i+1,j)
              case(PPRTS_REAR_FACE )
                verts(:,1) = xv(:,k  ,i  ,j)
                verts(:,2) = xv(:,k  ,i+1,j)
                verts(:,3) = xv(:,k+1,i+1,j)
                verts(:,4) = xv(:,k+1,i  ,j)
              case(PPRTS_FRONT_FACE)
                verts(:,1) = xv(:,k  ,i  ,j+1)
                verts(:,2) = xv(:,k  ,i+1,j+1)
                verts(:,3) = xv(:,k+1,i+1,j+1)
                verts(:,4) = xv(:,k+1,i  ,j+1)
              case default
                call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1)+1))
              end select

              do l = 1, size(verts,2)
                write (funit,*) verts([2,3,1],l) - [real(ireals) :: solver%atm%dx/2, solver%atm%dy/2, 0]
              enddo
            end associate
          enddo

          write (funit,*) '</DataItem>'
          write (funit,*) '</Geometry>'

          ! write data attributes

          ! planck
          if(allocated(B%planck)) then
            write (funit,*) '<Attribute Center="Cell" Name="planck">'
            write (funit,*) '<DataItem Format="XML" '//&
              & 'Dimensions="', size(B%iface), '">'
            write (funit,*) B%planck
            write (funit,*) '</DataItem>'
            write (funit,*) '</Attribute>'
          endif

          ! temp
          if(allocated(B%temp)) then
            write (funit,*) '<Attribute Center="Cell" Name="temp">'
            write (funit,*) '<DataItem Format="XML" '//&
              & 'Dimensions="', size(B%iface), '">'
            write (funit,*) B%temp
            write (funit,*) '</DataItem>'
            write (funit,*) '</Attribute>'
          endif

          ! albedo
          if(allocated(b%albedo)) then
            write (funit,*) '<Attribute Center="Cell" Name="albedo">'
            write (funit,*) '<DataItem Format="XML" '//&
              & 'Dimensions="', size(B%iface), '">'
            write (funit,*) B%albedo
            write (funit,*) '</DataItem>'
            write (funit,*) '</Attribute>'
          endif

          ! incoming
          if(allocated(B%incoming)) then
            write (funit,*) '<Attribute Center="Cell" Name="incoming">'
            write (funit,*) '<DataItem Format="XML" '//&
              & 'Dimensions="', size(B%iface), '">'
            write (funit,*) B%incoming
            write (funit,*) '</DataItem>'
            write (funit,*) '</Attribute>'
          endif

          ! outgoing
          if(allocated(B%outgoing)) then
            write (funit,*) '<Attribute Center="Cell" Name="outgoing">'
            write (funit,*) '<DataItem Format="XML" '//&
              & 'Dimensions="', size(B%iface), '">'
            write (funit,*) B%outgoing
            write (funit,*) '</DataItem>'
            write (funit,*) '</Attribute>'
          endif

          ! edir
          if(allocated(B%edir)) then
            write (funit,*) '<Attribute Center="Cell" Name="edir">'
            write (funit,*) '<DataItem Format="XML" '//&
              & 'Dimensions="', size(B%iface), '">'
            write (funit,*) B%edir
            write (funit,*) '</DataItem>'
            write (funit,*) '</Attribute>'
          endif

          write (funit,*) '</Grid>'

        end associate
      end subroutine
  end subroutine

  function buildingid2str(bid) result(str)
    integer(iintegers), intent(in) :: bid
    character(len=16) :: str
    select case(bid)
    case(PPRTS_TOP_FACE)
      str = "PPRTS_TOP_FACE"
    case(PPRTS_BOT_FACE)
      str = "PPRTS_BOT_FACE"
    case(PPRTS_LEFT_FACE)
      str = "PPRTS_LEFT_FACE"
    case(PPRTS_RIGHT_FACE)
      str = "PPRTS_RIGHT_FACE"
    case(PPRTS_REAR_FACE )
      str = "PPRTS_REAR_FACE"
    case(PPRTS_FRONT_FACE)
      str = "PPRTS_FRONT_FACE"
    case default
      call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(bid))
    end select
  end function
end module
