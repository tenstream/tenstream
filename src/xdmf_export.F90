module m_xdmf_export
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
    & get_arg

  use m_pprts_base, only: &
    & set_dmda_cell_coordinates, &
    & t_solver

  use m_buildings, only: &
    & t_pprts_buildings, &
    & PPRTS_TOP_FACE, &
    & PPRTS_BOT_FACE, &
    & PPRTS_FRONT_FACE, &
    & PPRTS_LEFT_FACE, &
    & PPRTS_REAR_FACE, &
    & PPRTS_RIGHT_FACE

  implicit none

  private
  public :: &
    & xdmf_pprts_buildings, &
    & xdmf_pprts_srfc_flux

contains

  !> @brief: dump the buildings information as xdmf
  !> @details: basename of the file will be expanded by .xmf postfix
  !> \n        if you are wondering, we use duplicates of vertices because it is easy
  subroutine xdmf_pprts_buildings(solver, buildings, fbasename, ierr, verbose)
    class(t_solver), intent(in) :: solver
    type(t_pprts_buildings), intent(inout) :: buildings
    character(len=*), intent(in) :: fbasename
    integer(mpiint), intent(out) :: ierr
    logical, optional, intent(in) :: verbose

    type(tDM) :: coordDA
    type(tVec) :: coordinates
    real(ireals), pointer, dimension(:, :, :, :) :: xv => null()
    real(ireals), pointer, dimension(:) :: xv1d => null()
    integer(iintegers) :: zs, zm, xs, xm, ys, ym

    integer(iintegers) :: m, idx(4), l, n
    real(ireals) :: verts(3, 4)

    character(len=default_str_len) :: fname
    integer :: funit
    logical :: file_exists
    integer(mpiint) :: irank, numnodes

    ierr = 0

    fname = trim(fbasename)//'.xmf'

    if (solver%myid .eq. 0 .and. get_arg(.false., verbose)) print *, 'Writing xmf for buildings to ', trim(fname)

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    inquire (file=fname, exist=file_exists)
    if (file_exists) then
      call CHKWARN(1_mpiint, "skipping output because file already exists: "//trim(fname))
      return
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    if (solver%myid .eq. 0) then
      open (newunit=funit, file=fname, status='new', action='write')
      call write_header()
      close (funit)
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    do irank = 0, numnodes - 1
      if (irank .eq. solver%myid .and. size(buildings%iface) .gt. 0) then
        open (newunit=funit, file=fname, status='old', action='write', position='append')
        call write_grid()
        close (funit)
      end if
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end do

    if (solver%myid .eq. 0) then
      open (newunit=funit, file=fname, status='old', action='write', position='append')
      call write_footer()
    end if
    close (funit)

  contains
    subroutine write_header()
      write (funit, FMT="(A)") '<?xml version="1.0" encoding="utf-8"?>'
      write (funit, *) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write (funit, *) '<Xdmf Version="3.0">'
      write (funit, *) '<Domain>'

      write (funit, *) '<Grid Name="BuildingsQuads" GridType="Collection">'
    end subroutine

    subroutine write_footer()
      write (funit, *) '</Grid>'
      write (funit, *) '</Domain>'
      write (funit, FMT="(A)") '</Xdmf>'
    end subroutine

    subroutine write_grid()
      associate (               &
          & B => buildings,    &
          & C => solver%C_diff)

        call DMGetCoordinateDM(C%da, coordDA, ierr); call CHKERR(ierr)
        call DMDAGetGhostCorners(coordDA, zs, xs, ys, zm, xm, ym, ierr); call CHKERR(ierr)

        call DMGetCoordinatesLocal(C%da, coordinates, ierr); call CHKERR(ierr)
        if (coordinates .eq. PETSC_NULL_VEC) then
          call set_dmda_cell_coordinates(solver, solver%atm, C%da, ierr)
          call DMGetCoordinatesLocal(C%da, coordinates, ierr); call CHKERR(ierr)
        end if
        call VecGetArrayF90(coordinates, xv1d, ierr); call CHKERR(ierr)
        xv(0:2, zs:zs + zm - 1, xs:xs + xm - 1, ys:ys + ym - 1) => xv1d

        write (funit, *) '<Grid Name="Quads'//toStr(solver%myid)//'">'
        write (funit, *) '<Topology TopologyType="Quadrilateral" NumberOfElements="'//toStr(size(B%iface))//'">'
        write (funit, *) '<DataItem Format="XML" DataType="Int" '//&
          & 'Dimensions="'//toStr(size(B%iface))//' 4">'

        n = 0
        do m = 1, size(B%iface)
          write (funit, *) (l, l=(m - 1) * 4, m * 4 - 1)
        end do

        write (funit, *) '</DataItem>'
        write (funit, *) '</Topology>'

        write (funit, *) '<Geometry GeometryType="XYZ">'
        write (funit, *) '<DataItem Format="XML" '//&
          & 'Dimensions="', size(B%iface), 4, 3, '">'

        do m = 1, size(B%iface)

          call ind_1d_to_nd(B%da_offsets, B%iface(m), idx)
          idx(2:4) = idx(2:4) - 1 + [C%zs, C%xs, C%ys]

          associate (k => idx(2), i => idx(3), j => idx(4))

            select case (idx(1))
            case (PPRTS_TOP_FACE)
              verts(:, 1) = xv(:, k, i, j)
              verts(:, 2) = xv(:, k, i + 1, j)
              verts(:, 3) = xv(:, k, i + 1, j + 1)
              verts(:, 4) = xv(:, k, i, j + 1)
            case (PPRTS_BOT_FACE)
              verts(:, 1) = xv(:, k + 1, i, j)
              verts(:, 2) = xv(:, k + 1, i + 1, j)
              verts(:, 3) = xv(:, k + 1, i + 1, j + 1)
              verts(:, 4) = xv(:, k + 1, i, j + 1)
            case (PPRTS_LEFT_FACE)
              verts(:, 1) = xv(:, k, i, j)
              verts(:, 2) = xv(:, k, i, j + 1)
              verts(:, 3) = xv(:, k + 1, i, j + 1)
              verts(:, 4) = xv(:, k + 1, i, j)
            case (PPRTS_RIGHT_FACE)
              verts(:, 1) = xv(:, k, i + 1, j)
              verts(:, 2) = xv(:, k, i + 1, j + 1)
              verts(:, 3) = xv(:, k + 1, i + 1, j + 1)
              verts(:, 4) = xv(:, k + 1, i + 1, j)
            case (PPRTS_REAR_FACE)
              verts(:, 1) = xv(:, k, i, j)
              verts(:, 2) = xv(:, k, i + 1, j)
              verts(:, 3) = xv(:, k + 1, i + 1, j)
              verts(:, 4) = xv(:, k + 1, i, j)
            case (PPRTS_FRONT_FACE)
              verts(:, 1) = xv(:, k, i, j + 1)
              verts(:, 2) = xv(:, k, i + 1, j + 1)
              verts(:, 3) = xv(:, k + 1, i + 1, j + 1)
              verts(:, 4) = xv(:, k + 1, i, j + 1)
            case default
              call CHKERR(1_mpiint, 'unknown building face_idx '//toStr(idx(1) + 1))
            end select

            do l = 1, size(verts, 2)
              write (funit, *) verts([2, 3, 1], l) - [real(ireals) :: solver%atm%dx, solver%atm%dy, 0]
            end do
          end associate
        end do

        write (funit, *) '</DataItem>'
        write (funit, *) '</Geometry>'

        ! write data attributes

        ! planck
        if (allocated(B%planck)) then
          write (funit, *) '<Attribute Center="Cell" Name="planck">'
          write (funit, *) '<DataItem Format="XML" '//&
            & 'Dimensions="', size(B%iface), '">'
          write (funit, *) B%planck
          write (funit, *) '</DataItem>'
          write (funit, *) '</Attribute>'
        end if

        ! temp
        if (allocated(B%temp)) then
          write (funit, *) '<Attribute Center="Cell" Name="temp">'
          write (funit, *) '<DataItem Format="XML" '//&
            & 'Dimensions="', size(B%iface), '">'
          write (funit, *) B%temp
          write (funit, *) '</DataItem>'
          write (funit, *) '</Attribute>'
        end if

        ! albedo
        if (allocated(b%albedo)) then
          write (funit, *) '<Attribute Center="Cell" Name="albedo">'
          write (funit, *) '<DataItem Format="XML" '//&
            & 'Dimensions="', size(B%iface), '">'
          write (funit, *) B%albedo
          write (funit, *) '</DataItem>'
          write (funit, *) '</Attribute>'
        end if

        ! incoming
        if (allocated(B%incoming)) then
          write (funit, *) '<Attribute Center="Cell" Name="incoming">'
          write (funit, *) '<DataItem Format="XML" '//&
            & 'Dimensions="', size(B%iface), '">'
          write (funit, *) B%incoming
          write (funit, *) '</DataItem>'
          write (funit, *) '</Attribute>'
        end if

        ! outgoing
        if (allocated(B%outgoing)) then
          write (funit, *) '<Attribute Center="Cell" Name="outgoing">'
          write (funit, *) '<DataItem Format="XML" '//&
            & 'Dimensions="', size(B%iface), '">'
          write (funit, *) B%outgoing
          write (funit, *) '</DataItem>'
          write (funit, *) '</Attribute>'
        end if

        ! edir
        if (allocated(B%edir)) then
          write (funit, *) '<Attribute Center="Cell" Name="edir">'
          write (funit, *) '<DataItem Format="XML" '//&
            & 'Dimensions="', size(B%iface), '">'
          write (funit, *) B%edir
          write (funit, *) '</DataItem>'
          write (funit, *) '</Attribute>'
        end if

        write (funit, *) '</Grid>'
      end associate
    end subroutine
  end subroutine xdmf_pprts_buildings

  !> @brief: dump the surface flux information as xdmf
  !> @details: basename of the file will be expanded by .xmf postfix
  subroutine xdmf_pprts_srfc_flux(solver, fbasename, edir, edn, eup, ierr, verbose)
    class(t_solver), intent(in) :: solver
    character(len=*), intent(in) :: fbasename
    real(ireals), allocatable, dimension(:, :, :), intent(in) :: edir, edn, eup
    integer(mpiint), intent(out) :: ierr
    logical, optional, intent(in) :: verbose

    type(tDM) :: coordDA
    type(tVec) :: coordinates
    real(ireals), pointer, dimension(:, :, :, :) :: xv => null()
    real(ireals), pointer, dimension(:) :: xv1d => null()
    integer(iintegers) :: zs, zm, xs, xm, ys, ym

    character(len=default_str_len) :: fname
    integer :: funit
    logical :: file_exists
    integer(mpiint) :: irank, numnodes

    ierr = 0

    fname = trim(fbasename)//'.xmf'

    if (solver%myid .eq. 0 .and. get_arg(.false., verbose)) print *, 'Writing xmf for surface fluxes to ', trim(fname)

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    inquire (file=fname, exist=file_exists)
    if (file_exists) then
      call CHKWARN(1_mpiint, "skipping output because file already exists: "//trim(fname))
      return
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    if (solver%myid .eq. 0) then
      open (newunit=funit, file=fname, status='new', action='write')
      call write_header()
      close (funit)
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    do irank = 0, numnodes - 1
      if (irank .eq. solver%myid) then
        open (newunit=funit, file=fname, status='old', action='write', position='append')
        call write_grid()
        close (funit)
      end if
      call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    end do

    if (solver%myid .eq. 0) then
      open (newunit=funit, file=fname, status='old', action='write', position='append')
      call write_footer()
    end if
    close (funit)

  contains
    subroutine write_header()
      write (funit, FMT="(A)") '<?xml version="1.0" encoding="utf-8"?>'
      write (funit, *) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
      write (funit, *) '<Xdmf Version="3.0">'
      write (funit, *) '<Domain>'

      ! surface ground mesh
      !write (funit,*)'<Grid Name="GroundMesh" GridType="Uniform">'
      !write (funit,*)'<Topology TopologyType="3DCORECTMesh" NumberOfElements="1 ',&
      !  & solver%C_one%glob_xm+1, solver%C_one%glob_ym+1,'"/>'
      !write (funit,*)'<Geometry GeometryType="ORIGIN_DXDYDZ">'
      !write (funit,*)'<DataStructure Name="Origin" Dimensions="3" Format="XML"> 0 0 0'
      !write (funit,*)'</DataStructure>'
      !write (funit,*)'<DataStructure Name="Origin" Dimensions="3" Format="XML">'
      !write (funit,*) solver%atm%dx, solver%atm%dy, 0
      !write (funit,*)'</DataStructure>'
      !write (funit,*)'</Geometry>'
      !write (funit,*)'</Grid>'

      write (funit, *) '<Grid Name="SurfaceMesh" GridType="Collection">'
    end subroutine

    subroutine write_footer()
      write (funit, *) '</Grid>'
      write (funit, *) '</Domain>'
      write (funit, FMT="(A)") '</Xdmf>'
    end subroutine

    subroutine write_grid()
      associate (               &
          & C => solver%C_diff)

        call DMGetCoordinateDM(C%da, coordDA, ierr); call CHKERR(ierr)
        call DMDAGetGhostCorners(coordDA, zs, xs, ys, zm, xm, ym, ierr); call CHKERR(ierr)

        call DMGetCoordinatesLocal(C%da, coordinates, ierr); call CHKERR(ierr)
        call VecGetArrayF90(coordinates, xv1d, ierr); call CHKERR(ierr)
        xv(0:2, zs:zs + zm - 1, xs:xs + xm - 1, ys:ys + ym - 1) => xv1d

        write (funit, *) '<Grid Name="GroundSubMesh'//toStr(solver%myid)//'">'
        write (funit, *) '<Topology TopologyType="3DCORECTMesh" &
          & NumberOfElements=" 1 ', size(edn, dim=3) + 1, size(edn, dim=2) + 1, '"/>'
        write (funit, *) '<Geometry GeometryType="ORIGIN_DXDYDZ">'
        write (funit, *) '<DataStructure Name="Origin" Dimensions="3" Format="XML">'
        write (funit, *) xv([1, 2, 0], zs + zm - 1, xs, ys)! - [real(ireals) :: solver%atm%dx/2, solver%atm%dy/2, 0]
        write (funit, *) '</DataStructure>'
        write (funit, *) '<DataStructure Name="Spacing" Dimensions="3" Format="XML">'
        write (funit, *) solver%atm%dx, solver%atm%dy, 0
        write (funit, *) '</DataStructure>'
        write (funit, *) '</Geometry>'

        ! write data attributes
        ! edir
        if (allocated(edir)) then
          write (funit, *) '<Attribute Center="Cell" Name="edir">'
          write (funit, *) '<DataItem Format="XML" Dimensions="', size(edir, dim=2), size(edir, dim=3), '">'
          write (funit, *) edir(size(edir, dim=1), :, :)
          write (funit, *) '</DataItem>', '</Attribute>'
        end if

        ! edn
        if (allocated(edn)) then
          write (funit, *) '<Attribute Center="Cell" Name="edn">'
          write (funit, *) '<DataItem Format="XML" Dimensions="', size(edn, dim=2), size(edn, dim=3), '">'
          write (funit, *) edn(size(edn, dim=1), :, :)
          write (funit, *) '</DataItem>', '</Attribute>'
        end if

        ! eup
        if (allocated(eup)) then
          write (funit, *) '<Attribute Center="Cell" Name="eup">'
          write (funit, *) '<DataItem Format="XML" Dimensions="', size(eup, dim=2), size(eup, dim=3), '">'
          write (funit, *) eup(size(eup, dim=1), :, :)
          write (funit, *) '</DataItem>', '</Attribute>'
        end if

        write (funit, *) '</Grid>'
      end associate
    end subroutine
  end subroutine xdmf_pprts_srfc_flux
end module
