module m_xdmf_export
  use m_data_parameters, only: &
    & iintegers, &
    & ireals, &
    & mpiint, &
    & default_str_len, &
    & i0, i1, i2

  use m_helper_functions, only: &
    & CHKERR, &
    & CHKWARN, &
    & toStr, &
    & ind_1d_to_nd, &
    & get_arg

  use m_pprts_base, only: &
    & atmk, &
    & t_solver

  use m_buildings, only: &
    & t_pprts_buildings

  use m_boxmc_geometry, only: &
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

    real(ireals), allocatable :: xv(:, :, :, :)
    integer(iintegers) :: zs, zm, xs, xm, ys, ym
    integer(iintegers) :: k, i, j

    integer(iintegers) :: m, idx(4), l, n
    real(ireals) :: verts(3, 4)

    character(len=default_str_len) :: fname
    integer :: funit
    logical :: file_exists
    integer(mpiint) :: irank, numnodes
    integer(mpiint) :: iter

    ierr = 0

    fname = trim(fbasename)//'.xmf'

    iter = 0
99  continue
    if (solver%myid .eq. 0 .and. get_arg(.false., verbose)) print *, 'Writing xmf for buildings to ', trim(fname)

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    inquire (file=fname, exist=file_exists)
    if (file_exists) then
      iter = iter + 1
      fname = trim(fbasename)//'.'//toStr(iter)//'.xmf'
      if (solver%myid .eq. 0 .and. get_arg(.false., verbose)) &
        & call CHKWARN(1_mpiint, "adding suffix to output because file already exists: "//trim(fname))
      goto 99
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    if (solver%myid .eq. 0) then
      open (newunit=funit, file=fname, status='new', action='write')
      call write_header()
      close (funit)
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    associate (C => solver%C_diff, atm => solver%atm)
      allocate (xv(i0:i2, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye))
      do j = C%gys, C%gye
        do i = C%gxs, C%gxe
          do k = C%zs, C%ze
            xv(i0, k, i, j) = atm%hhl(i0, atmk(atm, k), i, j)
            xv(i1, k, i, j) = (real(i, ireals) + 0.5_ireals) * atm%dx
            xv(i2, k, i, j) = (real(j, ireals) + 0.5_ireals) * atm%dy
          end do
        end do
      end do
      zs = lbound(xv, 2); xs = lbound(xv, 3); ys = lbound(xv, 4)
      zm = size(xv, 2); xm = size(xv, 3); ym = size(xv, 4)
    end associate

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
  subroutine xdmf_pprts_srfc_flux(solver, fbasename, edn, eup, ierr, edir, verbose)
    class(t_solver), intent(in) :: solver
    character(len=*), intent(in) :: fbasename
    real(ireals), allocatable, dimension(:, :, :), intent(in) :: edn, eup
    integer(mpiint), intent(out) :: ierr
    real(ireals), allocatable, dimension(:, :, :), intent(in), optional :: edir
    logical, optional, intent(in) :: verbose

    real(ireals), allocatable :: xv(:, :, :, :)
    integer(iintegers) :: zs, zm, xs, xm, ys, ym
    integer(iintegers) :: k, i, j

    character(len=default_str_len) :: fname
    integer :: funit
    logical :: file_exists
    integer(mpiint) :: irank, numnodes
    integer(mpiint) :: iter

    ierr = 0

    fname = trim(fbasename)//'.xmf'

    iter = 0
99  continue
    if (solver%myid .eq. 0 .and. get_arg(.false., verbose)) print *, 'Writing xmf for surface fluxes to ', trim(fname)

    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)
    inquire (file=fname, exist=file_exists)
    if (file_exists) then
      iter = iter + 1
      fname = trim(fbasename)//'.'//toStr(iter)//'.xmf'
      if (solver%myid .eq. 0 .and. get_arg(.false., verbose)) &
        & call CHKWARN(1_mpiint, "adding suffix to output because file already exists: "//trim(fname))
      goto 99
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    if (solver%myid .eq. 0) then
      open (newunit=funit, file=fname, status='new', action='write')
      call write_header()
      close (funit)
    end if
    call mpi_barrier(solver%comm, ierr); call CHKERR(ierr)

    call mpi_comm_size(solver%comm, numnodes, ierr); call CHKERR(ierr)

    associate (C => solver%C_diff, atm => solver%atm)
      allocate (xv(i0:i2, C%zs:C%ze, C%gxs:C%gxe, C%gys:C%gye))
      do j = C%gys, C%gye
        do i = C%gxs, C%gxe
          do k = C%zs, C%ze
            xv(i0, k, i, j) = atm%hhl(i0, atmk(atm, k), i, j)
            xv(i1, k, i, j) = (real(i, ireals) + 0.5_ireals) * atm%dx
            xv(i2, k, i, j) = (real(j, ireals) + 0.5_ireals) * atm%dy
          end do
        end do
      end do
      zs = lbound(xv, 2); xs = lbound(xv, 3); ys = lbound(xv, 4)
      zm = size(xv, 2); xm = size(xv, 3); ym = size(xv, 4)
    end associate

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
      integer(iintegers) :: p, q, Ndim2, Ndim3
      real(ireals) :: x0, y0, z0

      Ndim2 = size(edn, dim=2)
      Ndim3 = size(edn, dim=3)
      ! Surface height from ground level (k=zs is surface, k=ze is top of atm)
      ! Surface is at ze (k=ze, hhl=0); zs is the top of the local domain
      ! Use C%xs/C%ys (non-ghost) — ghost cells have uninitialized hhl
      z0 = xv(i0, zs + zm - 1, solver%C_diff%xs, solver%C_diff%ys)
      ! Node origin: left/bottom edge of first local cell
      x0 = real(solver%C_diff%xs, ireals) * solver%atm%dx
      y0 = real(solver%C_diff%ys, ireals) * solver%atm%dy

      write (funit, *) '<Grid Name="GroundSubMesh'//toStr(solver%myid)//'">'
      ! 2DSMesh with explicit XYZ nodes — unambiguous in both VisIt and ParaView
      write (funit, *) '<Topology TopologyType="2DSMesh" NumberOfElements="', Ndim2 + 1, Ndim3 + 1, '"/>'
      write (funit, *) '<Geometry GeometryType="XYZ">'
      write (funit, *) '<DataItem Format="XML" Dimensions="', (Ndim2 + 1) * (Ndim3 + 1), 3, '">'
      do q = 0, Ndim3
        do p = 0, Ndim2
          write (funit, *) x0 + real(p, ireals) * solver%atm%dx, &
            & y0 + real(q, ireals) * solver%atm%dy, z0
        end do
      end do
      write (funit, *) '</DataItem>'
      write (funit, *) '</Geometry>'

      ! write data attributes
      ! edir
      if (present(edir)) then
        if (allocated(edir)) then
          write (funit, *) '<Attribute Center="Cell" Name="edir">'
          write (funit, *) '<DataItem Format="XML" Dimensions="', Ndim2, Ndim3, '">'
          write (funit, *) edir(size(edir, dim=1), :, :)
          write (funit, *) '</DataItem>', '</Attribute>'
        end if
      end if

      ! edn
      if (allocated(edn)) then
        write (funit, *) '<Attribute Center="Cell" Name="edn">'
        write (funit, *) '<DataItem Format="XML" Dimensions="', Ndim2, Ndim3, '">'
        write (funit, *) edn(size(edn, dim=1), :, :)
        write (funit, *) '</DataItem>', '</Attribute>'
      end if

      ! eup
      if (allocated(eup)) then
        write (funit, *) '<Attribute Center="Cell" Name="eup">'
        write (funit, *) '<DataItem Format="XML" Dimensions="', Ndim2, Ndim3, '">'
        write (funit, *) eup(size(eup, dim=1), :, :)
        write (funit, *) '</DataItem>', '</Attribute>'
      end if

      write (funit, *) '</Grid>'
    end subroutine
  end subroutine xdmf_pprts_srfc_flux
end module
