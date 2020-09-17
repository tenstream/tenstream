module m_buildings
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: iintegers, ireals, mpiint

  use m_helper_functions, only: &
    & CHKERR, &
    & toStr, &
    & ind_1d_to_nd, ind_nd_to_1d, ndarray_offsets, &
    & get_arg, &
    & deallocate_allocatable

  implicit none

  private
  public :: &
    & t_pprts_buildings,                  &
    & t_plex_buildings,                   &
    & init_buildings,                     &
    & clone_buildings,                    &
    & destroy_buildings,                  &
    & faceidx_by_cell_plus_offset,        &
    & check_buildings_consistency,        &
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

  integer, parameter :: &
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
    if(allocated(buildings_clone)) ierr = 1
    call CHKERR(ierr, 'buildings_clone already allocated')

    allocate(buildings_clone)
    associate( B => buildings, C => buildings_clone )
      C%iface      => B%iface
      C%iface_data => B%iface_data
      C%iface_data%ref_count = C%iface_data%ref_count+1

      if(l_copy_data) then
        if(allocated(B%da_offsets)) then
          allocate(C%da_offsets(size(B%da_offsets)), source=B%da_offsets)
        endif
        if(allocated(B%albedo)) then
          allocate(C%albedo(size(B%albedo)), source=B%albedo)
        endif
        if(allocated(B%planck)) then
          allocate(C%planck(size(B%planck)), source=B%planck)
        endif
      endif
    end associate
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
end module


module m_examples_pprts_buildings
  use m_data_parameters, only : &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & zero, one, pi, i1, i2, default_str_len

  use m_helper_functions, only : &
    & CHKERR, &
    & toStr, cstr, &
    & spherical_2_cartesian, rotate_angle_z, &
    & meanval, &
    & is_inrange

  use m_pprts, only : init_pprts, &
    & set_optical_properties, solve_pprts, &
    & pprts_get_result, set_angles, &
    & gather_all_toZero

  use m_pprts_base, only: t_solver, &
    & allocate_pprts_solver_from_commandline, destroy_pprts, &
    & t_solver_1_2, t_solver_3_6, t_solver_3_10, t_solver_3_16, &
    & t_solver_8_10, t_solver_8_16, t_solver_8_18

  use m_tenstream_options, only: read_commandline_options

  use m_buildings, only: t_pprts_buildings, &
    & faceidx_by_cell_plus_offset, &
    & check_buildings_consistency, &
    & PPRTS_TOP_FACE, PPRTS_BOT_FACE, &
    & PPRTS_LEFT_FACE, PPRTS_RIGHT_FACE, &
    & PPRTS_REAR_FACE, PPRTS_FRONT_FACE, &
    & init_buildings

  use m_netcdfio, only: ncwrite

  implicit none

contains
  subroutine ex_pprts_buildings(comm, lverbose, &
      & lthermal, lsolar, &
      & Nx, Ny, Nlay, icollapse, &
      & glob_box_i, glob_box_j, glob_box_k, &
      & box_albedo, box_planck, &
      & dx, dy, dz, &
      & incSolar, phi0, theta0, &
      & albedo, dtau, w0, &
      & gedir, gedn, geup, gabso, &
      & buildings, &
      & outfile )

    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lverbose, lthermal, lsolar
    integer(iintegers), intent(in) :: Nx, Ny, Nlay ! global domain size
    integer(iintegers), intent(in) :: icollapse    ! collapse upper nr of layers into 1 layer
    integer(iintegers), intent(in) :: glob_box_i, glob_box_j, glob_box_k ! global index of a single cube building
    real(ireals), intent(in) :: box_albedo         ! albedo of building faces
    real(ireals), intent(in) :: box_planck         ! planck emission of building faces (only used if lthermal=.True.)
    real(ireals), intent(in) :: dx, dy, dz         ! grid spacing in [m]
    real(ireals), intent(in) :: incSolar           ! solar constant at TOA [W/m2]
    real(ireals), intent(in) :: phi0, theta0       ! sun azimuth(phi) and zenith(theta) angle
    real(ireals), intent(in) :: albedo, dtau, w0   ! surface albedo, vertically integrated optical depth and constant single scattering albedo
    real(ireals), allocatable, dimension(:,:,:), intent(out) :: gedir, gedn, geup, gabso ! global arrays on rank 0
    type(t_pprts_buildings), allocatable, intent(inout) :: buildings
    character(len=*), intent(in), optional :: outfile ! output file to dump flux results

    real(ireals) :: dz1d(Nlay)
    real(ireals) :: sundir(3)
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,plck
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver
    integer :: Nbuildings
    logical :: lhave_box

    character(len=default_str_len) :: groups(2)

    integer(iintegers) :: k, i
    integer(iintegers) :: box_k, box_i, box_j
    integer(mpiint) :: myid, numnodes, ierr

    dz1d = dz

    call init_mpi_data_parameters(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    sundir = spherical_2_cartesian(phi0, theta0)
    call init_pprts(comm, Nlay, Nx, Ny, dx,dy, sundir, solver, dz1d, collapseindex=icollapse)

    associate(Ca => solver%C_one_atm, C1 => solver%C_one)
      allocate(kabs(Ca%zm  , Ca%xm, Ca%ym ))
      allocate(ksca(Ca%zm  , Ca%xm, Ca%ym ))
      allocate(g   (Ca%zm  , Ca%xm, Ca%ym ))

      if(lthermal) then
        allocate(plck(Ca%zm+1, Ca%xm, Ca%ym ))
        plck(:,:,:) = 0
        plck(Ca%zm+1,:,:) = 0
      endif

      kabs = dtau*(one-w0)/(dz*Nlay)
      ksca = dtau*w0/(dz*Nlay)
      g    = zero

      box_k = glob_box_k - C1%zs
      box_i = glob_box_i - C1%xs
      box_j = glob_box_j - C1%ys

      if(lverbose) print *, myid, 'Have box:', &
        & is_inrange(glob_box_k, C1%zs+1, C1%ze+1), &
        & is_inrange(glob_box_i, C1%xs+1, C1%xe+1), &
        & is_inrange(glob_box_j, C1%ys+1, C1%ye+1)

      if( &
        !& .False. .and. &
        & is_inrange(glob_box_k, C1%zs+1, C1%ze+1).and. &
        & is_inrange(glob_box_i, C1%xs+1, C1%xe+1).and. &
        & is_inrange(glob_box_j, C1%ys+1, C1%ye+1)      ) then

        lhave_box = .True.
        Nbuildings = 6
      else
        lhave_box = .False.
        Nbuildings = 0
      endif

      call init_buildings(buildings, &
        & [integer(iintegers) :: 6, C1%zm, C1%xm,  C1%ym], &
        & Nbuildings, &
        & ierr); call CHKERR(ierr)

      if(lthermal) allocate(buildings%planck(Nbuildings))
      do i=1,Nbuildings
        buildings%iface(i) = faceidx_by_cell_plus_offset( &
          & buildings%da_offsets, &
          & box_k, &
          & box_i, &
          & box_j, i)
        buildings%albedo(i) = box_albedo
        if(lthermal) buildings%planck(i) = box_planck
      enddo

      call check_buildings_consistency(buildings, C1%zm, C1%xm, C1%ym, ierr); call CHKERR(ierr)

      if(lthermal) then
        call set_optical_properties(solver, albedo, kabs, ksca, g, plck)
      else
        call set_optical_properties(solver, albedo, kabs, ksca, g)
      endif
      call set_angles(solver, sundir)

      call solve_pprts(solver, &
        & lthermal=lthermal, &
        & lsolar=lsolar, &
        & edirTOA=incSolar, &
        & opt_buildings=buildings)

      call pprts_get_result(solver, fdn, fup, fdiv, fdir, opt_buildings=buildings)

      if(lsolar) then
        call gather_all_toZero(solver%C_one_atm1, fdir, gedir)
      endif
      call gather_all_toZero(solver%C_one_atm1, fdn, gedn)
      call gather_all_toZero(solver%C_one_atm1, fup, geup)
      call gather_all_toZero(solver%C_one_atm, fdiv, gabso)

      if(myid.eq.0_mpiint.and.present(outfile)) then
        groups(1) = trim(outfile)
        if(lsolar) then
          groups(2) = 'edir'; call ncwrite(groups, gedir, ierr); call CHKERR(ierr)
        endif
        groups(2) = 'edn' ; call ncwrite(groups, gedn , ierr); call CHKERR(ierr)
        groups(2) = 'eup' ; call ncwrite(groups, geup , ierr); call CHKERR(ierr)
        groups(2) = 'abso'; call ncwrite(groups, gabso, ierr); call CHKERR(ierr)
      endif

      if(lverbose .and. myid.eq.0) then
        print *,'y-slice: i='//toStr(box_i)
        i = box_i
        if(lsolar) then
          do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
            print *, k, cstr('edir'//toStr( gedir(k, i, :) ), 'red')
          enddo
        endif ! lsolar
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' edn'//toStr( gedn(k, i, :) ), 'green')
        enddo
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' eup'//toStr( geup(k, i, :) ), 'blue')
        enddo
        do k = 1+solver%C_one%zs, 1+solver%C_one%ze
          print *, k, cstr('abso'//toStr( gabso(k, i, :) ), 'purple')
        enddo

        print *,''
        print *,'x-slice: j='//toStr(box_j)
        i = box_j
        if(lsolar) then
          do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
            print *, k, cstr('edir'//toStr( gedir(k, :, i) ), 'red')
          enddo
        endif
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' edn'//toStr( gedn(k, :, i) ), 'green')
        enddo
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' eup'//toStr( geup(k, :, i) ), 'blue')
        enddo
        do k = 1+solver%C_one%zs, 1+solver%C_one%ze
          print *, k, cstr('abso'//toStr( gabso(k, :, i) ), 'purple')
        enddo

        print *,''
        if(lsolar) then
          do k = lbound(fdiv,1), ubound(fdiv, 1)
            print *, k, 'mean ', &
              & 'edir', meanval(gedir(k,:,:)), &
              & 'edn' , meanval(gedn(k,:,:)), &
              & 'eup' , meanval(geup(k,:,:)), &
              & 'abso', meanval(gabso(k,:,:))
          enddo
          k = ubound(fdir, 1)
          print *, k, 'mean ', &
            & 'edir', meanval(gedir(k,:,:)), &
            & 'edn' , meanval(gedn(k,:,:)), &
            & 'eup' , meanval(geup(k,:,:))
        else
          do k = lbound(fdiv,1), ubound(fdiv, 1)
            print *, k, &
              & 'mean ', &
              & 'edn', meanval(gedn(k,:,:)), &
              & 'eup', meanval(geup(k,:,:)), &
              & 'abso', meanval(gabso(k,:,:))
          enddo
          k = ubound(fdir, 1)
          print *, k, 'mean ', &
            & 'edn', meanval(gedn(k,:,:)), &
            & 'eup', meanval(geup(k,:,:))
        endif
      endif
    end associate
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    if(lverbose .and. lhave_box) then
      print *,''
      if(allocated(buildings%edir)) then
        do i=1, size(buildings%iface)
          print *, 'building_face', i, 'edir', buildings%edir(i), &
            & 'in/out', buildings%incoming(i), buildings%outgoing(i)
        enddo
      else
        do i=1, size(buildings%iface)
          print *, 'building_face', i, &
            & 'in/out', buildings%incoming(i), buildings%outgoing(i)
        enddo
      endif
    endif

    call destroy_pprts(solver, .False.)
  end subroutine

end module
