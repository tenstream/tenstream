module m_plex2rayli

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i8, zero, one, default_str_len

  use m_helper_functions, only: CHKERR, deg2rad, itoa

  use m_netcdfio, only: ncwrite

  use m_plex_grid, only: t_plexgrid, &
    TOAFACE

  use m_plex_rt_base, only : t_state_container_plexrt
  use m_f2c_rayli, only: rfft_wedgeF90, rpt_img_wedgeF90

  implicit none

  logical, parameter :: ldebug=.False.

  contains

  !> @brief wrapper for the rayli MonteCarlo solver
  !> @details solve the radiative transfer equation with RayLi, currently only works for single task mpi runs
  subroutine rayli_wrapper(plex, kabs, ksca, g, albedo, sundir, solution, plck)
    use iso_c_binding
    type(t_plexgrid), intent(inout) :: plex
    type(tVec), intent(in) :: kabs, ksca, g, albedo
    type(tVec), intent(in), optional :: plck
    real(ireals), intent(in) :: sundir(:)
    type(t_state_container_plexrt), intent(inout) :: solution

    real(ireals), pointer :: xksca(:), xkabs(:), xg(:)

    logical :: lthermal, lsolar

    integer(mpiint) :: comm, myid, numnodes, ierr

    integer(iintegers) :: fStart, fEnd, cStart, cEnd, vStart, vEnd
    integer(iintegers) :: icell, iface, ivert, voff, idof
    integer(iintegers), pointer :: trans_closure(:), faces_of_cell(:)

    type(tVec) :: coordinates
    real(ireals), pointer :: coords(:)
    type(tPetscSection) :: coord_section

    integer(c_size_t), allocatable :: verts_of_face(:,:)
    integer(c_size_t), allocatable :: faces_of_wedges(:,:)
    real(c_double),    allocatable :: vert_coords(:,:)
    real(c_float),    allocatable :: rkabs(:), rksca(:), rg(:)
    real(c_float),    allocatable :: ralbedo_on_faces(:)
    real(c_float)                 :: rsundir(3)
    real(c_float),    allocatable :: flx_through_faces_edir(:)
    real(c_float),    allocatable :: flx_through_faces_ediff(:)
    real(c_float),    allocatable :: abso_in_cells(:)

    real(ireals) :: diffuse_point_origin(3)
    integer(c_size_t) :: Nphotons, Nwedges, Nfaces, Nverts
    real(ireals) :: opt_photons

    integer(c_size_t) :: outer_id
    logical :: lcyclic, lflg
    integer(c_int) :: icyclic

    opt_photons = 100000
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      "-rayli_photons", opt_photons, lflg,ierr) ; call CHKERR(ierr)
    Nphotons = int(opt_photons, c_size_t)

    call PetscObjectGetComm(plex%dm, comm, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(numnodes.gt.1) call CHKERR(numnodes, "Rayli currently only in single rank computations available")

    lsolar = solution%lsolar_rad
    lthermal = .not.solution%lsolar_rad

    if(present(plck).or.lthermal) &
      call CHKERR(1_mpiint, "You provided planck stuff, I guess you want to use thermal radiation computations."// &
                            "However, Rayli currently only supports solar radiation.")

    diffuse_point_origin = 0; idof=3
    call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      '-rayli_diff_flx_origin',diffuse_point_origin, idof, lflg, ierr); call CHKERR(ierr)
    if(lflg) call CHKERR(int(idof-i3, mpiint), 'must provide exactly 3 values for rayli_diff_flx_origin')

    lcyclic=.False.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, &
      '-rayli_cylic_bc', lcyclic, lflg, ierr); call CHKERR(ierr)
    if(lcyclic) then
      icyclic = 1
    else
      icyclic = 0
    endif

    !DEBUG if(solution%lsolar_rad) then
    !DEBUG   call VecSet(solution%edir, zero, ierr); call CHKERR(ierr)
    !DEBUG endif
    !DEBUG call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

    if(solution%lsolar_rad .and. norm2(sundir).le.zero) then
      return
    endif

    call DMPlexGetDepthStratum(plex%dm, i3, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetDepthStratum(plex%dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
    call DMPlexGetDepthStratum(plex%dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices
    Nwedges = cEnd - cStart
    Nfaces = fEnd - fStart
    Nverts = vEnd - vStart
    !print *,'Rayli Nwedges', Nwedges, 'Nfaces', Nfaces, 'Nverts', Nverts

    outer_id = -1

    allocate(verts_of_face(4, Nfaces), &
             faces_of_wedges(5, Nwedges), &
             vert_coords(3, Nverts), &
             rkabs(Nwedges), &
             rksca(Nwedges), &
             rg   (Nwedges), &
             ralbedo_on_faces(Nfaces), &
             flx_through_faces_edir(Nfaces), &
             flx_through_faces_ediff(2*Nfaces), &
             abso_in_cells(Nwedges) )

    do icell = cStart, cEnd-1
      call DMPlexGetCone(plex%dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
      faces_of_wedges(:,icell-cStart+1) = faces_of_cell - fStart
      call DMPlexRestoreCone(plex%dm, icell, faces_of_cell, ierr); call CHKERR(ierr)
    enddo

    do iface = fStart, fEnd-1
      call DMPlexGetTransitiveClosure(plex%dm, iface, PETSC_TRUE, trans_closure, ierr); call CHKERR(ierr)
      select case (size(trans_closure))
      case(14) ! 3 edges
        verts_of_face(1:3, iface-fStart+1) = trans_closure(9:size(trans_closure):2) - vStart
        verts_of_face(4  , iface-fStart+1) = outer_id
      case(18) ! 4 edges
        verts_of_face(:, iface-fStart+1) = trans_closure(11:size(trans_closure):2) - vStart
      end select
      call DMPlexRestoreTransitiveClosure(plex%dm, iface, PETSC_TRUE, trans_closure, ierr); call CHKERR(ierr)
    enddo

    call DMGetCoordinateSection(plex%dm, coord_section, ierr); call CHKERR(ierr)
    call DMGetCoordinatesLocal(plex%dm, coordinates, ierr); call CHKERR(ierr)
    call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

    do ivert = vStart, vEnd-1
      call PetscSectionGetOffset(coord_section, ivert, voff, ierr); call CHKERR(ierr)
      vert_coords(:, ivert-vStart+1) = coords(voff+1:voff+3)
    enddo
    call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    rkabs(:) = real(xkabs, kind(rkabs))
    rksca(:) = real(xksca, kind(rksca))
    rg   (:) = real(xg,    kind(rg   ))
    !call delta_scale(rkabs, rksca, rg, max_g=0._c_double)

    call VecRestoreArrayReadF90(kabs, xkabs, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(ksca, xksca, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(g   , xg   , ierr); call CHKERR(ierr)

    call fill_albedo(plex, albedo, ralbedo_on_faces)

    rsundir = real(-sundir, kind(rsundir))

!    call take_snap(comm, &
!      Nphotons, Nwedges, Nfaces, Nverts, &
!      verts_of_face, faces_of_wedges, vert_coords, &
!      rkabs, rksca, rg, &
!      ralbedo_on_faces, &
!      rsundir, &
!      solution, ierr)
    if(ierr.eq.1) return

    ierr = rfft_wedgeF90(Nphotons, Nwedges, Nfaces, Nverts, icyclic, &
      verts_of_face, faces_of_wedges, vert_coords, &
      rkabs, rksca, rg, &
      ralbedo_on_faces, rsundir, real(diffuse_point_origin, c_float), &
      flx_through_faces_edir, flx_through_faces_ediff, abso_in_cells ); call CHKERR(ierr)

!    call rayli_get_result(plex, &
!          flx_through_faces_edir, &
!          flx_through_faces_ediff, &
!          abso_in_cells, &
!          solution)

    !Twostream solver returns fluxes as [W/m^2]
    solution%lWm2 = .True.
    ! and mark solution that it is up to date, otherwise abso would be overwritten
    solution%lchanged  = .False.

    end subroutine

!    subroutine rayli_get_result(plex, &
!        flx_through_faces_edir, &
!        flx_through_faces_ediff, &
!        abso_in_cells, &
!        solution)
!      type(t_plexgrid), intent(inout) :: plex
!      real(c_float), intent(in) :: flx_through_faces_edir(:)
!      real(c_float), intent(in) :: flx_through_faces_ediff(:)
!      real(c_float), intent(in) :: abso_in_cells(:)
!      type(t_state_container_plexrt), intent(inout) :: solution
!
!      type(tIS) :: toa_ids
!      integer(iintegers) :: i, k, ke1, Ncol, ridx, numDof, geom_offset
!      integer(iintegers) :: icell, iface, idof, voff
!      integer(iintegers) :: ofStart, ofEnd, cStart, cEnd
!      integer(iintegers), pointer :: xtoa_faces(:)
!      real(ireals), pointer :: xedir(:), xediff(:), xabso(:), geoms(:)
!      real(ireals) :: area
!      type(tPetscSection) :: edir_section, ediff_section, abso_section, geomSection
!      integer(mpiint) :: ierr
!
!      call DMGetSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
!      call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
!
!      ! --------------------------- Direct Radiation ------------------
!      call DMPlexGetDepthStratum(plex%dm, i2, ofStart, ofEnd, ierr); call CHKERR(ierr) ! vertices
!      if(solution%lsolar_rad) then
!        call DMGetSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
!        call VecGetArrayF90(solution%edir, xedir, ierr); call CHKERR(ierr)
!        do iface = ofStart, ofEnd-1
!          call PetscSectionGetFieldOffset(geomSection, iface, i2, geom_offset, ierr); call CHKERR(ierr)
!          area = geoms(i1+geom_offset)
!
!          call PetscSectionGetOffset(edir_section, iface, voff, ierr); call CHKERR(ierr)
!          call PetscSectionGetDof(edir_section, iface, numDof, ierr); call CHKERR(ierr)
!          do idof = 1, numDof
!            xedir(voff+idof) = real(abs(flx_through_faces_edir(iface-ofStart+1)), ireals)! / area
!          enddo
!        enddo
!        call VecRestoreArrayF90(solution%edir, xedir, ierr); call CHKERR(ierr)
!      endif
!
!      ! --------------------------- Diffuse Radiation -----------------
!      call DMGetSection(plex%ediff_dm, ediff_section, ierr); call CHKERR(ierr)
!      call VecGetArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
!
!      call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', &
!        TOAFACE, toa_ids, ierr); call CHKERR(ierr)
!      call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)
!
!      ke1 = plex%Nlay+1
!
!      call ISGetIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
!      ridx=1
!      do i = 1, size(xtoa_faces)
!        iface = xtoa_faces(i)
!        do k = 0, ke1-2
!          call PetscSectionGetFieldOffset(ediff_section, iface+k, i0, voff, ierr); call CHKERR(ierr)
!          call PetscSectionGetFieldOffset(geomSection, iface+k, i2, geom_offset, ierr); call CHKERR(ierr)
!          area = geoms(i1+geom_offset)
!
!          xediff(i1+voff) = real(abs( flx_through_faces_ediff(ridx  ) ), ireals)! / area
!          xediff(i2+voff) = real(abs( flx_through_faces_ediff(ridx+1) ), ireals)! / area
!
!          ridx = ridx+2
!        enddo
!
!        ! at the surface, the ordering of incoming/outgoing fluxes is reversed because of cellid_surface == -1
!        call PetscSectionGetFieldOffset(ediff_section, iface+ke1-1, i0, voff, ierr); call CHKERR(ierr)
!        call PetscSectionGetFieldOffset(geomSection, iface+ke1-1, i2, geom_offset, ierr); call CHKERR(ierr)
!        area = geoms(i1+geom_offset)
!
!        xediff(i2+voff) = real(abs( flx_through_faces_ediff(ridx  ) ), ireals)! / area
!        xediff(i1+voff) = real(abs( flx_through_faces_ediff(ridx+1) ), ireals)! / area
!        ridx = ridx+2
!      enddo
!
!      call ISRestoreIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
!      call VecRestoreArrayF90(solution%ediff, xediff, ierr); call CHKERR(ierr)
!
!      ! --------------------------- Absorption ------------------------
!      call DMGetSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)
!      call VecGetArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)
!      call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
!      do icell = cStart, cEnd-1
!        call PetscSectionGetFieldOffset(geomSection, icell, i2, geom_offset, ierr); call CHKERR(ierr) ! cell volume
!        call PetscSectionGetOffset(abso_section, icell, voff, ierr); call CHKERR(ierr)
!        xabso(i1+voff) = real(abso_in_cells(i1+icell-cStart), ireals) / geoms(i1+geom_offset)
!      enddo
!      call VecRestoreArrayF90(solution%abso, xabso, ierr); call CHKERR(ierr)
!      call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
!    end subroutine

    subroutine fill_albedo(plex, albedo, ralbedo_on_faces)
      type(t_plexgrid), intent(inout) :: plex
      type(tVec), intent(in) :: albedo
      real(c_float), intent(out) :: ralbedo_on_faces(:)

      type(tIS) :: toa_ids
      integer(iintegers), pointer :: xtoa_faces(:)
      real(ireals), pointer :: xalbedo(:)
      integer(iintegers) :: srfc_face, i, fStart, fEnd
      integer(mpiint) :: ierr

      call DMGetStratumIS(plex%geom_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
      call ISGetIndicesF90(toa_ids, xtoa_faces, ierr); call CHKERR(ierr)
      call DMPlexGetDepthStratum(plex%dm, i2, fStart, fEnd, ierr); call CHKERR(ierr) ! faces

      call VecGetArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)

      ralbedo_on_faces(:) = -one
      do i=1,size(xtoa_faces)
        srfc_face = xtoa_faces(i) + plex%Nlay
        ralbedo_on_faces(srfc_face-fStart+1) = real(xalbedo(i), kind(ralbedo_on_faces))
      enddo

      call VecRestoreArrayReadF90(albedo, xalbedo, ierr); call CHKERR(ierr)
    end subroutine

!    subroutine take_snap(comm, &
!        Nphotons, Nwedges, Nfaces, Nverts, &
!        verts_of_face, faces_of_wedges, vert_coords, &
!        rkabs, rksca, rg, &
!        ralbedo_on_faces, &
!        rsundir, &
!        solution, ierr)
!      integer(mpiint),   intent(in) :: comm
!      integer(c_size_t), intent(in) :: Nphotons, Nwedges, Nfaces, Nverts
!      integer(c_size_t), intent(in) :: verts_of_face(:,:)
!      integer(c_size_t), intent(in) :: faces_of_wedges(:,:)
!      real(c_double),    intent(in) :: vert_coords(:,:)
!      real(c_float),     intent(in) :: rkabs(:), rksca(:), rg(:)
!      real(c_float),     intent(in) :: ralbedo_on_faces(:)
!      real(c_float),     intent(in) :: rsundir(3)
!      type(t_state_container_plexrt), intent(in) :: solution
!      integer(mpiint), intent(out) :: ierr
!
!      character(len=default_str_len) :: snap_path, groups(2)
!      logical :: lflg
!      integer(c_size_t) :: Nx=400, Ny=300
!      real(c_float), allocatable :: img(:,:)
!      real(c_float) :: cam_loc(3), cam_viewing_dir(3), cam_up_vec(3)
!      real(c_float) :: fov_width, fov_height
!      real(ireals), dimension(3) :: visit_focus, visit_view_normal, visit_view_up
!      real(ireals) :: visit_view_angle, visit_image_zoom, visit_parallel_scale
!      integer(iintegers) :: narg
!
!      integer(mpiint) :: myid
!
!      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
!
!      call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-rayli_snapshot", &
!        snap_path, lflg,ierr) ; call CHKERR(ierr)
!      if(lflg) then
!        if(len_trim(snap_path).eq.0) snap_path = 'rayli_snaphots.nc'
!        if(myid.eq.0) print *,'Capturing scene to file: '//trim(snap_path), len_trim(snap_path)
!        Nx = 400
!        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-rayli_snap_Nx', narg, lflg, ierr); call CHKERR(ierr)
!        if(lflg) Nx = narg
!
!        Ny = 300
!        call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-rayli_snap_Ny', narg, lflg, ierr); call CHKERR(ierr)
!        if(lflg) Ny = narg
!
!        allocate(img(Nx, Ny))
!
!        visit_view_angle = 30._ireals
!        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-visit_view_angle', visit_view_angle, lflg, ierr); call CHKERR(ierr)
!
!        visit_image_zoom = 1._ireals
!        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-visit_image_zoom', visit_image_zoom, lflg, ierr); call CHKERR(ierr)
!
!        visit_parallel_scale = 1e5_ireals
!        call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-visit_parallel_scale', visit_parallel_scale, lflg, ierr); call CHKERR(ierr)
!
!        visit_focus = 0._ireals
!        narg=3
!        call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-visit_focus', visit_focus, narg, lflg, ierr); call CHKERR(ierr)
!        if(lflg) &
!          call CHKERR(int(narg-3, mpiint), 'wrong number of input, need to be given comma separated without spaces')
!
!        visit_view_normal = -rsundir
!        narg=3
!        call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-visit_view_normal', visit_view_normal, narg, lflg, ierr); call CHKERR(ierr)
!        if(lflg) &
!          call CHKERR(int(narg-3, mpiint), 'wrong number of input, need to be given comma separated without spaces')
!
!        visit_view_up = [0,0,1]
!        narg=3
!        call PetscOptionsGetRealArray(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER,&
!          '-visit_view_up', visit_view_up, narg, lflg, ierr); call CHKERR(ierr)
!        if(lflg) &
!          call CHKERR(int(narg-3, mpiint), 'wrong number of input, need to be given comma separated without spaces')
!
!        cam_viewing_dir = -real(visit_view_normal, kind(cam_viewing_dir))
!        cam_up_vec = real(visit_view_up, kind(cam_up_vec))
!        cam_loc = real(visit_focus + visit_view_normal * visit_parallel_scale &
!          / tan(deg2rad(visit_view_angle)/2), kind(cam_loc))
!
!        fov_width = 2 * real(tan(deg2rad(visit_view_angle)/2) / visit_image_zoom, kind(fov_width))
!        fov_height = fov_width * real(Ny, c_float) / real(Nx, c_float)
!
!        ierr = rpt_img_wedgeF90( Nx, Ny, &
!          Nphotons, Nwedges, Nfaces, Nverts, &
!          verts_of_face, faces_of_wedges, vert_coords, &
!          rkabs, rksca, rg, &
!          ralbedo_on_faces, &
!          rsundir, & ! DEBUG note the kabs/ksca/
!          cam_loc, cam_viewing_dir, cam_up_vec, &
!          fov_width, fov_height, &
!          img); call CHKERR(ierr)
!
!        groups(1) = trim(snap_path)
!        groups(2) = "rpt_img_"//itoa(solution%uid)
!        call ncwrite(groups, img, ierr); call CHKERR(ierr)
!        ierr = 1
!        return
!      endif
!    end subroutine
end module
