module m_plex_rt

#include "petsc/finclude/petsc.h"
  use petsc
  use m_helper_functions, only: CHKERR, determine_normal_direction, &
    angle_between_two_vec, rad2deg, deg2rad, &
    vec_proj_on_plane, cross_3d, norm, rotation_matrix_world_to_local_basis, &
    imp_bcast, approx
  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5,  &
    zero, one, pi

  use m_plex_grid, only: t_plexgrid, TOP_BOT_FACE, SIDE_FACE

  implicit none

  private

  public :: create_src_vec, solve_plex_rt, compute_edir_absorption, create_edir_mat, get_normal_of_first_TOA_face

  logical, parameter :: ldebug=.True.
  contains

    subroutine create_src_vec(dm, globalVec)
      type(tDM),allocatable :: dm
      type(tVec) :: globalVec
      integer(iintegers) :: i, voff
      type(tPetscSection) :: s

      type(tVec) :: localVec
      real(ireals), pointer :: xv(:)

      integer(iintegers) :: cStart, cEnd
      integer(iintegers) :: fStart, fEnd

      type(tIS) :: toa_ids

      integer(iintegers), pointer :: xx_v(:)

      integer(mpiint) :: ierr

      if(.not.allocated(dm)) stop 'called create_src_vec but face_dm is not allocated'

      call DMGetDefaultSection(dm, s, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(s, PETSC_NULL_SECTION, '-show_src_section', ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(dm, 0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      call DMPlexGetHeightStratum(dm, 1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges

      ! Now lets get vectors!
      call DMGetGlobalVector(dm, globalVec,ierr); call CHKERR(ierr)
      call PetscObjectSetName(globalVec, 'srcVecGlobal', ierr);call CHKERR(ierr)
      call VecSet(globalVec, zero, ierr); call CHKERR(ierr)

      call DMGetLocalVector(dm, localVec,ierr); call CHKERR(ierr)
      call PetscObjectSetName(localVec, 'srcVec', ierr);call CHKERR(ierr)
      call VecSet(localVec, zero, ierr); call CHKERR(ierr)

      call DMGetStratumIS(dm, 'TOA', i1, toa_ids, ierr); call CHKERR(ierr)
      if (toa_ids.eq.PETSC_NULL_IS) then ! dont have TOA points
      else
        call PetscObjectViewFromOptions(toa_ids, PETSC_NULL_IS, '-show_IS_TOA', ierr); call CHKERR(ierr)

        call ISGetIndicesF90(toa_ids, xx_v, ierr); call CHKERR(ierr)

        call VecGetArrayF90(localVec, xv, ierr); call CHKERR(ierr)

        !do i = 1, size(xx_v)
        do i = size(xx_v), size(xx_v)
          print *,'debug: setting only last source entry...'
          call PetscSectionGetOffset(s, xx_v(i), voff, ierr); call CHKERR(ierr)
          !print *,myid,'index:',i,xx_v(i),'off',voff+1, lbound(xv,1), ubound(xv,1)
          xv(voff+1) = 100 !i
        enddo
        call VecRestoreArrayF90(localVec, xv, ierr); call CHKERR(ierr)

        call ISRestoreIndicesF90(toa_ids, xx_v, ierr); call CHKERR(ierr)
      endif

      call PetscObjectViewFromOptions(localVec, PETSC_NULL_VEC, '-show_src_vec_local', ierr); call CHKERR(ierr)

      call DMLocalToGlobalBegin(dm, localVec, ADD_VALUES, globalVec, ierr); call CHKERR(ierr)
      call DMLocalToGlobalEnd(dm, localVec, ADD_VALUES, globalVec, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(globalVec, PETSC_NULL_VEC, '-show_src_vec_global', ierr); call CHKERR(ierr)

      call DMRestoreLocalVector(dm, localVec, ierr); call CHKERR(ierr)

      !call DMRestoreGlobalVector(dm, globalVec, ierr); call CHKERR(ierr) ! Dont destroy the output vec
    end subroutine

    subroutine solve_plex_rt(plex, b, A, x)
      type(t_plexgrid) :: plex
      type(tVec), intent(in) :: b
      type(tVec), intent(inout) :: x
      type(tMat), intent(in) :: A

      type(tKSP) :: ksp
      integer(mpiint) :: ierr

      call KSPCreate(plex%comm, ksp, ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp, A, A, ierr); CHKERRQ(ierr)

      call KSPSetFromOptions(ksp, ierr); CHKERRQ(ierr)
      call KSPSetUp(ksp, ierr); CHKERRQ(ierr)

      call VecDuplicate(b, x, ierr); CHKERRQ(ierr)
      call KSPSolve(ksp, b, x, ierr); CHKERRQ(ierr)
      call KSPDestroy(ksp, ierr); CHKERRQ(ierr)
    end subroutine


  subroutine compute_edir_absorption(plex, edir, sundir, abso)
    type(t_plexgrid), intent(inout) :: plex
    type(tVec), intent(in) :: edir
    real(ireals), intent(in) :: sundir(3)
    type(tVec), intent(out) :: abso

    type(tVec) :: local_edir

    real(ireals), pointer :: xedir(:), xabso(:)

    type(tPetscSection) :: abso_section, edir_section

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: icell, iface
    integer(iintegers),pointer :: faces_of_cell(:)
    integer(mpiint) :: ierr

    type(tPetscSection) :: geomSection
    real(ireals) :: cell_center(3), face_normal(3), face_center(3)
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
    integer(iintegers) :: geom_offset, abso_offset, edir_offset
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    real(ireals) :: angle_to_sun

    if(.not.allocated(plex%edir_dm) .or. .not.allocated(plex%abso_dm)) stop 'called compute_edir_absorption with a dm which is not allocated?'

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%edir_dm, edir_section, ierr); call CHKERR(ierr)
    call DMGetDefaultSection(plex%abso_dm, abso_section, ierr); call CHKERR(ierr)

    call DMPlexGetHeightStratum(plex%abso_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells

    ! Now lets get vectors!
    call DMGetGlobalVector(plex%abso_dm, abso,ierr); call CHKERR(ierr)
    call PetscObjectSetName(abso, 'absoVecGlobal', ierr);call CHKERR(ierr)
    call VecSet(abso, zero, ierr); call CHKERR(ierr)

    call DMGetLocalVector(plex%edir_dm, local_edir,ierr); call CHKERR(ierr)
    call VecSet(local_edir, zero, ierr); call CHKERR(ierr)
    call DMGlobalToLocalBegin(plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)
    call DMGlobalToLocalEnd  (plex%edir_dm, edir, INSERT_VALUES, local_edir, ierr); call CHKERR(ierr)

    call VecGetArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)
    call VecGetArrayF90(abso, xabso, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1
      call PetscSectionGetOffset(abso_section, icell, abso_offset, ierr); call CHKERR(ierr)

      call DMPlexGetCone(plex%abso_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(edir_section, faces_of_cell(iface), edir_offset, ierr); call CHKERR(ierr)

        call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
        cell_center = geoms(1+geom_offset:3+geom_offset)

        call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); call CHKERR(ierr)
        face_center = geoms(1+geom_offset: 3+geom_offset)
        face_normal = geoms(4+geom_offset: 6+geom_offset)

        ! Determine the inward normal vec for the face
        face_normal = face_normal * determine_normal_direction(face_normal, face_center, cell_center)

        ! Then determine if the face is src(in line with the sun vec) or if it is destination(contra sun direction)
        angle_to_sun = angle_between_two_vec(face_normal, sundir)

        if(angle_to_sun.lt.pi/2) then
          lsrc(iface) = .True.
          xabso(abso_offset+i1) = xabso(abso_offset+i1) + xedir(edir_offset+i1)
        else
          lsrc(iface) = .False.
          xabso(abso_offset+i1) = xabso(abso_offset+i1) - xedir(edir_offset+i1)
        endif
      enddo
      call DMPlexRestoreCone(plex%abso_dm, icell, faces_of_cell,ierr); call CHKERR(ierr)
    enddo

    call VecRestoreArrayF90(abso, xabso, ierr); call CHKERR(ierr)
    call VecRestoreArrayReadF90(local_edir, xedir, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(abso, PETSC_NULL_VEC, '-show_abso', ierr); call CHKERR(ierr)
  end subroutine

  subroutine create_edir_mat(plex, sundir, A)
    type(t_plexgrid), intent(inout) :: plex
    real(ireals), intent(in) :: sundir(3)
    type(tMat), intent(out) :: A

    type(tPetscSection) :: sec

    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(mpiint) :: ierr

    integer(iintegers), pointer :: faces_of_cell(:)
    integer(iintegers) :: iface, irow, icol, icell, idst

    integer(iintegers) :: irows(5)
    real(ireals) :: dir2dir_coeffs(5)

    type(tDMLabel) :: faceposlabel, zindexlabel

    type(tPetscSection) :: geomSection
    real(ireals), pointer :: geoms(:) ! pointer to coordinates vec

    real(ireals) :: zenith, azimuth

    integer(iintegers) :: base_face   ! index of face which is closest to sun angle, index regarding faces_of_cell
    integer(iintegers) :: left_face   ! index of face which is left/right of base face, index regarding faces_of_cell
    integer(iintegers) :: right_face  ! index of face which is left/right of base face, index regarding faces_of_cell
    integer(iintegers) :: upper_face  ! index of face which is top/bot of base face, index regarding faces_of_cell
    integer(iintegers) :: bottom_face ! index of face which is top/bot of base face, index regarding faces_of_cell

    real(ireals) :: dir2dir(5,5), S(5), T(5) ! Nface**2
    logical :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); call CHKERR(ierr)
    call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

    call DMGetDefaultSection(plex%edir_dm, sec, ierr); call CHKERR(ierr)
    call DMPlexGetHeightStratum(plex%edir_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
    call DMPlexGetHeightStratum(plex%edir_dm, i1, fStart, fEnd, ierr); call CHKERR(ierr) ! faces / edges
    call DMGetLabel(plex%edir_dm, "Face Position", faceposlabel, ierr); call CHKERR(ierr)
    call DMGetLabel(plex%edir_dm, "Vertical Index", zindexlabel, ierr); call CHKERR(ierr)

    call DMCreateMatrix(plex%edir_dm, A, ierr); call CHKERR(ierr)
    !call MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, ierr); call CHKERR(ierr)

    do icell = cStart, cEnd-1

      call DMPlexGetCone(plex%edir_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell

      call compute_local_wedge_ordering(icell, faces_of_cell, geomSection, geoms, faceposlabel, zindexlabel, sundir, &
        zenith, azimuth, upper_face, bottom_face, base_face, left_face, right_face, lsrc)

      do iface = 1, size(faces_of_cell) !DEBUG

        if(lsrc(iface)) then ! this is really a source

          !retrieve coeffs:
          if (iface.eq.upper_face) then
            call compute_dir2dir_coeff(i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.bottom_face) then
            call compute_dir2dir_coeff(5*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.base_face) then
            call compute_dir2dir_coeff(2*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.left_face) then
            call compute_dir2dir_coeff(3*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else if (iface.eq.right_face) then
            call compute_dir2dir_coeff(4*i1, rad2deg(azimuth), rad2deg(zenith), S, T)
          else
            stop 'iface is not in local wedge face numbering... something must have gone terribly wrong in compute_local_wedge_ordering'
          endif

          dir2dir([upper_face,base_face,left_face,right_face,bottom_face], iface) = T

          call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces

          do idst=1,size(faces_of_cell)
            call PetscSectionGetOffset(sec, faces_of_cell(idst), irow, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces
            irows(idst) = irow
            dir2dir_coeffs(idst) = -dir2dir(idst, iface)
          enddo

          call MatSetValuesLocal(A, size(faces_of_cell), irows, i1, [icol], dir2dir_coeffs, INSERT_VALUES, ierr); call CHKERR(ierr)
        endif
      enddo


      ! Set diagonal entries
      do iface = 1, size(faces_of_cell)
        call PetscSectionGetOffset(sec, faces_of_cell(iface), irow, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces
        call PetscSectionGetOffset(sec, faces_of_cell(iface), icol, ierr); call CHKERR(ierr) ! this is the offset of the neighboring faces
        call MatSetValuesLocal(A, i1, [irow], i1, [icol], [one] , INSERT_VALUES, ierr); call CHKERR(ierr)
        !print *,'Setting diagonal entry:', icell, iface, '::', faces_of_cell(iface), '->', irow, icol
      enddo

      call DMPlexRestoreCone(plex%edir_dm, icell, faces_of_cell,ierr); call CHKERR(ierr)
    enddo

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call PetscObjectViewFromOptions(A, PETSC_NULL_VEC, '-show_A', ierr); call CHKERR(ierr)

    call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

  end subroutine

  subroutine compute_local_wedge_ordering(icell, faces_of_cell, geomSection, geoms, faceposlabel, zindexlabel, sundir, &
      zenith, azimuth, upper_face, bottom_face, base_face, left_face, right_face, lsrc)
    integer(iintegers), intent(in) :: icell
    integer(iintegers), intent(in), pointer :: faces_of_cell(:)
    type(tPetscSection) :: geomSection
    real(ireals), intent(in), pointer :: geoms(:) ! pointer to coordinates vec
    type(tDMLabel), intent(in) :: faceposlabel, zindexlabel
    real(ireals), intent(in) :: sundir(3)
    real(ireals), intent(out) :: zenith, azimuth
    integer(iintegers), intent(out) :: upper_face, bottom_face, base_face, left_face, right_face
    logical, intent(out) :: lsrc(5) ! is src or destination of solar beam (5 faces in a wedge)

    integer(iintegers) :: iface, geom_offset, ifacepos, izindex(2)

    real(ireals) :: cell_center(3)
    real(ireals) :: face_normals(3,5), face_centers(3,5)
    real(ireals) :: side_faces_angles_to_sun(5), proj_angles_to_sun(3), proj_normal(3)
    real(ireals) :: e_x(3), e_y(3), e_z(3) ! unit vectors of local coord system in which we compute the transfer coefficients
    real(ireals) :: projected_sundir(3)

    integer(iintegers) :: side_faces(3), top_faces(2) ! indices in faces_of_cell which give the top/bot and side faces via labeling
    integer(iintegers) :: iside_faces, itop_faces, ibase_face ! indices to fill above arrays
    real(ireals) :: MrotWorld2Local(3,3) ! Rotation Matrix into the local wedge space, (ex, ey, ez)

    real(ireals) :: local_normal_left(3), local_normal_right(3) ! normal vectors in local wedgemc geometry. For even sided triangles, this is smth. like left: [.5, -.8] or right [-.5, -.8]

    integer(mpiint) :: ierr

    call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); call CHKERR(ierr)
    cell_center = geoms(1+geom_offset:3+geom_offset)

    do iface = 1, size(faces_of_cell)
      call PetscSectionGetOffset(geomSection, faces_of_cell(iface), geom_offset, ierr); call CHKERR(ierr)
      face_centers(:,iface) = geoms(1+geom_offset: 3+geom_offset)
      face_normals(:,iface) = geoms(4+geom_offset: 6+geom_offset)
      face_normals(:,iface) = face_normals(:,iface) * determine_normal_direction(face_normals(:,iface), face_centers(:, iface), cell_center)
      ! to find the face whose normal is closest to sun vector:
      side_faces_angles_to_sun(iface) = angle_between_two_vec(face_normals(:,iface), sundir) ! in [rad]
      lsrc(iface) = side_faces_angles_to_sun(iface).lt.pi/2-deg2rad(.1_ireals) ! dont propagate energy along edges where sun is only .1 degrees off
    enddo

    ! get numbering for 2 top/bot faces and 3 side faces which indicate the position in the faces_of_cell vec
    iside_faces = 1; itop_faces = 1
    do iface = 1, size(faces_of_cell)
      call DMLabelGetValue(faceposlabel, faces_of_cell(iface), ifacepos, ierr); call CHKERR(ierr)
      select case(ifacepos)
      case(TOP_BOT_FACE)
        top_faces(itop_faces) = iface
        itop_faces = itop_faces+1
      case(SIDE_FACE)
        side_faces(iside_faces) = iface
        iside_faces = iside_faces+1
      case default
        stop 'wrong or no side face label'
      end select
    enddo

    call DMLabelGetValue(zindexlabel, faces_of_cell(top_faces(1)), izindex(1), ierr); call CHKERR(ierr)
    call DMLabelGetValue(zindexlabel, faces_of_cell(top_faces(2)), izindex(2), ierr); call CHKERR(ierr)
    if(izindex(1).gt.izindex(2)) then
      upper_face = top_faces(2)
      bottom_face = top_faces(1)
    else
      upper_face = top_faces(1)
      bottom_face = top_faces(2)
    endif

    do iface=1,size(side_faces)
      if(all(approx(sundir, face_normals(:,upper_face)))) then
        proj_angles_to_sun(iface) = 0
      else
        proj_normal = vec_proj_on_plane(sundir, face_normals(:,upper_face))
        !print *,'vec_proj_on_plane',iface, '::', sundir, 'vs', face_normals(:,upper_face), '->', proj_normal
        proj_angles_to_sun(iface) = angle_between_two_vec(proj_normal, face_normals(:,side_faces(iface)))
      endif
    enddo

    ibase_face = minloc(proj_angles_to_sun,dim=1)
    base_face = side_faces(ibase_face)

    e_y = face_normals(:, base_face)   ! inward facing normal -> in local wedgemc coordinates
    e_z = -face_normals(:, upper_face) ! outward facing normal with respect to the top plate
    e_x = cross_3d(e_y, e_z)           ! in local wedge_coords, this is y=0 coordinate

    MrotWorld2Local = rotation_matrix_world_to_local_basis(e_x, e_y, e_z)

    left_face  = side_faces(modulo(ibase_face,size(side_faces))+1)
    right_face = side_faces(modulo(ibase_face+1,size(side_faces))+1)

    local_normal_left  = matmul(MrotWorld2Local, face_normals(:,left_face))
    local_normal_right = matmul(MrotWorld2Local, face_normals(:,right_face))

    if(local_normal_left(1).lt.local_normal_right(1)) then ! switch right and left face
      iface = right_face
      right_face = left_face
      left_face = iface
    endif

    ! Now we all the info for the local wedge calculations as we do em with the MonteCarlo raytracer

    zenith = angle_between_two_vec(sundir, -e_z)

    projected_sundir = vec_proj_on_plane(matmul(MrotWorld2Local, sundir), [zero,zero,one])
    if(norm(projected_sundir).le.epsilon(zero)) then
      azimuth = 0
    else
      projected_sundir = projected_sundir
      azimuth = angle_between_two_vec([zero,one,zero], projected_sundir) * sign(one, projected_sundir(1))
    endif

    if(ldebug .and. norm(face_centers(:,upper_face)) .le. norm(face_centers(:,bottom_face))) then ! we expect the first face to be the upper one
      print *,'norm upper_face ', norm(face_centers(:,upper_face))
      print *,'norm bottom_face', norm(face_centers(:,bottom_face))
      print *,'we expect the first face to be the upper one but found:',icell, faces_of_cell(1), faces_of_cell(2)
      stop 'create_edir_mat() :: wrong zindexlabel'
    endif

    if(azimuth.lt.-60 .or. azimuth.gt.60) stop 'local azimuth greater than 60 deg. something must have gone wrong with the base face selection!'
    end subroutine

    subroutine compute_dir2dir_coeff(src, phi, theta, S,T)
      use m_boxmc, only : t_boxmc,t_boxmc_wedge_5_5
      integer(iintegers), intent(in) :: src
      real(ireals), intent(in) :: phi, theta
      real(ireals), intent(out) :: S(5),T(5)

      type(t_boxmc_wedge_5_5) :: bmc_wedge_5_5
      real(ireals) :: bg(3), dx,dy,dz
      real(ireals) :: S_tol(5),T_tol(5)

      call bmc_wedge_5_5%init(PETSC_COMM_SELF)
      !print *,'computing coeffs for src/phi/theta',src,phi,theta

      bg  = [1e-3_ireals, zero, one/2 ]

      !phi   =  0
      !theta = 45

      dx = 100
      dy = dx
      dz = 200

      call bmc_wedge_5_5%get_coeff(PETSC_COMM_SELF, bg, src, .True., &
        phi, theta, dx, dy, dz, S, T, S_tol, T_tol, inp_atol=5e-2_ireals, inp_rtol=2e-1_ireals)
    end subroutine

    function get_normal_of_first_TOA_face(plex)
      type(t_plexgrid) :: plex
      real(ireals) :: get_normal_of_first_TOA_face(3)
      type(tPetscSection) :: geomSection
      real(ireals), pointer :: geoms(:) ! pointer to coordinates vec
      real(ireals) :: cell_center(3), face_center(3)
      real(ireals),allocatable :: face_normal(:)

      type(tIS) :: toa_ids
      integer(iintegers), pointer :: xitoa(:), cell_support(:)
      integer(iintegers) :: geom_offset, iface, icell

      integer(mpiint) :: myid, ierr

      call mpi_comm_rank(plex%comm, myid, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
        !call VecGetArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)
        !call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)
        call DMGetStratumIS(plex%geom_dm, 'TOA', i1, toa_ids, ierr); call CHKERR(ierr)

        if (toa_ids.eq.PETSC_NULL_IS) then ! dont have TOA points
          stop 'This didnt work, we tried to set the sundir according to first face on rank 0 but it seems he does not have TOA faces'
        else
          call ISGetIndicesF90(toa_ids, xitoa, ierr); call CHKERR(ierr)
          iface = xitoa(1) ! first face of TOA faces

          call VecGetArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)
          call PetscSectionGetOffset(geomSection, iface, geom_offset, ierr); call CHKERR(ierr)

          face_center = geoms(geom_offset+1:geom_offset+3)
          allocate(face_normal(3))
          face_normal = geoms(geom_offset+4:geom_offset+6)

          call DMPlexGetSupport(plex%geom_dm, iface, cell_support, ierr); CHKERRQ(ierr) ! support of face is cell
          icell = cell_support(1)
          call DMPlexRestoreSupport(plex%geom_dm, iface, cell_support, ierr); CHKERRQ(ierr) ! support of face is cell

          call PetscSectionGetOffset(geomSection, icell, geom_offset, ierr); CHKERRQ(ierr)
          cell_center = geoms(1+geom_offset:3+geom_offset)

          ! Determine the inward normal vec for the face
          face_normal = face_normal * determine_normal_direction(face_normal, face_center, cell_center)

          call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); call CHKERR(ierr)

          call ISRestoreIndicesF90(toa_ids, xitoa, ierr); call CHKERR(ierr)
        endif
      endif

      call imp_bcast(plex%comm, face_normal, 0_mpiint)
      get_normal_of_first_TOA_face = face_normal

    end function
end module
