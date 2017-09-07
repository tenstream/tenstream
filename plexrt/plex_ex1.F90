module m_mpi_plex_ex1

#include "petsc/finclude/petsc.h"
  use petsc
  use m_helper_functions, only: CHKERR, norm, imp_bcast, determine_normal_direction
  use m_data_parameters, only : ireals, iintegers, mpiint, &
    default_str_len, &
    i0, i1, i2, i3, i4, i5,  &
    zero, one,       &
    init_mpi_data_parameters, myid

  use m_icon_plexgrid, only: t_plexgrid, load_plex_from_file, &
                       compute_face_geometry, print_dmplex,   &
                       setup_edir_dmplex, setup_abso_dmplex,  &
                       compute_edir_absorption, create_edir_mat

  implicit none

  logical, parameter :: ldebug=.True.

  integer(mpiint) :: ierr

  contains

    subroutine create_src_vec(dm, globalVec)
      type(tDM),allocatable :: dm
      type(tVec) :: globalVec
      integer(iintegers) :: i, voff
      type(tPetscSection) :: s

      Vec :: localVec
      real(ireals), pointer :: xv(:)

      PetscInt :: cStart, cEnd
      PetscInt :: fStart, fEnd

      type(tIS) :: toa_ids

      PetscInt, pointer :: xx_v(:)

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


      if(myid.eq.0) then
        call DMGetDefaultSection(plex%geom_dm, geomSection, ierr); CHKERRQ(ierr)
        !call VecGetArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)
        !call VecRestoreArrayReadF90(plex%geomVec, geoms, ierr); CHKERRQ(ierr)
        call DMGetStratumIS(plex%geom_dm, 'TOA', i1, toa_ids, ierr); call CHKERR(ierr)

        if (toa_ids.eq.PETSC_NULL_IS) then ! dont have TOA points
          stop 'This didnt work, we tried to set the sunpos according to first face on rank 0 but it seems he does not have TOA faces'
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

    subroutine solve(plex, b, A, x)
      type(t_plexgrid) :: plex
      Vec, intent(in) :: b
      Vec, intent(inout) :: x
      Mat, intent(in) :: A

      KSP :: ksp

      call KSPCreate(plex%comm, ksp, ierr); CHKERRQ(ierr)
      call KSPSetOperators(ksp, A, A, ierr); CHKERRQ(ierr)

      call KSPSetFromOptions(ksp, ierr); CHKERRQ(ierr)
      call KSPSetUp(ksp, ierr); CHKERRQ(ierr)

      call VecDuplicate(b, x, ierr); CHKERRQ(ierr)
      call KSPSolve(ksp, b, x, ierr); CHKERRQ(ierr)
      call KSPDestroy(ksp, ierr); CHKERRQ(ierr)
    end subroutine

    subroutine plex_ex1(plex)
      type(t_plexgrid) :: plex
      type(tVec) :: b, abso, edir
      type(tMat) :: A

      call compute_face_geometry(plex)


      call setup_edir_dmplex(plex, plex%edir_dm)
      plex%sunpos = get_normal_of_first_TOA_face(plex)
      print *,plex%sunpos
      plex%sunpos = [0.65403449436133709,  0.13472543825908251,  2.74437083263075809]
      plex%sunpos = plex%sunpos / norm(plex%sunpos)
      print *,plex%sunpos

      call print_dmplex(plex%comm, plex%edir_dm)

      call create_src_vec(plex%edir_dm, b)

      call setup_abso_dmplex(plex, plex%abso_dm)

      call create_edir_mat(plex, A)

      call solve(plex, b, A, edir)

      call PetscObjectViewFromOptions(edir, PETSC_NULL_VEC, '-show_edir', ierr); call CHKERR(ierr)

      call compute_edir_absorption(plex, edir, abso)

    end subroutine
end module

program main
  use m_mpi_plex_ex1
  implicit none

  character(len=default_str_len) :: gridfile
  logical :: lflg

  type(t_plexgrid) :: plex


  call PetscInitialize(PETSC_NULL_CHARACTER,ierr); call CHKERR(ierr)
  call init_mpi_data_parameters(PETSC_COMM_WORLD)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-plex', gridfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) stop 'need to supply a plex filename... please call with -plex <fname_of_plexfile.h5>'

  call load_plex_from_file(PETSC_COMM_WORLD, gridfile, plex)

  call plex_ex1(plex)

  call PetscFinalize(ierr)
end program
