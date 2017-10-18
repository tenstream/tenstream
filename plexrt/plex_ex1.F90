module m_mpi_plex_ex1

#include "petsc/finclude/petsc.h"
  use petsc
  use m_helper_functions, only: CHKERR, norm, imp_bcast, determine_normal_direction, spherical_2_cartesian
  use m_data_parameters, only : ireals, iintegers, mpiint, &
    default_str_len, &
    i0, i1, i2, i3, i4, i5,  &
    zero, one,       &
    init_mpi_data_parameters, myid

  use m_plex_grid, only: t_plexgrid, load_plex_from_file, &
                       compute_face_geometry, print_dmplex,   &
                       setup_edir_dmplex, setup_abso_dmplex

  use m_plex_rt, only: create_edir_mat, create_src_vec, compute_edir_absorption

  implicit none

  logical, parameter :: ldebug=.True.

  integer(mpiint) :: ierr

  contains


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
      !call KSPSolve(ksp, b, x, ierr); CHKERRQ(ierr)
      call KSPDestroy(ksp, ierr); CHKERRQ(ierr)
    end subroutine

    subroutine plex_ex1(plex)
      type(t_plexgrid) :: plex
      type(tVec) :: b, abso, edir
      type(tMat) :: A

      real(ireals) :: sundir(3) ! cartesian direction of sun rays in a global reference system


      call compute_face_geometry(plex)


      call setup_edir_dmplex(plex, plex%edir_dm)
      sundir = get_normal_of_first_TOA_face(plex)
      print *,'get_normal_of_first_TOA_face', sundir
      !sundir = [-0.71184089224108049, -3.7794622710146053E-002, -0.70132311428301020] ! original zenith 0
      !sundir = [-0.71184089224108049, -3.7794622710146053E-002, -0.60132311428301020] ! zenith 10 azi 1
      !sundir = [-0.51184089224108049, +0.37794622710146053, -0.00132311428301020] ! zenith 63 azi 17

      !sundir = [-one/100,-one/10,-one]
      !sundir = [one,one,-one]
      !sundir = spherical_2_cartesian(180*one,60*one)
      !sundir = spherical_2_cartesian(90*one,60*one)

      sundir = sundir / norm(sundir)
      print *,'sundir',sundir
      !stop 'debug'

      call print_dmplex(plex%comm, plex%edir_dm)

      call create_src_vec(plex%edir_dm, b)

      call setup_abso_dmplex(plex, plex%abso_dm)

      call create_edir_mat(plex, sundir, A)

      call solve(plex, b, A, edir)

      call PetscObjectViewFromOptions(edir, PETSC_NULL_VEC, '-show_edir', ierr); call CHKERR(ierr)

      call compute_edir_absorption(plex, edir, sundir, abso)

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
