module m_gen_fish_plex
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, i0, i1, init_mpi_data_parameters

  use m_icon_plexgrid, only : TOP_BOT_FACE, SIDE_FACE

  implicit none


  PetscErrorCode :: ierr

  contains
    subroutine init()
      type(tDM) :: dmcell

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr); CHKERRQ(ierr)

      call init_mpi_data_parameters(PETSC_COMM_WORLD)

      call create_plex_serial(dmcell)

      call PetscObjectViewFromOptions(dmcell, PETSC_NULL_VEC, "-show_plex", ierr); CHKERRQ(ierr)

      call DMDestroy(dmcell, ierr);CHKERRQ(ierr)
      call PetscFinalize(ierr)
    end subroutine

    subroutine create_plex_serial(dm)
      DM :: dm

      call DMPlexCreate(PETSC_COMM_SELF, dm, ierr);CHKERRQ(ierr)
      call PetscObjectSetName(dm, 'testplex', ierr);CHKERRQ(ierr)
      call DMSetDimension(dm, 3, ierr);CHKERRQ(ierr)

      call DMPlexSetChart(dm, 0, 57, ierr); CHKERRQ(ierr)

      call setup_plex(dm)

    end subroutine
    subroutine setup_plex(dm)
      DM :: dm

      PetscInt :: pStart, pEnd
      PetscInt :: cStart, cEnd
      PetscInt :: fStart, fEnd
      PetscInt :: eStart, eEnd
      PetscInt :: vStart, vEnd

      call set_wedge_connectivity(dm)

      call DMPlexSymmetrize(dm, ierr); CHKERRQ(ierr)
      call DMPlexStratify(dm, ierr); CHKERRQ(ierr)

      call set_coords_serial(dm)

      call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
      call DMPlexGetHeightStratum(dm, 0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
      call DMPlexGetHeightStratum(dm, 1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
      call DMPlexGetDepthStratum (dm, 1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
      call DMPlexGetDepthStratum (dm, 0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices

      print *,'pStart,End :: ',pStart, pEnd
      print *,'cStart,End :: ',cStart, cEnd
      print *,'fStart,End :: ',fStart, fEnd
      print *,'eStart,End :: ',eStart, eEnd
      print *,'vStart,End :: ',vStart, vEnd

      call create_face_labels(dm)
    end subroutine

    subroutine set_wedge_connectivity(dm)
      DM :: dm
      integer :: i
      ! Preallocation
      ! Every cell has 5 faces
      do i=0,3
        call DMPlexSetConeSize(dm, i, 5, ierr); CHKERRQ(ierr)
      enddo

      ! Faces have 3 or 4 edges
      do i=4,7
        call DMPlexSetConeSize(dm, i, 3, ierr); CHKERRQ(ierr)
      enddo
      do i=17,20
        call DMPlexSetConeSize(dm, i, 3, ierr); CHKERRQ(ierr)
      enddo
      do i=8,16
        call DMPlexSetConeSize(dm, i, 4, ierr); CHKERRQ(ierr)
      enddo

      ! Edges have 2 vertices
      do i=21,44
        call DMPlexSetConeSize(dm, i, 2, ierr); CHKERRQ(ierr)
      enddo

      call DMSetUp(dm, ierr); CHKERRQ(ierr) ! Allocate space for cones

      ! Setup Connections
      call DMPlexSetCone(dm,  0, [ 4, 8, 9,10,17], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  1, [ 5,10,11,12,18], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  2, [ 6,12,13,14,19], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  3, [ 7,14,15,16,20], ierr); CHKERRQ(ierr)

      call DMPlexSetCone(dm,  4, [21,22,23], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  5, [23,24,25], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  6, [25,26,27], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  7, [27,28,29], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 17, [36,37,38], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 18, [38,39,40], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 19, [40,41,42], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 20, [42,43,44], ierr); CHKERRQ(ierr)

      call DMPlexSetCone(dm,  8, [21,30,31,36], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  9, [22,30,32,37], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 10, [23,31,32,38], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 11, [24,31,33,39], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 12, [25,32,33,40], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 13, [26,33,34,41], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 14, [27,32,34,42], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 15, [28,32,35,43], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 16, [29,34,35,44], ierr); CHKERRQ(ierr)

      call DMPlexSetCone(dm, 21, [45,46], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 22, [45,47], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 23, [46,47], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 24, [46,48], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 25, [47,48], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 26, [48,49], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 27, [47,49], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 28, [47,50], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 29, [49,50], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 30, [45,51], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 31, [46,52], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 32, [47,53], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 33, [48,54], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 34, [49,55], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 35, [50,56], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 36, [51,52], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 37, [51,53], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 38, [52,53], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 39, [52,54], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 40, [53,54], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 41, [54,55], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 42, [53,55], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 43, [53,56], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 44, [55,56], ierr); CHKERRQ(ierr)
    end subroutine

  subroutine set_coords_serial(dm)
    type(tDM) :: dm
    real(ireals), pointer:: coords(:)
    type(tVec)           :: coordinates
    integer(iintegers)   :: dimEmbed, coordSize, vStart, vEnd, pStart, pEnd, eStart, eEnd, voff, cStart, cEnd
    type(tPetscSection)  :: coordSection
    integer(iintegers)   :: i

    real(ireals), parameter :: dx=2, dy=2, dz=.5
    real(ireals), parameter :: ds=sqrt(dy**2 - (dx/2)**2)

    call DMGetCoordinateDim(dm, dimEmbed, ierr); CHKERRQ(ierr)

    call DMGetCoordinateSection(dm, coordSection, ierr); CHKERRQ(ierr)

    call PetscSectionSetNumFields(coordSection, 1, ierr); CHKERRQ(ierr)
    call PetscSectionSetUp(coordSection, ierr); CHKERRQ(ierr)
    call PetscSectionSetFieldComponents(coordSection, 0, dimEmbed, ierr); CHKERRQ(ierr)

    call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
    call DMPlexGetDepthStratum (dm, 0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices
    call DMPlexGetDepthStratum (dm, 1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
    call DMPlexGetHeightStratum (dm, 0, cStart, cEnd, ierr); CHKERRQ(ierr) ! edges

    call PetscSectionSetChart(coordSection, vStart, vEnd, ierr);CHKERRQ(ierr)

    do i = vStart, vEnd-1
      call PetscSectionSetDof(coordSection, i, dimEmbed, ierr); CHKERRQ(ierr)
      call PetscSectionSetFieldDof(coordSection, i, 0, dimEmbed, ierr); CHKERRQ(ierr)
    enddo

    call PetscSectionSetUp(coordSection, ierr); CHKERRQ(ierr)
    call PetscSectionGetStorageSize(coordSection, coordSize, ierr); CHKERRQ(ierr)

    call VecCreate(PETSC_COMM_SELF, coordinates, ierr); CHKERRQ(ierr)
    call VecSetSizes(coordinates, coordSize, PETSC_DETERMINE, ierr);CHKERRQ(ierr)
    call VecSetBlockSize(coordinates, dimEmbed, ierr);CHKERRQ(ierr)
    call VecSetType(coordinates, VECSTANDARD, ierr);CHKERRQ(ierr)

    call PetscObjectSetName(coordinates, "coordinates", ierr); CHKERRQ(ierr)

    call VecGetArrayF90(coordinates, coords, ierr); CHKERRQ(ierr)

    i = 45; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [0*dx,0*ds,0*dz]
    i = 46; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [2*dx,0*ds,0*dz]
    i = 47; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [1*dx,2*ds,0*dz]
    i = 48; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [3*dx,2*ds,0*dz]
    i = 49; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [2*dx,4*ds,0*dz]
    i = 50; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [0*dx,4*ds,0*dz]
    i = 51; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [0*dx,0*ds,1*dz]
    i = 52; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [2*dx,0*ds,1*dz]
    i = 53; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [1*dx,2*ds,1*dz]
    i = 54; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [3*dx,2*ds,1*dz]
    i = 55; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [2*dx,4*ds,1*dz]
    i = 56; call PetscSectionGetOffset(coordSection, i, voff, ierr); coords(voff+1:voff+3) = [0*dx,4*ds,1*dz]

    call VecRestoreArrayF90(coordinates, coords, ierr); CHKERRQ(ierr)
    call DMSetCoordinatesLocal(dm, coordinates, ierr);CHKERRQ(ierr)
    call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); CHKERRQ(ierr)
    call VecDestroy(coordinates, ierr);CHKERRQ(ierr)
  end subroutine

  subroutine create_face_labels(dm)
    DM, intent(inout) :: dm
    PetscInt    :: i, depth
    PetscSection :: s

    PetscInt :: pStart, pEnd
    PetscInt :: cStart, cEnd
    PetscInt :: fStart, fEnd
    PetscInt :: eStart, eEnd
    PetscInt :: vStart, vEnd

    type(tDMLabel) :: faceposlabel, zindexlabel, TOAlabel
    type(PetscInt) :: lval, zval

    call DMCreateLabel(dm, "Face Position", ierr); CHKERRQ(ierr)
    call DMCreateLabel(dm, "Vertical Index", ierr); CHKERRQ(ierr)
    call DMCreateLabel(dm, "TOA", ierr); CHKERRQ(ierr)

    call DMGetLabel(dm, "Face Position", faceposlabel, ierr); CHKERRQ(ierr)
    call DMGetLabel(dm, "Vertical Index", zindexlabel, ierr); CHKERRQ(ierr)
    call DMGetLabel(dm, "TOA", TOAlabel, ierr); CHKERRQ(ierr)

    call DMGetDefaultSection(dm, s, ierr); CHKERRQ(ierr)
    call DMPlexGetDepth(dm, depth, ierr); CHKERRQ(ierr)
    call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, 0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
    call DMPlexGetHeightStratum(dm, 1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, 1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
    call DMPlexGetDepthStratum (dm, 0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices
    print *,'PLEX GetChart', pStart, pEnd, ":: fStart, fEnd", fStart, fEnd

    do i=fstart, fEnd-1
      select case(i)
      case (4:7)
        lval = TOP_BOT_FACE
        zval = 2
      case (17:20)
        lval = TOP_BOT_FACE
        zval = 1
      case default
        lval = SIDE_FACE
        zval = 1
      end select

      call DMLabelSetValue(faceposlabel, i, lval, ierr); CHKERRQ(ierr)
      call DMLabelSetValue(zindexlabel, i, zval, ierr); CHKERRQ(ierr)

      if (zval.eq.1 .and. lval.eq.TOP_BOT_FACE) call DMLabelSetValue(TOAlabel, i, zval, ierr); CHKERRQ(ierr)
    enddo
    do i=fstart, fEnd-1
      !call DMLabelGetValue(faceposlabel, i, lval, ierr); CHKERRQ(ierr)
      call DMLabelGetValue(TOAlabel, i, lval, ierr); CHKERRQ(ierr)
    enddo
  end subroutine

  end module

  program main
    use m_gen_fish_plex

    call init()

  end program
