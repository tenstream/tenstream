module m_gen_fish_plex
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, init_mpi_data_parameters, &
    i0, i1, i2, i3, i4, i5

  implicit none

  integer(iintegers), parameter :: TOP_BOT_FACE=1, SIDE_FACE=2

  PetscErrorCode :: ierr

  contains
    subroutine init()
      type(tDM) :: dmcell
      character(len=*), parameter :: default_options = '&
        & -show_plex hdf5:fish.h5'

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr); CHKERRQ(ierr)
      call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr); CHKERRQ(ierr)

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
      call DMSetDimension(dm, i3, ierr);CHKERRQ(ierr)

      call DMPlexSetChart(dm, i0, 57_iintegers, ierr); CHKERRQ(ierr)

      call setup_plex(dm)

    end subroutine
    subroutine setup_plex(dm)
      DM :: dm

      integer(iintegers) :: pStart, pEnd
      integer(iintegers) :: cStart, cEnd
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: eStart, eEnd
      integer(iintegers) :: vStart, vEnd

      call set_wedge_connectivity(dm)

      call DMPlexSymmetrize(dm, ierr); CHKERRQ(ierr)
      call DMPlexStratify(dm, ierr); CHKERRQ(ierr)

      call set_coords_serial(dm)

      call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
      call DMPlexGetHeightStratum(dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
      call DMPlexGetHeightStratum(dm, i1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
      call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
      call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices

      print *,'pStart,End :: ',pStart, pEnd
      print *,'cStart,End :: ',cStart, cEnd
      print *,'fStart,End :: ',fStart, fEnd
      print *,'eStart,End :: ',eStart, eEnd
      print *,'vStart,End :: ',vStart, vEnd

      call create_face_labels(dm)
    end subroutine

    subroutine set_wedge_connectivity(dm)
      DM :: dm
      integer(iintegers) :: i
      ! Preallocation
      ! Every cell has 5 faces
      do i=0,3
        call DMPlexSetConeSize(dm, i, i5, ierr); CHKERRQ(ierr)
      enddo

      ! Faces have 3 or 4 edges
      do i=4,7
        call DMPlexSetConeSize(dm, i, i3, ierr); CHKERRQ(ierr)
      enddo
      do i=17,20
        call DMPlexSetConeSize(dm, i, i3, ierr); CHKERRQ(ierr)
      enddo
      do i=8,16
        call DMPlexSetConeSize(dm, i, i4, ierr); CHKERRQ(ierr)
      enddo

      ! Edges have 2 vertices
      do i=21,44
        call DMPlexSetConeSize(dm, i, i2, ierr); CHKERRQ(ierr)
      enddo

      call DMSetUp(dm, ierr); CHKERRQ(ierr) ! Allocate space for cones

      ! Setup Connections
      call DMPlexSetCone(dm,  0_iintegers, [ 4_iintegers, 8_iintegers, 9_iintegers,10_iintegers,17_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  1_iintegers, [ 5_iintegers,10_iintegers,11_iintegers,12_iintegers,18_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  2_iintegers, [ 6_iintegers,12_iintegers,13_iintegers,14_iintegers,19_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  3_iintegers, [ 7_iintegers,14_iintegers,15_iintegers,16_iintegers,20_iintegers], ierr); CHKERRQ(ierr)

      call DMPlexSetCone(dm,  4_iintegers, [21_iintegers,22_iintegers,23_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  5_iintegers, [23_iintegers,24_iintegers,25_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  6_iintegers, [25_iintegers,26_iintegers,27_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  7_iintegers, [27_iintegers,28_iintegers,29_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 17_iintegers, [36_iintegers,37_iintegers,38_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 18_iintegers, [38_iintegers,39_iintegers,40_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 19_iintegers, [40_iintegers,41_iintegers,42_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 20_iintegers, [42_iintegers,43_iintegers,44_iintegers], ierr); CHKERRQ(ierr)

      call DMPlexSetCone(dm,  8_iintegers, [21_iintegers,30_iintegers,31_iintegers,36_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm,  9_iintegers, [22_iintegers,30_iintegers,32_iintegers,37_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 10_iintegers, [23_iintegers,31_iintegers,32_iintegers,38_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 11_iintegers, [24_iintegers,31_iintegers,33_iintegers,39_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 12_iintegers, [25_iintegers,32_iintegers,33_iintegers,40_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 13_iintegers, [26_iintegers,33_iintegers,34_iintegers,41_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 14_iintegers, [27_iintegers,32_iintegers,34_iintegers,42_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 15_iintegers, [28_iintegers,32_iintegers,35_iintegers,43_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 16_iintegers, [29_iintegers,34_iintegers,35_iintegers,44_iintegers], ierr); CHKERRQ(ierr)

      call DMPlexSetCone(dm, 21_iintegers, [45_iintegers,46_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 22_iintegers, [45_iintegers,47_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 23_iintegers, [46_iintegers,47_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 24_iintegers, [46_iintegers,48_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 25_iintegers, [47_iintegers,48_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 26_iintegers, [48_iintegers,49_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 27_iintegers, [47_iintegers,49_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 28_iintegers, [47_iintegers,50_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 29_iintegers, [49_iintegers,50_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 30_iintegers, [45_iintegers,51_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 31_iintegers, [46_iintegers,52_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 32_iintegers, [47_iintegers,53_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 33_iintegers, [48_iintegers,54_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 34_iintegers, [49_iintegers,55_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 35_iintegers, [50_iintegers,56_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 36_iintegers, [51_iintegers,52_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 37_iintegers, [51_iintegers,53_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 38_iintegers, [52_iintegers,53_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 39_iintegers, [52_iintegers,54_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 40_iintegers, [53_iintegers,54_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 41_iintegers, [54_iintegers,55_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 42_iintegers, [53_iintegers,55_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 43_iintegers, [53_iintegers,56_iintegers], ierr); CHKERRQ(ierr)
      call DMPlexSetCone(dm, 44_iintegers, [55_iintegers,56_iintegers], ierr); CHKERRQ(ierr)
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

    call PetscSectionSetNumFields(coordSection, i1, ierr); CHKERRQ(ierr)
    call PetscSectionSetUp(coordSection, ierr); CHKERRQ(ierr)
    call PetscSectionSetFieldComponents(coordSection, i0, dimEmbed, ierr); CHKERRQ(ierr)

    call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
    call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices
    call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
    call DMPlexGetHeightStratum (dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr) ! edges

    call PetscSectionSetChart(coordSection, vStart, vEnd, ierr);CHKERRQ(ierr)

    do i = vStart, vEnd-1
      call PetscSectionSetDof(coordSection, i, dimEmbed, ierr); CHKERRQ(ierr)
      call PetscSectionSetFieldDof(coordSection, i, i0, dimEmbed, ierr); CHKERRQ(ierr)
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
    integer(iintegers) :: i, depth
    PetscSection :: s

    integer(iintegers) :: pStart, pEnd
    integer(iintegers) :: cStart, cEnd
    integer(iintegers) :: fStart, fEnd
    integer(iintegers) :: eStart, eEnd
    integer(iintegers) :: vStart, vEnd

    type(tDMLabel) :: faceposlabel, zindexlabel, TOAlabel
    integer(iintegers) :: lval, zval

    call DMCreateLabel(dm, "Face Position", ierr); CHKERRQ(ierr)
    call DMCreateLabel(dm, "Vertical Index", ierr); CHKERRQ(ierr)
    call DMCreateLabel(dm, "TOA", ierr); CHKERRQ(ierr)

    call DMGetLabel(dm, "Face Position", faceposlabel, ierr); CHKERRQ(ierr)
    call DMGetLabel(dm, "Vertical Index", zindexlabel, ierr); CHKERRQ(ierr)
    call DMGetLabel(dm, "TOA", TOAlabel, ierr); CHKERRQ(ierr)

    call DMGetDefaultSection(dm, s, ierr); CHKERRQ(ierr)
    call DMPlexGetDepth(dm, depth, ierr); CHKERRQ(ierr)
    call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
    call DMPlexGetHeightStratum(dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr) ! cells
    call DMPlexGetHeightStratum(dm, i1, fStart, fEnd, ierr); CHKERRQ(ierr) ! faces / edges
    call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
    call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices
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
