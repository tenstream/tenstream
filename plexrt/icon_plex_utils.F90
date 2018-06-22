module m_icon_plex_utils

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, i0, i1, i2, i3

  use m_helper_functions, only: chkerr, itoa

  use m_plex_grid, only: print_dmplex

  implicit none

  private
  public :: dmplex_2D_to_3D, create_2d_fish_plex

  logical, parameter :: ldebug=.True.

  contains

    subroutine dmplex_2D_to_3D(dm2d, hhl, dm3d)
      type(tDM), intent(in) :: dm2d
      real(ireals), intent(in) :: hhl(:) ! height levels of interfaces, those will be added to base height of 2D elements
      type(tDM), intent(out) :: dm3d

      integer(mpiint) :: comm, ierr

      call PetscObjectGetComm(dm2d, comm, ierr); call CHKERR(ierr)

      call print_dmplex(comm, dm2d)

    end subroutine

    ! Create a 2D Torus grid with Nx vertices horizontally and Ny rows of Vertices vertically
    subroutine create_2d_fish_plex(dm, Nx, Ny, lcyclic)
      type(tDM) :: dm, dmdist
      integer(iintegers), intent(in) :: Nx, Ny
      logical, intent(in) :: lcyclic

      integer(iintegers) :: chartsize, Nfaces, Nedges, Nvertices

      integer(iintegers) :: pStart, pEnd
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: eStart, eEnd
      integer(iintegers) :: vStart, vEnd

      integer(mpiint) :: myid, numnodes, ierr

      logical, parameter :: ldebug=.False.

      call mpi_comm_rank(PETSC_COMM_WORLD, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(PETSC_COMM_WORLD, numnodes, ierr); call CHKERR(ierr)

      call DMPlexCreate(PETSC_COMM_WORLD, dm, ierr); call CHKERR(ierr)
      call PetscObjectSetName(dm, 'testplex Nx'//itoa(Nx)//'_Ny'//itoa(Ny), ierr); call CHKERR(ierr)
      call DMSetDimension(dm, i2, ierr); call CHKERR(ierr)

      if(modulo(Nx,2).ne.0) call CHKERR(1_mpiint, 'Nx has to be even, e.g. 2,4,6...')
      if(modulo(Ny,2).eq.0) call CHKERR(1_mpiint, 'Nx has to be uneven, e.g. 3,5...')

      if(myid.eq.0) then
        if(lcyclic) then
          Nfaces = (Nx-1) * (Ny-1) * 2        ! Number of faces
          Nedges = (Nx-1) * (Ny-1)            ! Number of edges on full height (y = 0, ...)
          Nedges = Nedges + (Nx*2-2) * (Ny-1) ! Number of edges on half height (y = {} + .5)
          Nvertices = (Nx-1) * (Ny-1)         ! Number of vertices
        else
          Nfaces = (Nx-1) * (Ny-1) * 2        ! Number of faces
          Nedges = (Nx-1) * Ny                ! Number of edges on full height (y = 0, ...)
          Nedges = Nedges + (Nx*2-1) * (Ny-1) ! Number of edges on half height (y = {} + .5)
          Nvertices = Nx * Ny                 ! Number of vertices
        endif
      else
        Nfaces = 0
        Nedges = 0
        Nvertices = 0
      endif
      if(ldebug) then
        print *, myid, 'Nfaces', Nfaces
        print *, myid, 'Nedges', Nedges
        print *, myid, 'Nverts', Nvertices
      endif

      chartsize = Nfaces + Nedges + Nvertices
      call DMPlexSetChart(dm, i0, chartsize, ierr); call CHKERR(ierr)

      call set_wedge_connectivity(dm, Nx, Ny, Nfaces, Nedges, Nvertices, lcyclic)

      call DMPlexSymmetrize(dm, ierr); CHKERRQ(ierr)
      call DMPlexStratify(dm, ierr); CHKERRQ(ierr)

      call set_coords_serial(dm, Nx, Ny, lcyclic)

      call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(dm, i0, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
      call DMPlexGetHeightStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetHeightStratum(dm, i2, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

      if(myid.eq.0.and.ldebug) then
        print *,'pStart,End :: ',pStart, pEnd
        print *,'fStart,End :: ',fStart, fEnd
        print *,'eStart,End :: ',eStart, eEnd
        print *,'vStart,End :: ',vStart, vEnd
      endif

      call DMPlexSetAdjacencyUseCone(dm, PETSC_TRUE, ierr); call CHKERR(ierr)
      call DMPlexSetAdjacencyUseClosure(dm, PETSC_FALSE, ierr); call CHKERR(ierr)

      call DMPlexDistribute(dm, i0, PETSC_NULL_SF, dmdist, ierr); call CHKERR(ierr)
      if(dmdist.ne.PETSC_NULL_DM) then
        call DMDestroy(dm, ierr); call CHKERR(ierr)
        dm   = dmdist
      endif

      contains

        subroutine set_wedge_connectivity(dm, Nx, Ny, Nfaces, Nedges, Nvertices, lcyclic)
          type(tDM) :: dm
          integer(iintegers) :: Nx, Ny, Nfaces, Nedges, Nvertices
          logical, intent(in) :: lcyclic

          integer(iintegers) :: k, i, j, ioff, cone3(3), cone2(2)
          ! Preallocation
          ! Faces have 3 edges
          ioff = 0
          do k = 1, Nfaces
            call DMPlexSetConeSize(dm, ioff, i3, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
          ! Edges have 2 vertices
          do k = 1, Nedges
            call DMPlexSetConeSize(dm, ioff, i2, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo

          call DMSetUp(dm, ierr); call CHKERR(ierr) ! Allocate space for cones

          ioff = 0
          do k = 1, Nfaces
            j = (k-1) / (2*(Nx-1)) ! row of faces
            i = k-1 - j*(2*(Nx-1)) ! col of faces
            !print *,'faces',k,': i,j', i, j

            ! determine edges of a face
            if(modulo(i+modulo(j,2),2).eq.0) then ! this has a bot edge
              if(lcyclic) then
                cone3(1) = Nfaces + j*(Nx-1) + j*(Nx*2-2) + i/2  ! Nfaces offset + number of edges of full height + number of edges of half heights + i offset
                cone3(2) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-2) + i  ! left edge  ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
                cone3(3) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-2) + i +1 ! right edge ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
                if(i.eq.(Nx-1)*2-1) cone3(3) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-2)
              else
                cone3(1) = Nfaces + j*(Nx-1) + j*(Nx*2-1) + i/2  ! Nfaces offset + number of edges of full height + number of edges of half heights + i offset
                cone3(2) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-1) + i  ! left edge  ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
                cone3(3) = Nfaces + (j+1)*(Nx-1) + j*(Nx*2-1) + i +1 ! right edge ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
              endif
              if(ldebug) print *,'upward edge of face', k-1, ':', cone3
            else
              if(lcyclic) then
                cone3(1) = Nfaces + (j+1)*(Nx-1) + (j+1)*(Nx*2-2) + i/2
                cone3(2) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-2) + i
                cone3(3) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-2) + i +1
                if(i.eq.(Nx-1)*2-1) cone3(3) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-2)
                if(j.eq.Ny-2) cone3(1) = Nfaces + i/2
              else
                cone3(1) = Nfaces + (j+1)*(Nx-1) + (j+1)*(Nx*2-1) + i/2
                cone3(2) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-1) + i
                cone3(3) = Nfaces + (j+1)*(Nx-1) + (j+0)*(Nx*2-1) + i +1
              endif
              if(ldebug) print *,'                                                                            downward edge of face', k-1, ':', i, j,':', cone3
            endif

            call DMPlexSetCone(dm,  ioff, cone3, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
          ! Edges have 2 vertices
          do k = 1, Nedges
            if(lcyclic) then
              j = (k-1) / (Nx-1 + Nx*2-2) ! row of edge ! goes up to Ny
              i = k-1 - j*(Nx-1 + Nx*2-2) ! col of edge
              if(i.lt.Nx-1) then ! bottom edge
                cone2(1) = Nfaces + Nedges + j*(Nx-1) + i
                cone2(2) = Nfaces + Nedges + j*(Nx-1) + i +1
                if(i.eq.Nx-2) cone2(2) = Nfaces + Nedges + j*(Nx-1)
                if(ldebug) print *,'bot edge',k,ioff,': i,j', i, j, ':', cone2
              else ! sideward edge
                if(modulo(i+j,2).ne.0) then ! slash
                  cone2(1) = Nfaces + Nedges + j*(Nx-1) + (i-Nx+1)/2
                  cone2(2) = Nfaces + Nedges + (j+1)*(Nx-1) + (i-Nx+1)/2 + modulo(j,2)
                  if(i.eq.(Nx-1+Nx*2-2-1)) cone2(2) = Nfaces + Nedges + (j+1)*(Nx-1)
                  if(j.eq.Ny-2) then
                    cone2(2) = Nfaces + Nedges + (i-Nx+1)/2 + modulo(j,2)
                    if(i.eq.(Nx-1+Nx*2-2-1)) cone2(2) = Nfaces + Nedges
                  endif
                  if(ldebug) print *,'                                                                            slash side edge',k,ioff,': i,j', i, j, ':', cone2
                else ! backslash
                  cone2(1) = Nfaces + Nedges + j*(Nx-1) + (i-Nx+1)/2 + modulo(j+1,2)
                  cone2(2) = Nfaces + Nedges + (j+1)*(Nx-1) + (i-Nx+1)/2
                  if(i.eq.(Nx-1+Nx*2-2-1)) cone2(1) = Nfaces + Nedges + j*(Nx-1)
                  if(j.eq.Ny-2) cone2(2) = Nfaces + Nedges + (i-Nx+1)/2
                  if(ldebug) print *,'                                                                            backs side edge',k,ioff,': i,j', i, j, ':', cone2
                endif
              endif
            else
              j = (k-1) / (Nx-1 + Nx*2-1) ! row of edge ! goes up to Ny
              i = k-1 - j*(Nx-1 + Nx*2-1) ! col of edge
              if(i.lt.Nx-1) then ! bottom edge
                cone2(1) = Nfaces + Nedges + j*Nx + i
                cone2(2) = Nfaces + Nedges + j*Nx + i +1
                if(ldebug) print *,'bot edge',k,ioff,': i,j', i, j, ':', cone2
              else ! sideward edge
                if(modulo(i+j,2).ne.0) then ! slash
                  cone2(1) = Nfaces + Nedges + j*Nx + (i-Nx+1)/2
                  cone2(2) = Nfaces + Nedges + (j+1)*Nx + (i-Nx+1)/2 + modulo(j,2)
                  !print *,'                                                                            slash side edge',k,ioff,': i,j', i, j, ':', cone2
                else ! backslash
                  cone2(1) = Nfaces + Nedges + j*Nx + (i-Nx+1)/2 + modulo(j+1,2)
                  cone2(2) = Nfaces + Nedges + (j+1)*Nx + (i-Nx+1)/2
                  !print *,'                                                                            backs side edge',k,ioff,': i,j', i, j, ':', cone2
                endif
              endif
            endif
            call DMPlexSetCone(dm,  ioff, cone2, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
        end subroutine

        subroutine set_coords_serial(dm, Nx, Ny, lcyclic)
          type(tDM) :: dm
          integer(iintegers), intent(in) :: Nx, Ny
          logical, intent(in) :: lcyclic
          real(ireals), pointer:: coords(:)
          type(tVec)           :: coordinates
          integer(iintegers)   :: dimEmbed, coordSize, vStart, vEnd, pStart, pEnd, eStart, eEnd, voff, cStart, cEnd
          type(tPetscSection)  :: coordSection
          integer(iintegers)   :: iv, i, j

          real(ireals), parameter :: dx=1, dy=1
          real(ireals), parameter :: ds=sqrt(dy**2 - (dx/2)**2)
          real(ireals) :: x, y, z

          call DMGetCoordinateDim(dm, dimEmbed, ierr); CHKERRQ(ierr)
          dimEmbed = 3

          call DMGetCoordinateSection(dm, coordSection, ierr); CHKERRQ(ierr)

          call PetscSectionSetNumFields(coordSection, i1, ierr); CHKERRQ(ierr)
          call PetscSectionSetUp(coordSection, ierr); CHKERRQ(ierr)
          call PetscSectionSetFieldComponents(coordSection, i0, dimEmbed, ierr); CHKERRQ(ierr)

          call DMPlexGetChart(dm, pStart, pEnd, ierr); CHKERRQ(ierr)
          call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); CHKERRQ(ierr) ! vertices
          call DMPlexGetDepthStratum (dm, i1, eStart, eEnd, ierr); CHKERRQ(ierr) ! edges
          call DMPlexGetHeightStratum (dm, i0, cStart, cEnd, ierr); CHKERRQ(ierr) ! faces

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

          do iv = vStart, vEnd-1
            if(lcyclic) then
              j = (iv-vStart) / (Nx-1)
              i = (iv-vStart) - j*(Nx-1)
              z = sqrt(x**2 + y**2)
            else
              j = (iv-vStart) / Nx
              i = (iv-vStart) - j*Nx
              z = 0
            endif
            x = (i+.5*modulo(j,2))*dx
            y = j*ds
            if(ldebug) print *,'iv',iv,':', i, j,'=>', x, y, z

            call PetscSectionGetOffset(coordSection, iv, voff, ierr); coords(voff+1:voff+3) = [x, y, z]
          enddo

          call VecRestoreArrayF90(coordinates, coords, ierr); CHKERRQ(ierr)
          call DMSetCoordinatesLocal(dm, coordinates, ierr);CHKERRQ(ierr)
          call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); CHKERRQ(ierr)
          call VecDestroy(coordinates, ierr);CHKERRQ(ierr)
        end subroutine
      end subroutine
end module
