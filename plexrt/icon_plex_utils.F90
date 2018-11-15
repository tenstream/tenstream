module m_icon_plex_utils

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    i0, i1, i2, i3, i4, i5, zero, one, default_str_len, pi

  use m_helper_functions, only: chkerr, itoa, norm, get_arg, imp_bcast, deg2rad

  use m_plex_grid, only: t_plexgrid, print_dmplex, create_plex_section, TOAFACE, &
    get_horizontal_faces_around_vertex

  use m_netcdfio, only: ncload

  implicit none

  private
  public :: dmplex_2D_to_3D, create_2d_fish_plex, dump_ownership, &
    gen_2d_plex_from_icongridfile, icon_hdcp2_default_hhl, &
    icon_ncvec_to_plex, Nz_Ncol_vec_to_celldm1, celldm1_vec_to_Nz_Ncol, &
    celldm_veccopy, dm2d_vec_to_Nz_Ncol, date_to_julian_day, get_sun_vector

  !logical, parameter :: ldebug=.True.
  logical, parameter :: ldebug=.False.

  real(ireals), parameter :: icon_hdcp2_default_hhl(151) = &
    [21000.000   , 20617.792   , 20290.235   , 19982.347   , 19687.435   , 19402.393   , 19125.412   , 18855.302   , 18591.217   , 18332.527   , 18078.737   , 17829.452   , 17584.349   , 17343.154   , 17105.635   , 16871.593   , 16640.853   , 16413.261   , 16188.679   , 15966.985   , 15748.067   , 15531.826   , 15318.170   , 15107.015   , 14898.282   , 14691.901   , 14487.805   , 14285.932   , 14086.225   , 13888.630   , 13693.095   , 13499.575   , 13308.022   , 13118.397   , 12930.657   , 12744.766   , 12560.687   , 12378.386   , 12197.831   , 12018.990   , 11841.834   , 11666.334   , 11492.464   , 11320.198   , 11149.510   , 10980.378   , 10812.778   , 10646.689   , 10482.089   , 10318.958   , 10157.277   , 9997.027    , 9838.190    , 9680.749    , 9524.686    , 9369.987    , 9216.635    , 9064.616    , 8913.914    , 8764.517    , 8616.410    , 8469.580    , 8324.016    , 8179.704    , 8036.633    , 7894.793    , 7754.171    , 7614.758    , 7476.543    , 7339.516    , 7203.669    , 7068.992    , 6935.476    , 6803.113    , 6671.895    , 6541.813    , 6412.861    , 6285.031    , 6158.317    , 6032.711    , 5908.207    , 5784.800    , 5662.484    , 5541.252    , 5421.100    , 5302.024    , 5184.017    , 5067.077    , 4951.198    , 4836.377    , 4722.610    , 4609.895    , 4498.227    , 4387.605    , 4278.026    , 4169.487    , 4061.987    , 3955.525    , 3850.098    , 3745.707    , 3642.350    , 3540.027    , 3438.738    , 3338.483    , 3239.264    , 3141.081    , 3043.935    , 2947.829    , 2852.765    , 2758.746    , 2665.774    , 2573.854    , 2482.990    , 2393.185    , 2304.447    , 2216.779    , 2130.190    , 2044.686    , 1960.274    , 1876.965    , 1794.766    , 1713.690    , 1633.747    , 1554.949    , 1477.311    , 1400.848    , 1325.575    , 1251.512    , 1178.678    , 1107.095    , 1036.786    , 967.780     , 900.104     , 833.792     , 768.881     , 705.412     , 643.431     , 582.990     , 524.151     , 466.982     , 411.564     , 357.994     , 306.385     , 256.878     , 209.648     , 164.919     , 122.997     , 84.314      , 49.554      , 20.000      , 0.000 ]

  contains

    subroutine dmplex_2D_to_3D(dm2d, ke1, hhl, dm3d, zindex)
      type(tDM), intent(in) :: dm2d
      integer(iintegers), intent(in) :: ke1 ! number of levels for the 3D DMPlex
      real(ireals), intent(in) :: hhl(:) ! height levels of interfaces, those will be added to base height of 2D elements, either of shape(nlev) or shape(nlev*nverts)
      type(tDM), intent(out) :: dm3d
      ! vertical layer / level of cells/faces/edges/vertices , pStart..pEnd-1, fortran indexing, i.e. start with k=1
      integer(iintegers), allocatable, intent(out) :: zindex(:)

      integer(iintegers) :: p2dStart, p2dEnd
      integer(iintegers) :: f2dStart, f2dEnd
      integer(iintegers) :: e2dStart, e2dEnd
      integer(iintegers) :: v2dStart, v2dEnd
      integer(iintegers) :: ke, chartsize
      integer(iintegers) :: Nfaces2d, Nedges2d, Nverts2d
      integer(iintegers) :: Ncells, Nfaces, Nedges, Nverts
      integer(mpiint) :: comm, ierr

      ke = ke1-1

      call PetscObjectGetComm(dm2d, comm, ierr); call CHKERR(ierr)

      call print_dmplex(comm, dm2d)

      call DMPlexGetChart(dm2d, p2dStart, p2dEnd, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(dm2d, i0, f2dStart, f2dEnd, ierr); call CHKERR(ierr) ! faces
      call DMPlexGetHeightStratum(dm2d, i1, e2dStart, e2dEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetHeightStratum(dm2d, i2, v2dStart, v2dEnd, ierr); call CHKERR(ierr) ! vertices

      Nfaces2d = f2dEnd - f2dStart
      Nedges2d = e2dEnd - e2dStart
      Nverts2d = v2dEnd - v2dStart

      Ncells = Nfaces2d * ke
      Nfaces = Nfaces2d * ke1 + Nedges2d * ke
      Nedges = Nedges2d * ke1 + Nverts2d * ke
      Nverts = Nverts2d * ke1

      chartsize = Ncells + Nfaces + Nedges + Nverts

      if(ldebug) then
        print *,'Nlev    ', ke1
        print *,'size hhl', size(hhl)
        print *,'Nfaces2d', Nfaces2d
        print *,'Nedges2d', Nedges2d
        print *,'Nverts2d', Nverts2d

        print *,'Ncells3d', Ncells
        print *,'Nfaces3d', Nfaces
        print *,'Nedges3d', Nedges
        print *,'Nverts3d', Nverts
        print *,'Chartsize3d', chartsize
      endif

      call DMPlexCreate(comm, dm3d, ierr); call CHKERR(ierr)
      call DMSetDimension(dm3d, i3, ierr); call CHKERR(ierr)
      call DMPlexSetChart(dm3d, i0, chartsize, ierr); call CHKERR(ierr)

      call set_connectivity(dm2d, dm3d)

      call set_sf_graph(dm2d, dm3d)

      call set_coords(dm2d, dm3d)

      call DMSetFromOptions(dm3d, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(dm2d, PETSC_NULL_DM, "-show_iconplex_2d", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(dm3d, PETSC_NULL_DM, "-show_iconplex_3d", ierr); call CHKERR(ierr)

      if(ldebug) then
        call gen_test_mat(dm3d)
        call gen_test_mat(dm2d)
      endif
      !call CHKERR(1_mpiint, 'DEBUG')

    contains

      subroutine gen_test_mat(dm)
        type(tDM), intent(in) :: dm
        type(tPetscSection) :: sec
        type(tMat) :: A

        call create_plex_section(dm, 'face_test_section', i1, [i0], [i1], [i0], [i0], sec)
        call DMSetSection(dm, sec, ierr); call CHKERR(ierr)
        call PetscSectionDestroy(sec, ierr); call CHKERR(ierr)

        call DMCreateMatrix(dm, A, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(A, PETSC_NULL_MAT, "-show_dmplex_2D_to_3D_test_mat", ierr); call CHKERR(ierr)
        call MatDestroy(A, ierr); call CHKERR(ierr)
      end subroutine

      subroutine set_sf_graph(dm2d, dm3d)
        type(tDM), intent(in) :: dm2d
        type(tDM), intent(inout) :: dm3d

        integer(mpiint) :: myid, ierr

        type(tDM) :: dmsf2d
        type(tPetscSF) :: sf2d, sf3d

        type(tPetscSection) :: section_2d_to_3d
        type(tVec) :: lVec, gVec
        real(ireals), pointer :: xv(:)

        integer(iintegers), pointer :: pmyidx(:) ! list of my indices that we do not own
        type(PetscSFNode), pointer  :: premote(:) ! rank and remote idx of those points
        integer(iintegers), allocatable :: myidx(:)
        type(PetscSFNode), allocatable :: remote(:)
        integer(iintegers) :: nroots2d   ! number of root vertices on the current process (possible targets for leaves)
        integer(iintegers) :: nleaves2d  ! number of leaf vertices on the current process, references roots on any process
        integer(iintegers) :: nroots3d, nleaves3d
        integer(iintegers),allocatable :: ilocal_elements(:)
        type(PetscSFNode),allocatable :: iremote_elements(:)
        PetscCopyMode,parameter :: localmode=PETSC_COPY_VALUES, remotemode=PETSC_COPY_VALUES

        integer(iintegers) :: i, k, voff, ileaf, owner

        call DMClone(dm2d, dmsf2d, ierr); call CHKERR(ierr)

        call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

        call DMGetPointSF(dmsf2d, sf2d, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(sf2d, PETSC_NULL_SF, "-show_plex_sf2d", ierr); call CHKERR(ierr)

        call PetscSFGetGraph(sf2d, nroots2d, nleaves2d, pmyidx, premote, ierr); call CHKERR(ierr)
        if(nroots2d.eq.-i1) then
          nleaves2d = 0
          allocate(myidx(0))
          allocate(remote(0))
        else
          allocate(myidx(nleaves2d), source=pmyidx)
          allocate(remote(nleaves2d), source=premote)
        endif
        pmyidx => NULL(); premote=>NULL()
        call PetscSortIntWithArrayPair(nleaves2d, myidx, remote(:)%rank, remote(:)%index, ierr); call CHKERR(ierr)

        call create_plex_section(dmsf2d, 'plex_2d_to_3d_sf_graph_info', i1, &
          [i0], [i0], [ke1+ke], [ke1+ke], section_2d_to_3d)

        call DMSetSection(dmsf2d, section_2d_to_3d, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(section_2d_to_3d, PETSC_NULL_SECTION, &
          '-show_dm2d_section_2d_to_3d', ierr); call CHKERR(ierr)

        call DMGetLocalVector(dmsf2d, lVec, ierr); call CHKERR(ierr)
        call DMGetGlobalVector(dmsf2d, gVec, ierr); call CHKERR(ierr)
        call VecSet(lVec, zero, ierr); call CHKERR(ierr)

        if(ldebug) then
          print *, myid, 'nroots', nroots2d, 'shape myidx', shape(myidx)
          print *, myid, 'nleaves', nleaves2d, 'shape remote', shape(remote)
          print *, myid, 'myidx', myidx
          print *, myid, 'remote', remote
        endif

        ! Distribute Info for 2D edges,
        call VecGetArrayF90(lVec, xv, ierr); call CHKERR(ierr)
        do i = e2dStart, e2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.lt.zero) then ! only add my local idx number if it belongs to me
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              xv(i1+voff+k) = real(iedge_top_icon_2_plex(i, k), ireals)
            enddo
            do k = 0, ke-1
              xv(i1+ke1+voff+k) = real(iface_side_icon_2_plex(i, k), ireals)
            enddo
          endif
        enddo
        do i = v2dStart, v2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.lt.zero) then
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              xv(i1+voff+k) = real(ivertex_icon_2_plex(i, k), ireals)
            enddo
            do k = 0, ke-1
              xv(i1+ke1+voff+k) = real(iedge_side_icon_2_plex(i, k), ireals)
            enddo
          endif
        enddo
        call VecRestoreArrayF90(lVec, xv, ierr); call CHKERR(ierr)

        call VecSet(gVec, zero, ierr); call CHKERR(ierr)
        call DMLocalToGlobalBegin(dmsf2d, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
        call DMLocalToGlobalEnd(dmsf2d, lVec, ADD_VALUES, gVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalBegin(dmsf2d, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)
        call DMGlobalToLocalEnd(dmsf2d, gVec, INSERT_VALUES, lVec, ierr); call CHKERR(ierr)

        nroots3d = chartsize
        nleaves3d = max(i0, nleaves2d * (ke+ke1))
        allocate(ilocal_elements(nleaves3d))
        allocate(iremote_elements(nleaves3d))

        ileaf = 1
        call VecGetArrayF90(lVec, xv, ierr); call CHKERR(ierr)
        do i = e2dStart, e2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.ge.zero) then ! this is owned by someone else
            owner = remote(i1+voff)%rank
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              ilocal_elements(ileaf) = iedge_top_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = int(xv(i1+voff+k), iintegers)
              !if(ldebug) print *,myid,' 2dEdge top', i,'::', k,' local edge index', iedge_top_icon_2_plex(i, k), &
              !  'remote idx', xv(i1+voff+k)
              ileaf = ileaf+1
            enddo
            do k = 0, ke-1
              ilocal_elements(ileaf) = iface_side_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = int(xv(i1+ke1+voff+k), iintegers)
              !if(ldebug) print *,myid,' 2dEdge fac', i,'::', k,' local face index', ilocal_elements(ileaf), &
              !  'remote idx', xv(i1+ke1+voff+k)
              ileaf = ileaf+1
            enddo
          endif
        enddo
        do i = v2dStart, v2dEnd-1
          call PetscFindInt(i, nleaves2d, myidx, voff, ierr); call CHKERR(ierr)
          if(voff.ge.zero) then ! this is owned by someone else
            owner = remote(i1+voff)%rank
            call PetscSectionGetOffset(section_2d_to_3d, i, voff, ierr); call CHKERR(ierr)
            do k = 0, ke1-1
              ilocal_elements(ileaf) = ivertex_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = int(xv(i1+voff+k), iintegers)
              !if(ldebug) print *,myid,' 2dVert ver', i,'::', k,' local vert index', ilocal_elements(ileaf), &
              !  'remote idx', xv(i1+voff+k)
              ileaf = ileaf+1
            enddo
            do k = 0, ke-1
              ilocal_elements(ileaf) = iedge_side_icon_2_plex(i, k)
              iremote_elements(ileaf)%rank = owner
              iremote_elements(ileaf)%index = int(xv(i1+ke1+voff+k), iintegers)
              !if(ldebug) print *,myid,' 2dVert edg', i,'::', k,' local face index', ilocal_elements(ileaf), &
              !  'remote idx', xv(i1+ke1+voff+k)
              ileaf = ileaf+1
            enddo
          endif
        enddo
        call VecRestoreArrayF90(lVec, xv, ierr); call CHKERR(ierr)
        if(ldebug) call CHKERR(int(ileaf-1-nleaves3d, mpiint), 'Does not add up... something is wrong')

        call PetscObjectViewFromOptions(lVec, PETSC_NULL_VEC, "-show_dm2d_sf_vec", ierr); call CHKERR(ierr)
        call DMRestoreLocalVector(dmsf2d, lVec, ierr); call CHKERR(ierr)
        call DMRestoreGlobalVector(dmsf2d, gVec, ierr); call CHKERR(ierr)
        call PetscSectionDestroy(section_2d_to_3d, ierr); call CHKERR(ierr)
        call DMDestroy(dmsf2d, ierr); call CHKERR(ierr)

        call DMGetPointSF(dm3d, sf3d, ierr); call CHKERR(ierr)
        call PetscSFSetGraph(sf3d, nroots3d, nleaves3d, ilocal_elements, localmode, &
          iremote_elements, remotemode, ierr); call CHKERR(ierr)
        call PetscSFSetUp(sf3d, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(sf3d, PETSC_NULL_SF, "-show_plex_sf3d", ierr); call CHKERR(ierr)
      end subroutine

      subroutine set_connectivity(dm2d, dm3d)
        type(tDM), intent(in) :: dm2d
        type(tDM), intent(inout) :: dm3d

        integer(iintegers) :: i, j, k
        integer(iintegers) :: icell, iface, iedge, ivert
        integer(iintegers), pointer :: cone(:)

        integer(iintegers) :: edge3(3), edge4(4), faces(5), vert2(2)

        ! Preallocation
        k = 0 ! k is running index for elements in DAG

        ! Every cell has 5 faces
        do i = 1, Ncells
          call DMPlexSetConeSize(dm3d, k, i5, ierr); call CHKERR(ierr)
          k = k+1
        enddo

        ! top/bottom faces have 3 edges
        do i = 1, Nfaces2d * ke1
          call DMPlexSetConeSize(dm3d, k, i3, ierr); call CHKERR(ierr)
          k = k+1
        enddo

        ! side faces have 4 edges
        do i = 1, Nedges2d * ke
          call DMPlexSetConeSize(dm3d, k, i4, ierr); call CHKERR(ierr)
          k = k+1
        enddo

        ! Edges have 2 vertices
        do i = 1, Nedges
          call DMPlexSetConeSize(dm3d, k, i2, ierr); call CHKERR(ierr)
          k = k+1
        enddo
        call CHKERR(int((chartsize-Nverts)-k, mpiint), 'This does not add up, we forgot something?')

        call DMSetUp(dm3d, ierr); call CHKERR(ierr) ! Allocate space for cones

        allocate(zindex(i0:chartsize-i1))
        if(ldebug) then
          zindex(:) = -1
        endif

        ! Setup Connections
        ! First set five faces of cell
        do k = 0, ke-1
          do i = 0, Nfaces2d-1
            icell = icell_icon_2_plex(i, k)
            faces(1) = iface_top_icon_2_plex(i, k)
            faces(2) = iface_top_icon_2_plex(i, k+1)

            call DMPlexGetCone(dm2d, i, cone, ierr); call CHKERR(ierr) ! edges of face
            !if(ldebug) print *,'iface2d', i, 'has edges', cone
            do j=1,size(cone)
              faces(2+j) = iface_side_icon_2_plex(cone(j), k)
            enddo
            call DMPlexRestoreCone(dm2d, i, cone, ierr); call CHKERR(ierr)

            call DMPlexSetCone(dm3d, icell, faces, ierr); call CHKERR(ierr)
            !if(ldebug) print *,'icell',icell,'has faces:',faces
            zindex(icell) = k+1
          enddo
        enddo

        ! set edges of top/bot faces
        do k = 0, ke1-1
          do i = 0, Nfaces2d-1
            iface = iface_top_icon_2_plex(i, k)

            call DMPlexGetCone(dm2d, i, cone, ierr); call CHKERR(ierr) ! edges of face
            do j=1,size(cone)
              edge3(j) = iedge_top_icon_2_plex(cone(j), k)
            enddo
            call DMPlexRestoreCone(dm2d, i, cone, ierr); call CHKERR(ierr)

            call DMPlexSetCone(dm3d, iface, edge3, ierr); call CHKERR(ierr)
            !if(ldebug) print *,'iface2d', i, 'iface3d', iface, 'gets edges3', edge3
            zindex(iface) = k+1
          enddo
        enddo

        ! set edges of vertical faces
        do k = 0, ke-1
          do i = 0, Nedges2d-1
            iedge = e2dStart+i
            iface = iface_side_icon_2_plex(iedge, k)
            edge4(1) = iedge_top_icon_2_plex(iedge, k)
            edge4(2) = iedge_top_icon_2_plex(iedge, k+1)

            call DMPlexGetCone(dm2d, iedge, cone, ierr); call CHKERR(ierr) ! vertices of 2d edge
            do j=1,size(cone)
              edge4(2+j) = iedge_side_icon_2_plex(cone(j), k)
            enddo
            !if(ldebug) print *,'below_edge2d', iedge, 'iface3d', iface, 'gets edges', edge4, 'cone',cone
            call DMPlexRestoreCone(dm2d, iedge, cone, ierr); call CHKERR(ierr)
            call DMPlexSetCone(dm3d, iface, edge4, ierr); call CHKERR(ierr)
            zindex(iface) = k+1
          enddo
        enddo

        ! and then set the two vertices of edges in each level, i.e. vertices for edges in horizontal plane
        do k = 0, ke1-1
          do i = 0, Nedges2d-1
            iedge = e2dStart+i
            call DMPlexGetCone(dm2d, iedge, cone, ierr); call CHKERR(ierr) ! vertices of 2d edge
            do j=1,size(cone)
              vert2(j) = ivertex_icon_2_plex(cone(j), k)
            enddo
            call DMPlexRestoreCone(dm2d, iedge, cone, ierr); call CHKERR(ierr)
            iedge = iedge_top_icon_2_plex(iedge, k)
            !if(ldebug) print *,'iedge', iedge, 'gets verts', vert2
            call DMPlexSetCone(dm3d, iedge, vert2, ierr); call CHKERR(ierr)
            zindex(iedge) = k+1
          enddo
        enddo

        ! and then set the two vertices of edges in each layer, i.e. vertices at the end of vertical edges
        do k = 0, ke-1
          do i = 0, Nverts2d-1
            ivert = v2dStart+i
            iedge = iedge_side_icon_2_plex(ivert, k)
            vert2(1) = ivertex_icon_2_plex(ivert, k)
            vert2(2) = ivertex_icon_2_plex(ivert, k+1)

            !if(ldebug) print *,'edge', iedge, 'gets verts', vert2
            call DMPlexSetCone(dm3d, iedge, vert2, ierr); call CHKERR(ierr)
            zindex(iedge) = k+1
          enddo
        enddo

        do k = 0, ke1-1
          do i = 0, Nverts2d-1
            ivert = v2dStart+i
            ivert = ivertex_icon_2_plex(ivert, k)
            zindex(ivert) = k+1
          enddo
        enddo

        if(ldebug) then
          if(any(zindex.eq.-1)) then
            print *,'zindex',zindex
            call CHKERR(1_mpiint, 'Seems we forgot to set zindex for some elements!')
          endif
        endif

        call DMPlexSymmetrize(dm3d, ierr); call CHKERR(ierr)
        call DMPlexStratify(dm3d, ierr); call CHKERR(ierr)
      end subroutine

      subroutine set_coords(dm2d, dm3d)
        type(tDM), intent(in) :: dm2d
        type(tDM), intent(inout) :: dm3d

        real(ireals), pointer :: coords2d(:), coords3d(:)
        type(tVec)            :: vec_coord2d, vec_coord3d
        type(tPetscSection)   :: coordSection2d, coordSection3d
        integer(iintegers)    :: coordsize, voff2d, voff3d
        integer(iintegers)    :: v2dStart, v2dEnd
        integer(iintegers)    :: v3dStart, v3dEnd

        real(ireals) :: distance, inv_distance, vert_height
        integer(mpiint) :: ierr
        integer(iintegers) :: i, k, ivertex
        logical :: lpolar, lflg, lhave_3d_surface_heights

        if(size(hhl).eq.ke1) then
          lhave_3d_surface_heights = .False.
        else if(size(hhl).eq.ke1*Nverts2d) then
          lhave_3d_surface_heights = .True.
        else
          call CHKERR(1_mpiint, 'heights array has the wrong shape: is' &
            //itoa(size(hhl))//' should be '//itoa(ke1)//' or '//itoa(ke1*Nverts2d))
        endif
        !print *,'size(hhl)', size(hhl), ':', ke1*Nverts2d, ':', lhave_3d_surface_heights

        call DMPlexGetDepthStratum (dm3d, i0, v3dStart, v3dEnd, ierr); call CHKERR(ierr) ! 3D vertices
        call DMPlexGetDepthStratum (dm2d, i0, v2dStart, v2dEnd, ierr); call CHKERR(ierr) ! 2D vertices

        lpolar=.True.
        call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-polar_coords', lpolar, lflg, ierr); call CHKERR(ierr)

        ! Create Coordinate stuff for 3D DM
        call DMGetCoordinateSection(dm3d, coordSection3d, ierr); call CHKERR(ierr)
        call PetscSectionSetNumFields(coordSection3d, i1, ierr); call CHKERR(ierr) ! just 1 field for spherical coords
        call PetscSectionSetUp(coordSection3d, ierr); call CHKERR(ierr)
        call PetscSectionSetFieldComponents(coordSection3d, i0, i3, ierr); call CHKERR(ierr)

        call PetscSectionSetChart(coordSection3d, v3dStart, v3dEnd, ierr);call CHKERR(ierr)
        do i = v3dStart, v3dEnd-1
          call PetscSectionSetDof(coordSection3d, i, i3, ierr); call CHKERR(ierr)
          call PetscSectionSetFieldDof(coordSection3d, i, i0, i3, ierr); call CHKERR(ierr)
        enddo
        call PetscSectionSetUp(coordSection3d, ierr); call CHKERR(ierr)
        call PetscSectionGetStorageSize(coordSection3d, coordSize, ierr); call CHKERR(ierr)

        ! Create New Vec to hold 3D coordinates
        call VecCreate(PETSC_COMM_SELF, vec_coord3d, ierr); call CHKERR(ierr)
        call VecSetSizes(vec_coord3d, coordSize, PETSC_DETERMINE, ierr);call CHKERR(ierr)
        call VecSetBlockSize(vec_coord3d, i3, ierr);call CHKERR(ierr)
        call VecSetType(vec_coord3d, VECSTANDARD, ierr);call CHKERR(ierr)

        call PetscObjectSetName(vec_coord3d, "coordinates", ierr); call CHKERR(ierr)

        ! Fill 3D Coord Vec
        call VecGetArrayF90(vec_coord3d, coords3d, ierr); call CHKERR(ierr)

        ! Get Coordinates from 2D DMPLEX
        call DMGetCoordinateSection(dm2d, coordSection2d, ierr); call CHKERR(ierr)
        call DMGetCoordinatesLocal(dm2d, vec_coord2d, ierr); call CHKERR(ierr)
        call VecGetArrayReadF90(vec_coord2d, coords2d, ierr); call CHKERR(ierr)

        do i = v2dStart, v2dEnd-1
          call PetscSectionGetOffset(coordSection2d, i, voff2d, ierr); call CHKERR(ierr)

          do k = 0, ke1-1
            ivertex = ivertex_icon_2_plex(i, k)
            call PetscSectionGetOffset(coordSection3d, ivertex, voff3d, ierr); call CHKERR(ierr)

            if(lhave_3d_surface_heights) then
              vert_height = hhl(i1 + (i-v2dStart)*ke1 + k)
            else
              vert_height = hhl(i1+k)
            endif

            if(lpolar) then
              distance = norm(coords2d(voff2d+i1 : voff2d+i3))
              inv_distance = one / distance

              coords3d(voff3d+1:voff3d+3) = coords2d(voff2d+1:voff2d+3) * (distance + vert_height) * inv_distance
              !if(ldebug) print *,'Setting coord for 2d', i, '3d vert', ivertex,':', &
              !  coords2d(voff2d+1:voff2d+3), '=>', &
              !  coords3d(voff3d+1:voff3d+3), &
              !  '(',(distance + hhl(k+1)) * inv_distance,')'
            else
              coords3d(voff3d+1:voff3d+3) = coords2d(voff2d+1:voff2d+3) + [zero, zero, vert_height]
            endif
          enddo
        enddo

        call VecRestoreArrayReadF90(vec_coord2d, coords2d, ierr); call CHKERR(ierr)
        call VecRestoreArrayF90(vec_coord3d, coords3d, ierr); call CHKERR(ierr)

        call DMSetCoordinatesLocal(dm3d, vec_coord3d, ierr);call CHKERR(ierr)
        call PetscObjectViewFromOptions(vec_coord2d, PETSC_NULL_VEC, "-show_plex_coordinates2d", ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(vec_coord3d, PETSC_NULL_VEC, "-show_plex_coordinates3d", ierr); call CHKERR(ierr)
        call VecDestroy(vec_coord3d, ierr);call CHKERR(ierr)
      end subroutine

      !> @brief return the dmplex cell index for an  2d grid face index
      function icell_icon_2_plex(iface, k)
        integer(iintegers),intent(in) :: iface    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: icell_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        icell_icon_2_plex = k + iface*ke
      end function

      !> @brief return the dmplex face index for an 2d face index situated at the top of a cell
      function iface_top_icon_2_plex(iface, k)
        integer(iintegers),intent(in) :: iface    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iface_top_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells
        iface_top_icon_2_plex = offset + k + iface*ke1
      end function

      !> @brief return the dmplex face index for an 2d edge index situated at the side of a cell, i.e. below a certain edge
      function iface_side_icon_2_plex(iedge, k)
        integer(iintegers),intent(in) :: iedge    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iface_side_icon_2_plex   !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells + Nfaces2d*ke1
        iface_side_icon_2_plex = offset + k + (iedge-e2dStart)*ke
      end function

      !> @brief return the dmplex edge index for a given 2d edge index, i.e. the edges on the top/bot faces of cells
      function iedge_top_icon_2_plex(iedge, k)
        integer(iintegers),intent(in) :: iedge    !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iedge_top_icon_2_plex !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells + Nfaces
        iedge_top_icon_2_plex = offset + k + (iedge-e2dStart)*ke1
      end function

      !> @brief return the dmplex edge index for a given 2d vertex index, i.e. the edges on the side faces of cells
      function iedge_side_icon_2_plex(ivertex, k)
        integer(iintegers),intent(in) :: ivertex  !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: iedge_side_icon_2_plex !< @param[out] icell_icon_2_plex, the cell index in the dmplex, starts from 0 and goes to plex%cEnd
        integer(iintegers) :: offset
        offset = Ncells + Nfaces + Nedges2d*ke1
        iedge_side_icon_2_plex = offset + k + (ivertex-v2dStart)*ke
      end function

      !> @brief return the dmplex vertex index for a given 2d vertex index
      function ivertex_icon_2_plex(ivertex, k)
        integer(iintegers),intent(in) :: ivertex  !< @param[in] icell, starts with 1 up to Nfaces (size of icon base grid)
        integer(iintegers),intent(in) :: k        !< @param[in] k, vertical index
        integer(iintegers) :: ivertex_icon_2_plex !< @param[out] icell_icon_2_plex, the vertex index in the dmplex, starts from plex%vStart and goes to plex%vEnd
        integer(iintegers) :: offset
        offset    = Ncells + Nfaces + Nedges
        ivertex_icon_2_plex = offset + k + (ivertex-v2dStart)*ke1
      end function
    end subroutine

    subroutine dump_ownership(dm, cmd_string_dump_ownership, cmd_string_dump_plex)
      type(tDM), intent(in) :: dm
      character(len=*), intent(in) :: cmd_string_dump_ownership
      character(len=*), intent(in), optional :: cmd_string_dump_plex

      type(tPetscSF) :: sf
      integer(iintegers) :: nroots, nleaves
      integer(iintegers), pointer :: pmyidx(:) ! list of my indices that we do not own
      type(PetscSFNode), pointer  :: premote(:) ! rank and remote idx of those points

      type(tDM) :: owner_dm
      type(tPetscSection) :: sec
      type(tVec) :: gVec
      real(ireals), pointer :: xv(:)
      integer(iintegers), pointer :: faces_of_cell(:)

      integer(iintegers) :: cStart, cEnd
      integer(iintegers) :: i, icell, idx, voff
      integer(iintegers), allocatable :: myidx(:)
      type(PetscSFNode), allocatable :: remote(:)

      integer(mpiint) :: comm, myid, ierr

      call PetscObjectGetComm(dm, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call DMClone(dm, owner_dm, ierr); call CHKERR(ierr)
      if(present(cmd_string_dump_plex)) then
        call PetscObjectViewFromOptions(owner_dm, PETSC_NULL_DM, cmd_string_dump_plex, ierr); call CHKERR(ierr)
      endif

      call DMGetPointSF(owner_dm, sf, ierr); call CHKERR(ierr)
      call PetscSFGetGraph(sf, nroots, nleaves, pmyidx, premote, ierr); call CHKERR(ierr)
      call CHKERR(int(nleaves-size(pmyidx), mpiint), 'wrong size of array')
      call CHKERR(int(nleaves-size(premote), mpiint), 'wrong size of array')
      allocate(myidx(nleaves), source=pmyidx)
      allocate(remote(nleaves), source=premote)
      pmyidx => NULL()
      premote => NULL()
      call PetscSortIntWithArrayPair(nleaves, myidx, remote(:)%rank, remote(:)%index, ierr); call CHKERR(ierr)

      call create_plex_section(owner_dm, 'dmplex_ownership info', i1, &
        [i1], [i0], [i0], [i0], sec)
      call DMSetSection(owner_dm, sec, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(sec, ierr); call CHKERR(ierr)

      call DMGetSection(owner_dm, sec, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(owner_dm, gVec, ierr); call CHKERR(ierr)

      ! Dump Cell Ownership
      call PetscObjectSetName(gVec, 'ownership_cells', ierr);call CHKERR(ierr)
      call VecGetArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      xv(:) = myid*one
      call VecRestoreArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(gVec, PETSC_NULL_VEC, cmd_string_dump_ownership, ierr); call CHKERR(ierr)

      ! Dump Face Ownership
      call PetscObjectSetName(gVec, 'ownership_non_local_faces', ierr);call CHKERR(ierr)
      call DMPlexGetHeightStratum(owner_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      call VecGetArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      do icell = cStart, cEnd-1
        call PetscSectionGetOffset(sec, icell, voff, ierr); call CHKERR(ierr)
        xv(i1+voff) = zero
        call DMPlexGetCone(owner_dm, icell, faces_of_cell, ierr); call CHKERR(ierr) ! Get Faces of cell
        do i=1,size(faces_of_cell)
          call PetscFindInt(faces_of_cell(i), nleaves, myidx, idx, ierr); call CHKERR(ierr)
          if(idx.ge.i0) then ! not a local element
            xv(i1+voff) = xv(i1+voff) + i1
          endif
        enddo
      enddo
      call VecRestoreArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(gVec, PETSC_NULL_VEC, cmd_string_dump_ownership, ierr); call CHKERR(ierr)

      ! Dump Sum of cell Ownership
      call PetscObjectSetName(gVec, 'ownership_sum_cells', ierr);call CHKERR(ierr)
      call DMPlexGetHeightStratum(owner_dm, i0, cStart, cEnd, ierr); call CHKERR(ierr) ! cells
      call VecGetArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      xv(:) = real(cEnd-cStart-i1, ireals)
      call VecRestoreArrayF90(gVec, xv, ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(gVec, PETSC_NULL_VEC, cmd_string_dump_ownership, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(owner_dm, gVec, ierr); call CHKERR(ierr)
      call DMDestroy(owner_dm, ierr); call CHKERR(ierr)
    end subroutine

    ! Create a 2D Torus grid with Nx vertices horizontally and Ny rows of Vertices vertically
    subroutine create_2d_fish_plex(Nx, Ny, dm, dmserial)
      type(tDM), intent(inout) :: dm
      type(tDM), intent(inout), optional :: dmserial
      type(tDM) :: dmdist
      type(tPetscSF) :: psf
      integer(iintegers), intent(in) :: Nx, Ny

      integer(iintegers) :: chartsize, Nfaces, Nedges, Nvertices
      integer(iintegers), allocatable :: cellslist(:)
      double precision, allocatable :: vertexCoords(:)

      integer(iintegers) :: pStart, pEnd
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: eStart, eEnd
      integer(iintegers) :: vStart, vEnd

      integer(mpiint) :: myid, numnodes, ierr

      logical, parameter :: lfromcelllist=.False.

      call mpi_comm_rank(PETSC_COMM_WORLD, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(PETSC_COMM_WORLD, numnodes, ierr); call CHKERR(ierr)

      if(modulo(Nx,i2).ne.0) call CHKERR(1_mpiint, 'Nx has to be even, e.g. 2,4,6...')
      if(modulo(Ny,i2).eq.0) call CHKERR(1_mpiint, 'Ny has to be uneven, e.g. 3,5...')

      if(myid.eq.0) then
          Nfaces = Nx * (Ny-1)                ! Number of faces
          Nvertices = (Nx/2+1) * Ny           ! Number of vertices
          Nedges = Nvertices + Nfaces - 1
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

      if(lfromcelllist) then
        call gen_cellslist_array(Nx, Nfaces, Nvertices, cellslist, vertexCoords)
        !print *,Nfaces, Nvertices, 'cells',cellslist,'vcord',vertexCoords
        call DMPlexCreateFromCellList(PETSC_COMM_WORLD, i2, Nfaces, Nvertices, &
                          i3, PETSC_TRUE, cellslist, i3, vertexCoords, dm, ierr); call CHKERR(ierr)
      else
        call DMPlexCreate(PETSC_COMM_WORLD, dm, ierr); call CHKERR(ierr)
        call PetscObjectSetName(dm, 'testplex Nx'//itoa(Nx)//'_Ny'//itoa(Ny), ierr); call CHKERR(ierr)
        call DMSetDimension(dm, i2, ierr); call CHKERR(ierr)

        chartsize = Nfaces + Nedges + Nvertices
        call DMPlexSetChart(dm, i0, chartsize, ierr); call CHKERR(ierr)

        call set_wedge_connectivity(dm, Nx, Nfaces, Nedges)

        call DMPlexSymmetrize(dm, ierr); CHKERRQ(ierr)
        call DMPlexStratify(dm, ierr); CHKERRQ(ierr)

        call set_coords_serial(dm, Nx)
      endif
      call PetscObjectViewFromOptions(dm, PETSC_NULL_DM, "-show_2d_fish", ierr); call CHKERR(ierr)

      call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
      call DMPlexGetHeightStratum(dm, i0, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
      call DMPlexGetHeightStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
      call DMPlexGetHeightStratum(dm, i2, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

      if(ldebug) then
        print *,'pStart,End serial :: ',pStart, pEnd
        print *,'fStart,End serial :: ',fStart, fEnd
        print *,'eStart,End serial :: ',eStart, eEnd
        print *,'vStart,End serial :: ',vStart, vEnd
      endif

      call DMPlexSetAdjacencyUseCone(dm, PETSC_TRUE, ierr); call CHKERR(ierr)
      call DMPlexSetAdjacencyUseClosure(dm, PETSC_FALSE, ierr); call CHKERR(ierr)

      if(present(dmserial)) then
        call DMClone(dm, dmserial, ierr); call CHKERR(ierr)
      endif

      psf = PETSC_NULL_SF
      call DMPlexDistribute(dm, i0, psf, dmdist, ierr); call CHKERR(ierr)
      if(dmdist.ne.PETSC_NULL_DM) then
        call DMDestroy(dm, ierr); call CHKERR(ierr)
        dm = dmdist

        call DMPlexGetChart(dm, pStart, pEnd, ierr); call CHKERR(ierr)
        call DMPlexGetHeightStratum(dm, i0, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
        call DMPlexGetHeightStratum(dm, i1, eStart, eEnd, ierr); call CHKERR(ierr) ! edges
        call DMPlexGetHeightStratum(dm, i2, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices

        if(ldebug) then
          print *,myid,'pStart,End distributed:: ',pStart, pEnd
          print *,myid,'fStart,End distributed:: ',fStart, fEnd
          print *,myid,'eStart,End distributed:: ',eStart, eEnd
          print *,myid,'vStart,End distributed:: ',vStart, vEnd
        endif
      endif

      contains
        subroutine gen_cellslist_array(Nx, Nfaces, Nverts, cellslist, vertexCoords)
          integer(iintegers), intent(in) :: Nx, Nfaces, Nverts
          integer(iintegers), allocatable, intent(out) :: cellslist(:)
          double precision, allocatable, intent(out) :: vertexCoords(:)

          integer(iintegers) :: i, j, k
          integer(iintegers) :: v1, v2, v3
          integer(iintegers) :: vert_per_row

          real(ireals), parameter :: dx=1, dy=1, dz=100
          real(ireals), parameter :: ds=sqrt(dy**2 - (dx/2)**2)

          vert_per_row = Nx/2 + 1

          allocate(cellslist(Nfaces*i3))
          allocate(vertexCoords(Nverts*i3))
          do k = 0, Nfaces-1
            j = k / Nx ! row of faces
            i = k - j*Nx ! col of faces
            if(modulo(i+modulo(j,i2),i2).eq.0) then ! this has a bot edge
              v1 = j*vert_per_row + i/2
              v2 = j*vert_per_row + i/2 + i1
              v3 = (j+1)*vert_per_row + i/2 + modulo(j,i2)
              cellslist(i1+k*i3:(k+i1)*i3) = [v1, v2, v3]
              if(ldebug) print *,k, ':',i, j, 'has bot edge', cellslist(k*i3+i1:(k+i1)*i3)
            else
              v1 = (j+1)*vert_per_row + i/2
              v2 = (j+1)*vert_per_row + i/2 + i1
              v3 = (j)*vert_per_row + i/2 + i1 - modulo(j,i2)
              cellslist(k*i3+i1:(k+i1)*i3) = [v1, v2, v3]
              if(ldebug) print *,k, ':',i, j, 'has top edge', cellslist(k*i3+i1:(k+i1)*i3)
            endif
          enddo

          do k=0,Nverts-1
            j = k / Nx ! row of verts
            i = k - j*Nx ! col of verts
            vertexCoords(i1+k*i3:(k+i1)*i3) = [i*dx + real(modulo(j,i2),ireals)*dx/2, j*ds, dz]
            if(ldebug) print *,k,':',i,j,'vcoord', vertexCoords(i1+k*i3:(k+i1)*i3)
          enddo
          if(ldebug) then
            call mpi_barrier(PETSC_COMM_WORLD, ierr); call CHKERR(ierr)
          endif
        end subroutine

        subroutine set_wedge_connectivity(dm, Nx, Nfaces, Nedges)
          type(tDM) :: dm
          integer(iintegers) :: Nx, Nfaces, Nedges

          integer(iintegers) :: k, i, j, l, ioff, cone3(3), cone2(2)
          integer(iintegers) :: vert_per_row, edge_per_row, base_edges_per_row

          vert_per_row = Nx/2 + 1
          base_edges_per_row = Nx/2
          edge_per_row = base_edges_per_row + (Nx + 1)

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
          do k = 0, Nfaces-1
            j = k / Nx ! row of faces
            i = k - j*Nx ! col of faces
            !print *,'faces',k,': i,j', i, j

            ! determine edges of a face
            if(modulo(i+modulo(j,i2),i2).eq.0) then ! this has a bot edge
                cone3(1) = Nfaces + j*edge_per_row + i/2           ! Nfaces offset + number of edges of full height + number of edges of half heights + i offset
                cone3(2) = Nfaces + j*edge_per_row + base_edges_per_row + i      ! left edge  ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
                cone3(3) = Nfaces + j*edge_per_row + base_edges_per_row + i +1   ! right edge ! Nfaces offset + number of edges of full height + number of edges of half heights + Nedges full heigths on this row + i offset
              !if(ldebug) print *,'upward edge of face', ioff, ':', cone3
            else
                cone3(1) = Nfaces + (j+1)*edge_per_row + i/2
                cone3(2) = Nfaces + j*edge_per_row + base_edges_per_row + i
                cone3(3) = Nfaces + j*edge_per_row + base_edges_per_row + i +1
              !if(ldebug) print *,'downward edge of face', ioff, ':', i, j,':', cone3
            endif

            call DMPlexSetCone(dm,  ioff, cone3, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
          ! Edges have 2 vertices
          do k = 0, Nedges-1
              j = k / edge_per_row   ! row of edge ! goes up to Ny
              i = k - j*edge_per_row ! col of edge
              if(i.lt.Nx/2) then ! bottom edge
                cone2(1) = Nfaces + Nedges + j*vert_per_row + i
                cone2(2) = Nfaces + Nedges + j*vert_per_row + i + i1
                !if(ldebug) print *,k,'- edge',ioff,':',i,j,'            cone',cone2
              else ! sideward edge
                l = i - Nx/2 ! index for vertical edges, i.e. edges that bridge rows
                if(modulo(l+j,i2).eq.0) then ! slash
                  cone2(1) = Nfaces + Nedges + j*vert_per_row + l/2
                  cone2(2) = Nfaces + Nedges + (j+1)*vert_per_row + (l+1)/2
                  !if(ldebug) print *,k,'s edge',ioff,':',i,j,l,'cone',cone2
                else ! backslash
                  cone2(1) = Nfaces + Nedges + j*vert_per_row + (l+1)/2
                  cone2(2) = Nfaces + Nedges + (j+1)*vert_per_row + l/2
                  !if(ldebug) print *,k,'b edge',ioff,':',i,j,l,'cone',cone2
                endif
              endif
            call DMPlexSetCone(dm,  ioff, cone2, ierr); call CHKERR(ierr)
            ioff = ioff + 1
          enddo
        end subroutine

        subroutine set_coords_serial(dm, Nx)
          type(tDM) :: dm
          integer(iintegers), intent(in) :: Nx
          real(ireals), pointer:: coords(:)
          type(tVec)           :: coordinates
          integer(iintegers)   :: dimEmbed, coordSize, vStart, vEnd, pStart, pEnd, eStart, eEnd, voff, cStart, cEnd
          type(tPetscSection)  :: coordSection
          integer(iintegers)   :: iv, i, j

          real(ireals), parameter :: dx=1, dy=1
          real(ireals), parameter :: ds=sqrt(dy**2 - (dx/2)**2)
          real(ireals) :: x, y, z

          integer(iintegers) :: vert_per_row
          vert_per_row = Nx/2 + 1

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
              j = (iv-vStart) / vert_per_row
              i = (iv-vStart) - j*vert_per_row
              z = 100
            x = (real(modulo(j,i2),ireals)*.5_ireals + real(i, ireals))*dx
            y = j*ds
            !if(ldebug) print *,'iv',iv,':', i, j,'=>', x, y, z

            call PetscSectionGetOffset(coordSection, iv, voff, ierr); coords(voff+1:voff+3) = [x, y, z]
          enddo

          call VecRestoreArrayF90(coordinates, coords, ierr); CHKERRQ(ierr)
          call DMSetCoordinatesLocal(dm, coordinates, ierr);CHKERRQ(ierr)
          call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); CHKERRQ(ierr)
          call VecDestroy(coordinates, ierr);CHKERRQ(ierr)
        end subroutine
      end subroutine

    subroutine gen_2d_plex_from_icongridfile(comm, gridfile, dm, dmdist, migration_sf, cell_ao_2d)
      integer(mpiint), intent(in) :: comm
      character(len=*), intent(in) :: gridfile
      type(tDM), intent(out) :: dm, dmdist
      type(tPetscSF), intent(out) :: migration_sf
      AO, allocatable, intent(out), optional :: cell_ao_2d

      integer(iintegers) :: Nfaces, Nedges, Nverts ! number of entries in base icon grid
      integer(iintegers), allocatable :: cell_index(:)
      integer(iintegers), allocatable :: edge_index(:)
      integer(iintegers), allocatable :: vert_index(:)
      real(ireals)      , allocatable :: cell_elevation(:)       ! dim(local_cells)
      real(ireals)      , allocatable :: cartesian_x_vertices(:) ! dim(local_cells)
      real(ireals)      , allocatable :: cartesian_y_vertices(:) ! dim(local_cells)
      real(ireals)      , allocatable :: cartesian_z_vertices(:) ! dim(local_cells)

      integer(iintegers), allocatable, dimension(:,:) :: edges_of_cell ! edges of each cell, dim(local_cells, 3)
      integer(iintegers), allocatable, dimension(:,:) :: edge_verts    ! vertices at edge  , dim(local_edges, 2)

      integer(iintegers), allocatable :: dmplex_idx(:)

      character(len=default_str_len) :: varname(2)
      integer(iintegers) :: chartsize, i, k, icell, iedge, edge3(3), vert2(2)
      integer(mpiint) :: myid, ierr

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        if (ldebug) print *,'Reading Icon icongrid File:', trim(gridfile)
        varname(1) = trim(gridfile)

        varname(2) = 'cell_index'           ; call ncload(varname, cell_index   , ierr); call CHKERR(ierr);
        varname(2) = 'edge_index'           ; call ncload(varname, edge_index   , ierr); call CHKERR(ierr);
        varname(2) = 'vertex_index'         ; call ncload(varname, vert_index   , ierr); call CHKERR(ierr);
        varname(2) = 'edge_of_cell'         ; call ncload(varname, edges_of_cell, ierr); call CHKERR(ierr)
        varname(2) = 'edge_vertices'        ; call ncload(varname, edge_verts   , ierr); call CHKERR(ierr)
        varname(2) = 'cell_elevation'       ; call ncload(varname, cell_elevation,ierr) ; call CHKERR(ierr);
        varname(2) = 'cartesian_x_vertices' ; call ncload(varname, cartesian_x_vertices, ierr) ; call CHKERR(ierr);
        varname(2) = 'cartesian_y_vertices' ; call ncload(varname, cartesian_y_vertices, ierr) ; call CHKERR(ierr);
        varname(2) = 'cartesian_z_vertices' ; call ncload(varname, cartesian_z_vertices, ierr) ; call CHKERR(ierr);

        Nfaces = size(cell_index)
        Nedges = size(edge_index)
        Nverts = size(vert_index)
      else
        Nfaces = i0
        Nedges = i0
        Nverts = i0
      endif

      call DMPlexCreate(comm, dm, ierr);call CHKERR(ierr)
      call PetscObjectSetName(dm, 'DMPLEX_read_from_'//trim(gridfile), ierr);call CHKERR(ierr)
      call DMSetDimension(dm, i2, ierr);call CHKERR(ierr)

      chartsize = Nfaces + Nedges + Nverts
      call DMPlexSetChart(dm, i0, chartsize, ierr); call CHKERR(ierr)

      ! Preallocation
      ! Every cell has three edges
      k = 0
      do i = 1, Nfaces
        call DMPlexSetConeSize(dm, k, i3, ierr); call CHKERR(ierr)
        k = k+1
      enddo

      ! Edges have 2 vertices
      do i = 1, Nedges
        call DMPlexSetConeSize(dm, k, i2, ierr); call CHKERR(ierr)
        k = k+1
      enddo

      call DMSetUp(dm, ierr); call CHKERR(ierr) ! Allocate space for cones

      ! Setup Connections
      ! First set three edges of cell
      do i = 1, Nfaces
        icell = cell_index(i)
        edge3 = edges_of_cell(icell,:)

        call DMPlexSetCone(dm, icell-i1, Nfaces-i1 + edge3, ierr); call CHKERR(ierr)
      enddo

      !! and then set the two vertices of edge
      do i = 1, Nedges
        iedge = edge_index(i)
        vert2 = edge_verts(iedge,:)

        call DMPlexSetCone(dm, Nfaces-i1 + iedge, Nfaces+Nedges-i1 + vert2, ierr); call CHKERR(ierr)
      enddo

      call DMPlexSymmetrize(dm, ierr); call CHKERR(ierr)
      call DMPlexStratify(dm, ierr); call CHKERR(ierr)

      call set_coords()

      if(present(cell_ao_2d)) then
        allocate(cell_ao_2d)
        if(.not.allocated(cell_index)) allocate(cell_index(Nfaces))
        allocate(dmplex_idx(Nfaces), source=(/ (i, i=0,Nfaces-1) /))
        call AOCreateMapping(comm, Nfaces, cell_index-1, dmplex_idx, cell_ao_2d, ierr); call CHKERR(ierr)
        !call AOView(cell_ao_2d, PETSC_VIEWER_STDOUT_WORLD, ierr)
      endif

      call DMPlexSetAdjacencyUseCone(dm, PETSC_TRUE, ierr); call CHKERR(ierr)
      call DMPlexSetAdjacencyUseClosure(dm, PETSC_FALSE, ierr); call CHKERR(ierr)

      call DMPlexDistribute(dm, i0, migration_sf, dmdist, ierr); call CHKERR(ierr)

      if(dmdist.ne.PETSC_NULL_DM) then
        call PetscObjectViewFromOptions(migration_sf, PETSC_NULL_SF, &
          '-show_migration_sf', ierr); call CHKERR(ierr)
      else
        migration_sf = PETSC_NULL_SF
        call DMClone(dm, dmdist, ierr); call CHKERR(ierr)
      endif
      contains
        subroutine set_coords()
          real(ireals),pointer :: coords(:)
          type(tVec)           :: coordinates
          integer(iintegers)   :: coordSize, voff, ind
          PetscSection         :: coordSection
          integer(iintegers)   :: vStart, vEnd, i, j
          integer(iintegers), allocatable :: faces_around_vertex(:)

          real(ireals) :: cart_coord(3), vert_elevation
          real(ireals), parameter :: sphere_radius = 6371229._ireals

          call DMGetCoordinateSection(dm, coordSection, ierr); call CHKERR(ierr)

          call PetscSectionSetNumFields(coordSection, i1, ierr); call CHKERR(ierr)
          call PetscSectionSetUp(coordSection, ierr); call CHKERR(ierr)
          call PetscSectionSetFieldComponents(coordSection, i0, i3, ierr); call CHKERR(ierr)

          call DMPlexGetDepthStratum (dm, i0, vStart, vEnd, ierr); call CHKERR(ierr) ! vertices
          call PetscSectionSetChart(coordSection, vStart, vEnd, ierr);call CHKERR(ierr)

          do i = vStart, vEnd-i1
            call PetscSectionSetDof(coordSection, i, i3, ierr); call CHKERR(ierr)
            call PetscSectionSetFieldDof(coordSection, i, i0, i3, ierr); call CHKERR(ierr)
          enddo

          call PetscSectionSetUp(coordSection, ierr); call CHKERR(ierr)
          call PetscSectionGetStorageSize(coordSection, coordSize, ierr); call CHKERR(ierr)
          !print *,'Coord Section has size:', coordSize

          call VecCreate(comm, coordinates, ierr); call CHKERR(ierr)
          call VecSetSizes(coordinates, coordSize, PETSC_DETERMINE, ierr);call CHKERR(ierr)
          call VecSetBlockSize(coordinates, i3, ierr);call CHKERR(ierr)
          call VecSetType(coordinates, VECSTANDARD, ierr);call CHKERR(ierr)

          call PetscObjectSetName(coordinates, "coordinates", ierr); call CHKERR(ierr)

          call VecGetArrayF90(coordinates, coords, ierr); call CHKERR(ierr)
          !print *,'bounds coords:', lbound(coords), ubound(coords)

          ! set vertices as coordinates
          do i = i1, Nverts
            ind = vStart + i - i1
            call PetscSectionGetOffset(coordSection, ind, voff, ierr); call CHKERR(ierr)

            call get_horizontal_faces_around_vertex(dm, ind, faces_around_vertex)
            vert_elevation = zero
            do j = 1, size(faces_around_vertex)
              vert_elevation = vert_elevation + cell_elevation(i1+faces_around_vertex(j))
            enddo
            vert_elevation = vert_elevation / size(faces_around_vertex)

            cart_coord = [cartesian_x_vertices(i), cartesian_y_vertices(i), cartesian_z_vertices(i)]
            !cart_coord = [icongrid%cartesian_x_vertices(i), icongrid%cartesian_y_vertices(i)]
            coords(voff+i1 : voff+i3) = cart_coord * (sphere_radius + vert_elevation)
            !print *,'setting coords',cart_coord,'to',voff+i1 , voff+dimEmbed
          enddo

          !print *,'coords', shape(coords), '::', coords
          call VecRestoreArrayF90(coordinates, coords, ierr); call CHKERR(ierr)

          call DMSetCoordinatesLocal(dm, coordinates, ierr);call CHKERR(ierr)
          call PetscObjectViewFromOptions(coordinates, PETSC_NULL_VEC, "-show_plex_coordinates", ierr); call CHKERR(ierr)
          call VecDestroy(coordinates, ierr);call CHKERR(ierr)
        end subroutine
    end subroutine

    subroutine icon_ncvec_to_plex(dm2d_serial, dm2d_parallel, migration_sf, &
        filename, varname, gVec, timeidx, dm3d)
      type(tDM), intent(in) :: dm2d_serial
      type(tDM), intent(inout) :: dm2d_parallel
      type(tPetscSF), intent(in) :: migration_sf
      character(len=*), intent(in) :: filename, varname
      type(tVec), allocatable, intent(inout) :: gvec
      integer(iintegers), intent(in), optional :: timeidx
      type(tDM), intent(in), optional :: dm3d

      type(tDM) :: dm3dcopy
      type(tVec) :: rank0vec, newVec
      type(tPetscSection) :: rank0section, newSection, dm3dSection

      real(ireals), allocatable :: arr(:,:)
      real(ireals), allocatable :: arr3d(:,:,:)
      character(len=default_str_len) :: ncgroups(2)

      real(ireals), pointer :: xloc(:)

      integer(iintegers) :: i, k, voff, ke

      integer(mpiint) :: comm, myid, ierr

      if(allocated(gVec)) call CHKERR(1_mpiint, 'gVec already allocated')

      call PetscObjectGetComm(dm2d_parallel, comm, ierr); call CHKERR(ierr)
      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        ncgroups(1) = trim(filename)
        ncgroups(2) = trim(varname)
        if(ldebug) print *,'plex_grid::ncvar2d_to_globalvec : Loading file ', trim(ncgroups(1)), ':', trim(ncgroups(2))
        call ncload(ncgroups, arr, ierr)
        if(ierr.ne.0) then !try to load 3D data
          call ncload(ncgroups, arr3d, ierr); call CHKERR(ierr, 'Could not load Data from NetCDF')

          k = get_arg(i1, timeidx)
          if(k.lt.i1 .or. k.gt.ubound(arr3d,3)) &
            call CHKERR(1_mpiint, 'invalid time index, shape of arr:' &
            //itoa(size(arr3d,1))//itoa(size(arr3d,2))//itoa(size(arr3d,3)))

          allocate(arr(size(arr3d,1), size(arr3d,2)), source=arr3d(:,:,k))
          deallocate(arr3d)
        endif
        ke = size(arr, dim=2)
      endif ! rank 0
      call imp_bcast(comm, ke, 0_mpiint); call CHKERR(ierr)

      call create_plex_section(dm2d_serial, 'face_section', i1, [i0], [ke], [i0], [i0], rank0section)
      call DMSetSection(dm2d_serial, rank0section, ierr); call CHKERR(ierr)

      call DMGetGlobalVector(dm2d_serial, rank0vec, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        call VecGetArrayF90(rank0Vec, xloc,ierr); call CHKERR(ierr)
        do i = 0, size(arr, dim=1)-1
          call PetscSectionGetOffset(rank0Section, i, voff, ierr); call CHKERR(ierr)
          do k = 0, size(arr, dim=2)-1
            xloc(i1+voff+k) = arr(i1+i, ke-k)
          enddo
        enddo
        call VecRestoreArrayF90(rank0Vec, xloc,ierr) ;call CHKERR(ierr)
      endif ! rank 0

      call PetscSectionCreate(comm, newSection, ierr); call CHKERR(ierr)
      call VecCreate(PETSC_COMM_SELF, newVec, ierr); call CHKERR(ierr)

      if(migration_sf.ne.PETSC_NULL_SF) then
        call DMPlexDistributeField(dm2d_serial, migration_sf, rank0section, rank0vec, newsection, newVec, ierr); call CHKERR(ierr)
      else
        newsection = rank0section
        call PetscSectionClone(rank0section, newsection, ierr); call CHKERR(ierr)
        call VecDuplicate(rank0vec, newVec, ierr); call CHKERR(ierr)
        call VecCopy(rank0vec, newVec, ierr); call CHKERR(ierr)
      endif

      call DMRestoreGlobalVector(dm2d_serial, rank0vec, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(rank0section, ierr); call CHKERR(ierr)

      allocate(gVec)
      if(present(dm3d)) then
        call DMClone(dm3d, dm3dcopy, ierr); call CHKERR(ierr)
        call create_plex_section(dm3dcopy, 'cell_section', i1, [i1], [i0], [i0], [i0], dm3dSection)
        call DMSetSection(dm3dcopy, dm3dSection, ierr); call CHKERR(ierr)
        call DMCreateGlobalVector(dm3dcopy, gVec, ierr); call CHKERR(ierr)
        call PetscObjectSetName(gvec, 'VecfromNC_'//trim(varname), ierr);call CHKERR(ierr)
        call VecCopy(newVec, gVec, ierr); call CHKERR(ierr)
        call PetscSectionDestroy(dm3dSection, ierr); call CHKERR(ierr)
        call DMDestroy(dm3dcopy, ierr); call CHKERR(ierr)
      else
        call DMSetSection(dm2d_parallel, newsection, ierr); call CHKERR(ierr)
        call PetscObjectViewFromOptions(newsection, PETSC_NULL_SECTION, "-show_icon_ncvec_to_plex_section", ierr); call CHKERR(ierr)
        call VecDuplicate(newVec, gVec, ierr); call CHKERR(ierr)
        call PetscObjectSetName(gvec, 'VecfromNC_'//trim(varname), ierr);call CHKERR(ierr)
        call VecCopy(newVec, gVec, ierr); call CHKERR(ierr)
      endif
      call VecDestroy(newVec, ierr); call CHKERR(ierr)
      call PetscSectionDestroy(newsection, ierr); call CHKERR(ierr)
    end subroutine

    subroutine Nz_Ncol_vec_to_celldm1(plex, inp, vec)
      type(t_plexgrid), intent(in)  :: plex
      real(ireals), intent(in) :: inp(:,:)
      type(tVec), intent(inout) :: vec

      type(tPetscSection) :: sec
      real(ireals), pointer :: xv(:)
      type(tIS) :: toa_ids
      integer(iintegers), pointer :: xitoa_faces(:), cell_support(:)
      integer(iintegers) :: i, iface, Ncol, ke, vecsize, voff
      integer(mpiint) :: ierr

      if(.not.allocated(plex%cell1_dm)) call CHKERR(1_mpiint, 'plex%cell1_dm has to be allocated')

      call DMGetStratumIS(plex%cell1_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
      call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)
      call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
      ke = vecsize/Ncol

      call CHKERR(int(ke - (plex%Nlay), mpiint), &
        'vertical vec sizes do not match '//itoa(ke)//' vs '//itoa(plex%Nlay)//' => '//itoa(ke - (plex%Nlay)))

      call DMGetSection(plex%cell1_dm, sec, ierr); call CHKERR(ierr)

      call VecGetArrayF90(vec, xv, ierr); call CHKERR(ierr)

      call ISGetIndicesF90(toa_ids, xitoa_faces, ierr); call CHKERR(ierr)
      do i = 1, size(xitoa_faces)
        iface = xitoa_faces(i)
        call DMPlexGetSupport(plex%cell1_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
        call PetscSectionGetOffset(sec, cell_support(1), voff, ierr); call CHKERR(ierr)
        xv(i1+voff: voff+ke) = inp(:,i)
        call DMPlexRestoreSupport(plex%cell1_dm, iface, cell_support, ierr); call CHKERR(ierr) ! support of face is cell
      enddo
      call ISRestoreIndicesF90(toa_ids, xitoa_faces, ierr); call CHKERR(ierr)

      call VecRestoreArrayF90(vec, xv, ierr); call CHKERR(ierr)
    end subroutine

    subroutine celldm1_vec_to_Nz_Ncol(plex, vec, arr)
      type(t_plexgrid), intent(in)  :: plex
      type(tVec), intent(in) :: vec
      real(ireals), allocatable, intent(inout) :: arr(:,:)

      type(tPetscSection) :: sec
      real(ireals), pointer :: xv(:)
      type(tIS) :: toa_ids
      integer(iintegers), pointer :: xitoa_faces(:)
      integer(iintegers) :: i, iface, Ncol, ke, vecsize, voff
      integer(mpiint) :: ierr

      call DMGetStratumIS(plex%cell1_dm, 'DomainBoundary', TOAFACE, toa_ids, ierr); call CHKERR(ierr)
      call ISGetSize(toa_ids, Ncol, ierr); call CHKERR(ierr)
      call VecGetLocalSize(vec, vecsize, ierr); call CHKERR(ierr)
      ke = vecsize/Ncol

      call CHKERR(int(ke - plex%Nlay, mpiint), &
        'vertical vec sizes do not match '//itoa(ke)//' vs '//itoa(plex%Nlay))

      if(.not.allocated(arr)) then

        allocate(arr(ke, Ncol))
      else
        ke = min(vecsize/Ncol, size(arr,dim=1,kind=iintegers))
        if(.not.all(shape(arr).eq.[ke, Ncol])) then
          print *,'shape cell1dm', [ke, Ncol], 'shape out_arr', shape(arr)
          call CHKERR(1_mpiint, 'shape of out_arr does not conform to cell1dm sizes')
        endif
      endif

      call DMGetSection(plex%cell1_dm, sec, ierr); call CHKERR(ierr)

      call VecGetArrayReadF90(vec, xv, ierr); call CHKERR(ierr)

      call ISGetIndicesF90(toa_ids, xitoa_faces, ierr); call CHKERR(ierr)
      do i = 1, size(xitoa_faces)
        iface = xitoa_faces(i)
        call PetscSectionGetOffset(sec, iface, voff, ierr); call CHKERR(ierr)
        arr(i1:ke, i) = xv(i1+voff: voff+ke)
      enddo
      call ISRestoreIndicesF90(toa_ids, xitoa_faces, ierr); call CHKERR(ierr)

      call VecRestoreArrayReadF90(vec, xv, ierr); call CHKERR(ierr)
    end subroutine

    subroutine dm2d_vec_to_Nz_Ncol(dm, vec, arr)
      type(tDM), intent(in)  :: dm
      type(tVec), intent(in) :: vec
      real(ireals), allocatable, intent(inout) :: arr(:,:)

      type(tPetscSection) :: sec
      real(ireals), pointer :: xv(:)
      integer(iintegers) :: fStart, fEnd
      integer(iintegers) :: iface, Ncol, ke, voff
      integer(mpiint) :: ierr

      call DMPlexGetHeightStratum(dm, i0, fStart, fEnd, ierr); call CHKERR(ierr) ! faces
      call DMGetSection(dm, sec, ierr); call CHKERR(ierr)
      call PetscSectionGetDof(sec, fStart, ke, ierr); call CHKERR(ierr)
      Ncol = fEnd - fStart

      if(.not.allocated(arr)) then
        allocate(arr(ke, Ncol))
      else
        ke = min(ke, size(arr, dim=1, kind=iintegers))
        if(.not.all(shape(arr).eq.[ke, Ncol])) then
          print *,'shape cell1dm', [ke, Ncol], 'shape out_arr', shape(arr)
          call CHKERR(1_mpiint, 'shape of out_arr does not conform to cell1dm sizes')
        endif
      endif

      call VecGetArrayReadF90(vec, xv, ierr); call CHKERR(ierr)

      do iface = fStart, fEnd-1
        call PetscSectionGetOffset(sec, iface, voff, ierr); call CHKERR(ierr)
        arr(i1:ke, i1+iface-fStart) = xv(i1+voff: voff+ke)
      enddo

      call VecRestoreArrayReadF90(vec, xv, ierr); call CHKERR(ierr)
    end subroutine

    subroutine celldm_veccopy(vin, vout)
      type(tVec), intent(in) :: vin
      type(tVec), intent(inout) :: vout
      integer(iintegers) :: vsize1, vsize2
      integer(mpiint) :: ierr
      if(ldebug) then
        call VecGetLocalSize(vin , vsize1, ierr); call CHKERR(ierr)
        call VecGetLocalSize(vout, vsize2, ierr); call CHKERR(ierr)
        call CHKERR(int(vsize1-vsize2, mpiint), 'VecSizes of input and output vecs dont match!')
      endif
      call VecCopy(vin, vout, ierr); call CHKERR(ierr)
    end subroutine


    !Convert a date to Julian Day.
    !
    ! Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet',
    !    4th ed., Duffet-Smith and Zwart, 2011.
    !
    !Parameters
    !----------
    !year : int
    !    Year as integer. Years preceding 1 A.D. should be 0 or negative.
    !    The year before 1 A.D. is 0, 10 B.C. is year -9.
    !month : int
    !    Month as integer, Jan = 1, Feb. = 2, etc.
    !day : float
    !    Day, may contain fractional part.
    !Returns
    !-------
    !jd : float
    !    Julian Day
    !Examples
    !--------
    !Convert 6 a.m., February 17, 1985 to Julian Day
    !
    !>>> date_to_jd(1985,2,17.25)
    !2446113.75
    function date_to_julian_day(year,month,day) result(jd)
    integer(iintegers), intent(in) :: year, month
    real(ireals), intent(in) :: day
    real(ireals) :: jd

    integer(iintegers) :: yearp, monthp, A, B, C, D
    if(month.eq.1 .or. month.eq.2) then
      yearp = year - 1
      monthp = month + 12
    else
      yearp = year
      monthp = month
    endif

    ! this checks where we are in relation to October 15, 1582, the beginning
    ! of the Gregorian calendar.
    if ((year .lt. 1582) .or. &
      (year .eq. 1582 .and. month .lt. 10) .or. &
      (year .eq. 1582 .and. month .eq. 10 .and. day .lt. 15._ireals)) then
      ! before start of Gregorian calendar
      B = 0
    else
      ! after start of Gregorian calendar
      A = int(yearp / 100._ireals)
      B = 2 - A + int(A / 4._ireals)
    endif

    if(yearp .lt. 0) then
      C = int((365.25_ireals * yearp) - 0.75_ireals)
    else
      C = int(365.25_ireals * yearp)
    endif
    D = int(30.6001_ireals * (monthp + 1))
    jd = B + C + D + day + 1720994.5_ireals
  end function

  ! after the wiki page: https://en.wikipedia.org/wiki/Position_of_the_Sun
  function get_sun_vector(year, month, day) result(sundir)
    integer(iintegers), intent(in) :: year, month
    real(ireals), intent(in) :: day
    real(ireals) :: sundir(3)

    real(ireals) :: n, L, g, lambda, delta, eps

    n = date_to_julian_day(year,month,day) - 2451545.0_ireals
    L = 280.460_ireals + 0.9856474_ireals * n
    g = 357.528_ireals + 0.9856003_ireals * n

    lambda = L + 1.915_ireals * sin(deg2rad(g)) + 0.020_ireals * sin(deg2rad(2*g))

    eps = 23.439_ireals - 4e-7_ireals * n

    delta = asin(sin(deg2rad(eps)) * sin(deg2rad(lambda)))

    sundir = [ -cos(2*pi*day) * cos(delta), &
                sin(2*pi*day) * cos(delta), &
                sin(delta) ]
  end function

end module
