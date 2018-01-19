module m_icon_grid

#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload

  use m_data_parameters, only : ireals, iintegers, mpiint, default_str_len

  use m_helper_functions, only: get_arg, imp_bcast, chkerr, unique

  implicit none

  private
  public :: t_icongrid, distribute_icon_grid, read_icon_grid_file, bcast_icongrid, &
    decompose_icon_grid

  type :: t_icongrid
    integer(iintegers) :: Nfaces, Nedges, Nvertices ! number of entries in base icon grid
    integer(iintegers), allocatable, dimension(:) :: parNfaces, parNedges, parNvertices ! number of entries in base icon grid on each processor, dim=(0..numnodes-1)

    ! index of cells in parent grid :: cell_index=dim(Nfaces)
    integer(iintegers), allocatable :: cell_index(:)
    integer(iintegers), allocatable :: edge_index(:)
    integer(iintegers), allocatable :: vertex_index(:)
    real(ireals)      , allocatable :: cartesian_x_vertices(:) ! dim(local_cells)
    real(ireals)      , allocatable :: cartesian_y_vertices(:) ! dim(local_cells)
    real(ireals)      , allocatable :: cartesian_z_vertices(:) ! dim(local_cells)

    integer(iintegers), allocatable, dimension(:,:) :: vertex_of_cell    ! vertices of each cell            , dim(local_cells, 3)
    integer(iintegers), allocatable, dimension(:,:) :: edge_of_cell      ! edges of each cell               , dim(local_cells, 3)
    integer(iintegers), allocatable, dimension(:,:) :: edge_vertices     ! vertices at the end of each edge , dim(local_edges, 2)

    ! Following variables are defined on parent cell/edge/vertex sizes (dimension of parent grid)
    ! i.e. gives local index for e.g. a global cell number
    integer(iintegers), allocatable :: local_cell_index(:)
    integer(iintegers), allocatable :: local_edge_index(:)
    integer(iintegers), allocatable :: local_vertex_index(:)

    ! Remote indices are the local indices on neighbouring ranks (has dimensions of parent grid)
    integer(iintegers), allocatable :: remote_cell_index(:)
    integer(iintegers), allocatable :: remote_edge_index(:)
    integer(iintegers), allocatable :: remote_vertex_index(:)

    integer(iintegers), allocatable, dimension(:,:) :: adj_cell_of_edge  ! cells adjacent to an edge        , dim(parent_edges, 2)
    integer(iintegers), allocatable, dimension(:,:) :: cells_of_vertex   ! cells adjacent to an vertex      , dim(parent_vertices, 6)

    integer(iintegers), allocatable :: cellowner(:)   ! dim=(parent_Nfaces)
    integer(iintegers), allocatable :: edgeowner(:)   ! dim=(parent_Nedges)
    integer(iintegers), allocatable :: vertexowner(:) ! dim=(parent_Nvertices)
  end type

  logical, parameter :: ldebug=.True.
  integer(iintegers), parameter :: ICONULL=-1

  contains
    subroutine distribute_icon_grid(comm, icongrid, local_icongrid)
      MPI_Comm, intent(in) :: comm
      type(t_icongrid), allocatable, intent(in) :: icongrid
      type(t_icongrid), allocatable, intent(out) :: local_icongrid

      integer(mpiint) :: myid, numnodes, ierr

      integer(iintegers) :: i, j, ic, ie, iv, ilocal
      integer(iintegers) :: owner
      logical,allocatable :: ladj_cell(:)
      integer(iintegers),allocatable :: par_cnt(:), iowner(:), unique_owner(:)

      if(.not.allocated(icongrid)) stop 'distribute_icon_grid :: global icongrid not allocated!'
      allocate(local_icongrid)

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      call decompose_icon_grid(icongrid, numnodes, local_icongrid%cellowner, local_icongrid%edgeowner, local_icongrid%vertexowner)

      allocate(ladj_cell(0:numnodes-1))

      allocate(local_icongrid%parNfaces(0:numnodes-1))
      do owner=0,numnodes-1
        local_icongrid%parNfaces(owner) = count(local_icongrid%cellowner.eq.owner)
      enddo
      local_icongrid%Nfaces = local_icongrid%parNfaces(myid)

      ! Count edges/vertices for local grid
      allocate(local_icongrid%parNedges(0:numnodes-1), source=0)
      do i=1,icongrid%Nedges
        ie = icongrid%edge_index(i)
        ladj_cell = .False.
        do j=1,size(icongrid%adj_cell_of_edge(ie,:))
          ic = icongrid%adj_cell_of_edge(ie,j)
          if(ic.gt.0) then
            ladj_cell(local_icongrid%cellowner(ic)) = .True.
          endif
        enddo
        where(ladj_cell)
          local_icongrid%parNedges = local_icongrid%parNedges+1
        endwhere
      enddo
      local_icongrid%Nedges = local_icongrid%parNedges(myid)

      allocate(local_icongrid%parNvertices(0:numnodes-1), source=0)
      do i=1,icongrid%Nvertices
        iv = icongrid%vertex_index(i)
        ladj_cell = .False.
        do j=1,size(icongrid%cells_of_vertex(iv,:))
          ic = icongrid%cells_of_vertex(iv,j)
          if(ic.gt.0) then
            ladj_cell(local_icongrid%cellowner(ic)) = .True.
          endif
        enddo
        where(ladj_cell)
          local_icongrid%parNvertices = local_icongrid%parNvertices+1
        endwhere
      enddo
      local_icongrid%Nvertices = local_icongrid%parNvertices(myid)
      deallocate(ladj_cell)

      if(ldebug) then
        do i=0,numnodes-1
          if(i.eq.myid) then
            print *,myid,'distribute_icon_grid :: locally, I have ', local_icongrid%Nfaces, 'cells'
            print *,myid,'distribute_icon_grid :: locally, I have ', local_icongrid%Nedges, 'edges'
            print *,myid,'distribute_icon_grid :: locally, I have ', local_icongrid%Nvertices, 'vertices'
          endif
          call mpi_barrier(comm, ierr)
        enddo
      endif

      ! Now that we know the sizes of the local domains, lets setup the relationships between local and global grid
      allocate(local_icongrid%cell_index  (local_icongrid%Nfaces))
      allocate(local_icongrid%edge_index  (local_icongrid%Nedges))
      allocate(local_icongrid%vertex_index(local_icongrid%Nvertices))

      allocate(local_icongrid%local_cell_index  (icongrid%Nfaces), source=ICONULL)
      allocate(local_icongrid%local_edge_index  (icongrid%Nedges), source=ICONULL)
      allocate(local_icongrid%local_vertex_index(icongrid%Nvertices), source=ICONULL)

      allocate(local_icongrid%remote_cell_index  (icongrid%Nfaces), source=ICONULL)
      allocate(local_icongrid%remote_edge_index  (icongrid%Nedges), source=ICONULL)
      allocate(local_icongrid%remote_vertex_index(icongrid%Nvertices), source=ICONULL)

      allocate(par_cnt(0:numnodes-1))
      par_cnt = 0
      do i=1,icongrid%Nfaces
        ic = icongrid%cell_index(i)          ! parent cell index
        owner = local_icongrid%cellowner(ic)
        par_cnt(owner) = par_cnt(owner)+1
        if(owner.eq.myid) then
          local_icongrid%cell_index(par_cnt(owner)) = ic
          local_icongrid%local_cell_index(ic) = par_cnt(owner)
        else
          local_icongrid%remote_cell_index(ic) = par_cnt(owner)
        endif
      enddo

      do i=1,local_icongrid%Nfaces
        ic = local_icongrid%cell_index(i)
        ilocal = local_icongrid%local_cell_index(ic)
        print *,myid,i,'global cell', ic, ' to local ->', ilocal, 'remote', local_icongrid%remote_cell_index(ic)
      enddo

      par_cnt = 0
      allocate(iowner(2))
      do i=1,icongrid%Nedges
        ie = icongrid%edge_index(i)
        do j=1,size(icongrid%adj_cell_of_edge(ie,:))
          ic = icongrid%adj_cell_of_edge(ie,j)
          if(ic.gt.0) then
            iowner(j) = local_icongrid%cellowner(ic)
          else
            iowner(j) = ICONULL
          endif
        enddo
        unique_owner = unique(iowner)
        do j=1,size(unique_owner)
          owner = unique_owner(j)
          if(owner.eq.ICONULL) cycle
          par_cnt(owner) = par_cnt(owner) +1
          if(owner.eq.myid) then
            local_icongrid%edge_index(par_cnt(owner)) = ie
            local_icongrid%local_edge_index(ie) = par_cnt(owner)
          else
            local_icongrid%remote_edge_index(ie) = par_cnt(owner)
          endif
        enddo
      enddo
      deallocate(iowner)
      do i=1,local_icongrid%Nedges
        ie = local_icongrid%edge_index(i)
        ilocal = local_icongrid%local_edge_index(ie)
        print *,myid,i,'global edges', ie, ' to local ->', ilocal, 'remote', local_icongrid%remote_edge_index(ie)
      enddo

      allocate(local_icongrid%cartesian_x_vertices(local_icongrid%Nvertices))
      allocate(local_icongrid%cartesian_y_vertices(local_icongrid%Nvertices))
      allocate(local_icongrid%cartesian_z_vertices(local_icongrid%Nvertices))
      allocate(iowner(size(icongrid%cells_of_vertex, dim=2)))
      par_cnt = 0
      do i=1,icongrid%Nvertices
        iv = icongrid%vertex_index(i)
        do j=1,size(icongrid%cells_of_vertex(iv,:))
          ic = icongrid%cells_of_vertex(iv,j)
          if(ic.gt.0) then
            iowner(j) = local_icongrid%cellowner(ic)
          else
            iowner(j) = ICONULL
          endif
        enddo
        unique_owner = unique(iowner)
        do j=1,size(unique_owner)
          owner = unique_owner(j)
          if(owner.eq.ICONULL) cycle
          par_cnt(owner) =  par_cnt(owner) +1
          if(owner.eq.myid) then
            local_icongrid%vertex_index(par_cnt(owner)) = iv
            local_icongrid%local_vertex_index(iv) = par_cnt(owner)
            local_icongrid%cartesian_x_vertices(par_cnt(owner)) = icongrid%cartesian_x_vertices(iv)
            local_icongrid%cartesian_y_vertices(par_cnt(owner)) = icongrid%cartesian_y_vertices(iv)
            local_icongrid%cartesian_z_vertices(par_cnt(owner)) = icongrid%cartesian_z_vertices(iv)
          else
            local_icongrid%remote_vertex_index(iv) = par_cnt(owner)
          endif
        enddo
      enddo
      deallocate(iowner)

      ! Translate auxilliary connections
      allocate(local_icongrid%vertex_of_cell  (local_icongrid%Nfaces, size(icongrid%vertex_of_cell, dim=2)))
      allocate(local_icongrid%edge_of_cell    (local_icongrid%Nfaces, size(icongrid%edge_of_cell  , dim=2)))
      allocate(local_icongrid%edge_vertices   (local_icongrid%Nedges, size(icongrid%edge_vertices, dim=2)))
      allocate(local_icongrid%adj_cell_of_edge(local_icongrid%Nedges, size(icongrid%adj_cell_of_edge, dim=2)))

      do i=1,size(local_icongrid%cell_index)
        ic = local_icongrid%cell_index(i) ! cell index in parent grid
        do j=1,size(icongrid%edge_of_cell(ic,:))
          ie = icongrid%edge_of_cell(ic, j)  ! edge index in parent grid
          ilocal = local_icongrid%local_edge_index(ie) ! edge index locally
          local_icongrid%edge_of_cell(i,j) = ilocal
        enddo
        do j=1,size(icongrid%vertex_of_cell(ic,:))
          iv = icongrid%vertex_of_cell(ic, j)  ! vertex index in parent grid
          local_icongrid%vertex_of_cell(i,j) = local_icongrid%local_vertex_index(iv)
        enddo
      enddo

      do i=1,local_icongrid%Nedges
        ie = local_icongrid%edge_index(i)  ! edge index in parent grid
        do j=1,size(icongrid%edge_vertices(ie,:))
          iv = icongrid%edge_vertices(ie,j)  ! vertex index in parent grid
          local_icongrid%edge_vertices(i,j) = local_icongrid%local_vertex_index(iv)
        enddo
        do j=1,size(icongrid%adj_cell_of_edge(ie,:))
          ic = icongrid%adj_cell_of_edge(ie,j)  ! cell index in parent grid
          if(ic.gt.0) then
            local_icongrid%adj_cell_of_edge(i,j) = local_icongrid%local_cell_index(ic)
          else
            local_icongrid%adj_cell_of_edge(i,j) = ICONULL
          endif
        enddo
      enddo

      if(ldebug) then
        do owner=0,numnodes-1
          if(owner.eq.myid) then
            print *,myid,'distribute_icon_grid :: my cells global indices ', local_icongrid%cell_index
            print *,myid,'distribute_icon_grid :: my edge global indices  ', local_icongrid%edge_index
            print *,myid,'distribute_icon_grid :: my vertex global indices', local_icongrid%vertex_index

            do i=1,size(local_icongrid%cell_index)
              ic = local_icongrid%cell_index(i)       ! global cell index in parent grid
              j  = local_icongrid%local_cell_index(ic) ! local cell index
              print *,myid,'global 2 local cell', ic, '->', j, ' ::: ', i, '<-', ic
            enddo
            do i=1,size(local_icongrid%edge_index)
              ie = local_icongrid%edge_index(i)       ! global edge index in parent grid
              j  = local_icongrid%local_edge_index(ie) ! local edge index
              print *,myid,'global 2 local edge', ie, '->', j, ' ::: ', i, '<-', ie
            enddo

            do i=1,size(local_icongrid%cell_index)
              ic = local_icongrid%cell_index(i)        ! global cell index in parent grid
              j  = local_icongrid%local_cell_index(ic) ! local cell index
              print *,myid,'edge_of_cell', ic, ': local cell', i, '::',icongrid%edge_of_cell(ic,:),'->',local_icongrid%edge_of_cell(j,:)
            enddo
            do i=1,size(local_icongrid%cell_index)
              ic = local_icongrid%cell_index(i)       ! global cell index in parent grid
              j  = local_icongrid%local_cell_index(ic) ! local cell index
              print *,myid,'vertex_of_cell', ic, '::', i,'->',local_icongrid%vertex_of_cell(j,:)
            enddo
          endif
          call mpi_barrier(comm, ierr)
        enddo
        call mpi_barrier(comm, ierr)
      endif
    end subroutine
    subroutine read_icon_grid_file(fname, icongrid)
      character(len=*),intent(in) :: fname
      type(t_icongrid),intent(inout),allocatable :: icongrid
      character(default_str_len) :: varname(2)
      integer(mpiint) :: ierr

      if(allocated(icongrid)) then
        stop 'Icongrid already loaded...'
      endif

      allocate(icongrid)

      if (ldebug) print *,'Reading Icon icongrid File:', trim(fname)
      varname(1) = trim(fname)

      varname(2) = 'vertex_of_cell'       ; call ncload(varname, icongrid%vertex_of_cell  , ierr)     ; call CHKERR(ierr)
      varname(2) = 'edge_of_cell'         ; call ncload(varname, icongrid%edge_of_cell    , ierr)     ; call CHKERR(ierr)
      varname(2) = 'edge_vertices'        ; call ncload(varname, icongrid%edge_vertices   , ierr)     ; call CHKERR(ierr)
      varname(2) = 'adjacent_cell_of_edge'; call ncload(varname, icongrid%adj_cell_of_edge, ierr)     ; call CHKERR(ierr)
      varname(2) = 'cells_of_vertex'      ; call ncload(varname, icongrid%cells_of_vertex , ierr)     ; call CHKERR(ierr)
      varname(2) = 'cell_index'           ; call ncload(varname, icongrid%cell_index          , ierr) ; call CHKERR(ierr);
      varname(2) = 'edge_index'           ; call ncload(varname, icongrid%edge_index          , ierr) ; call CHKERR(ierr);
      varname(2) = 'vertex_index'         ; call ncload(varname, icongrid%vertex_index        , ierr) ; call CHKERR(ierr);
      varname(2) = 'cartesian_x_vertices' ; call ncload(varname, icongrid%cartesian_x_vertices, ierr) ; call CHKERR(ierr);
      varname(2) = 'cartesian_y_vertices' ; call ncload(varname, icongrid%cartesian_y_vertices, ierr) ; call CHKERR(ierr);
      varname(2) = 'cartesian_z_vertices' ; call ncload(varname, icongrid%cartesian_z_vertices, ierr) ; call CHKERR(ierr);

      icongrid%Nfaces = size(icongrid%cell_index)
      icongrid%Nedges = size(icongrid%edge_index)
      icongrid%Nvertices = size(icongrid%vertex_index)

      if (ldebug) then
        print *,'shape vertex of cell  ', shape(icongrid%vertex_of_cell),'::', &
          minval(icongrid%vertex_of_cell), maxval(icongrid%vertex_of_cell)

        print *,'shape edge of cell    ', shape(icongrid%edge_of_cell), '::', &
          minval(icongrid%edge_of_cell), maxval(icongrid%edge_of_cell)

        print *,'shape edge vertices   ', shape(icongrid%edge_vertices), '::', &
          minval(icongrid%edge_vertices), maxval(icongrid%edge_vertices)

        print *,'shape adj_cell_of_edge', shape(icongrid%adj_cell_of_edge), '::', &
          minval(icongrid%adj_cell_of_edge), maxval(icongrid%adj_cell_of_edge)

        print *,'shape cells_of_vertex ', shape(icongrid%cells_of_vertex), '::', &
          minval(icongrid%cells_of_vertex), maxval(icongrid%cells_of_vertex)

        print *,'shape cell_index      ', shape(icongrid%cell_index  )
        print *,'shape edge_index      ', shape(icongrid%edge_index  )
        print *,'shape vertex_index    ', shape(icongrid%vertex_index)
      endif
    end subroutine

    subroutine bcast_icongrid(comm, icongrid)
      MPI_Comm, intent(in) :: comm
      type(t_icongrid), allocatable, intent(inout) :: icongrid

      integer(mpiint) :: myid, numnodes, ierr

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        if(.not.allocated(icongrid)) stop 'bcast_icongrid :: icon grid has to be allocated on rank 0!'
      else
        if(allocated(icongrid)) stop 'bcast_icongrid :: icon grid should not be be allocated on rank > 0!'
        allocate(icongrid)
      endif

      call imp_bcast(comm, icongrid%Nfaces,    0)
      call imp_bcast(comm, icongrid%Nedges,    0)
      call imp_bcast(comm, icongrid%Nvertices, 0)
      call imp_bcast(comm, icongrid%cell_index          , 0)
      call imp_bcast(comm, icongrid%edge_index          , 0)
      call imp_bcast(comm, icongrid%vertex_index        , 0)
      call imp_bcast(comm, icongrid%cartesian_x_vertices, 0)
      call imp_bcast(comm, icongrid%cartesian_y_vertices, 0)
      call imp_bcast(comm, icongrid%cartesian_z_vertices, 0)
      call imp_bcast(comm, icongrid%vertex_of_cell      , 0)
      call imp_bcast(comm, icongrid%edge_of_cell        , 0)
      call imp_bcast(comm, icongrid%edge_vertices       , 0)
      call imp_bcast(comm, icongrid%adj_cell_of_edge    , 0)
      call imp_bcast(comm, icongrid%cells_of_vertex     , 0)
    end subroutine

    subroutine decompose_icon_grid(icongrid, numnodes, cellowner, edgeowner, vertexowner)
      type(t_icongrid),intent(in) :: icongrid
      integer(iintegers), intent(in) :: numnodes

      integer(iintegers),allocatable,intent(out) :: cellowner(:)   ! dim=(icongrid%Nfaces)
      integer(iintegers),allocatable,intent(out) :: edgeowner(:)   ! dim=(icongrid%Nedges)
      integer(iintegers),allocatable,intent(out) :: vertexowner(:) ! dim=(icongrid%Nvertices)

      integer(iintegers) :: cellsofedge(2,icongrid%Nedges),neighborcells(3,icongrid%Nfaces)
      integer(iintegers) :: area(numnodes), ichange
      integer(iintegers) :: i, icell, ie, iedge, iiter
      logical :: lwin, ladjacent_check

      if(.not. allocated(icongrid%cell_index)) stop 'decompose_icon_grid :: the grid is not loaded from an icon gridfile ... cant be decomposed with this method'

      allocate(cellowner(icongrid%Nfaces))
      allocate(edgeowner(icongrid%Nedges))
      allocate(vertexowner(icongrid%Nvertices))

      area = 0
      cellsofedge   = ICONULL
      neighborcells = ICONULL

      do icell=1,icongrid%Nfaces
        !print *,'icell', icell, ':: e', icongrid%edge_of_cell(icell,:)

        do ie=1,3
          iedge=icongrid%edge_of_cell(icell,ie)
          if(cellsofedge(1,iedge).eq.ICONULL) then
            cellsofedge(1,iedge) = icell
          else if(cellsofedge(2,iedge).eq.ICONULL) then
            cellsofedge(2,iedge) = icell
          else
            print *,cellsofedge(:,iedge)
            stop 'should not be here?'
          endif
        enddo
      enddo

      do iedge=1,icongrid%Nedges
        !print *,'edge',iedge,'::',cellsofedge(:,iedge)
        icell = cellsofedge(1,iedge)
        if(neighborcells(1, icell).eq.ICONULL) then
          neighborcells(1, icell) = cellsofedge(2,iedge)
        else if(neighborcells(2, icell).eq.ICONULL) then
          neighborcells(2, icell) = cellsofedge(2,iedge)
        else if(neighborcells(3, icell).eq.ICONULL) then
          neighborcells(3, icell) = cellsofedge(2,iedge)
        else
          print *,'neighbors',icell,'::',neighborcells(:,icell)
          stop 'should not be here?'
        endif
        icell = cellsofedge(2,iedge)
        if(icell.eq.ICONULL) cycle  ! only one edge here

        if(neighborcells(1, icell).eq.ICONULL) then
          neighborcells(1, icell) = cellsofedge(1,iedge)
        else if(neighborcells(2, icell).eq.ICONULL) then
          neighborcells(2, icell) = cellsofedge(1,iedge)
        else if(neighborcells(3, icell).eq.ICONULL) then
          neighborcells(3, icell) = cellsofedge(1,iedge)
        else
          print *,'neighbors',icell,'::',neighborcells(:,icell)
          stop 'should not be here?'
        endif
      enddo

      cellowner = ICONULL
      ! initial conquering of cells
      do i=1,numnodes
        icell = icongrid%Nfaces / numnodes * (i-1) +1
        call conquer_cell(icell, i, lwin)
      enddo

      do iiter=1,icongrid%Nfaces*10
        ichange=0
        do icell=1,icongrid%Nfaces
          if(cellowner(icell).ne.ICONULL) then
            do i=1,3
              if(neighborcells(i, icell).ne.ICONULL) then
                call conquer_cell(neighborcells(i, icell), cellowner(icell), lwin) ! conquer my neighbors
                if(lwin) ichange = ichange+1
                ladjacent_check = check_adjacency()
              endif
            enddo
          endif
        enddo
        ladjacent_check = check_adjacency()
        if(ldebug) print *,'iiter',iiter, ichange, ladjacent_check
        if(ichange.eq.0 .and. ladjacent_check) exit
      enddo

      !do icell=1,icongrid%Nfaces
      !  if(ldebug) print *,'icell',icell,'::', cellowner(icell)
      !enddo
      do i=1,numnodes
        if(ldebug) print *,'node',i, 'area', area(i)
      enddo

      ! convert from 1..numnodes to mpi rank numbering
      cellowner(:) = cellowner(:) -1

      call distribute_edges()
      call distribute_vertices()

      contains
        subroutine distribute_edges()
          integer(iintegers) :: i, ic, ie, new_edge_owner
          do ie=1,icongrid%Nedges
            new_edge_owner = ICONULL
            do i=1,ubound(icongrid%adj_cell_of_edge,2)
              ic = icongrid%adj_cell_of_edge(ie,i)
              if(ic.gt.0) then
                if(cellowner(ic) .ge. new_edge_owner) new_edge_owner=cellowner(ic)
              endif
            enddo
            edgeowner(ie) = new_edge_owner
          enddo
          !if(ldebug) then
          !  do ie=1,icongrid%Nedges
          !    print *,'Edge',ie,'owned by',edgeowner(ie)
          !  enddo
          !endif
        end subroutine

        subroutine distribute_vertices()
          integer(iintegers) :: i, ic, iv, new_vertex_owner
          do iv=1,icongrid%Nvertices
            new_vertex_owner = ICONULL
            do i=1,ubound(icongrid%cells_of_vertex,2)
              ic = icongrid%cells_of_vertex(iv,i)
              if(ic.gt.0) then
                if(cellowner(ic) .ge. new_vertex_owner) new_vertex_owner=cellowner(ic)
              endif
            enddo
            vertexowner(iv) = new_vertex_owner
          enddo
          !if(ldebug) then
          !  do iv=1,icongrid%Nvertices
          !    print *,'Vertex',iv,'owned by',vertexowner(iv)
          !  enddo
          !endif
        end subroutine

        logical function check_adjacency()
          integer(iintegers) :: owners(3), ineigh, ineighcell, icell, iowner
          logical :: lwin

          check_adjacency=.True.

          do icell=1,icongrid%Nfaces
            owners=ICONULL
            do ineigh=1,3
              ineighcell = neighborcells(ineigh, icell)
              if(ineighcell.ne.ICONULL) owners(ineigh) = cellowner(ineighcell)
            enddo
            if(.not.any(owners.eq.cellowner(icell))) then
              if(ldebug) print *,'adjacency not fullfilled! :( --',icell,'owned by',cellowner(icell),'neighbours',neighborcells(:, icell)
              check_adjacency=.False.
              do iowner=1,3
                if(owners(iowner).ne.ICONULL) call conquer_cell(icell, owners(iowner), lwin, lforce=.True.)
                exit
              enddo
            endif
          enddo
        end

        subroutine conquer_cell(icell, conqueror, lwin, lforce)
          integer(iintegers),intent(in) :: icell, conqueror
          integer(iintegers) :: old_owner, ineigh
          real(ireals) :: wgt(2)
          logical,intent(out) :: lwin
          logical,intent(in),optional :: lforce
          integer(iintegers) :: ineigh_cell, neigh_owners(3)
          integer(iintegers) :: owner_cellcnt(numnodes)
          logical :: ladj

          lwin=.False.
          wgt = 0

          old_owner = cellowner(icell)
          if( old_owner.eq.ICONULL )then
            lwin = .True.
          else
            if(old_owner.eq.conqueror) return
            wgt(1) = (real(area(old_owner)) / real(area(conqueror)) -1)/2
            if(area(old_owner).eq.1) then
              print *,conqueror,' wanted to take field ',icell,' :: but cant conquer last field of ',old_owner
              return
            endif

            !neigh_owners = ICONULL
            !do ineigh=1,3
            !  ineigh_cell = neighborcells(ineigh, icell)
            !  if(ineigh_cell.ne.ICONULL) neigh_owners(ineigh) = cellowner(ineigh_cell)
            !enddo
            call count_neighbour_owners(icell, owner_cellcnt)
            !print *,icell,'::',conqueror,':',owner_cellcnt, '=>', owner_cellcnt(conqueror)/real(sum(owner_cellcnt))

            wgt(2) = owner_cellcnt(conqueror)/real(sum(owner_cellcnt)) !real(count(neigh_owners.eq.conqueror) - count(neigh_owners.eq.ICONULL) -1.5)  ! 1 == conqueror -> negative; 2 == conqueror => positive; -1 owners dont count
            lwin = sum(wgt) .gt. .5
            !print *,conqueror,'neighbor cells:',neighborcells(:, icell),'::',neigh_owners, ':', count(neigh_owners.eq.conqueror), '=>', wgt(2)
          endif

          if( get_arg(.False., lforce) ) lwin = .True.

          if(lwin) then
            !print *,'conquering ',icell,old_owner,'->',conqueror
            if( old_owner .ne. ICONULL) area(old_owner) = area(old_owner) -1
            cellowner(icell) = conqueror
            area(conqueror) = area(conqueror) +1
          endif
        end subroutine

        subroutine count_neighbour_owners(icell, owner_cellcnt)
          integer(iintegers),intent(in) :: icell
          integer(iintegers),intent(out) :: owner_cellcnt(:) ! dim numnodes

          integer(iintegers) :: cell_list(3*3+3), ineigh, ineigh2, ineigh_cell, ineigh_cell2, iowners, owner

          owner_cellcnt = 0

          do ineigh=1,3
            ineigh_cell = neighborcells(ineigh, icell)
            if(ineigh_cell.ne.ICONULL) then
              owner = cellowner(ineigh_cell)
              if(owner.ne.ICONULL) then
                owner_cellcnt(owner) = owner_cellcnt(owner) + 2
                do ineigh2=1,3
                  ineigh_cell2 = neighborcells(ineigh2, ineigh_cell)
                  if(ineigh_cell2.ne.ICONULL .and. ineigh_cell2.ne.icell) then
                    owner = cellowner(ineigh_cell2)
                    if(owner.ne.ICONULL) owner_cellcnt(owner) = owner_cellcnt(owner) + 1
                  endif
                enddo
              endif
            endif
          enddo
        end subroutine

    end subroutine
    end module
