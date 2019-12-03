module m_icon_grid

#include "petsc/finclude/petsc.h"
  use petsc
  use m_netcdfIO, only: ncload

  use m_data_parameters, only : ireals, iintegers, mpiint, default_str_len, &
    i0, i1, i2

  use m_helper_functions, only: get_arg, imp_bcast, chkerr, unique, cumsum

  implicit none

  private
  public :: t_icongrid, distribute_icon_grid, read_icon_grid_file, bcast_icongrid, ICONULL

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

  logical, parameter :: ldebug=.False.
  integer(iintegers), parameter :: ICONULL=-1

  contains
    subroutine distribute_icon_grid(comm, icongrid, local_icongrid, cell_ao )
      integer(mpiint), intent(in) :: comm
      type(t_icongrid), allocatable, intent(in) :: icongrid
      type(t_icongrid), allocatable, intent(out) :: local_icongrid

      AO, intent(out) :: cell_ao

      integer(mpiint) :: myid, numnodes, ierr

      integer(iintegers) :: i, j, ic, ie, iv, ilocal
      integer(iintegers) :: owner
      logical,allocatable :: ladj_cell(:)
      integer(iintegers),allocatable :: par_cnt(:), iowner(:), unique_owner(:)

      allocate(unique_owner(0))

      if(.not.allocated(icongrid)) stop 'distribute_icon_grid :: global icongrid not allocated!'
      allocate(local_icongrid)

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(ldebug.and.myid.eq.0) print *,'distribute_icon_grid :: starting'
      call decompose_icon_grid_parmetis(comm, icongrid, &
        local_icongrid%cellowner, local_icongrid%edgeowner, local_icongrid%vertexowner, cell_ao)

      allocate(ladj_cell(0:numnodes-1))

      allocate(local_icongrid%parNfaces(0:numnodes-1))
      do owner=0,numnodes-1
        local_icongrid%parNfaces(owner) = count(local_icongrid%cellowner.eq.owner)
      enddo
      local_icongrid%Nfaces = local_icongrid%parNfaces(myid)

      ! Count edges/vertices for local grid
      allocate(local_icongrid%parNedges(0:numnodes-1), source=i0)
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

      allocate(local_icongrid%parNvertices(0:numnodes-1), source=i0)
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
        !print *,myid,i,'global cell', ic, ' to local ->', ilocal, 'remote', local_icongrid%remote_cell_index(ic)
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
        !print *,myid,i,'global edges', ie, ' to local ->', ilocal, 'remote', local_icongrid%remote_edge_index(ie)
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
            !print *,myid,'distribute_icon_grid :: my cells global indices ', local_icongrid%cell_index
            !print *,myid,'distribute_icon_grid :: my edge global indices  ', local_icongrid%edge_index
            !print *,myid,'distribute_icon_grid :: my vertex global indices', local_icongrid%vertex_index

            do i=1,size(local_icongrid%cell_index)
              ic = local_icongrid%cell_index(i)       ! global cell index in parent grid
              j  = local_icongrid%local_cell_index(ic) ! local cell index
              !print *,myid,'global 2 local cell', ic, '->', j, ' ::: ', i, '<-', ic
            enddo
            do i=1,size(local_icongrid%edge_index)
              ie = local_icongrid%edge_index(i)       ! global edge index in parent grid
              j  = local_icongrid%local_edge_index(ie) ! local edge index
              !print *,myid,'global 2 local edge', ie, '->', j, ' ::: ', i, '<-', ie
            enddo

            do i=1,size(local_icongrid%cell_index)
              ic = local_icongrid%cell_index(i)        ! global cell index in parent grid
              j  = local_icongrid%local_cell_index(ic) ! local cell index
              !print *,myid,'edge_of_cell', ic, ': local cell', i, '::',icongrid%edge_of_cell(ic,:),'->',local_icongrid%edge_of_cell(j,:)
            enddo
            do i=1,size(local_icongrid%cell_index)
              ic = local_icongrid%cell_index(i)       ! global cell index in parent grid
              j  = local_icongrid%local_cell_index(ic) ! local cell index
              !print *,myid,'vertex_of_cell', ic, '::', i,'->',local_icongrid%vertex_of_cell(j,:)
            enddo
          endif
          call mpi_barrier(comm, ierr)
        enddo
        call mpi_barrier(comm, ierr)
      endif
      if(ldebug.and.myid.eq.0) print *,'distribute_icon_grid :: finished'
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

    subroutine decompose_icon_grid_parmetis(comm, icongrid, cellowner, edgeowner, vertexowner, cell_ao)
      MPI_Comm, intent(in) :: comm
      type(t_icongrid),intent(in) :: icongrid

      integer(iintegers),allocatable,intent(out) :: cellowner(:)   ! dim=(icongrid%Nfaces)
      integer(iintegers),allocatable,intent(out) :: edgeowner(:)   ! dim=(icongrid%Nedges)
      integer(iintegers),allocatable,intent(out) :: vertexowner(:) ! dim=(icongrid%Nvertices)

      AO, intent(out) :: cell_ao

      integer(iintegers) :: i, icell, j, ivertex, offset

      integer(iintegers) :: ncells_local, istartcell
      integer(iintegers),allocatable :: ii(:), jj(:)

      integer(iintegers),allocatable :: cells_per_proc(:), cum_cells_per_proc(:)

      integer(mpiint) :: myid, numnodes, ierr

      type(tMat) :: mesh, dual
      type(tIS)  :: is, isg, is_my_icon_cells
      MatPartitioning :: part

      if(.not. allocated(icongrid%cell_index)) &
        call CHKERR(1_mpiint, 'decompose_icon_grid :: the grid is not loaded from an icon gridfile...'// &
          'cant be decomposed with this method')

      call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
      call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

      if(ldebug.and.myid.eq.0) print *,'decompose_icon_grid_parmetis :: starting'

      allocate(cells_per_proc(0:numnodes-1))
      allocate(cum_cells_per_proc(0:numnodes-1))

      allocate(cellowner(icongrid%Nfaces))
      allocate(edgeowner(icongrid%Nedges))
      allocate(vertexowner(icongrid%Nvertices))

      ncells_local = icongrid%Nfaces / numnodes
      if(myid.eq.numnodes-1) ncells_local = ncells_local + modulo(icongrid%Nfaces, int(numnodes, kind=iintegers))  ! last rank gets the rest

      allocate(ii(0:ncells_local), source=ICONULL)
      allocate(jj(0:ncells_local * 3 -1), source=ICONULL)

      istartcell = icongrid%Nfaces / numnodes * myid

      offset = 0
      do i=1,ncells_local
        ii(i-1) = offset
        icell = icongrid%cell_index(istartcell + i)
        do j=1,size(icongrid%vertex_of_cell, dim=2)
          ivertex = icongrid%vertex_of_cell(icell, j)
          jj(offset) = ivertex -1
          offset = offset +1
        enddo
      enddo
      ii(ncells_local) = offset

      call MatCreateMPIAdj(comm, ncells_local , icongrid%Nvertices, &
        ii, jj, PETSC_NULL_INTEGER, mesh, ierr); call CHKERR(ierr)

      call MatMeshToCellGraph(mesh, i2, dual, ierr);

      call MatPartitioningCreate(MPI_COMM_WORLD,part, ierr); call CHKERR(ierr)
      call MatPartitioningSetAdjacency(part,dual, ierr); call CHKERR(ierr)
      call MatPartitioningSetFromOptions(part, ierr); call CHKERR(ierr)

      call MatPartitioningApply(part,is, ierr); call CHKERR(ierr)
      call ISPartitioningToNumbering(is, isg, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(is, PETSC_NULL_IS, "-show_is", ierr); call CHKERR(ierr)
      call PetscObjectViewFromOptions(isg, PETSC_NULL_IS, "-show_isg", ierr); call CHKERR(ierr)

      call ISPartitioningCount(is, int(numnodes, kind=iintegers), cells_per_proc, ierr); call CHKERR(ierr)
      cum_cells_per_proc = cumsum(cells_per_proc)

      call ISInvertPermutation(isg, cells_per_proc(myid), is_my_icon_cells, ierr); call CHKERR(ierr)

      call PetscObjectViewFromOptions(is_my_icon_cells, PETSC_NULL_IS, "-show_is_my_icon_cells", ierr); call CHKERR(ierr)

      call AOCreateBasicIS(is_my_icon_cells,PETSC_NULL_IS, cell_ao, ierr); call CHKERR(ierr)

      do i=1,icongrid%Nfaces
        j = i-1
        call AOPetscToApplication(cell_ao, i1, j, ierr); call CHKERR(ierr)
        cellowner(j+1) = count(cum_cells_per_proc/i.eq.0)
        !print *,'Owner of cell:',i,'->',j,'@',count(cum_cells_per_proc/i.eq.0)
      enddo

      call ISDestroy(is, ierr);
      call ISDestroy(isg, ierr);
      call ISDestroy(is_my_icon_cells, ierr);
      call MatPartitioningDestroy(part, ierr);

      call MatDestroy(mesh, ierr);
      call MatDestroy(dual, ierr);

      call distribute_edges()
      call distribute_vertices()

      if(ldebug.and.myid.eq.0) print *,'decompose_icon_grid_parmetis :: finished'
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
        end subroutine
    end subroutine
end module
