module m_examples_pprts_buildings
  use m_data_parameters, only : &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & zero, one, pi, i1, i2, i6, default_str_len

  use m_helper_functions, only : &
    & CHKERR, &
    & toStr, cstr, &
    & spherical_2_cartesian, rotate_angle_z, &
    & meanval, &
    & is_inrange

  use m_pprts, only : init_pprts, &
    & set_optical_properties, solve_pprts, &
    & pprts_get_result, set_angles, &
    & gather_all_toZero

  use m_pprts_base, only: t_solver, &
    & allocate_pprts_solver_from_commandline, destroy_pprts, &
    & t_solver_1_2, t_solver_3_6, t_solver_3_10, t_solver_3_16, &
    & t_solver_8_10, t_solver_8_16, t_solver_8_18

  use m_tenstream_options, only: read_commandline_options

  use m_buildings, only: &
    & check_buildings_consistency, &
    & faceidx_by_cell_plus_offset, &
    & init_buildings, &
    & t_pprts_buildings, &
    & PPRTS_TOP_FACE, &
    & PPRTS_BOT_FACE, &
    & PPRTS_FRONT_FACE, &
    & PPRTS_LEFT_FACE, &
    & PPRTS_REAR_FACE, &
    & PPRTS_RIGHT_FACE

  use m_xdmf_export, only: &
    & xdmf_pprts_buildings

  use m_netcdfio, only: ncwrite

  implicit none

contains
  subroutine ex_pprts_buildings(            &
      & comm, lverbose,                     &
      & lthermal, lsolar,                   &
      & Nx, Ny, Nlay, icollapse,            &
      & glob_box_i, glob_box_j, glob_box_k, &
      & box_Ni, box_Nj, box_Nk,             &
      & box_albedo, box_planck,             &
      & dx, dy, dz,                         &
      & incSolar, phi0, theta0,             &
      & albedo, dtau, w0,                   &
      & gedir, gedn, geup, gabso,           &
      & buildings,                          &
      & outfile)

    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lverbose, lthermal, lsolar
    integer(iintegers), intent(in) :: Nx, Ny, Nlay ! global domain size
    integer(iintegers), intent(in) :: icollapse    ! collapse upper nr of layers into 1 layer
    integer(iintegers), intent(in) :: glob_box_i, glob_box_j, glob_box_k ! global index of a single cube building
    integer(iintegers), intent(in) :: box_Ni, box_Nj, box_Nk ! number of layers of the building
    real(ireals), intent(in) :: box_albedo         ! albedo of building faces
    real(ireals), intent(in) :: box_planck         ! planck emission of building faces (only used if lthermal=.True.)
    real(ireals), intent(in) :: dx, dy, dz         ! grid spacing in [m]
    real(ireals), intent(in) :: incSolar           ! solar constant at TOA [W/m2]
    real(ireals), intent(in) :: phi0, theta0       ! sun azimuth(phi) and zenith(theta) angle
    real(ireals), intent(in) :: albedo, dtau, w0   ! surface albedo, vertically integrated optical depth and constant single scattering albedo
    real(ireals), allocatable, dimension(:,:,:), intent(out) :: gedir, gedn, geup, gabso ! global arrays on rank 0
    type(t_pprts_buildings), allocatable, intent(inout) :: buildings
    character(len=*), intent(in), optional :: outfile ! output file to dump flux results

    real(ireals) :: dz1d(Nlay)
    real(ireals) :: sundir(3)
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,plck
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver
    integer(iintegers) :: Nbuildings
    logical :: lhave_box

    character(len=default_str_len) :: groups(2)

    integer(iintegers) :: k, i, j, faceid, building_idx
    integer(iintegers) :: box_k, box_i, box_j
    integer(mpiint) :: myid, numnodes, ierr

    dz1d = dz

    call init_mpi_data_parameters(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    sundir = spherical_2_cartesian(phi0, theta0)
    call init_pprts(comm, Nlay, Nx, Ny, dx,dy, sundir, solver, dz1d, collapseindex=icollapse)

    associate(Ca => solver%C_one_atm, C1 => solver%C_one)
      allocate(kabs(Ca%zm  , Ca%xm, Ca%ym ))
      allocate(ksca(Ca%zm  , Ca%xm, Ca%ym ))
      allocate(g   (Ca%zm  , Ca%xm, Ca%ym ))

      if(lthermal) then
        allocate(plck(Ca%zm+1, Ca%xm, Ca%ym ))
        plck(:,:,:) = 0
        plck(Ca%zm+1,:,:) = 0
      endif

      kabs = dtau*(one-w0)/dz/real(Nlay, ireals)
      ksca = dtau*w0/dz/real(Nlay, ireals)
      g    = zero

      box_k = glob_box_k - C1%zs
      box_i = glob_box_i - C1%xs
      box_j = glob_box_j - C1%ys

      do i=-box_Ni+1, box_Ni-1
        do j=-box_Nj+1, box_Nj-1
          if(lverbose) print *, myid, 'Have box:', &
            & i,j,':',&
            & is_inrange(glob_box_i+i, C1%xs+1, C1%xe+1), &
            & is_inrange(glob_box_j+j, C1%ys+1, C1%ye+1)
        enddo
      enddo

      lhave_box = .False.
      Nbuildings = 0
      do i=-box_Ni+1, box_Ni-1
        do j=-box_Nj+1, box_Nj-1
          if( &
            & is_inrange(glob_box_i+i, C1%xs+1, C1%xe+1).and. &
            & is_inrange(glob_box_j+j, C1%ys+1, C1%ye+1)      ) then

            lhave_box = .True.
            Nbuildings = Nbuildings + 6 * box_Nk
          endif
        enddo
      enddo
      print *,myid,': Nbuildings',Nbuildings

      call init_buildings(buildings, &
        & [integer(iintegers) :: 6, C1%zm, C1%xm,  C1%ym], &
        & Nbuildings, &
        & ierr); call CHKERR(ierr)

      if(lthermal) allocate(buildings%planck(Nbuildings))
      building_idx = 1
      do i=-box_Ni+1, box_Ni-1
        do j=-box_Nj+1, box_Nj-1
          if( &
            & is_inrange(glob_box_i+i, C1%xs+1, C1%xe+1).and. &
            & is_inrange(glob_box_j+j, C1%ys+1, C1%ye+1)      ) then
            do k = 0, box_Nk-1
              do faceid = 1, 6
                buildings%iface(building_idx) = faceidx_by_cell_plus_offset( &
                  & buildings%da_offsets, &
                  & box_k-k, &
                  & box_i+i, &
                  & box_j+j, faceid)
                buildings%albedo(building_idx) = box_albedo
                if(lthermal) buildings%planck(building_idx) = box_planck
                !print *, building_idx, 'building kij',box_k-k,box_i+i,box_j+j, faceid, &
                !  & 'iface', buildings%iface(building_idx), 'albedo', buildings%albedo(building_idx)
                building_idx = building_idx+1
              enddo
            enddo
          endif
        enddo
      enddo

      call check_buildings_consistency(buildings, C1%zm, C1%xm, C1%ym, ierr); call CHKERR(ierr)

      if(lthermal) then
        call set_optical_properties(solver, albedo, kabs, ksca, g, plck)
      else
        call set_optical_properties(solver, albedo, kabs, ksca, g)
      endif
      call set_angles(solver, sundir)

      call solve_pprts(solver, &
        & lthermal=lthermal, &
        & lsolar=lsolar, &
        & edirTOA=incSolar, &
        & opt_buildings=buildings)

      call pprts_get_result(solver, fdn, fup, fdiv, fdir, opt_buildings=buildings)

      if(lsolar) then
        call gather_all_toZero(solver%C_one_atm1, fdir, gedir)
      endif
      call gather_all_toZero(solver%C_one_atm1, fdn, gedn)
      call gather_all_toZero(solver%C_one_atm1, fup, geup)
      call gather_all_toZero(solver%C_one_atm, fdiv, gabso)

      if(myid.eq.0_mpiint.and.present(outfile)) then
        groups(1) = trim(outfile)
        if(lsolar) then
          groups(2) = 'edir'; call ncwrite(groups, gedir, ierr); call CHKERR(ierr)
        endif
        groups(2) = 'edn' ; call ncwrite(groups, gedn , ierr); call CHKERR(ierr)
        groups(2) = 'eup' ; call ncwrite(groups, geup , ierr); call CHKERR(ierr)
        groups(2) = 'abso'; call ncwrite(groups, gabso, ierr); call CHKERR(ierr)
      endif

      if(lverbose .and. myid.eq.0) then
        print *,'y-slice: i='//toStr(box_i)
        i = box_i
        if(lsolar) then
          do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
            print *, k, cstr('edir'//toStr( gedir(k, i, :) ), 'red')
          enddo
        endif ! lsolar
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' edn'//toStr( gedn(k, i, :) ), 'green')
        enddo
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' eup'//toStr( geup(k, i, :) ), 'blue')
        enddo
        do k = 1+solver%C_one%zs, 1+solver%C_one%ze
          print *, k, cstr('abso'//toStr( gabso(k, i, :) ), 'purple')
        enddo

        print *,''
        print *,'x-slice: j='//toStr(box_j)
        i = box_j
        if(lsolar) then
          do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
            print *, k, cstr('edir'//toStr( gedir(k, :, i) ), 'red')
          enddo
        endif
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' edn'//toStr( gedn(k, :, i) ), 'green')
        enddo
        do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
          print *, k, cstr(' eup'//toStr( geup(k, :, i) ), 'blue')
        enddo
        do k = 1+solver%C_one%zs, 1+solver%C_one%ze
          print *, k, cstr('abso'//toStr( gabso(k, :, i) ), 'purple')
        enddo

        print *,''
        if(lsolar) then
          do k = lbound(fdiv,1), ubound(fdiv, 1)
            print *, k, 'mean ', &
              & 'edir', meanval(gedir(k,:,:)), &
              & 'edn' , meanval(gedn(k,:,:)), &
              & 'eup' , meanval(geup(k,:,:)), &
              & 'abso', meanval(gabso(k,:,:))
          enddo
          k = ubound(fdir, 1)
          print *, k, 'mean ', &
            & 'edir', meanval(gedir(k,:,:)), &
            & 'edn' , meanval(gedn(k,:,:)), &
            & 'eup' , meanval(geup(k,:,:))
        else
          do k = lbound(fdiv,1), ubound(fdiv, 1)
            print *, k, &
              & 'mean ', &
              & 'edn', meanval(gedn(k,:,:)), &
              & 'eup', meanval(geup(k,:,:)), &
              & 'abso', meanval(gabso(k,:,:))
          enddo
          k = ubound(fdir, 1)
          print *, k, 'mean ', &
            & 'edn', meanval(gedn(k,:,:)), &
            & 'eup', meanval(geup(k,:,:))
        endif
      endif
    end associate
    call mpi_barrier(comm, ierr); call CHKERR(ierr)

    if(lverbose .and. lhave_box) then
      print *,''
      if(allocated(buildings%edir)) then
        do i=1, size(buildings%iface)
          print *, 'building_face', i, 'edir', buildings%edir(i), &
            & 'in/out', buildings%incoming(i), buildings%outgoing(i)
        enddo
      else
        do i=1, size(buildings%iface)
          print *, 'building_face', i, &
            & 'in/out', buildings%incoming(i), buildings%outgoing(i)
        enddo
      endif
    endif

    call destroy_pprts(solver, .False.)
  end subroutine

end module
