module m_pprts_buildings
  use m_data_parameters, only : &
    & init_mpi_data_parameters, &
    & finalize_mpi, &
    & iintegers, ireals, mpiint, &
    & zero, one, pi, i1, i2, default_str_len

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

  use m_buildings, only: t_pprts_buildings, &
    & faceidx_by_cell_plus_offset, &
    & check_buildings_consistency, &
    & PPRTS_TOP_FACE, PPRTS_BOT_FACE, &
    & PPRTS_LEFT_FACE, PPRTS_RIGHT_FACE, &
    & PPRTS_REAR_FACE, PPRTS_FRONT_FACE, &
    & init_buildings

  use m_netcdfio, only: ncwrite

  implicit none

contains
  subroutine pprts_buildings(comm, &
      & outfile, &
      & lthermal, lsolar, &
      & Nx, Ny, Nlay, icollapse, &
      & dx, dy, dz, phi0, theta0, &
      & albedo, dtau, w0)
    integer(mpiint), intent(in) :: comm
    character(len=*), intent(in) :: outfile
    logical, intent(in) :: lthermal, lsolar
    integer(iintegers), intent(in) :: Nx,Ny,Nlay, icollapse
    real(ireals), intent(in) :: dx, dy, dz
    real(ireals), intent(in) :: phi0, theta0
    real(ireals), intent(in) :: albedo, dtau, w0
    real(ireals),parameter :: incSolar = 1
    real(ireals) :: dz1d(Nlay)

    real(ireals) :: sundir(3)
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,plck
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv
    real(ireals),allocatable, dimension(:,:,:) :: gedir, gedn, geup, gabso ! global arrays which we will dump to netcdf

    class(t_solver), allocatable :: solver
    type(t_pprts_buildings), allocatable :: buildings
    integer :: Nbuildings
    logical :: lhave_box

    character(len=default_str_len) :: groups(2)

    integer(iintegers) :: k, i
    integer(iintegers) :: glob_box_k, glob_box_i, glob_box_j
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

      kabs = dtau*(one-w0)/(dz*Nlay)
      ksca = dtau*w0/(dz*Nlay)
      g    = zero

      glob_box_k = C1%zm-1 ! one above the surface-touching cell
      glob_box_i = 7 !int((1+Ca%xm) / 2.)
      glob_box_j = 7 !int((1+Ca%ym) / 2.)

      box_k = glob_box_k - C1%zs
      box_i = glob_box_i - C1%xs
      box_j = glob_box_j - C1%ys

      print *, myid, 'Have box:', &
        & is_inrange(glob_box_k, C1%zs+1, C1%ze), &
        & is_inrange(glob_box_i, C1%xs+1, C1%xe), &
        & is_inrange(glob_box_j, C1%ys+1, C1%ye)

      if( &
        !& .False. .and. &
        & is_inrange(glob_box_k, C1%zs+1, C1%ze).and. &
        & is_inrange(glob_box_i, C1%xs+1, C1%xe).and. &
        & is_inrange(glob_box_j, C1%ys+1, C1%ye)      ) then

        lhave_box = .True.
        Nbuildings = 6
      else
        lhave_box = .False.
        Nbuildings = 0
      endif

      call init_buildings(buildings, &
        & [integer(iintegers) :: 6, C1%zm, C1%xm,  C1%ym], &
        & Nbuildings, &
        & ierr); call CHKERR(ierr)

      if(lthermal) allocate(buildings%planck(Nbuildings))
      do i=1,Nbuildings
        buildings%iface(i) = faceidx_by_cell_plus_offset( &
          & buildings%da_offsets, &
          & box_k, &
          & box_i, &
          & box_j, i)
        buildings%albedo(i) = .1_ireals
        if(lthermal) buildings%planck(i) = 100
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

      groups(1) = trim(outfile)

      if(lsolar) then
        call gather_all_toZero(solver%C_one_atm1, fdir, gedir)
      endif
      call gather_all_toZero(solver%C_one_atm1, fdn, gedn)
      call gather_all_toZero(solver%C_one_atm1, fup, geup)
      call gather_all_toZero(solver%C_one_atm, fdiv, gabso)

      if(myid.eq.0) then
        if(lsolar) then
          groups(2) = 'edir'; call ncwrite(groups, gedir, ierr); call CHKERR(ierr)
        endif
        groups(2) = 'edn' ; call ncwrite(groups, gedn , ierr); call CHKERR(ierr)
        groups(2) = 'eup' ; call ncwrite(groups, geup , ierr); call CHKERR(ierr)
        groups(2) = 'abso'; call ncwrite(groups, gabso, ierr); call CHKERR(ierr)
      endif

      if(myid.eq.0) then
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

    if(lhave_box) then
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


program main
#include "petsc/finclude/petsc.h"
  use petsc
  use m_pprts_buildings
  use mpi, only : MPI_COMM_WORLD
  implicit none

  character(len=default_str_len) :: outfile
  integer(iintegers) :: Nx,Ny,Nlay, icollapse
  real(ireals) :: dx, dy, dz
  real(ireals) :: phi0, theta0
  real(ireals) :: Ag, dtau, w0

  character(len=10*default_str_len) :: rayli_options
  logical :: lflg, lverbose, lrayli_opts, lsolar, lthermal
  integer(mpiint) :: ierr

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call read_commandline_options(mpi_comm_world)

  call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-out', outfile, lflg, ierr); call CHKERR(ierr)
  if(.not.lflg) call CHKERR(1_mpiint, 'need to supply a output filename... please call with -out <output.nc>')

  lsolar = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-solar', &
    lsolar, lflg, ierr) ; call CHKERR(ierr)

  lthermal = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-thermal', &
    lthermal, lflg, ierr) ; call CHKERR(ierr)

  Nx=5; Ny=5; Nlay=3
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nlay, lflg, ierr); call CHKERR(ierr)

  icollapse=1
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-icollapse", icollapse, lflg, ierr); call CHKERR(ierr)

  dx = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
  dy = dx
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr); call CHKERR(ierr)
  dz = 100
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dz", dz, lflg, ierr); call CHKERR(ierr)

  phi0 = 180._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi", phi0, lflg, ierr); call CHKERR(ierr)
  theta0 = 0._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta", theta0, lflg, ierr); call CHKERR(ierr)

  Ag=0.1_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ag", Ag, lflg, ierr); call CHKERR(ierr)
  dtau=1._ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dtau", dtau, lflg, ierr); call CHKERR(ierr)
  w0=0.5_ireals
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-w0", w0, lflg, ierr); call CHKERR(ierr)

  lverbose = .True.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-verbose', &
    lverbose, lflg, ierr) ; call CHKERR(ierr)

  lrayli_opts = .False.
  call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-rayli_opts', &
    lrayli_opts, lflg, ierr) ; call CHKERR(ierr)

  if(lrayli_opts) then
    rayli_options=''
    rayli_options=trim(rayli_options)//' -pprts_use_rayli'
    rayli_options=trim(rayli_options)//' -rayli_diff_flx_origin 0,0,-inf'
    rayli_options=trim(rayli_options)//' -rayli_cyclic_bc'
    rayli_options=trim(rayli_options)//' -show_rayli_dm3d hdf5:dm.h5'

    rayli_options=trim(rayli_options)//' -rayli_snapshot'
    rayli_options=trim(rayli_options)//' -rayli_snap_Nx 256'
    rayli_options=trim(rayli_options)//' -rayli_snap_Ny 256'
    rayli_options=trim(rayli_options)//' -visit_image_zoom .75'
    rayli_options=trim(rayli_options)//' -visit_parallel_scale 291.5'
    rayli_options=trim(rayli_options)//' -visit_focus 300,300,0'
    rayli_options=trim(rayli_options)//' -visit_view_normal -0.2811249083446944,-0.7353472951470268,0.6166304739697339'
    rayli_options=trim(rayli_options)//' -visit_view_up 0.1878717450780742,0.5879401184069877,0.7867849925925738'

    if(lverbose) print *,'Adding rayli Petsc Options:', trim(rayli_options)
    call PetscOptionsInsertString(PETSC_NULL_OPTIONS, trim(rayli_options), ierr); call CHKERR(ierr)
  endif

  if(lsolar) then
    call pprts_buildings(mpi_comm_world, outfile, &
      & .False., lsolar, Nx, Ny, Nlay, icollapse, dx, dy, dz, phi0, theta0, Ag, dtau, w0)
  endif
  if(lthermal) then
    call pprts_buildings(mpi_comm_world, outfile, &
      & lsolar, .False., Nx, Ny, Nlay, icollapse, dx, dy, dz, phi0, theta0, Ag, dtau, w0)
  endif

  call finalize_mpi(ierr)
end program
