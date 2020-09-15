module m_pprts_buildings
  use m_data_parameters, only : init_mpi_data_parameters, &
    & iintegers, ireals, mpiint, &
    & zero, one, pi, i1, i2, default_str_len

  use m_helper_functions, only : &
    & CHKERR, &
    & toStr, cstr, &
    & spherical_2_cartesian, rotate_angle_z, &
    & meanval

  use m_pprts, only : init_pprts, &
    & set_optical_properties, solve_pprts, &
    & pprts_get_result, set_angles

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

  implicit none

contains
  subroutine pprts_buildings(comm, Nx, Ny, Nlay, icollapse, dx, dy, dz, phi0, theta0, albedo, dtau, w0)
    integer(iintegers), intent(in) :: Nx,Ny,Nlay, icollapse
    real(ireals), intent(in) :: dx, dy, dz
    real(ireals), intent(in) :: phi0, theta0
    real(ireals), intent(in) :: albedo, dtau, w0
    integer(mpiint), intent(in) :: comm
    real(ireals),parameter :: incSolar = 1
    real(ireals) :: dz1d(Nlay)

    real(ireals) :: sundir(3)
    real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g
    real(ireals),allocatable,dimension(:,:,:) :: fdir,fdn,fup,fdiv

    class(t_solver), allocatable :: solver
    type(t_pprts_buildings), allocatable :: buildings
    integer, parameter :: Nbuildings = 6

    integer(iintegers) :: k, i, box_k, box_i, box_j
    integer(mpiint) :: id, myid, numnodes, ierr

    dz1d = dz

    call init_mpi_data_parameters(comm)
    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    sundir = spherical_2_cartesian(phi0, theta0)
    call init_pprts(comm, Nlay, Nx, Ny, dx,dy, sundir, solver, dz1d, collapseindex=icollapse)

    allocate(kabs(solver%C_one_atm%zm , solver%C_one_atm%xm,  solver%C_one_atm%ym ))
    allocate(ksca(solver%C_one_atm%zm , solver%C_one_atm%xm,  solver%C_one_atm%ym ))
    allocate(g   (solver%C_one_atm%zm , solver%C_one_atm%xm,  solver%C_one_atm%ym ))

    kabs = dtau*(one-w0)/(dz*Nlay)
    ksca = dtau*w0/(dz*Nlay)
    g    = zero

    call init_buildings(buildings, &
      & [integer(iintegers) :: 6, solver%C_one%zm, solver%C_one%xm,  solver%C_one%ym], &
      & ierr); call CHKERR(ierr)

    box_k = 2 !int((1+solver%C_one%zm) / 2.)
    box_i = 3 !int((1+solver%C_one%xm) / 2.)
    box_j = 3 !int((1+solver%C_one%ym) / 2.)

    print *, 'Have box:', &
      & solver%C_one%ys, solver%C_one%ye, &
      & box_k.gt.solver%C_one%zs, box_k.le.solver%C_one%ze+1, &
      & box_i.gt.solver%C_one%xs, box_i.le.solver%C_one%xe+1, &
      & box_j.gt.solver%C_one%ys, box_j.le.solver%C_one%ye+1

    if( box_k.gt.solver%C_one%zs.and.box_k.le.solver%C_one%ze+1 .and. &
      !& .False. .and. &
      & box_i.gt.solver%C_one%xs.and.box_i.le.solver%C_one%xe+1 .and. &
      & box_j.gt.solver%C_one%ys.and.box_j.le.solver%C_one%ye+1 ) then
      allocate(buildings%albedo(Nbuildings), buildings%iface(Nbuildings))
      do i=1,6
        buildings%iface(i) = faceidx_by_cell_plus_offset( &
          & buildings%da_offsets, box_k, box_i, box_j, i)
        buildings%albedo(i) = .1_ireals !+ i/10._ireals
      enddo
    else
      allocate(buildings%albedo(0), buildings%iface(0))
    endif

    call check_buildings_consistency(buildings, solver%C_one%zm, solver%C_one%xm, solver%C_one%ym, ierr); call CHKERR(ierr)

    call set_optical_properties(solver, albedo, kabs, ksca, g)
    call set_angles(solver, sundir)

    call solve_pprts(solver, &
      & lthermal=.False., &
      & lsolar=.True., &
      & edirTOA=incSolar, &
      & opt_buildings=buildings)

    call pprts_get_result(solver, fdn, fup, fdiv, fdir, opt_buildings=buildings)

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

    do id = 0, numnodes-1
      if(id.eq.myid) then
        print *,''
        print *,cstr(' ***************** Rank '//toStr(myid), 'red')
        if(box_i.gt.solver%C_one%xs.and.box_i.le.solver%C_one%xe+1) then
          print *,''
          print *,'Direct y-slice:'
          do i = box_i, box_i
            do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
              print *, 'edir', k,i, toStr( fdir(k, i, :) )
            enddo
          enddo
        endif

        if(box_j.gt.solver%C_one%ys.and.box_j.le.solver%C_one%ye+1) then
          print *,''
          print *,'Direct x-slice:'
          do i = box_j, box_j
            do k = 1+solver%C_dir%zs, 1+solver%C_dir%ze
              print *, 'edir', k,i, toStr( fdir(k, :, i) )
            enddo
          enddo
          print *,''
        endif

        if(box_i.gt.solver%C_one%xs.and.box_i.le.solver%C_one%xe+1) then
          print *,''
          print *,'Diffuse y-slice:'
          do i = box_i, box_i
            do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
              print *, 'edn', k,i, toStr( fdn(k, i, :) ), &
                & cstr(' eup'//toStr( fup(k, i, :) ), 'blue')
            enddo
          enddo
        endif

        if(box_j.gt.solver%C_one%ys.and.box_j.le.solver%C_one%ye+1) then
          print *,''
          print *,'Diffuse x-slice:'
          do i = box_j, box_j
            do k = 1+solver%C_diff%zs, 1+solver%C_diff%ze
              print *, 'edn', k,i, toStr( fdn(k, :, i) ), &
                & cstr(' eup'//toStr( fup(k, :, i) ), 'green')
            enddo
          enddo
          print *,''
        endif

        if(allocated(buildings%edir)) then
          do i=1, size(buildings%iface)
            print *, 'building_face', i, 'edir', buildings%edir(i), &
              & 'in/out', buildings%incoming(i), buildings%outgoing(i)
          enddo
        endif
      endif

      do k = lbound(fdiv,1), ubound(fdiv, 1)
        print *, k, &
          & 'mean edir', meanval(fdir(k,:,:)), &
          & 'edn', meanval(fdn(k,:,:)), &
          & 'eup', meanval(fup(k,:,:)), &
          & 'abso', meanval(fdiv(k,:,:))
      enddo
      k = ubound(fdir, 1)
      print *, k, &
        & 'mean edir', meanval(fdir(k,:,:)), &
        & 'edn', meanval(fdn(k,:,:)), &
        & 'eup', meanval(fup(k,:,:))

      call mpi_barrier(comm, ierr); call CHKERR(ierr)
    enddo

    call destroy_pprts(solver, .True.)
  end subroutine

end module


program main
#include "petsc/finclude/petsc.h"
  use petsc
  use m_pprts_buildings
  use mpi, only : MPI_COMM_WORLD
  implicit none

  integer(iintegers) :: Nx,Ny,Nlay, icollapse
  real(ireals) :: dx, dy, dz
  real(ireals) :: phi0, theta0
  real(ireals) :: Ag, dtau, w0

  character(len=10*default_str_len) :: rayli_options
  logical :: lflg, lverbose, lrayli_opts
  integer(mpiint) :: ierr

  call mpi_init(ierr)
  call init_mpi_data_parameters(mpi_comm_world)
  call read_commandline_options(mpi_comm_world)

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

  call pprts_buildings(mpi_comm_world, Nx, Ny, Nlay, icollapse, dx, dy, dz, phi0, theta0, Ag, dtau, w0)

  call mpi_finalize(ierr)
end program
