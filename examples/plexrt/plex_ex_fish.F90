module m_examples_plex_ex_fish

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only : read_commandline_options

  use m_helper_functions, only: &
    & angle_between_two_vec, &
    & CHKERR, &
    & cstr, &
    & deallocate_allocatable, &
    & deg2rad, &
    & determine_normal_direction, &
    & imp_bcast, &
    & imp_reduce_mean, &
    & meanval, &
    & rad2deg, &
    & reverse, &
    & spherical_2_cartesian, &
    & toStr

  use m_data_parameters, only : ireals, iintegers, mpiint, &
    default_str_len, &
    i0, i1, i2, i3, i4, i5,  &
    zero, one,       &
    init_mpi_data_parameters

  use m_icon_grid, only: t_icongrid, read_icon_grid_file, &
    bcast_icongrid, distribute_icon_grid

  use m_plex_grid, only: t_plexgrid, &
    setup_edir_dmplex, setup_abso_dmplex, compute_face_geometry, &
    ncvar2d_to_globalvec, setup_plexgrid, setup_cell1_dmplex, &
    gen_test_mat, get_normal_of_first_toa_face, get_horizontal_faces_around_vertex, &
    atm_dz_to_vertex_heights

  use m_plex_rt_base, only: t_plex_solver, allocate_plexrt_solver_from_commandline

  use m_plex_rt, only: compute_face_geometry, &
    init_plex_rt_solver, run_plex_rt_solver, set_plex_rt_optprop, &
    plexrt_get_result, destroy_plexrt_solver

  use m_netcdfio, only : ncload, ncwrite

  use m_icon_plex_utils, only: create_2d_fish_plex, create_2d_regular_plex, &
    dmplex_2D_to_3D, dump_ownership, Nz_Ncol_vec_to_celldm1

  implicit none

  logical, parameter :: ldebug=.False.

contains

  subroutine plex_ex_fish(&
      & comm, lverbose, &
      & lthermal, lsolar, &
      & lregular_mesh, &
      & Nx, Ny, Nz, &
      & dz, Ag, sundir, &
      & dtau, w0, Bplck, &
      & edir, edn, eup, abso)

    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lverbose              ! verbose output of steps and results
    logical, intent(in) :: lthermal, lsolar      ! if solar and thermal computation should be done
    logical, intent(in) :: lregular_mesh         ! if false, we use a fish mesh, if true, use a rectilinear mesh
    integer(iintegers), intent(in) :: Nx, Ny, Nz ! domain size
    real(ireals), intent(in) :: dz, Ag           ! vertical height of wedges and surface albedo
    real(ireals), intent(in) :: sundir(3)        ! cartesian direction of sun rays in a global reference system
    real(ireals), intent(in) :: dtau, w0         ! vertically integrated optical depth and assym. parameter
    real(ireals), intent(in) :: Bplck            ! constant planck emission value
    real(ireals), allocatable, dimension(:,:), intent(out) :: edir, edn, eup, abso

    type(tDM) :: dm2d, dm2d_dist, dm3d
    real(ireals) :: hhl(Nz)

    type(t_plexgrid), allocatable :: plex
    integer(iintegers), allocatable :: zindex(:)
    class(t_plex_solver), allocatable :: solver

    real(ireals) :: first_normal(3), kabs, ksca
    real(ireals) :: medir, medn, meup, mabso
    integer(iintegers) :: i, k, icol
    integer(mpiint) :: myid, numnodes, ierr

    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)
    call mpi_comm_size(comm, numnodes, ierr); call CHKERR(ierr)

    if(lregular_mesh) then
      call create_2d_regular_plex(comm, Nx, Ny, dm2d, dm2d_dist, lverbose=lverbose)
    else
      call create_2d_fish_plex(comm, Nx, Ny, dm2d, dm2d_dist, lverbose=lverbose)
    endif

    hhl(1) = zero
    do k=2,Nz
      hhl(k) = hhl(k-1) + dz
    enddo
    hhl = reverse(hhl)
    if(lverbose.and.myid.eq.0) print *,'hhl', hhl

    call dmplex_2D_to_3D(dm2d_dist, Nz, hhl, dm3d, zindex, lverbose=lverbose, lpolar_coords=.False.)

    call setup_plexgrid(dm2d_dist, dm3d, Nz-1, zindex, plex, hhl)
    deallocate(zindex)
    call DMDestroy(dm2d, ierr); call CHKERR(ierr)
    call DMDestroy(dm2d_dist, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(plex%dm, PETSC_NULL_DM, "-show_plex", ierr); call CHKERR(ierr)

    if(lregular_mesh) then
      call allocate_plexrt_solver_from_commandline(solver, 'rectilinear_5_8')
    else
      call allocate_plexrt_solver_from_commandline(solver, '5_8')
    endif
    call init_plex_rt_solver(plex, solver)

    if(lverbose.and.myid.eq.0) then
      first_normal = get_normal_of_first_TOA_face(solver%plex)
      print *,'Initial sundirection = ', sundir, ': sza', angle_between_two_vec(sundir, first_normal), 'rad'
      print *,'Initial sundirection = ', sundir, ': sza', rad2deg(angle_between_two_vec(sundir, first_normal)),'deg'
    endif

    kabs = dtau * (one - w0)
    ksca = dtau * w0
    call set_plex_rt_optprop(solver, &
      & vert_integrated_kabs=dtau * (one - w0), &
      & vert_integrated_ksca=dtau * w0, &
      & lverbose=lverbose)

    if(.not.allocated(solver%albedo)) then
      allocate(solver%albedo)
      call DMCreateGlobalVector(solver%plex%srfc_boundary_dm, solver%albedo, ierr); call CHKERR(ierr)
    endif
    call VecSet(solver%albedo, Ag, ierr); call CHKERR(ierr)

    if(lthermal) then
      if(.not.allocated(solver%plck)) then
        allocate(solver%plck)
        call DMCreateGlobalVector(solver%plex%horizface1_dm, solver%plck, ierr); call CHKERR(ierr)
      endif
      call VecSet(solver%plck, Bplck, ierr); call CHKERR(ierr)
    endif

    call run_plex_rt_solver(solver, lthermal=lthermal, lsolar=lsolar, sundir=sundir)

    call plexrt_get_result(solver, edn, eup, abso, redir=edir)

    if(lverbose) then
      do i = 0, numnodes-1
        if(myid.eq.i) then
          do icol = 1, size(abso,dim=2)
            print *,''
            print *,cstr('Column ','blue')//cstr(toStr(icol),'red')//&
              cstr('  k  Edir              Edn              Eup              abso', 'blue')
            do k = 1, ubound(abso,1)
              print *,k, edir(k,icol), edn(k,icol), eup(k,icol), abso(k,icol)
            enddo
            print *,k, edir(k,icol), edn(k,icol), eup(k,icol)
          enddo

        endif
        call mpi_barrier(comm, ierr); call CHKERR(ierr)
      enddo ! ranks
      call mpi_barrier(comm, ierr); call CHKERR(ierr)

      if(myid.eq.0) then
        print *, ''
        print *, 'Averages'
      endif
      do k = 1, ubound(abso,1)
        call imp_reduce_mean(comm, edir(k,:), medir)
        call imp_reduce_mean(comm, edn (k,:), medn)
        call imp_reduce_mean(comm, eup (k,:), meup)
        call imp_reduce_mean(comm, abso(k,:), mabso)
        if(myid.eq.0) print *,k, medir, medn, meup, mabso
      enddo
      k=size(edn,1)
      call imp_reduce_mean(comm, edir(k,:), medir)
      call imp_reduce_mean(comm, edn (k,:), medn)
      call imp_reduce_mean(comm, eup (k,:), meup)
      if(myid.eq.0) print *,k, medir, medn, meup
    endif

    call destroy_plexrt_solver(solver, lfinalizepetsc=.False.)
  end subroutine

end module
