module m_examples_pprts_specint_tree
  ! Example to compute with additional, optional optical properties with specint.
  ! This allows to add other species with custom optical properties to scenes.
  ! For example necessary for aerosols or in this case vegetation.
  ! The optical properties will given from the vegetation module,
  ! which provides USGS spectral reflectance data.
  ! The clear sky scene has a tree like structure in the center of the domain
  !   - with a straight column of 'bark'
  !   - and a cubic radial canopy of 'leaf'

#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only: &
    & init_mpi_data_parameters, &
    & iintegers, ireals, mpiint, &
    & default_str_len

  use m_helper_functions, only: &
    & CHKERR, &
    & cstr, &
    & domain_decompose_2d_petsc, &
    & get_arg, &
    & is_inrange, &
    & linspace, &
    & spherical_2_cartesian, &
    & toStr

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only: t_solver, allocate_pprts_solver_from_commandline

  use m_pprts, only: gather_all_to_all

  ! main entry point for solver, and desctructor
  use m_specint_pprts, only: specint_pprts, specint_pprts_destroy

  ! tenstr_atm holds info about tracer and merges dynamics grid vars with background grids
  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, destroy_tenstr_atm

  use m_vegetation_optprop, only: &
    & init_vegetation_types_simple, &
    & init_vegetation_types_from_usgs, &
    & get_veg_type_id, &
    & get_albedo_for_range

  implicit none

contains
  subroutine ex_pprts_specint_tree(         &
      & specint,                            &
      & comm, lverbose,                     &
      & lthermal, lsolar,                   &
      & Nx, Ny, Nlay,                       &
      & dx, dy,                             &
      & atm_filename,                       &
      & phi0, theta0,                       &
      & Ag_solar, Ag_thermal,               &
      & Ntree_height,                       &
      & luse_usgs_db,                       &
      & gedir, gedn, geup, gabso,           &
      & local_dims, icollapse               &
      & )

    character(len=*), intent(in) :: specint           ! name of module to use for spectral integration
    integer(mpiint), intent(in) :: comm
    logical, intent(in) :: lverbose, lthermal, lsolar
    integer(iintegers), intent(in) :: Nx, Ny, Nlay   ! global domain size
    real(ireals), intent(in) :: dx, dy               ! grid spacing in [m]
    character(len=*), intent(in) :: atm_filename     ! e.g. 'afglus_100m.dat'
    real(ireals), intent(in) :: phi0, theta0         ! sun azimuth(phi) and zenith(theta) angle
    real(ireals), intent(in) :: Ag_solar, Ag_thermal ! surface albedo
    integer(iintegers), intent(in) :: Ntree_height   ! Nr of vertical layers of the tree
    logical, intent(in) :: luse_usgs_db ! enable to load USGS spectral material database. Otherwise use simple version

    real(ireals), allocatable, dimension(:, :, :), intent(out) :: gedir, gedn, geup, gabso
    integer(iintegers), intent(out), optional :: local_dims(:) ! local domain indices (zs, zm, xs, xm, ys, ym), dim(6)
    integer(iintegers), intent(in), optional :: icollapse

    integer(iintegers) :: nxp, nyp, xs, ys ! local domain size in x and y aswell as start indices
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    real(ireals) :: sundir(3) ! cartesian vector pointing from sun to earth

    real(ireals), allocatable, dimension(:, :, :), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), allocatable, dimension(:, :, :), target :: tlev ! Temperature on layer interfaces [K]

    real(ireals), allocatable, dimension(:, :, :) :: tree_tau_solar, tree_w0_solar, tree_tau_thermal

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:, :) :: pplev, ptlev

    real(ireals), allocatable, dimension(:, :, :) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    class(t_solver), allocatable :: solver
    type(t_tenstr_atm) :: atm

    integer(iintegers) :: k
    integer(mpiint) :: myid, ierr

    call init_mpi_data_parameters(comm)
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(solver, '3_10', ierr); call CHKERR(ierr)

    call domain_decompose_2d_petsc(comm, Nx, Ny, &
      & nxp, nyp, xs, ys, nxproc, nyproc, ierr); call CHKERR(ierr)

    sundir = spherical_2_cartesian(phi0, theta0)
    if (lverbose) print *, 'sundir:', sundir

    allocate (plev(Nlay + 1, nxp, nyp))
    allocate (tlev(Nlay + 1, nxp, nyp))

    ! Start with a dynamics grid ranging from 1000 hPa up to 900 hPa and a
    ! Temperature difference of 10K
    do k = 1, Nlay + 1
      plev(k, :, :) = linspace(k, [1e3_ireals, 900._ireals], Nlay + 1)
      tlev(k, :, :) = linspace(k, [288._ireals, 278._ireals], Nlay + 1)
    end do

    pplev(1:size(plev, 1), 1:size(plev, 2) * size(plev, 3)) => plev
    ptlev(1:size(tlev, 1), 1:size(tlev, 2) * size(tlev, 3)) => tlev

    call setup_tenstr_atm( &
      & comm, .false.,     &
      & atm_filename,      &
      & pplev, ptlev,      &
      & atm)

    ! only init grid structures
    call specint_pprts(        &
      & specint,               &
      & comm,                  &
      & solver, atm,           &
      & nxp, nyp, dx, dy,      &
      & sundir,                &
      & Ag_thermal, Ag_solar,  &
      & lthermal, lsolar,      &
      & edir, edn, eup, abso,  &
      & icollapse=get_arg(1_iintegers, icollapse), &
      & nxproc=nxproc,         &
      & nyproc=nyproc,         &
      & lonly_initialize=.true.)

    call build_tree(Ntree_height)

    if (present(local_dims)) then
      local_dims(:) = [ &
        & solver%C_one%zs, solver%C_one%zm, &
        & solver%C_one%xs, solver%C_one%xm, &
        & solver%C_one%ys, solver%C_one%ym]
    end if

    ! call solver
    call specint_pprts(                             &
      & specint,                                    &
      & comm,                                       &
      & solver, atm,                                &
      & nxp, nyp, dx, dy,                           &
      & sundir,                                     &
      & Ag_thermal, Ag_solar,                       &
      & lthermal, lsolar,                           &
      & edir, edn, eup, abso,                       &
      & icollapse=get_arg(1_iintegers, icollapse),  &
      & nxproc=nxproc,                              &
      & nyproc=nyproc,                              &
      & opt_tau_solar=tree_tau_solar,               &
      & opt_w0_solar=tree_w0_solar,                 &
      & opt_tau_thermal=tree_tau_thermal            &
      & )

    if (allocated(edir)) &
      & call gather_all_to_all(solver%C_one1, edir, gedir)
    call gather_all_to_all(solver%C_one1, edn, gedn)
    call gather_all_to_all(solver%C_one1, eup, geup)
    call gather_all_to_all(solver%C_one, abso, gabso)

    ! Tidy up
    call specint_pprts_destroy(specint, solver, lfinalizepetsc=.true., ierr=ierr)
    call destroy_tenstr_atm(atm)

  contains

    ! Set trunk in the global domain center
    ! with Ntree_height layers sitting on the surface
    ! And add a canopy with radius on the top
    subroutine build_tree(Ntree_height)
      integer(iintegers), intent(in) :: Ntree_height
      integer(iintegers) :: center(3) ! global indices (i,j), of center box of domain, and (3) the vertical center of the leaves
      integer(iintegers) :: i, j, k, vid
      integer(iintegers) :: li, lj, lk
      real(ireals) :: radius, albedo, lambda_start, lambda_end
      real(ireals) :: LAI ! leaf area index
      real(ireals), parameter :: LAI_leaf = 1, LAI_bark = .1_ireals

      if (luse_usgs_db) then
        call init_vegetation_types_from_usgs(comm, ierr, verbose=lverbose); call CHKERR(ierr)
      else
        call init_vegetation_types_simple(ierr); call CHKERR(ierr)
      end if

      if (Nlay .le. Ntree_height) &
        & call CHKERR(1_mpiint, 'cant fit tree in domains so small: need nz > '//toStr(Ntree_height)//&
          & ' but have nz='//toStr(Nlay))
      if (Nx .lt. Ntree_height + 1) &
        & call CHKERR(1_mpiint, 'cant fit tree in domains so small: need nx >= '//toStr(Ntree_height + 1)//&
          & ' but have nx='//toStr(Nx))
      if (Ny .lt. Ntree_height + 1) &
        & call CHKERR(1_mpiint, 'cant fit tree in domains so small: need ny >= '//toStr(Ntree_height + 1)//&
          & ' but have ny='//toStr(Ny))

      associate (C1 => solver%C_one)
        allocate ( &
          & tree_tau_solar(Ntree_height, C1%xm, C1%ym), &
          & tree_w0_solar(Ntree_height, C1%xm, C1%ym), &
          & tree_tau_thermal(Ntree_height, C1%xm, C1%ym), &
          & source=0._ireals)

        center(1) = int(real(C1%glob_xm + 1) / 2.) ! midpoint of domain
        center(2) = int(real(C1%glob_ym + 1) / 2.) ! midpoint of domain
        center(3) = C1%glob_zm - Ntree_height + 1

        ! Solar wavelengths
        lambda_start = 450
        lambda_end = 1200

        ! trunk
        i = center(1)
        j = center(2)
        do k = C1%glob_zm - Ntree_height + 1, C1%glob_zm ! lowermost layer
          if (have_box(k, i, j)) then
            li = i - solver%C_one%xs ! local Fortran index
            lj = j - solver%C_one%ys ! local Fortran index
            lk = C1%glob_zm - k + 1 ! local Fortran index starting at bot

            if (luse_usgs_db) then
              vid = get_veg_type_id('WhitebarkPine YNP-WB-1 frst AVIRISb RTGC')
            else
              vid = get_veg_type_id('bark')
            end if
            LAI = LAI_bark
            albedo = get_albedo_for_range(vid, lambda_start, lambda_end)

            tree_w0_solar(lk, li, lj) = albedo
            tree_tau_solar(lk, li, lj) = LAI
          end if
        end do

        ! canopy

        do j = center(2) - Ntree_height / 2, center(2) + Ntree_height / 2
          do i = center(1) - Ntree_height / 2, center(1) + Ntree_height / 2
            do k = C1%glob_zm - Ntree_height, C1%glob_zm
              radius = sqrt(real(i - center(1))**2 + real(j - center(2))**2 + 4 * real(k - center(3))**2)
              if (radius .le. real(Ntree_height, ireals) / 3._ireals .and. have_box(k, i, j)) then
                li = i - solver%C_one%xs
                lj = j - solver%C_one%ys
                lk = C1%glob_zm - k + 1

                if (luse_usgs_db) then
                  vid = get_veg_type_id('Oak Oak-Leaf-1 fresh         ASDFRa AREF')
                else
                  vid = get_veg_type_id('leaf')
                end if
                LAI = LAI_leaf
                albedo = get_albedo_for_range(vid, lambda_start, lambda_end)

                ! mixing bark and leaf albedi
                tree_w0_solar(lk, li, lj) = &
                  & (albedo * LAI + tree_w0_solar(lk, li, lj) * tree_tau_solar(lk, li, lj)) &
                  & / (LAI + tree_tau_solar(lk, li, lj))

                tree_tau_solar(lk, li, lj) = tree_tau_solar(lk, li, lj) + LAI
              end if
            end do
          end do
        end do

        ! Thermal part
        lambda_start = 4000
        lambda_end = 40000

        ! trunk
        i = center(1)
        j = center(2)
        do k = C1%glob_zm - Ntree_height + 1, C1%glob_zm ! lowermost layer
          if (have_box(k, i, j)) then
            li = i - solver%C_one%xs ! F_idx
            lj = j - solver%C_one%ys ! F_idx
            lk = C1%glob_zm - k + 1 ! F_idx starting at bot

            LAI = LAI_bark
            tree_tau_thermal(lk, li, lj) = LAI
          end if
        end do

        ! canopy

        do j = center(2) - Ntree_height / 2, center(2) + Ntree_height / 2
          do i = center(1) - Ntree_height / 2, center(1) + Ntree_height / 2
            do k = C1%glob_zm - Ntree_height, C1%glob_zm ! lowermost layer
              radius = sqrt(real(i - center(1))**2 + real(j - center(2))**2 + 4 * real(k - center(3))**2)
              if (radius .le. real(Ntree_height, ireals) / 3._ireals .and. have_box(k, i, j)) then
                li = i - solver%C_one%xs ! F_idx
                lj = j - solver%C_one%ys ! F_idx
                lk = C1%glob_zm - k + 1 ! F_idx starting at bot

                LAI = LAI_leaf

                tree_tau_thermal(lk, li, lj) = tree_tau_thermal(lk, li, lj) + LAI
              end if
            end do
          end do
        end do
      end associate
    end subroutine

    logical function have_box(gk, gi, gj)
      integer(iintegers), intent(in) :: gk, gi, gj ! global inidces
      have_box = all([ &
        & is_inrange(gk, solver%C_one%zs + 1, solver%C_one%ze + 1), &
        & is_inrange(gi, solver%C_one%xs + 1, solver%C_one%xe + 1), &
        & is_inrange(gj, solver%C_one%ys + 1, solver%C_one%ye + 1)])
    end function

  end subroutine
end module
