module m_example_pprts_rrtmg_hill

#include "petsc/finclude/petsc.h"
  use petsc
  use mpi

  ! Import datatype from the TenStream lib. Depending on how PETSC is
  ! compiled(single or double floats, or long ints), this will determine what
  ! the Tenstream uses.
  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint, zero, one, i1, default_str_len

  use m_helper_functions, only : reverse, linspace, CHKERR, meanval, itoa, spherical_2_cartesian

  use m_search, only: search_sorted_bisection

  ! Import specific solver type: 3_10 for example uses 3 streams direct, 10 streams for diffuse radiation
  use m_pprts_base, only : t_solver, allocate_pprts_solver_from_commandline
  use m_netcdfIO, only : ncwrite, set_global_attribute
  use m_petsc_helpers, only: getvecpointer, restorevecpointer, petscGlobalVecToZero, petscVecToF90, f90VecToPetsc

  ! main entry point for solver, and desctructor
  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, &
    destroy_tenstr_atm, print_tenstr_atm

  implicit none

contains

  real(ireals) function hill_pressure_deficiency(jglob, ny_glob, hill_dP, hill_shape) result(dP)
    integer(iintegers), intent(in) :: jglob, ny_glob
    real(ireals), intent(in) :: hill_dP, hill_shape
    dP = hill_dP / ( 1._ireals + ((real(jglob, ireals)-(real(ny_glob, ireals)-1._ireals)/2._ireals)/hill_shape)**2 )
  end function

  subroutine example_pprts_rrtmg_hill(comm, nxp, nyp, nzp, dx, dy, nxproc, nyproc, xstart, ystart)
    integer(mpiint), intent(in) :: comm
    integer(iintegers), intent(in) :: nxp, nyp, nzp  ! local domain size for each rank
    integer(iintegers), intent(in) :: nxproc(:), nyproc(:) ! local domain sizes on 2d cartesian decomposition
    integer(iintegers), intent(in) :: xStart, yStart ! start indices of local domains
    real(ireals), intent(in) :: dx, dy               ! horizontal grid spacing in [m]

    ! MPI variables and domain decomposition sizes
    integer(mpiint) :: numnodes, myid, ierr

    real(ireals) :: phi0, theta0 ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(ireals),parameter :: albedo_th=0, albedo_sol=.12 ! broadband ground albedo for solar and thermal spectrum

    real(ireals), dimension(nzp+1,xStart:xStart+nxp-1,yStart:yStart+nyp-1), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1,xStart:xStart+nxp-1,yStart:yStart+nyp-1), target :: tlev ! Temperature on layer interfaces [K]

    ! Layer values for the atmospheric constituents -- those are actually all
    ! optional and if not provided, will be taken from the background profile file (atm_filename)
    ! see interface of `tenstream_rrtmg()` for units
    ! real(ireals), dimension(nzp,nxp,nyp) :: h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr

    ! Liquid water cloud content [g/kg] and effective radius in micron
    real(ireals), dimension(nzp,xStart:xStart+nxp-1,yStart:yStart+nyp-1), target :: lwc, reliq

    ! Fluxes and absorption in [W/m2] and [W/m3] respectively.
    ! Dimensions will probably be bigger than the dynamics grid, i.e. will have
    ! the size of the merged grid. If you only want to use heating rates on the
    ! dynamics grid, use the lower layers, i.e.,
    !   edn(ubound(edn,1)-nlay_dynamics : ubound(edn,1) )
    ! or:
    !   abso(ubound(abso,1)-nlay_dynamics+1 : ubound(abso,1) )
    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(len=default_str_len) :: atm_filename ! ='afglus_100m.dat'

    !------------ Local vars ------------------
    integer(iintegers) :: j, k, nlev, icld(2)

    ! reshape pointer to convert i,j vecs to column vecs
    real(ireals), pointer, dimension(:,:) :: pplev, ptlev, plwc, preliq
    real(ireals), pointer, dimension(:,:,:) :: patmlwc, patmreliq

    logical,parameter :: ldebug=.True.
    logical :: lthermal, lsolar

    class(t_solver), allocatable :: pprts_solver
    type(t_tenstr_atm), target :: atm

    real(ireals) :: hill_dP, hill_shape, dP, cld_lwc

    integer(iintegers) :: cld_width
    real(ireals) :: cld_bot, cld_top
    logical :: lflg
    character(len=default_str_len) :: outpath(2)
    real(ireals), pointer :: hhl(:,:,:,:)=>null(), hhl1d(:)=>null()

    call MPI_COMM_SIZE(comm, numnodes, ierr)
    call MPI_COMM_RANK(comm, myid, ierr)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)

    ! Start with a dynamics grid ranging from 1000 hPa up to 250 hPa and a
    ! Temperature difference of 60K
    hill_dP = 100 ! [hPa]
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-hill_dP", hill_dP, lflg, ierr)
    hill_shape = 10
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-hill_shape", hill_shape, lflg, ierr) ! the bigger the flatter

    do j=yStart,yStart+nyp-1
      dp = hill_pressure_deficiency( j, sum(nyproc), hill_dP, hill_shape) &
          -hill_pressure_deficiency(i1, sum(nyproc), hill_dP, hill_shape)
      dP = max(0._ireals, dP)

      do k=1,nzp+1
        plev(k,:,j) = linspace(k, [1e3_ireals - dP, 500._ireals], nzp+1)
        tlev(k,:,j) = 250 !linspace(k, [290._ireals, 250._ireals], nzp+1)
      enddo
    enddo

    ! define a cloud, with liquid water content and effective radius 10 micron
    lwc = 0
    reliq = 10

    icld = -1
    cld_lwc = 0e-2
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-lwc", cld_lwc, lflg, ierr)
    cld_width = 5
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-cld_width", cld_width, lflg, ierr)
    cld_bot = 750._ireals
    cld_top = 700._ireals
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-cld_bot", cld_bot, lflg, ierr)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-cld_top", cld_top, lflg, ierr)

    do j=yStart,yStart+nyp-1
      if( abs(j - real(sum(nyproc)-1, ireals)/2).le.cld_width) then
        icld(1) = nint(search_sorted_bisection(plev(:,1,j), cld_bot))
        icld(2) = nint(search_sorted_bisection(plev(:,1,j), cld_top))

        lwc  (icld(1):icld(2), :, j) = cld_lwc
      endif
    enddo

    if(myid.eq.0 .and. ldebug) print *,'Setup Atmosphere...'
    pplev(1:size(plev,1),1:size(plev,2)*size(plev,3)) => plev
    ptlev(1:size(tlev,1),1:size(tlev,2)*size(tlev,3)) => tlev
    plwc (1:size(lwc ,1),1:size(lwc ,2)*size(lwc ,3)) => lwc
    preliq(1:size(reliq,1),1:size(reliq,2)*size(reliq,3)) => reliq

    atm_filename='afglus_100m.dat'
    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, '-atm_filename', &
      atm_filename, lflg, ierr); call CHKERR(ierr)

    phi0   = 180
    theta0 = 20
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-phi0", phi0, lflg, ierr); call CHKERR(ierr)
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-theta0", theta0, lflg, ierr); call CHKERR(ierr)

    lsolar = .True.
    lthermal = .True.
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-solar", lsolar, lflg,ierr); call CHKERR(ierr)
    call PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-thermal", lthermal, lflg,ierr); call CHKERR(ierr)

    call allocate_pprts_solver_from_commandline(pprts_solver, '3_10')

    call setup_tenstr_atm(comm, .False., atm_filename, &
      pplev, ptlev, atm, &
      d_lwc=plwc, d_reliq=preliq)

    ! Hack the background atmosphere so that there is no jump in temperatures,
    ! i.e. add the difference at the glue location to all background temps
    atm%tlev(atm%d_ke1+1:size(atm%tlev,1),:) = atm%tlev(atm%d_ke1+1:size(atm%tlev,1),:) &
      + meanval(atm%tlev(atm%d_ke1,:)-atm%tlev(atm%d_ke1+1,:))
    if(myid.eq.0) call print_tenstr_atm(atm)

    call pprts_rrtmg(comm, pprts_solver, atm, nxp, nyp, &
      dx, dy, spherical_2_cartesian(phi0, theta0),   &
      albedo_th, albedo_sol,  &
      lthermal, lsolar,       &
      edir, edn, eup, abso,   &
      nxproc=nxproc, nyproc=nyproc )

    nlev = ubound(edn,1)
    if(myid.eq.0) then
      if(ldebug) then
        do k=1,nlev
          if(allocated(edir)) then
            print *,k,'edir', edir(k,1,1), 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
          else
            print *,k, 'edn', edn(k,1,1), 'eup', eup(k,1,1), abso(min(nlev-1,k),1,1)
          endif
        enddo
      endif

      if(allocated(edir)) &
        print *,'surface :: direct flux', edir(nlev,1,1)
      print *,'surface :: downw flux ', edn (nlev,1,1)
      print *,'surface :: upward fl  ', eup (nlev,1,1)
      print *,'surface :: absorption ', abso(nlev-1,1,1)

      if(allocated(edir)) &
        print *,'TOA :: direct flux', edir(1,1,1)
      print *,'TOA :: downw flux ', edn (1,1,1)
      print *,'TOA :: upward fl  ', eup (1,1,1)
      print *,'TOA :: absorption ', abso(1,1,1)

      if(all(icld.ne.-1)) then
        if(allocated(edir)) then
          print *,'icloud :: direct flux  ', edir(nlev-icld  ,1,1)
          print *,'icloud+1 :: direct flux', edir(nlev-icld+1,1,1)
        endif
        print *,'icloud :: downw flux   ', edn (nlev-icld+1,1,1)
        print *,'icloud :: upward fl    ', eup (nlev-icld  ,1,1)
        print *,'icloud :: absorption   ', abso(nlev-icld  ,1,1)
      endif
    endif

    outpath(1) = 'out_pprts_hill.nc'
    call PetscOptionsGetString(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-out", outpath(1), lflg, ierr)

    associate(&
        & C  => pprts_solver%C_one,     &
        & C1 => pprts_solver%C_one1,    &
        & Ca => pprts_solver%C_one_atm, &
        & Ca1=> pprts_solver%C_one_atm1_box )

      call dump_vec(Ca%da, pprts_solver%atm%dz, 'dz')

      patmlwc  (Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%lwc
      patmreliq(Ca%zs:Ca%ze, Ca%xs:Ca%xe, Ca%ys:Ca%ye) => atm%reliq
      call dump_vec(Ca%da, reverse(patmlwc), 'lwc')
      call dump_vec(Ca%da, reverse(patmreliq), 'reliq')

      if(allocated(edir)) &
        call dump_vec(C1%da, edir, 'edir')
      call dump_vec(C1%da, edn , 'edn')
      call dump_vec(C1%da, eup , 'eup')
      call dump_vec(C%da , abso, 'abso')

      call getVecPointer(pprts_solver%atm%hhl, Ca1%da, hhl1d, hhl)
      call dump_vec(Ca1%da, hhl(0,Ca1%zs:Ca1%ze,Ca1%xs:Ca1%xe,Ca1%ys:Ca1%ye), 'hhl')
      call restoreVecPointer(pprts_solver%atm%hhl, hhl1d, hhl)

      if(myid.eq.0) then
        call set_global_attribute(outpath(1), 'Nx', C%glob_xm, ierr); call CHKERR(ierr)
        call set_global_attribute(outpath(1), 'Ny', C%glob_ym, ierr); call CHKERR(ierr)
        call set_global_attribute(outpath(1), 'Nz', nzp, ierr); call CHKERR(ierr)
        call set_global_attribute(outpath(1), 'dx', dx, ierr); call CHKERR(ierr)
        call set_global_attribute(outpath(1), 'dy', dy, ierr); call CHKERR(ierr)
        call set_global_attribute(outpath(1), 'phi0', phi0, ierr); call CHKERR(ierr)
        call set_global_attribute(outpath(1), 'theta0', theta0, ierr); call CHKERR(ierr)
        call set_global_attribute(outpath(1), 'Ag_solar', albedo_sol, ierr); call CHKERR(ierr)
      endif
    end associate

    ! Tidy up
    call destroy_pprts_rrtmg(pprts_solver, lfinalizepetsc=.True.)
    call destroy_tenstr_atm(atm)
  contains
    subroutine dump_vec(dm, arr, varname)
      type(tDM), intent(in) :: dm
      real(ireals), intent(in) :: arr(:,:,:)
      character(len=*), intent(in) :: varname
      type(tVec) :: gvec, lVec
      real(ireals), allocatable :: larr(:,:,:)
      integer(mpiint) :: ierr

      call DMGetGlobalVector(dm, gvec, ierr); call CHKERR(ierr)
      call f90VecToPetsc(arr, dm, gvec)
      call petscGlobalVecToZero(gvec, dm, lVec)
      if(myid.eq.0) then
        call petscVecToF90(lVec, dm, larr, opt_l_only_on_rank0=.True.)

        outpath(2) = trim(varname)
        call ncwrite(outpath, larr, ierr); call CHKERR(ierr)
      endif
      call VecDestroy(lVec, ierr); call CHKERR(ierr)

      call DMRestoreGlobalVector(dm, gvec, ierr); call CHKERR(ierr)
    end subroutine
  end subroutine
end module

program main
#include "petsc/finclude/petsc.h"
  use petsc
  use mpi
  use m_data_parameters, only : iintegers, mpiint, ireals
  use m_helper_functions, only : domain_decompose_2d_petsc, CHKERR
  use m_example_pprts_rrtmg_hill, only: example_pprts_rrtmg_hill

  implicit none

  integer(mpiint) :: ierr, myid
  integer(iintegers),allocatable :: nxproc(:), nyproc(:)
  integer(iintegers) :: Nx, Ny, Nz, nxp, nyp, xStart, yStart
  real(ireals)       :: dx, dy
  logical :: lflg

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myid, ierr)

  call PetscInitialize(PETSC_NULL_CHARACTER ,ierr); call CHKERR(ierr)

  Nx=3; Ny=32; Nz=10
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nx", Nx, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Ny", Ny, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-Nz", Nz, lflg, ierr); call CHKERR(ierr)

  dx = 500
  dy = dx
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dx", dx, lflg, ierr); call CHKERR(ierr)
  call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-dy", dy, lflg, ierr); call CHKERR(ierr)

  call domain_decompose_2d_petsc(mpi_comm_world, Nx, Ny, &
    & nxp, nyp, xStart, yStart, nxproc, nyproc, ierr); call CHKERR(ierr)

  if (myid.eq.0) then
    print *,'Running rrtm_lw_sw example with grid size:', Nx, Ny, Nz,&
      & '(decomp:',nxproc,nyproc,')'
  endif

  call example_pprts_rrtmg_hill(mpi_comm_world, nxp, nyp, Nz, dx, dy, nxproc, nyproc, xStart, yStart)

  call mpi_finalize(ierr)
end program
