module test_pprts_slope_correction

  use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, irealLUT, zero, one,  pi, pi64, mpiint, default_str_len, i1

#include "petsc/finclude/petsc.h"
  use petsc

  use m_tenstream_options, only: read_commandline_options
  use m_helper_functions, only: itoa, linspace, CHKERR, deg2rad, spherical_2_cartesian, cross_3d
  use pfunit_mod
  use m_pprts_base, only : t_solver, allocate_pprts_solver_from_commandline

  use m_dyn_atm_to_rrtmg, only: t_tenstr_atm, setup_tenstr_atm, &
    destroy_tenstr_atm, print_tenstr_atm

  use m_pprts_rrtmg, only : pprts_rrtmg, destroy_pprts_rrtmg


  implicit none


contains
  @before
  subroutine setup(this)
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine setup

  @after
  subroutine teardown(this)
    class (MpiTestMethod), intent(inout) :: this
    continue
  end subroutine teardown

  real(ireals) function hill_pressure_deficiency(jglob, ny_glob, hill_dP, hill_shape) result(dP)
    integer(iintegers), intent(in) :: jglob, ny_glob
    real(ireals), intent(in) :: hill_dP, hill_shape
    dP = hill_dP / ( 1._ireals + ((real(jglob, ireals)-(real(ny_glob, ireals)+1._ireals)/2._ireals)/hill_shape)**2 )
  end function

  @test(npes =[1])
  subroutine test_pprts_slope_correction_ex1(this)
    class (MpiTestMethod), intent(inout) :: this

    integer(iintegers) :: j, jglob, k, x_index, y_index, cols
    integer(iintegers),parameter :: nxp=5,nyp=5,nzp=5
    real(ireals),parameter :: dx=100,dy=dx
    real(ireals) :: phi0, theta0
    real(ireals),parameter :: atolerance=.1
    real(ireals), parameter :: two=2

    real(ireals),parameter :: albedo_th=0, albedo_sol=.3 ! broadband ground albedo for solar and thermal spectrum
    logical :: lflg
    real(ireals), dimension(nzp+1,nxp,nyp), target :: plev ! pressure on layer interfaces [hPa]
    real(ireals), dimension(nzp+1,nxp,nyp), target :: tlev ! Temperature on layer interfaces [K]

    real(ireals), pointer, dimension(:,:) :: pplev, ptlev
    integer(iintegers), allocatable :: nxproc(:), nyproc(:)

    real(ireals),allocatable, dimension(:,:,:) :: edir, edn, eup, abso ! [nlev_merged(-1), nxp, nyp]
    ! static: save statt parameter
    character(len=default_str_len),parameter :: atm_filename='atmosphere.dat'

    logical :: lthermal, lsolar

    integer(mpiint) :: numnodes, comm, myid, N_ranks_x, N_ranks_y, ierr
    class(t_solver), allocatable :: pprts_solver
    type(t_tenstr_atm) :: atm
    real(ireals), allocatable :: hill(:), nnormal( :, : )
    real(ireals), allocatable :: edir_target(:, :)
    real(ireals), parameter :: E0=1
    real(ireals) :: hill_dP, hill_shape, dP, sundir(3)
    real(ireals) :: gradient_x, gradient_y, normal(3)

    comm     = this%getMpiCommunicator()
    numnodes = this%getNumProcesses()
    myid     = this%getProcessRank()

    N_ranks_x = 1
    N_ranks_y = numnodes

    allocate(nxproc(N_ranks_x), source=nxp) ! dimension will determine how many ranks are used along the axis
    allocate(nyproc(N_ranks_y), source=nyp) ! values have to define the local domain sizes on each rank (here constant on all processes)

    ! Have to call init_mpi_data_parameters() to define datatypes
    call init_mpi_data_parameters(comm)


    hill_dP = 100 ! [hPa]
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-hill_dP", hill_dP, lflg, ierr); call CHKERR(ierr)
    hill_shape = 3
    call PetscOptionsGetReal(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-hill_shape", hill_shape, lflg, ierr); call CHKERR(ierr) ! the bigger the flatter

    do j=1,nyp
      jglob = j + nyp*myid
      dp = hill_pressure_deficiency(jglob, nyp*numnodes, hill_dP, hill_shape) &
          -hill_pressure_deficiency(   i1, nyp*numnodes, hill_dP, hill_shape)
      dP = max(0._ireals, dP)

      do k=1,nzp+1
        plev(k,:,j) = linspace(k, [1e3_ireals - dP, 500._ireals], nzp+1)
        tlev(k,:,j) = 250 !linspace(k, [290._ireals, 250._ireals], nzp+1)
      enddo
      !print *,'plev@j',j,jglob,':',plev(1,1,j)
    enddo

    pplev(1:size(plev,1),1:size(plev,2)*size(plev,3)) => plev
    ptlev(1:size(tlev,1),1:size(tlev,2)*size(tlev,3)) => tlev

    lsolar = .True.
    lthermal = .False.

    call allocate_pprts_solver_from_commandline(pprts_solver, '3_10', ierr); call CHKERR(ierr)

    call setup_tenstr_atm(comm, .False., atm_filename, &
      pplev, ptlev, atm)

    !call print_tenstr_atm(atm)

    allocate( edir_target( nxp, nyp ) )
    cols = size( atm%zt, 2 )
    allocate( hill( cols ) )
    hill = ( maxval( atm%zt ) - atm%zt( ubound( atm%zt, 1 ), : ) )
    allocate( nnormal( cols, 3 ) )

    do j=1, cols
      if ( j < nxp  + 1 .or. j > cols - nxp ) then

        if ( j == 1 .or. j == cols - nxp + 1 ) then
          gradient_x = ( hill( j + 1 ) - hill( j ) ) / ( two * dx )
        else if ( j == nxp .or. j == cols ) then
          gradient_x = ( hill( j ) - hill( j - 1 ) ) / ( two * dx )
        else
          gradient_x = ( hill( j + 1 ) - hill( j - 1 ) ) / ( two * dx )
        end if

        if ( j > cols - nxp ) then
          gradient_y = ( hill( j ) - hill( j - nxp ) ) / ( two * dy )
        else
          gradient_y = ( hill( j + nxp ) - hill( j ) ) / ( two * dy )
        end if

      else
        if ( mod( j, nxp ) == 0 ) then
          gradient_x = ( hill( j ) - hill( j - 1 ) ) / ( two * dx )
        else if ( mod( j - 1, nxp ) == 0 ) then
          gradient_x = ( hill( j + 1 ) - hill( j ) ) / ( two * dx )
        else
          gradient_x = ( hill( j + 1)  - hill( j - 1) ) / ( two * dx )
        end if
        gradient_y = ( hill( j + nxp ) - hill( j - nxp ) ) / ( two * dy )

      end if

      normal = cross_3d( [ dx, zero, dx * gradient_x ], [ zero, dy, dy * gradient_y ] )
      nnormal( j, : ) = normal / norm2(normal)
    enddo

    do k=1, 4
      theta0 = 30._ireals
      phi0 = 90._ireals * k

      print *, '________________________________________________________'
      print *, 'k', k, 'theta0', theta0, 'phi0', phi0

      sundir = spherical_2_cartesian( phi0, theta0 )

      call pprts_rrtmg(comm, pprts_solver, atm, nxp, nyp, &
        dx, dy, sundir,         &
        albedo_th, albedo_sol,  &
        lthermal, lsolar,       &
        edir, edn, eup, abso,   &
        nxproc=nxproc, nyproc=nyproc, opt_solar_constant=E0)

      print *, 'sundir', sundir
      do j=1, cols
        if ( mod( j, nxp ) == 0 ) then
          x_index = nxp
        else
          x_index = mod( j, nxp )
        end if

        if ( ( j - 1) / nxp  == 0 ) then
          y_index = 1
        else
          y_index =  ( j - 1 ) / nxp + 1
        end if
        edir_target( x_index, y_index ) = E0 * cos(deg2rad(theta0)) * dot_product( - sundir, nnormal( j, : ) ) / dot_product( - sundir, [ zero, zero, one ] )
      enddo

      do j=1, nyp
        print *, 'edir_target', j, edir_target(1, j)
      enddo

      do j=1, nyp
        print *, 'edir', j, edir(ubound(edir, 1), 1, j)
      enddo

      @assertEqual(edir_target, edir(ubound(edir, 1), :, :), atolerance, 'Should be equal')

    enddo

    end subroutine
end module
