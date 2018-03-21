!-------------------------------------------------------------------------
! This file is part of the tenstream solver.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright (C) 2010-2015  Fabian Jakub, <fabian@jakub.com>
!-------------------------------------------------------------------------

module m_f2c_tenstream

      use iso_c_binding

      use m_data_parameters, only : init_mpi_data_parameters, &
        iintegers, ireals, mpiint, default_str_len, &
        zero, one

      use m_tenstr_rrtmg, only : tenstream_rrtmg, destroy_tenstream_rrtmg

      use m_tenstream_options, only: read_commandline_options

      use m_helper_functions, only: imp_bcast, mean, CHKERR

      use m_pprts, only : init_pprts, t_solver, t_solver_8_10, t_solver_3_10, t_solver_3_6, t_solver_1_2, &
        set_global_optical_properties, solve_pprts, destroy_pprts, &
        pprts_get_result_toZero, t_coord

      use m_petsc_helpers, only : petscVecToF90, petscGlobalVecToZero, f90VecToPetsc

#include "petsc/finclude/petsc.h"
      use petsc
      implicit none

      private
      public :: f2c_tenstream_rrtmg, f2c_destroy_tenstream_rrtmg, &
        pprts_f2c_init, pprts_f2c_set_global_optical_properties, &
        pprts_f2c_solve, pprts_f2c_get_result, pprts_f2c_destroy

      class(t_solver), allocatable :: solver

#include "f2c_solver_ids.h"

contains

  subroutine f2c_tenstream_rrtmg(comm, Nz, Nx, Ny, dx, dy, &
      phi0, theta0, albedo_thermal, albedo_solar,          &
      c_atm_filename, c_lthermal, c_lsolar,                &
      Nz_merged, cptr_edir, cptr_edn, cptr_eup, cptr_abso, &
      d_plev, d_tlev, d_lwc, d_reliq, d_iwc, d_reice,      &
      nprocx, nxproc, nprocy, nyproc) bind(C)

    integer(c_int), value :: comm
    integer(c_int), intent(in) :: Nx, Ny, Nz                       ! local subdomain size
    real(c_double), intent(in) :: dx, dy                           ! horizontal grid spacing in [m]
    real(c_double), intent(in) :: phi0, theta0                     ! Sun's angles, azimuth phi(0=North, 90=East), zenith(0 high sun, 80=low sun)
    real(c_double), intent(in) :: albedo_thermal, albedo_solar     ! broadband ground albedo for solar and thermal spectrum

    ! Filename of background atmosphere file. ASCII file with columns:
    ! z(km)  p(hPa)  T(K)  air(cm-3)  o3(cm-3) o2(cm-3) h2o(cm-3)  co2(cm-3) no2(cm-3)
    character(kind=c_char,len=1), intent(in) :: c_atm_filename(*)

    ! Compute solar or thermal radiative transfer. Or compute both at once.
    integer(c_int), intent(in) :: c_lthermal, c_lsolar

    integer(c_int), intent(out) :: Nz_merged                       ! will determine the number of layers of the result
    type(c_ptr), intent(out) :: cptr_edir, cptr_edn, cptr_eup      ! fluxes edir, edn, eup have shape(Nz_merged+1, Nx, Ny)
    type(c_ptr), intent(out) :: cptr_abso                          ! abso just (Nz_merged, Nx, Ny)

    real(c_double), dimension(Nz+1, Nx, Ny), intent(in) :: d_plev  ! pressure on layer interfaces    [hPa]
    real(c_double), dimension(Nz+1, Nx, Ny), intent(in) :: d_tlev  ! Temperature on layer interfaces [K]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_lwc   ! liq water content               [g/kg]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_reliq ! effective radius                [micron]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_iwc   ! ice water content               [g/kg]
    real(c_double), dimension(Nz, Nx, Ny), intent(in)   :: d_reice ! ice effective radius            [micron]

    integer(c_int), intent(in) :: nprocx, nprocy                  ! number of processors in x and y
    integer(c_int), intent(in) :: nxproc(nprocx), nyproc(nprocy)  ! local size of subdomain along x and y


    ! From here on local variables
    real(ireals), allocatable, target, save, dimension(:,:,:) :: edir, edn, eup, abso
    real(c_double), allocatable, target, save, dimension(:,:,:) :: edir_dp, edn_dp, eup_dp, abso_dp

    character(default_str_len) :: atm_filename
    logical :: lthermal, lsolar

    call init_mpi_data_parameters(comm)

    atm_filename = c_to_f_string(c_atm_filename)

    lthermal = c_int_2_logical(c_lthermal)
    lsolar   = c_int_2_logical(c_lsolar)

    call tenstream_rrtmg(comm,            &
      real(dx, kind=ireals),              &
      real(dy,kind=ireals),               &
      real(phi0, kind=ireals),            &
      real(theta0, kind=ireals),          &
      real(albedo_thermal, kind=ireals),  &
      real(albedo_solar, kind=ireals),    &
      atm_filename,                       &
      lthermal, lsolar,                   &
      edir,edn,eup,abso,                  &
      real(d_plev, kind=ireals),          &
      real(d_tlev, kind=ireals),          &
      d_lwc=real(d_lwc, kind=ireals),     &
      d_reliq=real(d_reliq, kind=ireals), &
      d_iwc=real(d_iwc, kind=ireals),     &
      d_reice=real(d_reice, kind=ireals), &
      nxproc=int(nxproc, kind=iintegers), &
      nyproc=int(nyproc, kind=iintegers))

    cptr_edn = C_NULL_PTR
    cptr_eup = C_NULL_PTR
    cptr_abso= C_NULL_PTR
    cptr_edir= C_NULL_PTR

    if(kind(dx).ne.ireals) then
      if(allocated(edn)) then
        if(.not.allocated(edn_dp )) allocate(edn_dp (size(edn ,dim=1), size(edn ,dim=2), size(edn ,dim=3)))
        edn_dp  = real(edn , c_double)
        cptr_edn  = c_loc(edn_dp )
      endif
      if(allocated(eup)) then
        if(.not.allocated(eup_dp )) allocate(eup_dp (size(eup ,dim=1), size(eup ,dim=2), size(eup ,dim=3)))
        eup_dp  = real(eup , c_double)
        cptr_eup  = c_loc(eup_dp )
      endif
      if(allocated(abso)) then
        if(.not.allocated(abso_dp)) allocate(abso_dp(size(abso,dim=1), size(abso,dim=2), size(abso,dim=3)))
        abso_dp = real(abso, c_double)
        cptr_abso = c_loc(abso_dp)
      endif
      if(allocated(edir)) then
        if(.not.allocated(edir_dp)) allocate(edir_dp(size(edir,dim=1), size(edir,dim=2), size(edir,dim=3)))
        edir_dp = real(edir, c_double)
        cptr_edir = c_loc(edir_dp)
      endif
    else
      if(allocated(edn )) cptr_edn  = c_loc(edn )
      if(allocated(eup )) cptr_eup  = c_loc(eup )
      if(allocated(abso)) cptr_abso = c_loc(abso)
      if(allocated(edir)) cptr_edir = c_loc(edir)
    endif

    Nz_merged = int(size(abso,dim=1), c_int)
  end subroutine

  subroutine f2c_destroy_tenstream_rrtmg(c_lfinalizepetsc) bind(C) ! Tidy up the solver
    integer(c_int), intent(in) :: c_lfinalizepetsc  ! determines if we drop the Petsc Environment. If you called petsc initialize in the C program, say False, i.e. 0
    call destroy_tenstream_rrtmg(c_int_2_logical(c_lfinalizepetsc))
  end subroutine

  function c_to_f_string(s) result(str)
    use iso_c_binding
    character(kind=c_char,len=1), intent(in) :: s(*)
    character(len=:), allocatable :: str
    integer i, nchars
    i = 1
    do
      if (s(i) == c_null_char) exit
      i = i + 1
    end do
    nchars = i - 1  ! Exclude null character from Fortran string
    allocate(character(len=nchars) :: str)
    str = transfer(s(1:nchars), str)
  end function c_to_f_string

  logical function c_int_2_logical(i)
    integer(c_int) :: i
    if(i.eq.0) then
      c_int_2_logical = .False.
    else
      c_int_2_logical = .True.
    endif
  end function

  ! initialize tenstream environment
  ! all nodes in communicator have to call this
  ! but only the zeroth node has to have meaningful values for the arguments except the communicator
  ! all but hhl is overwritten on nonzero nodes
  subroutine pprts_f2c_init(comm, solver_id, Nz,Nx,Ny,dx,dy,hhl, phi0, theta0, collapseindex) bind(C)
    integer(c_int), value :: comm
    integer(c_int),intent(inout) :: solver_id, Nx, Ny, Nz
    real(c_double),intent(inout) :: dx,dy
    real(c_float), intent(in),dimension(Nz+1) :: hhl
    real(c_float), intent(inout) :: phi0,theta0
    integer(c_int),intent(inout) :: collapseindex

    integer(iintegers) :: osolver_id, oNx, oNy, oNz, ocollapseindex
    real(ireals) :: odx,ody,ophi0,otheta0
    real(ireals),allocatable :: ohhl(:)

    real(ireals),allocatable :: odz(:)
    integer(iintegers) :: k
    integer(mpiint) :: myid, ierr

    logical,save :: initialized=.False.

    if(initialized) then
      ierr = 0
      select type(solver)
        class is (t_solver_8_10)
          if(solver_id.ne.SOLVER_ID_PPRTS_8_10) ierr = 1
        class is (t_solver_3_10)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_10) ierr = 1
        class is (t_solver_3_6)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_6) ierr = 1
        class is (t_solver_1_2)
          if(solver_id.ne.SOLVER_ID_PPRTS_1_2) ierr = 1
      end select
      if(ierr.ne.0) call CHKERR(ierr, 'seems you changed the solver type id in between calls... you must destroy the solver first before you use a different kind of solver')
      return
    endif

    call init_mpi_data_parameters(comm)
    call read_commandline_options()
    call mpi_comm_rank(comm, myid, ierr); call CHKERR(ierr)

    if(myid.eq.0) then
      osolver_id = int(solver_id, kind=iintegers)
      oNx     = int(Nx, kind=iintegers)
      oNy     = int(Ny, kind=iintegers)
      oNz     = int(Nz, kind=iintegers)
      odx     = real(dx    ,kind=ireals)
      ody     = real(dy    ,kind=ireals)
      ophi0   = real(phi0  ,kind=ireals)
      otheta0 = real(theta0,kind=ireals)
      ocollapseindex = int(collapseindex, kind=iintegers)

      allocate( ohhl(size(hhl)) )
      ohhl = hhl
    endif

    call imp_bcast(comm, osolver_id, 0_mpiint)
    call imp_bcast(comm, oNx    ,0_mpiint)
    call imp_bcast(comm, oNy    ,0_mpiint)
    call imp_bcast(comm, oNz    ,0_mpiint)
    call imp_bcast(comm, odx    ,0_mpiint)
    call imp_bcast(comm, ody    ,0_mpiint)
    call imp_bcast(comm, ophi0  ,0_mpiint)
    call imp_bcast(comm, otheta0,0_mpiint)
    call imp_bcast(comm, ocollapseindex,0_mpiint)

    call imp_bcast(comm, ohhl,0_mpiint)

    ! and overwrite input values to propagate values back to caller...
    solver_id = int(osolver_id, kind=c_int)
    Nx     = int(oNx, kind=c_int)
    Ny     = int(oNy, kind=c_int)
    Nz     = int(oNz, kind=c_int)
    dx     = real(odx    , kind=c_double)
    dy     = real(ody    , kind=c_double)
    phi0   = real(ophi0  , kind=c_float)
    theta0 = real(otheta0, kind=c_float)
    collapseindex = int(ocollapseindex, kind=c_int)

    ! Now every process has the correct values
    print *,myid,'Initializing Tenstream environment from C Language :: domainshape',solver_id
    print *,myid,'Initializing Tenstream environment from C Language :: domainshape',oNx,oNy,oNz,'::',shape(ohhl)

    allocate(odz(oNz))
    do k=1,Nz
      odz(k) = ohhl(k) - ohhl(k+1)
    enddo

    select case(solver_id)
    case(SOLVER_ID_PPRTS_8_10)
      allocate(t_solver_8_10::solver)
    case(SOLVER_ID_PPRTS_3_10)
      allocate(t_solver_3_10::solver)
    case(SOLVER_ID_PPRTS_3_6)
      allocate(t_solver_3_6::solver)
    case(SOLVER_ID_PPRTS_1_2)
      allocate(t_solver_1_2::solver)
    end select

    call init_pprts(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, solver, dz1d=odz)

    initialized=.True.
  end subroutine

  subroutine pprts_f2c_set_global_optical_properties(Nz, Nx, Ny, albedo, kabs, ksca, g, planck) bind(c)
    integer(c_int), value :: Nx,Ny,Nz
    real(c_float),intent(in) :: albedo
    real(c_float),intent(in),dimension(Nz  ,Nx,Ny) :: kabs, ksca, g
    real(c_float),intent(in),dimension(Nz+1,Nx,Ny) :: planck

    real(ireals) :: oalbedo
    real(ireals),allocatable,dimension(:,:,:) :: okabs, oksca, og, oplanck

    if(solver%myid.eq.0) then
      oalbedo = real(albedo, kind=ireals)
      allocate( okabs  (Nz  ,Nx,Ny) ); okabs   = real(kabs, ireals)
      allocate( oksca  (Nz  ,Nx,Ny) ); oksca   = real(ksca, ireals)
      allocate( og     (Nz  ,Nx,Ny) ); og      = real(g, ireals)
      allocate( oplanck(Nz+1,Nx,Ny) ); oplanck = real(planck, ireals)

      if(any(oplanck.gt.zero)) then
        call set_global_optical_properties(solver, oalbedo, okabs, oksca, og, oplanck)
      else
        call set_global_optical_properties(solver, oalbedo, okabs, oksca, og)
      endif

      print *,'mean kabs  ',sum(okabs)  /size(okabs)
      print *,'mean ksca  ',sum(oksca)  /size(oksca)
      print *,'mean g     ',sum(og)     /size(og)
      print *,'mean planck',sum(oplanck)/size(oplanck)
    else !slave
      call set_global_optical_properties(solver)
    endif
  end subroutine

  subroutine pprts_f2c_solve(comm, edirTOA) bind(c)
    ! solve tenstream equations
    ! optical properties have had to be set and environment had to be initialized
    ! incoming solar radiation need only be set by zeroth node
    integer(c_int), value :: comm
    real(c_float), value :: edirTOA
    real(ireals) :: oedirTOA

    if(solver%myid.eq.0) oedirTOA = real(edirTOA, ireals)
    call imp_bcast(comm, oedirTOA, 0_mpiint)

    call solve_pprts(solver, oedirTOA)
  end subroutine

  subroutine pprts_f2c_destroy() bind(c)
    call destroy_pprts(solver, lfinalizepetsc=.False.)
    deallocate(solver)
  end subroutine

  subroutine pprts_f2c_get_result(Nz,Nx,Ny, res_edn, res_eup, res_abso, res_edir) bind(c)
    ! after solving equations -- retrieve the results for edir,edn,eup and absorption
    ! only zeroth node gets the results back.

    integer(c_int), value :: Nx,Ny,Nz
    real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edn
    real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_eup
    real(c_float),intent(out),dimension(Nz  ,Nx,Ny) :: res_abso
    real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edir

    real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso

    call pprts_get_result_toZero(solver,redn,reup,rabso,redir)

    if(solver%myid.eq.0) then
      res_edn  = real(redn (:, 1:Nx, 1:Ny), kind=c_float)
      res_eup  = real(reup (:, 1:Nx, 1:Ny), kind=c_float)
      res_abso = real(rabso(:, 1:Nx, 1:Ny), kind=c_float)
      res_edir = real(redir(:, 1:Nx, 1:Ny), kind=c_float)
      print *,'pprts_f2c_get_result result_edir first column', res_edir(:,1,1)
      print *,'pprts_f2c_get_result redir first column', redir(:,1,1)
    endif

  end subroutine
end module
