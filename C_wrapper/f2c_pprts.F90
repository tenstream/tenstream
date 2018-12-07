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

module m_f2c_pprts

      use iso_c_binding

      use m_data_parameters, only : init_mpi_data_parameters, &
        iintegers, ireals, mpiint, default_str_len, &
        zero, one

      use m_tenstream_options, only: read_commandline_options

      use m_helper_functions, only: imp_bcast, meanval, CHKERR

      use m_pprts_base, only : t_solver, t_solver_3_10, t_solver_3_6, t_solver_1_2, t_coord, &
        t_solver_8_10, t_solver_3_16, t_solver_8_16, t_solver_8_18
      use m_pprts, only : init_pprts, set_global_optical_properties, solve_pprts, destroy_pprts, &
        pprts_get_result_toZero

      use m_petsc_helpers, only : petscVecToF90, petscGlobalVecToZero, f90VecToPetsc

#include "petsc/finclude/petsc.h"
      use petsc
      implicit none

      private
      public :: pprts_f2c_init, pprts_f2c_set_global_optical_properties, &
        pprts_f2c_solve, pprts_f2c_get_result, pprts_f2c_destroy

      class(t_solver), allocatable :: solver

#include "f2c_solver_ids.h"

contains

  ! initialize pprts environment
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
        class is (t_solver_1_2)
          if(solver_id.ne.SOLVER_ID_PPRTS_1_2) ierr = 1
        class is (t_solver_3_6)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_6) ierr = 1
        class is (t_solver_3_10)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_10) ierr = 1
        class is (t_solver_8_10)
          if(solver_id.ne.SOLVER_ID_PPRTS_8_10) ierr = 1
        class is (t_solver_3_16)
          if(solver_id.ne.SOLVER_ID_PPRTS_3_16) ierr = 1
        class is (t_solver_8_16)
          if(solver_id.ne.SOLVER_ID_PPRTS_8_16) ierr = 1
        class is (t_solver_8_18)
          if(solver_id.ne.SOLVER_ID_PPRTS_8_18) ierr = 1
      end select
      if(ierr.ne.0) call CHKERR(ierr, 'seems you changed the solver type id in between calls... you must destroy the solver first before you use a different kind of solver')
      return
    endif

    call init_mpi_data_parameters(comm)
    call read_commandline_options(comm)
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
    print *,myid,'Initializing pprts environment from C Language :: solver_id', solver_id
    print *,myid,'Initializing pprts environment from C Language :: domainshape',oNx,oNy,oNz,'::',shape(ohhl)

    allocate(odz(oNz))
    do k=1,Nz
      odz(k) = ohhl(k) - ohhl(k+1)
    enddo

    select case(solver_id)
    case(SOLVER_ID_PPRTS_1_2)
      allocate(t_solver_1_2::solver)
    case(SOLVER_ID_PPRTS_3_6)
      allocate(t_solver_3_6::solver)
    case(SOLVER_ID_PPRTS_3_10)
      allocate(t_solver_3_10::solver)
    case(SOLVER_ID_PPRTS_8_10)
      allocate(t_solver_8_10::solver)
    case(SOLVER_ID_PPRTS_3_16)
      allocate(t_solver_3_16::solver)
    case(SOLVER_ID_PPRTS_8_16)
      allocate(t_solver_8_16::solver)
    case(SOLVER_ID_PPRTS_8_18)
      allocate(t_solver_8_18::solver)
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

      print *,'mean kabs  ',meanval(okabs)
      print *,'mean ksca  ',meanval(oksca)
      print *,'mean g     ',meanval(og)
      print *,'mean planck',meanval(oplanck)
    else !slave
      call set_global_optical_properties(solver)
    endif
  end subroutine

  subroutine pprts_f2c_solve(comm, edirTOA) bind(c)
    ! solve pprts equations
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
    call destroy_pprts(solver, lfinalizepetsc=.True.)
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
      print *,'pprts_f2c_get_result rabso first column', rabso(:,1,1)
    endif

  end subroutine
end module
