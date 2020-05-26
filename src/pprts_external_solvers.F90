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

module m_pprts_external_solvers
#include "petsc/finclude/petsc.h"
  use petsc

  use m_data_parameters, only : ireals, iintegers, mpiint, i0, i1, zero
  use m_helper_functions, only : CHKERR, CHKWARN, itoa
  use m_petsc_helpers, only : getVecPointer, restoreVecPointer

  use m_pprts_base, only : t_solver, t_state_container, atmk
  use m_schwarzschild, only: schwarzschild, B_eff
  use m_twostream, only: delta_eddington_twostream, adding_delta_eddington_twostream

  implicit none
  private
  public :: twostream, schwarz

  logical,parameter :: ldebug=.False.

contains

  !> @brief wrapper for the delta-eddington twostream solver
  !> @details solve the radiative transfer equation for infinite horizontal slabs
  subroutine twostream(solver, edirTOA, solution)
    class(t_solver), intent(inout) :: solver
    real(ireals),intent(in)       :: edirTOA
    type(t_state_container)       :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: xv_dir=>null(),xv_diff=>null()
    real(ireals),pointer,dimension(:) :: xv_dir1d=>null(),xv_diff1d=>null()
    integer(iintegers) :: i,j,src

    real(ireals),allocatable :: dtau(:),kext(:),w0(:),g(:),S(:),Edn(:),Eup(:)
    real(ireals) :: mu0,incSolar,fac
    integer(mpiint) :: ierr

    associate(atm         => solver%atm, &
        C_diff      => solver%C_diff, &
        C_dir       => solver%C_dir, &
        C_one_atm   => solver%C_one_atm, &
        C_one_atm1  => solver%C_one_atm1)

      if(solution%lsolar_rad) then
        call PetscObjectSetName(solution%edir,'twostream_edir_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
        call VecSet(solution%edir ,zero,ierr); call CHKERR(ierr)
      endif

      call PetscObjectSetName(solution%ediff,'twostream_ediff_vec uid='//itoa(solution%uid),ierr) ; call CHKERR(ierr)
      call VecSet(solution%ediff,zero,ierr); call CHKERR(ierr)

      allocate( dtau(C_one_atm%zm) )
      allocate( kext(C_one_atm%zm) )
      allocate(   w0(C_one_atm%zm) )
      allocate(    g(C_one_atm%zm) )

      if(solution%lsolar_rad) &
        call getVecPointer(solution%edir  ,C_dir%da  ,xv_dir1d , xv_dir)
      call getVecPointer(solution%ediff ,C_diff%da ,xv_diff1d, xv_diff)

      allocate( S  (C_one_atm1%zs:C_one_atm1%ze) )
      allocate( Eup(C_one_atm1%zs:C_one_atm1%ze) )
      allocate( Edn(C_one_atm1%zs:C_one_atm1%ze) )

      do j=C_one_atm%ys,C_one_atm%ye
        do i=C_one_atm%xs,C_one_atm%xe

          mu0 = solver%sun%costheta(C_one_atm1%zs,i,j)
          incSolar = edirTOA* mu0

          kext = atm%kabs(:,i,j) + atm%ksca(:,i,j)
          dtau = atm%dz(:,i,j)* kext
          w0   = atm%ksca(:,i,j) / max(kext, epsilon(kext))
          g    = atm%g(:,i,j)

          if(allocated(atm%planck) ) then
            if(allocated(atm%Bsrfc)) then
              call delta_eddington_twostream(dtau, w0, g, &
                mu0, incSolar, atm%albedo(i,j), &
                S, Edn, Eup, &
                planck=atm%planck(:,i,j), &
                planck_srfc=atm%Bsrfc(i,j))
            else
              call delta_eddington_twostream(dtau, w0, g, &
                mu0, incSolar, atm%albedo(i,j), &
                S, Edn, Eup, &
                planck=atm%planck(:,i,j) )
            endif
          else
            !
            ! call adding_delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup )
            !
            !TODO investigate if this one is really ok...
            ! I recently had valgrind errors in VecNorm after calling this:
            ! make -j ex_pprts_rrtm_lw_sw &&
            ! mpirun -np 1 -wdir ../examples/rrtm_lw_sw/ valgrind $(pwd)/bin/ex_pprts_rrtm_lw_sw -twostr_only
            call delta_eddington_twostream(dtau,w0,g,mu0,incSolar,atm%albedo(i,j), S,Edn,Eup )
          endif

          if(solution%lsolar_rad) then
            fac = real(solver%dirtop%area_divider, ireals) / real(solver%dirtop%streams, ireals)
            do src=i0,solver%dirtop%dof-1
              xv_dir(src,C_dir%zs+1:C_dir%ze,i,j) = S(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
              xv_dir(src,C_dir%zs           ,i,j) = S(C_one_atm1%zs) * fac
            enddo
          endif

          fac = real(solver%difftop%area_divider, ireals) / real(solver%difftop%streams, ireals)
          do src = 1, solver%difftop%dof
            if(solver%difftop%is_inward(src)) then
              xv_diff(src-1,C_diff%zs+1:C_diff%ze,i,j) = Edn(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
              xv_diff(src-1,C_diff%zs            ,i,j) = Edn(C_one_atm1%zs) * fac
            else
              xv_diff(src-1,C_diff%zs+1:C_diff%ze,i,j) = Eup(atmk(atm, C_one_atm1%zs)+1:C_one_atm1%ze) * fac
              xv_diff(src-1,C_diff%zs            ,i,j) = Eup(C_one_atm1%zs) * fac
            endif
          enddo
        enddo
      enddo

      if(solution%lsolar_rad) &
        call restoreVecPointer(solution%edir, xv_dir1d, xv_dir  )
      call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

      !Twostream solver returns fluxes as [W]
      solution%lWm2_dir  = .True.
      solution%lWm2_diff = .True.
      ! and mark solution that it is not up to date
      solution%lchanged  = .True.

      deallocate(S)
      deallocate(Edn)
      deallocate(Eup)

    end associate
  end subroutine

  !> @brief simple schwarzschild solver
  !> @details Wrapper for the schwarzschild solver for the radiative transfer equation
  !> \n The solver neglects the scattering term and just solves for lambert beerschen transport + emission
  !> \n This is the simplest radiation solver but quite accurate for thermal calculations
  subroutine schwarz(solver, solution)
    class(t_solver)         :: solver
    type(t_state_container) :: solution

    real(ireals),pointer,dimension(:,:,:,:) :: xv_diff=>null()
    real(ireals),pointer,dimension(:)       :: xv_diff1d=>null()
    integer(iintegers) :: i,j,k,idof
    integer(iintegers) :: Nmu, ak
    logical :: lflg

    real(ireals),allocatable :: dtau(:),Edn(:),Eup(:)
    integer(mpiint) :: ierr


    associate( &
        C_diff => solver%C_diff, &
        atm    => solver%atm,    &
        C_one  => solver%C_one,  &
        C_one1 => solver%C_one1)

      if(solution%lsolar_rad) call CHKERR(1_mpiint, 'Tried calling schwarschild solver for solar calculation -- stopping!')
      if(.not.allocated(atm%planck)) call CHKERR(1_mpiint, 'Tried calling schwarschild solver but no planck was given -- stopping!')

      call VecSet(solution%ediff, zero, ierr); call CHKERR(ierr)

      allocate(dtau(size(atm%dz,dim=1)))

      call getVecPointer(solution%ediff, C_diff%da, xv_diff1d, xv_diff)

      allocate( Eup(0:size(atm%dz,dim=1)) )
      allocate( Edn(0:size(atm%dz,dim=1)) )

      Nmu = 10
      call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER , &
        "-schwarzschild_Nmu" , Nmu, lflg , ierr) ;call CHKERR(ierr)

      if(solver%myid.eq.0 .and. ldebug) print *,' CALCULATING schwarzschild ::'

      do j = C_diff%ys, C_diff%ye
        do i = C_diff%xs, C_diff%xe

          dtau = atm%dz(:, i, j) * atm%kabs(:, i, j)

          if(allocated(atm%Bsrfc)) then
            call schwarzschild(Nmu, dtau, atm%albedo(i,j), Edn, Eup, &
              atm%planck(:, i, j), opt_srfc_emission=atm%Bsrfc(i,j))
          else
            call schwarzschild(Nmu, dtau, atm%albedo(i,j), Edn, Eup, &
              atm%planck(:, i, j))
          endif

          ! icollapse needs special case for TOA flx's
          do idof = 0, solver%difftop%dof-1
            if (solver%difftop%is_inward(i1+idof)) then ! Edn
              xv_diff(idof,C_diff%zs,i,j) = Edn(0)
            else ! Eup
              xv_diff(idof,C_diff%zs,i,j) = Eup(0)
            endif
          enddo

          ! rest of the atmosphere
          do k=C_diff%zs+1,C_diff%ze
            ak = atmk(atm,k)
            do idof = 0, solver%difftop%dof-1
              if (solver%difftop%is_inward(i1+idof)) then ! Edn
                xv_diff(idof,k,i,j) = Edn(ak)
              else ! Eup
                xv_diff(idof,k,i,j) = Eup(ak)
              endif
            enddo
          enddo
        enddo
      enddo
      xv_diff = xv_diff / real(solver%difftop%streams, ireals)

      call restoreVecPointer(solution%ediff, xv_diff1d, xv_diff )

      !Schwarzschild solver returns fluxes as [W/m^2]
      solution%lWm2_dir  = .True.
      solution%lWm2_diff = .True.
      ! and mark solution that it is not up to date
      solution%lchanged         = .True.

      deallocate(Edn)
      deallocate(Eup)

    end associate
  end subroutine

end module
