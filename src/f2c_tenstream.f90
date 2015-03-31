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

module f2c_tenstream

      use iso_c_binding

      use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,imp_comm,myid,mpierr,zero

      use m_tenstream, only : init_tenstream, set_global_optical_properties, solve_tenstream, destroy_tenstream,&
!                            edir,ediff,abso, &
                            tenstream_get_result, getvecpointer, restorevecpointer, &
                            t_coord,C_dir,C_diff,C_one

      use m_tenstream_options, only: read_commandline_options

      use m_helper_functions, only: imp_bcast,mean

      implicit none
#include "finclude/petsc.h90"

      PetscErrorCode :: ierr

contains

      subroutine tenstr_f2c_init(comm, Nx,Ny,Nz,dx,dy,hhl, phi0, theta0, albedo ) bind(C)                                
        ! initialize tenstream environment
        ! all nodes in communicator have to call this 
        ! but only the zeroth node has to have meaningful values for the arguments except the communicator
        ! all but hhl is overwritten on nonzero nodes

        integer(c_int), value :: comm
        integer(c_int),intent(inout) :: Nx,Ny,Nz                     
        real(c_double),intent(inout) :: dx,dy 
        real(c_float), intent(inout) :: phi0,theta0,albedo
        real(c_float), intent(in),dimension(Nz+1) :: hhl

        integer(iintegers) :: oNx,oNy,oNz
        real(ireals) :: odx,ody,ophi0,otheta0,oalbedo
        real(ireals),allocatable :: ohhl(:)

        real(ireals),allocatable :: odz(:)
        integer(iintegers) :: k

        logical,save :: initialized=.False.

        if(initialized) return

        call init_mpi_data_parameters(comm)
        call read_commandline_options()

        if(myid.eq.0) then

          oNx     = Nx
          oNy     = Ny
          oNz     = Nz
          odx     = dx
          ody     = dy
          ophi0   = phi0
          otheta0 = theta0
          oalbedo = albedo

          allocate( ohhl(size(hhl)) )
          ohhl = hhl
        endif

        call imp_bcast(oNx    ,0_mpiint,myid)
        call imp_bcast(oNy    ,0_mpiint,myid)
        call imp_bcast(oNz    ,0_mpiint,myid)
        call imp_bcast(odx    ,0_mpiint,myid)
        call imp_bcast(ody    ,0_mpiint,myid)
        call imp_bcast(ophi0  ,0_mpiint,myid)
        call imp_bcast(otheta0,0_mpiint,myid)
        call imp_bcast(oalbedo,0_mpiint,myid)

        call imp_bcast(ohhl,0_mpiint,myid)

        Nx     = oNx
        Ny     = oNy
        Nz     = oNz
        dx     = odx
        dy     = ody
        phi0   = ophi0
        theta0 = otheta0
        albedo = oalbedo

        ! Now every process has the correct values
        !        print *,myid,'Initializing Tenstream environment from C Language :: domainshape',oNx,oNy,oNz,'::',shape(ohhl)

        allocate(odz(oNz))
        do k=1,Nz
          odz(k) = ohhl(k) - ohhl(k+1)
        enddo

        call init_tenstream(imp_comm, oNx,oNy,oNz, odx,ody, ophi0, otheta0, oalbedo, dz1d=odz)

        initialized=.True.
      end subroutine                                             

      subroutine tenstr_f2c_set_global_optical_properties(Nx,Ny,Nz, kabs, ksca, g, planck) bind(c)
        integer(c_int), value :: Nx,Ny,Nz
        real(c_float),intent(in),dimension(Nx,Ny,Nz) :: kabs, ksca, g
        real(c_float),intent(in),dimension(Nx,Ny,Nz+1) :: planck

        real(ireals),allocatable,dimension(:,:,:) :: okabs, oksca, og, oplanck
        
        if(myid.eq.0) then
          allocate( okabs  (Nx,Ny,Nz) ); okabs   = kabs
          allocate( oksca  (Nx,Ny,Nz) ); oksca   = ksca
          allocate( og     (Nx,Ny,Nz) ); og      = g   
          allocate( oplanck(Nx,Ny,Nz+1) ); oplanck = planck

          if(any(oplanck.gt.zero)) then
            call set_global_optical_properties(okabs, oksca, og, oplanck)
          else
            call set_global_optical_properties(okabs, oksca, og)
          endif

!          print *,'mean kabs  ',sum(okabs)  /size(okabs)
!          print *,'mean ksca  ',sum(oksca)  /size(oksca)
!          print *,'mean g     ',sum(og)     /size(og)
!          print *,'mean planck',sum(oplanck)/size(oplanck)
        else !slave
          call set_global_optical_properties()
        endif
      end subroutine

      subroutine tenstr_f2c_solve(edirTOA) bind(c)
        ! solve tenstream equations
        ! optical properties have had to be set and environment had to be initialized
        ! incoming solar radiation need only be set by zeroth node
        real(c_float), value :: edirTOA
        real(ireals) :: oedirTOA

        if(myid.eq.0) oedirTOA = edirTOA
        call imp_bcast(oedirTOA    ,0_mpiint,myid)

        call solve_tenstream(oedirTOA)
      end subroutine

      subroutine tenstr_f2c_destroy() bind(c)
        call destroy_tenstream(lfinalizepetsc=.False.)
      end subroutine

      subroutine tenstr_f2c_get_result(Nx,Ny,Nz, res_edir,res_edn,res_eup,res_abso) bind(c)
        ! after solving equations -- retrieve the results for edir,edn,eup and absorption
        ! only zeroth node gets the results back.

        integer(c_int), value :: Nx,Ny,Nz
        real(c_float),intent(out),dimension(Nx,Ny,Nz+1) :: res_edir
        real(c_float),intent(out),dimension(Nx,Ny,Nz+1) :: res_edn
        real(c_float),intent(out),dimension(Nx,Ny,Nz+1) :: res_eup
        real(c_float),intent(out),dimension(Nx,Ny,Nz  ) :: res_abso
        real(ireals),allocatable,dimension(:,:,:,:) :: res
        real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso

        Vec :: vec
        real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso

        PetscScalar,pointer,dimension(:,:,:,:) :: xinp=>null()
        PetscScalar,pointer,dimension(:) :: xinp1d=>null()

        allocate( redir(C_dir%xs :C_dir%xe , C_dir%ys :C_dir%ye  , C_dir%zs :C_dir%ze) )
        allocate( redn (C_diff%xs :C_diff%xe , C_diff%ys :C_diff%ye  , C_diff%zs :C_diff%ze) )
        allocate( reup (C_diff%xs :C_diff%xe , C_diff%ys :C_diff%ye  , C_diff%zs :C_diff%ze) )
        allocate( rabso(C_one%xs :C_one%xe , C_one%ys :C_one%ye  , C_one%zs :C_one%ze) )

        call tenstream_get_result(redir,redn,reup,rabso)

        call DMCreateGlobalVector(C_diff%da,vec,ierr) ; CHKERRQ(ierr)
        call getVecPointer(vec ,C_diff ,xinp1d, xinp)
        xinp(0,:,:,:) = redir
        xinp(1,:,:,:) = redn
        xinp(2,:,:,:) = reup
        xinp(3,C_one%xs :C_one%xe , C_one%ys :C_one%ye  , C_one%zs :C_one%ze) = rabso
        call restoreVecPointer(vec ,C_diff ,xinp1d, xinp )

        call globalVec2Local(vec,C_diff,res)

        call VecDestroy(vec,ierr); CHKERRQ(ierr)

        if(myid.eq.0) then
          res_edir = res(1,:,:,:)
          res_edn  = res(2,:,:,:)
          res_eup  = res(3,:,:,:)
          res_abso = res(4,:,:,C_one%zs+1 :C_one%ze+1)
        endif

        if(myid.eq.0) then
          print *,'Retrieving results:',Nx,Ny,Nz+1
          print *,sum(res_edir)/size(res_edir)
          print *,sum(res_edn) /size(res_edn)
          print *,sum(res_eup) /size(res_eup)
          print *,sum(res_abso)/size(res_abso),'surf',res_abso(1,1,Nz),'toa',res_abso(1,1,1)
        endif

        if(myid.eq.0) deallocate(res)
     end subroutine

      subroutine globalVec2Local(vec,C,res)
        Vec :: vec
        real(ireals),allocatable :: res(:,:,:,:)
        type(t_coord) :: C
        
        Vec :: natural, local
        VecScatter :: scatter_context

        PetscScalar,Pointer :: xloc(:)

        if(allocated(res)) deallocate(res)
        if(myid.eq.0) allocate( res(C%dof,C%glob_xm,C%glob_ym,C%glob_zm) )

        call DMDACreateNaturalVector(C%da, natural, ierr); CHKERRQ(ierr)

        call DMDAGlobalToNaturalBegin(C%da,vec, INSERT_VALUES, natural, ierr); CHKERRQ(ierr)
        call DMDAGlobalToNaturalEnd  (C%da,vec, INSERT_VALUES, natural, ierr); CHKERRQ(ierr)

        call VecScatterCreateToZero(natural, scatter_context, local, ierr); CHKERRQ(ierr)

        call VecScatterBegin(scatter_context, natural, local, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
        call VecScatterEnd  (scatter_context, natural, local, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)

        call VecScatterDestroy(scatter_context, ierr); CHKERRQ(ierr)

        if(myid.eq.0) then
          call VecGetArrayF90(local,xloc,ierr) ;CHKERRQ(ierr)

          res = reshape( xloc, (/ C%dof,C%glob_xm,C%glob_ym,C%glob_zm /) )

          call VecRestoreArrayF90(local,xloc,ierr) ;CHKERRQ(ierr)
        endif

        call VecDestroy(local,ierr); CHKERRQ(ierr)
        call VecDestroy(natural,ierr); CHKERRQ(ierr)
      end subroutine
        
end module
