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

      use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,imp_comm,myid,mpierr,zero, i0

      use m_tenstream, only : init_tenstream, set_global_optical_properties, solve_tenstream, destroy_tenstream,&
                            tenstream_get_result, getvecpointer, restorevecpointer, &
                            t_coord,C_dir,C_diff,C_one,C_one_atm, C_one_atm1

      use m_tenstream_options, only: read_commandline_options

      use m_helper_functions, only: imp_bcast,mean, CHKERR

#include "petsc/finclude/petsc.h"
      use petsc

      implicit none

      integer(mpiint) :: ierr

contains

      subroutine tenstr_f2c_init(comm, Nz,Nx,Ny,dx,dy,hhl, phi0, theta0, collapseindex) bind(C)
        ! initialize tenstream environment
        ! all nodes in communicator have to call this
        ! but only the zeroth node has to have meaningful values for the arguments except the communicator
        ! all but hhl is overwritten on nonzero nodes

        integer(c_int), value :: comm
        integer(c_int),intent(inout) :: collapseindex
        integer(c_int),intent(inout) :: Nx,Ny,Nz
        real(c_double),intent(inout) :: dx,dy
        real(c_float), intent(inout) :: phi0,theta0
        real(c_float), intent(in),dimension(Nz+1) :: hhl

        integer(iintegers) :: oNx,oNy,oNz,ocollapseindex
        real(ireals) :: odx,ody,ophi0,otheta0
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
          ocollapseindex = collapseindex

          allocate( ohhl(size(hhl)) )
          ohhl = hhl
        endif

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
        Nx     = oNx
        Ny     = oNy
        Nz     = oNz
        dx     = odx
        dy     = ody
        phi0   = ophi0
        theta0 = otheta0
        collapseindex=ocollapseindex

        ! Now every process has the correct values
        !        print *,myid,'Initializing Tenstream environment from C Language :: domainshape',oNx,oNy,oNz,'::',shape(ohhl)

        allocate(odz(oNz))
        do k=1,Nz
          odz(k) = ohhl(k) - ohhl(k+1)
        enddo

        !call init_tenstream(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, dz1d=odz, collapseindex=ocollapseindex)
        call init_tenstream(comm, oNz,oNx,oNy, odx,ody, ophi0, otheta0, dz1d=odz)

        initialized=.True.
      end subroutine

      subroutine tenstr_f2c_set_global_optical_properties(Nz,Nx,Ny, albedo, kabs, ksca, g, planck) bind(c)
        integer(c_int), value :: Nx,Ny,Nz
        real(c_float),intent(in) :: albedo
        real(c_float),intent(in),dimension(Nz  ,Nx,Ny) :: kabs, ksca, g
        real(c_float),intent(in),dimension(Nz+1,Nx,Ny) :: planck

        real(ireals) :: oalbedo
        real(ireals),allocatable,dimension(:,:,:) :: okabs, oksca, og, oplanck
        
        oalbedo = albedo

        if(myid.eq.0) then
          allocate( okabs  (Nz  ,Nx,Ny) ); okabs   = kabs
          allocate( oksca  (Nz  ,Nx,Ny) ); oksca   = ksca
          allocate( og     (Nz  ,Nx,Ny) ); og      = g   
          allocate( oplanck(Nz+1,Nx,Ny) ); oplanck = planck

          if(any(oplanck.gt.zero)) then
            call set_global_optical_properties(oalbedo, okabs, oksca, og, oplanck)
          else
            call set_global_optical_properties(oalbedo, okabs, oksca, og)
          endif

          print *,'mean kabs  ',sum(okabs)  /size(okabs)
          print *,'mean ksca  ',sum(oksca)  /size(oksca)
          print *,'mean g     ',sum(og)     /size(og)
          print *,'mean planck',sum(oplanck)/size(oplanck)
        else !slave
          call set_global_optical_properties(oalbedo)
        endif
      end subroutine

      subroutine tenstr_f2c_solve(edirTOA) bind(c)
        ! solve tenstream equations
        ! optical properties have had to be set and environment had to be initialized
        ! incoming solar radiation need only be set by zeroth node
        real(c_float), value :: edirTOA
        real(ireals) :: oedirTOA

        if(myid.eq.0) oedirTOA = edirTOA
        call imp_bcast(imp_comm, oedirTOA, 0_mpiint)

        call solve_tenstream(oedirTOA)
      end subroutine

      subroutine tenstr_f2c_destroy() bind(c)
        call destroy_tenstream(lfinalizepetsc=.False.)
      end subroutine

      subroutine tenstr_f2c_get_result(Nz,Nx,Ny, res_edir,res_edn,res_eup,res_abso) bind(c)
        ! after solving equations -- retrieve the results for edir,edn,eup and absorption
        ! only zeroth node gets the results back.

        integer(c_int), value :: Nx,Ny,Nz
        real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edir
        real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_edn
        real(c_float),intent(out),dimension(Nz+1,Nx,Ny) :: res_eup
        real(c_float),intent(out),dimension(Nz  ,Nx,Ny) :: res_abso
        real(ireals),allocatable,dimension(:,:,:) :: res

        real(ireals),allocatable,dimension(:,:,:) :: redir,redn,reup,rabso


        allocate( redir(C_one_atm1%zs:C_one_atm1%ze, C_dir%xs :C_dir%xe , C_dir%ys :C_dir%ye   )); redir=0
        allocate( redn (C_one_atm1%zs:C_one_atm1%ze, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye  )); redn =0
        allocate( reup (C_one_atm1%zs:C_one_atm1%ze, C_diff%xs:C_diff%xe, C_diff%ys:C_diff%ye  )); reup =0
        allocate( rabso(C_one_atm%zs :C_one_atm%ze , C_one%xs :C_one%xe , C_one%ys :C_one%ye   )); rabso=0

        call tenstream_get_result(redir,redn,reup,rabso)

        call exchange_var(C_one_atm1, redir, res_edir)
        call exchange_var(C_one_atm1, redn , res_edn )
        call exchange_var(C_one_atm1, reup , res_eup )
        call exchange_var(C_one_atm , rabso, res_abso)

        if(myid.eq.0) then
          print *,'Retrieving results:',Nx,Ny,Nz+1
          print *,sum(res_edir)/size(res_edir)
          print *,sum(res_edn) /size(res_edn)
          print *,sum(res_eup) /size(res_eup)
          print *,sum(res_abso)/size(res_abso),'surf',res_abso(Nz,1,1),'toa',res_abso(1,1,1)
        endif

        contains
            subroutine exchange_var(C, inp, outp)
                type(t_coord),intent(in) :: C
                real(ireals),intent(in) :: inp(:,:,:) ! local array from get_result
                real(c_float),intent(out) :: outp(:,:,:) ! global sized array on rank 0

                real(ireals),allocatable :: tmp(:,:,:,:)


                Vec :: vec
                PetscScalar,pointer,dimension(:,:,:,:) :: xinp=>null()
                PetscScalar,pointer,dimension(:) :: xinp1d=>null()

                call DMGetGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)
                call getVecPointer(vec ,C ,xinp1d, xinp)
                xinp(i0,:,:,:) = inp
                call restoreVecPointer(vec ,C ,xinp1d, xinp )

                call globalVec2Local(vec,C,tmp)

                call DMRestoreGlobalVector(C%da,vec,ierr) ; call CHKERR(ierr)

                if(myid.eq.0) outp = tmp(lbound(tmp,1), &
                                         lbound(tmp,2):lbound(tmp,2)+size(outp,1)-1,&
                                         lbound(tmp,3):lbound(tmp,3)+size(outp,2)-1,&
                                         lbound(tmp,4):lbound(tmp,4)+size(outp,3)-1)
            end subroutine
     end subroutine

      subroutine globalVec2Local(vec,C,res)
        Vec :: vec
        real(ireals),allocatable :: res(:,:,:,:)
        type(t_coord) :: C
        
        Vec :: natural, local
        VecScatter :: scatter_context

        PetscScalar,Pointer :: xloc(:)

        if(allocated(res)) deallocate(res)
        if(myid.eq.0) allocate( res(C%dof,C%glob_zm,C%glob_xm,C%glob_ym) )

        call DMDACreateNaturalVector(C%da, natural, ierr); call CHKERR(ierr)

        call DMDAGlobalToNaturalBegin(C%da,vec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)
        call DMDAGlobalToNaturalEnd  (C%da,vec, INSERT_VALUES, natural, ierr); call CHKERR(ierr)

        call VecScatterCreateToZero(natural, scatter_context, local, ierr); call CHKERR(ierr)

        call VecScatterBegin(scatter_context, natural, local, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)
        call VecScatterEnd  (scatter_context, natural, local, INSERT_VALUES, SCATTER_FORWARD, ierr); call CHKERR(ierr)

        call VecScatterDestroy(scatter_context, ierr); call CHKERR(ierr)

        if(myid.eq.0) then
          call VecGetArrayF90(local,xloc,ierr) ;call CHKERR(ierr)

          res = reshape( xloc, (/ C%dof,C%glob_zm,C%glob_xm,C%glob_ym /) )

          call VecRestoreArrayF90(local,xloc,ierr) ;call CHKERR(ierr)
        endif

        call VecDestroy(local,ierr); call CHKERR(ierr)
        call VecDestroy(natural,ierr); call CHKERR(ierr)
      end subroutine
        
end module
