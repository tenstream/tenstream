module f2c_tenstream

      use iso_c_binding

      use m_data_parameters, only : init_mpi_data_parameters, iintegers, ireals, mpiint ,imp_comm,myid

      use m_tenstream, only : init_tenstream, set_optical_properties, solve_tenstream, destroy_tenstream,&
                            edir,ediff,abso, &
                            edir_twostr,ediff_twostr,abso_twostr,&
                            t_coord,C_dir,C_diff,C_one

      use m_tenstream_options, only: read_commandline_options

      use m_helper_functions, only: imp_bcast

      implicit none
#include "finclude/petsc.h90"

      PetscErrorCode :: ierr

contains

      subroutine f2c_init_tenstream(comm, Nx,Ny,Nz,dx,dy,hhl, phi0, theta0, albedo ) bind(C)                                
        integer(c_int), value :: comm, Nx,Ny,Nz                     
        real(c_double), value :: dx,dy, phi0,theta0,albedo
        real(c_double),intent(in),dimension(Nz+1) :: hhl

        integer(iintegers) :: oNx,oNy,oNz
        real(ireals) :: odx,ody,ophi0,otheta0,oalbedo
        real(ireals),allocatable :: ohhl(:)

        call init_mpi_data_parameters(comm)
        call read_commandline_options()

        print *,myid,'Initializing Tenstream environment from C Language :: domainshape',Nx,Ny,Nz,'::',shape(hhl)

        oNx = Nx
        oNy = Ny
        oNz = Nz
        odx = dx
        ody = dy
        ophi0  = phi0
        otheta0= theta0
        oalbedo= albedo

        if(myid.eq.0) then
          allocate( ohhl(size(hhl)) )
          ohhl = hhl
        endif

        call imp_bcast(ohhl,0_mpiint,myid)

        call init_tenstream(imp_comm, oNx,oNy,oNz, odx,ody,ohhl ,ophi0, otheta0, oalbedo)
        print *,'Initializing Tenstream environment from C Language ... done'
      end subroutine                                             

      subroutine f2c_set_optical_properties(Nx,Ny,Nz, kabs, ksca, g) bind(c)
        integer(c_int), value :: Nx,Ny,Nz
        real(c_double),intent(in),dimension(Nx,Ny,Nz) :: kabs, ksca, g

        real(ireals),allocatable,dimension(:,:,:) :: okabs, oksca, og
        
        if(myid.eq.0) then
          allocate( okabs(Nx,Ny,Nz) ); okabs = kabs
          allocate( oksca(Nx,Ny,Nz) ); oksca = ksca
          allocate( og   (Nx,Ny,Nz) ); og    = g   
        endif

        call imp_bcast(okabs,0_mpiint,myid)
        call imp_bcast(oksca,0_mpiint,myid)
        call imp_bcast(og   ,0_mpiint,myid)

        call set_optical_properties(okabs, oksca, og)
      end subroutine

      subroutine f2c_solve_tenstream(edirTOA) bind(c)
        real(c_double), value :: edirTOA
        real(ireals) :: oedirTOA

        oedirTOA = edirTOA
        call solve_tenstream(oedirTOA)
      end subroutine

      subroutine f2c_destroy_tenstream() bind(c)
        call destroy_tenstream()
      end subroutine

      subroutine get_result(Nx,Ny,Nz, res_edir,res_edn,res_eup,res_abso) bind(c)
        integer(c_int), value :: Nx,Ny,Nz
        real(c_double),intent(out),dimension(Nx,Ny,Nz+1) :: res_edir
        real(c_double),intent(out),dimension(Nx,Ny,Nz+1) :: res_edn
        real(c_double),intent(out),dimension(Nx,Ny,Nz+1) :: res_eup
        real(c_double),intent(out),dimension(Nx,Ny,Nz  ) :: res_abso
        real(ireals),allocatable,dimension(:,:,:,:) :: res

        call globalVec2Local(edir,C_dir,res)
        res_edir = sum(res(1:4,1:Nx,1:Ny,:),dim=1) *.25_ireals

        call globalVec2Local(ediff,C_diff,res)
        res_edn = res(2,1:Nx,1:Ny,:)
        res_eup = res(1,1:Nx,1:Ny,:)

        call globalVec2Local(abso,C_one,res)
        res_abso = res(1,1:Nx,1:Ny,:)

        deallocate(res)
     end subroutine

      subroutine globalVec2Local(vec,C,res)
        Vec :: vec
        real(ireals),allocatable :: res(:,:,:,:)
        type(t_coord) :: C
        
        Vec :: natural, local
        VecScatter :: scatter_context

        PetscScalar,Pointer :: xloc(:)

        if(allocated(res)) deallocate(res)
        allocate( res(C%dof,C%glob_xm,C%glob_ym,C%glob_zm) )

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
