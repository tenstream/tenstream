module m_petsc_ts

!      use m_gridtransform
      use m_netcdfio
      use m_data_parameters, only : init_mpi_data_parameters,mpiint,ireals,iintegers,&
          imp_comm,myid, &
          nil,zero,one
      use m_kato_data
      use m_helper_functions, only : imp_bcast

      use m_tenstream, only : init_tenstream, set_global_optical_properties, solve_tenstream, destroy_tenstream, &
                            b,edir,ediff,abso, &
                            edir_twostr,ediff_twostr,abso_twostr

      use m_tenstream_options, only : read_commandline_options,lwriteall, &
                                      ident_dx,ident_dy,options_phi,options_theta,ident,basepath,output_prefix,ltwostr

      implicit none
#include "finclude/petsc.h90"

      PetscErrorCode :: ierr
      integer(iintegers) :: iierr

      contains

subroutine vec_from_hdf5(v,err_code)
      Vec :: v
      character(10),parameter :: suffix='.h5'
      character(110) :: fname
      logical fexists
      PetscFileMode :: fmode
      character(100) :: vecname,s_theta0
      PetscScalar :: rvecsum

      PetscErrorCode :: err_code
      
      PetscViewer view
      PetscErrorCode ierr

      call PetscObjectGetName(v,vecname,ierr) ;CHKERRQ(ierr)

      write(s_theta0,FMT='(".",I0)' ) int(options_theta)
      
      fname = trim(basepath) // trim(output_prefix) // '.' // trim(ident) // trim(s_theta0) // trim(suffix)
      inquire(file=trim(fname), exist=fexists)
      
      if(fexists) then
        if(myid.eq.0)  print *,myid,'reading vector from hdf5 file ',trim(fname),' vecname: ',vecname
        fmode = FILE_MODE_READ

!        call PetscViewerHDF5Open(imp_comm,trim(fname),fmode, view, ierr) ; CHKERRQ(ierr)
        call VecLoad(v, view, err_code)                                          ; CHKERRQ(ierr)
        call PetscViewerDestroy(view,ierr)                                       ; CHKERRQ(ierr)

      else
        err_code = 1378 ! Random Number(thought by Fabian) to show that file does not exist

      endif

      if(err_code.ne.0) then !Somehow couldnt read vector - set it to zero.
        call VecSet(v,zero,ierr) ; CHKERRQ(ierr)
      endif

      call VecSum(v,rvecsum,ierr) ;CHKERRQ(ierr)
      if(myid.eq.0) print *,myid,'reading from hdf5 file done :: error:',err_code,'sum of vector',rvecsum
end subroutine 
subroutine vec_to_hdf5(v)
      Vec,intent(in) :: v
      character(10),parameter :: suffix='.h5'
      character(110) :: fname
      logical fexists
      PetscFileMode :: fmode
      character(100) :: vecname,s_theta0
      
      PetscViewer :: view

      call PetscObjectGetName(v,vecname,ierr) ;CHKERRQ(ierr)

      write(s_theta0,FMT='(".",I0)' ) int(options_theta)
      
      fname = trim(basepath) // trim(output_prefix) // '.' // trim(ident) // trim(s_theta0) // trim(suffix)
      inquire(file=trim(fname), exist=fexists)
      
      if(fexists) then
        if(myid.eq.0)  print *,myid,'appending vector to hdf5 file ',trim(fname),' vecname: ',vecname
        fmode = FILE_MODE_APPEND
      else 
        if(myid.eq.0)  print *,myid,'writing vector to hdf5 file ',trim(fname),' vecname: ',vecname
        fmode = FILE_MODE_WRITE
      endif

!      call PetscViewerHDF5Open(imp_comm,trim(fname),fmode, view, ierr) ;CHKERRQ(ierr)
      call VecView(v, view, ierr) ;CHKERRQ(ierr)
      call PetscViewerDestroy(view,ierr) ;CHKERRQ(ierr)

      if(myid.eq.0 ) print *,myid,'writing to hdf5 file done'
end subroutine

  subroutine load_optprop(kato,iq,kabs,ksca,g,hhl,planck)
      integer(iintegers),intent(in) :: kato,iq
      real(ireals),allocatable,dimension(:) :: hhl
      real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g,planck

!      character(300) :: skato,siq
      character(300),parameter :: uvspec_file='/usr/users/jakub/cosmodata/tenstream/scenes/scenes.h5'
      character(300) :: groups(10)

      ! Make sure everything is cleaned
      if(allocated(kabs) ) deallocate(kabs)
      if(allocated(ksca) ) deallocate(ksca)
      if(allocated(g   ) ) deallocate(g   )
      if(allocated(hhl ) ) deallocate(hhl )

      if(allocated(planck))deallocate(planck)


      ! run test if we dont know which ident to run
      if(ident.eq.'run_test') then
        call load_test_optprop(kato,iq,kabs,ksca,g,hhl)
        return
      endif

      groups(1) = uvspec_file
      groups(2) = ident

      !Actually loading data

      groups(3) = 'hhl'
      if(myid.eq.0) then
        call ncload(groups(1:3),hhl,iierr)

        if(ierr.ne.0) then
          print *,'Tried loading hhl from ',trim(groups(1)),' ::  ',trim(groups(2)),'  ::  ',trim(groups(3))
          stop 'Error occured loading hhl'
        endif

        write(groups(3),FMT='("kato",I0)') kato
        write(groups(4),FMT='("iq",I0)') iq
        groups(5) = 'kabs'
        call ncload(groups(1:5),kabs,iierr)
        groups(5) = 'ksca'
        call ncload(groups(1:5),ksca,iierr) 
        groups(5) = 'g'
        call ncload(groups(1:5),g,iierr) 

        if(size(ksca).ne.size(kabs).or.size(kabs).ne.size(g).or.size(hhl).ne.ubound(kabs,3)+1) then !.or.(ubound(hhl1d,1)-1.ne.ubound(kabs,3))) then
          print *,'ERROR : shapes of optical properties do not match!!!'
          print *,'shape(kabs)',shape(kabs)
          print *,'shape(ksca)',shape(ksca)
          print *,'shape(g)',shape(g)
          print *,'shape(hhl)',shape(hhl)

          call exit()
        endif

      endif

      if(options_theta.lt.zero) then 
        !TODO this is just for debug purposes... we
        !     do not maintain the calculation of planck emission for ident test
        !     cases.... if you want to calc thermal heating rates, please use the
        !     tenstream solver in libRadtran.
        allocate(planck(ubound(kabs,1),ubound(kabs,2),ubound(hhl,1) ))
        planck = 100
      endif

      call imp_bcast(hhl ,0_mpiint,myid)
      call imp_bcast(kabs,0_mpiint,myid)
      call imp_bcast(ksca,0_mpiint,myid)
      call imp_bcast(g   ,0_mpiint,myid)

    contains
      subroutine load_test_optprop(kato,iq,kabs,ksca,g,hhl)
          integer(iintegers),intent(in) :: kato,iq
          real(ireals),allocatable,dimension(:) :: hhl

          real(ireals),allocatable,dimension(:,:,:) :: kabs,ksca,g

          integer(iintegers),parameter :: glob_Nx=40,glob_Ny=40,glob_Nz=40, zTOA=glob_Nz*40
          integer(iintegers) :: k

          if(myid.eq.0) print *,myid,'Creating Optical Properties here instead of taking them from kato',kato,iq

          allocate(kabs(glob_Nx,glob_Ny,glob_Nz))
          allocate(ksca(glob_Nx,glob_Ny,glob_Nz))
          allocate(g   (glob_Nx,glob_Ny,glob_Nz))
          allocate(hhl (glob_Nz+1))

          hhl(ubound(hhl,1))=0
          do k=ubound(hhl,1),2,-1
            hhl(k-1)=hhl(k)+zTOA/glob_Nz
          enddo

          g=0.5
          kabs=1e-4
          ksca=1e-4

          if(myid.eq.0) ksca(1:glob_Nx:5,1:glob_Ny:5,5) = 5e-3
          !      ksca(1,1,3:5) = 1e-2

  end subroutine
  end subroutine

  subroutine init_integral_vecs(intedir,intediff,intabso,intedir_twostr, intediff_twostr, intabso_twostr, linit)
        Vec :: intedir,intediff,intabso
        Vec :: intedir_twostr, intediff_twostr, intabso_twostr
        logical :: linit

        character(100) :: vecname

        if(ltwostr) then
          call VecDuplicate(edir_twostr , intedir_twostr , ierr) ; CHKERRQ(ierr)
          write(vecname,FMT='("edir.twostr.",I0,".",I0)') int(options_phi),int(options_theta) 
          call PetscObjectSetName(intedir_twostr,vecname,ierr)          ; CHKERRQ(ierr)

          call VecDuplicate(ediff_twostr, intediff_twostr, ierr) ; CHKERRQ(ierr)
          write(vecname,FMT='("ediff.twostr.",I0,".",I0)') int(options_phi),int(options_theta) 
          call PetscObjectSetName(intediff_twostr,vecname,ierr)          ; CHKERRQ(ierr)

          call VecDuplicate(abso_twostr , intabso_twostr , ierr) ; CHKERRQ(ierr)
          write(vecname,FMT='("abso.twostr.",I0,".",I0)') int(options_phi),int(options_theta) 
          call PetscObjectSetName(intabso_twostr,vecname,ierr)          ; CHKERRQ(ierr)

          call VecSet(intedir_twostr ,zero,ierr) ;CHKERRQ(ierr)
          call VecSet(intediff_twostr,zero,ierr) ;CHKERRQ(ierr)
          call VecSet(intabso_twostr ,zero,ierr) ;CHKERRQ(ierr)
        endif

        call VecDuplicate(edir , intedir , ierr) ; CHKERRQ(ierr)
        write(vecname,FMT='("edir.",I0,".",I0)') int(options_phi),int(options_theta) 
        call PetscObjectSetName(intedir,vecname,ierr)          ; CHKERRQ(ierr)

        call VecDuplicate(ediff, intediff, ierr) ; CHKERRQ(ierr)
        write(vecname,FMT='("ediff.",I0,".",I0)') int(options_phi),int(options_theta) 
        call PetscObjectSetName(intediff,vecname,ierr)          ; CHKERRQ(ierr)

        call VecDuplicate(abso , intabso , ierr) ; CHKERRQ(ierr)
        write(vecname,FMT='("abso.",I0,".",I0)') int(options_phi),int(options_theta) 
        call PetscObjectSetName(intabso,vecname,ierr)          ; CHKERRQ(ierr)

        call VecSet(intedir ,zero,ierr) ;CHKERRQ(ierr)
        call VecSet(intediff,zero,ierr) ;CHKERRQ(ierr)
        call VecSet(intabso ,zero,ierr) ;CHKERRQ(ierr)

        linit = .True.

    end subroutine

  subroutine dump_vectors(kato,iq)
        integer(iintegers),intent(in) :: kato,iq
        character(100) :: vecname

        if(ltwostr) then
          write(vecname,FMT='("edir.twostr.",I0,".",I0,"-",I0,"-",I0)') int(options_phi),int(options_theta),kato,iq
          call PetscObjectSetName(edir_twostr,vecname,ierr) ; CHKERRQ(ierr)
          call vec_to_hdf5(edir_twostr)

          write(vecname,FMT='("ediff.twostr.",I0,".",I0,"-",I0,"-",I0)') int(options_phi),int(options_theta),kato,iq
          call PetscObjectSetName(ediff_twostr,vecname,ierr) ; CHKERRQ(ierr)
          call vec_to_hdf5(ediff_twostr)

          write(vecname,FMT='("abso.twostr.",I0,".",I0,"-",I0,"-",I0)') int(options_phi),int(options_theta),kato,iq
          call PetscObjectSetName(abso_twostr,vecname,ierr) ; CHKERRQ(ierr)
          call vec_to_hdf5(abso_twostr)
        endif

        write(vecname,FMT='("b.",I0,".",I0,"-",I0,"-",I0)') int(options_phi),int(options_theta),kato,iq
        call PetscObjectSetName(b,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(b)

        write(vecname,FMT='("edir.",I0,".",I0,"-",I0,"-",I0)') int(options_phi),int(options_theta),kato,iq
        call PetscObjectSetName(edir,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(edir)

        write(vecname,FMT='("ediff.",I0,".",I0,"-",I0,"-",I0)') int(options_phi),int(options_theta),kato,iq
        call PetscObjectSetName(ediff,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(ediff)

        write(vecname,FMT='("abso.",I0,".",I0,"-",I0,"-",I0)') int(options_phi),int(options_theta),kato,iq
        call PetscObjectSetName(abso,vecname,ierr) ; CHKERRQ(ierr)
        call vec_to_hdf5(abso)

    end subroutine

end module

program main
        use m_petsc_ts
        implicit none

        Vec :: intedir,intediff,intabso
        Vec :: intedir_twostr, intediff_twostr, intabso_twostr

        integer(iintegers) :: iq, kato, k
        integer(iintegers) :: dims(3)

        real(ireals),allocatable,dimension(:,:,:) :: global_kabs,global_ksca,global_g, global_planck
        real(ireals),allocatable :: hhl(:),dz(:)
        real(ireals),parameter :: albedo = 0.05

        logical :: linit_integral_vecs = .False.

        call PetscInitialize(PETSC_NULL_CHARACTER,ierr) ;CHKERRQ(ierr)
        call init_mpi_data_parameters(PETSC_COMM_WORLD)
        call read_commandline_options()

        call load_optprop(1_iintegers,0_iintegers, global_kabs,global_ksca,global_g, hhl, global_planck)

        dims = shape(global_kabs)
        allocate(dz(dims(3)))
        do k=1,dims(3)
          dz(k) = hhl(k)-hhl(k+1)
        enddo
        call init_tenstream(imp_comm, dims(1),dims(2),dims(3), ident_dx, ident_dy, options_phi,options_theta,albedo,dz1d=dz )

        do kato=1,32
!                  do kato=11,11
          do iq=0,kato_bands(kato)
            if(myid.eq.0) print *,'-----------------------------------------------------------------------------------------------------------------------------'
            if(myid.eq.0) print *,'-------------------------- Calculate ',trim(ident),' sza',options_theta,' kato',kato,'iq',iq
            if(myid.eq.0) print *,'-----------------------------------------------------------------------------------------------------------------------------'

            call load_optprop(kato,iq, global_kabs,global_ksca,global_g, hhl, global_planck)
            if(allocated(global_planck)) then
              call set_global_optical_properties(global_kabs, global_ksca, global_g, global_planck)
            else
              call set_global_optical_properties(global_kabs, global_ksca, global_g)
            endif

            call solve_tenstream( get_edirTOA(kato,iq,hhl(1)) )

            if(.not.linit_integral_vecs) call init_integral_vecs(intedir,intediff,intabso,intedir_twostr, intediff_twostr, intabso_twostr, linit_integral_vecs)

            if(ltwostr) then
              call VecAXPY(intedir_twostr ,one,edir_twostr ,ierr) ;CHKERRQ(ierr)
              call VecAXPY(intediff_twostr,one,ediff_twostr,ierr) ;CHKERRQ(ierr)
              call VecAXPY(intabso_twostr ,one,abso_twostr ,ierr) ;CHKERRQ(ierr)
            endif

            call VecAXPY(intedir ,one,edir ,ierr) ;CHKERRQ(ierr)
            call VecAXPY(intediff,one,ediff,ierr) ;CHKERRQ(ierr)
            call VecAXPY(intabso ,one,abso ,ierr) ;CHKERRQ(ierr)

            if(lwriteall) call dump_vectors(kato,iq)

          enddo
        enddo

        if(ltwostr) then
          call vec_to_hdf5(intabso_twostr)
          call vec_to_hdf5(intedir_twostr)
          call vec_to_hdf5(intediff_twostr)

          print *,'Cleanup Result vectors'
          call VecDestroy(intedir_twostr,ierr) ;CHKERRQ(ierr)
          call VecDestroy(intediff_twostr   ,ierr) ;CHKERRQ(ierr)
          call VecDestroy(intabso_twostr,ierr) ;CHKERRQ(ierr)
          call mpi_barrier(imp_comm,ierr)
        endif

        call vec_to_hdf5(intabso)
        call vec_to_hdf5(intedir)
        call vec_to_hdf5(intediff)

        call VecDestroy(intedir,ierr) ;CHKERRQ(ierr)
        call VecDestroy(intediff,ierr) ;CHKERRQ(ierr)
        call VecDestroy(intabso,ierr) ;CHKERRQ(ierr)

        call destroy_tenstream()
        call PetscFinalize(ierr) ;CHKERRQ(ierr)
end program
