module m_tenstream_options
      use m_data_parameters, only : ireals,iintegers,myid,numnodes,one,i0,imp_comm
      use m_optprop_parameters, only: lut_basename
      implicit none
#include "finclude/petsc.h90"

      logical :: ltwostr, & ! additionally calculate delta eddington twostream solution
        lwriteall, & ! write out each and every solution -- mainly good to debug solver
        luse_eddington, & ! use delta eddington coefficients for upper atmosphere, if False, we use boxmc 2-str coeffs
        luse_hdf5_guess, & ! try loading initial guess from file
        luse_twostr_guess ! use twostream solution as first guess
      real(ireals) :: twostr_ratio,&
            ident_dx,&
            ident_dy,&
            options_phi,&
            options_theta,&
            options_max_solution_err, options_max_solution_time
                    
      integer(iintegers) :: pert_xshift,pert_yshift

      character(len=300) :: ident,output_prefix
      character(len=300) :: basepath

      contains 
        subroutine show_options()
          print *,'------------------------------------------------------------------------------------------------------------------'  
          print *,'------------------------------------------------------------------------------------------------------------------'  
          print *,'Tenstream options:'
          print *,'-show_options         :: show this text                                                                           '  
          print *,'-ident <run_*>        :: load optical properties from hdf5 -read petsc_solver::load_optprop (default = run_test)  '
          print *,'-ident run_test       :: load optical properties from function in petsc_solver::load_test_optprop                 '
          print *,'-out                  :: output prefix (default = ts)                                                             '  
          print *,'-basepath             :: output directory (default = ./)                                                          '  
          print *,'-dx -dy               :: domain size in [m] (mandatory if running with -ident <run_*> )                           '  
          print *,'-phi -theta           :: solar azimuth and zenith angle (default = (180,0) == south,overhead sun)                 '  
          print *,'-writeall             :: dump intermediate results                                                                '  
          print *,'-twostr               :: calculate delta eddington twostream solution                                             ' 
          print *,'-hdf5_guess           :: if run earlier with -writeall can now use dumped solutions as initial guess              '  
          print *,'-twostr_guess         :: use delta eddington twostream solution as first guess                                    '  
          print *,'-twostr_ratio <limit> :: when aspect ratio (dz/dx) is smaller than <limit> then we use twostr_coeffs(default = 1.)'  
          print *,'-pert_xshift <i>      :: shift optical properties in x direction by <i> pixels                                    '  
          print *,'-pert_yshift <j>      :: shift optical properties in Y direction by <j> pixels                                    '  
          print *,'-max_solution_err [W] :: if max error of solution is estimated below this value, skip calculation                 '  
          print *,'-max_solution_time[s] :: if last update of solution is older, update irrespective of estimated error              '  
          print *,'-lut_basename         :: path to LUT table files -- default is local dir                                          '  
          print *,'------------------------------------------------------------------------------------------------------------------'  
          print *,'------------------------------------------------------------------------------------------------------------------'  
        end subroutine
        subroutine read_commandline_options()
          logical :: lflg=.False.,lflg_ident=.False.
          PetscErrorCode :: ierr
          logical :: lshow_options=.False.
          logical,save :: linit=.False.
          if(linit) return

          call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-show_options",lshow_options,lflg,ierr) ;CHKERRQ(ierr)
          if(lshow_options) then
            if(myid.eq.0) call show_options()
            call mpi_abort(imp_comm,ierr)
          endif

          call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-ident',ident,lflg_ident,ierr) ; CHKERRQ(ierr)
          if(lflg_ident.eqv.PETSC_FALSE) ident = 'run_test'

          call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-out',output_prefix,lflg,ierr) ; CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) output_prefix = 'ts'

          call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-basepath',basepath,lflg,ierr) ; CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) basepath = './'

          call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-dx",ident_dx, lflg,ierr)  ; CHKERRQ(ierr)
          if( (lflg.eqv.PETSC_FALSE) .and. (lflg_ident.eqv.PETSC_TRUE) ) then
            print *,'If we run with -ident, you need to specify "-dx" commandline option e.g. -dx 70'
            call mpi_abort(imp_comm,ierr)
          endif

          call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-max_solution_err",options_max_solution_err, lflg,ierr)  ; CHKERRQ(ierr)
          if (lflg.eqv.PETSC_FALSE ) options_max_solution_err=0.01
          call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-max_solution_time",options_max_solution_time, lflg,ierr)  ; CHKERRQ(ierr)
          if (lflg.eqv.PETSC_FALSE ) options_max_solution_time=60

          call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-dy",ident_dy, lflg,ierr)  ; CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) ident_dy = ident_dx

          options_phi=180; options_theta=-1
          call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-phi"  ,options_phi, lflg,ierr)     ; CHKERRQ(ierr)
          call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-theta",options_theta, lflg,ierr) ; CHKERRQ(ierr)

          call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-eddington",luse_eddington,lflg,ierr) ;CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) luse_eddington = .True.

          call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-writeall",lwriteall,lflg,ierr) ;CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) lwriteall = .False.

          call PetscOptionsGetBool(PETSC_NULL_CHARACTER , "-twostr" , ltwostr , lflg , ierr) ;CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) ltwostr = .False.

          call PetscOptionsGetBool(PETSC_NULL_CHARACTER , "-hdf5_guess"   , luse_hdf5_guess   , lflg , ierr) ;CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) luse_hdf5_guess = .False.

          call PetscOptionsGetBool(PETSC_NULL_CHARACTER , "-twostr_guess" , luse_twostr_guess , lflg , ierr) ;CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) luse_twostr_guess = .False.
          if(luse_twostr_guess) ltwostr = .True.

          if(luse_twostr_guess.and.luse_hdf5_guess) then
            print *,'cant use twostr_guess .AND. hdf5_guess at the same time'
            call mpi_abort(imp_comm,ierr)
          endif

          call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-twostr_ratio",twostr_ratio, lflg,ierr)  ; CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) twostr_ratio=.5_ireals

          call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-pert_xshift",pert_xshift, lflg,ierr) ; CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) pert_xshift=0
          call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-pert_yshift",pert_yshift, lflg,ierr) ; CHKERRQ(ierr)
          if(lflg.eqv.PETSC_FALSE) pert_yshift=0

          call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-lut_basename',lut_basename,lflg,ierr) ; CHKERRQ(ierr)

          if(myid.eq.0) then
            print *,'********************************************************************'
            print *,'***   Running Job: ',ident
            print *,'***   nr. of Nodes:',numnodes
            print *,'***   dx,dy        ',ident_dx,ident_dy
            print *,'***   eddington    ',luse_eddington
            print *,'***   writeall     ',lwriteall
            print *,'***   twostr       ',ltwostr
            print *,'***   twostr_guess ',luse_twostr_guess
            print *,'***   hdf5_guess   ',luse_hdf5_guess
            print *,'***   twostr_ratio ',twostr_ratio
            print *,'***   out          ',trim(output_prefix)
            print *,'***   solar azimuth',options_phi
            print *,'***   solar zenith ',options_theta
            print *,'***   size_of ireal/iintegers',sizeof(one),sizeof(i0)
            print *,'***   max_solution_err       ',options_max_solution_err
            print *,'***   max_solution_time      ',options_max_solution_time
            print *,'***   lut_basename           ',trim(lut_basename)
            print *,'********************************************************************'
            print *,''
          endif

          call mpi_barrier(imp_comm,ierr)

          linit=.True.
      end subroutine
end module
