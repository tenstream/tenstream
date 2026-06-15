program main
#ifdef HAVE_PETSC
#include "petsc/finclude/petsc.h"
  use petsc
#endif
  use mpi
  use m_tenstream_options, only: read_commandline_options
  use m_data_parameters, only: mpiint, default_str_len
  use m_example_pprts_specint_lw_sw_from_dump, only: ex_pprts_specint_lw_sw_from_dump
  use m_helper_functions, only: CHKERR, get_petsc_opt

  implicit none

  integer(mpiint) :: ierr, comm, myid
  logical :: lflg
  character(len=default_str_len) :: inpfile, outfile, specint

  comm = MPI_COMM_WORLD
  call mpi_init(ierr)
  call mpi_comm_rank(comm, myid, ierr)

#ifdef HAVE_PETSC
  call PetscInitialize('', ierr)
#endif

  call read_commandline_options(comm)

  specint = 'no_default_set'
  call get_petsc_opt('', "-specint", specint, lflg, ierr); call CHKERR(ierr)

  call get_petsc_opt('', '-inp', inpfile, lflg, ierr); call CHKERR(ierr)
  if (.not. lflg) call CHKERR(1_mpiint, 'need to supply a input filename... please call with -inp <input.nc>')

  outfile = ''
  call get_petsc_opt('', '-out', outfile, lflg, ierr); call CHKERR(ierr)

  call ex_pprts_specint_lw_sw_from_dump(specint, comm, inpfile, outfile)

  call mpi_finalize(ierr)
end program
