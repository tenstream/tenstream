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

program main

#include "petsc/finclude/petsc.h"
      use petsc

      use m_data_parameters, only: mpiint, ireals, init_mpi_data_parameters
      use m_helper_functions, only: CHKERR
      use mpi
      use m_optprop_LUT, only : t_optprop_LUT_3_6
      use m_tenstream_options, only : read_commandline_options

      integer(mpiint) :: myid,comm

      type(t_optprop_LUT_3_6) :: OPP

      PetscErrorCode :: ierr

      call mpi_init(ierr)
      comm = MPI_COMM_WORLD
      call mpi_comm_rank(comm,myid,ierr)
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr) ;call CHKERR(ierr)

      call init_mpi_data_parameters(MPI_COMM_WORLD)

      call read_commandline_options(comm)

      call OPP%init(comm)

      call mpi_finalize(ierr)
end program
