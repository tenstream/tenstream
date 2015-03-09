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
#include "finclude/petscdef.h"
      use petsc
      use m_data_parameters, only: mpiint, ireals, init_mpi_data_parameters
      use mpi
      use m_optprop_LUT, only : t_optprop_LUT_8_10
      use m_tenstream_options, only : read_commandline_options

      integer(mpiint) :: myid,comm

      character(len=80) :: arg
      real(ireals) :: dx,user_sza
      real(ireals) :: azis(2),szas(6)

      type(t_optprop_LUT_8_10) :: OPP

      PetscErrorCode :: ierr
      integer :: i

      call mpi_init(ierr)
      comm = MPI_COMM_WORLD
      call mpi_comm_rank(comm,myid,ierr)
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr) ;CHKERRQ(ierr)

      call init_mpi_data_parameters(MPI_COMM_WORLD)

      call read_commandline_options()

      azis = [0,90]
      szas = [-1,0,20,40,60,80]

      dx=-1
      do i=1,10
        call get_command_argument(i, arg)
        if(len_trim(arg) .gt. 0) then

          if(arg.eq.'-dx') then
            call get_command_argument(i+1, arg)
            read (arg,*) dx
          endif

          if(arg.eq.'-sza') then 
            call get_command_argument(i+1, arg)
            read (arg,*) user_sza
            szas=user_sza
          endif
        endif
      enddo

      if(dx.eq.-1) stop 'Need to supply option dx to create Lookuptables... stopping' 

      print *,'calculating coeffs for dx',dx,'szas',szas
      call OPP%init(dx,dx,azis,szas,comm)
      print *,'loaded 8_10 coeffs for dx',dx,'szas',szas,'azis',azis

      call mpi_finalize(ierr)
end program
