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

      use m_data_parameters, only: mpiint, init_mpi_data_parameters
      use m_helper_functions, only: CHKERR
      use m_optprop_LUT, only : t_optprop_LUT, &
        t_optprop_LUT_1_2,  &
        t_optprop_LUT_3_6,  &
        t_optprop_LUT_3_10, &
        t_optprop_LUT_3_16, &
        t_optprop_LUT_8_10, &
        t_optprop_LUT_8_12, &
        t_optprop_LUT_8_16, &
        t_optprop_LUT_8_18, &
        t_optprop_LUT_wedge_5_8, &
        t_optprop_LUT_rectilinear_wedge_5_8, &
        t_optprop_LUT_wedge_18_8

      use m_tenstream_options, only : read_commandline_options

      integer(mpiint) :: myid,comm

      character(len=80) :: arg
      class(t_optprop_LUT), allocatable :: OPP

      PetscErrorCode :: ierr

      call mpi_init(ierr)
      comm = MPI_COMM_WORLD
      call mpi_comm_rank(comm,myid,ierr)
      call PetscInitialize(PETSC_NULL_CHARACTER ,ierr) ;call CHKERR(ierr)

      call init_mpi_data_parameters(MPI_COMM_WORLD)

      call read_commandline_options(comm)

      call get_command_argument(1, arg)

      select case(arg)
      case ('1_2')
        allocate(t_optprop_LUT_1_2::OPP)

      case ('3_6')
        allocate(t_optprop_LUT_3_6::OPP)

      case ('3_10')
        allocate(t_optprop_LUT_3_10::OPP)

      case ('3_16')
        allocate(t_optprop_LUT_3_16::OPP)

      case ('8_10')
        allocate(t_optprop_LUT_8_10::OPP)

      case ('8_12')
        allocate(t_optprop_LUT_8_12::OPP)

      case ('8_16')
        allocate(t_optprop_LUT_8_16::OPP)

      case ('8_18')
        allocate(t_optprop_LUT_8_18::OPP)

      case ('wedge_5_8')
        allocate(t_optprop_LUT_wedge_5_8::OPP)

      case ('rectilinear_wedge_5_8')
        allocate(t_optprop_LUT_rectilinear_wedge_5_8::OPP)

      case ('wedge_18_8')
        allocate(t_optprop_LUT_wedge_18_8::OPP)

      case default
        print *,'error, have to provide solver type as argument, e.g. call with'
        print *,'createLUT 1_2'
        print *,'createLUT 3_6'
        print *,'createLUT 3_10'
        print *,'createLUT 3_16'
        print *,'createLUT 8_10'
        print *,'createLUT 8_12'
        print *,'createLUT 8_16'
        print *,'createLUT 8_18'
        print *,'createLUT wedge_5_8'
        print *,'createLUT rectilinear_wedge_5_8'
        print *,'createLUT wedge_18_8'
        stop
      end select

      call OPP%init(comm, skip_load_LUT=.False.)

      call mpi_finalize(ierr)
end program
