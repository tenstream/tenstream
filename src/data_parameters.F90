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

module m_data_parameters
#include <petsc/finclude/petsc.h>
      use petsc

      use iso_fortran_env, only: INT32, INT64, REAL32, REAL64, REAL128
      use ieee_arithmetic, only: ieee_support_nan, ieee_quiet_nan, ieee_value, &
        & ieee_support_inf, ieee_positive_inf, ieee_negative_inf
      implicit none

      private
      public pi, pi_irealLUT, pi_ireal_params, pi32, pi64, pi128,&
             clight, nil, zero, one,                             &
             i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,inil,         &
             iintegers,mpiint,                                   &
             ireals,ireal128,ireal_dp,irealLUT,                  &
             ireal_params, irealbmc, irealeddington,             &
             nan32, nan64, nan,                                  &
             neginf, inf,                                        &
             imp_iinteger,imp_int4, imp_int8,                    &
             imp_ireals,imp_real_dp,imp_irealLUT,                &
             imp_logical, imp_character,                         &
             imp_REAL32, imp_REAL64,                             &
             init_mpi_data_parameters, finalize_mpi,             &
             default_str_len,                                    &
             EXP_MINVAL, EXP_MAXVAL, EXP_MINVAL128, EXP_MAXVAL128

      integer :: mpiint_dummy
      PetscInt :: petscint_dummy
      PetscReal :: petscreal_dummy

      integer,parameter :: &
          default_str_len = 512,            &
          iintegers = kind(petscint_dummy), &
          irealLUT = REAL32,                &
          ireals = kind(petscreal_dummy),   &
          ireal_params = REAL64, &
          ireal_dp     = REAL64, &
          ireal128     = REAL128, &
!         ireal128 = selected_real_kind(33, 4931), &
          irealbmc     = REAL64, &
          irealeddington = REAL64, &
          mpiint = kind(mpiint_dummy)

      real(ireals),parameter :: pi=3.141592653589793_ireals, clight=299792458._ireals, nil=-9999._ireals
      real(irealLUT),parameter :: pi_irealLUT=3.141592653589793_irealLUT
      real(REAL32),parameter :: pi32=3.141592653589793_REAL32
      real(REAL64),parameter :: pi64=3.141592653589793_REAL64
      real(ireal_params),parameter :: pi_ireal_params=3.141592653589793_ireal_params
      real(REAL128),parameter :: pi128=4 * ATAN (1._REAL128)
      real(ireals),parameter :: zero=0, one=1

      integer(iintegers) ,parameter :: i0=0,i1=1,i2=2,i3=3,i4=4,&
                                       i5=5,i6=6,i7=7,i8=8,i9=9,&
                                       i10=10,i11=11,inil=-9999_iintegers

      real(ireals), parameter :: EXP_MINVAL=epsilon(EXP_MINVAL), EXP_MAXVAL=-log(epsilon(EXP_MAXVAL))
      real(ireal128), parameter :: EXP_MINVAL128=epsilon(EXP_MINVAL), EXP_MAXVAL128=-log(epsilon(EXP_MAXVAL))

      real(real32) :: nan32 =  transfer(-4194304_int32, 1._real32)
      real(real64) :: nan64 =  transfer(-2251799813685248_int64, 1._real64)

      integer(mpiint) :: imp_irealLUT, imp_ireals, imp_real_dp, imp_REAL32, imp_REAL64
      integer(mpiint) :: imp_logical, imp_character
      integer(mpiint) :: imp_iinteger, imp_int4, imp_int8
      real(ireals) :: nan, inf, neginf

contains
subroutine init_mpi_data_parameters(comm)
  integer(mpiint),intent(in) :: comm
  integer(mpiint) :: dtsize, ierr, myid, numnodes, mpierr
  logical :: lmpi_is_initialized, lpetsc_is_initialized, lallsame
  logical :: file_exists

  call mpi_initialized( lmpi_is_initialized, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  if(.not.lmpi_is_initialized) call mpi_init(mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_COMM_RANK( comm, myid, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_Comm_size( comm, numnodes, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)


  call MPI_SIZEOF(i0, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, dtsize, imp_iinteger, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_SIZEOF(1_4, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, dtsize, imp_int4, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_SIZEOF(1_8, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, dtsize, imp_int8, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_SIZEOF(one, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, dtsize, imp_ireals, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_SIZEOF(1._irealLUT, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, dtsize, imp_irealLUT, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_SIZEOF(1._REAL32, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, dtsize, imp_real32, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_SIZEOF(1._REAL64, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, dtsize, imp_real64, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  call MPI_SIZEOF(1._ireal_dp, dtsize, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, dtsize, imp_real_dp, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  imp_logical = mpi_logical
  imp_character = mpi_character

  if(ireal128.lt.i0) then
    if(myid.eq.0) print *,'128 bit reals not supported :( -- '// &
      'you can switch to double precision instead -- '// &
      'beware that the twostream coefficients may not be stable -- '// &
      'please edit data_parameters'
    call mpi_abort(comm, 1_mpiint, ierr)
  endif

  !if(ieee_support_nan(nan32)) nan32=ieee_value(1._real32, ieee_quiet_nan)
  !if(ieee_support_nan(nan64)) nan64=ieee_value(1._real64, ieee_quiet_nan)
  if(ieee_support_inf(inf)) then
    inf   = ieee_value(1._ireals, ieee_positive_inf)
    neginf= ieee_value(1._ireals, ieee_negative_inf)
  else
    inf   = huge(1._ireals)
    neginf= -huge(1._ireals)
  endif

!  if(myid.eq.0) print *,myid,'init_mpi_data_parameters :: imp_int',imp_int,' :: imp_real',imp_real,'epsilon(real)',epsilon(one)
!  print *,'init_mpi_data_parameters :: MPI_INTEGER',MPI_INTEGER,' :: MPI_DOUBLE_PRECISION',MPI_DOUBLE_PRECISION,' :: MPI_REAL',MPI_REAL
  select case(kind(nan))
  case (kind(nan32))
    nan = real(nan32, kind(nan))
  case (kind(nan64))
    nan = real(nan64, kind(nan))
  case default
    nan = real(-99999, kind(nan))
  end select

  call PetscInitialized(lpetsc_is_initialized, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  lallsame = mpi_logical_all_same(comm,lpetsc_is_initialized)
  if(.not.lallsame) then
    print *,myid,'the provided communicator does not agree on lpetsc_is_initialized',lpetsc_is_initialized
    call mpi_abort(comm, 1_mpiint, ierr)
  endif

  PETSC_COMM_WORLD = comm
  if(.not.lpetsc_is_initialized) then
    inquire(file='tenstream.options', exist=file_exists)
    if(file_exists) then
      call PetscInitialize('tenstream.options', mpierr)
    else
      call PetscInitialize(PETSC_NULL_CHARACTER, mpierr)
    endif
    if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  endif
end subroutine

subroutine finalize_mpi(comm, lfinalize_mpi, lfinalize_petsc)
  integer(mpiint),intent(in) :: comm
  logical, intent(in) :: lfinalize_mpi, lfinalize_petsc
  integer(mpiint) :: ierr, mpierr
  logical :: lmpi_is_initialized, lpetsc_is_initialized

  call mpi_initialized( lmpi_is_initialized, mpierr)
  if(.not.lmpi_is_initialized) return ! if we dont even have mpi, petsc cant live either

  call PetscInitialized(lpetsc_is_initialized, mpierr)
  if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)

  if(lpetsc_is_initialized) then
    if(lfinalize_petsc) then
      call PetscFinalize(mpierr)
      if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
    endif
  endif

  if(lfinalize_mpi) then
    if(lmpi_is_initialized) call MPI_Finalize(mpierr)
    if(mpierr.ne.0) call mpi_abort(comm, mpierr, ierr)
  endif
end subroutine

!duplicate of helper function. otherwise get circular dependency
function mpi_logical_all_same(comm,lval) result(lsame)
  integer(mpiint),intent(in) :: comm
  logical :: lsame
  logical,intent(in) :: lval
  integer(mpiint) :: i, isum, commsize, ierr
  call mpi_comm_size(comm, commsize, ierr)
  if(commsize.gt.1_mpiint) then
    if(lval) then
      i = 1
    else
      i = 0
    endif
    call mpi_allreduce(i, isum, 1_mpiint, imp_int4, MPI_SUM, comm, ierr)
    if(lval) then
      lsame = isum.eq.commsize
    else
      lsame = isum.eq.0
    endif
  else
    lsame = .True.
  endif
end function
end module
