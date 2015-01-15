module m_data_parameters

#ifdef _XLF
        use mpi 
#else
        use mpi ,only:mpi_sizeof, mpi_type_match_size 
#endif

        implicit none
#include "finclude/petsc.h90"


      private
      public pi,clight,nil,zero,one,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,inil, &
             iintegers,ireals,ireal128,ireal_dp,                                &
             mpiint,imp_int,imp_real,imp_logical,imp_comm,myid,numnodes,mpierr, &
             init_mpi_data_parameters

!      integer,parameter :: &
!      ireals=selected_real_kind(8,100), &
!      ireals=selected_real_kind(16,200), &
!      iintegers=selected_int_kind(12)

      integer :: mpiint_dummy
      PetscInt :: petscint_dummy
      PetscReal :: petscreal_dummy

      integer,parameter :: &
          iintegers = kind(petscint_dummy), &
          ireals = kind(petscreal_dummy),   &
!          ireal128 = selected_real_kind(33, 4931), &
          ireal128 = selected_real_kind(15, 307), &
          ireal_dp = selected_real_kind(15, 307), &
          mpiint = kind(mpiint_dummy)

      real(ireals),parameter :: pi=3.141592653589793_ireals, clight=299792458._ireals, nil=-9999._ireals
      real(ireals),parameter :: zero=0._ireals, one=1._ireals
      integer(iintegers) ,parameter :: i0=0,i1=1,i2=2,i3=3,i4=4,i5=5,i6=6,i7=7,i8=8,i9=9,i10=10,i11=11,inil=-9999_iintegers


      integer(mpiint) :: imp_int, imp_real, imp_logical, imp_comm
      integer(mpiint) :: myid,numnodes,mpierr

contains 
subroutine init_mpi_data_parameters(comm)
  integer,intent(in) :: comm
  integer :: size,ierror

  imp_comm = comm

  call MPI_COMM_RANK( imp_comm, myid, mpierr )
  call MPI_Comm_size( imp_comm, numnodes, mpierr)

  call MPI_SIZEOF(i0, size, ierror)    
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER, size, imp_int, ierror)

  call MPI_SIZEOF(one, size, ierror)    
  call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_REAL, size, imp_real, ierror)

  imp_logical = mpi_logical

  if(ireal128.lt.i0) then
    if(myid.eq.0) print *,'128 bit reals not supported :( -- you can switch to double precision instead -- beware that the twostream coefficients may not be stable -- please edit data_parameters'
  endif

  if(myid.eq.0) print *,myid,'init_mpi_data_parameters :: imp_int',imp_int,' :: imp_real',imp_real,'epsilon(real)',epsilon(one)
!  print *,'init_mpi_data_parameters :: MPI_INTEGER',MPI_INTEGER,' :: MPI_DOUBLE_PRECISION',MPI_DOUBLE_PRECISION,' :: MPI_REAL',MPI_REAL
end subroutine
end module
