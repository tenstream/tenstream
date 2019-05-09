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

module m_twostream

#ifdef _XLF
  use ieee_arithmetic
#define isnan ieee_is_nan
#endif

#include "petsc/finclude/petsc.h"
  use petsc

  use iso_fortran_env, only: REAL32, REAL64

  use m_data_parameters, only: ireals,iintegers,mpiint,zero,one,pi,i0,i1,i2
  use m_eddington, only: eddington_coeff_zdun
  use m_helper_functions, only : delta_scale_optprop, CHKERR, itoa
  implicit none

  private
  public delta_eddington_twostream, petsc_delta_eddington_twostream, adding_delta_eddington_twostream

  logical, parameter :: ldebug=.False.

contains

  subroutine delta_eddington_twostream(dtau,w0,g,mu0,incSolar,albedo, S,Edn,Eup, planck)
    real(ireals),intent(in),dimension(:) :: dtau,w0,g
    real(ireals),intent(in) :: albedo,mu0,incSolar
    real(ireals),dimension(:),intent(out):: S,Edn,Eup
    real(ireals),dimension(:),intent(in),optional :: planck

    real(ireals),dimension(size(dtau)) :: a11,a12,a13,a23,a33,g1,g2

    integer(iintegers) :: i,j,k,ke,ke1,bi
    real(ireals) :: R,T, emis, b0, b1
    real(ireals),allocatable :: AB (:,:)
    real(ireals),allocatable :: B (:,:)
    integer(iintegers), allocatable :: IPIV(:)
    integer(iintegers) :: N, KLU,  KL, KU, NRHS, LDAB, LDB, INFO

    ke = size(dtau)
    ke1 = ke+1
    N = int(2*ke1, kind(N))
    KL =2
    KU =2
    NRHS=1
    LDAB = 2*KL+KU+1
    LDB  = N
    KLU = KL+KU+1

    do k=1,ke
      call eddington_coeff_zdun (dtau(k), w0(k),g(k), mu0,a11(k),a12(k),a13(k),a23(k),a33(k), g1(k),g2(k) )
      if(ldebug) then
        if(any(isnan( [a11(k),a12(k),a13(k),a23(k),a33(k),g1(k),g2(k)] )) ) then
          print *,'eddington',k,' :: ',dtau(k), w0(k),g(k), mu0,'::',a11(k),a12(k),a13(k),a23(k),a33(k),'::',g1(k),g2(k)
          call exit()
        endif
      endif
    enddo

    if(mu0.gt.zero) then
      S(1) = incSolar ! irradiance on tilted plane
      do k=1,ke
        S(k+1) = S(k) * a33(k)
      enddo
    else
      S = zero
    endif

    allocate( IPIV(N))
    allocate( AB (LDAB,N), source=zero)
    allocate( B (LDB,NRHS))

    ! Setup solar src vector
    do k=1,ke
      B(2*k-1,1) = S(k) * a13(k) ! Eup
      B(2*k+2,1) = S(k) * a23(k) ! Edn
    enddo
    B(2,1) = zero ! no Edn at TOA
    B(2*ke1-1,1) = S(ke1) * albedo

    ! Setup thermal src vector
    if(present(planck) ) then
      do k=1,ke
        emis = max(zero, min(one, one-a11(k)-a12(k))) * pi
        b0 = emis*planck(k  ) + (one-emis)*planck(k+1)
        b1 = emis*planck(k+1) + (one-emis)*planck(k  )
        B(2*k-1,1) = B(2*k-1,1) + emis * b0
        B(2*k+2,1) = B(2*k+2,1) + emis * b1
      enddo
      B(2*ke1-1,1) = B(2*ke1-1,1) + planck(ke1)*(one-albedo)*pi
    endif

    !diagonal entries
    do i=1,N
      j=i
      bi= KLU+i-j
      AB( bi,j ) = one
    enddo

    do k=1,ke
      T = a11(k)
      R = a12(k)

      ! setting Eup coeffs
      i=2*k-1 ; j=i+2 ; bi= KLU+i-j ; AB(bi,j) = -T
      i=2*k-1 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = -R

      ! setting Edn coeffs
      i=2*(k+1) ; j=i-2 ; bi= KLU+i-j ; AB(bi,j) = -T
      i=2*(k+1) ; j=i-1 ; bi= KLU+i-j ; AB(bi,j) = -R
    enddo
    i=2*ke1-1 ; j=i+1 ; bi= KLU+i-j ; AB(bi,j) = -albedo ! Eup at surface is Edn*albedo

    if(ldebug) then
      do k=1,ke
        if(any(isnan( [a11(k),a12(k),a13(k),a23(k),a33(k)] )) ) then
          print *,'eddington coefficients',k,' source',B(2*k-1,1),B(2*k,1), 'eddington',a11(k),a12(k),' :: ',a13(k),a23(k), ' :: ',a33(k), ' :: ',dtau(k),w0(k),g(k),mu0,'S=',S(k)
          call exit()
        endif
      enddo
    endif

    INFO=-1
    if(ireals.eq.REAL32) then !single_precision
      call SGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
    else if(ireals.eq.REAL64) then !double_precision
      call DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
    else
      call CHKERR(-5_mpiint, 'Dont know which LAPACK routine to call for real kind'//itoa(ireals))
    endif

    call CHKERR(int(INFO, mpiint), 'Error in twostream calculation - lapack returned Error')

    ! retrieve result from solver
    do k=1,ke1
      Eup(k) = B(2*k-1,NRHS) ! Eup
      Edn(k) = B(2*k,NRHS) ! Edn
    enddo
    if(ldebug) then
      do k=1,ke1
        if(any(isnan( [Eup(k),Edn(k)] )) ) &
          print *,'setting value for Eup,Edn',k,' LAPACK entries',B(2*k-1,1),B(2*k,1),'Eup/dn', Eup(k),Edn(k),'IPIV',IPIV(2*k-1:2*k)
      enddo
    endif

  end subroutine

  subroutine petsc_delta_eddington_twostream(dtau_in,w0_in,g_in,mu0,incSolar,albedo, S,Edn,Eup, planck)
    real(ireals),intent(in),dimension(:) :: dtau_in,w0_in,g_in
    real(ireals),intent(in) :: albedo,mu0,incSolar
    real(ireals),dimension(:),intent(out):: S,Edn,Eup
    real(ireals),dimension(:),intent(in),optional :: planck

    character(len=*), parameter :: default_options=''// &
      ' -twostream_ksp_type preonly'// &
      ' -twostream_pc_type lu'// &
      !' -twostream_ksp_view_mat'// &
      !' -twostream_ksp_view'// &
      !' -twostream_show_b'// &
      !' -twostream_show_x'// &
      !' -twostream_ksp_view_mat draw'// &
      ''

    real(ireals),dimension(size(dtau_in)) :: dtau,w0,g
    real(ireals),dimension(size(dtau_in)) :: a11,a12,a13,a23,a33,g1,g2

    integer(iintegers) :: i,j,k, ke, ke1, N
    integer(mpiint) :: ierr

    type(tMat), allocatable, save :: A
    type(tVec), allocatable, save :: b, x
    type(tKSP), allocatable, save :: ksp

    type(tVec) :: vEdn, vEup

    S=zero
    Edn=zero
    Eup=zero

    ke = size(dtau)
    ke1 = ke+1
    N = int(2*ke1, kind(N))

    dtau = dtau_in
    w0   = w0_in
    g    = g_in

    do k=1,ke
      call eddington_coeff_zdun (dtau(k), w0(k),g(k), mu0,a11(k),a12(k),a13(k),a23(k),a33(k), g1(k),g2(k) )
      if(ldebug) then
        if(any(isnan( [a11(k),a12(k),a13(k),a23(k),a33(k),g1(k),g2(k)] )) ) then
          print *,'eddington',k,' :: ',dtau(k), w0(k),g(k), mu0,'::',a11(k),a12(k),a13(k),a23(k),a33(k),'::',g1(k),g2(k)
          call exit()
        endif
      endif
    enddo

    S(1) = incSolar ! irradiance on tilted plane
    do k=1,ke
      S(k+1) = S(k) * a33(k)
    enddo

    if(.not.allocated(A)) then
      allocate(A)
      call MatCreateSeqAIJ(PETSC_COMM_SELF, N, N, 5_iintegers, PETSC_NULL_INTEGER, A, ierr); call CHKERR(ierr)
      call PetscObjectSetName(A, 'Twostream Mat', ierr);call CHKERR(ierr)
    endif
    if(.not.allocated(b)) then
      allocate(b)
      call VecCreateSeq(PETSC_COMM_SELF, N, b, ierr); call CHKERR(ierr)
      call VecSetBlockSize(b, i2, ierr); call CHKERR(ierr)
      call PetscObjectSetName(b, 'Twostream SRC Vec', ierr);call CHKERR(ierr)
    endif
    if(.not.allocated(x)) then
      allocate(x)
      call VecDuplicate(b, x, ierr); call CHKERR(ierr)
      call VecSetBlockSize(x, i2, ierr); call CHKERR(ierr)
      call PetscObjectSetName(x, 'Twostream Solution Vec', ierr);call CHKERR(ierr)
    endif

    if(.not.allocated(ksp)) then
      allocate(ksp)
      call KSPCreate(PETSC_COMM_SELF, ksp, ierr); call CHKERR(ierr)
      call KSPSetOptionsPrefix(ksp, "twostream_", ierr); call CHKERR(ierr)
      call PetscOptionsInsertString(PETSC_NULL_OPTIONS, default_options, ierr); call CHKERR(ierr)
      call KSPSetFromOptions(ksp, ierr); call CHKERR(ierr)
    endif


    ! Setup solar src vector
    do k=1,ke
      call VecSetValues(b, i2, [2*k-2, 2*k+1], [S(k) * a13(k), S(k) * a23(k)], INSERT_VALUES,ierr); call CHKERR(ierr)
    enddo
    call VecSetValues(b, i2, [i1, 2*ke1-2], [zero, S(ke1) * albedo], INSERT_VALUES,ierr); call CHKERR(ierr)

    call VecAssemblyBegin(b, ierr); call CHKERR(ierr)
    call VecAssemblyEnd(b, ierr); call CHKERR(ierr)

    ! Setup thermal src vector
    if(present(planck) ) then
      do k=1,ke
        call VecSetValues(b, i2, &
          [2*k-2, 2*k+1], &
          [(one-a11(k)-a12(k)) * planck(k) *pi, (one-a11(k)-a12(k)) * planck(k) *pi], &
          ADD_VALUES, ierr)
        call CHKERR(ierr)
      enddo
    endif

    call VecAssemblyBegin(b, ierr); call CHKERR(ierr)
    call VecAssemblyEnd(b, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(b, PETSC_NULL_VEC, '-twostream_show_b', ierr); call CHKERR(ierr)

    ! Set diagonal
    do i = 0, N-1
      j = i
      call MatSetValues(A, i1, i, i1, j, one, INSERT_VALUES, ierr); call CHKERR(ierr)
    enddo
    do k=1,ke
      ! setting Eup coeffs
      i = 2*k-2
      call MatSetValues(A, i1, i, i2, [i+2, i+1], [-a11(k), -a12(k)], INSERT_VALUES, ierr); call CHKERR(ierr)

      ! setting Edn coeffs
      i = 2*(k+1)-1
      call MatSetValues(A, i1, i, i2, [i-2, i-1], [-a11(k), -a12(k)], INSERT_VALUES, ierr); call CHKERR(ierr)
    enddo
    i = 2*ke1-2
    call MatSetValues(A, i1, i, i1, [i+1], [-albedo], INSERT_VALUES, ierr); call CHKERR(ierr)

    call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)
    call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr); call CHKERR(ierr)

    call KSPSetOperators(ksp, A, A, ierr); call CHKERR(ierr)
    call KSPSetFromOptions(ksp, ierr); call CHKERR(ierr)

    call KSPSolve(ksp, b, x, ierr); call CHKERR(ierr)

    call PetscObjectViewFromOptions(x, PETSC_NULL_VEC, '-twostream_show_x', ierr); call CHKERR(ierr)

    call VecCreateSeqWithArray(PETSC_COMM_SELF, i1, ke1, Eup, vEup, ierr); call CHKERR(ierr)
    call VecCreateSeqWithArray(PETSC_COMM_SELF, i1, ke1, Edn, vEdn, ierr); call CHKERR(ierr)

    call VecStrideGather(x, i0 , vEup, INSERT_VALUES, ierr); call CHKERR(ierr)
    call VecStrideGather(x, i1 , vEdn, INSERT_VALUES, ierr); call CHKERR(ierr)

    call  VecDestroy(vEdn, ierr); call CHKERR(ierr)
    call  VecDestroy(vEup, ierr); call CHKERR(ierr)
  end subroutine

  subroutine adding_delta_eddington_twostream(dtau, omega0, g, mu0, S0, Ag, Edir, Edn, Eup)
    real(ireals),intent(in),dimension(:) :: dtau, omega0, g
    real(ireals),intent(in) :: Ag, mu0, S0
    real(ireals),dimension(:),intent(out):: Edir, Edn, Eup

    integer(iintegers) :: k, ke, ke1

    real(ireals),dimension(size(dtau)) :: a11,a12,a13,a23,a33,g1,g2
    real(ireals),dimension(size(dtau)) :: T, R, Tdir, Sdir

    ke = size(dtau)
    ke1 = ke+1

    do k=1,ke
      call eddington_coeff_zdun (dtau(k), omega0(k),g(k), mu0, a11(k), a12(k), a13(k), a23(k), a33(k), g1(k), g2(k) )
      !call eddington_v2(dtau(k), omega0(k),g(k), mu0, a11(k), a12(k), a13(k), a23(k), a33(k))
      if(ldebug) then
        if(any(isnan( [a11(k),a12(k),a13(k),a23(k),a33(k),g1(k),g2(k)] )) ) then
          print *,'eddington',k,' :: ',dtau(k), omega0(k),g(k), mu0,'::',a11(k),a12(k),a13(k),a23(k),a33(k),'::',g1(k),g2(k)
          call exit()
        endif
      endif
    enddo

    Edir(1) = S0

    R(1) = a12(1)
    T(1) = a11(1)
    Tdir(1) = a33(1)
    Sdir(1) = a23(1)

    do k = 1, ke-1
      R(k+1) = a12(k+1) + (R(k) * a11(k+1) * a11(k+1)) / (one - R(k) * a12(k+1))
      T(k+1) = T(k) * a11(k+1) / (one - R(k) * a12(k+1))

      Tdir(k+1) = Tdir(k) * a33(k+1)
      Sdir(k+1) = (a11(k+1) * Sdir(k) + Tdir(k) * a13(k+1) * R(k) * a11(k+1)) / (one - R(k) * a12(k+1)) + Tdir(k) * a23(k+1)
    enddo

    do k = ke1, 2, -1
      Edir(k) = Tdir(k-1) * Edir(1)
    enddo

    Edn(ke1) = (Sdir(ke) + Tdir(ke) * R(ke) * Ag) / (one - R(ke) * Ag) * Edir(1)
    Eup(ke1) = Ag * (Edn(ke1) + Edir(ke1))

    do k = ke1, 3, -1
      Edn(k-1) = (R(k-2) * a11(k-1) * Eup(k) + Edir(1) * Sdir(k-2) + Edir(k-1) * a13(k-1) * R(k-2))  /  (one - R(k-2) * a12(k-1))
      Eup(k-1) = (a11(k-1) * Eup(k) + Edir(1) * Sdir(k-2) * a12(k-1) + Edir(k-1) * a13(k-1)) / (one - R(k-2) * a12(k-1))
    enddo
    Eup(1) = a11(1) * Eup(2) + a13(1) * Edir(1)


  end subroutine
!  subroutine eddington_bernhard_v2(dtau, g, omega0, mu0, &
!      t, r, rdir, sdir, tdir)
!    real(ireals),intent(in) :: dtau, g, omega0, mu0
!    real(ireals),intent(out) :: t, r, rdir, sdir, tdir
!    real(ireals) :: alpha1, alpha2, alpha3, alpha4, alpha5, alpha6
!    real(ireals) :: a11, a12, a13, a23, a33
!    real(ireals) :: lambda, b, A, denom
!
!    alpha1= (one-omega0)+0.75_ireals*(one-omega0*g)
!    alpha2=-(one-omega0)+0.75_ireals*(one-omega0*g)
!
!    lambda=sqrt(alpha1*alpha1-alpha2*alpha2)
!
!    A=one/(alpha2/(alpha1-lambda)*exp(lambda*dtau)-alpha2/(alpha1+lambda)*exp(-lambda*dtau))
!
!    a11=A*2.0_ireals*lambda/alpha2
!    a12=A*(exp(lambda*dtau)-exp(-lambda*dtau))
!
!    b=0.5-0.75*g*mu0
!    alpha3=-omega0*b
!    alpha4=omega0*(one-b)
!    denom = (one/mu0/mu0-lambda*lambda)
!    alpha5=((alpha1-one/mu0)*alpha3-alpha2*alpha4)/denom
!    alpha6=(alpha2*alpha3-(alpha1+one/mu0)*alpha4)/denom
!
!    a33=exp(-dtau/mu0)
!
!    a13=alpha5*(one-(a11)*(a33))-alpha6*(a12)
!    a23=-(a12)*alpha5*(a33)+alpha6*((a33)-(a11))
!
!    t    = a11
!    r    = a12
!    tdir = a33
!    rdir = a13 / mu0;
!    sdir = a23 / mu0;
!  end subroutine
end module
