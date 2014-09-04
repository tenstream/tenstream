module m_fitLUT
      use m_helper_functions, only: rmse
      use m_data_parameters, only: ireals,iintegers, init_mpi_data_parameters,numnodes,myid,zero,one,imp_comm

      implicit none
#include "finclude/petsc.h90"

      PetscErrorCode :: ierr

    contains
      subroutine iterative_fit(Norder, test_inp,test_out, inp, out, coeff, luse_coeff)
        integer(iintegers),intent(in) :: Norder
        real(ireals),intent(in) :: inp(:,:),out(:)
        real(ireals),intent(in) :: test_inp(:,:),test_out(:)
        real(ireals),intent(out),allocatable :: coeff(:)
        logical,intent(inout) :: luse_coeff(:)

        integer(iintegers) :: j,k,jmin
        Mat :: A
        Vec :: vcoeff
        Vec :: vinp,vout
        Vec :: vtest_inp,vtest_out
        Vec :: vtmp_coeff
        KSP :: ksp

        real(ireals) :: err,tmp_err
        PetscScalar, pointer :: xcoeff(:)

        logical :: ltmp_coeff(size(luse_coeff) )
        luse_coeff=.False.
        err=huge(err)

        call setup_general(test_inp,test_out,vtest_inp,vtest_out,inp,out,vinp,vout)

        jmin=0
        do k=1,10

          call clean_small_coeffs(coeff,luse_coeff)

          call setup_iter_step(size(inp), count(luse_coeff)+1, vcoeff,vtmp_coeff, A, ksp)

          do j=1,size(luse_coeff)
            ltmp_coeff = luse_coeff
            if( ltmp_coeff(j) ) cycle ! already included
            ltmp_coeff(j) = .True.

            call fill_A(size(inp), ltmp_coeff, vinp, A)
            call solve_lsqr(vout, vtmp_coeff, A, ksp)
            call KSPGetResidualNorm(ksp,tmp_err,ierr)

            call VecGetArrayF90(vtmp_coeff,xcoeff,ierr)
            print *,'error was ',err,'is now ',tmp_err,'at outer loop',k,' iteration',j,'coeff',xcoeff
            call VecRestoreArrayF90(vtmp_coeff,xcoeff,ierr)


            if( tmp_err.lt.err) then ! .and. mincoeff.gt. 1e-8_ireals) then ! if error was reduced -- we also want the coefficients to be big -- small coefficients are usually not stable?
              err=tmp_err
              jmin = j
              call VecCopy(vtmp_coeff,vcoeff,ierr)
            endif

          enddo

          if( luse_coeff(jmin) ) exit ! if the best coeff was already in use, that means we couldnt find a better one, we can safely stop the search...

          luse_coeff(jmin) =.True. 

          if(allocated(coeff) ) deallocate(coeff)
          allocate(coeff(count(luse_coeff)))

          call VecGetArrayF90(vcoeff,xcoeff,ierr)
          coeff = xcoeff
          call VecRestoreArrayF90(vcoeff,xcoeff,ierr)

          call KSPDestroy(ksp,ierr)
          call MatDestroy(A,ierr)
          call VecDestroy(vcoeff,ierr)
          call VecDestroy(vtmp_coeff,ierr)

          if( check_if_converged(err,coeff) ) exit  
        enddo

        print *,'will use coeffs:',luse_coeff,' :: ', coeff

      end subroutine
      function check_if_converged(err,coeff)
        logical :: check_if_converged
!        Vec :: vcoeff
        real(ireals) :: err,coeff(:)
!        integer(iintegers) :: minloc
!        real(ireals) :: mincoeff

        check_if_converged=.True.

        if( err.gt.1e-4_ireals ) check_if_converged=.False.
!        call VecMin(vcoeff, PETSC_NULL_INTEGER, mincoeff, ierr )
        if( minval(coeff).lt.1e-8_ireals) check_if_converged=.False.
        
      end
      subroutine clean_small_coeffs(coeff,luse_coeff)
        logical,intent(inout) :: luse_coeff(:)
        real(ireals),allocatable :: coeff(:)
        real(ireals),allocatable :: tmp(:)
        integer(iintegers) :: i,itmp,icoeff,coeff_cnt
        real(ireals),parameter :: limit=1e-8

        coeff_cnt = count(coeff.gt.limit)

        if(coeff_cnt.ne.count(luse_coeff)) then
          allocate(tmp(size(coeff)),source=coeff)
          deallocate(coeff); allocate(coeff(coeff_cnt))
          itmp=1
          icoeff=1
          do i=1,size(luse_coeff)
            if(.not.luse_coeff(i) ) cycle
            if(tmp(itmp).gt.limit ) then 
              coeff(icoeff) = tmp(itmp)
              icoeff=icoeff+1
            else
              luse_coeff(i)=.False.
            endif
            itmp=itmp+1
          enddo
          print *,'will use coeff:',luse_coeff
        endif

      end subroutine

      subroutine solve_lsqr(vout, vcoeff, A, ksp)
        Mat :: A
        KSP :: ksp
        Vec :: vout,vcoeff

        ! Set KSP solver
        call KSPSetOperators(ksp, A, A, ierr)
        call KSPSetType(ksp, KSPLSQR, ierr)
        call KSPSetFromOptions(ksp, ierr)
        call KSPSetUp(ksp, ierr)
        call KSPSolve(ksp, vout, vcoeff, ierr)
      end subroutine

      subroutine setup_general(test_inp,test_out,vtest_inp,vtest_out,inp,out,vinp,vout)
        real(ireals),intent(in) :: test_inp(:,:),test_out(:)
        Vec :: vtest_inp, vtest_out
        real(ireals),intent(in) :: inp(:,:),out(:)
        Vec :: vinp, vout

        vinp      = createVec(size(inp),reshape( inp , (/size(inp)/) ) )
        vout      = createVec(size(out),out)
        vtest_inp = createVec(size(test_inp),reshape( test_inp, (/size(test_inp)/)) )
        vtest_out = createVec(size(test_out),test_out)
      end subroutine

      subroutine fill_A(Nsample,luse_coeff,vinp, A)
        integer(iintegers) :: Nsample
        logical,intent(in) :: luse_coeff(:)
        Mat :: A
        Vec :: vinp
        integer(iintegers) :: icol,i,j
        PetscScalar, pointer :: xinp(:)

        call VecGetArrayF90(vinp,xinp,ierr)

        icol=0
        do j=1,size(luse_coeff)
          if(.not.luse_coeff(j)) cycle
          do i=1,Nsample
            call MatSetValues(A, 1_iintegers, [i-1], 1_iintegers, [icol], xinp(i)**(j-1), INSERT_VALUES, ierr)
          enddo
          icol=icol+1
        enddo

        call MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY, ierr)
        call VecRestoreArrayF90(vinp,xinp,ierr)
        call MatAssemblyEnd  ( A, MAT_FINAL_ASSEMBLY, ierr)
      end subroutine

      subroutine setup_iter_step(Nrows,Ncols,vcoeff,vtmp_coeff,A,ksp)
        integer(iintegers),intent(in) :: Nrows,Ncols
        Vec :: vcoeff,vtmp_coeff
        Mat :: A
        KSP :: ksp
        call KSPCreate(imp_comm, ksp, ierr)
        call init_Mat(A,Nrows, Ncols )
        vcoeff = createVec(Ncols)
        call VecDuplicate(vcoeff,vtmp_coeff,ierr)
      end subroutine

      subroutine init_Mat(A,Nrow,Ncol)
        integer(iintegers),intent(in) :: Nrow,Ncol
        Mat :: A
        print *,'Creating Matrix with ',Nrow,'rows and ',Ncol,'columns'
        call MatCreate(imp_comm, A, ierr)
        call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, Nrow, Ncol, ierr)
        call MatSetUp(A, ierr)
      end subroutine

      function poly1d(Norder, luse_coeff, coeff, inp)
        real(ireals) :: poly1d
        integer(iintegers),intent(in) :: Norder
        logical,intent(in) :: luse_coeff(:)
        real(ireals),intent(in) :: coeff(:), inp(:)
        real(ireals) :: x(size(coeff)), exponent(1)

        integer(iintegers) :: j,icoeff

        poly1d = 0

        icoeff=1
        do j=1,size(luse_coeff)
          if(.not.luse_coeff(j) ) cycle
          exponent = exponent_poly_1d(Norder,j)
          x(icoeff) = inp(1)**exponent(1)
          icoeff=icoeff+1
        enddo
!        poly1d = poly1d + coeff(icoeff)* x**(j-1) 
        poly1d = dot_product(x, coeff) 
      end function

      function createVec(N,src)
        Vec :: createVec
        integer(iintegers),intent(in) :: N
        real(ireals),optional :: src(:)

        PetscScalar, pointer :: xx_v(:)

        call VecCreate(PETSC_COMM_WORLD,createVec,ierr);
        call VecSetFromOptions(createVec,ierr)
        call VecSetSizes(createVec,PETSC_DECIDE,N,ierr);

        if(present(src)) then
          call VecGetArrayF90(createVec,xx_v,ierr)
          xx_v = src
          call VecRestoreArrayF90(createVec,xx_v,ierr)
        endif
      end function

      function test_func( X, dim, order)
        real(ireals) :: test_func
        integer(iintegers) :: dim,order
        real(ireals) :: X(dim)
        integer(iintegers) :: i,j
        real(ireals) :: t

!        t=0
!        do i=1,order+1
!!          do j=1,order+1-i
!!            do k=1,order-j
!              t = t + p( X(1),i ) * p( X(2),j ) !* p( X(3),k )
!              print *,t,'adding:',p( X(1),i ), p( X(2),j ),':',p( X(1),i ) * p( X(2),j )
!!            enddo
!!          enddo
!        enddo

        test_func = X(1)**2 + 5
      end function
      function p(x,e)
        real(ireals) :: p
        real(ireals),intent(in) :: x
        integer(iintegers),intent(in) :: e
        p=x**e
      end function
      function index_poly_1d(Norder,exponent)
        integer(iintegers) :: index_poly_1d
        integer(iintegers),intent(in) :: exponent(1),Norder
        index_poly_1d = exponent(1) +1_iintegers
      end function
      function exponent_poly_1d(Norder,index)
        integer(iintegers) :: exponent_poly_1d(1)
        integer(iintegers),intent(in) :: index,Norder
        exponent_poly_1d(1) = modulo(index,Norder+1_iintegers) -1
      end function

      function index_poly_2d(Norder,exponent)
        integer(iintegers) :: index_poly_2d
        integer(iintegers),intent(in) :: exponent(2),Norder
        index_poly_2d = exponent(1)*(Norder+1_iintegers) + exponent(2) +1_iintegers
      end function
      function exponent_poly_2d(Norder,index)
        integer(iintegers) :: exponent_poly_2d(2)
        integer(iintegers),intent(in) :: index,Norder
        exponent_poly_2d(2) = modulo(index,Norder+1_iintegers)
        exponent_poly_2d(1) = ( index-exponent_poly_2d(2) ) / ( Norder+1 )
      end function

end module

program main
      use m_fitLUT 
      implicit none

      integer(iintegers),parameter :: Ndim=1,Nsample=15,Ntest=Nsample*100         ! Size of original data (here created by test_func)
      integer(iintegers),parameter :: Norder=6                                    ! Order of the polynomial that we want to fit
      real(ireals) :: inp(Ndim,Nsample), out(Nsample)                             ! inp and out are the coordinates and the result of test_func
      real(ireals) :: test_inp(Ndim,Ntest), test_out(Ntest)                       ! inp and out are the coordinates and the result of test_func
      real(ireals) :: fit(Nsample)                                                ! result of fit function                       

      real(ireals),allocatable :: coeff(:)                                             ! coefficients for polynomial
      logical :: luse_coeff(Norder**Ndim)

      integer(iintegers) :: i

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr) ;CHKERRQ(ierr)
      call init_mpi_data_parameters(PETSC_COMM_WORLD)

      do i=1,Nsample
        inp(:,i) = i
        out(i) = test_func( inp(:,i), Ndim, 1_iintegers )
      enddo

      do i=1,Ntest
        test_inp(:,i) = i * Nsample/Ntest
        test_out(i) = test_func( test_inp(:,i), Ndim, 1_iintegers)
      enddo

      call iterative_fit(Norder, test_inp,test_out,inp,out,coeff,luse_coeff)

!      print *,'coeff of fit: ',coeff
      do i=1,Nsample
        fit(i) = poly1d(Norder, luse_coeff, coeff, inp(:,i) )
        print *,' inp ',nint(inp(:,i)),' out ',nint(out(i)),' :: fit ',nint(fit(i))
      enddo


      call PetscFinalize(ierr) ;CHKERRQ(ierr)
end program
