module m_fitLUT
      use m_helper_functions, only: rmse
      use m_data_parameters, only: ireals,iintegers, init_mpi_data_parameters,numnodes,myid,zero,one,imp_comm
      use m_arrayio

      use m_optprop_LUT, only : t_optprop_LUT_8_10

      implicit none
#include "finclude/petsc.h90"

      PetscErrorCode :: ierr

type context
      integer(iintegers) :: xs,xe ! start row and end row of local vectors
      integer(iintegers) :: gN,lN ! global and local size
end type

      type(context) :: train_inp_ctx, train_out_ctx, test_inp_ctx, test_out_ctx, coeff_ctx

    contains
      subroutine iterative_fit(Norder,Ndim, test_inp,test_out, inp, out, coeff, luse_coeff)
        integer(iintegers),intent(in) :: Norder,Ndim
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
        do k=1,size(luse_coeff)*Norder**2

          call clean_small_coeffs(coeff,luse_coeff)

          call setup_iter_step(size(out), count(luse_coeff)+1, vcoeff,vtmp_coeff, A, ksp)

          jmin=0
          do j=1,size(luse_coeff)
            ltmp_coeff = luse_coeff
            if( ltmp_coeff(j) ) cycle ! already included
            ltmp_coeff(j) = .True.

            call fill_A(Norder,Ndim, size(out), ltmp_coeff, vinp, A)
            call solve_lsqr(vout, vtmp_coeff, A, ksp)
            call KSPGetResidualNorm(ksp,tmp_err,ierr)

            if( tmp_err.lt.err) then ! .and. mincoeff.gt. 1e-8_ireals) then ! if error was reduced -- we also want the coefficients to be big -- small coefficients are usually not stable?
              err=tmp_err
              jmin = j
              call VecCopy(vtmp_coeff,vcoeff,ierr)

              call VecGetArrayF90(vtmp_coeff,xcoeff,ierr)
              print *,'iteration',k,' :: selecting new coeff:',jmin,'with err',err/size(out),'::',xcoeff
              call VecRestoreArrayF90(vtmp_coeff,xcoeff,ierr)
            endif

          enddo

          if( jmin.eq.0 ) exit ! didnt find any coeff that could make it better
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

        print *,'Final coeffs:',luse_coeff,' :: ', coeff

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
          
        if(.not. allocated(coeff) ) return
        coeff_cnt = count(abs(coeff).gt.limit)

        if(coeff_cnt.ne.count(luse_coeff)) then
          allocate(tmp(size(coeff)),source=coeff)
          deallocate(coeff); allocate(coeff(coeff_cnt))
          itmp=1
          icoeff=1
          do i=1,size(luse_coeff)
            if(.not.luse_coeff(i) ) cycle
            if(abs(tmp(itmp)).gt.limit ) then 
              coeff(icoeff) = tmp(itmp)
              icoeff=icoeff+1
            else
              luse_coeff(i)=.False.
            endif
            itmp=itmp+1
          enddo
          print *,'coeffs after clean_small_coeffs:',luse_coeff
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

!        vinp      = createVec(size(inp),reshape( inp , (/size(inp)/) ) )
!        vout      = createVec(size(out),out)
!        vtest_inp = createVec(size(test_inp),reshape( test_inp, (/size(test_inp)/)) )
!        vtest_out = createVec(size(test_out),test_out)
      end subroutine

      subroutine fill_A(Norder,Ndim,Nsample,luse_coeff,vinp, A)
        integer(iintegers) :: Norder,Ndim,Nsample
        logical,intent(in) :: luse_coeff(:)
        Mat :: A
        Vec :: vinp
        integer(iintegers) :: icol,i,j
        PetscScalar, pointer :: xinp(:)
        real(ireals) :: v,e(2)

        call VecGetArrayF90(vinp,xinp,ierr)

        icol=0
        do j=1,size(luse_coeff)
          if(.not.luse_coeff(j)) cycle
          do i=1,Nsample
            e = exponent_poly_2d(Norder,j)
            v = xinp( Ndim*(i-1)+1 )**e(1) * xinp(Ndim*(i-1)+2)**e(2)
            call MatSetValues(A, 1_iintegers, [i-1], 1_iintegers, [icol], v, INSERT_VALUES, ierr)
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
!        vcoeff = createVec(Ncols)
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

      function poly2d(Norder, coeff_ind, coeff, inp)
        real(ireals) :: poly2d
        integer(iintegers),intent(in) :: Norder
        integer(iintegers),intent(in) :: coeff_ind(:)
        real(ireals),intent(in) :: coeff(:), inp(:)

        real(ireals) :: x(size(coeff)), exponent(2)
        integer(iintegers) :: j

        do j=1,size( coeff )
          exponent = exponent_poly_2d( Norder,coeff_ind(j) )
          x(j) = inp(1)**exponent(1) * inp(2)**exponent(2)
          print *,'poly2d for coeff',j,'got e:',exponent,'inp', inp
        enddo
        poly2d = dot_product(x, coeff) 
      end function
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
        poly1d = dot_product(x, coeff) 
      end function

      function createVec(lN,gN,src)
        Vec :: createVec
        integer(iintegers),intent(in) :: gN,lN
        real(ireals),optional :: src(:)

        PetscScalar, pointer :: xx_v(:)

        call VecCreate(PETSC_COMM_WORLD,createVec,ierr);
        call VecSetFromOptions(createVec,ierr)
        call VecSetSizes(createVec,lN,gN,ierr);

        if(present(src)) then
          call VecGetArrayF90(createVec,xx_v,ierr)
          xx_v = src
          call VecRestoreArrayF90(createVec,xx_v,ierr)
        endif
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
!        print *,'exp poly 2d ::',Norder,index,'==>',exponent_poly_2d
      end function

      subroutine fitLUT(dx,azimuth,zenith)
        real(ireals),intent(in) :: dx,azimuth,zenith
      
        type(t_optprop_LUT_8_10) :: OPP  

        call OPP%init(dx,dx,[azimuth],[zenith],imp_comm)
      end subroutine
      subroutine setup_context(v,ctx)
        Vec,intent(in) :: v
        type(context),intent(out) :: ctx
        call VecGetOwnershipRange(v, ctx%xs, ctx%xe, ierr)
        call VecGetSize(v,ctx%gN,ierr)
        call VecGetLocalSize(v,ctx%lN,ierr)
        print *,myid,'Setting up context with local indices:',ctx%xs,ctx%xe,' and sizes', ctx%lN,ctx%gN
      end subroutine
end module

program main
      use m_fitLUT 
      implicit none

      integer(iintegers),parameter :: Ndim=2,Nsample=5,Ntest=Nsample*2           ! Size of original data (here created by test_func)
      integer(iintegers),parameter :: Norder=2                                    ! Order of the polynomial that we want to fit
      Vec :: inp,out
      Vec :: test_inp,test_out
      Vec :: fit_out
      Vec :: coeff,coeff_ind
      
      integer(iintegers) :: i

      character(len=300) :: groups(3)

      PetscRandom :: rctx
      PetscScalar, pointer :: xv_coeff_ind(:)
      PetscScalar, pointer :: xv_coeff(:),xv_inp(:),xv_out(:)

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr) ;CHKERRQ(ierr)
      call init_mpi_data_parameters(PETSC_COMM_WORLD)

!      call fitLUT(40._ireals,zero,zero)
      
      out = createVec(PETSC_DECIDE,Nsample)      
      call setup_context(out, train_out_ctx)

      inp = createVec(train_out_ctx%lN*Ndim ,Nsample*Ndim)
      call setup_context(inp, train_inp_ctx)

      test_out = createVec(PETSC_DECIDE, Ntest)      
      call setup_context(test_out, test_out_ctx)

      test_inp = createVec(test_out_ctx%lN*Ndim, Ntest*Ndim)
      call setup_context(test_inp, test_inp_ctx)

      fit_out  = createVec(PETSC_DECIDE, Ntest)      

      coeff     = createVec(Norder**Ndim,PETSC_DECIDE)
      call setup_context(coeff, coeff_ctx)
      coeff_ind = createVec(Norder**Ndim, PETSC_DECIDE)

      call PetscRandomCreate(PETSC_COMM_WORLD,rctx,ierr)
      call PetscRandomSetFromOptions(rctx,ierr)
      call VecSetRandom(coeff,rctx,ierr)
      call VecSet(coeff,one,ierr)

      call VecGetArrayF90(coeff_ind, xv_coeff_ind, ierr)
      call VecGetArrayF90(coeff, xv_coeff, ierr)
      do i=1,size(xv_coeff_ind)
        xv_coeff_ind(i) = i
      enddo

      !Setup initial input and output for training dataset
      call VecGetArrayF90(  inp, xv_inp,   ierr)
      call VecGetArrayF90(  out, xv_out,   ierr)

      do i=1,train_out_ctx%lN
        xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim ) = [ one*i+train_out_ctx%xs, one ]
        xv_out(i) = poly2d(Norder, int(xv_coeff_ind), xv_coeff, xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim  ) ) !*(one+((xv_out(i)-.5)*40/100_ireals)) 
        print *,myid,'setting train_input for',i,'::',(i-1)*Ndim + 1 , (i-1)*Ndim + Ndim ,' ::: ',xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim  ) 
      enddo
      call VecRestoreArrayF90(  inp, xv_inp,   ierr)
      call VecRestoreArrayF90(  out, xv_out,   ierr)


      !Additionally setup test function
!      call VecGetArrayF90(  test_inp, xv_inp,   ierr)
!      call VecGetArrayF90(  test_out, xv_out,   ierr)
!      do i=1,test_out_ctx%lN
!        xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim ) = [ one*i+test_out_ctx%xs, one ] * [ one*Nsample/Ntest, one ] ! double sampling along x
!        xv_out(i) = poly2d(Norder, int(xv_coeff_ind), xv_coeff, xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim  ) ) !*(one+((xv_out(i)-.5)*40/100_ireals)) 
!      enddo
!      call VecRestoreArrayF90(  inp, xv_inp,   ierr)
!      call VecRestoreArrayF90(  out, xv_out,   ierr)


      call VecRestoreArrayF90(coeff_ind, xv_coeff_ind, ierr)
      call VecRestoreArrayF90(coeff, xv_coeff, ierr)

      call VecView(inp, PETSC_VIEWER_STDOUT_WORLD ,ierr)
      call VecView(out, PETSC_VIEWER_STDOUT_WORLD ,ierr)
!      call VecView(test_inp, PETSC_VIEWER_STDOUT_WORLD ,ierr)
!      call VecView(test_out, PETSC_VIEWER_STDOUT_WORLD ,ierr)

      if(myid.eq.0) print *,'Initial coefficients, used for test dataset'
      call VecView(coeff, PETSC_VIEWER_STDOUT_WORLD ,ierr)
      call VecView(coeff_ind, PETSC_VIEWER_STDOUT_WORLD ,ierr)

!      call iterative_fit(Norder,Ndim, test_inp,test_out,inp,out,coeff,luse_coeff)
!      print *,'coeff of fit: ',coeff
!
!      do i=1,Ntest
!        fit(i) = poly2d(Norder, luse_coeff, coeff, test_inp(:,i) )
!!        print *,' inp ',nint(test_inp(:,i)),' out ',nint(test_out(i)),' :: fit ',nint(fit(i))
!      enddo
!
!      print *,'RMSE:',rmse(test_out,fit)
!
!      groups(1) = "./fit.out"; groups(2)="/fit"
!      groups(3)="out";   call h5write(groups,fit,ierr)
!      groups(3)="x";     call h5write(groups,test_inp(1,:),ierr)
!      groups(3)="y";     call h5write(groups,test_inp(2,:),ierr)
!      groups(3)='coeff'; call h5write(groups,coeff,ierr)
!
!      groups(1) = "./fit.out"; groups(2)="/test"
!      groups(3)="out";   call h5write(groups,test_out,ierr)
!      groups(3)="x";     call h5write(groups,test_inp(1,:),ierr)
!      groups(3)="y";     call h5write(groups,test_inp(2,:),ierr)
!
!      groups(1) = "./fit.out"; groups(2)="/train"
!      groups(3)="out";   call h5write(groups,out,ierr)
!      groups(3)="x";     call h5write(groups,inp(1,:),ierr)
!      groups(3)="y";     call h5write(groups,inp(2,:),ierr)

      call PetscFinalize(ierr) ;CHKERRQ(ierr)
end program
