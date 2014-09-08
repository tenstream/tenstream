module m_fitLUT
      use m_helper_functions, only: rmse
      use m_data_parameters, only: ireals,iintegers, init_mpi_data_parameters,numnodes,myid,zero,one,imp_comm
      use m_arrayio

      use m_optprop_LUT, only : t_optprop_LUT_8_10

      implicit none
#include "finclude/petsc.h90"

      PetscErrorCode :: ierr
      PetscViewer :: view

type Matrix
      integer(iintegers) :: xs,xe ! start row and end row of local portion
      integer(iintegers) :: lNrow,lNcol ! local size
      integer(iintegers) :: gNrow,gNcol ! global size
      Mat :: A
end type
type Vector
      integer(iintegers) :: xs,xe ! start row and end row of local vectors
      integer(iintegers) :: gN,lN ! global and local size
      Vec :: v
end type
      integer(iintegers) :: poly_func
      integer(iintegers) :: poly_order
      real(ireals) :: err_ratio

    contains
      subroutine iterative_fit(Norder,Ndim, train_inp, train_out, luse_coeff, coeff)
        integer(iintegers),intent(in) :: Norder,Ndim
        type(Vector) :: train_inp, train_out, coeff
        logical,intent(out) :: luse_coeff(:)

        logical :: ltmp_coeff(size(luse_coeff))

        type(Matrix) :: A
        type(Vector) :: tmp_coeff,tmp_out

        KSP :: ksp

        real(ireals) :: err=huge(err),tmp_err=huge(err),last_err

        integer(iintegers) :: k,j, jmin

      if(myid.eq.0) print *,'Maximum nr of polynomials for dimension/order',Ndim,Norder,' is',polygonal_number(Ndim,Norder)

        call VecSet(coeff%v         , zero)

        call createVec(coeff%lN     , coeff%gN     , tmp_coeff)
        call createVec(train_out%lN , train_out%gN , tmp_out)

        call KSPCreate(imp_comm, ksp, ierr)

        luse_coeff=.False. ! start without any coefficients
        call init_Mat(A,train_out%lN,coeff%lN, train_out%gN,coeff%gN)

        k=0
        do ! do this a while, inifinite would be nice, can we find a suitable exit condition?
          k=k+1

          jmin=0; 
          last_err = err
          do j=1,size(luse_coeff)
            ltmp_coeff = luse_coeff
            if( ltmp_coeff(j) ) cycle ! already included
            ltmp_coeff(j) = .True.

            call MatZeroEntries(A%A,ierr)
            call fill_A(Norder,Ndim, ltmp_coeff,train_inp, A)

            call solve_lsqr(train_out, tmp_coeff , A, ksp)
            call KSPGetResidualNorm(ksp,tmp_err,ierr)

!            call MatMult(A%A,coeff%v, tmp_out%v,ierr)
!            call VecView( tmp_out%v,PETSC_VIEWER_STDOUT_WORLD,ierr)

            if(myid.eq.0) print *,'k=',k,' :: testing new coeff:',j,'(',size(luse_coeff),')','old err',err,err/train_inp%gN,'with new err',tmp_err,tmp_err/train_inp%gN ,'::',tmp_err/err
            if( tmp_err.lt.err) then ! .and. mincoeff.gt. 1e-8_ireals) then ! if error was reduced -- we also want the coefficients to be big -- small coefficients are usually not stable?
              err=tmp_err
              jmin = j
              call VecSet(coeff%v,zero,ierr)
              call VecAXPY(coeff%v, one, tmp_coeff%v,ierr)
              call VecCopy(tmp_coeff%v,coeff%v,ierr)

!              call VecView( coeff%v,PETSC_VIEWER_STDOUT_WORLD,ierr)
            endif

          enddo

          if( err/last_err .gt. err_ratio ) then
            if(myid.eq.0) print *,'exiting iterative fit routine :: last coefficient we added reduce the error by less than would like to ',err,last_err,' ratio is:',err/last_err
            exit
          endif
          if( jmin.eq.0 ) then
            if(myid.eq.0) print *,'exiting iterative fit routine :: didnt find any coeff that could make it better'
            exit
          endif

          luse_coeff(jmin) =.True. 
          if(myid.eq.0) print *,'iteration',k,' :: selecting new coeff:',jmin,'with err',err,err/train_inp%gN  !,'::',luse_coeff

          if( check_if_converged(err/train_inp%gN) ) then
            if(myid.eq.0) print *,'exiting iterative fit routine :: we converged'
            exit
          endif
        enddo

        call MatDestroy(A%A,ierr)
        call KSPDestroy(ksp,ierr)
        call VecDestroy(tmp_coeff%v,ierr)

        if(myid.eq.0) print *,'Final coeffs:',luse_coeff
        call VecView( coeff%v,PETSC_VIEWER_STDOUT_WORLD,ierr)

      end subroutine
      function check_if_converged(err)
        logical :: check_if_converged
        real(ireals) :: err

        check_if_converged=.True.

        if( err.gt.1e-6_ireals ) check_if_converged=.False.
!        call VecMin(vcoeff, PETSC_NULL_INTEGER, mincoeff, ierr )
!        if( minval(coeff).lt.1e-8_ireals) check_if_converged=.False.
        
      end

      subroutine solve_lsqr(vout, vcoeff, A, ksp)
        type(Matrix) :: A
        type(Vector) :: vout,vcoeff
        KSP :: ksp

        ! Set KSP solver
        call KSPSetOperators(ksp, A%A, A%A, ierr)
        call KSPSetType(ksp, KSPLSQR, ierr)
        call KSPSetFromOptions(ksp, ierr)
        call KSPSetUp(ksp, ierr)
        call KSPSolve(ksp, vout%v, vcoeff%v, ierr)
      end subroutine

      subroutine fill_A(Norder,Ndim, luse_coeff,inp, A)
        integer(iintegers) :: Norder,Ndim
        type(Matrix) :: A
        type(Vector) :: inp
        logical,intent(in) :: luse_coeff(A%gNcol)

        integer(iintegers) :: icol,irow
        integer(iintegers) :: col_ind(A%gNcol)
        real(ireals) :: xvals(A%gNcol)
        real(ireals) :: xvals_inp(Ndim)

        PetscScalar, pointer :: xinp(:)

        col_ind = [ (icol-1, icol=1,A%gNcol)  ] 

        call VecGetArrayF90(inp%v,xinp,ierr)
        do irow=A%xs,A%xe-1
          xvals_inp = xinp( (irow-A%xs)*Ndim+1 : (irow-A%xs)*Ndim + Ndim )
          xvals = poly(Norder,Ndim, luse_coeff, xvals_inp)

!          print *,irow,':',irow-A%xs,'xvals_inp',xvals_inp,'xvals',xvals,'A%gNcol',A%gNcol,'col_ind',col_ind
!          if(any(isnan(xvals)) ) stop 'fill_A :: NaN in xvals -- bug!'
          call MatSetValues(A%A, 1_iintegers, [irow], A%gNcol, col_ind, xvals, INSERT_VALUES, ierr)
        enddo
        call VecRestoreArrayF90(inp%v,xinp,ierr)

        call MatAssemblyBegin( A%A, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEnd  ( A%A, MAT_FINAL_ASSEMBLY, ierr)
      end subroutine

      subroutine init_Mat(A,lNrow,lNcol, gNrow,gNcol)
        integer(iintegers),intent(in) :: lNrow,lNcol, gNrow,gNcol
        type(Matrix) :: A
        call MatCreate(imp_comm, A%A, ierr)
        call MatSetFromOptions(A%A,ierr)
        call MatSetSizes(A%A, lNrow,lNcol, gNrow,gNcol, ierr)
        call MatMPIAIJSetPreallocation(A%A, lNcol,PETSC_NULL_INTEGER, gNcol-lNcol,PETSC_NULL_INTEGER,ierr);CHKERRQ(ierr)
        call MatSeqAIJSetPreallocation(A%A, lNcol, PETSC_NULL_INTEGER, ierr) ;CHKERRQ(ierr)
        call MatSetUp(A%A, ierr)
        call MatGetOwnershipRange(A%A,A%xs,A%xe,ierr)
        call MatGetLocalSize(A%A,A%lNrow,A%lNcol,ierr)
        call MatGetSize(A%A, A%gNrow,A%gNcol,ierr)
!        print *,myid,'local portion is from:',A%xs,A%xe
!        print *,myid,'local rows and columns :',A%lNrow,A%lNcol
      end subroutine
 
      subroutine createVec(lN,gN,V)
        type(Vector) :: V
        integer(iintegers),intent(in) :: gN,lN

        call VecCreate(PETSC_COMM_WORLD,V%v,ierr);
        call VecSetFromOptions(V%v,ierr)
        call VecSetSizes(V%v,lN,gN,ierr);

        call VecGetOwnershipRange(V%v, V%xs, V%xe, ierr)
        call VecGetSize(V%v,V%gN,ierr)
        call VecGetLocalSize(V%v,V%lN,ierr)
        
        !print *,myid,'Setting up context with local indices:',V%xs,V%xe,' and sizes', V%lN,V%gN
      end subroutine

       function poly(Norder,Ndim, luse_coeff, inp, print_func,coeff)
        logical,intent(in) :: luse_coeff(:)
        real(ireals) :: poly(size(luse_coeff))
        integer(iintegers),intent(in) :: Norder,Ndim
        real(ireals),intent(in) :: inp(Ndim)

        character(len=*),optional :: print_func
        real(ireals),optional :: coeff(size(luse_coeff))

        integer(iintegers) :: j,xd,yd,zd,kd

        if(present(print_func) ) write( *,FMT='(A)', advance='no') trim(print_func)
        j=1
        select case (Ndim)
        case (4)
          do xd=0,Norder
            do yd=0,Norder-xd
              if(xd+yd.gt.Norder) cycle
              do zd=0,Norder-yd
                if(xd+yd+zd.gt.Norder) cycle
                do kd=0,Norder-zd
                  if(xd+yd+zd+kd.gt.Norder) cycle

!                    print *,'poly:',xd,yd,zd,kd,'::',j

                  if(luse_coeff(j) ) then

                    if(present(print_func) ) then
                      write( *,FMT='( "+",E14.7,"* ( x**",I0,"* y**",I0,"* z**",I0,"* k**",I0," ) ")',advance='no') coeff(j),xd,yd,zd,kd
                    else
                      poly(j) = p(inp(1),xd) * p(inp(2),yd) * p(inp(3),zd) * p(inp(4),kd)
                    endif

                  else
                    poly(j) = PETSC_NULL_REAL
                  endif
                  j=j+1
                enddo
              enddo
            enddo
          enddo

        case (2)
          do xd=0,Norder-1
            do yd=0,Norder-1-xd
              if(xd+yd.ge.Norder) cycle
              if(luse_coeff(j) ) then
                poly(j) = p(inp(1),xd) * p(inp(2),yd)
              else
                poly(j) = PETSC_NULL_REAL
              endif
              j=j+1
            enddo
          enddo

        case default
          print *,'did not implement poly function for ',Ndim,'dimensions'
          call mpi_abort(imp_comm,ierr)
        end select

        if(present(print_func) ) print *,''
      end function 
      pure function p(v,e) 
        real(ireals) :: p
        real(ireals),intent(in) :: v
        integer(iintegers),intent(in) :: e
        select case(poly_func)
        case(1)
          p=v**e
        case(2)
          p = legendre(v,e)
        case(3)
          p = exp(v*e)
        case default
!          stop 'dont know poly_func option'
          p = legendre(v,e)
        end select
      end function

      pure function legendre(x,n)
        real(ireals) :: legendre
        real(ireals),intent(in) :: x
        integer(iintegers),intent(in) :: n
        real(ireals),dimension(n+1) :: pn,pd
      !http://de.wikipedia.org/wiki/Legendre-Polynom
!        legendre=0
!        do k=0,n/2
!          legendre = legendre + (one)**k * factorial(2*n-2*k) /factorial(n-k)/factorial(n-2*k)/factorial(k)/2**n * x**(n-2*k)
!        enddo
        call LPN(n,2._ireals*x-one, pn,pd)
        legendre=pn(n+1)
      end function
      pure SUBROUTINE LPN(N,X,PN,PD)
!
!       ===============================================
!       Purpose: Compute Legendre polynomials Pn(x)
!                and their derivatives Pn'(x)
!       Input :  x --- Argument of Pn(x)
!                n --- Degree of Pn(x) ( n = 0,1,...)
!       Output:  PN(n) --- Pn(x)
!                PD(n) --- Pn'(x)
!       ===============================================
!
        real(ireals),intent(in) :: x
        integer(iintegers),intent(in) :: N
        real(ireals),intent(out),DIMENSION(0:N) :: pn,pd
        integer(iintegers) :: k
        real(ireals) :: p0,p1,pf
        PN(0)=one
        PN(1)=X
        PD(0)=zero
        PD(1)=one
        P0=one
        P1=X
        DO 10 K=2,N
           PF=(2._ireals*K-one)/K*X*P1-(K-one)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.one) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           P0=P1
10         P1=PF
        RETURN
        END

      pure recursive function factorial(n)  result(fact)
        integer(iintegers) :: fact
        integer(iintegers), intent(in) :: n
        if (n == 0) then
          fact = 1
        else
          fact = n * factorial(n-1)
        endif 
      end function          
      pure function polygonal_number(dim,vars)
        !http://en.wikipedia.org/wiki/Polygonal_number 
        integer(iintegers) :: polygonal_number
        integer(iintegers),intent(in) :: dim,vars 
        real(ireals) :: tmp
        integer(iintegers) :: i
        tmp=1
        do i=1,vars
          tmp = tmp * ((dim+vars)+one-i)/i
!          tmp = tmp + (s-2_ireals)* ( i*(i-1) *.5_ireals ) +i
        enddo
        polygonal_number = nint(tmp)
      end function
      subroutine write_poly_fit(inp,out,fit,fname,prefix)
        type(Vector) :: inp,out,fit
        character(len=*) :: fname,prefix
        character(300) :: vecname

        PetscFileMode :: fmode
        logical :: fexists
        
        if(myid.eq.0) print *,'**************************** Writing Solution **********************************'
        write( vecname,FMT='(A,".poly",I0,".order",I0,".")' ) trim(prefix), poly_func, poly_order

        inquire(file=trim(fname), exist=fexists)
        if(fexists) then
          if(myid.eq.0)  print *,myid,'appending vector to hdf5 file ',trim(fname),' vecname: ',vecname
          fmode = FILE_MODE_APPEND
        else 
          if(myid.eq.0)  print *,myid,'writing vector to hdf5 file ',trim(fname),' vecname: ',vecname
          fmode = FILE_MODE_WRITE
        endif

        call PetscViewerHDF5Open(imp_comm,fname, fmode, view, ierr) ;CHKERRQ(ierr)
        call PetscObjectSetName(inp%v, trim(vecname)//'inp',ierr) ; CHKERRQ(ierr);   call VecView(inp%v, view, ierr) ;CHKERRQ(ierr)
        call PetscObjectSetName(out%v, trim(vecname)//'out',ierr) ; CHKERRQ(ierr);   call VecView(out%v, view, ierr) ;CHKERRQ(ierr)
        call PetscObjectSetName(fit%v, trim(vecname)//'fit',ierr) ; CHKERRQ(ierr);   call VecView(fit%v, view, ierr) ;CHKERRQ(ierr)
        call PetscViewerDestroy(view,ierr) ;CHKERRQ(ierr)
      end subroutine

      subroutine fitLUT(dx,azimuth,zenith)
        real(ireals),intent(in) :: dx,azimuth,zenith

!        integer(iintegers),parameter :: Norder=8
        type(t_optprop_LUT_8_10) :: OPP  
        integer(iintegers) :: Ndim,Nsample,Ncoeff

        type(Vector) :: train_inp, train_out
        type(Vector) :: fit_out,coeff,glob_coeff
        type(Matrix) :: A
        VecScatter :: scatter_ctx

        logical,allocatable :: luse_coeff(:)

        PetscScalar, pointer :: xv_inp(:),xv_out(:),xv_coeff(:)

        real(ireals),allocatable :: coords(:),coords_norm(:,:), dummy_inp(:), dummy_coeff(:)

        integer(iintegers) :: i,j,idz,ikabs,iksca,ig,icoeff

        call OPP%init(dx,dx,[azimuth],[zenith],imp_comm)

        Ndim = size(shape(OPP%diffLUT%S%c))-1
        Nsample = OPP%Ndz*OPP%Nkabs*OPP%Nksca*OPP%Ng
        Ncoeff = ubound(OPP%diffLUT%S%c,1)
        allocate(coords(Ndim))
        allocate(coords_norm(2,Ndim))

        coords_norm(:,1) = [ minval( OPP%diffLUT%pspace%dz   ), maxval( OPP%diffLUT%pspace%dz   ) ]
        coords_norm(:,2) = [ minval( OPP%diffLUT%pspace%kabs ), maxval( OPP%diffLUT%pspace%kabs ) ]
        coords_norm(:,3) = [ minval( OPP%diffLUT%pspace%ksca ), maxval( OPP%diffLUT%pspace%ksca ) ]
        coords_norm(:,4) = [ minval( OPP%diffLUT%pspace%g    ), maxval( OPP%diffLUT%pspace%g    ) ]

!        coords_norm(2,:) = one/ ( coords_norm(2,:)-coords_norm(1,:) )
!        coords_norm(2,:) = one/ ( coords_norm(2,:) )

        if(myid.eq.0) then
          print *,'LUT has',Ndim,'dimensions and ',Nsample,'entries and',Ncoeff,'coefficients'
          print *,'LUT dimensions :: dz  ',OPP%diffLUT%pspace%dz
          print *,'LUT dimensions :: kabs',OPP%diffLUT%pspace%kabs
          print *,'LUT dimensions :: ksca',OPP%diffLUT%pspace%ksca
          print *,'LUT dimensions :: g   ',OPP%diffLUT%pspace%g 
        endif

        call createVec(PETSC_DECIDE,Nsample,train_out)      
        call createVec(train_out%lN*Ndim ,Nsample*Ndim, train_inp)

        !Setup initial input and output for training dataset
        call VecGetArrayF90(  train_inp%v, xv_inp,   ierr)
        call VecGetArrayF90(  train_out%v, xv_out,   ierr)

        icoeff=2

        i=1
        j=-1 ! j counts the global entries... i counts only when we reach the local part
        do idz=1,OPP%Ndz
          do ikabs=1,OPP%Nkabs
            do iksca=1,OPP%Nksca
              do ig=1,OPP%Ng
                j=j+1
                if(j.lt.train_out%xs.or.j.ge.train_out%xe) cycle
                coords(1) = OPP%diffLUT%pspace%dz(idz)     / OPP%Ndz   !idz  !
                coords(2) = OPP%diffLUT%pspace%kabs(ikabs) / OPP%Nkabs !ikabs!
                coords(3) = OPP%diffLUT%pspace%ksca(iksca) / OPP%Nksca !iksca!
                coords(4) = OPP%diffLUT%pspace%g(ig)       / OPP%Ng    !ig   !

!                coords = (coords-coords_norm(1,:) ) * coords_norm(2,:)
!                coords = (coords) * coords_norm(2,:)

                xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim ) = coords

                xv_out(i) = OPP%diffLUT%S%c(icoeff,idz,ikabs,iksca,ig)
                i=i+1
              enddo
            enddo
          enddo
        enddo
        if(myid.eq.0) then
          print *,'we have set up till index',i,'::',(i-1)*Ndim + Ndim
          print *,'setting trainings out',xv_out(1:20)
        endif

        call VecRestoreArrayF90(  train_inp%v, xv_inp,   ierr)
        call VecRestoreArrayF90(  train_out%v, xv_out,   ierr)

        ! Find a fit for training set
        call createVec(PETSC_DECIDE, polygonal_number(Ndim,poly_order), coeff)
        allocate(luse_coeff(polygonal_number(Ndim,poly_order)) )

        call iterative_fit(poly_order,Ndim, train_inp, train_out, luse_coeff, coeff)

        ! Sample the result on same grid

        call createVec(PETSC_DECIDE, Nsample, fit_out)      

        call init_Mat(A,fit_out%lN,coeff%lN, fit_out%gN,coeff%gN)

        if(myid.eq.0) print *,'************************************ Solution **********************************'
        call fill_A(poly_order,Ndim, luse_coeff,train_inp, A)
        call MatMult(A%A,coeff%v,fit_out%v,ierr)

        call write_poly_fit(train_inp,train_out,fit_out,'./fitLUT.h5','train')

        call VecScatterCreateToZero(coeff%v, scatter_ctx, glob_coeff%v, ierr)
        call VecScatterBegin(scatter_ctx,coeff%v,glob_coeff%v,INSERT_VALUES,SCATTER_FORWARD,ierr);
        call VecScatterEnd(scatter_ctx,coeff%v,glob_coeff%v,INSERT_VALUES,SCATTER_FORWARD,ierr);
        call VecScatterDestroy(scatter_ctx,ierr);

        if(myid.eq.0) then
          call VecGetArrayF90(  glob_coeff%v, xv_coeff,   ierr)
          print *,'coeffs are now on comm 0:',xv_coeff

          allocate(dummy_inp(train_inp%gN))
          allocate(dummy_coeff(coeff%gN))
          dummy_coeff = poly(poly_order,Ndim, luse_coeff, dummy_inp, 'c1=', xv_coeff)
          call VecRestoreArrayF90(  glob_coeff%v, xv_coeff,   ierr)
        endif

        call VecDestroy(glob_coeff%v,ierr);

      end subroutine
!      subroutine create_benchmark(train_inp,train_out)
!      type(Vector) :: train_inp, train_out
!
!      integer(iintegers),parameter :: Ndim=2,Nsample=100      ! Size of original data (here created by test_func)
!      integer(iintegers),parameter :: Norder=4               ! Order of the polynomial that we want to fit
!
!      type(Vector) ::  coeff!, fit_out
!      type(Matrix) :: A
!      logical,allocatable :: luse_coeff(:)
!      
!      integer(iintegers) :: i
!      real(ireals) :: rnd
!
!      PetscScalar, pointer :: xv_inp(:),xv_out(:)
!
!      PetscViewer :: view
!
!      if(myid.eq.0) print *,'Maximum nr of polynomials for dimension/order',Ndim,Norder,' is',polygonal_number(Ndim,Norder)
!
!      call createVec(PETSC_DECIDE,Nsample,train_out)      
!      call createVec(train_out%lN*Ndim ,Nsample*Ndim, train_inp)
!
!      call createVec(PETSC_DECIDE, polygonal_number(Ndim,Norder), coeff)
!      allocate(luse_coeff(polygonal_number(Ndim,Norder)) )
!
!      call init_Mat(A,train_out%lN,coeff%lN, train_out%gN,coeff%gN)
!
!      call VecSet(coeff%v,one,ierr)
!
!      do i=1,size(luse_coeff)
!        call random_number(rnd)
!        luse_coeff(i) = .True.
!        if(rnd.lt..5_ireals) luse_coeff(i) = .False.
!      enddo
!
!      !Setup initial input and output for training dataset
!      call VecGetArrayF90(  train_inp%v, xv_inp,   ierr)
!      do i=1,train_out%lN
!        xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim ) = [ one*i+train_out%xs,  (one*i+train_out%xs)**3 ]/ [ train_out%gN, train_out%gN**3 ]
!      enddo
!      call VecRestoreArrayF90(  train_inp%v, xv_inp,   ierr)
!
!      if(myid.eq.0) print *,'************************************ Benchmark *********************************'
!      call fill_A(Norder,Ndim, luse_coeff,train_inp, A)
!      call MatMult(A%A,coeff%v,train_out%v,ierr)
!
!      call VecGetArrayF90(  train_out%v, xv_out,   ierr)
!      do i=1,train_out%lN
!        call random_number(rnd)
!        xv_out(i) = xv_out(i) * ( one + .3_ireals*(rnd-.5_ireals) )
!      enddo
!      call VecRestoreArrayF90(  train_out%v, xv_out,   ierr)
!
!      if(myid.eq.0) print *,'************************************ Benchmark *********************************'
!      end subroutine

      subroutine get_cmd_line_options()
        logical :: lflg=.False.
        call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-poly",poly_func, lflg,ierr)  ; CHKERRQ(ierr)
        if(lflg.eqv.PETSC_FALSE) poly_func = 2

        call PetscOptionsGetInt(PETSC_NULL_CHARACTER,"-order",poly_order, lflg,ierr)  ; CHKERRQ(ierr)
        if(lflg.eqv.PETSC_FALSE) poly_order = 4

        call PetscOptionsGetReal(PETSC_NULL_CHARACTER,"-err_ratio",err_ratio, lflg,ierr)  ; CHKERRQ(ierr)
        if(lflg.eqv.PETSC_FALSE) err_ratio = .99

        print *,' ************* poly_func is now ',poly_func,' *************'
        print *,' ************* poly_order is now',poly_order,' *************'
        print *,' ************* error ratio is : ',err_ratio,' *************'
      end subroutine
end module

program main
      use m_fitLUT 
      implicit none

!      integer(iintegers),parameter :: Ndim=2,Nsample=100          ! Size of original data (here created by test_func)
!      integer(iintegers),parameter :: Norder=2                    ! Order of the polynomial that we want to fit

!      type(Vector) :: train_inp, train_out, coeff, fit_out
!
!      type(Matrix) :: A
!      logical,allocatable :: luse_coeff(:)
!      
!      integer(iintegers) :: i
!      real(ireals) :: rnd
!
!      character(len=300) :: groups(3)
!
!      PetscRandom :: rctx
!      PetscScalar, pointer :: xv_inp(:),xv_out(:)

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr) ;CHKERRQ(ierr)
      call init_mpi_data_parameters(PETSC_COMM_WORLD)
      call get_cmd_line_options()

      call fitLUT(40._ireals,zero,zero)

      call PetscFinalize(ierr) ;CHKERRQ(ierr);  return

!
!
!
!
!      call create_benchmark(train_inp,train_out)
!
!      ! Find a fit for benchmark training set
!      call createVec(PETSC_DECIDE, polygonal_number(Ndim,Norder), coeff)
!      allocate(luse_coeff(polygonal_number(Ndim,Norder)) )
!
!      call iterative_fit(Norder,Ndim, train_inp, train_out, luse_coeff, coeff)
!     
!      ! Sample the result on a fine grid
!      call createVec(PETSC_DECIDE, Nsample, fit_out)      
!
!      call init_Mat(A,fit_out%lN,coeff%lN, fit_out%gN,coeff%gN)
!
!      !Setup input grid 
!      call VecGetArrayF90(  train_inp%v, xv_inp,   ierr)
!      do i=1,fit_out%lN
!        xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim ) = [ one*i+fit_out%xs,  (one*i+fit_out%xs)**3 ]/ [ fit_out%gN, fit_out%gN**3 ]
!      enddo
!      call VecRestoreArrayF90(  train_inp%v, xv_inp,   ierr)
!
!      if(myid.eq.0) print *,'************************************ Solution **********************************'
!      call fill_A(Norder,Ndim, luse_coeff,train_inp, A)
!      call MatMult(A%A,coeff%v,fit_out%v,ierr)
!!      call VecView(train_out%v, PETSC_VIEWER_STDOUT_WORLD ,ierr)
!      if(myid.eq.0) print *,'************************************ Solution **********************************'
!      call VecAXPY(fit_out%v,-one,train_out%v,ierr) ! a-b
!      call VecNorm(fit_out%v, NORM_2 , rnd,ierr)
!      print *,'rmse is :',rnd/Nsample
!
!!      call PetscViewerHDF5Open(imp_comm,'./fitLUT.h5', FILE_MODE_WRITE, view, ierr) ;CHKERRQ(ierr)
!!      call PetscObjectSetName(train_inp%v, 'train_inp',ierr)          ; CHKERRQ(ierr);   call VecView(train_inp%v, view, ierr) ;CHKERRQ(ierr)
!!      call PetscObjectSetName(train_out%v, 'train_out',ierr)          ; CHKERRQ(ierr);   call VecView(train_out%v, view, ierr) ;CHKERRQ(ierr)
!!      call PetscObjectSetName(fit_out%v, 'fit_out',ierr)          ; CHKERRQ(ierr);   call VecView(fit_out%v, view, ierr) ;CHKERRQ(ierr)
!!      call PetscViewerDestroy(view,ierr) ;CHKERRQ(ierr)
!
!
!      call PetscFinalize(ierr) ;CHKERRQ(ierr)
end program
