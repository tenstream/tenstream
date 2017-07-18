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

module m_poly_fitLUT

#ifdef _XLF
      use ieee_arithmetic
#define isnan ieee_is_nan
#endif

#include "petsc/finclude/petsc.h"
      use petsc

      use m_helper_functions, only: rmse, CHKERR
      use m_data_parameters, only: ireals,iintegers, &
        init_mpi_data_parameters, numnodes, myid,    &
        zero, one, imp_comm, mpiint, nil,            &
        default_str_len
      use m_netcdfio

      use m_optprop_LUT, only : t_optprop_LUT_8_10

      implicit none

      integer(mpiint) :: ierr
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
      logical :: lfull_poly

    contains
      subroutine iterative_fit(Norder,Ndim, train_inp, train_out, luse_coeff, coeff)
        integer(iintegers),intent(in) :: Norder,Ndim
        type(Vector) :: train_inp, train_out, coeff
        logical,intent(out) :: luse_coeff(:)

        logical :: ltmp_coeff(size(luse_coeff))

        type(Matrix) :: A
        type(Vector) :: tmp_coeff,tmp_out

        KSP :: ksp

        real(ireals) :: err=huge(err),tmp_err=huge(err),last_err,avg,rmse
        logical :: learlyexit

        integer(iintegers) :: k,j, jmin

        if(myid.eq.0) print *,'Maximum nr of polynomials for dimension/order',Ndim,Norder,' is',polygonal_number(Ndim,Norder)

        call VecSum(train_out%v, avg, ierr); call CHKERR(ierr)
        avg = avg/train_out%gN

        call VecSet(coeff%v         , zero, ierr)

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
          learlyexit=.False.
          do j=1,size(luse_coeff)
            ltmp_coeff = luse_coeff
            if( ltmp_coeff(j) ) cycle ! already included
            ltmp_coeff(j) = .True.

            if(lfull_poly) then
              ltmp_coeff = .True.
              luse_coeff = .True.
              learlyexit = .True.
            endif

            call MatZeroEntries(A%A,ierr)
            call fill_A(Norder,Ndim, ltmp_coeff,train_inp, A)

            call solve_lsqr(train_out, tmp_coeff , A, ksp)
            call KSPGetResidualNorm(ksp,tmp_err,ierr)


            if( tmp_err/err .lt. .99_ireals ) then !enhancement big --> exit directly and do not search if we could find another, better one....
              if(myid.eq.0) print *,'adding coeff',j,' because error was reduced by big factor:', tmp_err/err
              luse_coeff(j) = .True.
            endif

            if( tmp_err.lt.err ) then 
              err=tmp_err
              jmin = j
              call VecCopy(tmp_coeff%v,coeff%v,ierr)

              ! calculate root mean squared error
              call MatMult(A%A,coeff%v,tmp_out%v, ierr)
              call VecAXPY(tmp_out%v, -one, train_out%v, ierr) 
              call VecScale(tmp_out%v, sqrt( one/tmp_out%gN ), ierr) 
              call VecNorm(tmp_out%v, NORM_2, rmse, ierr)

              if( rmse/avg .lt. 20e-2_ireals ) learlyexit=.True. ! relative error smaller
              if(myid.eq.0) print *,'k=',k,'N',count(luse_coeff),' :: testing new coeff:',j,'(',size(luse_coeff),')','old err',err,err/train_inp%gN,'with new err',tmp_err,tmp_err/train_inp%gN ,'::',rmse/avg*100,'%'
            endif

            if(learlyexit)  exit 
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
          if(myid.eq.0) print *,'iteration',k,' :: after searching through all coeffs, adding new coeff:',jmin,'with err',err,err/train_inp%gN  !,'::',luse_coeff

          ! calculate root mean squared error
          call MatMult(A%A,coeff%v,tmp_out%v, ierr)
          call VecAXPY(tmp_out%v, -one, train_out%v, ierr) 
          call VecScale(tmp_out%v, sqrt( one/tmp_out%gN ), ierr) 
          call VecNorm(tmp_out%v, NORM_2, rmse, ierr)


          if(myid.eq.0) print *,'relative error is:', rmse/avg*100._ireals,'%'
          if( rmse/avg .lt. 50e-2_ireals ) then ! relative error smaller 50%
            if(myid.eq.0) print *,'exiting iterative fit routine :: we converged with relative RMSE:',rmse/avg*100._ireals,'%'
            exit
          endif
        enddo

        call MatDestroy(A%A,ierr)
        call KSPDestroy(ksp,ierr)
        call VecDestroy(tmp_coeff%v,ierr)

        if(myid.eq.0) print *,'Final coeffs:',luse_coeff
!        call VecView( coeff%v,PETSC_VIEWER_STDOUT_WORLD,ierr)

    end subroutine

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

            if(any(isnan(xvals)).or.any(isnan(xvals_inp)) ) then
              print *,irow,':',irow-A%xs,'xvals_inp',xvals_inp,'xvals',xvals,'A%gNcol',A%gNcol,'col_ind',col_ind
              stop 'fill_A :: NaN in xvals -- bug!'
            endif
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
          call MatMPIAIJSetPreallocation(A%A, lNcol,PETSC_NULL_INTEGER, gNcol-lNcol,PETSC_NULL_INTEGER,ierr);call CHKERR(ierr)
          call MatSeqAIJSetPreallocation(A%A, lNcol, PETSC_NULL_INTEGER, ierr) ;call CHKERR(ierr)
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

          integer(iintegers) :: j,xd,yd,zd,kd,ld,md

          if(present(print_func) ) write( *,FMT='(A)', advance='no') trim(print_func)
          j=1
          select case (Ndim)
          case (6)
            do xd=0,Norder
              do yd=0,Norder-xd
                if(xd+yd.gt.Norder) cycle
                do zd=0,Norder-yd
                  if(xd+yd+zd.gt.Norder) cycle
                  do kd=0,Norder-zd
                    if(xd+yd+zd+kd.gt.Norder) cycle
                    do ld=0,Norder-kd
                      if(xd+yd+zd+kd+ld.gt.Norder) cycle
                      do md=0,Norder-ld
                        if(xd+yd+zd+kd+ld+md.gt.Norder) cycle

!                    print *,'poly:',xd,yd,zd,kd,'::',j
                        if(luse_coeff(j) ) then

                          if(present(print_func) ) then
                            write( *,FMT='( "+",E14.7,"* ( x**",I0,"* y**",I0,"* z**",I0,"* k**",I0,  "* l**",I0,  "* m**",I0, " ) ")',advance='no') coeff(j),xd,yd,zd,kd,ld,md
                          else
                            poly(j) = p(inp(1),xd) * p(inp(2),yd) * p(inp(3),zd) * p(inp(4),kd) * p(inp(5),ld) * p(inp(6),md)
                          endif

                        else
                          poly(j) = nil
                        endif
                        j=j+1
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
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
                      poly(j) = nil
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
                  poly(j) = nil
                endif
                j=j+1
              enddo
            enddo

          case default
            print *,'did not implement poly function for ',Ndim,'dimensions'
            call CHKERR(1_mpiint)
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
!            p = legendre(2._ireals*v-1,e)
            p = legendre(v,e)
          case(3)
            p = laguerre(v,e)
          case(4)
            p = chebyshev(v,e)
          case default
!          stop 'dont know poly_func option'
            p = legendre(2._ireals*v-1,e)
          end select
        end function

     pure function legendre_func ( x, n )
      integer (iintegers),intent(in) :: n
      real (ireals),intent(in) :: x
      real (ireals) :: legendre_func
      real(ireals) :: v(0:n)

      integer (iintegers) :: j

      v(0) = one
      legendre_func = one
      if(n.eq.0) return

      v(1) = x

      do j = 2, n
        v(j) = ( 2._ireals*(j-1)*x * v(j-1) - (j-1)*v(j-2) )/j
      end do
      legendre_func = v(n)
    end function
     pure function chebyshev ( x, n )
      integer (iintegers),intent(in) :: n
      real (ireals),intent(in) :: x
      real (ireals) :: chebyshev
      real(ireals) :: v(0:n)

      integer (iintegers) :: j

      v(0) = one
      chebyshev = one
      if(n.eq.0) return

      v(1) = x

      do j = 2, n
        v(j) = 2._ireals * v(j-1) - v(j-2)
      end do
      chebyshev = v(n)
    end function

        pure function legendre(x,n)
          real(ireals) :: legendre
          real(ireals),intent(in) :: x
          integer(iintegers),intent(in) :: n
          real(ireals),dimension(0:n) :: pn,pd
          !http://de.wikipedia.org/wiki/Legendre-Polynom
!        legendre=0
!        do k=0,n/2
!          legendre = legendre + (one)**k * factorial(2*n-2*k) /factorial(n-k)/factorial(n-2*k)/factorial(k)/2**n * x**(n-2*k)
!        enddo
          call LPN(n,2._ireals*x-one, pn,pd)
!        call LPN(n,x, pn,pd)
          legendre=pn(n)
        end function

      pure function laguerre(x,n)
        real(ireals) :: laguerre
        real(ireals),intent(in) :: x
        integer(iintegers),intent(in) :: n
        real (ireals) v(0:n)
        call laguerre_polynomial ( n, exp(-x), v )
        laguerre=v(n)
      end function
     pure subroutine laguerre_polynomial ( n, x, v )
      implicit none

      integer (iintegers),intent(in) :: n
      real (ireals),intent(in) :: x
      real (ireals),intent(out) :: v(0:n)

      integer (iintegers) :: j

      v(0) = one
      if(n.eq.0) return

      v(1) = one - x

      do j = 2, n
        v(j) = ( ( 2._ireals*(j-1) - x )*v(j-1) - (j-1)*v(j-2) ) / j
      end do
    end subroutine
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
          PD(0)=zero
          if(N.eq.0) return
          PN(1)=X
          PD(1)=one
          if(N.eq.1) return
          P0=one
          P1=X
          DO K=2,N
            PF=(2._ireals*K-one)/K*X*P1-(K-one)/K*P0
            PN(K)=PF
            IF (ABS(X).EQ.one) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
            ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
            ENDIF
            P0=P1
            P1=PF
          enddo
          RETURN
          END SUBROUTINE

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
          subroutine write_poly_fit(inp,out,fname,prefix)
            type(Vector),intent(in) :: inp,out
            character(len=*) :: fname,prefix
            character(default_str_len) :: vecname

            PetscFileMode :: fmode
            logical :: fexists

            if(myid.eq.0) print *,'**************************** Writing Solution **********************************'
            write( vecname,FMT='(A,".poly",I0,".order",I0,".full",L1,".")' ) trim(prefix), poly_func, poly_order,lfull_poly

            inquire(file=trim(fname), exist=fexists)
            if(fexists) then
              if(myid.eq.0)  print *,myid,'appending vector to hdf5 file ',trim(fname),' vecname: ',vecname
              fmode = FILE_MODE_APPEND
            else 
              if(myid.eq.0)  print *,myid,'writing vector to hdf5 file ',trim(fname),' vecname: ',vecname
              fmode = FILE_MODE_WRITE
            endif

!            call PetscViewerHDF5Open(imp_comm,fname, fmode, view, ierr) ;call CHKERR(ierr)
            call PetscObjectSetName(inp%v, trim(vecname)//'inp',ierr) ; call CHKERR(ierr);   call VecView(inp%v, view, ierr) ;call CHKERR(ierr)
            call PetscObjectSetName(out%v, trim(vecname)//'out',ierr) ; call CHKERR(ierr);   call VecView(out%v, view, ierr) ;call CHKERR(ierr)
            call PetscViewerDestroy(view,ierr) ;call CHKERR(ierr)
          end subroutine

          subroutine fitdirLUT(azimuth,zenith)
            real(ireals),intent(in) :: azimuth,zenith

            type(t_optprop_LUT_8_10) :: OPP
            integer(iintegers) :: Ndim,Nsample,Ncoeff

            type(Vector) :: train_inp, train_out
            type(Vector) :: fit_out,coeff,glob_coeff
            type(Matrix) :: A
            VecScatter :: scatter_ctx

            logical,allocatable :: luse_coeff(:)

            PetscScalar, pointer :: xv_inp(:),xv_out(:),xv_coeff(:)

            real(ireals),allocatable :: coords(:), dummy_inp(:), dummy_coeff(:)

            integer(iintegers) :: i,j,itaux,itauz,iw0,ig,iphi,itheta,icoeff

            call OPP%init([azimuth],[zenith],imp_comm)

            icoeff=2

            Ndim = 6
            Nsample = 0
            do iphi=1,OPP%Nphi
              do itheta=1,OPP%Ntheta
                if(allocated(OPP%dirLUT%S(iphi,itheta)%c) ) then
                  Ncoeff = ubound(OPP%dirLUT%S(iphi,itheta)%c,1)
                  Nsample = Nsample + size(OPP%dirLUT%S(iphi,itheta)%c(icoeff,:,:,:,:) )
                endif
              enddo
            enddo
            allocate(coords(Ndim))

            if(myid.eq.0) then
              print *,'LUT has',Ndim,'dimensions and ',Nsample,'entries and',Ncoeff,'coefficients'
              print *,'LUT dimensions :: tau  ',OPP%dirLUT%pspace%tau
              print *,'LUT dimensions :: w0   ',OPP%dirLUT%pspace%w0
              print *,'LUT dimensions :: g    ',OPP%dirLUT%pspace%g
              print *,'LUT dimensions :: phi  ',OPP%dirLUT%pspace%phi
              print *,'LUT dimensions :: theta',OPP%dirLUT%pspace%theta
            endif

            call createVec(PETSC_DECIDE,Nsample,train_out)
            call createVec(train_out%lN*Ndim ,Nsample*Ndim, train_inp)

            !Setup initial input and output for training dataset
            call VecGetArrayF90(  train_inp%v, xv_inp,   ierr)
            call VecGetArrayF90(  train_out%v, xv_out,   ierr)

            i=1
            j=-1 ! j counts the global entries... i counts only when we reach the local part
            do iphi=1,OPP%Nphi
              do itheta=1,OPP%Ntheta
                if( .not. allocated(OPP%dirLUT%S(iphi,itheta)%c) ) cycle
                do iw0=1,OPP%Nw0
                  do itauz=1,OPP%Ntau
                    do itaux=1,OPP%Ntau
                      do ig=1,OPP%Ng
                        j=j+1
                        if(j.lt.train_out%xs.or.j.ge.train_out%xe) cycle
                        coords(1) = (itaux -one)/ (OPP%Ntau  -one)
                        coords(2) = (itauz -one)/ (OPP%Ntau  -one)
                        coords(3) = (iw0   -one)/ (OPP%Nw0   -one)
                        coords(4) = (ig    -one)/ (OPP%Ng    -one)
                        coords(4) = (iphi  -one)/ (OPP%Nphi  -one)
                        coords(4) = (itheta-one)/ (OPP%Ntheta-one)

                        xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim ) = coords

                        xv_out(i) = OPP%dirLUT%S(iphi,itheta)%c(icoeff,itaux,itauz,iw0,ig)
                        i=i+1
                      enddo
                    enddo
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

            call write_poly_fit(train_inp,train_out,'./fitLUT.h5','train.dir.S')

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
          subroutine run_fit4d(x,y,z,k, luse_coeff,coeff,prefix)
            real(ireals),intent(in),dimension(:) :: x,y,z,k
            logical,intent(in) :: luse_coeff(:)
            type(Vector),intent(in) :: coeff
            character(len=*),intent(in) :: prefix

            integer(iintegers),parameter :: Ndim=4
            integer(iintegers)           :: Nsample

            type(Vector) :: inp,out
            type(Matrix) :: A

!            type(Vector) :: glob_coeff
!            VecScatter :: scatter_ctx
!            PetscScalar, pointer :: xv_coeff(:)
!            real(ireals),allocatable :: dummy_inp(:), dummy_coeff(:)

            Nsample = size(x)*size(y)*size(z)*size(k)

            call createVec(PETSC_DECIDE, Nsample      , out)
            call createVec(out%lN*Ndim , Nsample*Ndim , inp)
            call create_grid4d(x,y,z,k,inp)

            ! Sample the result on grid
            call init_Mat(A,out%lN,coeff%lN, out%gN,coeff%gN)

            call fill_A(poly_order,Ndim, luse_coeff,inp, A)
            call MatMult(A%A,coeff%v,out%v,ierr)

            call MatDestroy(A%A,ierr)

            call write_poly_fit(inp,out,'./fitLUT.h5',prefix)

!            call VecScatterCreateToZero(coeff%v, scatter_ctx, glob_coeff%v, ierr)
!            call VecScatterBegin(scatter_ctx,coeff%v,glob_coeff%v,INSERT_VALUES,SCATTER_FORWARD,ierr);
!            call VecScatterEnd(scatter_ctx,coeff%v,glob_coeff%v,INSERT_VALUES,SCATTER_FORWARD,ierr);
!            call VecScatterDestroy(scatter_ctx,ierr);
!
!            if(myid.eq.0) then
!              call VecGetArrayF90(  glob_coeff%v, xv_coeff,   ierr)
!              print *,'coeffs are now on comm 0:',xv_coeff
!
!              allocate(dummy_inp(train_inp%gN))
!              allocate(dummy_coeff(coeff%gN))
!              dummy_coeff = poly(poly_order,Ndim, luse_coeff, dummy_inp, 'c1=', xv_coeff)
!              call VecRestoreArrayF90(  glob_coeff%v, xv_coeff,   ierr)
!            endif
!            call VecDestroy(glob_coeff%v,ierr);

            call VecDestroy(inp%v,ierr)
            call VecDestroy(out%v,ierr)
          end subroutine
          subroutine create_grid4d(x,y,z,k,inp, reverse)
            real(ireals),intent(in),dimension(:) :: x,y,z,k
            type(Vector),intent(inout) :: inp
            logical,optional  :: reverse

            integer(iintegers),parameter :: Ndim=4
            integer(iintegers) :: i,j,ix,iy,iz,ik
            real(ireals) :: coords(Ndim),norm(2,Ndim) ! Ndim entries for shifted coordinates and norm are 2 entries for value normalization
            PetscScalar, pointer :: xv_inp(:)

            norm(:,1) = [ minval(x), one/(maxval(x)-minval(x)) ]
            norm(:,2) = [ minval(y), one/(maxval(y)-minval(y)) ]
            norm(:,3) = [ minval(z), one/(maxval(z)-minval(z)) ]
            norm(:,4) = [ minval(k), one/(maxval(k)-minval(k)) ]
!            if( present(reverse) ) then
!              if(reverse) then !TODO this section does not make any sense at all.... need to get a coffee.... initial idea was to somehow retransform from [0,1] to original input space
!                norm(:,1) = [ minval(x), one/(maxval(x)-minval(x)) ]
!                norm(:,2) = [ minval(y), one/(maxval(y)-minval(y)) ]
!                norm(:,3) = [ minval(z), one/(maxval(z)-minval(z)) ]
!                norm(:,4) = [ minval(k), one/(maxval(k)-minval(k)) ]
!              endif 
!            endif
!
!            if( poly_func .eq. 3) then !Laguerre polynomials dont need normalization?
!                norm(1,:) = zero
!                norm(2,:) = one
!            endif

            if(myid.eq.0) then
              print *,'grid has',Ndim,'dimensions'
              print *,'grid dimensions :: x',x
              print *,'grid dimensions :: y',y
              print *,'grid dimensions :: z',z
              print *,'grid dimensions :: k',k
            endif
            !Setup initial input and output for training dataset
            call VecGetArrayF90(  inp%v, xv_inp,   ierr)

            i=1
            j=0 ! j counts the global entries... i counts only when we reach the local part
            do ik=1,size(k)
              do iz=1,size(z)
                do iy=1,size(y)
                  do ix=1,size(x)
                    j=j+1
                    if( (j-1)*Ndim .lt. inp%xs .or. (j-1)*Ndim+Ndim .ge. inp%xe) cycle
                    coords(1) = (x(ix)-norm(1,1)  )*norm(2,1) ! linear transformation into orthogonal space of legendre polynomials [0,1]!(ix-one) / (size(x)-one)  !
                    coords(2) = (y(iy)-norm(1,2)  )*norm(2,2)                                                                            !(iy-one) / (size(y)-one)  !
                    coords(3) = (z(iz)-norm(1,3)  )*norm(2,3)                                                                            !(iz-one) / (size(z)-one)  !
                    coords(4) = (k(ik)-norm(1,4)  )*norm(2,4)                                                                            !(ik-one) / (size(k)-one)  !

                    xv_inp( (i-1)*Ndim + 1 : (i-1)*Ndim + Ndim ) = coords
                    i=i+1
                  enddo
                enddo
              enddo
            enddo
            call VecRestoreArrayF90(  inp%v, xv_inp,   ierr)
          end subroutine
          subroutine fit4d(x,y,z,k, val, luse_coeff, coeff)
            real(ireals),intent(in),dimension(:) :: x,y,z,k
            real(ireals),intent(in) :: val(:,:,:,:)
            logical,allocatable,intent(out) :: luse_coeff(:)

            integer(iintegers),parameter :: Ndim=4
            integer(iintegers)           :: Nsample

            type(Vector) :: train_inp, train_out, coeff


            PetscScalar, pointer :: xv_out(:)

            integer(iintegers) :: i,j,ix,iy,iz,ik

            Nsample = size(val)

            call createVec(PETSC_DECIDE      , Nsample      , train_out)
            call createVec(train_out%lN*Ndim , Nsample*Ndim , train_inp)
            call create_grid4d(x,y,z,k,train_inp)

            !Setup initial input and output for training dataset
            call VecGetArrayF90(  train_out%v, xv_out,   ierr)

            i=1
            j=0 ! j counts the global entries... i counts only when we reach the local part
            do ik=1,size(k)
              do iz=1,size(z)
                do iy=1,size(y)
                  do ix=1,size(x)
                    j=j+1
                    if(j.lt.train_out%xs+1.or.j.gt.train_out%xe) cycle
                    xv_out(i) = val(ix ,iy ,iz ,ik)
!                    xv_out(i) = ( xv_out(i) - norm(1) ) * norm(2)  ! now normalized to [ 0, 1]
!                    xv_out(i) = 2._ireals*xv_out(i) -one           ! now normalized to [-1, 1]

                    i=i+1
                  enddo
                enddo
              enddo
            enddo

            call VecRestoreArrayF90(  train_out%v, xv_out,   ierr)

            ! Find a fit for training set
            call createVec(PETSC_DECIDE, polygonal_number(Ndim,poly_order), coeff)
            allocate(luse_coeff(polygonal_number(Ndim,poly_order)) )

            call iterative_fit(poly_order,Ndim, train_inp, train_out, luse_coeff, coeff)

            call write_poly_fit(train_inp,train_out,'./fitLUT.h5','train')

            call VecDestroy(train_inp%v,ierr)
            call VecDestroy(train_out%v,ierr)
          end subroutine
          subroutine fitLUT(azimuth,zenith)
            real(ireals),intent(in) :: azimuth,zenith

            type(t_optprop_LUT_8_10) :: OPP

            type(Vector) :: coeff
            logical,allocatable :: luse_coeff(:)

            integer(iintegers) :: icoeff,i
            real(ireals),allocatable,dimension(:) :: x,y,z,k
            real(ireals),allocatable,dimension(:) :: xx,yy,zz,kk

            call OPP%init([azimuth],[zenith],imp_comm)

            icoeff = 2

            allocate( x(size(OPP%diffLUT%pspace%tau )),source=(OPP%diffLUT%pspace%tau )**.1_ireals )
            allocate( y(size(OPP%diffLUT%pspace%tau )),source=(OPP%diffLUT%pspace%tau )**.1_ireals )
            allocate( z(size(OPP%diffLUT%pspace%w0  )),source=(OPP%diffLUT%pspace%w0  )**.1_ireals )
            allocate( k(size(OPP%diffLUT%pspace%g   )),source=OPP%diffLUT%pspace%g )

            allocate( xx(2*size(x)-1 ) ) ; xx(1:size(xx):2) = x;  xx(2:size(xx):2) = [ ( (xx(i-1)+xx(i+1))/2  , i=2,size(xx),2 ) ]
            allocate( yy(2*size(y)-1 ) ) ; yy(1:size(yy):2) = y;  yy(2:size(yy):2) = [ ( (yy(i-1)+yy(i+1))/2  , i=2,size(yy),2 ) ]
            allocate( zz(2*size(z)-1 ) ) ; zz(1:size(zz):2) = z;  zz(2:size(zz):2) = [ ( (zz(i-1)+zz(i+1))/2  , i=2,size(zz),2 ) ]
            allocate( kk(2*size(k)-1 ) ) ; kk(1:size(kk):2) = k;  kk(2:size(kk):2) = [ ( (kk(i-1)+kk(i+1))/2  , i=2,size(kk),2 ) ]


            call fit4d    (x,y,z,k, OPP%diffLUT%S%c(icoeff,:,:,:,:), luse_coeff, coeff)
            call run_fit4d(x,y,z,k, luse_coeff,coeff,'fit')

            call run_fit4d(xx,yy,zz,kk, luse_coeff,coeff,'dfit')

          end subroutine

          subroutine get_cmd_line_options()
            logical :: lflg=.False.
            call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-poly",poly_func, lflg,ierr)  ; call CHKERR(ierr)
            if(lflg.eqv.PETSC_FALSE) poly_func = 2

            call PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-order",poly_order, lflg,ierr)  ; call CHKERR(ierr)
            if(lflg.eqv.PETSC_FALSE) poly_order = 4

            call PetscOptionsGetReal(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-err_ratio",err_ratio, lflg,ierr)  ; call CHKERR(ierr)
            if(lflg.eqv.PETSC_FALSE) err_ratio = .9999

            call PetscOptionsGetBool(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-full_poly",lfull_poly,lflg,ierr) ;call CHKERR(ierr)
            if(lflg.eqv.PETSC_FALSE) lfull_poly = .False.

            print *,' ************* poly_func is now ',poly_func,' *************'
            print *,' ************* poly_order is now',poly_order,' *************'
            print *,' ************* error ratio is : ',err_ratio,' *************'
          end subroutine
end module

program main
          use m_poly_fitLUT
          implicit none

          call PetscInitialize(PETSC_NULL_CHARACTER,ierr) ;call CHKERR(ierr)
          call init_mpi_data_parameters(PETSC_COMM_WORLD)
          call get_cmd_line_options()

          call fitLUT(zero,zero)

          call PetscFinalize(ierr) ;call CHKERR(ierr) 

end program
