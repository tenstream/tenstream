module gridtransform
      use iso_c_binding
      implicit none

      type gridinfo
              double precision,allocatable,dimension(:) :: dx,dy,dz
              double precision,allocatable,dimension(:) :: x,y,z
              integer(c_int) :: Nx,Ny,Nz
      end type

      type(gridinfo) :: oldgrid,newgrid

      double precision :: solution_dx=-1,solution_dy=-1,solution_dz=-1
!      double precision,parameter :: solution_dx=250,solution_dy=250,solution_dz=200

      contains

      subroutine init_grid_transform(odx,ody,odz,lsame)
        double precision,dimension(:),intent(in) :: odx,ody,odz
        logical,intent(in) :: lsame

        integer(c_int) i,j,k

        if(allocated(newgrid%z)) return

        oldgrid%Nx=size(odx)
        oldgrid%Ny=size(ody)
        oldgrid%Nz=size(odz)

        allocate(oldgrid%dx(oldgrid%Nx)) ; oldgrid%dx(:) = odx
        allocate(oldgrid%dy(oldgrid%Ny)) ; oldgrid%dy(:) = ody
        allocate(oldgrid%dz(oldgrid%Nz)) ; oldgrid%dz(:) = odz

        allocate(oldgrid%x(lbound(oldgrid%dx,1):ubound(oldgrid%dx,1)+1)) ; oldgrid%x(lbound(oldgrid%x,1)) = 0
        allocate(oldgrid%y(lbound(oldgrid%dy,1):ubound(oldgrid%dy,1)+1)) ; oldgrid%y(lbound(oldgrid%y,1)) = 0
        allocate(oldgrid%z(lbound(oldgrid%dz,1):ubound(oldgrid%dz,1)+1)) ; oldgrid%z(ubound(oldgrid%z,1)) = 0

        do i=1,ubound(oldgrid%x,1)-1
           oldgrid%x(i+1) = oldgrid%x(i) + oldgrid%dx(i)
        enddo
        do j=1,ubound(oldgrid%y,1)-1
           oldgrid%y(j+1) = oldgrid%y(j) + oldgrid%dy(j)
        enddo
        do k=ubound(oldgrid%z,1),lbound(oldgrid%z,1)+1,-1
           oldgrid%z(k-1) = oldgrid%z(k) + oldgrid%dz(k-1)
        enddo

        if(lsame) then
                newgrid%Nx = oldgrid%Nx
                newgrid%Ny = oldgrid%Ny
                newgrid%Nz = oldgrid%Nz
        else
                newgrid%Nx = int(sum(odx)/solution_dx)
                newgrid%Ny = int(sum(ody)/solution_dy)
                newgrid%Nz = int(sum(odz)/solution_dz)
        endif


        allocate(newgrid%dx(newgrid%Nx)) 
        allocate(newgrid%dy(newgrid%Ny)) 
        allocate(newgrid%dz(newgrid%Nz)) 
        if(lsame) then
                newgrid%dx=oldgrid%dx
                newgrid%dy=oldgrid%dy
                newgrid%dz=oldgrid%dz
        else
                newgrid%dx=solution_dx
                newgrid%dy=solution_dy
                newgrid%dz=solution_dz
        endif

        allocate(newgrid%x(newgrid%Nx+1)) ;newgrid%x(1)=0
        allocate(newgrid%y(newgrid%Ny+1)) ;newgrid%y(1)=0
        allocate(newgrid%z(newgrid%Nz+1)) ;newgrid%z(newgrid%Nz+1)=0

        do i=1,ubound(newgrid%x,1)-1
           newgrid%x(i+1) = newgrid%x(i) + newgrid%dx(i)
        enddo
        do j=1,ubound(newgrid%y,1)-1
           newgrid%y(j+1) = newgrid%y(j) + newgrid%dy(j)
        enddo
        do k=ubound(newgrid%z,1),lbound(newgrid%z,1)+1,-1
           newgrid%z(k-1) = newgrid%z(k) + newgrid%dz(k-1)
        enddo

        print *,'Old grid had',oldgrid%Nx,oldgrid%Ny,oldgrid%Nz,'boxes'
        print *,'New grid has',newgrid%Nx,newgrid%Ny,newgrid%Nz,'boxes'
!
!        print *,'-------------------------------'
!        print *,'dx',oldgrid%dx,' :: ',newgrid%dx
!        print *,'dy',oldgrid%dy,' :: ',newgrid%dy
!        print *,'dz',oldgrid%dz,' :: ',newgrid%dz
!        print *,'-------------------------------'
!        print *,'x',oldgrid%x,' :: ',newgrid%x
!        print *,'y',oldgrid%y,' :: ',newgrid%y
!        print *,'z',oldgrid%z,' :: ',newgrid%z
      end subroutine

      subroutine grid_old_to_new(A)
        double precision,allocatable,dimension(:,:,:),intent(inout) :: A
        double precision,allocatable,dimension(:,:,:) :: tmp

        integer(c_int) i,j,k,oi,oj,ok
        double precision :: val

        allocate(tmp(size(newgrid%dx),size(newgrid%dy),size(newgrid%dz) ))
        tmp = A
        deallocate(A)
        allocate(A(size(newgrid%dx),size(newgrid%dy),size(newgrid%dz)) )

        do i=1,ubound(A,1)
                val = newgrid%x(i)+newgrid%dx(i)/2
                oi = search_sorted_bisec(oldgrid%x,val) -1

                do j=1,ubound(A,2)
                        val = newgrid%y(j)+newgrid%dy(j)/2
                        oj = search_sorted_bisec(oldgrid%y,val) -1

                        do k=1,ubound(A,3)
                                val = newgrid%z(k)-newgrid%dz(k)/2
                                ok = search_sorted_bisec(oldgrid%z,val)

                                A(i,j,k) = tmp(oi,oj,ok)
                        enddo
                enddo
        enddo
      end subroutine
      subroutine grid_new_to_old(A)
        double precision,allocatable,dimension(:,:,:),intent(inout) :: A
        double precision,dimension(lbound(A,1):ubound(A,1),lbound(A,2):ubound(A,2),lbound(A,3):ubound(A,3)) :: tmp

        integer(c_int) i,j,k,oi,oj,ok
        double precision :: val

        tmp = A
        deallocate(A)
        allocate(A(size(oldgrid%dx),size(oldgrid%dy),size(oldgrid%dz)) )

        do i=lbound(oldgrid%dx,1),ubound(oldgrid%dx,1)
                val = oldgrid%x(i)+oldgrid%dx(i)/2
                oi = search_sorted_bisec(newgrid%x,val) -1

                do j=lbound(oldgrid%dy,1),ubound(oldgrid%dy,1)
                        val = oldgrid%y(j)+oldgrid%dy(j)/2
                        oj = search_sorted_bisec(newgrid%y,val) -1

                        do k=lbound(oldgrid%dz,1),ubound(oldgrid%dz,1)
                                val = oldgrid%z(k)-oldgrid%dz(k)/2
                                ok = search_sorted_bisec(newgrid%z,val)

                                A(i,j,k) = tmp(oi,oj,ok)
                        enddo
                enddo
        enddo
      end subroutine

function search_sorted_bisec(arr,val) ! return index i where arr(i) .gt. val
! For example o
!  height   100.00000       90.000000       80.000000       70.000000       60.000000       50.000000       40.000000       30.000000       20.000000       10.000000    
!  lvl           9               8           7           6           5           4           3           2           1
!  search   110.00000     return index           1 which is value   100.00000     at lvl           9
!  search   20.100000     return index           8 which is value   30.000000     at lvl           2
!  search   20.000000     return index           9 which is value   20.000000     at lvl           1
!  search   12.500000     return index           9 which is value   20.000000     at lvl           1
!  search   10.000000     return index           9 which is value   20.000000     at lvl           1

  integer(c_int) :: search_sorted_bisec
  double precision,dimension(:),intent(in) :: arr
  double precision,intent(in)              :: val
  integer(c_int) :: i,j,k
  i=lbound(arr,1)
  j=ubound(arr,1)
  if(arr(i).gt.arr(j)) then
    do
      k=(i+j)/2
      if (val > arr(k)) then
        j=k
      else
        i = k
      endif
      if (i+1 >= j) then ! only single or tuple left
         search_sorted_bisec=i
         exit
      endif
    end do
  else
    do
      k=(i+j)/2
      if (val < arr(k)) then
        j=k
      else
        i = k
      endif
      if (i+1 >= j) then ! only single or tuple left
         search_sorted_bisec=i+1
         exit
      endif
    end do
  endif
end function

end module

!program main
!      use gridtransform
!      implicit none
!
!      integer(c_int),parameter :: Nx=1,Ny=1,Nz=3
!
!      double precision,allocatable :: &
!      kabs(:,:,:) ,&
!      dx(:),dy(:),dz(:)
!
!      allocate(kabs(Nx,Ny,Nz)) ; kabs = 0 ; kabs(1,1,3)=1
!      allocate(dx(Nx) ) ; dx=2700
!      allocate(dy(Ny) ) ; dy=2700
!      allocate(dz(Nz) ) ; dz=200
!
!      print *,'This is gridtransform main'
!
!      call init_grid_transform(dx,dy,dz)
!
!      print *,'old kabs',kabs
!      call grid_old_to_new(kabs)
!      print *,'new kabs',kabs
!
!end program
