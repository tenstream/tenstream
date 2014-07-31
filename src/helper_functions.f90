module m_helper_functions
      use m_data_parameters,only : iintegers,ireals,pi,one,imp_real,imp_int,mpiint,imp_comm

      implicit none

      public imp_bcast

      interface imp_bcast
        module procedure imp_bcast_real_1d,imp_bcast_real_3d
      end interface

      integer(mpiint) :: mpierr

      contains

      pure function norm(v)
        real(ireals) :: norm
        real(ireals),intent(in) :: v(:)
        norm = sqrt(dot_product(v,v))
      end function

      elemental function deg2rad(deg)
        real(ireals) :: deg2rad
        real(ireals),intent(in) :: deg
        deg2rad = deg *pi/180._ireals
      end function

      pure function rmse(a,b)
          real(ireals) :: rmse(2)
          real(ireals),intent(in) :: a(:),b(:)
          rmse(1) = sqrt( mean( (a-b)**2 ) )
          rmse(2) = rmse(1)/max( mean(b), epsilon(rmse) )
      end function

      pure function mean(arr)
          real(ireals) :: mean
          real(ireals),intent(in) :: arr(:)
          mean = sum(arr)/size(arr)
      end function

      elemental logical function approx(a,b,precision)
        real(ireals),intent(in) :: a,b
        real(ireals),intent(in),optional :: precision
        real(ireals) :: factor
        if(present(precision) ) then
          factor = precision
        else
          factor = 10*epsilon(b)
        endif
        if( a.le.b+factor .and. a.ge.b-factor ) then
          approx = .True.
        else
          approx = .False.
        endif
      end function

      subroutine  imp_bcast_real_1d(arr,sendid,myid)
          real(ireals),allocatable,intent(inout) :: arr(:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot

          if(sendid.eq.myid) Ntot = size(arr)
          call mpi_bcast(Ntot,1,imp_int,sendid,imp_comm,mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot) )
          call mpi_bcast(arr,size(arr),imp_real,sendid,imp_comm,mpierr)
      end subroutine
      subroutine  imp_bcast_real_3d(arr,sendid,myid)
          real(ireals),allocatable,intent(inout) :: arr(:,:,:)
          integer(mpiint),intent(in) :: sendid,myid

          integer(iintegers) :: Ntot(3)

          if(sendid.eq.myid) Ntot = shape(arr)
          call mpi_bcast(Ntot,3,imp_int,sendid,imp_comm,mpierr)

          if(myid.ne.sendid) allocate( arr(Ntot(1), Ntot(2), Ntot(3) ) )
          call mpi_bcast(arr,size(arr),imp_real,sendid,imp_comm,mpierr)
      end subroutine

elemental subroutine delta_scale( kabs,ksca,g ) 
          real(ireals),intent(inout) :: kabs,ksca,g ! kabs, ksca, g
          real(ireals) :: dtau, w0
          dtau = kabs+ksca
          w0   = ksca/dtau
          g    = g
          call delta_scale_optprop( dtau, w0, g)
          kabs= dtau * (one-w0)
          ksca= dtau * w0
      end subroutine
elemental subroutine delta_scale_optprop( dtau, w0, g) 
          real(ireals),intent(inout) :: dtau,g,w0
          real(ireals) :: f

          f = g**2
!          f = g
          dtau = dtau * ( one - w0 * f )
          g    = ( g - f ) / ( one - f )
          w0   = w0 * ( one - f ) / ( one - f * w0 )
      end subroutine

      end module
