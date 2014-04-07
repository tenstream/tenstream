module helper_functions
      use data_parameters,only : ireals,pi,one

      implicit none
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
          factor = epsilon(b)
        endif
        if( a.le.b+factor .and. a.ge.b-factor ) then
          approx = .True.
        else
          approx = .False.
        endif
      end function


      end module
