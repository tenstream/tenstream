module helper_functions
      use data_parameters,only : ireals,pi

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

      elemental logical function approx(a,b)
        real(ireals),intent(in) :: a,b
        if( a.le.b+1e3*epsilon(b) .and. a.ge.b-1e3*epsilon(b) ) then
          approx = .True.
        else
          approx = .False.
        endif
      end function


      end module
