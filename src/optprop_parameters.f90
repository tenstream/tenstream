module optprop_parameters
      use data_parameters,only : ireals,iintegers
      implicit none
      
      logical,parameter :: delta_scale=.True.

      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.9848_ireals ! .9848 = 10 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.8660_ireals ! .8660 = 30 degrees delta scaling
!      real(ireals),parameter :: stddev_rtol=5e-3_ireals
!      real(ireals),parameter :: stddev_rtol=2e-3_ireals
      real(ireals),parameter :: stddev_rtol=1e-2_ireals

!      integer(iintegers) ,parameter :: Ndz=20, Nkabs=40, Nksca=40, Ng=15, Nphi=2, Ntheta=10, interp_mode=2
      integer(iintegers) ,parameter :: Ndz=2, Nkabs=50, Nksca=50, Ng=15, Nphi=2, Ntheta=10, interp_mode=2
!      integer(iintegers) ,parameter :: Ndz=5, Nkabs=5, Nksca=5, Ng=5, Nphi=2, Ntheta=10, interp_mode=2 
!      integer(iintegers) ,parameter :: Ndz=2, Nkabs=2, Nksca=2, Ng=2, Nphi=2, Ntheta=10, interp_mode=2 

end module
