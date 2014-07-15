module optprop_parameters
      use data_parameters,only : ireals,iintegers
      implicit none
      
!      logical,parameter :: delta_scale=.True.
      logical,parameter :: delta_scale=.False.

      real(ireals),parameter :: delta_scale_truncate=1.0_ireals ! .9962 = 5 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.9848_ireals ! .9848 = 10 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.8660_ireals ! .8660 = 30 degrees delta scaling
      real(ireals),parameter :: stddev_rtol=5e-3_ireals
!      real(ireals),parameter :: stddev_rtol=2e-3_ireals
!      real(ireals),parameter :: stddev_rtol=1e-2_ireals

!      integer(iintegers) ,parameter :: Ndz=20, Nkabs=40, Nksca=40, Ng=15, Nphi=2, Ntheta=10, interp_mode=2
      integer(iintegers) ,parameter :: Ndz_8_10=3, Nkabs_8_10=80, Nksca_8_10=80, Ng_8_10=15, Nphi_8_10=10, Ntheta_8_10=10, interp_mode_8_10=2
!      integer(iintegers) ,parameter :: Ndz=5, Nkabs=5, Nksca=5, Ng=5, Nphi=2, Ntheta=10, interp_mode=2 
!      integer(iintegers) ,parameter :: Ndz=2, Nkabs=2, Nksca=2, Ng=2, Nphi=2, Ntheta=10, interp_mode=2 

      integer(iintegers) ,parameter :: Ndz_1_2=30, Nkabs_1_2=30, Nksca_1_2=30, Ng_1_2=15, Nphi_1_2=10, Ntheta_1_2=10, interp_mode_1_2=2
end module
