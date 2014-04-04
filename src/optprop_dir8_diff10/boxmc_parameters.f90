module boxmc_parameters_8_10
      use data_parameters,only : ireals,iintegers
      implicit none
      
      integer(iintegers),parameter :: dir_streams=8, diff_streams=10
      logical,parameter :: delta_scale=.True.

      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling

      integer(iintegers) ,parameter :: Ndz=20, Nkabs=40, Nksca=40, Ng=15, Nphi=2, Ntheta=10, interp_mode=2 
!      integer(iintegers) ,parameter :: Ndz=5, Nkabs=5, Nksca=5, Ng=5, Nphi=2, Ntheta=10, interp_mode=2 
!      integer(iintegers) ,parameter :: Ndz=2, Nkabs=2, Nksca=2, Ng=2, Nphi=2, Ntheta=10, interp_mode=2 

end module
