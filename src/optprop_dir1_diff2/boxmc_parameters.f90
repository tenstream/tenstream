module boxmc_parameters_1_2
      use data_parameters,only : ireals,iintegers
      implicit none
      
      integer(iintegers),parameter :: dir_streams=1, diff_streams=2
      logical,parameter :: delta_scale=.True.

      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling

      integer(iintegers) ,parameter :: Ndz=10, Nkabs=40, Nksca=40, Ng=10, Nphi=2, Ntheta=10, interp_mode=2 
!      integer(iintegers) ,parameter :: Ndz=5, Nkabs=5, Nksca=5, Ng=5, Nphi=2, Ntheta=10, interp_mode=2 

end module
