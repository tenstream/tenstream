module boxmc_parameters
      use data_parameters,only : ireals,iintegers
      implicit none
      
      integer(iintegers),parameter :: dir_streams=1, diff_streams=2
      logical,parameter :: delta_scale=.True.

      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling

      integer(iintegers) ,parameter :: Ndz=20, Nkabs=30, Nksca=30, Ng=20, Nphi=2, Ntheta=10, interp_mode=2 

end module
