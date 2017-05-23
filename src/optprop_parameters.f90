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

module m_optprop_parameters
      use m_data_parameters,only : ireals, iintegers, default_str_len
      implicit none

      !> \page optprop_parameters Parameters concerning the transport coefficients
      !! You should define some parameters about how and where to find the Look
      !! Up Tables foir the transport coefficients
      !> Have a look at the options in m_optprop_parameters::

      !-----------------------------------------
      !> Define the path to the Lookuptables
      !!
      !!  This has to be a reachable path for rank 0,
      !!  At MIM in Munich please set to
      !!  '/home/opt/cosmo_tica_lib/tenstream/optpropLUT/LUT'
      !!  At ZMAW in Hamburg please set to
      !!  '/scratch/mpi/mpiaes/m300362/tenstream_LUT/LUT'
      !-----------------------------------------

      character(default_str_len) :: lut_basename='./LUT' !'/scratch/mpi/mpiaes/m300362/tenstream_LUT/LUT'

      !-----------------------------------------
      !> Define the mode to calculate coeffs   -
      !!-----------------------------------------
      !! 0 :: retrieve directly from Lookuptable
      !! 1 :: retrieve from ANN (experimental) ::
      !!        this assumes you actually have a
      !!        ANN that can be queried....
      !!
      integer(iintegers) :: coeff_mode = 0

      !-----------------------------------------
      !- Define the size of the Lookuptables:  -
      !-----------------------------------------
      !
      ! You should not need to change this... but feel free to play around...
      ! interp_mode 1 == nearest neighbour interpolation
      ! interp_mode 2 == linear interpolation, nearest neighbour in solar angles
      ! interp_mode 3 == linear interpolation, nearest neighbour in solar azimuth
      ! interp_mode 4 == linear interpolation in all dimensions

      integer(iintegers) ,parameter :: Ndiff_8_10=10, Ndir_8_10=8, interp_mode_8_10=2

      integer(iintegers) ,parameter :: Ndiff_1_2=2, Ndir_1_2=1, interp_mode_1_2=2

      integer(iintegers) ,parameter :: Nphi=10

      ! integer(iintegers) ,parameter :: Ntau=40, Nw0=10, Ng=2, Ntheta=19
      ! real(ireals), parameter :: preset_tau(1)  =[0]
      ! real(ireals), parameter :: preset_w0(1)   =[0]
      ! real(ireals), parameter :: preset_g(1)    =[0]
      ! real(ireals), parameter :: preset_theta(1)=[0]


      ! We may also pre-compute the dimensions for the LUT using eddington coeffs as proxy for good values
      !     -- see python script: ''eddington_to_LUT.py''
      ! This way, the output of the script has to be put here... and stuff above commented out
			logical, parameter :: use_prescribed_LUT_dims=.True.

!      integer(iintegers), parameter :: Ntau=102, Nw0=10, Ng=3, Ntheta=19
!			real(ireals), parameter :: preset_tau(102) = [1e-10,9.01101825167e-08,0.00011747647831,0.000459246031624,0.00102551144299,0.00181774047313,0.0028363788343,0.00408142349113,0.00555913729201,0.00726636713699,0.00920009297672,0.0113863053351,0.0137887340035,0.016460536223,0.0193725050106,0.0225314004723,0.0259452600421,0.0296199706257,0.0335584593541,0.0377596330187,0.042216976843,0.0470178466654,0.0521145731098,0.0574715388266,0.0631079017627,0.0691810406145,0.0754928415292,0.0821576289514,0.0892558929205,0.0965281440447,0.104479776842,0.112580022708,0.121315425255,0.13031049739,0.139936937786,0.14989146143,0.16054791759,0.171530409588,0.183390600761,0.195474424896,0.208759361254,0.222044297612,0.237021757979,0.252018597476,0.268652951636,0.285744317294,0.304295186079,0.324007074608,0.344867727746,0.367953056571,0.391787710049,0.419384285043,0.447466251539,0.481433648365,0.516630578919,0.560407048941,0.611045323274,0.672590979942,0.785438933925,0.918783237441,1.01081970932,1.09279951321,1.17318537381,1.25137398969,1.32948728738,1.40929582745,1.48974842634,1.57150971865,1.65781955673,1.74412939482,1.83554449403,1.93122786266,2.02691123129,2.13294581641,2.24210089986,2.35125598332,2.47854359759,2.60634899337,2.73997183085,2.8932442289,3.04651662695,3.22054928104,3.40843738858,3.60271194545,3.8375836044,4.07245526336,4.36284253396,4.66142220682,5.03114165197,5.4160118427,5.9174258949,6.48042592743,7.13917536419,8.00693271984,9.09659402734,10.5288734794,12.5111466027,15.6698358077,21.2104809095,35.4308085013,100.0,100.0]
!			real(ireals), parameter :: preset_w0(10) = [0.0,0.309611541089,0.54233191301,0.707203158386,0.820526619394,0.894993804547,0.939905081045,0.961866428699,0.980503546279,0.99999]
!			real(ireals), parameter :: preset_g(3) = [0.0,0.267018506789,0.5]
!			real(ireals), parameter :: preset_theta(19) = [0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0]

      integer(iintegers), parameter :: Ntau=22, Nw0=15, Ng=3, Ntheta=19
      real(ireals), parameter :: preset_tau(22) = [1e-10,9.01101825167e-08,0.00307841992306,0.0123699155746,0.0282161842775,0.0512938546048,0.0825084072368,0.1236142154,0.17714839488,0.247282753424,0.342437693132,0.485009163821,0.833265411384,1.29248625163,1.72141627953,2.23635589547,2.91744513386,3.92411526823,5.70630418871,10.203209312,100.0,100.0]
      real(ireals), parameter :: preset_w0(15) = [0.0,0.206223071171,0.384160332835,0.528227775139,0.643710500558,0.735455592396,0.806918561817,0.861885907586,0.902875161822,0.932179798039,0.949561345869,0.963194640071,0.975148542423,0.987462994037,0.99999]
      real(ireals), parameter :: preset_g(3) = [0.0,0.267018506789,0.5]
      real(ireals), parameter :: preset_theta(19) = [0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0]


      !-----------------------------------------
      !- Define precision of coefficients      -
      !-----------------------------------------
      ! absolute tolerance and relatice tolerance have to be reached for every
      ! coefficient

!      real(ireals),parameter :: stddev_atol=1e-2_ireals
!      real(ireals),parameter :: stddev_atol=5e-3_ireals
      real(ireals),parameter :: stddev_atol=1e-3_ireals
!      real(ireals),parameter :: stddev_atol=1e-4_ireals
!      real(ireals),parameter :: stddev_atol=5e-6_ireals

      real(ireals),parameter :: stddev_rtol=1e-2_ireals
!      real(ireals),parameter :: stddev_rtol=1e-3_ireals
!      real(ireals),parameter :: stddev_rtol=1e-4_ireals

      ! Do some sanity checks on coefficients -- only disable if you are sure
      ! what to expect.
!      logical,parameter :: ldebug_optprop=.False.
      logical,parameter :: ldebug_optprop=.True.

      ! Use delta scaling on optical properties? -- this significantly reduces
      ! the size of the lookuptables.
      logical,parameter :: ldelta_scale=.True.

      ! Treat direct2diffuse radiation in a cone around solar angle as direct
      ! radiation.
!      real(ireals),parameter :: delta_scale_truncate=.9848_ireals ! .9848 = 10 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling
      real(ireals),parameter :: delta_scale_truncate=1.0_ireals   !1.     = 0 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.8660_ireals ! .8660 = 30 degrees delta scaling


end module
