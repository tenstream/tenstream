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

      logical, parameter :: luse_memory_map=.True.

      !-----------------------------------------
      !- Define the size of the Lookuptables:  -
      !-----------------------------------------
      !
      ! You should not need to change this... but feel free to play around...
      ! interp_mode 1 == nearest neighbour interpolation
      ! interp_mode 2 == linear interpolation

      integer(iintegers), parameter :: interp_mode_8_10=2
      integer(iintegers), parameter :: interp_mode_1_2=2
      integer(iintegers), parameter :: interp_mode_3_6=2
      integer(iintegers), parameter :: interp_mode_3_10=2
      integer(iintegers), parameter :: interp_mode_wedge_5_8=2

      ! We may also pre-compute the dimensions for the LUT using eddington coeffs as proxy for good values
      !     -- see python script: ''eddington_to_LUT.py''
      ! This way, the output of the script has to be put here... and stuff above commented out

      !real(ireals), parameter :: preset_aspect10(10) = [0.01, 0.56, 1.12, 1.67, 2.23, 2.78, 3.34, 3.89, 4.45, 5.0]
      !real(ireals), parameter :: preset_aspect21(21) = [0.01, 0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2.0, 2.25,  2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.]
      real(ireals), parameter :: preset_aspect23(23) = [0.02, 0.032, 0.042, 0.056, 0.075, 0.1, 0.133, 0.178, 0.237, 0.316, 0.422, 0.562, 0.75, 1., 1.25, 1.562, 1.953, 2.441, 3.052, 3.815, 4.768, 5.96, 7.451]

      !real(ireals), parameter :: preset_tau10(10) = [1e-10,0.0137388296978,0.0573711123567,0.139789935148,0.285526816981,0.610418869282,1.57085127911,2.73919277554,5.41500325409,100.0]
      !real(ireals), parameter :: preset_tau15(15) = [1e-10,9.60528972228e-06,0.00023849409901,0.00166198651975,0.00673267911551,0.0202911873538,0.0511713576613,0.11623968405,0.258968664233,0.568687940597,1.2527762873,2.64735111809,5.53090672088,14.5233716971,100.0]
      !real(ireals), parameter :: preset_tau21(21) = [1e-10,9.01101825167e-08,0.00307841992306,0.0123699155746,0.0282161842775,0.0512938546048,0.0825084072368,0.1236142154,0.17714839488,0.247282753424,0.342437693132,0.485009163821,0.833265411384,1.29248625163,1.72141627953,2.23635589547,2.91744513386,3.92411526823,5.70630418871,10.203209312,100.0]
      !real(ireals), parameter :: preset_tau23(23) = [1e-10,9.01101825167e-08,1e-7,1e-5,0.00307841992306,0.0123699155746,0.0282161842775,0.0512938546048,0.0825084072368,0.1236142154,0.17714839488,0.247282753424,0.342437693132,0.485009163821,0.833265411384,1.29248625163,1.72141627953,2.23635589547,2.91744513386,3.92411526823,5.70630418871,10.203209312,100.0]
      real(ireals), parameter :: preset_tau31(31) = [1e-10,3.62266272998e-07,7.04565803675e-06,4.47545500233e-05,0.000172126759821,0.000495994753047,0.00119161313679,0.00251026980343,0.00480799264297,0.00856221891924,0.0143961482731,0.0231530284254,0.0358868239775,0.0541358315379,0.079959118223,0.11623968405,0.167882053841,0.246414427244,0.350199325489,0.502459974196,0.759082408765,1.08083180518,1.5415157991,2.19832932733,3.04549626819,4.27145477454,6.16953841432,9.43719309835,15.7335501106,29.5819342206,100.0]





      !real(ireals), parameter :: preset_w08(8) = [0.0,0.384011887456,0.643717368476,0.806922184345,0.902876665378,0.949561805714,0.975148725847,0.99999]
      real(ireals), parameter :: preset_w010(10) = [0.0,0.146844960107,0.299864348265,0.441659869071,0.571826424536,0.689506880372,0.793582370219,0.882074852633,0.94968540101,0.99999]
      !real(ireals), parameter :: preset_w015(15) = [0.0,0.206223071171,0.384160332835,0.528227775139,0.643710500558,0.735455592396,0.806918561817,0.861885907586,0.902875161822,0.932179798039,0.949561345869,0.963194640071,0.975148542423,0.987462994037,0.99999]
      real(ireals), parameter :: preset_w020(20) = [0.0,0.152960717624,0.295085090042,0.416951893959,0.521358613652,0.610087211908,0.684967634054,0.747886390181,0.800286677013,0.84336972609,0.878674797098,0.906377786525,0.928097831502,0.943463164595,0.954135786554,0.963824066888,0.972632134967,0.981529289348,0.990759644674,0.99999]


      real(ireals), parameter :: preset_g1(1) = [0.0]
      real(ireals), parameter :: preset_g3(3) = [0.0,0.267018506789,0.5]

      !-----------------------------------------
      !- Define precision of coefficients      -
      !-----------------------------------------
      ! absolute tolerance and relatice tolerance have to be reached for every
      ! coefficient

!      real(ireals),parameter :: stddev_atol=1e-2_ireals
!      real(ireals),parameter :: stddev_atol=1e-3_ireals
!      real(ireals),parameter :: stddev_atol=5e-4_ireals
      real(ireals),parameter :: stddev_atol=5e-4_ireals
!      real(ireals),parameter :: stddev_atol=5e-6_ireals

      real(ireals),parameter :: stddev_rtol=5e-1_ireals
!      real(ireals),parameter :: stddev_rtol=1e-3_ireals

      ! Do some sanity checks on coefficients -- only disable if you are sure
      ! what to expect.
      logical,parameter :: ldebug_optprop=.False.
!      logical,parameter :: ldebug_optprop=.True.

      ! Use delta scaling on optical properties? -- this significantly reduces
      ! the size of the lookuptables.
      logical,parameter :: ldelta_scale=.True.

      ! Treat direct2diffuse radiation in a cone around solar angle as direct
      ! radiation.
!      real(ireals),parameter :: delta_scale_truncate=.9848_ireals ! .9848 = 10 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.9962_ireals ! .9962 = 5 degrees delta scaling
      real(ireals),parameter :: delta_scale_truncate=1.0_ireals   !1.     = 0 degrees delta scaling
!      real(ireals),parameter :: delta_scale_truncate=.8660_ireals ! .8660 = 30 degrees delta scaling


      ! Used to signal that all Angles possible should be loaded when initializing the LUT object -- pass this as azi and zenith
      integer(iintegers), parameter :: OPP_LUT_ALL_ANGLES=361

end module
