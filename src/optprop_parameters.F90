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
      use m_data_parameters,only : irealLUT, iintegers, default_str_len
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

      character(default_str_len) :: lut_basename='./LUT'

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
      integer(iintegers), parameter :: LUT_MAX_DIM=8

      ! We pre-compute the dimensions for the LUT using eddington coeffs as proxy for good values
      !     -- see python script: ''eddington_to_LUT.py''
      real(irealLUT), parameter :: param_eps = 1e-4
      real, parameter :: peps_r = real(param_eps)

      real(irealLUT), parameter :: preset_param_phi6(6) = &
        [-2.      , -1.-peps_r, -1.+peps_r, &
        +1.-peps_r, +1.+peps_r, 2.]

      real(irealLUT), parameter :: preset_param_phi11(11) = &
        [-2.      , -1.5      , &
        -1.-peps_r, -1.+peps_r, &
        -0.5      , 0.        , 0.5, &
        +1.-peps_r, +1.+peps_r, &
        1.5       , 2.]

      real(irealLUT), parameter :: preset_param_phi19(19) = [ &
        -2.       , -1.75     , -1.5 , -1.25, &
        -1.-peps_r, -1.+peps_r, &
        -0.75     , -0.5      , -0.25, 0.   , 0.25, 0.5, 0.75, &
        +1.-peps_r, +1.+peps_r, &
        1.25      , 1.5       , 1.75 , 2.]

      real(irealLUT), parameter :: preset_param_phi83(83) = [ &
       -2.        , -1.95     , -1.9      , -1.85     , -1.8 , -1.75, -1.7 , -1.65, -1.6 , &
       -1.55      , -1.5      , -1.45     , -1.4      , -1.35, -1.3 , -1.25, -1.2 , -1.15, &
       -1.1       , -1.05     , -1.-peps_r, -1.+peps_r, &
       -0.95      , -0.9      , -0.85     , -0.8      , -0.75, -0.7 , &
       -0.65      , -0.6      , -0.55     , -0.5      , -0.45, -0.4 , -0.35, -0.3 , -0.25, &
       -0.2       , -0.15     , -0.1      , -0.05     , 0.   , 0.05 , 0.1  , 0.15 , 0.2  , &
        0.25      , 0.3       , 0.35      , 0.4       , 0.45 , 0.5  , 0.55 , 0.6  , 0.65 , &
        0.7       , 0.75      , 0.8       , 0.85      , 0.9  , 0.95 , &
        +1.-peps_r, +1.+peps_r, 1.05      , 1.1       , &
        1.15      , 1.2       , 1.25      , 1.3       , 1.35 , 1.4  , 1.45 , 1.5  , 1.55 , &
        1.6       , 1.65      , 1.7       , 1.75      , 1.8  , 1.85 , 1.9  , 1.95 , 2. ]


      real(irealLUT), parameter :: preset_param_theta4(4) = [ &
        -1., -peps_r, +peps_r, 1.]

      real(irealLUT), parameter :: preset_param_theta13(13) = [ &
        -1., -peps_r, +peps_r, &
        .1 , .2     , .3     , .4, .5, .6, .7, .8, .9, 1.]


      real(irealLUT), parameter :: preset_aspect5(5) = [.2, .5, 1., 2., 5.]
      real(irealLUT), parameter :: preset_aspect7(7) = &
        [ 0.03125, 0.0625, 0.125, 0.25, 0.5, 1., 2.]
      real(irealLUT), parameter :: preset_aspect11(11) = [ &
        0.125     , 0.16493849, 0.21763764, 0.28717459, 0.37892914, &
        0.5       , 0.65975396, 0.87055056, 1.14869835, 1.51571657, &
        2.        ]
      real(irealLUT), parameter :: preset_aspect13(13) = [ &
        0.02 , 0.042, 0.075, 0.133 , 0.237, 0.422, 0.75, 1., 1.25, &
        1.953, 3.052, 4.768, 7.451]
      real(irealLUT), parameter :: preset_aspect17(17) = [ &
        0.02 , 0.032, 0.042, 0.056, 0.075, 0.1 , 0.133, 0.178, 0.237, &
        0.316, 0.422, 0.562, 0.75 , 1.   , 1.25, 1.562, 2.]
      real(irealLUT), parameter :: preset_aspect18(18) = [ &
        0.02 , 0.032, 0.042, 0.056, 0.075, 0.1 , 0.133, 0.178, 0.237, &
        0.316, 0.422, 0.562, 0.75 , 1.   , 1.25, 1.562, 2.   , 3.]
      real(irealLUT), parameter :: preset_aspect23(23) = [ &
        0.02 , 0.032, 0.042, 0.056, 0.075 , 0.1 , 0.133, 0.178, 0.237, &
        0.316, 0.422, 0.562, 0.75 , 1.    , 1.25, 1.562, 1.953, 2.441, &
        3.052, 3.815, 4.768, 5.96 , 7.451]

      real(irealLUT), parameter :: preset_tau15(15) = &
        [1e-10          , 9.60528972228e-06, 0.00023849409901, 0.00166198651975, &
        0.00673267911551, 0.0202911873538  , 0.051171357661  , 0.116239684     , &
        0.258968664233  , 0.568687940597   , 1.2527762873    , 2.64735111809   , &
        5.53090672088   , 14.5233716971    , 100.0]

      real(irealLUT), parameter :: preset_tau20(20) = &
        [1e-10          , 2.33773213401e-06, 5.40185638224e-05, 0.000365962943669, &
        0.00145415861897, 0.00431514105527 , 0.0105306225135  , &
        0.0225104907999 , 0.044534085216   , 0.0835690735283  , &
        0.152160041198  , 0.271322429414   , 0.492503225042   , &
        0.91860742252   , 1.60959133986    , 2.79337830498    , 4.89077663742    , &
        9.35922562367   , 21.643468069     , 100.0]

      real(irealLUT), parameter :: preset_tau31(31) = &
        [1e-10           , 3.62266272998e-07, 7.04565803675e-06, 4.47545500233e-05, &
        0.000172126759821, 0.000495994753047, 0.00119161313679 , &
        0.00251026980343 , 0.00480799264297 , 0.00856221891924 , &
        0.0143961482731  , 0.0231530284254  , 0.0358868239775  , &
        0.0541358315379  , 0.079959118223   , 0.11623968405    , 0.167882053841   , &
        0.246414427244   , 0.350199325489   , 0.502459974196   , 0.759082408765   , &
        1.08083180518    , 1.5415157991     , 2.19832932733    , 3.04549626819    , &
        4.27145477454    , 6.16953841432    , 9.43719309835    , 15.7335501106    , &
        29.5819342206    , 100.0]

      real(irealLUT), parameter :: preset_w010(10) = &
        [0.0          , 0.146844960107, 0.299864348265, 0.441659869071, 0.571826424536, &
        0.689506880372, 0.793582370219, 0.882074852633, 0.94968540101 , 0.99999]

      real(irealLUT), parameter :: preset_w020(20) = &
        [0.0          , 0.152960717624, 0.295085090042, 0.416951893959, 0.521358613652, &
        0.610087211908, 0.684967634054, 0.747886390181, 0.800286677013, &
        0.84336972609 , 0.878674797098, 0.906377786525, 0.928097831502, &
        0.943463164595, 0.954135786554, 0.963824066888, 0.972632134967, &
        0.981529289348, 0.990759644674, 0.99999]

      real(irealLUT), parameter :: preset_w021(21) = &
        [0.0          , 0.152960717624, 0.295085090042, 0.416951893959, 0.521358613652, &
        0.610087211908, 0.684967634054, 0.747886390181, 0.800286677013, &
        0.84336972609 , 0.878674797098, 0.906377786525, 0.928097831502, &
        0.943463164595, 0.954135786554, 0.963824066888, 0.972632134967, &
        0.981529289348, 0.990759644674, 0.99999, 1.0]

      real(irealLUT), parameter :: preset_g2(2) = [0.0, 0.5]
      real(irealLUT), parameter :: preset_g3(3) = [0.0, 0.2670, 0.5]
      real(irealLUT), parameter :: preset_g4(4) = [0.0, 0.1936, 0.4258, 0.65]
      real(irealLUT), parameter :: preset_g6(6) = [0.0, 0.2424, 0.4137, 0.5717, 0.7144, 0.85]

      !-----------------------------------------
      !- Define precision of coefficients      -
      !-----------------------------------------
      ! absolute tolerance and relatice tolerance have to be reached for every
      ! coefficient

!      real(irealLUT),parameter :: stddev_atol=1e-2_irealLUT
!      real(irealLUT),parameter :: stddev_atol=5e-3_irealLUT
      real(irealLUT) :: stddev_atol=1e-3_irealLUT
!      real(irealLUT),parameter :: stddev_atol=2e-4_irealLUT
!      real(irealLUT),parameter :: stddev_atol=5e-6_irealLUT

      real(irealLUT) :: stddev_rtol=2e-1_irealLUT
!      real(irealLUT),parameter :: stddev_rtol=1e-3_irealLUT

      ! Do some sanity checks on coefficients -- only disable if you are sure
      ! what to expect.
      ! logical,parameter :: ldebug_optprop=.False.
      logical,parameter :: ldebug_optprop=.True.

      ! Use delta scaling on optical properties? -- this significantly reduces
      ! the size of the lookuptables.
      logical,parameter :: ldelta_scale=.True.

      ! Treat direct2diffuse radiation in a cone around solar angle as direct
      ! radiation.
!      real(irealLUT),parameter :: delta_scale_truncate=.9848_irealLUT ! .9848 = 10 degrees delta scaling
!      real(irealLUT),parameter :: delta_scale_truncate=.9962_irealLUT ! .9962 = 5 degrees delta scaling
      real(irealLUT),parameter :: delta_scale_truncate=1.0_irealLUT   !1.     = 0 degrees delta scaling
!      real(irealLUT),parameter :: delta_scale_truncate=.8660_irealLUT ! .8660 = 30 degrees delta scaling

      ! spherical correction for wedge computations,
      ! this is tuned towards earth radius and average dx = 1000m sized elements
      real(irealLUT), parameter :: wedge_sphere_radius = 6371e3_irealLUT / 1000._irealLUT

      real(irealLUT), parameter :: LUT_dump_interval=3600 ! dump the LUT every 60 minutes
      real(irealLUT), parameter :: LUT_max_create_jobtime=3600*3 ! after 3hrs, cancel the createLUT jobs in any case
end module
