# Install script for directory: /home/c/Carolin.Klinger/tenstream/tests

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/c/Carolin.Klinger/tenstream/build/install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/c/Carolin.Klinger/tenstream/build/tests/eddington/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/interpolation/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/plexrt_error_growth_tracking/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/pprts_error_growth_tracking/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/reorder_mpi_comm/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_LUT_3_10/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_LUT_8_10/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_LUT_wedge_18_8/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_LUT_wedge_5_8/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc_3_10/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc_3_6/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc_3_6_tau_scaling/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc_8_10/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc_8_12/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc_8_16/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_boxmc_8_18/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_convolution/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_fish_plex/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_helper_functions/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_helper_functions_dp/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_icon_plex_utils/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_mmap/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_netcdfio/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_petsc_scatterToZero/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_petsc_sort/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_plex_solution_vecscale/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_pprts_rrtm_icollapse/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_pprts_symmetry/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_rrtm_lw_Bsrfc/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_rrtm_lw_sw/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_schwarzschild/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_search/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_tenstr_atm/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_twostr/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_wedge_boxmc/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_wedge_boxmc_18_8/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_wedge_boxmc_5_8/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_wedge_boxmc_5_8_spherical/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/tests/test_wedge_param_phi/cmake_install.cmake")

endif()

