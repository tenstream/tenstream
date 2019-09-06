# Install script for directory: /home/c/Carolin.Klinger/tenstream/examples

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
  include("/home/c/Carolin.Klinger/tenstream/build/examples/box_cld/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/libRadtran_cld_file/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/plexrt/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/pprts/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/pprts_ex1/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/pprts_hill/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/python/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/rrtm_lw_sw/cmake_install.cmake")
  include("/home/c/Carolin.Klinger/tenstream/build/examples/wetterstein_ts/cmake_install.cmake")

endif()

