message(STATUS "Installing pfUnit  --  Installation taken from https://github.com/bast/pfunit-demo/blob/master/CMakeLists.txt")

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)



message(STATUS "Installing fypp - Fortran Preprocessor")
#ExternalProject_Add(fypp
#  GIT_REPOSITORY https://github.com/aradi/fypp.git
#  SOURCE_DIR "${PROJECT_BINARY_DIR}/external/fypp"
#  CONFIGURE_COMMAND ""
#  BUILD_COMMAND ""
#  INSTALL_COMMAND ""
#)
set(FYPP_EXE ${PROJECT_SOURCE_DIR}/external/fypp CACHE STRING "Python Script Location for Fortran PreProcessor")

