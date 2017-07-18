# UBUNTU CMAKE config 
#
# example Config script for compile under ubuntu
#
# Make sure you have set PETSC_ARCH and PETSC_DIR environment variables

message(STATUS "### Using ubuntu.cmake config ###")

set(CMAKE_C_COMPILER   "mpicc")
set(CMAKE_Fortran_COMPILER   "/usr/bin/mpif90")
set(Fortran_COMPILER_WRAPPER "/usr/bin/mpif90")

set(USER_C_FLAGS               "-cpp -W -Wall -Wuninitialized --std=c99") 
set(USER_Fortran_FLAGS         "-cpp -ffree-line-length-none -W -Wall -g")
set(USER_Fortran_FLAGS_RELEASE "-fno-backtrace -fno-range-check -O3") 
set(USER_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=nan -W -Wall -Wuninitialized -g -pg -fcheck=all -fbounds-check -pedantic -Wsurprising -fno-range-check")
