# Default GCC
#
# example Config script for a gcc linux environment 08.Jan.2015
#
# You need to provide paths for the NETCDF installation
# and additionally, you need to provide an environment variable so that we can find PETSC
# Make sure you have set PETSC_ARCH and PETSC_DIR
# If unclear, see the PETSC installation instructions and read the README for further hints.
#

message(STATUS "### USING GCC CONFIG ###")

set(USER_C_FLAGS               "-cpp -W -Wall -Wuninitialized --std=c99")
set(USER_Fortran_FLAGS         "-cpp -W -Wall -Wuninitialized -g")
set(USER_Fortran_FLAGS_RELEASE "-fno-backtrace -fno-range-check -O3")
set(USER_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=nan -pg -fcheck=all -fbounds-check -pedantic -Wsurprising")
