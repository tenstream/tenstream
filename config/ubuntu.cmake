# UBUNTU CMAKE config
#
# example Config script for compile under ubuntu
#
# Make sure you have set PETSC_ARCH and PETSC_DIR environment variables

message(STATUS "### Using ubuntu.cmake config ###")

set(USER_C_FLAGS               "-W -Wall -Wuninitialized --std=c99")
set(USER_Fortran_FLAGS         "-W -Wall -g")
set(USER_Fortran_FLAGS_RELEASE "-fno-backtrace -fno-range-check -O3")
set(USER_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=nan -Wall -Werror -Wuninitialized -g -pg -fcheck=all -fbounds-check -pedantic -Wsurprising -ffpe-trap=invalid")
