# Default GCC
#
# includes default compile flags for gcc type compilers

message(STATUS "### USING GCC CONFIG ###")

set(USER_C_FLAGS               "-W -Wall")
set(USER_Fortran_FLAGS         "-W -Wall --std=legacy" )
set(USER_Fortran_FLAGS_RELEASE "-fno-backtrace -fno-range-check -O3 -mtune=native")
set(USER_Fortran_FLAGS_DEBUG   "-Og -ggdb -fbacktrace -fcheck=all -fbounds-check -Wsurprising -Wuninitialized -ffpe-trap=invalid")
