# Default GCC
#
# Config script for the linux environment at the Meteorological Institute Munich 08.Jan.2015
#
# There may be a maintained(Fabian Jakub) petsc library available through
#
#  module load petsc
#
# if so, you might get away and use this config file here with
# `cmake <tenstream_root_dir> -DSYST:STRING=lmu_mim`

set(USER_C_FLAGS               "-W -std=c99")
set(USER_Fortran_FLAGS         "-g")
set(USER_Fortran_FLAGS_RELEASE "-fno-range-check -O3 -mtune=native")
set(USER_Fortran_FLAGS_DEBUG   "-fbacktrace -finit-real=snan -Werror -Wall -pedantic -g -pg -fcheck=all -fbounds-check -Wuninitialized -ffpe-trap=invalid ")
