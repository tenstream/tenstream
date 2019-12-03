# Default GCC
#
# Config script for the linux environment at the Meteorological Institute Munich 08.Jan.2015
#
# and use this config file with `cmake <tenstream_root_dir> -DSYST:STRING=lmu_mim`

set(USER_C_FLAGS       "-Wall -std=c99")
set(USER_Fortran_FLAGS "-g -cpp -sox -no-wrap-margin -mkl=sequential -warn all")
set(USER_Fortran_FLAGS_RELEASE "-xCORE-AVX2 -march=native -Ofast -ftz -pc64 -fp-model fast=2 -no-prec-div -no-prec-sqrt -fast-transcendentals")
set(USER_Fortran_FLAGS_DEBUG "-traceback -extend_source -g -fp-model strict -ftrapuv -warn all -warn errors -fpe0 -O2 -g -check all -check nopointers -check noarg_temp_created ")
