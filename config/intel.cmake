# Default intel compiler flags
set(USER_C_FLAGS       "-W -Wall -mkl=sequential")
set(USER_Fortran_FLAGS "-g -sox -no-wrap-margin -mkl=sequential -warn all,noexternal")
set(USER_Fortran_FLAGS_RELEASE "-mtune=native -O3 -ftz -pc64 -fp-model fast=2 -no-prec-sqrt -fast-transcendentals")
set(USER_Fortran_FLAGS_DEBUG "-traceback -extend_source -g -fp-model strict -ftrapuv -warn errors -fpe0 -O2 -g -check all -check nopointers -check noarg_temp_created ")
