# Default intel compiler flags
set(USER_C_FLAGS       "-W -Wall")
set(USER_Fortran_FLAGS "-g -sox -no-wrap-margin -mkl=sequential -warn all")
set(USER_Fortran_FLAGS_RELEASE "-mtune=native -Ofast -ftz -pc64 -fp-model fast=2 -no-prec-div -no-prec-sqrt -fast-transcendentals")
set(USER_Fortran_FLAGS_DEBUG "-traceback -extend_source -g -fp-model strict -ftrapuv -warn all -warn errors -fpe0 -O2 -g -check all -check nopointers -check noarg_temp_created ")
