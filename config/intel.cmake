# Default intel compiler flags
set(USER_C_FLAGS       "-W -Wall -qmkl=sequential")
set(USER_Fortran_FLAGS "-g -no-wrap-margin -qmkl=sequential -warn all,noexternal")
set(USER_Fortran_FLAGS_RELEASE "-mtune=native -O3 -ftz -fp-model fast=2")
set(USER_Fortran_FLAGS_DEBUG "-traceback -g -fp-model strict -ftrapuv -warn errors -fpe0 -O2 -g -check all -check nopointers -check noarg_temp_created ")
