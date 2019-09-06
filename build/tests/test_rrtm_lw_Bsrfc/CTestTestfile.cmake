# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/test_rrtm_lw_Bsrfc
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/test_rrtm_lw_Bsrfc
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_test_rrtm_lw_Bsrfc "mpirun" "-n" "2" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_test_rrtm_lw_Bsrfc")
set_tests_properties(pfunit_test_rrtm_lw_Bsrfc PROPERTIES  WORKING_DIRECTORY "/home/c/Carolin.Klinger/tenstream/tests/test_rrtm_lw_Bsrfc")
