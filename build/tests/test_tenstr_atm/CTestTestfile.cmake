# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/test_tenstr_atm
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/test_tenstr_atm
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_test_tenstr_atm "mpirun" "-n" "2" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_test_tenstr_atm")
set_tests_properties(pfunit_test_tenstr_atm PROPERTIES  WORKING_DIRECTORY "/home/c/Carolin.Klinger/tenstream/tests/test_tenstr_atm")
