# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/plexrt_error_growth_tracking
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/plexrt_error_growth_tracking
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_plexrt_error_growth_tracking "mpirun" "-n" "2" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_plexrt_error_growth_tracking")
set_tests_properties(pfunit_plexrt_error_growth_tracking PROPERTIES  WORKING_DIRECTORY "/home/c/Carolin.Klinger/tenstream/tests/plexrt_error_growth_tracking")
