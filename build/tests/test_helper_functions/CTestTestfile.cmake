# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/test_helper_functions
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/test_helper_functions
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_test_helper_functions "mpirun" "-n" "2" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_test_helper_functions")
