# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/test_petsc_scatterToZero
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/test_petsc_scatterToZero
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_test_petsc_scatterToZero "mpirun" "-n" "3" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_test_petsc_scatterToZero")
set_tests_properties(pfunit_test_petsc_scatterToZero PROPERTIES  WORKING_DIRECTORY "/home/c/Carolin.Klinger/tenstream/tests/test_petsc_scatterToZero")
