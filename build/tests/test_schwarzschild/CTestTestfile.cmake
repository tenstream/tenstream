# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/test_schwarzschild
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/test_schwarzschild
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_test_schwarzschild "mpirun" "-n" "1" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_test_schwarzschild")
