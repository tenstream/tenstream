# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/eddington
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/eddington
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_eddington "mpirun" "-n" "1" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_eddington")
