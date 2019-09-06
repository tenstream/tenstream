# CMake generated Testfile for 
# Source directory: /home/c/Carolin.Klinger/tenstream/tests/reorder_mpi_comm
# Build directory: /home/c/Carolin.Klinger/tenstream/build/tests/reorder_mpi_comm
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pfunit_reorder_mpi_comm "mpirun" "-n" "9" "/home/c/Carolin.Klinger/tenstream/build/bin/pfunit_reorder_mpi_comm")
