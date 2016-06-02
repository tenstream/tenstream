#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

add_executable( tenstream_example_sw tenstream_example.f90 )
target_link_libraries( tenstream_example_sw tenstr_rrtm_sw tenstream)
