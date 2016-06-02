#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

#add_executable( rrtm_example_lw fpda_example_lw.c )
#target_link_libraries( rrtm_example_lw m tenstr_rrtm_lw )
#
#add_executable( rrtm_example_sw fpda_example_sw.c )
#target_link_libraries( rrtm_example_sw m tenstr_rrtm_sw )

add_executable( tenstream_example_sw tenstream_example.f90 )
target_link_libraries( tenstream_example_sw tenstr_rrtm_sw tenstream)


#install(TARGETS rrtm_example_lw DESTINATION bin)
#install(TARGETS rrtm_example_sw DESTINATION bin)
