#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file( GLOB files *.f90 )

SET_SOURCE_FILES_PROPERTIES(rrtmg_lw_k_g.f90 PROPERTIES COMPILE_FLAGS -O0)

add_library( fpda_rrtm_lw STATIC ${files} )
target_link_libraries( fpda_rrtm_lw )

add_executable( fpda_example_lw fpda_example.F90 )
target_link_libraries( fpda_example_lw fpda_rrtm_lw )

install(TARGETS fpda_rrtm_lw   DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
