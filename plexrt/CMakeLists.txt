include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file( GLOB files *.F* )

add_library( plexrt STATIC ${files} )
target_link_libraries( plexrt tenstream)

add_executable( gen_plex_from_icon gen_plex_from_icon.F90 )
target_link_libraries( gen_plex_from_icon plexrt tenstream )

install(TARGETS gen_plex_from_icon DESTINATION bin)

install(TARGETS plexrt DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
