include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( general_sources
  "plex_grid.F90"
  )

add_library( plexrt ${general_sources} )
target_link_libraries( plexrt tenstream)

add_executable( gen_plex_from_icon gen_plex_from_icon.F90 )
target_link_libraries( gen_plex_from_icon plexrt tenstream )

install(TARGETS gen_plex_from_icon DESTINATION bin)


add_executable( plex_ex1 plex_ex1.F90 )
target_link_libraries( plex_ex1 plexrt tenstream )

install(TARGETS plex_ex1 DESTINATION bin)

install(TARGETS plexrt DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
