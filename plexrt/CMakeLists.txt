include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( general_sources
  "plex_grid.F90"
  "plex_rt.F90"
  "plex2rayli.F90"
  "icon_grid.F90"
  "icon_plex_utils.F90"
  )

add_library( plexrt ${general_sources} )
target_link_libraries(plexrt pprts_base ${LIB} ${EXTLIB})

add_executable( gen_plex_from_icon gen_plex_from_icon.F90 )
target_link_libraries( gen_plex_from_icon plexrt pprts_base ${LIB} ${EXTLIB} )
install(TARGETS gen_plex_from_icon DESTINATION bin)

add_executable( gen_fish_plex gen_fish_plex.F90 )
target_link_libraries( gen_fish_plex plexrt pprts_base)
install(TARGETS gen_fish_plex DESTINATION bin)

install(TARGETS plexrt DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
