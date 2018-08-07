include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( general_sources
  "plex_grid.F90"
  "plex_rt.F90"
  "icon_grid.F90"
  "icon_plex_utils.F90"
  )

add_library( plexrt ${general_sources} )
target_link_libraries( plexrt pprts_base ${LIB} ${EXTLIB} )


add_executable( gen_plex_from_icon gen_plex_from_icon.F90 )
target_link_libraries( gen_plex_from_icon plexrt pprts_base ${LIB} ${EXTLIB} )
install(TARGETS gen_plex_from_icon DESTINATION bin)

file( GLOB plex_examples plex_ex*.F90 )
foreach ( plex_ex ${plex_examples} )
  GET_FILENAME_COMPONENT(basename "${plex_ex}" NAME_WE )

  add_executable("${basename}" "${basename}.F90" )
  target_link_libraries( "${basename}" plexrt pprts_base)

  install(TARGETS ${basename} RUNTIME DESTINATION bin)
endforeach ()

add_executable( gen_fish_plex gen_fish_plex.F90 )
target_link_libraries( gen_fish_plex plexrt pprts_base)
install(TARGETS gen_fish_plex DESTINATION bin)

file( GLOB streams createLUT_*.f90 )
foreach ( stream ${streams} )
  GET_FILENAME_COMPONENT(basename "${stream}" NAME_WE )

  add_executable("${basename}" "${basename}.f90" )
  target_link_libraries( "${basename}" plexrt pprts_base )

  install(TARGETS ${basename} RUNTIME DESTINATION bin)
endforeach ()

install(TARGETS plexrt DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
