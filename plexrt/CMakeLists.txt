include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")

set( general_sources
  "plex_rt_base.F90"
  "plex_rt.F90"
  "plex_grid.F90"
  "plexrt_external_solvers.F90"
  "plex2rayli.F90"
  "icon_grid.F90"
  "icon_plex_utils.F90"
  "nca_multi_tri.F90"
  )

add_library( plexrt_objs OBJECT ${general_sources} )
add_dependencies( plexrt_objs
  pprts_base
  tenstr_disort
  f2c_rayli)

add_library( plexrt $<TARGET_OBJECTS:plexrt_objs> )
target_link_libraries(plexrt
  pprts_base
  tenstr_disort
  f2c_rayli
  ${EXTLIB})

# example executables to generate meshes
add_executable( gen_plex_from_icon gen_plex_from_icon.F90 )
target_link_libraries( gen_plex_from_icon plexrt ${EXTLIB})

add_executable( gen_fish_plex gen_fish_plex.F90 )
target_link_libraries( gen_fish_plex plexrt ${EXTLIB})

add_executable( gen_regular_plex gen_regular_plex.F90 )
target_link_libraries( gen_regular_plex plexrt ${EXTLIB})
