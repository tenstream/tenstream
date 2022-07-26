include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( files
  "specint_pprts.F90"
  "specint_plexrt.F90"
  )

add_library( tenstr_specint_objs OBJECT ${files} )
add_dependencies(tenstr_specint_objs
  pprts
  pprts_base
  pprts_rrtmg
  repwvl_pprts
  plexrt
  plexrt_rrtmg
  repwvl_plexrt
  )

add_library( specint ${files} $<TARGET_OBJECTS:tenstr_specint_objs> )

TARGET_LINK_LIBRARIES( specint
  pprts
  pprts_base
  tenstr_disort
  pprts_rrtmg
  repwvl_pprts
  plexrt
  plexrt_rrtmg
  repwvl_plexrt
  ${EXTLIB} ${USERLIB} ${LIBS})
