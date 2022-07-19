include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( files
  "specint_pprts.F90"
  )

add_library( tenstr_specint_objs OBJECT ${files} )
add_dependencies(tenstr_specint_objs
  pprts
  pprts_base
  pprts_rrtmg
  repwvl_pprts
  )

add_library( specint_pprts "specint_pprts.F90" $<TARGET_OBJECTS:tenstr_specint_objs> )

TARGET_LINK_LIBRARIES( specint_pprts
  pprts
  pprts_base
  tenstr_disort
  pprts_rrtmg
  repwvl_pprts
  ${EXTLIB} ${USERLIB} ${LIBS})
