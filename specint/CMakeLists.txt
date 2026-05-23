include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( files
  "specint_pprts.F90"
  )
if(ENABLE_PLEXRT)
  list(APPEND files "specint_plexrt.F90")
endif()

add_library( tenstr_specint_objs OBJECT ${files} )
add_dependencies(tenstr_specint_objs
  pprts
  pprts_base
  pprts_rrtmg
  repwvl_pprts
  ecckd_pprts
  )
if(ENABLE_PLEXRT)
  add_dependencies(tenstr_specint_objs
    plexrt
    plexrt_rrtmg
    repwvl_plexrt
    )
endif()

add_library( specint $<TARGET_OBJECTS:tenstr_specint_objs> )
if(ENABLE_PLEXRT)
  TARGET_LINK_LIBRARIES( specint
    pprts
    pprts_base
    tenstr_disort
    pprts_rrtmg
    repwvl_pprts
    ecckd_pprts
    plexrt
    plexrt_rrtmg
    repwvl_plexrt
    ${EXTLIB} ${USERLIB} ${LIBS})
else()
  TARGET_LINK_LIBRARIES( specint
    pprts
    pprts_base
    tenstr_disort
    pprts_rrtmg
    repwvl_pprts
    ecckd_pprts
    ${EXTLIB} ${USERLIB} ${LIBS})
endif()
