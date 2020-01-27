#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")

set( files
  "tenstr_rrtmg_base.F90"
  "dyn_atm_to_rrtmg.F90"
  "optprop_rrtmg.F90"
  )

add_library( tenstr_rrtmg_objs OBJECT ${files} )
add_dependencies(tenstr_rrtmg_objs
  pprts_rrtm_lw
  pprts_rrtm_sw
  plexrt
  pprts
  pprts_base
  tenstr_disort
  )


# PPRTS
add_library( pprts_rrtmg "pprts_rrtmg.F90" $<TARGET_OBJECTS:tenstr_rrtmg_objs> )

TARGET_LINK_LIBRARIES( pprts_rrtmg
  pprts_rrtm_lw
  pprts_rrtm_sw
  pprts
  pprts_base
  tenstr_disort
  ${EXTLIB} ${USERLIB} ${LIBS})


# PlexRT
add_library( plexrt_rrtmg "plexrt_rrtmg.F90" $<TARGET_OBJECTS:tenstr_rrtmg_objs> )

TARGET_LINK_LIBRARIES( plexrt_rrtmg
  pprts_rrtm_lw
  pprts_rrtm_sw
  plexrt
  f2c_twomax
  pprts_base
  tenstr_disort
  ${EXTLIB} ${USERLIB} ${LIBS})
