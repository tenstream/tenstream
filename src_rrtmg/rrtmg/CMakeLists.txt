#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( files
  "dyn_atm_to_rrtmg.f90"
  "optprop_rrtmg.f90"
  "pprts_rrtmg.f90"
  )

add_library( pprts_rrtmg STATIC ${files} )
TARGET_LINK_LIBRARIES( pprts_rrtmg
  pprts_rrtm_lw pprts_rrtm_sw pprts pprts_base)

install(TARGETS pprts_rrtmg DESTINATION lib)

set( files
  "dyn_atm_to_rrtmg.f90"
  "optprop_rrtmg.f90"
  "plexrt_rrtmg.f90"
  )

add_library( plexrt_rrtmg STATIC ${files} )
TARGET_LINK_LIBRARIES( plexrt_rrtmg
  pprts_rrtm_lw pprts_rrtm_sw pprts_base plexrt)

install(TARGETS plexrt_rrtmg DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
