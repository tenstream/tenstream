#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file( GLOB files *.f* )

add_library( pprts_rrtmg STATIC ${files} )
TARGET_LINK_LIBRARIES( pprts_rrtmg
  pprts_rrtm_lw pprts_rrtm_sw plexrt pprts pprts_base)

install(TARGETS pprts_rrtmg DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
