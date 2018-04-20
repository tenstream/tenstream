#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file( GLOB files *.f* )

add_library( tenstr_rrtmg STATIC ${files} )
TARGET_LINK_LIBRARIES( tenstr_rrtmg
  tenstr_rrtm_lw tenstr_rrtm_sw tenstream pprts pprts_base)

install(TARGETS tenstr_rrtmg DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
