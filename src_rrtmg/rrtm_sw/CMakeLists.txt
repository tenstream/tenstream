#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file( GLOB files *.f90 )

SET_SOURCE_FILES_PROPERTIES(rrtmg_sw_k_g.f90 PROPERTIES COMPILE_FLAGS -O0)
SET_SOURCE_FILES_PROPERTIES(mcica_random_numbers.f90 PROPERTIES COMPILE_FLAGS -fno-range-check)

add_library( tenstr_rrtm_sw STATIC ${files} )
TARGET_LINK_LIBRARIES( tenstr_rrtm_sw tenstream)

install(TARGETS tenstr_rrtm_sw   DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
