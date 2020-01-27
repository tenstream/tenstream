#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")

file( GLOB files *.f* )

# Dont use Werror for files that are not ours
string(COMPARE EQUAL ${CMAKE_Fortran_COMPILER_ID} "Intel" _is_Intel)
if(${_is_Intel})
  set(WNOERRORFLG "-warn noerrors -warn nounused")
else()
  set(WNOERRORFLG "-w -Wno-error -Wno-pedantic")
endif(${_is_Intel})
file(GLOB rrtm_files *.f*)

foreach(rrtm_file ${rrtm_files})
  SET_SOURCE_FILES_PROPERTIES(${rrtm_file} PROPERTIES COMPILE_FLAGS "${WNOERRORFLG}")
endforeach()

SET_SOURCE_FILES_PROPERTIES(rrtmg_sw_k_g.f90 PROPERTIES COMPILE_FLAGS "-O0 ${WNOERRORFLG}")
SET_SOURCE_FILES_PROPERTIES(mcica_random_numbers.f90 PROPERTIES COMPILE_FLAGS "-fno-range-check ${WNOERRORFLG}")

add_library( pprts_rrtm_sw_objs   OBJECT ${files} )
add_dependencies( pprts_rrtm_sw_objs pprts_objs )
add_library( pprts_rrtm_sw $<TARGET_OBJECTS:pprts_rrtm_sw_objs> )
TARGET_LINK_LIBRARIES( pprts_rrtm_sw pprts)
