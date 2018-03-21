#message( "Entering tenstream src dir with EXTLIB==${EXTLIB} ")

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file( GLOB files *.f* )

# Dont use Werror for files that are not ours
set(WNOERRORFLG "-w -Wno-error -Wno-pedantic")
file(GLOB rrtm_files *.f* )
file(GLOB tenstr_files tenstr_*.f* )

foreach(tenstr_file ${tenstr_files})
  list(REMOVE_ITEM rrtm_files EXCLUDE ${tenstr_file})
endforeach()
foreach(rrtm_file ${rrtm_files})
  SET_SOURCE_FILES_PROPERTIES(${rrtm_file} PROPERTIES COMPILE_FLAGS "${WNOERRORFLG}")
endforeach()

SET_SOURCE_FILES_PROPERTIES(rrtmg_sw_k_g.f90 PROPERTIES COMPILE_FLAGS "-O0 ${WNOERRORFLG}")
SET_SOURCE_FILES_PROPERTIES(mcica_random_numbers.f90 PROPERTIES COMPILE_FLAGS "-fno-range-check ${WNOERRORFLG}")

add_library( tenstr_rrtm_sw STATIC ${files} )
TARGET_LINK_LIBRARIES( tenstr_rrtm_sw tenstream)

install(TARGETS tenstr_rrtm_sw   DESTINATION lib)
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")
