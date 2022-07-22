if(NOT ENABLE_PLEXRT)
  return()
endif(NOT ENABLE_PLEXRT)

#get_filename_component(CDIR ${CMAKE_CURRENT_LIST_DIR} NAME)
#message( STATUS "Creating example for directory: CDIR=${CDIR}")
#
#file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.[fF]90")
#
#foreach(_file ${files})
#  get_filename_component (_example_NAME ${_file} NAME   )
#  get_filename_component (_example_WE   ${_file} NAME_WE)
#  get_filename_component (_example_EXT  ${_file} EXT    )
#
#  add_executable("ex_${_example_WE}" "${_example_NAME}")
#  target_link_libraries("ex_${_example_WE}" plexrt_rrtmg pprts_rrtm_lw pprts_rrtm_sw plexrt pprts_base f2c_rayli)
#endforeach()
get_filename_component(CDIR ${CMAKE_CURRENT_LIST_DIR} NAME)
message( STATUS "Creating example for directory: CDIR=${CDIR}")

# While those files that start with "ex_" contain the main() functions and are used to generate executables
file(GLOB exe_files "${CMAKE_CURRENT_SOURCE_DIR}/ex_*.[fF]90")

# We build a lib out of all example source files, this is so that we can use the examples in tests
file(GLOB src_files "${CMAKE_CURRENT_SOURCE_DIR}/*.[fF]90")
list(REMOVE_ITEM src_files ${exe_files})

add_library(plexrt_examples ${src_files})
target_link_libraries(plexrt_examples ${TENSTR_SINGLE_LIB} ${EXTLIB})

foreach(_file ${exe_files})
  get_filename_component (_example_NAME ${_file} NAME   )
  get_filename_component (_example_WE   ${_file} NAME_WE)
  get_filename_component (_example_EXT  ${_file} EXT    )

  add_executable("${_example_WE}" "${_example_NAME}")
  target_link_libraries("${_example_WE}" plexrt_examples)
endforeach()
