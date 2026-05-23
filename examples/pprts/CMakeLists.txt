get_filename_component(CDIR ${CMAKE_CURRENT_LIST_DIR} NAME)
message( STATUS "Creating example for directory: CDIR=${CDIR}")

# While those files that start with "ex_" contain the main() functions and are used to generate executables
file(GLOB exe_files "${CMAKE_CURRENT_SOURCE_DIR}/ex_*.[fF]90")

# We build a lib out of all example source files, this is so that we can use the examples in tests
file(GLOB src_files "${CMAKE_CURRENT_SOURCE_DIR}/*.[fF]90")
list(REMOVE_ITEM src_files ${exe_files})

add_library(pprts_examples_objs OBJECT ${src_files})
add_dependencies(pprts_examples_objs pprts pprts_base)

add_library(pprts_examples $<TARGET_OBJECTS:pprts_examples_objs>)
target_link_libraries(pprts_examples ${TENSTR_SINGLE_LIB} ${EXTLIB})

foreach(_file ${exe_files})
  get_filename_component (_example_NAME ${_file} NAME   )
  get_filename_component (_example_WE   ${_file} NAME_WE)
  get_filename_component (_example_EXT  ${_file} EXT    )

  add_library("${_example_WE}_objs" OBJECT "${_example_NAME}")
  add_dependencies("${_example_WE}_objs" pprts_examples_objs)

  add_executable("${_example_WE}" EXCLUDE_FROM_ALL $<TARGET_OBJECTS:${_example_WE}_objs>)
  target_link_libraries("${_example_WE}" pprts_examples)
endforeach()
