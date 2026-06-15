get_filename_component(CDIR ${CMAKE_CURRENT_LIST_DIR} NAME)
message( STATUS "Creating example for directory: CDIR=${CDIR}")

file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.[fF]90")

foreach(_file ${files})
  get_filename_component (_example_NAME ${_file} NAME   )
  get_filename_component (_example_WE   ${_file} NAME_WE)
  get_filename_component (_example_EXT  ${_file} EXT    )

  add_library("ex_${_example_WE}_objs" OBJECT "${_example_NAME}")
  add_dependencies("ex_${_example_WE}_objs" ${TENSTR_LIBS})

  add_executable("ex_${_example_WE}" EXCLUDE_FROM_ALL $<TARGET_OBJECTS:ex_${_example_WE}_objs>)
  target_link_libraries("ex_${_example_WE}" ${TENSTR_SINGLE_LIB} ${EXTLIB})
  add_dependencies(examples "ex_${_example_WE}")
endforeach()
