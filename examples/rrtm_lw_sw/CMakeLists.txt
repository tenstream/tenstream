get_filename_component(CDIR ${CMAKE_CURRENT_LIST_DIR} NAME)
message( STATUS "Creating example for directory: CDIR=${CDIR}")

file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/*.[fF]90")
message(STATUS "Found following fortran files: ${files}")

foreach(_file ${files})
  get_filename_component (_example_NAME ${_file} NAME   )
  get_filename_component (_example_WE   ${_file} NAME_WE)
  get_filename_component (_example_EXT  ${_file} EXT    )

  add_executable("ex_${_example_WE}" "${_example_NAME}")
  target_link_libraries("ex_${_example_WE}" tenstr_rrtmg tenstr_rrtm_sw tenstr_rrtm_lw tenstream)
endforeach()
