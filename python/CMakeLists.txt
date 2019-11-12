message(STATUS "Creating python interfaces with ${F2PY_PATH}")

set(DEPENDENCY_LIST tenstream)

list(APPEND wrappers
  "py_boxmc.f90"
  "py_eddington.f90"
  "py_param_phi.f90"
  )

list(APPEND rrtm_wrappers
  "py_rrtm_lw_sw.f90"
  )

set(lib_line "")
foreach ( dep ${DEPENDENCY_LIST} )
  list(APPEND lib_line "-l${dep}")
endforeach()
foreach ( dep ${DEPENDENCY_LIST} )
  list(APPEND lib_line "-l${dep}")
endforeach()
foreach(extlib ${EXTLIB})
  list(APPEND lib_line "${extlib}")
endforeach()
message(STATUS "Linker line: ${lib_line}")

get_property(include_dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
list(APPEND include_dirs "${CMAKE_Fortran_MODULE_DIRECTORY}")
set(inc_line "")
foreach ( inc_dir ${include_dirs} )
  list(APPEND inc_line "-I${inc_dir}")
endforeach()
message(STATUS "Include line: ${inc_line}")

foreach ( wrapper ${wrappers} )
  get_filename_component(basename "${wrapper}" NAME_WE )
  get_filename_component(wrappername "${wrapper}" NAME )

  add_custom_command(OUTPUT ${basename}.so
    COMMAND ${F2PY_PATH}
    --f90exec=${CMAKE_Fortran_COMPILER}
    -m ${basename}
    -c ${CMAKE_CURRENT_SOURCE_DIR}/${wrapper}
    -L${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    ${lib_line}
    ${inc_line}
    -I${CMAKE_Fortran_MODULE_DIRECTORY}
    DEPENDS ${DEPENDENCY_LIST} ${CMAKE_CURRENT_SOURCE_DIR}/${wrapper}
    COMMENT "Building f2py target"
    VERBATIM)

  add_custom_target(do_python_interface_${basename} DEPENDS ${basename}.so)
  add_dependencies(python_interface do_python_interface_${basename})
endforeach ()

#add_executable(test_f2py "py_boxmc.f90" )
#target_link_libraries( test_f2py pprts_base)
