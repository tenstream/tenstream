include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

file(COPY ./
  DESTINATION ${CMAKE_BINARY_DIR}/include/
  FILES_MATCHING PATTERN "*.h" )

add_library( f2c_twomax twomaxrnd.c f2c_twomax.F90)
add_dependencies( f2c_twomax pprts_base )
target_link_libraries(f2c_twomax pprts_base)
