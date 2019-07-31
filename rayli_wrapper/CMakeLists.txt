set(RAYLI_DIR "$ENV{RAYLI_DIR}")

if(RAYLI_DIR)
  message(STATUS "Configuring with RayliDIR: ${RAYLI_DIR}")
  set(HAVE_RAYLI TRUE)

  include_directories(${RAYLI_DIR}/include)
  find_library(RAYLI_LIB rayli_c_wrapper ${RAYLI_DIR} ${RAYLI_DIR}/lib )
  find_library(RAYLI_LIB2 rayli_common    ${RAYLI_DIR} ${RAYLI_DIR}/lib )
  list(APPEND RAYLI_LIB ${RAYLI_LIB2} stdc++ pthread m)
  message(STATUS "Found RAYLI LIB: ${RAYLI_LIB}")
endif(RAYLI_DIR)

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

if(HAVE_RAYLI)
  SET_SOURCE_FILES_PROPERTIES(f2c_rayli.F90
    PROPERTIES COMPILE_DEFINITIONS __HAVE_RAYLI__)

  add_executable(single_wedge_2_rayli "single_wedge_2_rayli.c" )
  target_link_libraries(single_wedge_2_rayli ${RAYLI_LIB})
endif()

add_library( f2c_rayli_objs   OBJECT "f2c_rayli.F90")
add_dependencies( f2c_rayli_objs pprts_base )
add_library( f2c_rayli  $<TARGET_OBJECTS:f2c_rayli_objs> )
target_link_libraries(f2c_rayli pprts_base ${RAYLI_LIB})
