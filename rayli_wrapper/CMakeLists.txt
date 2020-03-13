include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set(depends pprts_base)

if(RAYLI_DIR)
  message(STATUS " Have Rayli and suppose we have RayliLib: ${RAYLI_LIBS}")
  #list(APPEND depends rayli)

  SET_SOURCE_FILES_PROPERTIES(f2c_rayli.F90
    PROPERTIES COMPILE_DEFINITIONS HAVE_RAYLI)
  add_executable(single_wedge_2_rayli "single_wedge_2_rayli.c" )
  target_link_libraries(single_wedge_2_rayli ${RAYLI_LIBS} ${EXTLIB})
  include_directories(${RAYLI_DIR}/include)
endif()

add_library( f2c_rayli_objs OBJECT "f2c_rayli.F90")
add_dependencies( f2c_rayli_objs ${depends})
add_library( f2c_rayli  $<TARGET_OBJECTS:f2c_rayli_objs> )
target_link_libraries(f2c_rayli pprts_base ${RAYLI_LIBS} ${EXTLIB})
