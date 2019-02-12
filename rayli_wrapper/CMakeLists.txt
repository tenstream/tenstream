set(RAYLI_DIR "$ENV{RAYLI_DIR}")

if(RAYLI_DIR)
  message(STATUS "Configuring with RayliDIR: ${RAYLI_DIR}")
  set(HAVE_RAYLI TRUE)

  include_directories(${RAYLI_DIR}/include)
  find_library(RAYLI_LIB rayli_c_wrapper ${RAYLI_DIR} ${RAYLI_DIR}/lib )
  find_library(RAYLI_LIB2 rayli_common    ${RAYLI_DIR} ${RAYLI_DIR}/lib )
  list(APPEND RAYLI_LIB ${RAYLI_LIB2})
  message(STATUS "Found RAYLI LIB: ${RAYLI_LIB}")

  add_executable(tenstr2rayli "tenstr2rayli.c" )
  target_link_libraries(tenstr2rayli ${RAYLI_LIB} stdc++ m)
endif(RAYLI_DIR)

