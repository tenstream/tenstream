list( APPEND test_cmake_SRCS test_cmake.f90 )
#
#list( APPEND solver_sources *.f90 )

file (GLOB_RECURSE Fortran_sourcefiles "*.?90")
#message( "Sources for project ${APP_NAME} are :: ${Fortran_sourcefiles}" )

add_executable( "${APP_NAME}" "${test_cmake_SRCS}" )
