cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

enable_language(Fortran)

enable_testing()

macro(subdirlist result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child})
        list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR} )

foreach(subdir ${subdirs})
    add_subdirectory(${subdir})
endforeach()
