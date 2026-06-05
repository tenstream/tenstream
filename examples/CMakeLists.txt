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

# These example directories use PETSc DMDA/PLEX types directly in their
# public API and cannot compile without PETSc headers.
set(PETSC_ONLY_EXAMPLES
  geometric_coeffs_dir2diff_hacks
  geometric_coeffs_num_intervals
  plexrt
  wetterstein_ts
)

add_custom_target(examples)

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR})

foreach(subdir ${subdirs})
  if(NOT WITH_PETSC AND ${subdir} IN_LIST PETSC_ONLY_EXAMPLES)
    message(STATUS "Skipping example ${subdir} (requires WITH_PETSC=ON)")
  else()
    add_subdirectory(${subdir})
  endif()
endforeach()
