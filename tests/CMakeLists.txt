include_directories(${PFUNIT_DIR}/mod)
include_directories(${PFUNIT_DIR}/include)

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

# Tests that use PETSc DMDA/PLEX types directly in their source files.
set(PETSC_ONLY_TESTS
  test_petsc_dm_global_to_shared_mem
  test_petsc_scatterToZero
  test_petsc_sort
  test_pprts_rayli
)

# Tests that require PETSc PLEX / unstructured grid solver.
set(PLEXRT_ONLY_TESTS
  plexrt_error_growth_tracking
  test_fish_plex
  test_icon_plex_utils
  test_plex_solution_vecscale
  test_plexrt_fish
  test_plexrt_nca
  test_plexrt_nr_iterations_after_restart
  test_plexrt_rayli
  test_plexrt_rrtmg_lw_sw
)

subdirlist(subdirs ${CMAKE_CURRENT_SOURCE_DIR})

foreach(subdir ${subdirs})
  if(NOT WITH_PETSC AND ${subdir} IN_LIST PETSC_ONLY_TESTS)
    message(STATUS "Skipping test ${subdir} (requires WITH_PETSC=ON)")
  elseif(NOT ENABLE_PLEXRT AND ${subdir} IN_LIST PLEXRT_ONLY_TESTS)
    message(STATUS "Skipping test ${subdir} (requires ENABLE_PLEXRT=ON)")
  else()
    add_subdirectory(${subdir})
  endif()
endforeach()
