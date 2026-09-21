# Bundle the objects of several static libraries into one archive.
#
# Invoked as a build-time script:
#   cmake -DAR=<ar> -DOBJ_DIR=<dir> -DOUTPUT=<lib.a> -P combine_archives.cmake
#
# Each input archive has already been unpacked into its own subdirectory of
# OBJ_DIR, because unpacking them all side by side is lossy: several component
# libraries produce objects with identical basenames from *different* sources
# (rrtm_sw/parkind.f90 vs rrtm_lw/parkind.f90, and likewise
# mcica_random_numbers.f90), and a flat unpack silently dropped one of each.
#
# The flip side is that some sources are genuinely compiled into more than one
# component library (repwvl_pprts and repwvl_plexrt both build fu_ice.F90,
# mie_tables.F90, ...), which would then land in the archive twice and give
# duplicate symbols to anyone linking with --whole-archive. Those copies are
# byte-identical, so deduplicate on content hash: same content means the same
# object, different content means two objects that both have to be kept.

file(GLOB_RECURSE _objs "${OBJ_DIR}/*.o")
list(SORT _objs)

set(_keep "")
foreach(_obj ${_objs})
  file(MD5 "${_obj}" _hash)
  if(NOT DEFINED _seen_${_hash})
    set(_seen_${_hash} 1)
    list(APPEND _keep "${_obj}")
  endif()
endforeach()

list(LENGTH _objs _n_all)
list(LENGTH _keep _n_keep)
message(STATUS "combine_archives: ${_n_keep} objects (${_n_all} unpacked, duplicates dropped) -> ${OUTPUT}")

file(REMOVE "${OUTPUT}")
execute_process(
  COMMAND ${AR} -qc "${OUTPUT}" ${_keep}
  RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "combine_archives: ${AR} failed with ${_rc}")
endif()
