file(REMOVE_RECURSE
  "../lib/libpprts_base.pdb"
  "../lib/libpprts_base.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang C Fortran)
  include(CMakeFiles/pprts_base.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
