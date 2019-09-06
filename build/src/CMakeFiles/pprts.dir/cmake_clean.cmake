file(REMOVE_RECURSE
  "../lib/libpprts.pdb"
  "../lib/libpprts.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/pprts.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
