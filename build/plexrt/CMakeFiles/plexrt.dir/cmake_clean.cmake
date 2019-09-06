file(REMOVE_RECURSE
  "../lib/libplexrt.pdb"
  "../lib/libplexrt.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/plexrt.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
