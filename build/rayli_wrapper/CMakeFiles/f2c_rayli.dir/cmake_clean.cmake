file(REMOVE_RECURSE
  "../lib/libf2c_rayli.pdb"
  "../lib/libf2c_rayli.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/f2c_rayli.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
