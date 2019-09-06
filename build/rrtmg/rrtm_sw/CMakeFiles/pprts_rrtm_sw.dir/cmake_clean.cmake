file(REMOVE_RECURSE
  "../../lib/libpprts_rrtm_sw.pdb"
  "../../lib/libpprts_rrtm_sw.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/pprts_rrtm_sw.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
