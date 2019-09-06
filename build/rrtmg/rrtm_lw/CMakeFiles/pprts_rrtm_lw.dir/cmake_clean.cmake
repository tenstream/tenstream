file(REMOVE_RECURSE
  "../../lib/libpprts_rrtm_lw.pdb"
  "../../lib/libpprts_rrtm_lw.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/pprts_rrtm_lw.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
