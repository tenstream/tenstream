include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")
install(DIRECTORY ./ DESTINATION include/ FILES_MATCHING PATTERN "*.h")

set( files
  "mie_tables.F90"
  "fu_ice.F90"
  "rayleigh.F90"
  "repwvl_base.F90"
  "repwvl_pprts.F90"
  )

add_library( tenstr_repwvl_objs OBJECT ${files} )
add_dependencies(tenstr_repwvl_objs
  plexrt
  pprts
  pprts_base
  pprts_rrtmg
  )

# PPRTS
add_library( repwvl_pprts "repwvl_pprts.F90" $<TARGET_OBJECTS:tenstr_repwvl_objs> )

TARGET_LINK_LIBRARIES( repwvl_pprts
  pprts
  pprts_base
  tenstr_disort
  pprts_rrtmg
  ${EXTLIB} ${USERLIB} ${LIBS})


add_executable(compute_training_data "compute_training_data.F90" )
target_link_libraries( compute_training_data repwvl_pprts ${EXTLIB})
