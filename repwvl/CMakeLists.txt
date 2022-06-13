include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

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


add_executable(compute_repwvl_training_data "compute_repwvl_training_data.F90" )
target_link_libraries( compute_repwvl_training_data repwvl_pprts ${EXTLIB})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/fu.ice.general.nc" ${TENSTREAM_SHARE_DIR}/fu.ice.general.nc   COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/mie.wc.table.nc"   ${TENSTREAM_SHARE_DIR}/mie.wc.table.nc     COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/repwvl_thermal.lut" ${TENSTREAM_SHARE_DIR}/repwvl_thermal.lut COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/repwvl_solar.lut"   ${TENSTREAM_SHARE_DIR}/repwvl_solar.lut   COPYONLY)
