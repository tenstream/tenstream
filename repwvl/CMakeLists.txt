include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( base_files
  "mie_tables.F90"
  "fu_ice.F90"
  "rayleigh.F90"
  "repwvl_base.F90"
  "repwvl_optprop.F90"
  )

add_library( tenstr_repwvl_objs OBJECT ${base_files} )
add_dependencies(tenstr_repwvl_objs
  pprts_base
  pprts_rrtmg
  )
if(ENABLE_PLEXRT)
  add_dependencies(tenstr_repwvl_objs
    plexrt
    pprts
    )
endif()

# PPRTS
add_library( repwvl_pprts_objs OBJECT "repwvl_pprts.F90" )
add_dependencies( repwvl_pprts_objs tenstr_repwvl_objs pprts_base pprts_rrtmg )

add_library( repwvl_pprts $<TARGET_OBJECTS:repwvl_pprts_objs> $<TARGET_OBJECTS:tenstr_repwvl_objs> )
TARGET_LINK_LIBRARIES( repwvl_pprts
  pprts_base
  tenstr_disort
  pprts_rrtmg
  ${EXTLIB} ${USERLIB} ${LIBS})

if(ENABLE_PLEXRT)
  # PlexRT
  add_library( repwvl_plexrt_objs OBJECT "repwvl_plexrt.F90" )
  add_dependencies( repwvl_plexrt_objs tenstr_repwvl_objs plexrt plexrt_rrtmg )

  add_library( repwvl_plexrt $<TARGET_OBJECTS:repwvl_plexrt_objs> $<TARGET_OBJECTS:tenstr_repwvl_objs> )
  TARGET_LINK_LIBRARIES( repwvl_plexrt
    plexrt
    plexrt_rrtmg
    ${EXTLIB} ${USERLIB} ${LIBS})
endif()

if(WITH_PETSC)
  add_executable(compute_repwvl_training_data "compute_repwvl_training_data.F90" )
  target_link_libraries( compute_repwvl_training_data ${TENSTR_SINGLE_LIB} ${EXTLIB})
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/fu.ice.general.nc" ${TENSTREAM_SHARE_DIR}/fu.ice.general.nc   COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/mie.wc.table.nc"   ${TENSTREAM_SHARE_DIR}/mie.wc.table.nc     COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/repwvl_thermal.lut" ${TENSTREAM_SHARE_DIR}/repwvl_thermal.lut COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/repwvl_solar.lut"   ${TENSTREAM_SHARE_DIR}/repwvl_solar.lut   COPYONLY)
