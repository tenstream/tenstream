include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

set( files
  "ecckd_base.F90"
  "ecckd_pprts.F90"
  "ecckd_optprop.F90"
  )

add_library( tenstr_ecckd_objs OBJECT ${files} )
add_dependencies(tenstr_ecckd_objs
  plexrt
  pprts
  pprts_base
  pprts_rrtmg
  repwvl_pprts
  )

# PPRTS
add_library( ecckd_pprts "ecckd_pprts.F90" $<TARGET_OBJECTS:tenstr_ecckd_objs> )

TARGET_LINK_LIBRARIES( ecckd_pprts
  pprts
  pprts_base
  tenstr_disort
  pprts_rrtmg
  repwvl_pprts
  ${EXTLIB} ${USERLIB} ${LIBS})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/ecckd-1.0_lw_climate_fsck-16_ckd-definition.nc" ${TENSTREAM_SHARE_DIR}/ecckd-1.0_lw_climate_fsck-16_ckd-definition.nc COPYONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/data/ecckd-1.0_sw_climate_rgb-16_ckd-definition.nc"  ${TENSTREAM_SHARE_DIR}/ecckd-1.0_sw_climate_rgb-16_ckd-definition.nc  COPYONLY)
