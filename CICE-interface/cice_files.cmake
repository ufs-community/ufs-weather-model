list(APPEND cice_shared_files
  #Shared List:
  CICE/cicecore/shared/ice_arrays_column.F90
  CICE/cicecore/shared/ice_calendar.F90
  CICE/cicecore/shared/ice_constants.F90
  CICE/cicecore/shared/ice_distribution.F90
  CICE/cicecore/shared/ice_domain_size.F90
  CICE/cicecore/shared/ice_fileunits.F90
  CICE/cicecore/shared/ice_init_column.F90
  CICE/cicecore/shared/ice_kinds_mod.F90
  CICE/cicecore/shared/ice_restart_column.F90
  CICE/cicecore/shared/ice_restart_shared.F90
  CICE/cicecore/shared/ice_spacecurve.F90

  #Analysis List:
  CICE/cicecore/cicedyn/analysis/ice_diagnostics.F90
  CICE/cicecore/cicedyn/analysis/ice_diagnostics_bgc.F90
  CICE/cicecore/cicedyn/analysis/ice_history.F90
  CICE/cicecore/cicedyn/analysis/ice_history_bgc.F90
  CICE/cicecore/cicedyn/analysis/ice_history_drag.F90
  CICE/cicecore/cicedyn/analysis/ice_history_fsd.F90
  CICE/cicecore/cicedyn/analysis/ice_history_mechred.F90
  CICE/cicecore/cicedyn/analysis/ice_history_pond.F90
  CICE/cicecore/cicedyn/analysis/ice_history_shared.F90
  CICE/cicecore/cicedyn/analysis/ice_history_snow.F90

  #Dynamics List:
  CICE/cicecore/cicedyn/dynamics/ice_dyn_eap.F90
  CICE/cicecore/cicedyn/dynamics/ice_dyn_evp.F90
  CICE/cicecore/cicedyn/dynamics/ice_dyn_evp_1d.F90
  CICE/cicecore/cicedyn/dynamics/ice_dyn_shared.F90
  CICE/cicecore/cicedyn/dynamics/ice_dyn_vp.F90
  CICE/cicecore/cicedyn/dynamics/ice_transport_driver.F90
  CICE/cicecore/cicedyn/dynamics/ice_transport_remap.F90

  #General List:
  CICE/cicecore/cicedyn/general/ice_flux.F90
  CICE/cicecore/cicedyn/general/ice_flux_bgc.F90
  CICE/cicecore/cicedyn/general/ice_forcing.F90
  CICE/cicecore/cicedyn/general/ice_forcing_bgc.F90
  CICE/cicecore/cicedyn/general/ice_init.F90
  CICE/cicecore/cicedyn/general/ice_state.F90
  CICE/cicecore/cicedyn/general/ice_step_mod.F90

  #Infrastructure List
  CICE/cicecore/cicedyn/infrastructure/ice_blocks.F90
  CICE/cicecore/cicedyn/infrastructure/ice_domain.F90
  CICE/cicecore/cicedyn/infrastructure/ice_grid.F90
  CICE/cicecore/cicedyn/infrastructure/ice_read_write.F90
  CICE/cicecore/cicedyn/infrastructure/ice_restart_driver.F90
  CICE/cicecore/cicedyn/infrastructure/ice_restoring.F90
  CICE/cicecore/cicedyn/infrastructure/ice_memusage.F90
  CICE/cicecore/cicedyn/infrastructure/ice_memusage_gptl.c)


#Icepack List:
list(APPEND icepack_files
  CICE/icepack/columnphysics/icepack_aerosol.F90
  CICE/icepack/columnphysics/icepack_age.F90
  CICE/icepack/columnphysics/icepack_algae.F90
  CICE/icepack/columnphysics/icepack_atmo.F90
  CICE/icepack/columnphysics/icepack_brine.F90
  CICE/icepack/columnphysics/icepack_firstyear.F90
  CICE/icepack/columnphysics/icepack_flux.F90
  CICE/icepack/columnphysics/icepack_fsd.F90
  CICE/icepack/columnphysics/icepack_intfc.F90
  CICE/icepack/columnphysics/icepack_isotope.F90
  CICE/icepack/columnphysics/icepack_itd.F90
  CICE/icepack/columnphysics/icepack_kinds.F90
  CICE/icepack/columnphysics/icepack_mechred.F90
  CICE/icepack/columnphysics/icepack_meltpond_lvl.F90
  CICE/icepack/columnphysics/icepack_meltpond_topo.F90
  CICE/icepack/columnphysics/icepack_mushy_physics.F90
  CICE/icepack/columnphysics/icepack_ocean.F90
  CICE/icepack/columnphysics/icepack_orbital.F90
  CICE/icepack/columnphysics/icepack_parameters.F90
  CICE/icepack/columnphysics/icepack_shortwave.F90
  CICE/icepack/columnphysics/icepack_snow.F90
  CICE/icepack/columnphysics/icepack_therm_bl99.F90
  CICE/icepack/columnphysics/icepack_therm_itd.F90
  CICE/icepack/columnphysics/icepack_therm_mushy.F90
  CICE/icepack/columnphysics/icepack_therm_shared.F90
  CICE/icepack/columnphysics/icepack_therm_vertical.F90
  CICE/icepack/columnphysics/icepack_tracers.F90
  CICE/icepack/columnphysics/icepack_warnings.F90
  CICE/icepack/columnphysics/icepack_wavefracspec.F90
  CICE/icepack/columnphysics/icepack_zbgc.F90
  CICE/icepack/columnphysics/icepack_zbgc_shared.F90
  CICE/icepack/columnphysics/icepack_zsalinity.F90
)

list(APPEND cice_shared_files_c
  CICE/cicecore/cicedyn/infrastructure/ice_shr_reprosum86.c
)

#-- Using MPI
list(APPEND cice_mpi_comm_files
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_boundary.F90
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_broadcast.F90
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_communicate.F90
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_exit.F90
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_gather_scatter.F90
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_global_reductions.F90
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_reprosum.F90
  CICE/cicecore/cicedyn/infrastructure/comm/mpi/ice_timers.F90
)

#-- Using Serial
list(APPEND cice_serial_comm_files
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_boundary.F90
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_broadcast.F90
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_communicate.F90
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_exit.F90
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_gather_scatter.F90
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_global_reductions.F90
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_reprosum.F90
  CICE/cicecore/cicedyn/infrastructure/comm/serial/ice_timers.F90
)

#-- Using binary IO
list(APPEND cice_binary_io_files
  CICE/cicecore/cicedyn/infrastructure/io/io_binary/ice_history_write.F90
  CICE/cicecore/cicedyn/infrastructure/io/io_binary/ice_restart.F90
)

#-- Using NetCDF IO
list(APPEND cice_netcdf_io_files
  CICE/cicecore/cicedyn/infrastructure/io/io_netcdf/ice_history_write.F90
  CICE/cicecore/cicedyn/infrastructure/io/io_netcdf/ice_restart.F90
)

#PIO2 I/O List:
list(APPEND cice_pio2_io_files
  CICE/cicecore/cicedyn/infrastructure/io/io_pio2/ice_history_write.F90
  CICE/cicecore/cicedyn/infrastructure/io/io_pio2/ice_pio.F90
  CICE/cicecore/cicedyn/infrastructure/io/io_pio2/ice_restart.F90
)

#-- Using NUOPC CMEPS driver
list(APPEND cice_nuopc_cmeps_driver_files
  CICE/cicecore/drivers/nuopc/cmeps/CICE_FinalMod.F90
  CICE/cicecore/drivers/nuopc/cmeps/CICE_InitMod.F90
  CICE/cicecore/drivers/nuopc/cmeps/CICE_RunMod.F90
  CICE/cicecore/drivers/nuopc/cmeps/cice_wrapper_mod.F90
  CICE/cicecore/drivers/nuopc/cmeps/ice_comp_nuopc.F90
  CICE/cicecore/drivers/nuopc/cmeps/ice_import_export.F90
  CICE/cicecore/drivers/nuopc/cmeps/ice_prescribed_mod.F90
  CICE/cicecore/drivers/nuopc/cmeps/ice_scam.F90
  CICE/cicecore/drivers/nuopc/cmeps/ice_shr_methods.F90
  CICE/cicecore/drivers/nuopc/cmeps/ice_mesh_mod.F90
)
