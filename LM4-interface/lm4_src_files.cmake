list(APPEND lm4_src_files
     LM4/land_data.F90
     LM4/transitions/transitions.F90
     LM4/snow/snow_tile.F90
     LM4/snow/snow.F90
     LM4/canopy_air/canopy_air.F90
     LM4/canopy_air/cana_tile.F90
     LM4/lake/lake_tile.F90
     LM4/lake/lake.F90
     LM4/river/river_physics.F90
     LM4/river/river.F90
     LM4/river/river_type.F90
     LM4/vegetation/vegn_radiation.F90
     LM4/vegetation/vegn_tile.F90
     LM4/vegetation/vegetation.F90
     LM4/vegetation/vegn_cohort.F90
     LM4/vegetation/read_remap_cohort_data_new.inc
     LM4/vegetation/vegn_static_override.F90
     LM4/vegetation/vegn_cohort_io.F90
     LM4/vegetation/vegn_harvesting.F90
     LM4/vegetation/vegn_data.F90
     LM4/vegetation/vegn_photosynthesis.F90
     LM4/vegetation/vegn_disturbance.F90
     LM4/vegetation/vegn_dynamics.F90
     LM4/shared/land_utils.F90
     LM4/shared/land_tile_io.F90
     LM4/shared/land_tile_diag_sel.F90
     LM4/shared/sphum.F90
     LM4/shared/land_tile_diag_buff.F90
     LM4/shared/table_printer.F90
     LM4/shared/land_io.F90
     LM4/shared/land_debug.F90
     LM4/shared/land_tile_diag.F90
     LM4/shared/debug.inc
     LM4/shared/version_variable.inc
     LM4/shared/land_numerics.F90
     LM4/shared/sat_vapor_pres/sat_vapor_pres.F90
     LM4/shared/sat_vapor_pres/sat_vapor_pres_k.F90		
     LM4/land_tracers/land_tracer_driver.F90
     LM4/land_tracers/land_dust.F90
     LM4/land_tracers/land_tracers.F90
     LM4/predefined_tiles/tiling_input.F90
     LM4/predefined_tiles/tiling_input_types.F90
     LM4/land_constants.F90
     LM4/land_model.F90
     LM4/soil/soil.F90
     LM4/soil/hillslope_tile.F90
     LM4/soil/soil_carbon.F90
     LM4/soil/hillslope.F90
     LM4/soil/uptake.F90
     LM4/soil/hillslope_hydrology.F90
     LM4/soil/soil_tile.F90
     LM4/land_tile.F90
     LM4/land_version.inc
     LM4/topo_rough/topo_rough.F90
     LM4/glacier/glac_tile.F90
     LM4/glacier/glacier.F90
     LM4/land_chksum.F90

     LM4/nuopc_cap/lm4_cap.F90
     LM4/nuopc_cap/lm4_driver.F90
     LM4/nuopc_cap/lm4_type.F90
     LM4/nuopc_cap/lnd_import_export.F90
     LM4/nuopc_cap/nuopc_shr_methods.F90
     LM4/nuopc_cap/shr_file_mod.F90
     LM4/nuopc_cap/shr_kind_mod.F90
     LM4/nuopc_cap/shr_log_mod.F90
     LM4/nuopc_cap/shr_sys_mod.F90
     LM4/nuopc_cap/domain_create.F90
     LM4/nuopc_cap/proc_bounds.F90

    )

# TODO: seperate cap files into seperate list
# list(APPEND lm4_nuopc_src_files
#      nuopc_cap/lm4_cap.F90
# )