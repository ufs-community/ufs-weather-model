#!/usr/bin/env python

# CCPP prebuild config for GFDL Finite-Volume Cubed-Sphere Model (FV3) v1


###############################################################################
# Definitions                                                                 #
###############################################################################

HOST_MODEL_IDENTIFIER = "FV3"

# Add all files with metadata tables on the host model side,
# relative to basedir = top-level directory of host model
VARIABLE_DEFINITION_FILES = [
    'ccpp/physics/physics/machine.F',
    'FV3/gfsphysics/CCPP_layer/CCPP_typedefs.F90',
    'FV3/gfsphysics/GFS_layer/GFS_typedefs.F90',
    'FV3/gfsphysics/CCPP_layer/CCPP_data.F90',
    ]

# Add all physics scheme dependencies relative to basedir - note that the CCPP
# rules stipulate that dependencies are not shared between the schemes!
SCHEME_FILES_DEPENDENCIES = [
    'ccpp/physics/physics/GFDL_parse_tracers.F90',
    'ccpp/physics/physics/aer_cloud.F',
    'ccpp/physics/physics/aerclm_def.F',
    'ccpp/physics/physics/aerinterp.F90',
    'ccpp/physics/physics/calpreciptype.f90',
    'ccpp/physics/physics/cldwat2m_micro.F',
    'ccpp/physics/physics/cldmacro.F',
    'ccpp/physics/physics/date_def.f',
    'ccpp/physics/physics/funcphys.f90',
    'ccpp/physics/physics/gcycle.F90',
    'ccpp/physics/physics/gfs_phy_tracer_config.F',
    'ccpp/physics/physics/gocart_tracer_config_stub.f',
    'ccpp/physics/physics/h2o_def.f',
    'ccpp/physics/physics/h2ointerp.f90',
    'ccpp/physics/physics/iccn_def.F',
    'ccpp/physics/physics/iccninterp.F90',
    'ccpp/physics/physics/iounitdef.f',
    'ccpp/physics/physics/machine.F',
    'ccpp/physics/physics/mersenne_twister.f',
    'ccpp/physics/physics/mfpbl.f',
    'ccpp/physics/physics/micro_mg_utils.F90',
    'ccpp/physics/physics/micro_mg2_0.F90',
    'ccpp/physics/physics/micro_mg3_0.F90',
    'ccpp/physics/physics/module_bfmicrophysics.f',
    'ccpp/physics/physics/multi_gases.F90',
    'ccpp/physics/physics/module_gfdl_cloud_microphys.F90',
    'ccpp/physics/physics/module_nst_model.f90',
    'ccpp/physics/physics/module_nst_parameters.f90',
    'ccpp/physics/physics/module_nst_water_prop.f90',
    'ccpp/physics/physics/module_mp_radar.F90',
    'ccpp/physics/physics/module_mp_thompson.F90',
    'ccpp/physics/physics/module_mp_thompson_make_number_concentrations.F90',
    'ccpp/physics/physics/module_bl_mynn.F90',
    'ccpp/physics/physics/module_sf_mynn.F90',
    'ccpp/physics/physics/namelist_soilveg.f',
    'ccpp/physics/physics/mfpblt.f',
    'ccpp/physics/physics/mfscu.f',
    'ccpp/physics/physics/num_parthds.F',
    'ccpp/physics/physics/ozne_def.f',
    'ccpp/physics/physics/ozinterp.f90',
    'ccpp/physics/physics/physcons.F90',
    'ccpp/physics/physics/physparam.f',
    'ccpp/physics/physics/radcons.f90',
    'ccpp/physics/physics/radiation_aerosols.f',
    'ccpp/physics/physics/radiation_astronomy.f',
    'ccpp/physics/physics/radiation_clouds.f',
    'ccpp/physics/physics/radiation_gases.f',
    'ccpp/physics/physics/radiation_surface.f',
    'ccpp/physics/physics/radlw_datatb.f',
    'ccpp/physics/physics/radlw_param.f',
    'ccpp/physics/physics/radsw_datatb.f',
    'ccpp/physics/physics/radsw_param.f',
    'ccpp/physics/physics/sfcsub.F',
    'ccpp/physics/physics/sflx.f',
    'ccpp/physics/physics/set_soilveg.f',
    'ccpp/physics/physics/surface_perturbation.F90',
    'ccpp/physics/physics/cu_gf_deep.F90',
    'ccpp/physics/physics/cu_gf_sh.F90',
    'ccpp/physics/physics/tridi.f',
    'ccpp/physics/physics/wv_saturation.F',
    'ccpp/physics/physics/module_sf_ruclsm.F90',
    'ccpp/physics/physics/namelist_soilveg_ruc.F90',
    'ccpp/physics/physics/set_soilveg_ruc.F90',
    'ccpp/physics/physics/module_soil_pre.F90',
    # stochastic physics
    'ccpp/physics/stochastic_physics/spectral_layout.f',
    'ccpp/physics/stochastic_physics/stochy_gg_def.f',
    'ccpp/physics/stochastic_physics/stochy_resol_def.f',
    'ccpp/physics/stochastic_physics/stochy_layout_lag.f',
    'ccpp/physics/stochastic_physics/four_to_grid_stochy.F',
    'ccpp/physics/stochastic_physics/glats_stochy.f',
    'ccpp/physics/stochastic_physics/sumfln_stochy.f',
    'ccpp/physics/stochastic_physics/gozrineo_stochy.f',
    'ccpp/physics/stochastic_physics/get_ls_node_stochy.f',
    'ccpp/physics/stochastic_physics/get_lats_node_a_stochy.f',
    'ccpp/physics/stochastic_physics/setlats_a_stochy.f',
    'ccpp/physics/stochastic_physics/setlats_lag_stochy.f',
    'ccpp/physics/stochastic_physics/epslon_stochy.f',
    'ccpp/physics/stochastic_physics/getcon_lag_stochy.f',
    'ccpp/physics/stochastic_physics/pln2eo_stochy.f',
    'ccpp/physics/stochastic_physics/dozeuv_stochy.f',
    'ccpp/physics/stochastic_physics/dezouv_stochy.f',
    'ccpp/physics/stochastic_physics/getcon_spectral.F90',
    'ccpp/physics/stochastic_physics/stochy_namelist_def.F90',
    'ccpp/physics/stochastic_physics/compns_stochy.F90',
    'ccpp/physics/stochastic_physics/stochy_internal_state_mod.F90',
    'ccpp/physics/stochastic_physics/stochy_patterngenerator.F90',
    'ccpp/physics/stochastic_physics/stochy_data_mod.F90',
    'ccpp/physics/stochastic_physics/get_stochy_pattern.F90',
    'ccpp/physics/stochastic_physics/initialize_spectral_mod.F90',
    'ccpp/physics/stochastic_physics/stochy_ccpp.F90',
    # derived data type definitions
    'FV3/gfsphysics/GFS_layer/GFS_typedefs.F90',
    'FV3/gfsphysics/CCPP_layer/CCPP_typedefs.F90',
    ]

# Add all physics scheme files relative to basedir
SCHEME_FILES = {
    # Relative path to source (from where ccpp_prebuild.py is called) : [ list of physics sets in which scheme may be called ];
    # current restrictions are that each scheme can only belong to one physics set, and all schemes within one group in the
    # suite definition file have to belong to the same physics set
    'ccpp/physics/physics/GFS_DCNV_generic.F90'              : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_MP_generic.F90'                : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_PBL_generic.F90'               : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_SCNV_generic.F90'              : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_debug.F90'                     : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_phys_time_vary.fv3.F90'        : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_rad_time_vary.fv3.F90'         : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_rrtmg_post.F90'                : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_rrtmg_pre.F90'                 : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_rrtmg_setup.F90'               : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_stochastics.F90'               : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_suite_interstitial.F90'        : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_surface_generic.F90'           : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_surface_composites.F90'        : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_surface_loop_control.F90'      : [ 'slow_physics' ],
    'ccpp/physics/physics/GFS_time_vary_pre.fv3.F90'         : [ 'slow_physics' ],
    'ccpp/physics/physics/cnvc90.f'                          : [ 'slow_physics' ],
    'ccpp/physics/physics/cs_conv.F90'                       : [ 'slow_physics' ],
    'ccpp/physics/physics/cs_conv_aw_adj.F90'                : [ 'slow_physics' ],
    'ccpp/physics/physics/cu_ntiedtke_pre.F90'               : [ 'slow_physics' ],
    'ccpp/physics/physics/cu_ntiedtke.F90'                   : [ 'slow_physics' ],
    'ccpp/physics/physics/cu_ntiedtke_post.F90'              : [ 'slow_physics' ],
    'ccpp/physics/physics/dcyc2.f'                           : [ 'slow_physics' ],
    'ccpp/physics/physics/gcm_shoc.F90'                      : [ 'slow_physics' ],
    'ccpp/physics/physics/get_prs_fv3.F90'                   : [ 'slow_physics' ],
    'ccpp/physics/physics/gfdl_cloud_microphys.F90'          : [ 'slow_physics' ],
    'ccpp/physics/physics/gfdl_fv_sat_adj.F90'               : [ 'fast_physics' ],
    'ccpp/physics/physics/gscond.f'                          : [ 'slow_physics' ],
    'ccpp/physics/physics/gwdc.f'                            : [ 'slow_physics' ],
    'ccpp/physics/physics/gwdps.f'                           : [ 'slow_physics' ],
    'ccpp/physics/physics/h2ophys.f'                         : [ 'slow_physics' ],
    'ccpp/physics/physics/samfdeepcnv.f'                     : [ 'slow_physics' ],
    'ccpp/physics/physics/samfshalcnv.f'                     : [ 'slow_physics' ],
    'ccpp/physics/physics/maximum_hourly_diagnostics.F90'    : [ 'slow_physics' ],
    'ccpp/physics/physics/m_micro.F90'                       : [ 'slow_physics' ],
    'ccpp/physics/physics/m_micro_interstitial.F90'          : [ 'slow_physics' ],
    'ccpp/physics/physics/cu_gf_driver_pre.F90'              : [ 'slow_physics' ],
    'ccpp/physics/physics/cu_gf_driver.F90'                  : [ 'slow_physics' ],
    'ccpp/physics/physics/cu_gf_driver_post.F90'             : [ 'slow_physics' ],
    'ccpp/physics/physics/moninedmf.f'                       : [ 'slow_physics' ],
    'ccpp/physics/physics/moninshoc.f'                       : [ 'slow_physics' ],
    'ccpp/physics/physics/satmedmfvdif.F'                    : [ 'slow_physics' ],
    'ccpp/physics/physics/shinhongvdif.F90'                  : [ 'slow_physics' ],
    'ccpp/physics/physics/ysuvdif.F90'                       : [ 'slow_physics' ],
    'ccpp/physics/physics/module_MYNNPBL_wrapper.F90'        : [ 'slow_physics' ],
    'ccpp/physics/physics/module_MYNNSFC_wrapper.F90'        : [ 'slow_physics' ],
    'ccpp/physics/physics/module_MYNNrad_pre.F90'            : [ 'slow_physics' ],
    'ccpp/physics/physics/module_MYNNrad_post.F90'           : [ 'slow_physics' ],
    'ccpp/physics/physics/mp_thompson_pre.F90'               : [ 'slow_physics' ],
    'ccpp/physics/physics/mp_thompson.F90'                   : [ 'slow_physics' ],
    'ccpp/physics/physics/mp_thompson_post.F90'              : [ 'slow_physics' ],
    'ccpp/physics/physics/ozphys.f'                          : [ 'slow_physics' ],
    'ccpp/physics/physics/ozphys_2015.f'                     : [ 'slow_physics' ],
    'ccpp/physics/physics/precpd.f'                          : [ 'slow_physics' ],
    'ccpp/physics/physics/radlw_main.f'                      : [ 'slow_physics' ],
    'ccpp/physics/physics/radsw_main.f'                      : [ 'slow_physics' ],
    'ccpp/physics/physics/rayleigh_damp.f'                   : [ 'slow_physics' ],
    'ccpp/physics/physics/rrtmg_lw_post.F90'                 : [ 'slow_physics' ],
    'ccpp/physics/physics/rrtmg_lw_pre.F90'                  : [ 'slow_physics' ],
    'ccpp/physics/physics/rrtmg_sw_post.F90'                 : [ 'slow_physics' ],
    'ccpp/physics/physics/rrtmg_sw_pre.F90'                  : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_diag.f'                        : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_diag_post.F90'                 : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_drv_ruc.F90'                   : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_diff.f'                        : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_drv.f'                         : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_nst.f'                         : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_ocean.F'                       : [ 'slow_physics' ],
    'ccpp/physics/physics/sfc_sice.f'                        : [ 'slow_physics' ],
    # stochastic physics
    'ccpp/physics/stochastic_physics/stochastic_physics.F90' : [ 'slow_physics' ],
    # for testing the <init> and <finalize> sections
    'ccpp/physics/physics/GFS_suite_init_finalize_test.F90'  : [ 'slow_physics' ],
    }

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE = 'ccpp/physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE = 'ccpp/physics/CCPP_SCHEMES.cmake'

# CCPP host cap in which to insert the ccpp_field_add statements;
# determines the directory to place ccpp_{modules,fields}.inc
TARGET_FILES = [
    'FV3/atmos_cubed_sphere/driver/fvGFS/atmosphere.F90',
    'FV3/CCPP_layer/CCPP_Driver.F90',
    ]

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE = 'ccpp/physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE = 'ccpp/physics/CCPP_CAPS.cmake'

# Directory where to put all auto-generated physics caps
CAPS_DIR = 'ccpp/physics/physics'

# Directory where the suite definition files are stored
SUITES_DIR = 'ccpp/suites'

# Optional arguments - only required for schemes that use
# optional arguments. ccpp_prebuild.py will throw an exception
# if it encounters a scheme subroutine with optional arguments
# if no entry is made here. Possible values are: 'all', 'none',
# or a list of standard_names: [ 'var1', 'var3' ].
OPTIONAL_ARGUMENTS = {
    'rrtmg_sw' : {
        'rrtmg_sw_run' : [
            'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step',
            'components_of_surface_downward_shortwave_fluxes',
            'cloud_liquid_water_path',
            'mean_effective_radius_for_liquid_cloud',
            'cloud_ice_water_path',
            'mean_effective_radius_for_ice_cloud',
            'cloud_rain_water_path',
            'mean_effective_radius_for_rain_drop',
            'cloud_snow_water_path',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'rrtmg_lw' : {
        'rrtmg_lw_run' : [
            'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step',
            'cloud_liquid_water_path',
            'mean_effective_radius_for_liquid_cloud',
            'cloud_ice_water_path',
            'mean_effective_radius_for_ice_cloud',
            'cloud_rain_water_path',
            'mean_effective_radius_for_rain_drop',
            'cloud_snow_water_path',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'mp_thompson' : {
        'mp_thompson_init' : [
            'water_friendly_aerosol_number_concentration',
            'ice_friendly_aerosol_number_concentration',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            ],
        'mp_thompson_run' : [
            'cloud_droplet_number_concentration_updated_by_physics',
            'water_friendly_aerosol_number_concentration_updated_by_physics',
            'ice_friendly_aerosol_number_concentration_updated_by_physics',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            'mean_effective_radius_for_liquid_cloud',
            'mean_effective_radius_for_ice_cloud',
            'mean_effective_radius_for_snow_flake',
            ],
        },
    'mp_thompson_pre' : {
        'mp_thompson_pre_run' : [
            'cloud_droplet_number_concentration_updated_by_physics',
            'water_friendly_aerosol_number_concentration_updated_by_physics',
            'ice_friendly_aerosol_number_concentration_updated_by_physics',
            'tendency_of_water_friendly_aerosols_at_surface',
            'tendency_of_ice_friendly_aerosols_at_surface',
            ],
        },
    #'subroutine_name_1' : 'all',
    #'subroutine_name_2' : 'none',
    #'subroutine_name_2' : [ 'var1', 'var3'],
    }

# Names of Fortran include files in the host model cap (do not change);
# both files will be written to the directory of each target file, only
# used by the dynamic builds
MODULE_INCLUDE_FILE = 'ccpp_modules_{set}.inc'
FIELDS_INCLUDE_FILE = 'ccpp_fields_{set}.inc'

# Directory where to write static API to
STATIC_API_DIR = 'FV3/gfsphysics/CCPP_layer'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = 'ccpp/physics/CCPP_VARIABLES_FV3.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = 'ccpp/framework/doc/DevelopersGuide/CCPP_VARIABLES_FV3.tex'


###############################################################################
# Template code to generate include files                                     #
###############################################################################

# Name of the CCPP data structure in the host model cap;
# in the case of FV3, this is a 2-dimensional array with
# the number of blocks as the first and the number of
# OpenMP threads as the second dimension; nb is the loop
# index for the current block, nt for the current thread.
# Internally, the model uses an associate construct to
# reference cdata(nb,nt) with cdata (recommended).
CCPP_DATA_STRUCTURE = 'cdata'
