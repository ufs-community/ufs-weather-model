#!/usr/bin/env python

# CCPP prebuild config for GFDL Finite-Volume Cubed-Sphere Model (FV3) v0


###############################################################################
# Definitions                                                                 #
###############################################################################

# Add all files with metadata tables on the host model side,
# relative to basedir = top-level directory of host model
VARIABLE_DEFINITION_FILES = [
    'FV3/gfsphysics/GFS_layer/GFS_typedefs.F90',
    'FV3/gfsphysics/physics/physcons.f90',
    ]

# Add all physics scheme dependencies relative to basedir - note that these are all violations
# of the CCPP requirement to not use any external modules except Fortran standard modules!
SCHEME_FILES_DEPENDENCIES = [
    'FV3/gfsphysics/physics/GFDL_parse_tracers.F90',
    'FV3/gfsphysics/physics/GFS_phys_time_vary.fv3.f90',
    'FV3/gfsphysics/physics/GFS_rad_time_vary.fv3.f90',
    'FV3/gfsphysics/physics/GFS_radupdate.f90',
    'FV3/gfsphysics/physics/date_def.f',
    'FV3/gfsphysics/physics/funcphys.f90',
    'FV3/gfsphysics/physics/gfs_phy_tracer_config.f',
    'FV3/gfsphysics/physics/gocart_tracer_config_stub.f',
    'FV3/gfsphysics/physics/gcycle.f90',
    'FV3/gfsphysics/physics/h2o_def.f',
    'FV3/gfsphysics/physics/h2ointerp.f90',
    'FV3/gfsphysics/physics/iounitdef.f',
    'FV3/gfsphysics/physics/machine.F',
    'FV3/gfsphysics/physics/mersenne_twister.f',
    'FV3/gfsphysics/physics/mfpbl.f',
    'FV3/gfsphysics/physics/module_bfmicrophysics.f',
    'FV3/gfsphysics/physics/module_nst_model.f90',
    'FV3/gfsphysics/physics/module_nst_parameters.f90',
    'FV3/gfsphysics/physics/module_nst_water_prop.f90',
    'FV3/gfsphysics/physics/namelist_soilveg.f',
    'FV3/gfsphysics/physics/num_parthds.F',
    'FV3/gfsphysics/physics/ozinterp.f90',
    'FV3/gfsphysics/physics/ozne_def.f',
    'FV3/gfsphysics/physics/physcons.f90',
    'FV3/gfsphysics/physics/physparam.f',
    'FV3/gfsphysics/physics/radcons.f90',
    'FV3/gfsphysics/physics/rad_initialize.f',
    'FV3/gfsphysics/physics/radiation_aerosols.f',
    'FV3/gfsphysics/physics/radiation_astronomy.f',
    'FV3/gfsphysics/physics/radiation_clouds.f',
    'FV3/gfsphysics/physics/radiation_gases.f',
    'FV3/gfsphysics/physics/radiation_surface.f',
    'FV3/gfsphysics/physics/radlw_datatb.f',
    'FV3/gfsphysics/physics/radlw_param.f',
    'FV3/gfsphysics/physics/radsw_datatb.f',
    'FV3/gfsphysics/physics/radsw_param.f',
    'FV3/gfsphysics/physics/rascnvv2.f',
    'FV3/gfsphysics/physics/set_soilveg.f',
    'FV3/gfsphysics/physics/sfcsub.F',
    'FV3/gfsphysics/physics/sflx.f',
    'FV3/gfsphysics/physics/tridi.f',
    'FV3/gfsphysics/physics/wam_f107_kp_mod.f90',
    'FV3/gfsphysics/GFS_layer/GFS_abstraction_layer.F90',
    'FV3/gfsphysics/GFS_layer/GFS_diagnostics.F90',
    'FV3/gfsphysics/GFS_layer/GFS_driver.F90',
    'FV3/gfsphysics/GFS_layer/GFS_radiation_driver.F90',
    'FV3/gfsphysics/GFS_layer/GFS_restart.F90',
    'FV3/gfsphysics/GFS_layer/GFS_typedefs.F90',
    'FV3/gfsphysics/IPD_layer/IPD_CCPP_driver.F90',
    'FV3/gfsphysics/IPD_layer/IPD_driver.F90',
    'FV3/gfsphysics/IPD_layer/IPD_driver_cap.F90',
    'FV3/gfsphysics/IPD_layer/IPD_typedefs.F90',
    ]

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    'FV3/gfsphysics/physics/GFS_DCNV_generic.f90',
    'FV3/gfsphysics/physics/GFS_MP_generic_post.f90',
    'FV3/gfsphysics/physics/GFS_MP_generic_pre.f90',
    'FV3/gfsphysics/physics/GFS_PBL_generic.f90',
    'FV3/gfsphysics/physics/GFS_SCNV_generic.f90',
    'FV3/gfsphysics/physics/GFS_calpreciptype.f90',
    'FV3/gfsphysics/physics/GFS_debug.F90',
    'FV3/gfsphysics/physics/GFS_rrtmg_post.F90',
    'FV3/gfsphysics/physics/GFS_rrtmg_pre.F90',
    'FV3/gfsphysics/physics/GFS_stochastics.f90',
    'FV3/gfsphysics/physics/GFS_suite_interstitial.ccpp.F90',
    'FV3/gfsphysics/physics/GFS_surface_generic.f90',
    'FV3/gfsphysics/physics/GFS_surface_loop_control.f',
    'FV3/gfsphysics/physics/GFS_zhao_carr_pre.f90',
    'FV3/gfsphysics/physics/cnvc90.f',
    'FV3/gfsphysics/physics/dcyc2.f',
    'FV3/gfsphysics/physics/get_prs_fv3.f90',
    'FV3/gfsphysics/physics/gscond.f',
    'FV3/gfsphysics/physics/gwdc.f',
    'FV3/gfsphysics/physics/gwdps.f',
    'FV3/gfsphysics/physics/memcheck.F90',
    'FV3/gfsphysics/physics/mfdeepcnv.f',
    'FV3/gfsphysics/physics/mfshalcnv.f',
    'FV3/gfsphysics/physics/moninedmf.f',
    'FV3/gfsphysics/physics/ozphys.f',
    'FV3/gfsphysics/physics/precpd.f',
    'FV3/gfsphysics/physics/radlw_main.f',
    'FV3/gfsphysics/physics/radsw_main.f',
    'FV3/gfsphysics/physics/rayleigh_damp.f',
    'FV3/gfsphysics/physics/rrtmg_lw_post.F90',
    'FV3/gfsphysics/physics/rrtmg_lw_pre.F90',
    'FV3/gfsphysics/physics/rrtmg_sw_post.F90',
    'FV3/gfsphysics/physics/rrtmg_sw_pre.F90',
    'FV3/gfsphysics/physics/sfc_diag.f',
    'FV3/gfsphysics/physics/sfc_diff.f',
    'FV3/gfsphysics/physics/sfc_drv.f',
    'FV3/gfsphysics/physics/sfc_nst.f',
    'FV3/gfsphysics/physics/sfc_sice.f',
    ]

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE = 'FV3/gfsphysics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE = 'FV3/gfsphysics/CCPP_SCHEMES.cmake'

# CCPP host cap in which to insert the ccpp_field_add statements;
# determines the directory to place ccpp_{modules,fields}.inc
TARGET_FILES = [
    'FV3/gfsphysics/IPD_layer/IPD_CCPP_Driver.F90',
    ]

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE = 'FV3/gfsphysics/CCPP_CAPS.mk'
CAPS_CMAKEFILE = 'FV3/gfsphysics/CCPP_CAPS.cmake'

# Directory where to put all auto-generated physics caps
CAPS_DIR = 'FV3/gfsphysics/physics'

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
    #'subroutine_name_1' : 'all',
    #'subroutine_name_2' : 'none',
    #'subroutine_name_2' : [ 'var1', 'var3'],
    }

# Names of Fortran include files in the host model cap (do not change);
# both files will be written to the directory of each target file
MODULE_INCLUDE_FILE = 'ccpp_modules.inc'
FIELDS_INCLUDE_FILE = 'ccpp_fields.inc'

# HTML document containing the model-defined CCPP variables
HTML_VARTABLE_FILE = 'FV3/gfsphysics/CCPP_VARIABLES_FV3.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = 'ccpp-framework/doc/DevelopersGuide/CCPP_VARIABLES_FV3.tex'


###############################################################################
# Template code to generate include files                                     #
###############################################################################

# Name of the CCPP data structure in the host model cap;
# in the case of FV3, this is a 2-dimensional array with
# the number of blocks as the first and the number of
# OpenMP threads as the second dimension; nb is the loop
# index for the current block, nt for the current thread
CCPP_DATA_STRUCTURE = 'cdata_block(nb,nt)'

# Modules to load for auto-generated ccpp_field_add code
# in the host model cap (e.g. error handling)
MODULE_USE_TEMPLATE_HOST_CAP = \
'''
use ccpp_errors, only: ccpp_error
'''

# Modules to load for auto-generated ccpp_field_get code
# in the physics scheme cap (e.g. derived data types)
MODULE_USE_TEMPLATE_SCHEME_CAP = \
'''
       use machine, only: kind_phys
       use module_radlw_parameters, only: sfcflw_type, topflw_type
       use module_radsw_parameters, only: cmpfsw_type, sfcfsw_type, topfsw_type
       use GFS_typedefs, only: GFS_statein_type,  GFS_stateout_type,    &
                               GFS_sfcprop_type,  GFS_coupling_type,    &
                               GFS_control_type,  GFS_grid_type,        &
                               GFS_tbd_type,      GFS_cldprop_type,     &
                               GFS_radtend_type,  GFS_diag_type,        &
                               GFS_interstitial_type
'''
