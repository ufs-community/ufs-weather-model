#!/usr/bin/env python

# CCPP prebuild config for GMTB Single Column Model (SCM) v2.0


###############################################################################
# Definitions                                                                 #
###############################################################################

# Add all files with metadata tables on the host model side,
# relative to basedir = top-level directory of host model
VARIABLE_DEFINITION_FILES = [
    'scm/src/gmtb_scm_type_defs.f90',
    'scm/src/gmtb_scm_physical_constants.f90'
    ]

# Can be empty, since all physics schemes and their
# dependencies are hardcoded in CMakeLists in
# ccpp-physics - to fix, c.f. FV3 v1
SCHEME_FILES_DEPENDENCIES = []

# Add all physics scheme files relative to basedir
SCHEME_FILES = [
    'ccpp-physics/GFS_layer/GFS_initialize_scm.F90',
    'ccpp-physics/GFS_layer/GFS_finalize_scm.F90',
    'ccpp-physics/physics/GFS_DCNV_generic.f90',
    'ccpp-physics/physics/GFS_MP_generic_post.f90',
    'ccpp-physics/physics/GFS_MP_generic_pre.f90',
    'ccpp-physics/physics/GFS_PBL_generic.f90',
    'ccpp-physics/physics/GFS_SCNV_generic.f90',
    'ccpp-physics/physics/GFS_calpreciptype.f90',
    'ccpp-physics/physics/GFS_phys_time_vary.scm.f90',
    'ccpp-physics/physics/GFS_rad_time_vary.scm.f90',
    'ccpp-physics/physics/GFS_rrtmg_post.F90',
    'ccpp-physics/physics/GFS_rrtmg_pre.F90',
    'ccpp-physics/physics/GFS_suite_interstitial.ccpp.F90',
    'ccpp-physics/physics/GFS_surface_generic.f90',
    'ccpp-physics/physics/GFS_surface_loop_control.f',
    'ccpp-physics/physics/GFS_zhao_carr_pre.f90',
    'ccpp-physics/physics/cnvc90.f',
    'ccpp-physics/physics/dcyc2.f',
    'ccpp-physics/physics/get_prs_fv3.f90',
    'ccpp-physics/physics/gscond.f',
    'ccpp-physics/physics/gwdc.f',
    'ccpp-physics/physics/gwdps.f',
    'ccpp-physics/physics/mfdeepcnv.f',
    'ccpp-physics/physics/mfshalcnv.f',
    'ccpp-physics/physics/moninedmf.f',
    'ccpp-physics/physics/ozphys.f',
    'ccpp-physics/physics/precpd.f',
    'ccpp-physics/physics/radlw_main.f',
    'ccpp-physics/physics/radsw_main.f',
    'ccpp-physics/physics/rayleigh_damp.f',
    'ccpp-physics/physics/rrtmg_lw_post.F90',
    'ccpp-physics/physics/rrtmg_lw_pre.F90',
    'ccpp-physics/physics/rrtmg_sw_post.F90',
    'ccpp-physics/physics/rrtmg_sw_pre.F90',
    'ccpp-physics/physics/sfc_diag.f',
    'ccpp-physics/physics/sfc_diff.f',
    'ccpp-physics/physics/sfc_drv.f',
    'ccpp-physics/physics/sfc_nst.f',
    'ccpp-physics/physics/sfc_sice.f',
    ]

# Auto-generated makefile/cmakefile snippets that contain all schemes
SCHEMES_MAKEFILE = 'ccpp-physics/CCPP_SCHEMES.mk'
SCHEMES_CMAKEFILE = 'ccpp-physics/CCPP_SCHEMES.cmake'

# CCPP host cap in which to insert the ccpp_field_add statements;
# determines the directory to place ccpp_{modules,fields}.inc
TARGET_FILES = [
    'scm/src/gmtb_scm.f90',
    ]

# Auto-generated makefile/cmakefile snippets that contain all caps
CAPS_MAKEFILE = 'ccpp-physics/CCPP_CAPS.mk'
CAPS_CMAKEFILE = 'ccpp-physics/CCPP_CAPS.cmake'

# Directory where to put all auto-generated physics caps
CAPS_DIR = 'ccpp-physics/physics'

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
HTML_VARTABLE_FILE = 'ccpp-physics/CCPP_VARIABLES_SCM.html'

# LaTeX document containing the provided vs requested CCPP variables
LATEX_VARTABLE_FILE = 'ccpp-framework/doc/DevelopersGuide/CCPP_VARIABLES_SCM.tex'


###############################################################################
# Template code to generate include files                                     #
###############################################################################

# Name of the CCPP data structure in the host model cap;
# in the case of SCM, this is a vector with loop index i
CCPP_DATA_STRUCTURE = 'cdata(i)'

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
                               GFS_sfcprop_type,                        &
                               GFS_coupling_type, GFS_control_type,     &
                               GFS_grid_type,     GFS_tbd_type,         &
                               GFS_cldprop_type,  GFS_radtend_type,     &
                               GFS_diag_type,     GFS_interstitial_type,&
                               GFS_init_type
'''
