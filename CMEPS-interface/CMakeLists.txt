###############################################################################
### CMEPS compiler flags
###############################################################################

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace -ffree-line-length-none")
  #Set CPP defintions for the ufs/cdeps_share source component of the cmeps target library
  list(APPEND CDEPS_SHARE_DEFS "CPRGNU")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -assume realloc_lhs")
  if(FASTER)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise -assume buffered_stdout -fno-alias -align all -debug minimal -qoverride-limits -ftz -no-ip")
  else()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  endif()
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  set(CMAKE_Fortran_LINK_FLAGS "")
  #Set CPP defintions for the ufs/cdeps_share source component of the cmeps target library
  list(APPEND CDEPS_SHARE_DEFS "CPRINTEL")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

###############################################################################
### Source files
###############################################################################

list(APPEND _ufs_util_files
  CMEPS/ufs/flux_atmocn_driver_mod.F90
  CMEPS/ufs/flux_atmocn_bulk_mod.F90
  CMEPS/ufs/ufs_kind_mod.F90
  CMEPS/ufs/ufs_flux_mod.F90
  CMEPS/ufs/ufs_const_mod.F90
  CMEPS/ufs/perf_mod.F90
  CMEPS/ufs/wtracers_mod.F90
  CMEPS/ufs/glc_elevclass_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/shr_orb_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/shr_const_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/shr_abort_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/shr_strconvert_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/shr_log_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/shr_sys_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/shr_kind_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/ufs/cdeps_share/shr_assert_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/ufs/cdeps_share/shr_infnan_mod.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/CDEPS/share/nuopc_shr_methods.F90
  ${PROJECT_SOURCE_DIR}/CDEPS-interface/ufs/cdeps_share/shr_is_restart_fh_mod.F90)

list(APPEND _mediator_files
  CMEPS/mediator/med_phases_restart_mod.F90
  CMEPS/mediator/med_map_mod.F90
  CMEPS/mediator/med_methods_mod.F90
  CMEPS/mediator/med_phases_prep_ice_mod.F90
  CMEPS/mediator/med_phases_history_mod.F90
  CMEPS/mediator/med_phases_prep_glc_mod.F90
  CMEPS/mediator/med_internalstate_mod.F90
  CMEPS/mediator/med_phases_profile_mod.F90
  CMEPS/mediator/esmFldsExchange_hafs_mod.F90
  CMEPS/mediator/med_phases_prep_rof_mod.F90
  CMEPS/mediator/esmFldsExchange_cesm_mod.F90
  CMEPS/mediator/med_merge_mod.F90
  CMEPS/mediator/med_constants_mod.F90
  CMEPS/mediator/med_kind_mod.F90
  CMEPS/mediator/esmFldsExchange_ufs_mod.F90
  CMEPS/mediator/med_phases_prep_lnd_mod.F90
  CMEPS/mediator/med_phases_prep_atm_mod.F90
  CMEPS/mediator/med_phases_prep_ocn_mod.F90
  CMEPS/mediator/esmFlds.F90
  CMEPS/mediator/med.F90
  CMEPS/mediator/med_phases_ocnalb_mod.F90
  CMEPS/mediator/med_phases_prep_wav_mod.F90
  CMEPS/mediator/med_utils_mod.F90
  CMEPS/mediator/med_fraction_mod.F90
  CMEPS/mediator/med_phases_aofluxes_mod.F90
  CMEPS/mediator/med_io_mod.F90
  CMEPS/mediator/med_diag_mod.F90
  CMEPS/mediator/med_phases_post_atm_mod.F90
  CMEPS/mediator/med_phases_post_ice_mod.F90
  CMEPS/mediator/med_phases_post_ocn_mod.F90
  CMEPS/mediator/med_phases_post_lnd_mod.F90
  CMEPS/mediator/med_phases_post_glc_mod.F90
  CMEPS/mediator/med_phases_post_rof_mod.F90
  CMEPS/mediator/med_phases_post_wav_mod.F90
  CMEPS/mediator/med_ufs_trace_wrapper.F90
  CMEPS/mediator/med_field_info_mod.F90)

if(CDEPS_INLINE)
  list(APPEND _mediator_files CMEPS/mediator/med_phases_cdeps_mod.F90)
  add_definitions(-DCDEPS_INLINE)
endif()

###############################################################################
### Host model for CCPP
###############################################################################

if(CMEPS_AOFLUX)
  list(APPEND _ufs_util_files
    CMEPS/ufs/flux_atmocn_ccpp_mod.F90)

  ### Query CCPP path
  get_filename_component(UFSATM_PATH "../UFSATM" REALPATH)

  ### If UFSATM is available build CCPP host model under CMEPS too
  if(EXISTS ${UFSATM_PATH})
    ### Print out information
    message(STATUS "UFSATM_PATH = ${UFSATM_PATH}")
    message(STATUS "Compiling CMEPS with CCPP support.")

    #------------------------------------------------------------------------------
    # BEGIN call to CCPP Framework
    #------------------------------------------------------------------------------
    set(PARENT_DIR "UFSATM")
    include(${UFSATM_PATH}/cmake/ccpp_capgen.cmake)

    #------------------------------------------------------------------------------
    # Physics scheme metadata/fortran files (CCPP Physics)
    #------------------------------------------------------------------------------
    set(SCHEME_FILES_DIR "${CMAKE_SOURCE_DIR}/UFSATM/ccpp/physics/physics")

    set(SCHEME_METADATA_FILES
       "${SCHEME_FILES_DIR}/SFC_Models/Ocean/UFS/sfc_ocean.meta"
       "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diff.meta"
       "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part1.meta"
       "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_loop_control_part2.meta"
       "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_pre.meta"
       "${SCHEME_FILES_DIR}/Interstitials/UFS_SCM_NEPTUNE/GFS_surface_composites_post.meta"
       "${SCHEME_FILES_DIR}/SFC_Layer/UFS/sfc_diag.meta"
     )

    # Create list of scheme source files from metadata files
    ccpp_source_files(SCHEME_FORTRAN_FILES ${SCHEME_METADATA_FILES})
    message(STATUS "${SCHEME_METADATA_FILES} --> ${SCHEME_FORTRAN_FILES}")

    #------------------------------------------------------------------------------
    # Host-model metadata/fortran files (CMEPS)
    #------------------------------------------------------------------------------
    set(HOST "cmeps")
    set(HOST_FILES_DIR "${CMAKE_SOURCE_DIR}/CMEPS-interface/CMEPS/ufs/ccpp")

    set(HOST_METADATA_FILES
       "${HOST_FILES_DIR}/data/MED_typedefs.meta"
       "${HOST_FILES_DIR}/data/MED_data.meta"
       "${HOST_FILES_DIR}/driver/med_ccpp_driver.meta"
       "${SCHEME_FILES_DIR}/hooks/physcons.meta"
    )

    # Create list of host source files from metadata files
    ccpp_source_files(HOST_FORTRAN_FILES ${HOST_METADATA_FILES})
    message(STATUS "${HOST_METADATA_FILES} --> ${HOST_FORTRAN_FILES}")

    #------------------------------------------------------------------------------
    # Define Suite Definition Files (CMEPS)
    #------------------------------------------------------------------------------
    set(SUITE_FILES_DIR "${CMAKE_SOURCE_DIR}/CMEPS-interface/CMEPS/ufs/ccpp/suites")

    set(SUITE_FILES
       "${SUITE_FILES_DIR}/suite_FV3_sfc_ocean.xml"
    )

    #------------------------------------------------------------------------------
    # Set output directory for CCPP files.
    #------------------------------------------------------------------------------
    #set(OUTPUT_ROOT "${CMAKE_SOURCE_DIR}/CMEPS-interface/CMEPS/ufs/ccpp")
    set(OUTPUT_ROOT "${CMAKE_CURRENT_BINARY_DIR}/ccpp")

    set(CCPP_VERBOSITY "2" CACHE STRING "Verbosity level of output (default: 0)")

    #------------------------------------------------------------------------------
    # Run ccpp_validator
    #------------------------------------------------------------------------------
    ccpp_validator(VERBOSITY      ${CCPP_VERBOSITY}
                   SOURCE_FILES   ${SCHEME_FORTRAN_FILES}
                   METADATA_FILES ${SCHEME_METADATA_FILES}
                   TYPE "SCHEME")

    ccpp_validator(VERBOSITY      ${CCPP_VERBOSITY}
                   SOURCE_FILES   ${HOST_FORTRAN_FILES}
	           METADATA_FILES ${HOST_METADATA_FILES}
                   TYPE "HOST")

    #------------------------------------------------------------------------------
    # Call Capgen
    #------------------------------------------------------------------------------
    ccpp_capgen(VERBOSITY    ${CCPP_VERBOSITY}
                HOSTFILES    ${HOST_METADATA_FILES}
                SCHEMEFILES  ${SCHEME_METADATA_FILES}
                SUITES       ${SUITE_FILES}
                HOST_NAME    ${HOST}
                OUTPUT_ROOT "${OUTPUT_ROOT}"
                KIND_SPECS  "${CCPP_KINDS}")

    #------------------------------------------------------------------------------
    # Retrieve the list of Fortran files required for test host from datatable.xml;
    # this includes capgen-generated files and dependencies inferred from metadata
    #------------------------------------------------------------------------------
    ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
                  REPORT_NAME "--dependencies")
    set(CAPGEN_DEPENDENCIES ${CCPP_FILES})
    ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
                   REPORT_NAME "--scheme-files")
    set(SCHEME_FORTRAN_FILES ${CCPP_FILES})
    ccpp_datafile(DATATABLE "${OUTPUT_ROOT}/datatable.xml"
                  REPORT_NAME "--capgen-files")
    set(CAPGEN_FILES ${CCPP_FILES})

    message(STATUS "List of capgen dependencies: ${CAPGEN_DEPENDENCIES}")
    message(STATUS "List of filtered scheme files: ${SCHEME_FORTRAN_FILES}")
    message(STATUS "List of capgen-generated files: ${CAPGEN_FILES}")

    message(STATUS "CCPP Framework step complete")

    #------------------------------------------------------------------------------
    # END call to CCPP Framework
    #------------------------------------------------------------------------------

    ### Add extra driver file to the source
    list(APPEND _ufs_util_files
      CMEPS/ufs/ufs_io_mod.F90)

    ### Add special CPP flag for CMEPS/CCPP support
    add_definitions(-DUFS_AOFLUX)
  endif()
else()
  list(APPEND _ufs_util_files
    CMEPS/ufs/flux_atmocn_ccpp_stub.F90)
endif()

###############################################################################
### Libraries
###############################################################################

### Create target library
add_library(cmeps STATIC ${_ufs_util_files} ${_mediator_files} ${SCHEME_FORTRAN_FILES} ${CAPGEN_DEPENDENCIES} ${SCHEME_FORTRAN_FILES} ${HOST_FORTRAN_FILES} ${CAPGEN_FILES})

### Set PUBLIC interfaces on the library
set_target_properties(cmeps PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(cmeps PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                        $<INSTALL_INTERFACE:mod>)
target_compile_definitions(cmeps PRIVATE ${CDEPS_SHARE_DEFS})
target_link_libraries(cmeps PUBLIC ESMF::ESMF
                                   PIO::PIO_Fortran)
if(UFS_TRACING)
  target_link_libraries(cmeps PUBLIC ufs_tracing)
endif()
if(OpenMP_Fortran_FOUND)
  target_link_libraries(cmeps PRIVATE OpenMP::OpenMP_Fortran)
endif()

### Set dependency for CDEPS if inline capability is requested
if(CDEPS_INLINE)
  add_dependencies(cmeps cdeps::cdeps)
  target_link_libraries(cmeps PUBLIC cdeps::cdeps)
endif()

###############################################################################
### Install
###############################################################################

install(
  TARGETS cmeps
  EXPORT  cmeps-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/CMEPS)

install(EXPORT      cmeps-config
        DESTINATION lib/cmake)
