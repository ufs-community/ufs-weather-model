# CMEPS compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-none")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  set(CMAKE_Fortran_LINK_FLAGS "")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

list(APPEND _nems_util_files
  CMEPS/ufs/flux_atmocn_mod.F90
  CMEPS/ufs/glc_elevclass_mod.F90
  CMEPS/ufs/ufs_kind_mod.F90
  CMEPS/ufs/perf_mod.F90
  CMEPS/ufs/ufs_const_mod.F90)

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
  CMEPS/mediator/esmFldsExchange_nems_mod.F90
  CMEPS/mediator/med_phases_prep_lnd_mod.F90
  CMEPS/mediator/med_phases_prep_atm_mod.F90
  CMEPS/mediator/med_phases_prep_ocn_mod.F90
  CMEPS/mediator/esmFlds.F90
  CMEPS/mediator/med.F90
  CMEPS/mediator/med_time_mod.F90
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
  CMEPS/mediator/med_phases_post_wav_mod.F90)

### Build required files for CCPP aoflux calculation

### Query FV3 path
get_filename_component(FV3_PATH "../FV3" REALPATH)

### If FV3 is available then use CMEPS as CCPP host model
if(EXISTS ${FV3_PATH})
  ### Print out information
  message("FV3_PATH = ${FV3_PATH}")
  message("Compiling CMEPS with CCPP support.")

  ### Set environment variable for FV3 directory, used by config file
  set(ENV{FV3_PATH} ${FV3_PATH})

  ### Parse CMAKE_FLAGS variable to get CCPP suites
  set(CMAKE_FLAGS_LIST $ENV{CMAKE_FLAGS}) 
  separate_arguments(CMAKE_FLAGS_LIST)
  foreach(myflag ${CMAKE_FLAGS_LIST})
    if(${myflag} MATCHES "-DCCPP_SUITES*")
      string(REGEX REPLACE "-DCCPP_SUITES=" "" CCPP_SUITES ${myflag})
      set(_ccpp_suites_arg "--suites=${CCPP_SUITES}")
    endif()
  endforeach()

  ### Trigger ccpp_prebuild.py script to connect host model with CCPP physics
  execute_process(
                  COMMAND ${Python_EXECUTABLE} 
                  "${FV3_PATH}/ccpp/framework/scripts/ccpp_prebuild.py"
                  "--config=CMEPS/ufs/ccpp/config/ccpp_prebuild_config.py"
                  "--builddir=${CMAKE_CURRENT_BINARY_DIR}" ${_ccpp_suites_arg}
                  "--debug"
                  "--verbose"
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.out
                  ERROR_FILE  ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.err
                  RESULT_VARIABLE return_code)

  ### Check return code from CCPP prebuild.py
  if(return_code EQUAL 0)
     message (STATUS "CCPP prebuild step completed successfully")
  else(return_code EQUAL 0)
     message (FATAL_ERROR "CCPP prebuild step failed, check ccpp_prebuild.out/ccpp_prebuild.err")
  endif(return_code EQUAL 0)

  list(APPEND _nems_util_files
    CMEPS/ufs/flux_atmocn_ccpp_mod.F90)
  #add_definitions(-DUFS_AOFLUX)
endif()

### Create target library and set PUBLIC interfaces on the library
add_library(cmeps STATIC ${_nems_util_files} ${_mediator_files} ${_ccpp_physics})
set_target_properties(cmeps PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(cmeps PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                        $<INSTALL_INTERFACE:mod>)
target_link_libraries(cmeps PUBLIC esmf
                                   PIO::PIO_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(cmeps PRIVATE OpenMP::OpenMP_Fortran)
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
