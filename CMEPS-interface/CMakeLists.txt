# CMEPS compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv")
  set(CMAKE_Fortran_LINK_FLAGS "")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

list(APPEND _nems_util_files
  CMEPS/nems/util/shr_abort_mod.F90
  CMEPS/nems/util/shr_log_mod.F90
  CMEPS/nems/util/shr_pio_mod.F90
  CMEPS/nems/util/shr_sys_mod.F90
  CMEPS/nems/util/shr_flux_mod.F90
  CMEPS/nems/util/shr_mpi_mod.F90
  CMEPS/nems/util/glc_elevclass_mod.F90
  CMEPS/nems/util/shr_mem_mod.F90
  CMEPS/nems/util/shr_kind_mod.F90
  CMEPS/nems/util/perf_mod.F90
  CMEPS/nems/util/shr_const_mod.F90)

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
  CMEPS/mediator/med_io_mod.F90)

list(APPEND _mediator_defs ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR}
                           ESMF_VERSION_MINOR=${ESMF_VERSION_MINOR}
                           INTERNAL_PIO_INIT)

### Use common object library for building target library
add_library(nems_util_obj OBJECT ${_nems_util_files})
set_target_properties(nems_util_obj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(nems_util_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CMEPS/nems/util>)
target_link_libraries(nems_util_obj PRIVATE esmf
                                            PIO::PIO_C PIO::PIO_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(nems_util_obj PRIVATE OpenMP::OpenMP_Fortran)
endif()

add_library(mediator_obj OBJECT ${_mediator_files})
set_target_properties(mediator_obj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(mediator_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CMEPS/mediator>)
target_compile_definitions(mediator_obj PRIVATE "${_mediator_defs}")
target_link_libraries(mediator_obj PRIVATE nems_util_obj
                                           esmf
                                           PIO::PIO_C PIO::PIO_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(mediator_obj PRIVATE OpenMP::OpenMP_Fortran)
endif()

### Create target library and set PUBLIC interfaces on the library
add_library(cmeps STATIC $<TARGET_OBJECTS:nems_util_obj>
                         $<TARGET_OBJECTS:mediator_obj>)
target_compile_definitions(cmeps PUBLIC "${_mediator_defs}")
target_include_directories(cmeps PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                        $<INSTALL_INTERFACE:mod>)
target_link_libraries(cmeps PUBLIC esmf
                                   PIO::PIO_C PIO::PIO_Fortran)
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
