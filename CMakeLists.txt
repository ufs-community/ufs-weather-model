cmake_minimum_required(VERSION 2.8.12)

foreach(env_var IN ITEMS
         CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_Fortran_COMPILER
         CMAKE_Platform)
  if ( NOT DEFINED ENV{${env_var}})
    message(FATAL_ERROR "${env_var} is not defined")
  endif()
endforeach()

set(CMAKE_C_COMPILER $ENV{CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER $ENV{CMAKE_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER $ENV{CMAKE_Fortran_COMPILER})
set(CMAKE_Platform $ENV{CMAKE_Platform})

project(NEMSfv3gfs C CXX Fortran)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(MPI REQUIRED)
find_package(ESMF REQUIRED)

include(${PROJECT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)

if(NOT DEFINED PHYS)
  set(PHYS gfs)
endif()
message("")
message("Selected physics package: ${PHYS}")
message("")

add_definitions(-Duse_libMPI)
add_definitions(-Duse_netCDF)
add_definitions(-Duse_WRTCOMP)
add_definitions(-DSPMD)
add_definitions(-DUSE_LOG_DIAG_FIELD_INFO)
#add_definitions(-Duse_LARGEFILE)
add_definitions(-DUSE_GFSL63)
add_definitions(-DGFS_PHYS)
add_definitions(-DMOIST_CAPPA)
add_definitions(-DUSE_COND)
add_definitions(-DNEW_TAUCTMAX)
add_definitions(-DINTERNAL_FILE_NML)

###############################################################################
### FMS
###############################################################################
add_library(
    fms

    FMS/amip_interp/amip_interp.F90
    FMS/astronomy/astronomy.F90
    FMS/axis_utils/axis_utils.F90
    FMS/block_control/block_control.F90
    FMS/column_diagnostics/column_diagnostics.F90
    FMS/constants/constants.F90
    FMS/coupler/atmos_ocean_fluxes.F90
    FMS/coupler/coupler_types.F90
    FMS/coupler/ensemble_manager.F90
    FMS/data_override/data_override.F90
    FMS/diag_manager/diag_axis.F90
    FMS/diag_manager/diag_data.F90
    FMS/diag_manager/diag_grid.F90
    FMS/diag_manager/diag_manager.F90
    FMS/diag_manager/diag_manifest.F90
    FMS/diag_manager/diag_output.F90
    FMS/diag_manager/diag_table.F90
    FMS/diag_manager/diag_util.F90
    FMS/drifters/cloud_interpolator.F90
    FMS/drifters/drifters_comm.F90
    FMS/drifters/drifters_core.F90
    FMS/drifters/drifters.F90
    FMS/drifters/drifters_input.F90
    FMS/drifters/drifters_io.F90
    FMS/drifters/quicksort.F90
    FMS/exchange/stock_constants.F90
    FMS/exchange/xgrid.F90
    FMS/fft/fft99.F90
    FMS/fft/fft.F90
    FMS/field_manager/field_manager.F90
    FMS/field_manager/fm_util.F90
    FMS/fms/fms.F90
    FMS/fms/fms_io.F90
    FMS/horiz_interp/horiz_interp_bicubic.F90
    FMS/horiz_interp/horiz_interp_bilinear.F90
    FMS/horiz_interp/horiz_interp_conserve.F90
    FMS/horiz_interp/horiz_interp.F90
    FMS/horiz_interp/horiz_interp_spherical.F90
    FMS/horiz_interp/horiz_interp_type.F90
    FMS/interpolator/interpolator.F90
    FMS/memutils/memutils.F90
    FMS/mosaic/gradient.F90
    FMS/mosaic/grid.F90
    FMS/mosaic/mosaic.F90
    FMS/mpp/mpp_data.F90
    FMS/mpp/mpp_domains.F90
    FMS/mpp/mpp_efp.F90
    FMS/mpp/mpp.F90
    FMS/mpp/mpp_io.F90
    FMS/mpp/mpp_memutils.F90
    FMS/mpp/mpp_parameter.F90
    FMS/mpp/mpp_pset.F90
    FMS/mpp/mpp_utilities.F90
    FMS/oda_tools/oda_core_ecda.F90
    FMS/oda_tools/oda_core.F90
    FMS/oda_tools/oda_types.F90
    FMS/oda_tools/write_ocean_data.F90
    FMS/oda_tools/xbt_drop_rate_adjust.f90
    FMS/platform/platform.F90
    FMS/random_numbers/MersenneTwister.F90
    FMS/random_numbers/random_numbers.F90
    FMS/sat_vapor_pres/sat_vapor_pres.F90
    FMS/sat_vapor_pres/sat_vapor_pres_k.F90
    FMS/station_data/station_data.F90
    FMS/time_interp/time_interp_external.F90
    FMS/time_interp/time_interp.F90
    FMS/time_manager/get_cal_time.F90
    FMS/time_manager/time_manager.F90
    FMS/topography/gaussian_topog.F90
    FMS/topography/topography.F90
    FMS/tracer_manager/tracer_manager.F90
    FMS/tridiagonal/tridiagonal.F90
    FMS/memutils/memuse.c
    FMS/mosaic/create_xgrid.c
    FMS/mosaic/gradient_c2l.c
    FMS/mosaic/interp.c
    FMS/mosaic/mosaic_util.c
    FMS/mosaic/read_mosaic.c
    FMS/mpp/affinity.c
    FMS/mpp/nsclock.c
    FMS/mpp/threadloc.c
)
# stupid cmake can not figure out dependency of fft.F90 on fft99.F90 because 'use fft99_mod' is inside ifdefs
set_property(SOURCE FMS/fft/fft.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-DSGICRAY=0 -DNAGFFT=0")

target_include_directories(fms PUBLIC ${PROJECT_SOURCE_DIR}/FMS/include FMS/include FMS/fms FMS/mpp/include)
target_include_directories(fms PRIVATE ${NETCDF_INC_DIR})

set_target_properties(fms PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS)
set_target_properties(fms PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)
target_include_directories(fms PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)

###############################################################################
### FV3
###############################################################################
add_subdirectory(FV3)

###############################################################################
### WW3
###############################################################################
if(WW3)
  if(${CMAKE_Platform} STREQUAL "theia.intel")
    set(WW3_COMP "theia")
  else()
    message(FATAL_ERROR "Platform ${CMAKE_Platform} is not supported in WW3")
  endif()
  message("Build WW3:")
  message("   run: ${CMAKE_BUILD_TOOL} WW3_PARCOMPN=4 WW3_COMP=${WW3_COMP} ww3_nems")
  message("   in:  ${PROJECT_SOURCE_DIR}/WW3/model/esmf")
  message("")

  add_custom_target(ww3_nems
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/WW3/model/esmf
    COMMAND ${CMAKE_BUILD_TOOL} WW3_PARCOMPN=4 WW3_COMP=${WW3_COMP} ww3_nems > make.log 2>&1
  )
endif()

###############################################################################
### NEMS
###############################################################################
configure_file(${PROJECT_SOURCE_DIR}/NEMS/src/ESMFVersionDefine_ESMF_NUOPC.h
               ${PROJECT_BINARY_DIR}/NEMS/src/ESMFVersionDefine.h
               COPYONLY)
# until we remove "../ESMFVersionDefine.h" from ENS_Cpl/ENS_CplComp_ESMFMod_STUB.F90
configure_file(${PROJECT_SOURCE_DIR}/NEMS/src/ESMFVersionDefine_ESMF_NUOPC.h
               ${PROJECT_BINARY_DIR}/NEMS/ESMFVersionDefine.h
               COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/NEMS/src/ESMFConvenienceMacros.h
               ${PROJECT_BINARY_DIR}/NEMS/src/ESMFConvenienceMacros.h
               COPYONLY)

add_executable(
    NEMS.exe

    NEMS/src/MAIN_NEMS.F90
    NEMS/src/module_NEMS_UTILS.F90
    NEMS/src/module_MEDIATOR_methods.F90
    NEMS/src/module_MEDIATOR.F90
    NEMS/src/module_MEDIATOR_SpaceWeather.F90
    NEMS/src/module_EARTH_INTERNAL_STATE.F90
    NEMS/src/module_EARTH_GRID_COMP.F90
    NEMS/src/module_NEMS_INTERNAL_STATE.F90
    NEMS/src/module_NEMS_GRID_COMP.F90
    NEMS/src/module_NEMS_Rusage.F90
    NEMS/src/nems_c_rusage.c
    NEMS/src/ENS_Cpl/ENS_CplComp_ESMFMod_STUB.F90
)

target_compile_definitions(NEMS.exe PRIVATE -DESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})
target_compile_definitions(NEMS.exe PRIVATE -DFRONT_FMS)
target_compile_definitions(NEMS.exe PRIVATE -DFRONT_FV3=fv3gfs_cap_mod)

target_include_directories(NEMS.exe PRIVATE ${ESMF_MOD} )
target_include_directories(NEMS.exe PRIVATE ${PROJECT_BINARY_DIR}/NEMS/src)

set_target_properties(NEMS.exe PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/NEMS/mod)

if(WW3)
    target_compile_definitions(NEMS.exe PRIVATE -DFRONT_WW3=WMESMFMD)
    set_target_properties(NEMS.exe PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/WW3/model/mod)
    set(WW3_LIBS ${PROJECT_SOURCE_DIR}/WW3/model/obj/libww3_multi_esmf.a)
    add_dependencies(NEMS.exe ww3_nems)
endif()

if(CCPP)
    set(CCPP_LIBRARIES ccppdriver ccppphys ccpp)
else()
    set(IPD_LIBRARIES ipd)
endif()

target_link_libraries(NEMS.exe
                      fv3cap
                      fv3core
                      io
                      ${IPD_LIBRARIES}
                      ${PHYS}physics
                      ${CCPP_LIBRARIES}
                      fv3cpl
                      stochastic_physics
                      fms
                      ${WW3_LIBS}
                      ${NCEP_LIBS} ${ESMF_LIBS} ${NETCDF_LIBS})

###############################################################################
### done
###############################################################################
