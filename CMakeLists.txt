cmake_minimum_required(VERSION 3.15)

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

project(ufs-weather-model
        VERSION 1.0
        LANGUAGES C CXX Fortran)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

message("")
message("Setting configuration for $ENV{CMAKE_Platform}")
message("")
include(${PROJECT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)

get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("C       compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} (${C_COMPILER_NAME})")
message("CXX     compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CXX_COMPILER_NAME})")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${Fortran_COMPILER_NAME})")
message("")

include(cmake/${CMAKE_Fortran_COMPILER_ID}.cmake)
message("")

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

find_package(NetCDF REQUIRED C Fortran)
find_package(MPI REQUIRED)
find_package(ESMF MODULE REQUIRED)

find_package(bacio REQUIRED)
find_package(w3emc REQUIRED)
find_package(sigio REQUIRED)
find_package(sp REQUIRED)
find_package(nemsio REQUIRED)
find_package(g2 REQUIRED)
find_package(g2tmpl REQUIRED)

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

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0.0)
  message(FATAL_ERROR "GNU Compiler >= 9 is required")
endif()

if(CCPP)

    find_package(Python 3 QUIET COMPONENTS Interpreter)
    if (NOT Python_Interpreter_FOUND)
      find_package(Python 2.7 QUIET REQUIRED COMPONENTS Interpreter)
    endif()
    message("Found Python: ${Python_EXECUTABLE}")

    if(DEFINED CCPP_SUITES)
        message("Calling CCPP code generator (ccpp_prebuild.py) for SUITES = ${CCPP_SUITES}")
        execute_process(COMMAND FV3/ccpp/framework/scripts/ccpp_prebuild.py
                                "--config=FV3/ccpp/config/ccpp_prebuild_config.py"
                                "--suites=${CCPP_SUITES}"
                                "--builddir=${PROJECT_BINARY_DIR}/FV3"
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                        OUTPUT_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.out
                        ERROR_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.err
                        RESULT_VARIABLE RC
                       )
    else()
        message("Calling CCPP code generator (ccpp_prebuild.py) ...")
        execute_process(COMMAND FV3/ccpp/framework/scripts/ccpp_prebuild.py
                                "--config=FV3/ccpp/config/ccpp_prebuild_config.py"
                                "--builddir=${PROJECT_BINARY_DIR}/FV3"
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                        OUTPUT_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.out
                        ERROR_FILE ${PROJECT_BINARY_DIR}/ccpp_prebuild.err
                        RESULT_VARIABLE RC
                       )
    endif()

    # Check return code from ccpp_prebuild.py
    if(RC EQUAL 0)
        message("")
    else()
        message(FATAL_ERROR "An error occured while running ccpp_prebuild.py, check ${PROJECT_BINARY_DIR}/ccpp_prebuild.{out,err}")
    endif()

    # this should not be necessary; including CCPP_*.cmake here and passing
    # SCHEMES, CAPS and TYPEDEFS via environment variables to CCPP build.
    # CCPP should be able to directly include those three .cmake files.

    include(${PROJECT_BINARY_DIR}/FV3/ccpp/physics/CCPP_SCHEMES.cmake)
    include(${PROJECT_BINARY_DIR}/FV3/ccpp/physics/CCPP_CAPS.cmake)
    include(${PROJECT_BINARY_DIR}/FV3/ccpp/physics/CCPP_TYPEDEFS.cmake)

    set(ENV{CCPP_SCHEMES} "${SCHEMES}")
    set(ENV{CCPP_CAPS} "${CAPS}")
    set(ENV{CCPP_TYPEDEFS} "${TYPEDEFS}")

endif()

###############################################################################
### FMS
###############################################################################
add_library(
    fms

    FMS/amip_interp/amip_interp.F90
    FMS/astronomy/astronomy.F90
    FMS/axis_utils/axis_utils.F90
    FMS/axis_utils/axis_utils2.F90
    FMS/block_control/block_control.F90
    FMS/column_diagnostics/column_diagnostics.F90
    FMS/constants/constants.F90
    FMS/coupler/atmos_ocean_fluxes.F90
    FMS/coupler/coupler_types.F90
    FMS/coupler/ensemble_manager.F90
    FMS/data_override/data_override.F90
    FMS/diag_integral/diag_integral.F90
    FMS/diag_manager/diag_axis.F90
    FMS/diag_manager/diag_data.F90
    FMS/diag_manager/diag_grid.F90
    FMS/diag_manager/diag_manager.F90
    FMS/diag_manager/diag_output.F90
    FMS/diag_manager/diag_table.F90
    FMS/diag_manager/diag_util.F90
    FMS/drifters/cloud_interpolator.F90
    FMS/drifters/drifters.F90
    FMS/drifters/drifters_comm.F90
    FMS/drifters/drifters_core.F90
    FMS/drifters/drifters_input.F90
    FMS/drifters/drifters_io.F90
    FMS/drifters/quicksort.F90
    FMS/exchange/stock_constants.F90
    FMS/exchange/xgrid.F90
    FMS/field_manager/field_manager.F90
    FMS/field_manager/fm_util.F90
    FMS/fms/fms.F90
    FMS/fms/fms_io.F90
    FMS/fms2_io/blackboxio.F90
    FMS/fms2_io/fms2_io.F90
    FMS/fms2_io/fms_io_utils.F90
    FMS/fms2_io/fms_netcdf_domain_io.F90
    FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90
    FMS/fms2_io/netcdf_io.F90
    FMS/horiz_interp/horiz_interp.F90
    FMS/horiz_interp/horiz_interp_bicubic.F90
    FMS/horiz_interp/horiz_interp_bilinear.F90
    FMS/horiz_interp/horiz_interp_conserve.F90
    FMS/horiz_interp/horiz_interp_spherical.F90
    FMS/horiz_interp/horiz_interp_type.F90
    FMS/interpolator/interpolator.F90
    FMS/memutils/memutils.F90
    FMS/monin_obukhov/monin_obukhov.F90
    FMS/monin_obukhov/monin_obukhov_kernel.F90
    FMS/mosaic/gradient.F90
    FMS/mosaic/grid.F90
    FMS/mosaic/mosaic.F90
    FMS/mosaic/mosaic2.F90
    FMS/mpp/mpp.F90
    FMS/mpp/mpp_data.F90
    FMS/mpp/mpp_domains.F90
    FMS/mpp/mpp_efp.F90
    FMS/mpp/mpp_io.F90
    FMS/mpp/mpp_memutils.F90
    FMS/mpp/mpp_parameter.F90
    FMS/mpp/mpp_pset.F90
    FMS/mpp/mpp_utilities.F90
    FMS/platform/platform.F90
    FMS/random_numbers/MersenneTwister.F90
    FMS/random_numbers/random_numbers.F90
    FMS/sat_vapor_pres/sat_vapor_pres.F90
    FMS/sat_vapor_pres/sat_vapor_pres_k.F90
    FMS/time_interp/time_interp.F90
    FMS/time_interp/time_interp_external.F90
    FMS/time_interp/time_interp_external2.F90
    FMS/time_manager/get_cal_time.F90
    FMS/time_manager/time_manager.F90
    FMS/topography/gaussian_topog.F90
    FMS/topography/topography.F90
    FMS/tracer_manager/tracer_manager.F90
    FMS/tridiagonal/tridiagonal.F90

    FMS/mosaic/create_xgrid.c
    FMS/mosaic/gradient_c2l.c
    FMS/mosaic/interp.c
    FMS/mosaic/mosaic_util.c
    FMS/mosaic/read_mosaic.c
    FMS/mpp/mpp_memuse.c
    FMS/mpp/nsclock.c
    FMS/mpp/threadloc.c
)
# stupid cmake can not figure out dependency of fft.F90 on fft99.F90 because 'use fft99_mod' is inside ifdefs
set_property(SOURCE FMS/fft/fft.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-DSGICRAY=0 -DNAGFFT=0")

target_include_directories(fms PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/include>
                                      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/fms>
                                      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/fms2_io/include>
                                      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/mpp/include>)
target_link_libraries(fms NetCDF::NetCDF_Fortran)
set_target_properties(fms PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS)
set_target_properties(fms PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)
target_include_directories(fms PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/FMS/mod>)

###############################################################################
### FV3
###############################################################################
add_subdirectory(FV3)

###############################################################################
### WW3
###############################################################################
if(WW3)
  set(WW3_COMP ${CMAKE_Platform})
  if(${CMAKE_Platform} STREQUAL "hera.intel")
    set(WW3_COMP "hera")
  endif()
  if(${CMAKE_Platform} STREQUAL "orion.intel")
    set(WW3_COMP "orion")
  endif()
  message("Build WW3:")
  message("   run: ${CMAKE_BUILD_TOOL} WW3_PARCOMPN=4 WW3_COMP=${WW3_COMP} ww3_nemslibonly")
  message("   in:  ${PROJECT_SOURCE_DIR}/WW3/model/esmf")
  message("")

  add_custom_target(ww3_nems
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/WW3/model/esmf
    COMMAND ${CMAKE_BUILD_TOOL} WW3_PARCOMPN=4 WW3_COMP=${WW3_COMP} ww3_nemslibonly > make.log 2>&1
  )

    set(WW3_LIBS ${PROJECT_SOURCE_DIR}/WW3/model/obj/libww3_multi_esmf.a)
endif()

###############################################################################
### NEMS (i.e. ufs-weather-model-lib)
###############################################################################
add_library(
    ufs-weather-model-lib

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
add_dependencies(ufs-weather-model-lib fv3cap)
target_link_libraries(ufs-weather-model-lib PRIVATE esmf)

target_compile_definitions(ufs-weather-model-lib PRIVATE -DESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})
target_compile_definitions(ufs-weather-model-lib PRIVATE -DFRONT_FMS)
target_compile_definitions(ufs-weather-model-lib PRIVATE -DFRONT_FV3=fv3gfs_cap_mod)

target_include_directories(ufs-weather-model-lib PRIVATE ${PROJECT_BINARY_DIR}/NEMS/src)
target_include_directories(ufs-weather-model-lib PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/FV3/mod)
target_include_directories(ufs-weather-model-lib PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)

if(WW3)
    target_compile_definitions(ufs-weather-model-lib PRIVATE -DFRONT_WW3=WMESMFMD)
    target_include_directories(ufs-weather-model-lib PRIVATE ${PROJECT_SOURCE_DIR}/WW3/model/mod)
    set_target_properties(ufs-weather-model-lib PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/WW3/model/mod)
    add_dependencies(ufs-weather-model-lib ww3_nems)
endif()

set_target_properties(ufs-weather-model-lib PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
target_include_directories(ufs-weather-model-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/FV3/mod>
                                       $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/mod>)

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

###############################################################################
### ufs-weather-model
###############################################################################
add_executable(
    NEMS.exe

    NEMS/src/MAIN_NEMS.F90
)
add_dependencies(NEMS.exe ufs-weather-model-lib)

target_compile_definitions(NEMS.exe PRIVATE -DESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})
target_include_directories(NEMS.exe PRIVATE ${PROJECT_BINARY_DIR}/NEMS/src)

set_target_properties(NEMS.exe PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)

if(WW3)
    target_compile_definitions(NEMS.exe PRIVATE -DFRONT_WW3=WMESMFMD)
    set_target_properties(NEMS.exe PROPERTIES Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/WW3/model/mod)
    add_dependencies(NEMS.exe ww3_nems)
endif()

if(CCPP)
    set(CCPP_LIBRARIES ccppdriver ccppphys ccpp)
endif()

target_link_libraries(NEMS.exe
                      ufs-weather-model-lib
                      fv3cap
                      fv3core
                      io
                      ipd
                      gfsphysics
                      ${CCPP_LIBRARIES}
                      cpl
                      stochastic_physics
                      fms
                      ${WW3_LIBS}
                      $ENV{POST_LIB}
                      nemsio::nemsio g2::g2_4 g2tmpl::g2tmpl w3emc::w3emc_d sp::sp_d
                      esmf NetCDF::NetCDF_Fortran ${MKL_LIB})

target_link_libraries(ufs-weather-model-lib INTERFACE
                      fv3cap
                      fv3core
                      io
                      ipd
                      gfsphysics
                      ${CCPP_LIBRARIES}
                      cpl
                      stochastic_physics
                      fms
                      ${WW3_LIBS}
                      $ENV{POST_LIB}
                      nemsio::nemsio g2::g2_4 g2tmpl::g2tmpl w3emc::w3emc_d sp::sp_d
                      esmf NetCDF::NetCDF_Fortran ${MKL_LIB})

###############################################################################
### Install
###############################################################################

install(
  TARGETS ufs-weather-model-lib fv3cap fv3core io ipd gfsphysics ${CCPP_LIBRARIES} cpl stochastic_physics fms
  EXPORT ufs-weather-model-lib-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT ufs-weather-model-lib-config
  DESTINATION lib/cmake
)

###############################################################################
### done
###############################################################################
