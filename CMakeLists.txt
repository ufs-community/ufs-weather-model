cmake_minimum_required(VERSION 3.15)

foreach(env_var IN ITEMS CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_Fortran_COMPILER)
  if(NOT DEFINED ENV{${env_var}})
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

set(32BIT           OFF CACHE BOOL "Enable 32BIT (single precision arithmetic in dycore)")
set(AVX2            ON  CACHE BOOL "Enable AVX2 instruction set")
set(SIMDMULTIARCH   OFF CACHE BOOL "Enable multi-target SIMD instruction sets")
set(CCPP            ON  CACHE BOOL "Enable CCPP")
set(DEBUG           OFF CACHE BOOL "Enable DEBUG mode")
set(INLINE_POST     OFF CACHE BOOL "Enable inline post")
set(MULTI_GASES     OFF CACHE BOOL "Enable MULTI_GASES")
set(OPENMP          ON  CACHE BOOL "Enable OpenMP threading")
set(PARALLEL_NETCDF OFF CACHE BOOL "Enable parallel NetCDF")
set(QUAD_PRECISION  ON  CACHE BOOL "Enable QUAD_PRECISION (for certain grid metric terms in dycore)")
set(REPRO           OFF CACHE BOOL "Enable REPRO mode")
set(WW3             OFF CACHE BOOL "Enable WW3")

if(CMAKE_Platform)
  message("")
  if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)
    message("Setting configuration for ${CMAKE_Platform}")
    include(${PROJECT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)
  else()
    message("Platform '${CMAKE_Platform}' configuration file does not exist")
  endif()
endif()

message("")
message("32BIT ............ ${32BIT}")
message("AVX2 ............. ${AVX2}")
message("SIMDMULTIARCH      ${SIMDMULTIARCH}")
message("CCPP ............. ${CCPP}")
message("DEBUG ............ ${DEBUG}")
message("INLINE_POST ...... ${INLINE_POST}")
message("MULTI_GASES ...... ${MULTI_GASES}")
message("OPENMP ........... ${OPENMP}")
message("PARALLEL_NETCDF .. ${PARALLEL_NETCDF}")
message("QUAD_PRECISION ... ${QUAD_PRECISION}")
message("REPRO ............ ${REPRO}")
message("WW3 .............. ${WW3}")
message("")

get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("C       compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} (${C_COMPILER_NAME})")
message("CXX     compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CXX_COMPILER_NAME})")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${Fortran_COMPILER_NAME})")
message("")

include(cmake/${CMAKE_Fortran_COMPILER_ID}.cmake)
message("")

find_package(NetCDF REQUIRED C Fortran)
find_package(MPI REQUIRED)
find_package(ESMF MODULE REQUIRED)

find_package(bacio REQUIRED)
find_package(nemsio REQUIRED)
find_package(sp REQUIRED)
find_package(w3emc REQUIRED)
find_package(w3nco REQUIRED)
if(INLINE_POST)
  find_package(nceppost REQUIRED)
endif()
# See https://github.com/NOAA-EMC/NCEPLIBS-nemsio/pull/22
target_link_libraries(nemsio::nemsio INTERFACE w3emc::w3emc_d bacio::bacio_4)

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
include(fms_files.cmake)
add_library(fms ${fms_src_files})
target_include_directories(fms PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/include>
                                      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/fms>
                                      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/fms2_io/include>
                                      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/FMS/mpp/include>)
target_link_libraries(fms PRIVATE NetCDF::NetCDF_Fortran)
set_target_properties(fms PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS)
set_target_properties(fms PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)
target_include_directories(fms PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/FMS/mod>)

###############################################################################
### stochastic_physics
###############################################################################
add_subdirectory(stochastic_physics)

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
    COMMAND ${CMAKE_BUILD_TOOL} WW3_PARCOMPN=4 WW3_COMP=${WW3_COMP} ww3_nemslibonly > ${PROJECT_BINARY_DIR}/ww3_make.log 2>&1
  )

    set(WW3_LIBS ${PROJECT_SOURCE_DIR}/WW3/model/obj/libww3_multi_esmf.a)
endif()

###############################################################################
### NEMS (ufs-weather-model-lib)
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
add_dependencies(ufs-weather-model-lib fv3atm)
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

target_link_libraries(ufs-weather-model-lib PUBLIC fv3atm)

###############################################################################
### UFS executable
###############################################################################
add_executable(NEMS.exe
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

target_link_libraries(NEMS.exe
                      ufs-weather-model-lib
                      ${WW3_LIBS}
                      esmf NetCDF::NetCDF_Fortran ${MKL_LIB})

###############################################################################
### Install
###############################################################################
install(
  TARGETS ufs-weather-model-lib fms
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
