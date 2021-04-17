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

project(ufs
        VERSION 1.0
        LANGUAGES C CXX Fortran)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules/Modules)

set(32BIT           OFF CACHE BOOL "Enable 32BIT (single precision arithmetic in dycore)")
set(AVX2            ON  CACHE BOOL "Enable AVX2 instruction set")
set(SIMDMULTIARCH   OFF CACHE BOOL "Enable multi-target SIMD instruction sets")
set(DEBUG           OFF CACHE BOOL "Enable DEBUG mode")
set(DEBUG_LINKMPI   ON  CACHE BOOL "Enable linkmpi option when DEBUG mode is on")
set(INLINE_POST     OFF CACHE BOOL "Enable inline post")
set(MULTI_GASES     OFF CACHE BOOL "Enable MULTI_GASES")
set(OPENMP          ON  CACHE BOOL "Enable OpenMP threading")
set(PARALLEL_NETCDF OFF CACHE BOOL "Enable parallel NetCDF")
set(REPRO           OFF CACHE BOOL "Enable REPRO mode")
set(JEDI_DRIVER     OFF CACHE BOOL "Enable JEDI as top level driver")

# Component Options
set(FMS             OFF CACHE BOOL "Enable FMS")
set(FV3             OFF CACHE BOOL "Enable FV3")
set(MOM6            OFF CACHE BOOL "Enable MOM6")
set(CICE6           OFF CACHE BOOL "Enable CICE6")
set(WW3             OFF CACHE BOOL "Enable WW3")
set(STOCH_PHYS      OFF CACHE BOOL "Enable Stochastic Physics")
set(NEMSdatm        OFF CACHE BOOL "Enable NEMSdatm")
set(CMEPS           OFF CACHE BOOL "Enable CMEPS")
set(CDEPS           OFF CACHE BOOL "Enable CDEPS")

# Valid applications and choices
list(APPEND VALID_APPS ATM ATMW S2S S2SW DATM DATM_NEMS)
set(APP NONE CACHE BOOL "Application Name")
if(NOT (APP IN_LIST VALID_APPS))
  message(FATAL_ERROR "${APP} is not a valid application.\nValid Applications are: ${VALID_APPS}")
endif()

###############################################################################
### Set CMAKE_BUILD_TYPE for DEBUG mode
###############################################################################
if(DEBUG)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Set type of build to Debug." FORCE)
else()
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Set type of build to Release." FORCE)
endif()

if(CMAKE_Platform)
  message("")
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)
    message("Setting configuration for ${CMAKE_Platform}")
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)
  else()
    message("Platform '${CMAKE_Platform}' configuration file does not exist")
  endif()

  # DH* 20210208 temporary workaround for FMS issues on Cheyenne with GNU 9.1.0
  if (CMAKE_Platform MATCHES "cheyenne.gnu")
    message("Applying patch cheyenne_gnu_fms_mpp_util_mpi_inc.patch")
    execute_process(COMMAND patch -N -p0 INPUT_FILE cheyenne_gnu_fms_mpp_util_mpi_inc.patch
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/cheyenne_gnu_fms_mpp_util_mpi_inc.patch.out
                    ERROR_FILE  ${CMAKE_CURRENT_BINARY_DIR}/cheyenne_gnu_fms_mpp_util_mpi_inc.patch.err
                    RESULT_VARIABLE RC)
    # *DH 20210208
  endif()
endif()

message("")
message("AVX2 ............. ${AVX2}")
message("SIMDMULTIARCH ...  ${SIMDMULTIARCH}")
message("DEBUG ............ ${DEBUG}")
message("DEBUG_LINKMPI .... ${DEBUG_LINKMPI}")
message("INLINE_POST ...... ${INLINE_POST}")
message("MULTI_GASES ...... ${MULTI_GASES}")
message("OPENMP ........... ${OPENMP}")
message("PARALLEL_NETCDF .. ${PARALLEL_NETCDF}")
message("REPRO ............ ${REPRO}")
message("")

# Configure selected application specific components
include(cmake/configure_apps.cmake)

# Components
message("")
message("FMS .............. ${FMS}")
message("FV3 .............. ${FV3}")
message("MOM6 ............. ${MOM6}")
message("CICE6 ............ ${CICE6}")
message("WW3 .............. ${WW3}")
message("STOCH_PHYS ....... ${STOCH_PHYS}")
message("NEMSdatm ......... ${NEMSdatm}")
message("CDEPS ............ ${CDEPS}")
message("CMEPS ............ ${CMEPS}")

# Options that applications reset
message("32BIT ............ ${32BIT}")
message("")

get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message("C       compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION} (${C_COMPILER_NAME})")
message("CXX     compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} (${CXX_COMPILER_NAME})")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION} (${Fortran_COMPILER_NAME})")
message("")

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0.0)
  message(FATAL_ERROR "GNU Compiler >= 9 is required")
endif()

include(cmake/${CMAKE_Fortran_COMPILER_ID}.cmake)
message("")

###############################################################################
### Find Dependencies
###############################################################################
find_package(MPI REQUIRED)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()
find_package(NetCDF REQUIRED C Fortran)
find_package(ESMF MODULE REQUIRED)
if(CMEPS)
  find_package(PIO REQUIRED COMPONENTS C Fortran STATIC)
endif()

find_package(bacio REQUIRED)
find_package(nemsio REQUIRED)
find_package(sp REQUIRED)
find_package(w3emc REQUIRED)
find_package(w3nco REQUIRED)
if(INLINE_POST)
  find_package(upp REQUIRED)
endif()

# See https://github.com/NOAA-EMC/NCEPLIBS-nemsio/pull/22
target_link_libraries(nemsio::nemsio INTERFACE w3emc::w3emc_d bacio::bacio_4)

# Configure Python
find_package(Python 3 REQUIRED COMPONENTS Interpreter)
message("Found Python: ${Python_EXECUTABLE}")

###############################################################################
### FMS
###############################################################################
if(FMS)
  set(GFS_PHYS ON CACHE BOOL "Enable GFS Physics")
  if(NOT 32BIT)
    set(64BIT ON CACHE BOOL "Enable 64-bit")
  endif()

  add_subdirectory(FMS)

  if(32BIT)
    add_library(fms ALIAS fms_r4)
  else()
    add_library(fms ALIAS fms_r8)
  endif()
endif(FMS)

###############################################################################
### stochastic_physics
###############################################################################
if(STOCH_PHYS)
  add_subdirectory(stochastic_physics)
endif()

###############################################################################
### Atmosphere Components [FV3, MPAS?]
###############################################################################
if(FV3)
  add_subdirectory(FV3)
endif()

###############################################################################
### Wave components [WW3]
###############################################################################
if(WW3)
  if(CMAKE_Platform)
    if(${CMAKE_Platform} STREQUAL "wcoss2")
      set(WW3_COMP "wcoss_cray")
    elseif(${CMAKE_Platform} MATCHES "linux*" OR ${CMAKE_Platform} MATCHES "macosx*")
      if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(WW3_COMP "gnu")
      elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(WW3_COMP "intel")
      endif()
    else()
      set(WW3_COMP ${CMAKE_Platform})
    endif()
  endif()
  message("Build WW3:")
  message("   run: ${CMAKE_BUILD_TOOL} WW3_PARCOMPN=4 WW3_COMP=${WW3_COMP} ww3_nemslibonly")
  message("   in:  ${CMAKE_CURRENT_SOURCE_DIR}/WW3/model/esmf")
  message("")

  add_custom_target(ww3_nems
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/WW3/model/esmf
    COMMAND ${CMAKE_BUILD_TOOL} WW3_PARCOMPN=4 WW3_COMP=${WW3_COMP} ww3_nemslibonly > ${CMAKE_CURRENT_BINARY_DIR}/ww3_make.log 2>&1)

  set(WW3_LIBS ${CMAKE_CURRENT_SOURCE_DIR}/WW3/model/obj/libww3_multi_esmf.a)
endif()

###############################################################################
### Marine Components [MOM6, CICE6]
###############################################################################
if(MOM6)
  add_subdirectory(MOM6-interface)
endif()

if(CICE6)
  add_subdirectory(CICE-interface)
endif()

###############################################################################
### Mediator Components [CMEPS]
###############################################################################
if(CMEPS)
  add_subdirectory(CMEPS-interface)
endif()

###############################################################################
### Data Components [NEMSdatm, CDEPS]
###############################################################################
if(NEMSdatm)
  add_subdirectory(DATM)
endif()

if(CDEPS)
  add_subdirectory(CDEPS-interface)
endif()

###############################################################################
### UFS Library [NEMS]
###############################################################################
list(APPEND _nems_srcs NEMS/src/module_NEMS_UTILS.F90
                       NEMS/src/module_EARTH_INTERNAL_STATE.F90
                       NEMS/src/module_EARTH_GRID_COMP.F90
                       NEMS/src/module_NEMS_INTERNAL_STATE.F90
                       NEMS/src/module_NEMS_GRID_COMP.F90
                       NEMS/src/module_NEMS_Rusage.F90
                       NEMS/src/nems_c_rusage.c)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/NEMS/src/ESMFVersionDefine_ESMF_NUOPC.h
               ${CMAKE_CURRENT_BINARY_DIR}/NEMS/src/ESMFVersionDefine.h
               COPYONLY)
# until we remove "../ESMFVersionDefine.h" from ENS_Cpl/ENS_CplComp_ESMFMod_STUB.F90
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/NEMS/src/ESMFVersionDefine_ESMF_NUOPC.h
               ${CMAKE_CURRENT_BINARY_DIR}/NEMS/ESMFVersionDefine.h
               COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/NEMS/src/ESMFConvenienceMacros.h
               ${CMAKE_CURRENT_BINARY_DIR}/NEMS/src/ESMFConvenienceMacros.h
               COPYONLY)

add_library(ufs ${_nems_srcs})
set_target_properties(ufs PROPERTIES Fortran_MODULE_DIRECTORY
                                     ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(ufs PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/NEMS/src)
target_include_directories(ufs INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                         $<INSTALL_INTERFACE:mod>)

list(APPEND _ufs_defs_private ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})

if(JEDI_DRIVER)
  list(APPEND _ufs_defs_private JEDI_DRIVER=ON)
endif()

if(WW3)
  add_dependencies(ufs ww3_nems)
  list(APPEND _ufs_defs_private FRONT_WW3=WMESMFMD)
  target_include_directories(ufs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/WW3/model/mod)
  list(APPEND _ufs_libs_public ${WW3_LIBS})
endif()

list(APPEND _ufs_libs_public esmf)

if(STOCH_PHYS)
  list(APPEND _ufs_libs_public stochastic_physics)
endif()

if(FMS)
  add_dependencies(ufs fms)
  list(APPEND _ufs_defs_private FRONT_FMS)
endif()

if(NEMSdatm)
  add_dependencies(ufs datatm)
  list(APPEND _ufs_defs_private FRONT_DATM=datm)
  list(APPEND _ufs_libs_public datatm)
endif()

if(FV3)
  add_dependencies(ufs fv3atm)
  list(APPEND _ufs_defs_private FRONT_FV3=fv3gfs_cap_mod)
  list(APPEND _ufs_libs_public fv3atm)
endif()

if(MOM6)
  add_dependencies(ufs mom6)
  list(APPEND _ufs_defs_private FRONT_MOM6=mom_cap_mod)
  list(APPEND _ufs_libs_public mom6)
endif()

if(CICE6)
  add_dependencies(ufs cice)
  list(APPEND _ufs_defs_private FRONT_CICE6=ice_comp_nuopc)
  list(APPEND _ufs_libs_public cice)
endif()

if(CMEPS)
  add_dependencies(ufs cmeps)
  list(APPEND _ufs_defs_private CMEPS FRONT_CMEPS=MED)
  list(APPEND _ufs_libs_public cmeps)
endif()

if(CDEPS)
  add_dependencies(ufs cdeps::datm)
endif()

target_compile_definitions(ufs PRIVATE "${_ufs_defs_private}")
target_link_libraries(ufs PUBLIC "${_ufs_libs_public}")

###############################################################################
### UFS executable
###############################################################################
add_executable(ufs_model NEMS/src/MAIN_NEMS.F90)
add_dependencies(ufs_model ufs)

list(APPEND _ufs_model_defs_private ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})

set_target_properties(ufs_model PROPERTIES Fortran_MODULE_DIRECTORY
                                           ${CMAKE_CURRENT_BINARY_DIR}/mod)

target_include_directories(ufs_model PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/NEMS/src)

if(WW3)
  list(APPEND _ufs_model_defs_private FRONT_WW3=WMESMFMD)
  target_include_directories(ufs_model PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/WW3/model/mod)
endif()

target_compile_definitions(ufs_model PRIVATE "${_ufs_model_defs_private}")

if(NEMSdatm OR CDEPS)
  target_link_libraries(ufs_model PUBLIC ufs w3nco::w3nco_d)
endif()

target_link_libraries(ufs_model PRIVATE ufs
                                        esmf
                                        NetCDF::NetCDF_Fortran)

###############################################################################
### Install
###############################################################################
install(TARGETS ufs
        EXPORT  ufs-config
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(EXPORT      ufs-config
        DESTINATION lib/cmake)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

###############################################################################
### Done
###############################################################################
