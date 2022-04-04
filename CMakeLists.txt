cmake_minimum_required(VERSION 3.19)

project(ufs
        VERSION 1.0
        LANGUAGES C CXX Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU" AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9.0.0)
  message(FATAL_ERROR "GNU Compiler >= 9 is required")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules/Modules)

###############################################################################
### App and Component Options
###############################################################################

# Valid applications and choices
list(APPEND VALID_APPS ATM ATMAERO ATMAQ ATMW S2S S2SA S2SW S2SWA HAFS HAFSW HAFS-ALL NG-GODAS)
set(APP NONE CACHE BOOL "Application Name")
if(NOT (APP IN_LIST VALID_APPS))
  message(FATAL_ERROR "${APP} is not a valid application.\nValid Applications are: ${VALID_APPS}")
endif()

set(FMS             OFF CACHE BOOL "Enable FMS")
set(FV3             OFF CACHE BOOL "Enable FV3")
set(AQM             OFF CACHE BOOL "Enable AQM")
set(UFS_GOCART      OFF CACHE BOOL "Enable GOCART")
set(MOM6            OFF CACHE BOOL "Enable MOM6")
set(HYCOM           OFF CACHE BOOL "Enable HYCOM")
set(CICE6           OFF CACHE BOOL "Enable CICE6")
set(WW3             OFF CACHE BOOL "Enable WW3")
set(STOCH_PHYS      OFF CACHE BOOL "Enable Stochastic Physics")
set(CMEPS           OFF CACHE BOOL "Enable CMEPS")
set(CDEPS           OFF CACHE BOOL "Enable CDEPS")

# Configure selected application specific components
message("")
include(cmake/configure_apps.cmake)

message("")
message("FMS .............. ${FMS}")
message("FV3 .............. ${FV3}")
message("AQM .............. ${AQM}")
message("GOCART ........... ${UFS_GOCART}")
message("MOM6 ............. ${MOM6}")
message("HYCOM ............ ${HYCOM}")
message("CICE6 ............ ${CICE6}")
message("WW3 .............. ${WW3}")
message("STOCH_PHYS ....... ${STOCH_PHYS}")
message("CDEPS ............ ${CDEPS}")
message("CMEPS ............ ${CMEPS}")

###############################################################################
### Build Options
###############################################################################
set(32BIT           OFF CACHE BOOL "Enable 32BIT (single precision arithmetic in dycore)")
set(AVX2            ON  CACHE BOOL "Enable AVX2 instruction set")
set(AVX             OFF CACHE BOOL "Enable AVX-I instruction set")
set(SIMDMULTIARCH   OFF CACHE BOOL "Enable multi-target SIMD instruction sets")
set(DEBUG           OFF CACHE BOOL "Enable DEBUG mode")
set(DEBUG_LINKMPI   ON  CACHE BOOL "Enable linkmpi option when DEBUG mode is on")
set(INLINE_POST     OFF CACHE BOOL "Enable inline post")
set(MULTI_GASES     OFF CACHE BOOL "Enable MULTI_GASES")
set(MOVING_NEST     OFF CACHE BOOL "Enable moving nest code")
set(OPENMP          ON  CACHE BOOL "Enable OpenMP threading")
set(PARALLEL_NETCDF OFF CACHE BOOL "Enable parallel NetCDF")
set(REPRO           OFF CACHE BOOL "Enable REPRO mode")
set(JEDI_DRIVER     OFF CACHE BOOL "Enable JEDI as top level driver")

set(CMAKE_Platform $ENV{CMAKE_Platform})
if(CMAKE_Platform)
  message("")
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)
    message("Setting configuration for ${CMAKE_Platform}")
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)
  else()
    message("Platform '${CMAKE_Platform}' configuration file does not exist")
  endif()
endif()

message("")
message("32BIT ............ ${32BIT}")
message("AVX2 ............. ${AVX2}")
message("SIMDMULTIARCH .... ${SIMDMULTIARCH}")
message("DEBUG ............ ${DEBUG}")
message("DEBUG_LINKMPI .... ${DEBUG_LINKMPI}")
message("INLINE_POST ...... ${INLINE_POST}")
message("MULTI_GASES ...... ${MULTI_GASES}")
message("MOVING_NEST ...... ${MOVING_NEST}")
message("OPENMP ........... ${OPENMP}")
message("PARALLEL_NETCDF .. ${PARALLEL_NETCDF}")
message("REPRO ............ ${REPRO}")
message("JEDI_DRIVER ...... ${JEDI_DRIVER}")
message("")

###############################################################################
### Set CMAKE_BUILD_TYPE for DEBUG mode
###############################################################################
if(DEBUG)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Set type of build to Debug." FORCE)
else()
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Set type of build to Release." FORCE)
endif()

include(cmake/${CMAKE_Fortran_COMPILER_ID}.cmake)

###############################################################################
### Find Dependencies
###############################################################################
find_package(MPI REQUIRED)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()
find_package(NetCDF REQUIRED C Fortran)
find_package(ESMF MODULE REQUIRED)
if(FMS)
  find_package(FMS REQUIRED COMPONENTS R4 R8)
  if(32BIT)
    add_library(fms ALIAS FMS::fms_r4)
  else()
    add_library(fms ALIAS FMS::fms_r8)
  endif()
endif()
if(CMEPS)
  find_package(PIO REQUIRED COMPONENTS C Fortran STATIC)
endif()

find_package(bacio REQUIRED)
find_package(sp REQUIRED)
find_package(w3nco REQUIRED)

# Configure Python
find_package(Python 3.6 REQUIRED COMPONENTS Interpreter)
message("Found Python: ${Python_EXECUTABLE}")

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
### AQM
###############################################################################
if(AQM)
  add_subdirectory(AQM)
endif()

###############################################################################
### GOCART
###############################################################################
if(UFS_GOCART)
  add_subdirectory(GOCART)
endif()

###############################################################################
### Wave components [WW3]
###############################################################################
if(WW3)
  # Set options for WW3 build
  set(SWITCH "multi_esmf" CACHE STRING "Default to switch_multi_esmf")
  set(MULTI_ESMF ON CACHE BOOL "Enable ESMF for WW3")
  add_subdirectory(WW3)
endif()

###############################################################################
### Marine Components [MOM6, HYCOM, CICE6]
###############################################################################
if(MOM6)
  add_subdirectory(MOM6-interface)
endif()

if(HYCOM)
  add_subdirectory(HYCOM-interface)
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
### Data Components [CDEPS]
###############################################################################
if(CDEPS)
  add_subdirectory(CDEPS-interface)
endif()

###############################################################################
### UFS Library
###############################################################################
add_library(ufs driver/EARTH_GRID_COMP.F90)
set_target_properties(ufs PROPERTIES Fortran_MODULE_DIRECTORY
                                     ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(ufs INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                         $<INSTALL_INTERFACE:mod>)

if(JEDI_DRIVER)
  list(APPEND _ufs_defs_private JEDI_DRIVER=ON)
endif()

if(AQM)
  add_dependencies(ufs aqm)
  list(APPEND _ufs_defs_private FRONT_AQM=aqm)
  list(APPEND _ufs_libs_public aqm)
endif()

if(UFS_GOCART)
  add_dependencies(ufs UFS_Aerosols)
  list(APPEND _ufs_defs_private FRONT_GOCART=Aerosol_Cap)
  list(APPEND _ufs_libs_public UFS_Aerosols)
endif()

if(WW3)
  list(APPEND _ufs_defs_private FRONT_WW3=WMESMFMD)
  target_link_libraries(ufs PUBLIC WW3::WW3)
endif()

list(APPEND _ufs_libs_public esmf)

if(STOCH_PHYS)
  list(APPEND _ufs_libs_public stochastic_physics)
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

if(HYCOM)
  list(APPEND _ufs_defs_private FRONT_HYCOM=HYCOM_Mod)
  add_dependencies(ufs hycom)
  target_link_libraries(ufs PUBLIC hycom)
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
  add_dependencies(ufs cdeps::cdeps)
  list(APPEND _ufs_defs_private FRONT_CDEPS_DATM=cdeps_datm_comp)
  list(APPEND _ufs_defs_private FRONT_CDEPS_DOCN=cdeps_docn_comp)
  target_link_libraries(ufs PUBLIC cdeps::cdeps)
endif()

target_compile_definitions(ufs PRIVATE "${_ufs_defs_private}")
target_link_libraries(ufs PUBLIC "${_ufs_libs_public}")

###############################################################################
### UFS executable
###############################################################################
add_executable(ufs_model driver/UFS.F90)
add_dependencies(ufs_model ufs)
target_link_libraries(ufs_model ufs esmf w3nco::w3nco_d)

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
