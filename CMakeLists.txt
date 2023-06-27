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
list(APPEND VALID_APPS ATM ATMAERO ATMAQ ATMW ATMWM ATML LND S2S S2SA S2SW S2SWA S2SWAL HAFS HAFSW HAFS-ALL NG-GODAS)
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
set(NOAHMP          OFF CACHE BOOL "Enable NOAHMP")

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
message("NOAHMP ........... ${NOAHMP}")

###############################################################################
### Build Options
###############################################################################
set(32BIT           OFF CACHE BOOL "Enable 32BIT (single precision arithmetic in dycore and fast physics)")
set(CCPP_32BIT      OFF CACHE BOOL "Enable CCPP_32BIT (single precision arithmetic in slow physics)")
set(AVX2            ON  CACHE BOOL "Enable AVX2 instruction set")
set(AVX             OFF CACHE BOOL "Enable AVX-I instruction set")
set(SIMDMULTIARCH   OFF CACHE BOOL "Enable multi-target SIMD instruction sets")
set(DEBUG           OFF CACHE BOOL "Enable DEBUG mode")
set(DISABLE_FMA     OFF CACHE BOOL "Disable Fused Multiply-Add instructions (workaround needed for AMD EPYC)" FORCE)
set(INLINE_POST     ON  CACHE BOOL "Enable inline post")
set(MULTI_GASES     OFF CACHE BOOL "Enable MULTI_GASES")
set(MOVING_NEST     OFF CACHE BOOL "Enable moving nest code")
set(OPENMP          ON  CACHE BOOL "Enable OpenMP threading")
set(PARALLEL_NETCDF OFF CACHE BOOL "Enable parallel NetCDF")
set(JEDI_DRIVER     OFF CACHE BOOL "Enable JEDI as top level driver")
set(CMEPS_AOFLUX    OFF CACHE BOOL "Enable atmosphere-ocean flux calculation in mediator")
set(PDLIB           OFF CACHE BOOL "Enable Domain Decomposition in WW3 via PDLIB")

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

option(BUILD_WITH_IFI "Build NCEPpost with In-Flight Icing (IFI) library if present" OFF)
option(REQUIRE_IFI "Abort if libIFI is not found ; enables BUILD_WITH_IFI=ON" OFF)
option(INTERNAL_IFI "Compile with IFI inside the executable, instead of a library" OFF)

if(INTERNAL_IFI)
  message("Building with libIFI inside executable. Do not distribute executable outside of WCOSS2.")
  set(BUILD_WITH_IFI OFF)
  set(REQUIRE_IFI OFF)
endif()

if(REQUIRE_IFI)
  set(BUILD_WITH_IFI ON)
endif()

if(BUILD_WITH_IFI)
  if(REQUIRE_IFI)
    find_package(IFI CONFIG REQUIRED)
  else()
    find_package(IFI CONFIG)
  endif()
endif()

message("")
message("32BIT ............ ${32BIT}")
message("CCPP_32BIT ....... ${CCPP_32BIT}")
message("AVX2 ............. ${AVX2}")
message("SIMDMULTIARCH .... ${SIMDMULTIARCH}")
message("DEBUG ............ ${DEBUG}")
message("INLINE_POST ...... ${INLINE_POST}")
message("MULTI_GASES ...... ${MULTI_GASES}")
message("MOVING_NEST ...... ${MOVING_NEST}")
message("OPENMP ........... ${OPENMP}")
message("PARALLEL_NETCDF .. ${PARALLEL_NETCDF}")
message("JEDI_DRIVER ...... ${JEDI_DRIVER}")
message("CMEPS_AOFLUX ..... ${CMEPS_AOFLUX}")

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

find_package(NetCDF 4.7.4 REQUIRED C Fortran)
find_package(ESMF 8.3.0 MODULE REQUIRED)
if(FMS)
  find_package(FMS 2022.04 REQUIRED COMPONENTS R4 R8)
  if(APP MATCHES "^(HAFSW)$")
    add_library(fms ALIAS FMS::fms_r4)
  elseif (APP MATCHES "^(S2S|S2SA|S2SW|S2SWA|S2SWAL|NG-GODAS)$")
    add_library(fms ALIAS FMS::fms_r8)
  endif()
  if(APP MATCHES "^(ATM|ATMAERO|ATMAQ|ATMWM|ATMW|ATML|LND|HAFS|HAFS-ALL)$")
    if(32BIT)
      add_library(fms ALIAS FMS::fms_r4)
    else()
      add_library(fms ALIAS FMS::fms_r8)
    endif()
  endif()
endif()

if(CMEPS)
  find_package(PIO 2.5.3 REQUIRED COMPONENTS C Fortran)
endif()

find_package(bacio 2.4.0 REQUIRED)
find_package(sp 2.3.3 REQUIRED)
find_package(w3emc 2.9.2 REQUIRED)

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
  if(APP MATCHES "^(ATMWM)$")
    set(SWITCH "multi_esmf" CACHE STRING "ESMF cap")
    set(UFS_CAP "MULTI_ESMF" CACHE STRING "Build with MULTI_ESMF cap")
    list(APPEND _ufs_defs_private FRONT_WW3=WMESMFMD)
  else()
    if(PDLIB)
      set(SWITCH "meshcap_pdlib" CACHE STRING "NUOPC mesh cap")
    else()
      set(SWITCH "meshcap" CACHE STRING "NUOPC mesh cap")
    endif()
    set(UFS_CAP "NUOPC_MESH" CACHE STRING "Build with NUOPC_MESH cap")
    list(APPEND _ufs_defs_private FRONT_WW3=wav_comp_nuopc)
  endif()
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
### Land Components [NOAHMP]
###############################################################################
if(NOAHMP)
  add_subdirectory(NOAHMP-interface)
endif()

###############################################################################
### UFS Library
###############################################################################
add_library(ufs driver/UFSDriver.F90)
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
  target_link_libraries(ufs PUBLIC WW3::WW3)
endif()

list(APPEND _ufs_libs_public esmf)

if(STOCH_PHYS)
  list(APPEND _ufs_libs_public stochastic_physics)
endif()

if(FV3)
  add_dependencies(ufs fv3atm)
  list(APPEND _ufs_defs_private FRONT_FV3=fv3atm_cap_mod)
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

if(NOAHMP)
  add_dependencies(ufs noahmp)
  list(APPEND _ufs_defs_private FRONT_NOAHMP=lnd_comp_nuopc)
  list(APPEND _ufs_libs_public noahmp)
endif()

target_compile_definitions(ufs PRIVATE "${_ufs_defs_private}")
target_link_libraries(ufs PUBLIC "${_ufs_libs_public}")

###############################################################################
### UFS executable
###############################################################################
add_executable(ufs_model driver/UFS.F90)
add_dependencies(ufs_model ufs)
target_link_libraries(ufs_model ufs esmf w3emc::w3emc_d)
set_target_properties(ufs_model PROPERTIES LINKER_LANGUAGE Fortran)

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
