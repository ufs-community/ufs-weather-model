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
set(QUAD_PRECISION  ON  CACHE BOOL "Enable QUAD_PRECISION (for certain grid metric terms in dycore)")
set(REPRO           OFF CACHE BOOL "Enable REPRO mode")
set(WW3             OFF CACHE BOOL "Enable WW3")
set(S2S             OFF CACHE BOOL "Enable S2S")
set(JEDI_DRIVER     OFF CACHE BOOL "Enable JEDI as top level driver")
set(DATM_NEMS       OFF CACHE BOOL "Enable NEMS Data Atmosphere")
set(CDEPS           OFF CACHE BOOL "Enable CDEPS")
set(DATM            OFF CACHE BOOL "Enable CDEPS Data Atmosphere")

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
message("32BIT ............ ${32BIT}")
message("AVX2 ............. ${AVX2}")
message("SIMDMULTIARCH ...  ${SIMDMULTIARCH}")
message("DEBUG ............ ${DEBUG}")
message("DEBUG_LINKMPI .... ${DEBUG_LINKMPI}")
message("INLINE_POST ...... ${INLINE_POST}")
message("MULTI_GASES ...... ${MULTI_GASES}")
message("OPENMP ........... ${OPENMP}")
message("PARALLEL_NETCDF .. ${PARALLEL_NETCDF}")
message("QUAD_PRECISION ... ${QUAD_PRECISION}")
message("REPRO ............ ${REPRO}")
message("WW3 .............. ${WW3}")
message("S2S .............. ${S2S}")
message("DATM_NEMS......... ${DATM_NEMS}")
message("CDEPS ............ ${CDEPS}")
message("DATM ............. ${DATM}")
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
if(S2S OR CDEPS)
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
### Checks
###############################################################################
if(DATM AND NOT S2S)
   message(FATAL_ERROR "DATM=ON and S2S=OFF are incompatible, ABORT!")
endif()

if(S2S AND 32BIT)
   message(FATAL_ERROR "S2S=ON and 32BIT=ON are incompatible, ABORT!")
endif()

if(DATM_NEMS AND DATM)
   message(FATAL_ERROR "DATM_NEMS=ON and DATM=ON are incompatible, ABORT!")
endif()

###############################################################################
### FMS
###############################################################################
include(fms_files.cmake)
add_library(fms ${fms_src_files} ${fms_headers})

list(APPEND _fms_defs_public use_libMPI
                             use_netCDF
                             GFS_PHYS
                             INTERNAL_FILE_NML)
if(QUAD_PRECISION)
  list(APPEND _fms_defs_public ENABLE_QUAD_PRECISION)
endif()

# check gettid
include(CheckFunctionExists)
check_function_exists(gettid HAVE_GETTID)
if(HAVE_GETTID)
  list(APPEND _fms_defs_public HAVE_GETTID)
endif()

target_compile_definitions(fms PUBLIC "${_fms_defs_public}")

if(32BIT)
  list(APPEND _fms_defs_private OVERLOAD_R4
                                OVERLOAD_R8)
endif()

target_compile_definitions(fms PRIVATE "${_fms_defs_private}")

target_include_directories(fms PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FMS/include>
                                      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FMS/fms>
                                      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FMS/fms2_io/include>
                                      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/FMS/mod>
                                      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FMS/mpp/include>)
target_include_directories(fms INTERFACE
                                      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
                                      $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/mod>)
target_link_libraries(fms PUBLIC MPI::MPI_Fortran
                                 NetCDF::NetCDF_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(fms PRIVATE OpenMP::OpenMP_Fortran)
endif()
set_target_properties(fms PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS)
set_target_properties(fms PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)
set_target_properties(fms PROPERTIES PUBLIC_HEADER "${fms_headers}" )

###############################################################################
### stochastic_physics
###############################################################################
if(NOT DATM_NEMS AND NOT DATM)
  add_subdirectory(stochastic_physics)
endif()

###############################################################################
### FV3 or Data Atmosphere
###############################################################################
if(DATM_NEMS)
  add_subdirectory(DATM)
elseif(NOT DATM)
  add_subdirectory(FV3)
endif()

###############################################################################
### WW3
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
### S2S Components [CMEPS, MOM6, CICE]
###############################################################################
if(S2S)
  add_subdirectory(CMEPS-interface)
  add_subdirectory(MOM6-interface)
  add_subdirectory(CICE-interface)
endif()

###############################################################################
### CDEPS Component
###############################################################################
if(CDEPS)
  add_subdirectory(CDEPS)
endif()

### UFS Library [NEMS]
###############################################################################
list(APPEND _nems_srcs NEMS/src/module_NEMS_UTILS.F90
                       NEMS/src/module_MEDIATOR_methods.F90
                       NEMS/src/module_MEDIATOR.F90
                       NEMS/src/module_MEDIATOR_SpaceWeather.F90
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

list(APPEND _ufs_defs_private ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR}
                              FRONT_FMS)
add_dependencies(ufs fms)

if(DATM_NEMS)
  list(APPEND _ufs_defs_private FRONT_DATM=datm)
  add_dependencies(ufs datatm)
elseif(NOT DATM)
  list(APPEND _ufs_defs_private FRONT_FV3=fv3gfs_cap_mod)
  add_dependencies(ufs fv3atm)
endif()

if(JEDI_DRIVER)
  list(APPEND _ufs_defs_private JEDI_DRIVER=ON)
endif()
set_target_properties(ufs PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(ufs PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/NEMS/src)

if(WW3)
  list(APPEND _ufs_defs_private FRONT_WW3=WMESMFMD)
  target_include_directories(ufs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/WW3/model/mod)
  add_dependencies(ufs ww3_nems)
  target_link_libraries(ufs PUBLIC ${WW3_LIBS})
endif()

target_include_directories(ufs INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                         $<INSTALL_INTERFACE:mod>)

target_link_libraries(ufs PUBLIC esmf)

if(DATM_NEMS)
  target_link_libraries(ufs PUBLIC datatm)
elseif(DATM)
  list(APPEND _ufs_defs_private CDEPS
                                FRONT_CDEPS_DATM=atm_comp_nuopc)
  add_dependencies(ufs datm)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/CDEPS/datm)
  target_link_libraries(ufs PUBLIC datm)
else()
  target_link_libraries(ufs PUBLIC fv3atm)
endif()

if(S2S)
  list(APPEND _ufs_defs_private FRONT_MOM6=mom_cap_mod
                                FRONT_CICE6=ice_comp_nuopc
                                CMEPS
                                FRONT_CMEPS=MED)
  add_dependencies(ufs mom6 cice cmeps)
  target_link_libraries(ufs PUBLIC mom6
                                   cice
                                   cmeps)
endif()

if(CDEPS)
  add_dependencies(ufs streams dshr cdeps_share FoX_dom)
  target_link_libraries(ufs PUBLIC streams dshr cdeps_share FoX_dom)
endif()

target_compile_definitions(ufs PRIVATE "${_ufs_defs_private}")

###############################################################################
### UFS executable
###############################################################################
add_executable(ufs_model NEMS/src/MAIN_NEMS.F90)

list(APPEND _ufs_model_defs_private ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})

set_target_properties(ufs_model PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

target_include_directories(ufs_model PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/NEMS/src)

if(WW3)
  list(APPEND _ufs_model_defs_private FRONT_WW3=WMESMFMD)
  target_include_directories(ufs_model PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/WW3/model/mod)
  add_dependencies(ufs_model ww3_nems)
endif()

target_compile_definitions(ufs_model PRIVATE "${_ufs_model_defs_private}")

if(DATM_NEMS OR CDEPS)
  target_link_libraries(ufs_model PUBLIC ufs w3nco::w3nco_d)
endif()
target_link_libraries(ufs_model PRIVATE ufs
                                        esmf
                                        NetCDF::NetCDF_Fortran
                                        ${MKL_LIB})

###############################################################################
### Install
###############################################################################
install(
  TARGETS fms
  EXPORT fms-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  PUBLIC_HEADER DESTINATION include )

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT fms-config
  DESTINATION lib/cmake
)
install(
  TARGETS ufs
  EXPORT ufs-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT ufs-config
  DESTINATION lib/cmake
)
###############################################################################
### Done
###############################################################################
