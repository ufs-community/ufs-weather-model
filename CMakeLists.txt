cmake_minimum_required(VERSION 2.8.12)

# disable in source builds to protect against overwriting other Makefiles
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
if( NOT DEFINED ENV{CC} )
  find_path( crayComp "ftn" )
  find_path( wcossIntel "mpfort" )
  find_path( intelComp "ifort" )
  find_path( pgiComp "pgf90" )
  if( crayComp )
    message("Setting CrayLinuxEnvironment")
    set(CMAKE_SYSTEM_NAME "CrayLinuxEnvironment")
    set(CMAKE_C_COMPILER "${crayComp}/cc")
    set(CMAKE_CXX_COMPILER "${crayComp}/CC")
    set(CMAKE_Fortran_COMPILER "${crayComp}/ftn")
  endif()
  if( intelComp )
    set(ENV{CC} "icc")
    set(ENV{CXX} "icpc")
    set(ENV{FC} "ifort")
  endif()
  if( wcossIntel )
    message("Setting env for wcoss intel")
    set(ENV{CC} "mpcc")
    set(ENV{CXX} "mpCC")
    set(ENV{FC} "mpfort")
  endif()
  if( pgiComp )
    set(ENV{CC} "pgcc")
    set(ENV{CXX} "pgCC")
    set(ENV{FC} "pgf90")
  endif()
endif()

project(NEMSfv3gfs C CXX Fortran)
if( DEFINED ENV{CMAKE_Platform} )
  set(CMAKE_Platform $ENV{CMAKE_Platform})
  string(TOLOWER ${CMAKE_Fortran_COMPILER_ID} COMPILER_TYPE)  
  include(${PROJECT_SOURCE_DIR}/cmake/configure_${CMAKE_Platform}.cmake)
else()
  string(TOLOWER ${CMAKE_SYSTEM_NAME} OS)
  string(TOLOWER ${CMAKE_Fortran_COMPILER_ID} COMPILER_TYPE)  
  set(CMAKE_Platform configure_${OS}.${COMPILER_TYPE}.cmake)
  message("platform will be ${CMAKE_Platform}")
endif()
if( DEFINED NCEPLIBS_INSTALL_DIR )
  find_path( LIBXML2_INCLUDE_DIR NAMES libxml2 HINTS ${NCEPLIBS_INSTALL_DIR}/include)
  find_path( LIBXML2_LIB_DIR NAMES libxml2.so HINTS ${NCEPLIBS_INSTALL_DIR}/lib)
  message("xml are ${LIBXML2_INCLUDE_DIR} and ${LIBXML2_LIB_DIR} ")
endif()
if(DEFINED ${ESMFMKFILE} )
  message("setting esmf mk file via env variable")
  set(ENV{ESMFMKFILE} ${ESMFMKFILE})
else()
   if( DEFINED NCEPLIBS_INSTALL_DIR )
      message("setting esmf mk file to point to nceplibs install")
      set(ENV{ESMFMKFILE} ${NCEPLIBS_INSTALL_DIR}/lib/esmf.mk)
   endif()  
endif()  
if(DEFINED ENV{NETCDF})
  message("netcdf is defined")
  set(NETCDF_DIR ${NETCDF})
else()
   message("netcdf is NOT defined")
   if( NCEPLIBS_INSTALL_DIR )
      message("setting NETCDF dir to point to nceplibs install")
      set(NETCDF_DIR ${NCEPLIBS_INSTALL_DIR})
      set(ENV{NETCDF_DIR} ${NCEPLIBS_INSTALL_DIR})
      message("HEY!! before find netcdf ${NETCDF_DIR}")
   endif()  
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/../cmake/Modules/")

message("HEY system name is ${CMAKE_SYSTEM_NAME}")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(MPI REQUIRED)
#set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
#set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
#set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
find_package(ESMF REQUIRED)
#set(ESMF_MOD ${ESMF_F90COMPILEPATHS})
#set(ESMF_LIB ${ESMF_LIBSDIR}/libesmf.so)

set(NETCDF_F90 YES)
set(NETCDF_C YES)
find_package(NetCDF REQUIRED )
find_package(BACIO)
find_package(NEMSIO)
find_package(SP)
find_package(W3EMC)
find_package(W3NCO)


if(NOT DEFINED CCPP_SUITES)
  set(CCPP_SUITES FV3_GFS_2017_gfdlmp)
endif()
find_package( PythonInterp 2.7 REQUIRED )
find_package( PythonLibs 2.7 REQUIRED )

#set( prebuildcommand  ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/FV3/ccpp/framework/scripts/ccpp_prebuild.py --config=./FV3/ccpp/config/ccpp_prebuild_config.py --static --suites=${CCPP_SUITES} --builddir=${CMAKE_CURRENT_BINARY_DIR}/FV3 >  ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.log)
set( prebuildcommand ${PROJECT_SOURCE_DIR}/genCCPP.sh ${PROJECT_BINARY_DIR} )
message("${prebuildcommand}")
execute_process(COMMAND ${prebuildcommand} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE STDOUT )
message("preprocess output was ${STDOUT}")

#set( CCPPLOG ${CMAKE_CURRENT_BINARY_DIR}/ccpp_prebuild.log )
#ADD_CUSTOM_COMMAND( OUTPUT ${CCPPLOG}
#    PRE_BUILD
#    COMMAND ${prebuildcommand}
#    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#    COMMENT "HEY! running ccpp_prebuild"
#    DEPENDS ${CCPPLOG}
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/FV3/ccpp/framework/scripts/ccpp_prebuild.py
#)
#add_custom_target(ccpp_prebuild ALL DEPENDS ${CCPPLOG} )
set(CCPP ON)
set(STATIC ON)
set(SUITES ${CCPP_SUITES})
if(NOT DEFINED PHYS)
  set(PHYS gfs)
endif()
message("")
message("Selected physics package: ${PHYS}")
message("")

if(${COMPILER_TYPE} STREQUAL "gnu")
  include(${CMAKE_SOURCE_DIR}/cmake/GNU.cmake)
  setGNU()
endif()
if(${COMPILER_TYPE} STREQUAL "intel")
  include(${CMAKE_SOURCE_DIR}/cmake/Intel.cmake)
  setIntel()
endif()
message("HEY fort flags are ${CMAKE_Fortran_FLAGS}")
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
    FMS/affinity/affinity.c
    FMS/mpp/nsclock.c
    FMS/mpp/threadloc.c
)
#add_dependencies(fms run )
#add_dependencies( fms ccpp_prebuild )
# stupid cmake can not figure out dependency of fft.F90 on fft99.F90 because 'use fft99_mod' is inside ifdefs
set_property(SOURCE FMS/fft/fft.F90 APPEND_STRING PROPERTY COMPILE_FLAGS "-DSGICRAY=0 -DNAGFFT=0")

target_include_directories(fms PUBLIC ${PROJECT_SOURCE_DIR}/FMS/include FMS/include FMS/fms FMS/mpp/include)
#set(NETCDF_INC_DIR ${NETCDF_INCLUDES})
message("HEY, netcdf includes are ${NETCDF_INCLUDES_F90} and ${NETCDF_INCLUDES}")
#target_include_directories(fms PRIVATE ${NETCDF_INC_DIR})
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(fms PRIVATE ${NETCDF_INCLUDES} ${NETCDF_INCLUDES_F90} ${MPI_Fortran_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} "./" ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})
link_directories(${MPI_Fortran_LIBRARIES} ${ARCHIVE_OUTPUT_PATH} )

set_target_properties(fms PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS)
set_target_properties(fms PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)
target_include_directories(fms PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/FMS/mod)

set(delim "=")
file(READ ${PROJECT_BINARY_DIR}//FV3/ccpp/physics/CCPP_SCHEMES.sh CCPP_SCHEMES_FILE)
string(FIND "${CCPP_SCHEMES_FILE}" "${delim}" SPLIT_POS)
string(LENGTH "${CCPP_SCHEMES_FILE}" FILE_LEN)
MATH(EXPR SPLIT_POS "${SPLIT_POS}+2")
MATH(EXPR LINE_LEN "${FILE_LEN}-${SPLIT_POS} -2")
string(SUBSTRING "${CCPP_SCHEMES_FILE}" ${SPLIT_POS} ${LINE_LEN} CCPP_FILES)
set(ENV{CCPP_SCHEMES} "${CCPP_FILES}")
file(READ ${PROJECT_BINARY_DIR}//FV3/ccpp/physics/CCPP_CAPS.sh CCPP_CAPS_FILE)
string(FIND "${CCPP_CAPS_FILE}" "${delim}" SPLIT_POS)
string(LENGTH "${CCPP_CAPS_FILE}" FILE_LEN)
MATH(EXPR SPLIT_POS "${SPLIT_POS}+2")
MATH(EXPR LINE_LEN "${FILE_LEN}-${SPLIT_POS} -2")
string(SUBSTRING "${CCPP_CAPS_FILE}" ${SPLIT_POS} ${LINE_LEN} CCPP_FILES)
set(ENV{CCPP_CAPS} "${CCPP_FILES}")
file(READ ${PROJECT_BINARY_DIR}//FV3/ccpp/physics/CCPP_STATIC_API.sh CCPP_API_FILE)
string(FIND "${CCPP_API_FILE}" "${delim}" SPLIT_POS)
string(LENGTH "${CCPP_API_FILE}" FILE_LEN)
MATH(EXPR SPLIT_POS "${SPLIT_POS}+2")
MATH(EXPR LINE_LEN "${FILE_LEN}-${SPLIT_POS} -2")
string(SUBSTRING "${CCPP_API_FILE}" ${SPLIT_POS} ${LINE_LEN} CCPP_API)
set(ENV{CCPP_STATIC_API} "${CCPP_API}")
set(ENV{CCPP_INCLUDE_DIRS} "${NETCDF_INCLUDES}")
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
target_include_directories(NEMS.exe PRIVATE ${MPI_Fortran_INCLUDE_DIRS} )
target_include_directories(NEMS.exe PRIVATE ${NCEPLIBS_INSTALL_DIR}/include )
target_include_directories(NEMS.exe PRIVATE ${NCEPLIBS_INSTALL_DIR}/mod )
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
set(NCEP_LIBS ${NEMSIO_LIB};${BACIO_LIB4};${SP_LIBd};${W3EMC_LIBd};${W3NCO_LIBd})
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
                      ${NEMSIO_LIB} 
                      ${NCEP_LIBS} ${ESMF_LIB} ${NETCDF_LIBRARIES} ${NETCDF_LIBRARIES_F90} ${MPI_Fortran_LIBRARIES} ${MKL_LIB})

###############################################################################
### done
###############################################################################
