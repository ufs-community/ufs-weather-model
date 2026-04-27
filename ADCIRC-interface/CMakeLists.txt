###############################################################################
### Fortran compiler flags
###############################################################################

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Waliasing -fcray-pointer -fconvert=big-endian -fno-range-check -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
  set(Fortran_LINELENGTH_FLAG "-ffree-line-length-none -ffixed-line-length-none")
  set(Fortran_COMPILER_SPECIFIC_FLAG "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -sox")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  set(CMAKE_Fortran_LINK_FLAGS "")
  set(Fortran_LINELENGTH_FLAG "-132")
  set(Fortran_COMPILER_SPECIFIC_FLAG "-mcmodel=medium")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ADCIRC_Extra_FORTRAN_FLAGS}")

###############################################################################
### File lists and macros
###############################################################################

include("adcirc_files.cmake")
include("cmake/macros.cmake")

###############################################################################
### Options
###############################################################################

# User option to build configuration
# Default configuration is MPI parallel version of ADCIRC (PADCIRC)
set(ADCIRC_CONFIG PADCIRC CACHE STRING "Build configuration for ADCIRC")
set_property(CACHE ADCIRC_CONFIG PROPERTY STRINGS ADCIRC ADCSWAN SWAN PADCIRC PADCSWAN PUNSWAN)

if(ADCIRC_CONFIG STREQUAL "ADCIRC")
  message(STATUS "Build the serial version of ADCIRC")
  set(src_files ${adcirc_src_files})
elseif(ADCIRC_CONFIG STREQUAL "ADCSWAN")
  message(STATUS "Build the serial version of ADCIRC+SWAN")
  set(src_files ${adcswan1_src_files} ${adcswan2_src_files})
elseif(ADCIRC_CONFIG STREQUAL "SWAN")
  message(STATUS "Build the serial version of SWAN")
elseif(ADCIRC_CONFIG STREQUAL "PADCIRC")
  message(STATUS "Build the MPI parallel version of ADCIRC")
  set(src_files ${padcirc_src_files})
elseif(ADCIRC_CONFIG STREQUAL "PADCSWAN")
  message(STATUS "Build the MPI parallel version of ADCIRC+SWAN")
  set(src_files ${padcswan1_src_files} ${padcswan_src_files}) 
elseif(ADCIRC_CONFIG STREQUAL "PUNSWAN")
  message(STATUS "Build the MPI parallel version of ADCIRC+UNSTRUCTURED SWAN")
else()
  message(FATAL_ERROR "${ADCIRC_CONFIG} configuration is not supported! 
    Options: ADCIRC|ADCSWAN|SWAN|PADCIRC|PADCSWAN|PUNSWAN")
endif()

# User option for coupling
set(COUPLED OFF CACHE BOOL "Build NUOPC cap for coupling.")
if(COUPLED)
  if(ADCIRC_CONFIG STREQUAL "PADCIRC" OR ADCIRC_CONFIG STREQUAL "PADCSWAN" OR ADCIRC_CONFIG STREQUAL "PUNSWAN")
    set(src_files ${src_files} ${nuopc_src_files})
  endif()
endif()

# User option for 8-byte real numbers (or 4-byte if turned off)
set(PRECISION_8BYTE ON CACHE BOOL "Use 8-byte real numbers instead of 4-byte.")
if(PRECISION_8BYTE)
  set(PRECISION_FLAG "REAL8")
else(PRECISION_8BYTE)
  set(PRECISION_FLAG "REAL4")
endif(PRECISION_8BYTE)

# User option to build adcprep
if(MPI)
  set(BUILD_ADCPREP OFF CACHE BOOL "Build the MPI parallel ADCIRC preprocessor")
endif(MPI)

# User option to build utilities
set(BUILD_UTILITIES OFF CACHE BOOL "Build the ADCIRC utility programs")

###############################################################################
### Set variables
###############################################################################

# Machine specific
if(SUN)
  set(MACHINE_FLAG "CMACHSUN")
elseif(SGI)
  set(MACHINE_FLAG "SGI")
elseif(CRAY)
  set(MACHINE_FLAG "CRAY")
elseif(CRAYX1)
  set(MACHINE_FLAG "CRAYX1")
elseif(UNIX)
  set(MACHINE_FLAG "LINUX")
elseif(CYGWIN)
  set(MACHINE_FLAG "LINUX")
elseif(WIN32)
  set(MACHINE_FLAG "WINDOWS")
elseif(APPLE)
  set(MACHINE_FLAG "LINUX")
endif(SUN)

if(SUN_MACHINE)
  set(MACHINE_FLAG "${MACHINE_FLAG} -CMACHSUN")
endif(SUN_MACHINE)

if(VECTOR_COMPUTER)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} CVEC")
else(VECTOR_COMPUTER)
  set(ADCIRC_OPTION_FLAGS "${ADCIRC_OPTION_FLAGS} CSCA")
endif(VECTOR_COMPUTER)
string(STRIP ${ADCIRC_OPTION_FLAGS} ADCIRC_OPTION_FLAGS)

# SWAN specific
set(SWAN_FLAG "CSWAN")
set(PREP_SWAN_FLAG "ADCSWAN")

# MPI specific
set(ADCIRC_MPI_FLAG "CMPI")

# Combined 
set(LOCAL_COMPILER_DEFINITIONS "${ADCIRC_OPTION_FLAGS} ${PRECISION_FLAG} ${MACHINE_FLAG}")
separate_arguments(LOCAL_COMPILER_DEFINITIONS)

###############################################################################
### Find component specific dependencies
###############################################################################

if(ADCIRC_CONFIG STREQUAL "ADCSWAN" OR ADCIRC_CONFIG STREQUAL "SWAN" OR
   ADCIRC_CONFIG STREQUAL "PADCSWAN" OR ADCIRC_CONFIG STREQUAL "PUNSWAN") 
  find_package(Perl)
endif()

###############################################################################
### Libraries
###############################################################################

# Metis
if(BUILD_ADCPREP)
  find_package(METIS)
  if(NOT METIS_FOUND)
    add_library(libmetis STATIC ${metis_src_files})
    target_include_directories(libmetis PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod> $<INSTALL_INTERFACE:mod>)
    set_target_properties(libmetis PROPERTIES EXCLUDE_FROM_ALL TRUE)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      if(${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 10 OR ${CMAKE_C_COMPILER_VERSION} VERSION_EQUAL 10)
        set_target_properties(libmetis PROPERTIES COMPILE_FLAGS "-Wno-implicit-function-declaration")
      endif()
    endif()
  endif()
endif(BUILD_ADCPREP)

# Version
add_library(version OBJECT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/ADCIRC/scripts/adcirc_version.py)
  set(GEN_VER_CMD ./adcirc_version.py --create-version-file --directory ${CMAKE_CURRENT_SOURCE_DIR}/ADCIRC)
else()
  set(GEN_VER_CMD ./adcircVersion.sh ${CMAKE_CURRENT_SOURCE_DIR}/ADCIRC)
endif()
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
  COMMAND ${GEN_VER_CMD}
  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/ADCIRC/version.F
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/version_cmake.F
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ADCIRC/scripts
  COMMENT "Generating ADCIRC version...")
set_target_properties(version PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
set_target_properties(version PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}")
set_target_properties(version PROPERTIES EXCLUDE_FROM_ALL TRUE)

# Mkdir
add_library(mkdir OBJECT ${CMAKE_CURRENT_SOURCE_DIR}/ADCIRC/prep/mkdir.c)
target_include_directories(mkdir PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod> $<INSTALL_INTERFACE:mod>)
set_target_properties(mkdir PROPERTIES EXCLUDE_FROM_ALL TRUE)

###############################################################################
### Application 
###############################################################################

# ADCIRC
add_library(adcirc STATIC ${src_files})
set_target_properties(adcirc PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(adcirc PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod> $<INSTALL_INTERFACE:mod>)
target_include_directories(adcirc PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/mod)
set_target_properties(adcirc PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}")
target_compile_definitions(adcirc PRIVATE ${LOCAL_COMPILER_DEFINITIONS})
target_link_libraries(adcirc PUBLIC ESMF::ESMF $<TARGET_OBJECTS:version> $<TARGET_OBJECTS:mkdir>)
# mpi support
if(ADCIRC_CONFIG STREQUAL "PADCIRC" OR ADCIRC_CONFIG STREQUAL "PADCSWAN" OR ADCIRC_CONFIG STREQUAL "PUNSWAN")
  add_mpi(adcirc)
endif()
# netcdf support
add_netcdf(adcirc)

# SWAN
if(ADCIRC_CONFIG STREQUAL "PADCSWAN")
  add_library(swan1p STATIC "${swan1parallel_src_files}")
  add_library(swan2p STATIC "${swan2parallel_src_files}")
  add_custom_command(
    OUTPUT "${swan1parallel_src_files}" "${swan2parallel_src_files}"
    COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source
    COMMAND ${PERL} switch.pl -pun -adcirc -unix -outdir
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/swan_parallel_source *.ftn *.ftn90
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ADCIRC/thirdparty/swan
    COMMENT "Generating parallel unstructured SWAN sources...")
  set_target_properties(swan1p PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
  set_target_properties(swan2p PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
  add_library(padcswan1 STATIC ${padcswan1_src_files})
endif()

###############################################################################
### Install utilities 
###############################################################################

# Utilities
if(BUILD_UTILITIES)
  foreach(UTIL_FILE ${util_src_files})
    get_filename_component(UTIL_FILE_EXE ${UTIL_FILE} NAME_WE)
    add_executable(${UTIL_FILE_EXE} ${UTIL_FILE})
    set_target_properties(${UTIL_FILE_EXE} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
    install(TARGETS ${UTIL_FILE_EXE} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  endforeach()
endif(BUILD_UTILITIES)

# ADCPREP
if(BUILD_ADCPREP)
  # main
  add_executable(adcprep ${adcprep_src_files})
  set_target_properties(adcprep PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod2)
  set_target_properties(adcprep PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}")
  target_compile_definitions(adcprep PRIVATE ${LOCAL_COMPILER_DEFINITIONS})
  target_include_directories(adcprep PRIVATE prep)
  install(TARGETS adcprep RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  # mkdir, version and metis
  if(METIS_FOUND)
    target_link_libraries(adcprep PUBLIC metis)
    target_link_directories(adcprep PUBLIC ${METIS_LIBRARY_DIR})
    target_include_directories(adcprep PUBLIC ${METIS_INCLUDE_DIR})
  else()
    target_link_libraries(adcprep PUBLIC libmetis)
    add_dependencies(adcprep libmetis)
  endif()
  target_link_libraries(adcprep PUBLIC mkdir version)
  add_dependencies(adcprep mkdir version)

  # netcdf support
  add_netcdf(adcprep)
endif(BUILD_ADCPREP)

###############################################################################
### Install
###############################################################################

install(
  TARGETS adcirc
  EXPORT  adcirc-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/ADCIRC)

install(EXPORT      adcirc-config
        DESTINATION lib/cmake)
