###############################################################################
### Fortran compiler flags
###############################################################################

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Waliasing -fcray-pointer -fconvert=big-endian -fno-range-check -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-512")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
  set(C_PREPROCESS_FLAG "-cpp")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -sox")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  set(CMAKE_Fortran_LINK_FLAGS "")
  set(C_PREPROCESS_FLAG "-cpp")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

###############################################################################
### Add module directory
###############################################################################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/SCHISM-interface/SCHISM/cmake/modules)

###############################################################################
### Options
###############################################################################

option(USE_CICE "Build SCHISM with sea-ice coupling support" OFF)

###############################################################################
### Set variables
###############################################################################

set(BLD_STANDALONE "OFF")

if(USE_CICE)
  add_definitions(-DUSE_CICE)
endif()

if(NOT NO_PARMETIS)
  # find parmetis
  find_package(ParMETIS) 
  if(ParMETIS_FOUND)
    set(PARMETIS_DIR $ENV{PARMETIS_ROOT})
    set(PARMETIS_VER "ParMetis-${PARMETIS_VERSION}")
    # find metis
    find_package(METIS)
    message(STATUS "METIS_LIBRARY_DIR: ${METIS_LIBRARY_DIR}")
    message(STATUS "METIS_INCLUDE_DIR: ${METIS_INCLUDE_DIR}")
  else()
    set(PARMETIS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/SCHISM/src/ParMetis-4.0.3")
    set(PARMETIS_VER ParMetis-4.0.3 CACHE STRING "Version and relative path to Parmetis library inside /src.") 
    add_subdirectory(SCHISM/src/${PARMETIS_VER})
  endif()
  message(STATUS "PARMETIS_DIR: ${PARMETIS_DIR}")
  message(STATUS "PARMETIS_VER: ${PARMETIS_VER}")
endif()

if(NetCDF_Fortran_FOUND)
  get_filename_component(_root ${NetCDF_Fortran_INCLUDE_DIR}/.. ABSOLUTE)
  set(NetCDF_FORTRAN_DIR ${_root})
  include_directories(${_root}/include)
endif()

if(NetCDF_C_FOUND)
  get_filename_component(_root ${NetCDF_C_INCLUDE_DIR}/.. ABSOLUTE)
  set(NetCDF_C_DIR ${_root})
  include_directories(${_root}/include)
endif()

set(NetCDFLIBS "${NetCDF_Fortran_LIBRARIES} ${NetCDF_C_LIBRARIES}")

message(STATUS "NetCDF_FORTRAN_DIR: ${NetCDF_FORTRAN_DIR}")
message(STATUS "NetCDF_C_DIR: ${NetCDF_C_DIR}")
message(STATUS "NetCDFLIBS: ${NetCDFLIBS}")

if(NetCDF_PARALLEL)
  set(HDF5_LIBRARIES ${NetCDF_C_LIBRARIES})
endif() 

###############################################################################
### File lists and macros
###############################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/SCHISM/cmake/modules")
add_subdirectory(SCHISM/src)

list(APPEND schism_nuopc_src_files
  SCHISM-ESMF/src/schism/schism_nuopc_util.F90
  SCHISM-ESMF/src/schism/schism_nuopc_cap.F90
  SCHISM-ESMF/src/schism/schism_esmf_util.F90
  SCHISM-ESMF/src/schism/schism_esmf_cap.F90
  SCHISM-ESMF/src/schism/schism_bmi.F90
)

###############################################################################
### Application
###############################################################################

# SCHISM
if(USE_WWM)
  add_library(schism STATIC ${schism_nuopc_src_files} $<TARGET_OBJECTS:core> $<TARGET_OBJECTS:hydro> $<TARGET_OBJECTS:wwmIII>)
else(USE_WWM)
  add_library(schism STATIC ${schism_nuopc_src_files} $<TARGET_OBJECTS:core> $<TARGET_OBJECTS:hydro>)
endif(USE_WWM)
set_target_properties(schism PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(schism PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod> $<INSTALL_INTERFACE:mod>)
target_link_libraries(schism PUBLIC ESMF::ESMF)
if(NOT NO_PARMETIS)
  list(APPEND schism_libs parmetis metis)
endif()
target_link_libraries(schism PUBLIC ${schism_libs})
if(ParMETIS_FOUND)
  target_link_directories(schism PUBLIC ${PARMETIS_DIR}/lib)
endif()
if(METIS_FOUND)
  message(INFO ${METIS_LIBRARY_DIR})
  target_link_directories(schism PUBLIC ${METIS_LIBRARY_DIR})
  target_include_directories(schism PUBLIC ${METIS_INCLUDE_DIR})
endif()

# Sub-modules
list(APPEND schism_modules
  core
  hydro) 

foreach(ITEM ${schism_modules})
  set_target_properties(${ITEM} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
  target_include_directories(${ITEM} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
endforeach(ITEM)

###############################################################################
### Install
###############################################################################

list(APPEND schism_targets schism)
if(NOT NO_PARMETIS)
  if (NOT ParMETIS_FOUND)
    list(APPEND schism_targets metis parmetis)
  endif()
endif()

install(
  TARGETS ${schism_targets}
  EXPORT  schism-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
  DESTINATION ${CMAKE_INSTALL_PREFIX}/SCHISM)

install(
  EXPORT schism-config
  DESTINATION lib/cmake)
