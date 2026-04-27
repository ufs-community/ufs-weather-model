###############################################################################
### Fortran compiler flags
###############################################################################

message( STATUS "---------------------------------------------------------" )
message( STATUS "ROMS-UFS Configuration:" )
message( STATUS "---------------------------------------------------------" )

if( $ENV{FC} MATCHES "mpiifort" )
  set( COMM "intel" )
endif()

message( STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}" )
message( STATUS "COMM = ${COMM}" )

###############################################################################
### Set variables
###############################################################################

#------------------------------------------------------------------------------
# Specify ROMS specific variables

if( NOT ROMS_APP AND NOT ROMS_APP_DIR )
  message( STATUS "Using ROMS default application" )

  set( ROMS_APP IRENE )
  set( ROMS_APP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/APP/${ROMS_APP} )

  set( ROMS_UFS_APPLICATION     ${ROMS_APP} )
  set( ROMS_UFS_APPLICATION_DIR ${ROMS_APP_DIR} )

elseif( ROMS_APP AND NOT ROMS_APP_DIR )
  message( ERROR "If you set ROMS_APP you must also set ROMS_APP_DIR to
                  the absolute path containing your application configuration." )
else()
  set( ROMS_UFS_APPLICATION     ${ROMS_APP} )
  set( ROMS_UFS_APPLICATION_DIR ${ROMS_APP_DIR} )
endif()

message( STATUS "ROMS_UFS_APPLICATION = ${ROMS_UFS_APPLICATION}" )
message( STATUS "ROMS_UFS_APPLICATION_DIR = ${ROMS_UFS_APPLICATION_DIR}" )

# Allow the UFS and ROMS build systems to set ROMS source code and
# application header directories. ROMS build system always overwrites.

set( ROMS_SRC_DIR    ${CMAKE_CURRENT_SOURCE_DIR}/ROMS CACHE
     PATH "ROMS Source Directory" )
set( MY_HEADER_DIR   ${ROMS_UFS_APPLICATION_DIR} CACHE
     PATH "ROMS Application Header Directory" )

set( LIBTYPE STATIC )
set( ROMS_EXECUTABLE OFF )

#------------------------------------------------------------------------------
# Set netcdf variables

if( NetCDF_Fortran_FOUND )
  set( ENV{NETCDF_INCDIR} "${NetCDF_Fortran_INCLUDE_DIRS}" )
  set( ENV{NETCDF_LIBDIR} "${NetCDF_Fortran_LIBRARIES}" )
endif()

#------------------------------------------------------------------------------
# Set build directory

set( BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} )

###############################################################################
### Set configuration specific variables
###############################################################################

# ROMS CPP option to use UFS coupling.

add_definitions(-DCMEPS)
add_definitions(-DMPI)

message( STATUS "ROMS_SRC_DIR = ${ROMS_SRC_DIR}" )
message( STATUS "MY_HEADER_DIR = ${MY_HEADER_DIR}" )

# Add MY_CPP_FLAGS from the build script

if( MY_CPP_FLAGS )
  message( STATUS "ROMS Adding MY_CPP_FLAGS: ${MY_CPP_FLAGS}" )
  foreach( flag ${MY_CPP_FLAGS} )
    add_compile_definitions( ${flag} )
  endforeach()
endif()

###############################################################################
### Add include directories and ROMS project 
###############################################################################

if ( ESMF_FOUND )
  message ( STATUS "ESMF_F90COMPILEPATHS = ${ESMF_F90COMPILEPATHS}" )
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_F90COMPILEPATHS}")
endif()

add_subdirectory( ${ROMS_SRC_DIR} ${CMAKE_CURRENT_BINARY_DIR}/ROMS )

target_include_directories(ROMS_static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/ROMS/modules>
                                        $<INSTALL_INTERFACE:mod>)

###############################################################################
### Install
###############################################################################

install(
  TARGETS ROMS_static 
  EXPORT  roms-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library )

install( DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/ROMS/modules DESTINATION ${CMAKE_INSTALL_PREFIX}/ROMS )

install( EXPORT      roms-config
         DESTINATION lib/cmake )
