# Set default project to unknown
if(NOT PROJECT)
  message(STATUS "Setting CCPP project to 'unknown' as none was specified.")
  set(PROJECT "Unknown")
endif (NOT PROJECT)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

#------------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.11)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
    project(check VERSION 0.0.1)
else(POLICY CMP0048)
    project(check)
    set(PROJECT_VERSION 0.0.1)
    set(PROJECT_VERSION_MAJOR 0)
    set(PROJECT_VERSION_MINOR 0)
    set(PROJECT_VERSION_PATCH 1)
endif(POLICY CMP0048)

#------------------------------------------------------------------------------
set(PACKAGE "check")
set(AUTHORS  "Timothy Brown" "Dom Heinzeller")
string(TIMESTAMP YEAR "%Y")

#------------------------------------------------------------------------------
# Enable Fortran
enable_language(Fortran)

#------------------------------------------------------------------------------
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake")

#------------------------------------------------------------------------------
# By default we want a shared library
option(BUILD_SHARED_LIBS "Build a shared library" ON)

#------------------------------------------------------------------------------
# Add the CCPP include/module directory and libraries, currently depends on build
# see FV3_current_trunk/ccpp/CMakeLists.txt on how to set CCPP_INCLUDE_DIRS etc.
if (PROJECT STREQUAL "CCPP-FV3")
  # Add the CCPP include/module directory
  set(CCPP_INCLUDE_DIRS "" CACHE FILEPATH "Path to ccpp includes")
  set_property(DIRECTORY PROPERTY INCLUDE_DIRECTORIES ${CCPP_INCLUDE_DIRS})
  # Add the CCPP library
  set(CCPP_LIB_DIRS "" CACHE FILEPATH "Path to ccpp library")
  link_directories(${CCPP_LIB_DIRS})
  list(APPEND LIBS "ccpp")
else (PROJECT STREQUAL "CCPP-FV3")
  # Add the CCPP include/module directory
  INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}/../../../../src")
  # Add the CCPP library
  LINK_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}/../../../../src")
  list(APPEND LIBS "ccpp")
endif (PROJECT STREQUAL "CCPP-FV3")

#------------------------------------------------------------------------------
# Set the sources
set(SOURCES
    check_test.f90
    check_noop.f90
)

# DH the mechanism to create caps from XML is no longer supported; todo:
# replace with a test that parses a simple metadata table in check_test.f90
# and creates check_test_cap.f90 - also needs a variable definition file.
##------------------------------------------------------------------------------
## Add the auto-generated caps
#add_custom_command(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/check_test_cap.f90
#    DEPENDS ${CCPP_MKCAP}
#    COMMAND ${CCPP_MKCAP} -o ${CMAKE_CURRENT_BINARY_DIR}/check_test_cap.f90
#                             ${CMAKE_CURRENT_SOURCE_DIR}/scheme.xml
#)
#list(APPEND SOURCES ${CMAKE_CURRENT_BINARY_DIR}/check_test_cap.f90)

#------------------------------------------------------------------------------
# The Fortran compiler/linker flag inserted by cmake to create shared libraries
# with the Intel compiler is deprecated (-i_dynamic), correct here.
# CMAKE_Fortran_COMPILER_ID = {"Intel", "PGI", "GNU", "Clang", "MSVC", ...}
if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
    string(REPLACE "-i_dynamic" "-shared-intel"
           CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
           "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
endif()

# Guard for undefined/empty CMAKE_Fortran_FLAGS
set(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS}")

add_library(check ${SOURCES})
target_link_libraries(check LINK_PUBLIC ${LIBS})
set_target_properties(check PROPERTIES VERSION ${PROJECT_VERSION}
                                       SOVERSION ${PROJECT_VERSION_MAJOR}
                                       COMPILE_FLAGS ${CMAKE_Fortran_FLAGS}
                                       LINK_FLAGS ${CMAKE_Fortran_FLAGS})
