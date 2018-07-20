cmake_minimum_required(VERSION 2.8.11)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

PROJECT(CCPP-FV3)
set(PROJECT "CCPP-FV3")

ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Coverage")
endif()

#------------------------------------------------------------------------------
# Set include directories for make, also set CCPP_... for external projects
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/framework/src)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/physics)
SET(CCPP_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/framework/src;${CMAKE_BINARY_DIR}/physics")
SET(CCPP_LIB_DIRS "${CMAKE_BINARY_DIR}/framework/src;${CMAKE_BINARY_DIR}/physics")

#------------------------------------------------------------------------------
# Generate Compiler flags for C/CXX/Fortran
set (CMAKE_C_FLAGS       "-fPIC")
set (CMAKE_CXX_FLAGS     "-fPIC")
set (CMAKE_Fortran_FLAGS "-fPIC")

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  ADD_DEFINITIONS(-DDEBUG)
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS} -O0 -g")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS} -O0 -g")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -O0 -g")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS} -O3")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS} -O3")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -O3")
endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")

# Set compiler options
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fbacktrace -ffree-line-length-none -cpp")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fcheck=all -ffpe-trap=zero -finit-real=nan")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fno-range-check -fbacktrace")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mfree -traceback -Mpreprocess")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mipa=noconst -Mfptrap -Mbounds")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mfprelaxed=div -Mfprelaxed=sqrt")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -align array64byte -traceback -qno-opt-dynamic-align -fpp")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fpe0 -debug extended -ftrapuv -warn")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ftz -nowarn -debug minimal")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
else (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  message (FATAL_ERROR "CCPP has only been compiled with gfortran, pgf90 and ifort.")
endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

#------------------------------------------------------------------------------
# The Fortran compiler/linker flag inserted by cmake to create shared libraries
# with the Intel compiler is deprecated (-i_dynamic), correct here.
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
  string(REPLACE "-i_dynamic" "-shared-intel"
         CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS
         "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS}")
  string(REPLACE "-i_dynamic" "-shared-intel"
         CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS
         "${CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS}")
endif()

#------------------------------------------------------------------------------
# CMake Modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/framework/cmake")

#------------------------------------------------------------------------------
# Set MPI flags for C/C++/Fortran preprocessor
if (MPI)
  ADD_DEFINITIONS(-DMPI)
endif (MPI)

#------------------------------------------------------------------------------
# Set OpenMP flags for C/C++/Fortran compiler and preprocessor
if (OPENMP)
  include(detect_openmp)
  detect_openmp()
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
  message(STATUS "Enable OpenMP support for C/C++/Fortran compiler")
  ADD_DEFINITIONS(-DOPENMP)
else (OPENMP)
  message (STATUS "Disable OpenMP support for C/C++/Fortran compiler")
endif (OPENMP)

#------------------------------------------------------------------------------
# Build CCPP framework and physics

SET(CCPP_FRAMEWORK_SRC ${CMAKE_SOURCE_DIR}/framework)
ADD_SUBDIRECTORY(${CCPP_FRAMEWORK_SRC} ${CMAKE_BINARY_DIR}/framework)

# Add directory containing compiled CCPP Fortran modules to included directories for physics
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${CMAKE_BINARY_DIR}/framework/src")

SET(CCPP_PHYSICS_SRC ${CMAKE_SOURCE_DIR}/physics)
ADD_SUBDIRECTORY(${CCPP_PHYSICS_SRC} ${CMAKE_BINARY_DIR}/physics)
