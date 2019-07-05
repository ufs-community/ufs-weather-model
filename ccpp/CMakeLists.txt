cmake_minimum_required(VERSION 2.8.11)

# Use rpaths on MacOSX
set(CMAKE_MACOSX_RPATH 1)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif(POLICY CMP0048)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif(POLICY CMP0042)

PROJECT(CCPP-FV3)
set(PROJECT "CCPP-FV3")

ENABLE_LANGUAGE(Fortran)

include(CMakeForceCompiler)

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
# Generate Compiler flags for C/CXX/Fortran - set to match NEMSfv3gfs flags
# for DEBUG, REPRO and PROD mode for all three compilers

# Set basic compiler optimization levels, depending on build type
if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  ADD_DEFINITIONS(-DDEBUG)
  set (CMAKE_C_FLAGS_DEBUG       "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_CXX_FLAGS_DEBUG     "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fPIC" CACHE STRING "" FORCE)
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS_DEBUG}")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS_DEBUG}")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS_DEBUG}")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_C_FLAGS_RELEASE       "-O1 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_CXX_FLAGS_RELEASE     "-O1 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-O1 -fPIC" CACHE STRING "" FORCE)
  else (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_CXX_FLAGS_RELEASE     "-O2 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS_RELEASE}")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS_RELEASE}")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS_RELEASE}")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Bitforbit")
  if (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_C_FLAGS_RELEASE       "-O1 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_CXX_FLAGS_RELEASE     "-O1 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-O1 -fPIC" CACHE STRING "" FORCE)
  else (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
    set (CMAKE_C_FLAGS_RELEASE       "-O2 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_CXX_FLAGS_RELEASE     "-O2 -fPIC" CACHE STRING "" FORCE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fPIC" CACHE STRING "" FORCE)
  endif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set (CMAKE_C_FLAGS        "${CMAKE_C_FLAGS_RELEASE}")
  set (CMAKE_CXX_FLAGS      "${CMAKE_CXX_FLAGS_RELEASE}")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS_RELEASE}")
endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")

# Set compiler-specific options, depending on build type; some optimization flags may only be used
# for certain files (controlled for example in ccpp-physics' CMakeLists.txt). These must be set as
# CMAKE_Fortran_FLAGS_OPT, which by default is empty. Flags for controlling the default precision
# are set as CMAKE_Fortran_FLAGS_DEFAULT_PREC.
set (CMAKE_Fortran_FLAGS_OPT "")
if (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fcray-pointer -ffree-line-length-none -fno-range-check")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fbacktrace -cpp")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -ffpe-trap=invalid,zero,overflow -fcheck=bounds -fbacktrace -fno-range-check")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (CMAKE_Fortran_FLAGS_DEFAULT_PREC "-fdefault-real-8 -fdefault-double-8")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mcray=pointer -Kieee -byteswapio -Mfprelaxed=div -Mfprelaxed=sqrt")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -traceback")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -Mpreprocess")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -tp=x64 -Ktrap=fp -Mbounds")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Bitforbit")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -tp=x64")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -tp=haswell")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (CMAKE_Fortran_FLAGS_DEFAULT_PREC "-i4 -r8")
elseif (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -nowarn -sox -align array64byte")
  set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align -fpp -save-temps")
  if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ftrapuv -traceback")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftrapuv -traceback")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fp-stack-check")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -fstack-protector-all -fpe0 -traceback -debug -ftrapuv")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Bitforbit")
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -debug minimal -fp-model consistent -qoverride-limits -g -traceback")
  elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    # Specify aggressive optimization flags (to be overwritten for individual files in ccpp-physics' CMakeLists.txt)
    if (SIMDMULTIARCH)
      set (CMAKE_Fortran_FLAGS_OPT "-no-prec-div -no-prec-sqrt -axSSE4.2,AVX,CORE-AVX2,CORE-AVX512")
    else (SIMDMULTIARCH)
      set (CMAKE_Fortran_FLAGS_OPT "-no-prec-div -no-prec-sqrt -xCORE-AVX2")
    endif (SIMDMULTIARCH)
    set (CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -debug minimal -fp-model source -qoverride-limits -qopt-prefetch=3")
  endif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
  set (CMAKE_Fortran_FLAGS_DEFAULT_PREC "-i4 -real-size 64")
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
# Detect OS and set preprocessor flags
if(APPLE)
  ADD_DEFINITIONS(-DMACOSX)
elseif(UNIX)
  ADD_DEFINITIONS(-DLINUX)
else (APPLE)
  message (FATAL_ERROR "Unsupported platform, only Linux and MacOSX are supported at this time.")
endif(APPLE)

#------------------------------------------------------------------------------
# Add -DCCPP preprocessor flag (needed to preprocess GFS_typedefs.F90 from FV3)
ADD_DEFINITIONS(-DCCPP)

#------------------------------------------------------------------------------
# Add host-model specific preprocessor flag (needed for some physics schemes)
ADD_DEFINITIONS(-DFV3)

#------------------------------------------------------------------------------
# Set MPI flags for C/C++/Fortran preprocessor
if (MPI)
  ADD_DEFINITIONS(-DMPI)
endif (MPI)

#------------------------------------------------------------------------------
# Set flags for static building/linking of CCPP
if (STATIC)
  ADD_DEFINITIONS(-DSTATIC)
endif (STATIC)

#------------------------------------------------------------------------------
# Set flags for CCPP transition (selective reduction of optimization, Intel15)
if (TRANSITION)
  ADD_DEFINITIONS(-DTRANSITION)
endif (TRANSITION)

#------------------------------------------------------------------------------
# Set flag for 32bit dynamics build
if (DYN32)
  message(STATUS "Dynamics compiled with 32-bit option, adjust fv_sat_adj types")
  ADD_DEFINITIONS(-DOVERLOAD_R4)
endif (DYN32)

#------------------------------------------------------------------------------
# If a dynamic build is requested, set TARGET_SUPPORTS_SHARED_LIBS to TRUE; on
# some systems, cmake doesn't set this even though the system supports shared
# libraries (e.g. Cray systems)
if(NOT STATIC)
SET_PROPERTY(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
endif(NOT STATIC)

#------------------------------------------------------------------------------
# Add model-specific flags for C/C++/Fortran preprocessor
ADD_DEFINITIONS(-DMOIST_CAPPA -DUSE_COND -DNEMS_GSM)
ADD_DEFINITIONS(-DINTERNAL_FILE_NML)

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
# Set SIONlib flags for preprocessor, compiler and linker (if defined)
if (SIONLIB)
  execute_process(COMMAND ${SIONLIB}/bin/sionconfig --mpi --cflags --f90
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE SIONlib_Fortran_FLAGS
                  ERROR_VARIABLE stderr
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (NOT exit_code EQUAL "0")
     message (FATAL_ERROR "Execution of '${SIONLIB}/bin/sionconfig --mpi --cflags --f90' failed")
  endif (NOT exit_code EQUAL "0")
  execute_process(COMMAND ${SIONLIB}/bin/sionconfig --mpi --libs --f90
                  RESULT_VARIABLE exit_code
                  OUTPUT_VARIABLE SIONlib_Fortran_LIBS
                  ERROR_VARIABLE stderr
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if (NOT exit_code EQUAL "0")
     message (FATAL_ERROR "Execution of '${SIONLIB}/bin/sionconfig --mpi --libs --f90' failed")
  endif (NOT exit_code EQUAL "0")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${SIONlib_Fortran_FLAGS} ${SIONlib_Fortran_LIBS}")
  ADD_DEFINITIONS(-DSION)
  message (STATUS "Enable SIONlib support")
else (SIONLIB)
  message (STATUS "Disable SIONlib support")
endif (SIONLIB)

#------------------------------------------------------------------------------
# Set Intel MKL flags for preprocessor, compiler and linker (if defined)
if(MKL_DIR)
  set (MKL_INC "-m64 -I${MKL_DIR}/include")
  set (MKL_LIB "-L${MKL_DIR}/lib -Wl,-rpath,${MKL_DIR}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_INC} ${MKL_LIB}")
  ADD_DEFINITIONS(-DMKL)
  message (STATUS "Enable Intel MKL support")
else(MKL_DIR)
  message (STATUS "Disable Intel MKL support")
endif(MKL_DIR)

#------------------------------------------------------------------------------
# Set netCDF flags for preprocessor, compiler and linker (if defined)
if(NETCDF_DIR)
  set (NETCDF_INC "-I${NETCDF_DIR}/include")
  set (NETCDF_LIB "-L${NETCDF_DIR}/lib -lnetcdff -lnetcdf")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_INC} ${NETCDF_LIB}")
  ADD_DEFINITIONS(-DNETCDF)
  message (STATUS "Enable netCDF support")
else(NETCDF_DIR)
  message (STATUS "Disable netCDF support")
endif(NETCDF_DIR)

#------------------------------------------------------------------------------
# Set ESMF flags for preprocessor, compiler and linker (if defined)
if(ESMF_LIB_DIR)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -L${ESMF_LIB_DIR} -lesmf")
  message (STATUS "Enable ESMF support")
else(ESMF_LIB_DIR)
  message (STATUS "Disable ESMF support")
endif(ESMF_LIB_DIR)

#------------------------------------------------------------------------------
# Set multi-gases flags for preprocessor (if defined)
if(MULTI_GASES)
  message (STATUS "Enable multi-gases physics")
  ADD_DEFINITIONS(-DMULTI_GASES)
else(MULTI_GASES)
  message (STATUS "Disable multi-gases physics")
endif(MULTI_GASES)

#------------------------------------------------------------------------------
# Build CCPP framework and physics

SET(CCPP_FRAMEWORK_SRC ${CMAKE_SOURCE_DIR}/framework)
ADD_SUBDIRECTORY(${CCPP_FRAMEWORK_SRC} ${CMAKE_BINARY_DIR}/framework)

# Add directory containing compiled CCPP Fortran modules to included directories for physics
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I${CMAKE_BINARY_DIR}/framework/src")

SET(CCPP_PHYSICS_SRC ${CMAKE_SOURCE_DIR}/physics)
ADD_SUBDIRECTORY(${CCPP_PHYSICS_SRC} ${CMAKE_BINARY_DIR}/physics)
