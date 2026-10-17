cmake_minimum_required(VERSION 3.15)
project(ufs_share LANGUAGES Fortran)

# ----------------------------------------------------------------------------
# 1. Initialization & Dependencies
# ----------------------------------------------------------------------------
include(CTest)

if(DEFINED ENV{ESMFMKFILE})
  get_filename_component(TEMP_DIR $ENV{ESMFMKFILE} DIRECTORY)
  get_filename_component(ESMF_ROOT ${TEMP_DIR} DIRECTORY)
  list(APPEND CMAKE_MODULE_PATH ${ESMF_ROOT}/cmake)
endif()
find_package(ESMF REQUIRED)

# ----------------------------------------------------------------------------
# 2. Build the Production Library
# ----------------------------------------------------------------------------
add_library(ufs_share qc_coupling_mod.F90)

# PUBLIC linking ensures anything that links to qc_coupling also gets ESMF
target_link_libraries(ufs_share PUBLIC ESMF::ESMF)

# Put .mod files in the current binary dir
set_target_properties(ufs_share PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# share the .mod directory with downstream targets
target_include_directories(ufs_share PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# ----------------------------------------------------------------------------
# 3. Conditionally add unit tests (default is ON)
# ----------------------------------------------------------------------------
if(BUILD_TESTING)
    add_subdirectory(unit_tests)
endif()
