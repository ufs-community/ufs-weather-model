# Set default project to unknown
if(NOT PROJECT)
  message(STATUS "Setting CCPP project to 'unknown' as none was specified.")
  set(PROJECT "Unknown")
endif (NOT PROJECT)

#------------------------------------------------------------------------------
# Set the sources
set(SOURCES_C
    ccpp_dl.h
    ccpp_dl.c
    ccpp_fields_idx.h
    ccpp_fields_idx.c
    ccpp_memory.h
    ccpp_memory.c
    ccpp_utils.h
    ccpp_utils.c
    ccpp_xml.h
    ccpp_xml.c
)
set(SOURCES_F90
    ccpp.F90
    ccpp_dl.F90
    ccpp_errors.F90
    ccpp_fcall.F90
    ccpp_fields.F90
    ccpp_memory.F90
    ccpp_strings.F90
    ccpp_suite.F90
    ccpp_types.F90
    ccpp_xml.F90
    ccpp_api.F90
)
# Generate list of Fortran modules from defined sources
foreach(source_f90 ${SOURCES_F90})
    string(REGEX REPLACE ".F90" ".mod" module_f90 ${source_f90})
    list(APPEND MODULES_F90 ${CMAKE_CURRENT_BINARY_DIR}/${module_f90})
endforeach()

#------------------------------------------------------------------------------
# Find libXML2
find_package(LibXml2 REQUIRED)
if (LIBXML2_FOUND)
    include_directories(${LIBXML2_INCLUDE_DIR})
    list(APPEND LIBS ${LIBXML2_LIBRARIES})
endif(LIBXML2_FOUND)

#------------------------------------------------------------------------------
# CMake Modules
# Set the CMake module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

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

#------------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                 "MinSizeRel" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------
# Pass debug/release flag to Fortran files for preprocessor
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDEBUG")
else(CMAKE_BUILD_TYPE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DRELEASE")
endif()

#------------------------------------------------------------------------------
# Add the toplevel source directory to our include directoies (for .h)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add the toplevel binary directory to our include directoies (for .mod)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Set a cached variable containing the includes, so schemes can use them
set(${PACKAGE}_INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}$<SEMICOLON>${CMAKE_CURRENT_BINARY_DIR}"
     CACHE FILEPATH "${PACKAGE} include directories")
set(${PACKAGE}_LIB_DIRS
    "${CMAKE_CURRENT_BINARY_DIR}"
     CACHE FILEPATH "${PACKAGE} library directories")

#------------------------------------------------------------------------------
# Add the tests
add_subdirectory(tests)

#------------------------------------------------------------------------------
# Define the executable and what to link
add_library(ccpp ${SOURCES_C} ${SOURCES_F90})
target_link_libraries(ccpp LINK_PUBLIC ${LIBS} ${CMAKE_DL_LIBS})
set_target_properties(ccpp PROPERTIES VERSION ${PROJECT_VERSION}
                                      SOVERSION ${PROJECT_VERSION_MAJOR}
                                      LINK_FLAGS ${CMAKE_Fortran_FLAGS})

#------------------------------------------------------------------------------
# Installation
#

# Find all the C headers and Fortran modules
file(GLOB HEADERS_C
     "${CMAKE_CURRENT_SOURCE_DIR}/ccpp*.h"
)

if (PROJECT STREQUAL "CCPP-FV3")
  target_include_directories(ccpp PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include>
  )
else (PROJECT STREQUAL "CCPP-FV3")
  target_include_directories(ccpp PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
      $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
  )
endif (PROJECT STREQUAL "CCPP-FV3")

# Define where to install the library
install(TARGETS ccpp
        EXPORT ccpp-targets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION lib
)

# Export our configuration
install(EXPORT ccpp-targets
        FILE ccpp-config.cmake
        DESTINATION lib/cmake
)

if (PROJECT STREQUAL "CCPP-FV3")
  # Define where to install the C headers and Fortran modules
  install(FILES ${HEADERS_C} DESTINATION include)
  install(FILES ${MODULES_F90} DESTINATION include)
else (PROJECT STREQUAL "CCPP-FV3")
  # Define where to install the C headers and Fortran modules
  install(FILES ${HEADERS_C} DESTINATION include/${PROJECT_NAME})
  install(FILES ${MODULES_F90} DESTINATION include/${PROJECT_NAME})
endif (PROJECT STREQUAL "CCPP-FV3")
