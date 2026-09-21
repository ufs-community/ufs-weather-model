# =============================================================================
# MOM6 Cap OutputLog Unit Tests Build Configuration
# =============================================================================
if(NOT PROJECT_NAME)
  cmake_minimum_required(VERSION 3.19)
  project(MOM6_UnitTests LANGUAGES Fortran)
  enable_testing()
  find_package(ESMF REQUIRED)
endif()

add_library(mom6_test_utils STATIC test_utils.F90)
target_link_libraries(mom6_test_utils PRIVATE mom6 ESMF::ESMF)
if(MOM6_INC_DIR)
  target_include_directories(mom6_test_utils PRIVATE ${MOM6_INC_DIR})
endif()

add_library(nc_fixture_mod STATIC nc_fixture_mod.F90)
target_link_libraries(nc_fixture_mod PRIVATE mom6 ESMF::ESMF NetCDF::NetCDF_Fortran MPI::MPI_Fortran)
if(MOM6_INC_DIR)
  target_include_directories(nc_fixture_mod PRIVATE ${MOM6_INC_DIR})
endif()

add_library(mom6_test_helpers STATIC test_helpers.F90)
target_link_libraries(mom6_test_helpers PRIVATE mom6 ESMF::ESMF mom6_test_utils nc_fixture_mod)
if(MOM6_INC_DIR)
  target_include_directories(mom6_test_helpers PRIVATE ${MOM6_INC_DIR})
endif()
# -----------------------------------------------------------------------------
# Registers the executable, its common link libraries / include dir /
# UNIT_TESTING compile definition, and its ctest entry in one call, so the
# four test targets below don't each repeat the same five CMake calls.
# -----------------------------------------------------------------------------
function(add_mom6_unit_test exe_name)
  cmake_parse_arguments(ARG "USE_MPI" "CTEST_NAME;NPROCS" "SOURCES;EXTRA_LIBS" ${ARGN})

  if(NOT ARG_CTEST_NAME)
    set(ARG_CTEST_NAME ${exe_name})
  endif()
  if(NOT ARG_NPROCS)
    set(ARG_NPROCS 2)
  endif()

  add_executable(${exe_name} ${ARG_SOURCES})
  target_link_libraries(${exe_name} PRIVATE mom6 ESMF::ESMF ${ARG_EXTRA_LIBS})
  target_compile_definitions(${exe_name} PRIVATE UNIT_TESTING)
  if(MOM6_INC_DIR)
    target_include_directories(${exe_name} PRIVATE ${MOM6_INC_DIR})
  endif()

  if(ARG_USE_MPI)
    add_test(NAME ${ARG_CTEST_NAME}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ARG_NPROCS} $<TARGET_FILE:${exe_name}>)
  else()
    add_test(NAME ${ARG_CTEST_NAME} COMMAND ${exe_name})
  endif()

  # Umbrella target: every test executable registered through this function
  # becomes a dependency automatically.
  if(NOT TARGET mom6_unit_tests)
    add_custom_target(mom6_unit_tests)
  endif()
  add_dependencies(mom6_unit_tests ${exe_name})
endfunction()

# readnml test
add_mom6_unit_test(test_outputlog_readnml
  CTEST_NAME MOM6_Outputlog_Readnml
  SOURCES test_outputlog_readnml.F90
  EXTRA_LIBS mom6_test_utils)

# alarminit test
add_mom6_unit_test(test_outputlog_alarm
  CTEST_NAME MOM6_Outputlog_Alarm
  SOURCES test_outputlog_alarminit.F90
  EXTRA_LIBS mom6_test_utils)

# fixture completion test
add_mom6_unit_test(test_outputlog_completion
 CTEST_NAME MOM6_Outputlog_Completion
 SOURCES test_outputlog_completion.F90
 EXTRA_LIBS nc_fixture_mod mom6_test_helpers NetCDF::NetCDF_Fortran MPI::MPI_Fortran
 USE_MPI)

# track_freqn test
add_mom6_unit_test(test_outputlog_freqn
   CTEST_NAME MOM6_Outputlog_Freqn
   SOURCES test_outputlog_freqn.F90
   EXTRA_LIBS mom6_test_utils nc_fixture_mod mom6_test_helpers NetCDF::NetCDF_Fortran MPI::MPI_Fortran
   USE_MPI)

# track_restn test
add_mom6_unit_test(test_outputlog_restn
  CTEST_NAME MOM6_Outputlog_Restn
  SOURCES test_outputlog_restn.F90
  EXTRA_LIBS mom6_test_utils mom6_test_helpers NetCDF::NetCDF_Fortran MPI::MPI_Fortran
  USE_MPI)

if(POLICY CMP0144)
    cmake_policy(SET CMP0144 NEW)
endif()

find_package(Doxygen)

if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/../docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(feature_docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating NUOPC cap and unit-test documentation"
        VERBATIM)
endif()
