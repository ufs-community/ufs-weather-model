#------------------------------------------------------------------------------
# Add source files
# NUOPC cap
list(APPEND _noahmp_cap_files noahmp/drivers/nuopc/lnd_comp_kind.F90
                              noahmp/drivers/nuopc/lnd_comp_types.F90
                              noahmp/drivers/nuopc/lnd_comp_shr.F90
                              noahmp/drivers/nuopc/lnd_comp_io.F90
                              noahmp/drivers/nuopc/lnd_comp_domain.F90
                              noahmp/drivers/nuopc/lnd_comp_import_export.F90
                              noahmp/drivers/nuopc/lnd_comp_nuopc.F90
                              noahmp/drivers/nuopc/lnd_comp_driver.F90)

# CCPP interface 
list(APPEND _noahmp_ccpp_files noahmp/drivers/ccpp/noahmpdrv.F90
                               noahmp/drivers/ccpp/sfc_diff.f
                               noahmp/drivers/ccpp/machine.F
                               noahmp/drivers/ccpp/noahmp_tables.f90
                               noahmp/drivers/ccpp/namelist_soilveg.f
                               noahmp/drivers/ccpp/set_soilveg.f
                               noahmp/drivers/ccpp/funcphys.f90
                               noahmp/drivers/ccpp/physcons.F90)

# NoahMP
list(APPEND _noahmp_files noahmp/src/module_sf_noahmplsm.F90
                          noahmp/src/module_sf_noahmp_glacier.F90)

#------------------------------------------------------------------------------
# Set CCPP flags for C/C++/Fortran preprocessor 
add_definitions(-DCCPP)

#------------------------------------------------------------------------------
# Set MPI flags for C/C++/Fortran preprocessor
if(MPI)
  add_definitions(-DMPI)
endif()

#------------------------------------------------------------------------------
# Set flag for 32bit dynamics build
if(32BIT)
  message(STATUS "Compile CCPP slow physics with 64-bit precision, fast physics with 32-bit precision")
  add_definitions(-DOVERLOAD_R4)
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-real-size 64 -no-prec-div -no-prec-sqrt")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_PHYSICS  "-fdefault-real-8 -fdefault-double-8")
  endif()
  set(CMAKE_Fortran_FLAGS_DYNAMICS "")
else()
  message(STATUS "Compile CCPP physics with 64-bit precision")
  remove_definitions(-DOVERLOAD_R8)
  remove_definitions(-DOVERLOAD_R4)
  set(CMAKE_Fortran_FLAGS_PHYSICS  "")
  set(CMAKE_Fortran_FLAGS_DYNAMICS "")
endif()

#------------------------------------------------------------------------------
# Add model-specific flags for C/C++/Fortran preprocessor
add_definitions(-DMOIST_CAPPA -DUSE_COND -DNEMS_GSM)
add_definitions(-DINTERNAL_FILE_NML)

#------------------------------------------------------------------------------
# NOAHMP
add_library(noahmp STATIC ${_noahmp_cap_files} ${_noahmp_ccpp_files} ${_noahmp_files})
set_target_properties(noahmp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(noahmp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(noahmp PUBLIC esmf)

###############################################################################
### Install
###############################################################################

install(
  TARGETS noahmp
  EXPORT  noahmp-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/NOAHMP)

install(EXPORT      noahmp-config
        DESTINATION lib/cmake)
