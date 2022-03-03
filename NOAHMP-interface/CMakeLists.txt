# Compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fdefault-real-8 -fdefault-double-8")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(FORTRAN_PREPROCESS_FLAG "-x f95-cpp-input")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model source -real-size 64")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv")
  set(FORTRAN_PREPROCESS_FLAG "-fpp")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

# Add source files for NUOPC cap
list(APPEND _noahmp_cap_files noahmp/drivers/nuopc/lnd_comp_kind.F90
                              noahmp/drivers/nuopc/lnd_comp_types.F90
                              noahmp/drivers/nuopc/lnd_comp_shr.F90
                              noahmp/drivers/nuopc/lnd_comp_io.F90
                              noahmp/drivers/nuopc/lnd_comp_domain.F90
                              noahmp/drivers/nuopc/lnd_comp_import_export.F90
                              noahmp/drivers/nuopc/lnd_comp_nuopc.F90
                              noahmp/drivers/nuopc/lnd_comp_driver.F90)

# Add source files for CCPP interface 
list(APPEND _noahmp_ccpp_files noahmp/drivers/ccpp/sfc_noahmp_drv.F90
                               noahmp/drivers/ccpp/sfc_diff.f
                               noahmp/drivers/ccpp/machine.F
                               noahmp/drivers/ccpp/noahmp_tables.f90
                               noahmp/drivers/ccpp/namelist_soilveg.f
                               noahmp/drivers/ccpp/set_soilveg.f
                               noahmp/drivers/ccpp/funcphys.f90
                               noahmp/drivers/ccpp/physcons.F90)

# Add source files for NoahMP
list(APPEND _noahmp_files noahmp/src/module_sf_noahmplsm.f90
                          noahmp/src/module_sf_noahmp_glacier.f90)
add_definitions(-DCCPP)

# NOAHMP
add_library(noahmp STATIC ${_noahmp_cap_files} ${_noahmp_ccpp_files} ${_noahmp_files})
set_target_properties(noahmp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
set_target_properties(noahmp PROPERTIES COMPILE_FLAGS ${FORTRAN_PREPROCESS_FLAG})
target_include_directories(noahmp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(noahmp PUBLIC esmf fms)

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
