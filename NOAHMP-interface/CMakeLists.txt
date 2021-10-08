# Compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(FORTRAN_PREPROCESS_FLAG "-x f95-cpp-input")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv")
  set(FORTRAN_PREPROCESS_FLAG "-fpp")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

list(APPEND _noahmp_cap_files NOAHMP/drivers/nuopc/lnd_comp_kind.F90
                              NOAHMP/drivers/nuopc/lnd_comp_types.F90
                              NOAHMP/drivers/nuopc/lnd_comp_shr.F90
                              NOAHMP/drivers/nuopc/lnd_import_export.F90
                              NOAHMP/drivers/nuopc/lnd_comp_nuopc.F90
                              NOAHMP/drivers/nuopc/lnd_comp_driver.F90)

get_filename_component(FV3_PATH "../FV3" REALPATH)
if(EXISTS ${FV3_PATH})
  message("FV3_PATH = ${FV3_PATH}")
  message("Compiling NOAHMP using CCPP source code")
  list(APPEND _ccpp_physics
    ${FV3_PATH}/ccpp/physics/physics/machine.F
    ${FV3_PATH}/ccpp/physics/physics/funcphys.f90
    ${FV3_PATH}/ccpp/physics/physics/physcons.F90
    ${FV3_PATH}/ccpp/physics/physics/set_soilveg.f
    ${FV3_PATH}/ccpp/physics/physics/namelist_soilveg.f
    ${FV3_PATH}/ccpp/physics/physics/module_sf_noahmplsm.f90
    ${FV3_PATH}/ccpp/physics/physics/module_sf_noahmp_glacier.f90
    ${FV3_PATH}/ccpp/physics/physics/sfc_noahmp_drv.F90)
  add_definitions(-DCCPP)
endif()

# NOAHMP
add_library(noahmp STATIC ${_noahmp_cap_files} ${_ccpp_physics})
set_target_properties(noahmp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
set_target_properties(noahmp PROPERTIES COMPILE_FLAGS ${FORTRAN_PREPROCESS_FLAG})
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
