### CICE Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fconvert=big-endian -ffree-line-length-none ")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -FR -convert big_endian -assume byterecl -ftz -align array64byte")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -warn -warn noerrors -fpe0 -ftrapuv")
  set(CMAKE_Fortran_LINK_FLAGS "")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

### CICE C compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS "-g")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
  set(CMAKE_C_LINK_FLAGS "")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set( CMAKE_C_FLAGS "-g -traceback")
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHOST")
  set( CMAKE_C_FLAGS_RELEASE "-O2 -fp-model precise")
  set( CMAKE_C_FLAGS_DEBUG "-O0 -ftrapuv")
  set( CMAKE_C_LINK_FLAGS "")
elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(CMAKE_C_FLAGS "-g")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  set( CMAKE_C_FLAGS_RELEASE "-O3" )
  set( CMAKE_C_FLAGS_DEBUG "-O0" )
  set( CMAKE_C_LINK_FLAGS "" )
else()
  message(WARNING "C compiler with ID ${CMAKE_C_COMPILER_ID} will be used with CMake default options")
endif()

# Too many files to list, so include them via this file
include("cice_files.cmake")

# Collect source files for library
list(APPEND lib_src_files
  ${cice_shared_files}
  ${cice_shared_files_c}
  ${icepack_files}
  ${cice_mpi_comm_files}
  ${cice_nuopc_cmeps_driver_files})

# NetCDF and PIO are both required in CICE
if(NOT NetCDF_Fortran_FOUND OR NOT PIO_Fortran_FOUND)
  message(FATAL_ERROR "NetCDF or PIO required")
else()
  list(APPEND lib_src_files ${cice_netcdf_io_files})
  list(APPEND lib_src_files ${cice_pio2_io_files})
endif()

add_library(cice ${lib_src_files})
set_target_properties(cice PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(cice INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                          $<INSTALL_INTERFACE:mod>)
target_compile_definitions(cice PRIVATE FORTRANUNDERSCORE
                                        coupled
                                        USE_NETCDF)
target_link_libraries(cice PUBLIC esmf
                                  NetCDF::NetCDF_Fortran
                                  PIO::PIO_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(cice PRIVATE OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### Install
###############################################################################

install(
  TARGETS cice
  EXPORT cice-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/CICE)

install(EXPORT      cice-config
        DESTINATION lib/cmake)
