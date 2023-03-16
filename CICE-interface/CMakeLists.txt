### CICE Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fconvert=big-endian -ffree-line-length-none ")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -Wall -Wextra -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -FR -convert big_endian -assume byterecl -ftz -align array64byte")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created -init=snan,arrays")
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
  set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
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

# Configuration Options
set(CICE_IO "NetCDF" CACHE STRING "CICE OPTIONS: Choose IO options.")
set_property(CACHE CICE_IO PROPERTY STRINGS "NetCDF" "PIO" "Binary")

# Too many files to list, so include them via this file
include("cice_files.cmake")

# Collect source files for library
list(APPEND lib_src_files
  ${cice_shared_files}
  ${cice_shared_files_c}
  ${icepack_files}
  ${cice_mpi_comm_files}
  ${cice_nuopc_cmeps_driver_files})

list(APPEND _cice_defs FORTRANUNDERSCORE
                       coupled)

# Select IO source files based on CICE_IO
if(CICE_IO MATCHES "NetCDF")
  list(APPEND lib_src_files ${cice_netcdf_io_files})
  list(APPEND _cice_defs USE_NETCDF)
elseif(CICE_IO MATCHES "PIO")
  list(APPEND lib_src_files ${cice_pio2_io_files})
  list(APPEND _cice_defs USE_NETCDF)
elseif(CICE_IO MATCHES "Binary")
  list(APPEND lib_src_files ${cice_binary_io_files})
endif()

### Create target library and set PUBLIC interfaces on the library
add_library(cice STATIC ${lib_src_files})
set_target_properties(cice PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_compile_definitions(cice PUBLIC "${_cice_defs}")
target_include_directories(cice PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(cice PUBLIC esmf)
if(CICE_IO MATCHES "^(NetCDF|PIO)$")
  target_link_libraries(cice PUBLIC NetCDF::NetCDF_Fortran)
  if(CICE_IO MATCHES "PIO")
    target_link_libraries(cice PUBLIC PIO::PIO_Fortran)
  endif()
endif()
if(OpenMP_Fortran_FOUND)
  target_link_libraries(cice PRIVATE OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### Install
###############################################################################

install(
  TARGETS cice
  EXPORT  cice-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/CICE)

install(EXPORT      cice-config
        DESTINATION lib/cmake)
