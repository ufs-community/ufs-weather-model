get_filename_component(hycom_dir "${CMAKE_CURRENT_SOURCE_DIR}/HYCOM" ABSOLUTE)

### HYCOM Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Waliasing -fcray-pointer -fconvert=big-endian -ffree-line-length-none -fno-range-check -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback -mcmodel=small -r8 ${HYCOM_Extra_FORTRAN_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model source -warn nogeneral")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all -fpe0 -ftrapuv -link_mpi=dbg -init=snan,arrays")
  set(CMAKE_Fortran_LINK_FLAGS "-V ${CMAKE_Fortran_FLAGS} -static-intel")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

### HYCOM C compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS "-g -fbacktrace")
  set(CMAKE_C_FLAGS_RELEASE "-O3")
  set(CMAKE_C_FLAGS_DEBUG "-O0")
  set(CMAKE_C_LINK_FLAGS "")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS "-g -traceback -mcmodel=small ${HYCOM_Extra_C_FLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "-O")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -ftrapuv")
  set(CMAKE_C_LINK_FLAGS "-V ${CMAKE_C_FLAGS} -static-intel")
else()
  message(WARNING "C compiler with ID ${CMAKE_C_COMPILER_ID} will be used with CMake default options")
endif()

list(APPEND _hycom_defs IA32
                        REAL8
                        MPI
                        SERIAL_IO
                        ENDIAN_IO
                        NAN2003
                        TIME
                        RELO
                        EOS_SIG2
                        EOS_17T
                        ESPC_COUPLE
)

# User option to build HYCOM offline executable
set(HYCOMOFFLINE OFF CACHE BOOL "Build HYCOM offline")

message("Build HYCOM:")
message("   in:  ${hycom_dir}")
message("   HYCOMOFFLINE: ${HYCOMOFFLINE}")
message("")

# Too many files to list, so include them via this file
include("hycom_files.cmake")

### Use common object library for building target library
add_library(hycom_obj OBJECT ${_hycom_src_files})
set_target_properties(hycom_obj PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/HYCOM)
set_target_properties(hycom_obj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(hycom_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_include_directories(hycom_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/HYCOM>)
target_compile_definitions(hycom_obj PRIVATE "${_hycom_defs}")


### Use NUOPC object library for building target library
add_library(hycom_nuopc_obj OBJECT ${_hycom_nuopc_src_files})
add_dependencies(hycom_nuopc_obj hycom_obj)
set_target_properties(hycom_nuopc_obj PROPERTIES ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/HYCOM)
set_target_properties(hycom_nuopc_obj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(hycom_nuopc_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_include_directories(hycom_nuopc_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/HYCOM>
                                                   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/HYCOM/NUOPC>)
target_compile_definitions(hycom_nuopc_obj PRIVATE "${_hycom_defs}")
target_link_libraries(hycom_nuopc_obj PRIVATE hycom_obj
                                              esmf
                                              NetCDF::NetCDF_Fortran)

### Create target library and set PUBLIC interfaces on the library
add_library(hycom STATIC $<TARGET_OBJECTS:hycom_obj>
                         $<TARGET_OBJECTS:hycom_nuopc_obj>)
target_include_directories(hycom PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                        $<INSTALL_INTERFACE:mod>)
target_include_directories(hycom PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/HYCOM>
                                        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/HYCOM/NUOPC>)
target_compile_definitions(hycom PRIVATE "${_hycom_defs}")
target_link_libraries(hycom PUBLIC esmf
                                   NetCDF::NetCDF_Fortran)

### Create HYCOM offline executable
if(HYCOMOFFLINE)
  message("Building HYCOM offline executable")
  add_executable(hycomoffline ${_hycom_offline_src_files})
  add_dependencies(hycomoffline hycom_obj)
  set_target_properties(hycomoffline PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_offline)
  target_include_directories(hycomoffline PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
  target_include_directories(hycomoffline PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_offline>)
  target_include_directories(hycomoffline PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/HYCOM>)

  target_link_libraries(hycomoffline PRIVATE hycom_obj)
endif()

###############################################################################
### Install
###############################################################################

install(
  TARGETS hycom
  EXPORT  hycom-config
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/HYCOM)

install(EXPORT      hycom-config
        DESTINATION lib/cmake)

if(HYCOMOFFLINE)
  install(TARGETS hycomoffline RUNTIME DESTINATION bin)
endif()
