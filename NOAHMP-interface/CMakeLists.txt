# Compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -ffree-form")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs -free")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

# Too many files to list, so include them via this file
include("noahmp_files.cmake")

# NOAHMP
add_library(noahmp STATIC ${_noahmp_src_files})
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
