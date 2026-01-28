### MOM6 Fortran compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace")
  list(APPEND CDEPS_SHARE_DEFS "CPRGNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Waliasing -fcray-pointer -fconvert=big-endian -ffree-line-length-none -fno-range-check -fbacktrace")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl -sox")
  if(FASTER)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise -assume buffered_stdout -fno-alias -align all -debug minimal -qoverride-limits -ftz -no-ip")
  else()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model source")
  endif()
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  list(APPEND CDEPS_SHARE_DEFS "CPRINTEL")
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r8")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-alias -auto -safe-cray-ptr -ftz -assume byterecl")
  if(FASTER)
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fp-model precise -assume buffered_stdout -fno-alias -align all -debug minimal -qoverride-limits -ftz")
  else()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal -fp-model precise")
  endif()
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  set(CMAKE_Fortran_LINK_FLAGS "")
  list(APPEND CDEPS_SHARE_DEFS "CPRINTEL")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MOM6_Extra_FORTRAN_FLAGS}")

# User option to build MOM6 solo executable
set(MOM6SOLO OFF CACHE BOOL "Build MOM6 solo")

# Too many files to list, so include them via this file
include("mom6_files.cmake")

### Use common object library for building target library
add_library(mom6_obj OBJECT ${mom6_src_files})
set_target_properties(mom6_obj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(mom6_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_compile_definitions(mom6_obj PRIVATE ${CDEPS_SHARE_DEFS} "DISABLE_FoX")
target_include_directories(mom6_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/CDEPS-interface/mod>)
add_dependencies(mom6_obj cdeps)

if(REGIONAL_MOM6)
  target_include_directories(mom6_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_symmetric>
                                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
else()
  target_include_directories(mom6_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_nonsymmetric>
                                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
endif()
target_link_libraries(mom6_obj PRIVATE fms
                                       ESMF::ESMF
                                       stochastic_physics
                                       NetCDF::NetCDF_Fortran
                                       cdeps::cdeps)
# OpenMP is disabled in MOM6
#if(OpenMP_Fortran_FOUND)
# target_link_libraries(mom6_obj PRIVATE OpenMP::OpenMP_Fortran)
#endif()

add_library(mom6_nuopc_obj OBJECT ${mom6_nuopc_src_files})
set_target_properties(mom6_nuopc_obj PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(mom6_nuopc_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_include_directories(mom6_nuopc_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/CDEPS-interface/mod>)
add_dependencies(mom6_nuopc_obj cdeps)
if(REGIONAL_MOM6)
  target_include_directories(mom6_nuopc_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_symmetric>
                                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
else()
  target_include_directories(mom6_nuopc_obj PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_nonsymmetric>
                                                    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
endif()
target_link_libraries(mom6_nuopc_obj PRIVATE mom6_obj
                                             fms
                                             stochastic_physics
                                             ESMF::ESMF)
if(UFS_TRACING)
  target_link_libraries(mom6_nuopc_obj PRIVATE ufs_tracing)
endif()
# OpenMP is disabled in MOM6
#if(OpenMP_Fortran_FOUND)
# target_link_libraries(mom6_nuopc_obj PRIVATE OpenMP::OpenMP_Fortran)
#endif()

### Create target library and set PUBLIC interfaces on the library
add_library(mom6 STATIC $<TARGET_OBJECTS:mom6_obj>
                        $<TARGET_OBJECTS:mom6_nuopc_obj>)
target_include_directories(mom6 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_include_directories(mom6 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_nonsymmetric>
                                       $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
target_link_libraries(mom6 PUBLIC fms
                                  stochastic_physics
                                  ESMF::ESMF
                                  NetCDF::NetCDF_Fortran
                                  cdeps::cdeps)
if(UFS_TRACING)
  target_link_libraries(mom6 PRIVATE ufs_tracing)
endif()
# OpenMP is disabled in MOM6
#if(OpenMP_Fortran_FOUND)
#  target_link_libraries(mom6 PRIVATE OpenMP::OpenMP_Fortran)
#endif()

### Create standalone MOM6 executable
if(MOM6SOLO)
  message("Building MOM6 standalone executable")
  add_executable(mom6solo ${mom6_solo_src_files})
  add_dependencies(mom6solo mom6_obj)
  set_target_properties(mom6solo PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_solo)
  target_include_directories(mom6solo PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
  target_include_directories(mom6solo PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_solo>)
  target_include_directories(mom6solo PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/config_src/memory/dynamic_nonsymmetric>
                                              $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/MOM6/src/framework>)
  target_link_libraries(mom6solo PRIVATE mom6_obj
                                         fms
                                         stochastic_physics
                                         ESMF::ESMF
                                         NetCDF::NetCDF_Fortran)
  # OpenMP is disabled in MOM6
  #if(OpenMP_Fortran_FOUND)
  #  target_link_libraries(mom6solo PRIVATE OpenMP::OpenMP_Fortran)
  #endif()
endif()

###############################################################################
### Install
###############################################################################

install(
  TARGETS mom6
  EXPORT  mom6-config
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/MOM6)

install(EXPORT      mom6-config
        DESTINATION lib/cmake)

if(MOM6SOLO)
  install(TARGETS mom6solo RUNTIME DESTINATION bin)
endif()
