# CDEPS compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace -ffree-line-length-none")
  #Set CPP defintions for the ufs/cdeps_share target library
  list(APPEND CDEPS_SHARE_DEFS "CPRGNU")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv -init=snan,arrays")
  #Set CPP defintions for the ufs/cdeps_share target library
  list(APPEND CDEPS_SHARE_DEFS "CPRINTEL")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

# Too many files to list, so include them via this file
include("cdeps_files.cmake")

# CDEPS/share
add_library(share OBJECT ${cdeps_share_files} ${ufs_cdeps_share_files})
set_target_properties(share PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(share PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ufs/cdeps_share)
target_include_directories(share PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_compile_definitions(share PRIVATE ${CDEPS_SHARE_DEFS})
target_link_libraries(share PRIVATE esmf)

# CDEPS/streams
add_library(streams OBJECT ${cdeps_streams_files})
add_dependencies(streams share)
set_target_properties(streams PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(streams PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_compile_definitions(streams PUBLIC "DISABLE_FoX")
target_link_libraries(streams PRIVATE share)
target_link_libraries(streams PRIVATE esmf PIO::PIO_Fortran)

# CDEPS/dshr
add_library(dshr OBJECT ${cdeps_dshr_files})
add_dependencies(dshr share streams)
set_target_properties(dshr PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dshr PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(dshr PRIVATE share streams)
target_link_libraries(dshr PRIVATE esmf PIO::PIO_Fortran)

# CDEPS/datm
add_library(datm OBJECT ${cdeps_datm_files})
add_dependencies(datm share streams dshr)
set_target_properties(datm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(datm PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(datm PRIVATE dshr share streams)
target_link_libraries(datm PRIVATE esmf PIO::PIO_Fortran)

# CDEPS/dice
add_library(dice OBJECT ${cdeps_dice_files})
add_dependencies(dice share streams dshr)
set_target_properties(dice PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dice PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(dice PRIVATE dshr share streams)
target_link_libraries(dice PRIVATE esmf)

# CDEPS/dlnd
add_library(dlnd OBJECT ${cdeps_dlnd_files})
add_dependencies(dlnd share streams dshr)
set_target_properties(dlnd PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dlnd PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(dlnd PRIVATE dshr share streams)
target_link_libraries(dlnd PRIVATE esmf)

# CDEPS/docn
add_library(docn OBJECT ${cdeps_docn_files})
add_dependencies(docn share streams dshr)
set_target_properties(docn PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(docn PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(docn PRIVATE dshr share streams)
target_link_libraries(docn PRIVATE esmf PIO::PIO_Fortran)

# CDEPS/drof
add_library(drof OBJECT ${cdeps_drof_files})
add_dependencies(drof share streams dshr)
set_target_properties(drof PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(drof PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(drof PRIVATE dshr share streams)
target_link_libraries(drof PRIVATE esmf)

# CDEPS/dwav
add_library(dwav OBJECT ${cdeps_dwav_files})
add_dependencies(dwav share streams dshr)
set_target_properties(dwav PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dwav PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)
target_link_libraries(dwav PRIVATE dshr share streams)
target_link_libraries(dwav PRIVATE esmf)

# CDEPS
add_library(cdeps STATIC $<TARGET_OBJECTS:share>
                         $<TARGET_OBJECTS:streams>
                         $<TARGET_OBJECTS:dshr>
                         $<TARGET_OBJECTS:datm>
                         $<TARGET_OBJECTS:dice>
                         $<TARGET_OBJECTS:dlnd>
                         $<TARGET_OBJECTS:docn>
                         $<TARGET_OBJECTS:drof>
                         $<TARGET_OBJECTS:dwav>)
add_library(cdeps::cdeps ALIAS cdeps)
target_include_directories(cdeps PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                        $<INSTALL_INTERFACE:mod>)
target_link_libraries(cdeps PUBLIC esmf PIO::PIO_Fortran)
list(APPEND LIB_TARGETS cdeps)

###############################################################################
### Install
###############################################################################

install(
  TARGETS   ${LIB_TARGETS}
  EXPORT    CDEPSExports
  LIBRARY   DESTINATION lib
  ARCHIVE   DESTINATION lib
  COMPONENT Library)

export(EXPORT    CDEPSExports
       NAMESPACE cdeps::
       FILE      cdeps-targets.cmake)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/CDEPS)

install(EXPORT      CDEPSExports
        NAMESPACE   cdeps::
        FILE        cdeps-targets.cmake
        DESTINATION lib/cmake)
