# CDEPS compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -ffree-line-length-none")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )
  set(CMAKE_Fortran_LINK_FLAGS "")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "-g -traceback")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -assume realloc_lhs")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -fp-model precise")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check -check noarg_temp_created -check nopointer -fpe0 -ftrapuv")
  set(CMAKE_Fortran_LINK_FLAGS "")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

# Too many files to list, so include them via this file
include("cdeps_files.cmake")

# CDEPS/share
add_library(share STATIC ${cdeps_share_files} ${ufs_cdeps_share_files})
add_library(cdeps::share ALIAS share)
set_target_properties(share PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(share PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                        $<INSTALL_INTERFACE:mod>)
target_link_libraries(share PUBLIC esmf)
list(APPEND LIB_TARGETS share)

# CDEPS/streams
add_library(streams STATIC ${cdeps_streams_files})
add_library(cdeps::streams ALIAS streams)
add_dependencies(streams cdeps::share)
set_target_properties(streams PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(streams PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                          $<INSTALL_INTERFACE:mod>)
target_compile_definitions(streams PUBLIC "DISABLE_FoX")
target_link_libraries(streams PUBLIC cdeps::share)
target_link_libraries(streams PUBLIC esmf PIO::PIO_Fortran)
list(APPEND LIB_TARGETS streams)

# CDEPS/dshr
add_library(dshr STATIC ${cdeps_dshr_files})
add_library(cdeps::dshr ALIAS dshr)
add_dependencies(dshr cdeps::share cdeps::streams)
set_target_properties(dshr PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dshr PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(dshr PUBLIC cdeps::share cdeps::streams)
target_link_libraries(dshr PUBLIC esmf PIO::PIO_Fortran)
list(APPEND LIB_TARGETS dshr)

# CDEPS/datm
add_library(datm STATIC ${cdeps_datm_files})
add_library(cdeps::datm ALIAS datm)
add_dependencies(datm cdeps::share cdeps::streams cdeps::dshr)
set_target_properties(datm PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(datm PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(datm PUBLIC cdeps::dshr cdeps::share cdeps::streams)
target_link_libraries(datm PUBLIC esmf PIO::PIO_Fortran)
list(APPEND LIB_TARGETS datm)

# CDEPS/dice
add_library(dice STATIC ${cdeps_dice_files})
add_library(cdeps::dice ALIAS dice)
add_dependencies(dice cdeps::share cdeps::streams cdeps::dshr)
set_target_properties(dice PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dice PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(dice PUBLIC cdeps::dshr cdeps::share cdeps::streams)
target_link_libraries(dice PUBLIC esmf)
list(APPEND LIB_TARGETS dice)

# CDEPS/dlnd
add_library(dlnd STATIC ${cdeps_dlnd_files})
add_library(cdeps::dlnd ALIAS dlnd)
add_dependencies(dlnd cdeps::share cdeps::streams cdeps::dshr)
set_target_properties(dlnd PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dlnd PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(dlnd PUBLIC cdeps::dshr cdeps::share cdeps::streams)
target_link_libraries(dlnd PUBLIC esmf)
list(APPEND LIB_TARGETS dlnd)

# CDEPS/docn
add_library(docn STATIC ${cdeps_docn_files})
add_library(cdeps::docn ALIAS docn)
add_dependencies(docn cdeps::share cdeps::streams cdeps::dshr)
set_target_properties(docn PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(docn PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(docn PUBLIC cdeps::dshr cdeps::share cdeps::streams)
target_link_libraries(docn PUBLIC esmf)
list(APPEND LIB_TARGETS docn)

# CDEPS/drof
add_library(drof STATIC ${cdeps_drof_files})
add_library(cdeps::drof ALIAS drof)
add_dependencies(drof cdeps::share cdeps::streams cdeps::dshr)
set_target_properties(drof PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(drof PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(drof PUBLIC cdeps::dshr cdeps::share cdeps::streams)
target_link_libraries(drof PUBLIC esmf)
list(APPEND LIB_TARGETS drof)

# CDEPS/dwav
add_library(dwav STATIC ${cdeps_dwav_files})
add_library(cdeps::dwav ALIAS dwav)
add_dependencies(dwav cdeps::share cdeps::streams cdeps::dshr)
set_target_properties(dwav PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(dwav PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(dwav PUBLIC cdeps::dshr cdeps::share cdeps::streams)
target_link_libraries(dwav PUBLIC esmf)
list(APPEND LIB_TARGETS dwav)

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
