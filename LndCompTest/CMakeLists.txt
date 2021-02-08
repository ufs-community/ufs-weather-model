#set(SRCS lnd_comp_nuopc.F90)

### Create target library and set PUBLIC interfaces on the library
add_library(lnd_comp STATIC
	    clm_varctl.F90
	    # decompInitMod.F90
	    # domainMod.F90
	    lnd_comp_nuopc.F90
	    lnd_import_export.F90
	    lnd_set_decomp_and_domain.F90
	    nuopc_shr_methods.F90
	    shr_file_mod.F90
	    shr_kind_mod.F90
	    shr_log_mod.F90
	    shr_sys_mod.F90
	    spmdMod.F90
	    )
set_target_properties(lnd_comp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_compile_definitions(lnd_comp PUBLIC "${_lnd_comp_defs}")
target_include_directories(lnd_comp PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                       $<INSTALL_INTERFACE:mod>)
target_link_libraries(lnd_comp PUBLIC esmf)

#target_link_libraries(lnd_comp PUBLIC esmf)

# test this:
target_include_directories(lnd_comp PRIVATE ${CMAKE_BINARY_DIR}/FV3/ccpp/physics)

###############################################################################
### Install
###############################################################################

# install(TARGETS lnd_comp)

install(
  TARGETS lnd_comp
  EXPORT  lnd_comp-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  COMPONENT Library)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX}/LND_COMP)

install(EXPORT      lnd_comp-config
        DESTINATION lib/cmake)

