if(32BIT)
  message ("Force 64 bits in stochastic_physics")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -real-size 64")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif()
endif()

list(APPEND _stoch_phys_srcs
  kinddef.F90
  mpi_wrapper.F90
  halo_exchange.fv3.F90
  plumes.F90
  mersenne_twister.F90
  random_numbers.F90
  stochy_namelist_def.F90
  compns_stochy.F90
  stochy_internal_state_mod.F90
  spectral_transforms.F90
  stochastic_physics.F90
  stochy_patterngenerator.F90
  stochy_data_mod.F90
  get_stochy_pattern.F90
  cellular_automata_global.F90
  cellular_automata_sgs.F90
  update_ca.F90
  lndp_apply_perts.F90
)

add_library(stochastic_physics ${_stoch_phys_srcs})
add_dependencies(stochastic_physics fms)

if(32BIT)
  target_compile_definitions(stochastic_physics PRIVATE INTERNAL_FILE_NML OVERLOAD_R4)
else()
  target_compile_definitions(stochastic_physics PRIVATE INTERNAL_FILE_NML)  
endif()

set_target_properties(stochastic_physics PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(stochastic_physics INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
                                                        $<INSTALL_INTERFACE:mod>)

target_link_libraries(stochastic_physics PUBLIC fms)
if(OPENMP)
  target_link_libraries(stochastic_physics PRIVATE OpenMP::OpenMP_Fortran)
endif()

###############################################################################
### Install
###############################################################################
install(
  TARGETS stochastic_physics
  EXPORT  stochastic_physics-config
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_PREFIX})

install(EXPORT stochastic_physics-config
    DESTINATION lib/cmake)
