set(PARALLEL_NETCDF ON  CACHE BOOL "Enable parallel NetCDF" FORCE)
set(DEBUG_LINKMPI   OFF CACHE BOOL "Enable linkmpi option when DEBUG mode is on" FORCE)
set(AVX2            OFF CACHE BOOL "Enable AVX2 instruction set" FORCE)

if(FASTER)
  set(DISABLE_FMA   ON  CACHE BOOL "Disable Fused Multiply-Add instructions (workaround needed for AMD EPYC)" FORCE)
endif()
