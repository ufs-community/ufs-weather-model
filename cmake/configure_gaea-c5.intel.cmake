set(PARALLEL_NETCDF ON  CACHE BOOL "Enable parallel NetCDF" FORCE)
set(MOM6_Extra_FORTRAN_FLAGS "-xsse2")
set(HYCOM_Extra_FORTRAN_FLAGS "-xSSE4.2")
set(HYCOM_Extra_C_FLAGS "-xSSE4.2")