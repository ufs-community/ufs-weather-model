mpif90 -C -traceback -o compare_output -I${NETCDF}/include compare_output.F90 -L${NETCDF}/lib -lnetcdff -lnetcdf -L${HDF5_LIBRARIES} -lhdf5_hl -lhdf5 -L${ZLIB_LIBRARIES} -lz
