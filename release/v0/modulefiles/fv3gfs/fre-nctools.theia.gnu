module load gcc/6.2.0

# mvapich2-2.2, compiled with gnu/6.2.0
prepend-path PATH            /scratch4/BMC/gmtb/mvapich2-2.2/gnu-6.2.0/bin
prepend-path LD_LIBRARY_PATH /scratch4/BMC/gmtb/mvapich2-2.2/gnu-6.2.0/lib
setenv MV2_ENABLE_AFFINITY 0

# netcdf-4.5.0, compiled with gnu/6.2.0 and mvapich2-2.2, and its dependencies
prepend-path PATH            /scratch4/BMC/gmtb/zlib-1.2.11/gnu-6.2.0/bin
prepend-path PATH            /scratch4/BMC/gmtb/szip-2.1.1/gnu-6.2.0/bin
prepend-path PATH            /scratch4/BMC/gmtb/hdf5-1.8.20/gnu-6.2.0/mvapich2-2.2/bin
prepend-path PATH            /scratch4/BMC/gmtb/parallel-netcdf-1.8.1/gnu-6.2.0/mvapich2-2.2/bin
prepend-path PATH            /scratch4/BMC/gmtb/netcdf-4.5.0/gnu-6.2.0/mvapich2-2.2/bin
prepend-path LD_LIBRARY_PATH /scratch4/BMC/gmtb/zlib-1.2.11/gnu-6.2.0/lib
prepend-path LD_LIBRARY_PATH /scratch4/BMC/gmtb/szip-2.1.1/gnu-6.2.0/lib
prepend-path LD_LIBRARY_PATH /scratch4/BMC/gmtb/hdf5-1.8.20/gnu-6.2.0/mvapich2-2.2/lib
prepend-path LD_LIBRARY_PATH /scratch4/BMC/gmtb/parallel-netcdf-1.8.1/gnu-6.2.0/mvapich2-2.2/lib
prepend-path LD_LIBRARY_PATH /scratch4/BMC/gmtb/netcdf-4.5.0/gnu-6.2.0/mvapich2-2.2/lib
setenv NETCDF /scratch4/BMC/gmtb/netcdf-4.5.0/gnu-6.2.0/mvapich2-2.2
