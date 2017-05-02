module load fre/bronx-10
module load intel_compilers/11.1.073 mpich2/1.2.1p1
setenv SITE pan

alias make make HDF5_HOME=/usr/local/hdf5-1.8.8_optimized NETCDF_HOME=/usr/local/netcdf-4.2_optimized -f fre-nctools.mk


