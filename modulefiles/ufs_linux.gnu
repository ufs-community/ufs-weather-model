#!/bin/bash

#%Module######################################################################
##
##    UFS prerequisites: Linux (tested: Ubuntu) with gcc/gfortran compilers

echo "Setting environment variables for UFS Model on Linux with gcc/gfortran"

##
## load programming environment: compiler, flags, paths
##
export CC=${CC:-mpicc}
export CXX=${CXX:-mpicxx}
export F77=${F77:-mpif77}
export F90=${F90:-mpif90}
export FC=${FC:-mpif90}

##
## set up variables for ../cmake/configure_linux.gnu.cmake
##
export CMAKE_Platform=linux.gnu

##
## use own NetCDF library
##
export NETCDF=${NETCDF:-/home/builder/opt}

##
## use pre-compiled EMSF library for above compiler / MPI combination
##
export ESMFMKFILE=${ESMFMKFILE:-/home/builder/opt/lib/esmf.mk}

##
## NCEP libraries (need to download and build manually, see doc/README_{UBUNTU,CENTOS,...}.txt and https://github.com/NCAR/NCEPlibs)
##
export NCEPLIBS_DIR=${NCEPLIBS_DIR:-/home/builder/opt}
export CMAKE_PREFIX_PATH=${NCEPLIBS_DIR}
