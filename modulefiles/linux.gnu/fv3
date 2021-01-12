#!/bin/bash

#%Module######################################################################
##
##    UFS prerequisites: Linux (tested: Ubuntu) with gcc/gfortran compilers

echo "Setting environment variables for NEMSfv3gfs on Linux with gcc/gfortran"

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
export CMAKE_C_COMPILER=${CC}
export CMAKE_CXX_COMPILER=${CXX}
export CMAKE_Fortran_COMPILER=${FC}

##
## use own NetCDF library
##
export NETCDF=${NETCDF:-/home/builder/opt}

##
## use SIONlib library if installed and environment variable is set
##
SIONLIB=${SIONLIB:-}
if [ ! "x$SIONLIB" == "x" ]; then
  echo "Use SIONlib installation in ${SIONLIB}"
  export SIONLIB_INC="-I${SIONLIB}/include -I${SIONLIB}/include/mod_64"
  export SIONLIB_LIB="-L${SIONLIB}/lib -lsionmpi_f90_64 -lsionser_f90_64 -lsionmpi_64 -lsiongen_64 -lsionser_64 -lsioncom_64 -lsioncom_64_lock_none"
fi

##
## use pre-compiled EMSF library for above compiler / MPI combination
##
export ESMFMKFILE=${ESMFMKFILE:-/home/builder/opt/lib/esmf.mk}

##
## NCEP libraries (need to download and build manually, see doc/README_{UBUNTU,CENTOS,...}.txt and https://github.com/NCAR/NCEPlibs)
##
export NCEPLIBS_DIR=${NCEPLIBS_DIR:-/home/builder/opt}
export CMAKE_PREFIX_PATH=${NCEPLIBS_DIR}
