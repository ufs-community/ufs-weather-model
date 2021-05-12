#!/bin/bash

#%Module######################################################################
##
##    UFS prerequisites: MACOSX with clang/gfortran compilers

echo "Setting environment variables for NEMSfv3gfs on MACOSX with gcc/gfortran or clang/gfortran"

##
## load programming environment: compiler, flags, paths
##
export CC=${MPICC:-mpicc}
export CXX=${MPICXX:-mpicxx}
export F77=${MPIF77:-mpif77}
export F90=${MPIF90:-mpif90}
export FC=${MPIFORT:-mpifort}
export CPP=${CPP:-"${F90} -E -x f95-cpp-input"}
export MPICC=${MPICC:-mpicc}
export MPIF90=${MPIF90:-mpif90}

##
## load cmake
##
export CMAKE_Platform=macosx.gnu

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
## Intel MKL library
##
#export MKL_DIR=${MKL_DIR:-/opt/intel/compilers_and_libraries_2019.4.233/mac/mkl}
#export MKL_INC="-m64 -I${MKL_DIR}/include"
#export MKL_LIB="-L${MKL_DIR}/lib -Wl,-rpath,${MKL_DIR}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
