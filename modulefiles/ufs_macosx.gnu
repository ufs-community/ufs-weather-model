#!/bin/bash

#%Module######################################################################
##
##    UFS prerequisites: macOS with gcc/gfortran or clang/gfortran compilers

echo "Setting environment variables for UFS Model on macOS with gcc/gfortran or clang/gfortran"

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
