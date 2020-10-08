#!/bin/bash
set -eu

if [[ $(uname -s) == Darwin ]]; then
  readonly UFS_MODEL_DIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly UFS_MODEL_DIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi

export CMAKE_C_COMPILER=${CMAKE_C_COMPILER:-${CC:-mpicc}}
export CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER:-${CXX:-mpicxx}}
export CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER:-${FC:-mpif90}}

export NETCDF=${NETCDF:-${NETCDF_ROOT:?"Please set NETCDF/NETCDF_ROOT environment variable"}}
export ESMFMKFILE=${ESMFMKFILE:?"Please set ESMFMKFILE environment variable"}

BUILD_DIR=${BUILD_DIR:-${UFS_MODEL_DIR}/build}
mkdir -p ${BUILD_DIR}

[[ -n "${CCPP_SUITES:-""}" ]] && CMAKE_FLAGS+=" -DCCPP_SUITES=${CCPP_SUITES}"
CMAKE_FLAGS+=" -DNETCDF_DIR=${NETCDF}"

cd ${BUILD_DIR}
cmake ${UFS_MODEL_DIR} ${CMAKE_FLAGS}
make -j ${BUILD_JOBS:-4} VERBOSE=${BUILD_VERBOSE:-0}
