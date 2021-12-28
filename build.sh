#!/bin/bash
set -eu

if [[ $(uname -s) == Darwin ]]; then
  readonly UFS_MODEL_DIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly UFS_MODEL_DIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi

export CC=${CC:-mpicc}
export CXX=${CXX:-mpicxx}
export FC=${FC:-mpif90}

BUILD_DIR=${BUILD_DIR:-${UFS_MODEL_DIR}/build}
mkdir -p ${BUILD_DIR}

cd ${BUILD_DIR}
cmake ${UFS_MODEL_DIR} ${CMAKE_FLAGS}
# Turn off OpenMP threading for parallel builds
# to avoid exhausting the number of user processes
OMP_NUM_THREADS=1 make -j ${BUILD_JOBS:-4} VERBOSE=${BUILD_VERBOSE:-}
