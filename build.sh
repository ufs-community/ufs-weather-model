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

export ESMFMKFILE=${ESMFMKFILE:?"Please set ESMFMKFILE environment variable"}

BUILD_DIR=${BUILD_DIR:-${UFS_MODEL_DIR}/build}
mkdir -p ${BUILD_DIR}

[[ -n "${CCPP_SUITES:-""}" ]] && CMAKE_FLAGS+=" -DCCPP_SUITES=${CCPP_SUITES}"

cd ${BUILD_DIR}
cmake ${UFS_MODEL_DIR} ${CMAKE_FLAGS}
make -j ${BUILD_JOBS:-4} VERBOSE=${BUILD_VERBOSE:-}
