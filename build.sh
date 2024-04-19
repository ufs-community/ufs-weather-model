#!/bin/bash
set -eu
uname_s=$(uname -s)
if [[ ${uname_s} == Darwin ]]; then
  UFS_MODEL_DIR=$(greadlink -f -n "${BASH_SOURCE[0]}")
  UFS_MODEL_DIR=$(dirname "${UFS_MODEL_DIR}")
  UFS_MODEL_DIR=$(cd "${UFS_MODEL_DIR}" && pwd -P)
else
  UFS_MODEL_DIR=$(readlink -f -n "${BASH_SOURCE[0]}")
  UFS_MODEL_DIR=$(dirname "${UFS_MODEL_DIR}")
  UFS_MODEL_DIR=$(cd "${UFS_MODEL_DIR}" && pwd -P)
fi
echo "UFS MODEL DIR: ${UFS_MODEL_DIR}"
readonly UFS_MODEL_DIR

export CC=${CC:-mpicc}
export CXX=${CXX:-mpicxx}
export FC=${FC:-mpif90}

BUILD_DIR=${BUILD_DIR:-${UFS_MODEL_DIR}/build}
mkdir -p "${BUILD_DIR}"

cd "${BUILD_DIR}"
ARR_CMAKE_FLAGS=()
for i in ${CMAKE_FLAGS}; do ARR_CMAKE_FLAGS+=("${i}") ; done
cmake "${UFS_MODEL_DIR}" "${ARR_CMAKE_FLAGS[@]}"
# Turn off OpenMP threading for parallel builds
# to avoid exhausting the number of user processes
OMP_NUM_THREADS=1 make -j "${BUILD_JOBS:-4}" "VERBOSE=${BUILD_VERBOSE:-1}"