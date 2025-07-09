#!/bin/bash
set -eu

SCRIPT_REALPATH=$(realpath "${BASH_SOURCE[0]}")
UFS_MODEL_DIR=$(dirname "${SCRIPT_REALPATH}")
readonly UFS_MODEL_DIR
echo "UFS MODEL DIR: ${UFS_MODEL_DIR}"

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
OMP_NUM_THREADS=1 make -j "${BUILD_JOBS:-4}" "VERBOSE=${BUILD_VERBOSE:-}"
