#!/bin/bash
set -eu

MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)

export CMAKE_Platform=${CMAKE_Platform:?"Please set the CMAKE_Platform environment variable, e.g. [macosx.gnu|linux.gnu|linux.intel|hera.intel|...]"}
export CMAKE_C_COMPILER=${CMAKE_C_COMPILER:-mpicc}
export CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER:-mpicxx}
export CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER:-mpif90}

export NETCDF=${NETCDF:?"Please set NETCDF environment variable"}
export ESMFMKFILE=${ESMFMKFILE:?"Please set ESMFMKFILE environment variable"}

BUILD_DIR=${MYDIR}/build
rm -rf ${BUILD_DIR}
mkdir ${BUILD_DIR}

CCPP_SUITES="${CCPP_SUITES:-FV3_GFS_2017_gfdlmp}"
CMAKE_FLAGS+=" -DCCPP=ON -DCCPP_SUITES=${CCPP_SUITES} -DNETCDF_DIR=${NETCDF}"
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Debug -DDEBUG=YES -DWW3=ON"

cd ${BUILD_DIR}
cmake .. ${CMAKE_FLAGS}
make -j ${BUILD_JOBS:-4}
cp NEMS.exe ${MYDIR}/ufs_weather_model
