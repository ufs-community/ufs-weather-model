#!/bin/bash
set -eu

MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)

export CMAKE_Platform=${CMAKE_Platform:?"Please set the CMAKE_Platform environment variable, e.g. [macosx.gnu|linux.gnu|linux.intel|hera.intel|...]"}
export CMAKE_C_COMPILER=${CMAKE_C_COMPILER:-mpicc}
export CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER:-mpicxx}
export CMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER:-mpif90}

export BACIO_LIB4=${BACIO_LIB4:?"Please set BACIO_LIB4 environment variable"}
export NEMSIO_INC=${NEMSIO_INC:?"Please set NEMSIO_INC environment variable"}
export NEMSIO_LIB=${NEMSIO_LIB:?"Please set NEMSIO_LIB environment variable"}
export SP_LIBd=${SP_LIBd:?"Please set SP_LIBd environment variable"}
export W3EMC_LIBd=${W3EMC_LIBd:?"Please set W3EMC_LIBd environment variable"}
export W3NCO_LIBd=${W3NCO_LIBd:?"Please set W3NCO_LIBd environment variable"}
export NETCDF=${NETCDF:?"Please set NETCDF environment variable"}
export ESMFMKFILE=${ESMFMKFILE:?"Please set ESMFMKFILE environment variable"}

BUILD_DIR=${MYDIR}/build
rm -rf ${BUILD_DIR}
mkdir ${BUILD_DIR}

CCPP_SUITES="${CCPP_SUITES:-FV3_GFS_v15p2}"

set +e
echo "Calling ccpp_prebuild.py ..."
./FV3/ccpp/framework/scripts/ccpp_prebuild.py \
    --config=FV3/ccpp/config/ccpp_prebuild_config.py \
    --static \
    --suites=${CCPP_SUITES} \
    --builddir=${BUILD_DIR}/FV3 > ${BUILD_DIR}/ccpp_prebuild.log 2>&1
if [ $? -ne 0 ]; then
  echo "Call to ccpp_prebuild.py failed, check ${BUILD_DIR}/ccpp_prebuild.log"
  exit 1
fi
set -e

source ${BUILD_DIR}/FV3/ccpp/physics/CCPP_SCHEMES.sh
source ${BUILD_DIR}/FV3/ccpp/physics/CCPP_CAPS.sh
source ${BUILD_DIR}/FV3/ccpp/physics/CCPP_STATIC_API.sh

CMAKE_FLAGS+=" -DCCPP=ON -DSTATIC=ON -DSUITES=${CCPP_SUITES} -DNETCDF_DIR=${NETCDF}"

cd ${BUILD_DIR}
cmake .. ${CMAKE_FLAGS}
make -j ${BUILD_JOBS:-4}
cp -v NEMS.exe ${MYDIR}/ufs_weather_model
