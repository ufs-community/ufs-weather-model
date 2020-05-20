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
CMAKE_FLAGS+=" -DCCPP=ON -DSUITES=${CCPP_SUITES} -DNETCDF_DIR=${NETCDF}"
CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Debug -DDEBUG=YES"

#NCEPLIBS=/scratch2/NCEPDEV/fv3-cam/Dusan.Jovic/nceplibs/local
#CMAKE_FLAGS+=" -DCMAKE_PREFIX_PATH=${NCEPLIBS}/w3emc_2.5.0;${NCEPLIBS}/w3nco_2.1.0;${NCEPLIBS}/sigio_2.2.0;${NCEPLIBS}/nemsio_2.3.0;${NCEPLIBS}/crtm_2.3.0;${NCEPLIBS}/g2_3.2.0;${NCEPLIBS}/g2tmpl_1.7.0"

#CMAKE_FLAGS+=" -DCMAKE_PREFIX_PATH=/scratch2/NCEPDEV/stmp1/Rahul.Mahajan/opt/intel-18.0.5.274/impi-2018.0.4/nceplibs-ufs/1.0.0"

CMAKE_FLAGS+=" -DCMAKE_PREFIX_PATH=/scratch2/NCEPDEV/fv3-cam/Dusan.Jovic/nceplibs/install"

cd ${BUILD_DIR}
cmake .. ${CMAKE_FLAGS}
make -j ${BUILD_JOBS:-4}
cp NEMS.exe ${MYDIR}/ufs_weather_model
