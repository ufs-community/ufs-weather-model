#!/bin/bash

set +x
set -eu

# List of valid/tested machines
VALID_MACHINES=( theia.intel theia.gnu theia.pgi cheyenne.intel cheyenne.gnu cheyenne.pgi macosx.gnu )

###################################################################################################

function usage   {
  echo "Usage: "
  echo "build_ccpp.sh MACHINE_ID [ MAKE_OPT ] [ clean_before ] [ clean_after ]"
  echo "    Where: MACHINE      [required] can be : ${VALID_MACHINES[@]}"
  echo "           CCPP_DIR     [required] is the target installation directory for CCPP"
  echo "           MAKE_OPT     [optional] can be any of the NEMSfv3gfs MAKE_OPT options; used:"
  echo "                                   DEBUG=Y/N  (default N)"
  echo "                                   OPENMP=Y/N (default N)"
  echo "           clean_before [optional] can be 'YES' (default) or 'NO'"
  echo "           clean_after  [optional] can be 'YES' (default) or 'NO'"
  exit 1
}

function checkvalid {
# Ensure value ($2) of variable ($1) is contained in list of validvalues ($3)
  if [ $# -lt 3 ]; then
    echo $FUNCNAME requires at least 3 arguments: stopping
    exit 1
  fi

  var_name=$1 && shift
  input_val=$1 && shift
  valid_vars=($@)

  for x in ${valid_vars[@]}; do
    if [ "$input_val" == "$x" ]; then
      echo "${var_name}=${input_val} is valid."
      return
    fi
  done
  echo "ERROR: ${var_name}=${input_val} is invalid."
  usage
  exit 1
}

function trim {
    local var="$1"
    # remove leading whitespace characters
    var="${var#"${var%%[![:space:]]*}"}"
    # remove trailing whitespace characters
    var="${var%"${var##*[![:space:]]}"}"   
    echo -n "$var"
}

###################################################################################################

# Check and process command line arguments

if [[ $# -lt 2 ]]; then usage; fi

readonly MACHINE_ID=$1
readonly CCPP_DIR=$2
readonly MAKE_OPT=${3:-}
readonly clean_before=${4:-YES}
readonly clean_after=${5:-YES}

checkvalid MACHINE_ID $MACHINE_ID ${VALID_MACHINES[@]}

# Generate CCPP cmake flags from MAKE_OPT

CCPP_CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=${CCPP_DIR} -DNCEPLIBS_DIR=${NCEPLIBS_DIR} -DMPI=ON"
CCPP_MAKE_FLAGS=""
if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
  CCPP_MAKE_FLAGS="${CCPP_MAKE_FLAGS} VERBOSE=1"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Release"
fi
if [[ "${MAKE_OPT}" == *"OPENMP=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DOPENMP=ON"
else
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DOPENMP=OFF"
fi
CCPP_CMAKE_FLAGS=$(trim "${CCPP_CMAKE_FLAGS}")
CCPP_MAKE_FLAGS=$(trim "${CCPP_MAKE_FLAGS}")

# Build and install CCPP

echo "Building CCPP with options '${CCPP_CMAKE_FLAGS}' ..."
PATH_CCPP=${PWD}
PATH_CCPP_BUILD=${PWD}/build
if [ $clean_before = YES ]; then
    rm -fr ${PATH_CCPP_BUILD}
fi
mkdir -p ${PATH_CCPP_BUILD}
cd ${PATH_CCPP_BUILD}
cmake ${CCPP_CMAKE_FLAGS} ${PATH_CCPP}
make ${CCPP_MAKE_FLAGS}
make ${CCPP_MAKE_FLAGS} install

if [ $clean_after = YES ]; then
    rm -fr ${PATH_CCPP_BUILD}
fi
