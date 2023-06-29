#!/bin/bash
set -eux

function trim {
    local var="$1"
    # remove leading whitespace characters
    var="${var#"${var%%[![:space:]]*}"}"
    # remove trailing whitespace characters
    var="${var%"${var##*[![:space:]]}"}"
    echo -n "$var"
}

SECONDS=0

if [[ $(uname -s) == Darwin ]]; then
  readonly MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
else
  readonly MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi

# ----------------------------------------------------------------------
# Parse arguments.

readonly ARGC=$#

if [[ $ARGC -lt 2 ]]; then
  echo "Usage: $0 MACHINE_ID [ MAKE_OPT [ COMPILE_NR ] [ RT_COMPILER ] [ clean_before ] [ clean_after ] ]"
  echo Valid MACHINE_IDs:
  echo $( ls -1 ../cmake/configure_* | sed s:.*configure_::g | sed s:\.cmake:: ) | fold -sw72
  exit 1
else
  MACHINE_ID=$1
  MAKE_OPT=${2:-}
  COMPILE_NR=${3:+_$3}
  RT_COMPILER=${4:-intel}
  clean_before=${5:-YES}
  clean_after=${6:-YES}
fi

BUILD_NAME=fv3${COMPILE_NR}

PATHTR=${PATHTR:-$( cd ${MYDIR}/.. && pwd )}
BUILD_DIR=$(pwd)/build_${BUILD_NAME}

# ----------------------------------------------------------------------
# Make sure we have reasonable number of threads.

if [[ $MACHINE_ID == cheyenne ]]; then
    BUILD_JOBS=${BUILD_JOBS:-3}
fi

BUILD_JOBS=${BUILD_JOBS:-8}

hostname

set +x
if [[ $MACHINE_ID == macosx ]] || [[ $MACHINE_ID == linux ]]; then
  source $PATHTR/modulefiles/ufs_${MACHINE_ID}.${RT_COMPILER}
else
  # Activate lua environment for gaea
  if [[ $MACHINE_ID == gaea ]]; then
    source /lustre/f2/dev/role.epic/contrib/Lmod_init.sh
  fi
  # Load fv3 module
  module purge
  module use $PATHTR/modulefiles
  modulefile="ufs_${MACHINE_ID}.${RT_COMPILER}"
  module load $modulefile
  module list
fi
set -x

echo "Compiling ${MAKE_OPT} into $BUILD_NAME.exe on $MACHINE_ID"

# set CMAKE_FLAGS based on $MAKE_OPT

CMAKE_FLAGS=$MAKE_OPT
CMAKE_FLAGS+=" -DMPI=ON"

if [[ ${MAKE_OPT} == *-DDEBUG=ON* ]]; then
  CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Debug"
else
  CMAKE_FLAGS+=" -DCMAKE_BUILD_TYPE=Release"
  if [[ ${MACHINE_ID} == jet ]] && [[ ${RT_COMPILER} == intel ]]; then
    CMAKE_FLAGS+=" -DSIMDMULTIARCH=ON"
  fi
fi

# Check if suites argument is provided or not
set +ex
TEST=$( echo $MAKE_OPT | grep -e "-DCCPP_SUITES=" )
if [[ $? -eq 0 ]]; then
  SUITES=$( echo $MAKE_OPT | sed 's/.*-DCCPP_SUITES=//' | sed 's/ .*//' )
  echo "Compiling suites ${SUITES}"
fi
set -ex

# Valid applications

if [[ "${MAKE_OPT}" == *"-DAPP=S2S"* ]]; then
    CMAKE_FLAGS+=" -DMOM6SOLO=ON"
fi

if [[ "${MAKE_OPT}" == *"-DAPP=NG-GODAS"* ]]; then
    CMAKE_FLAGS+=" -DMOM6SOLO=ON"
fi

CMAKE_FLAGS=$(trim "${CMAKE_FLAGS}")
echo "CMAKE_FLAGS = ${CMAKE_FLAGS}"

if [ $clean_before = YES ] ; then
  rm -rf ${BUILD_DIR}
fi

export BUILD_VERBOSE=1
export BUILD_DIR
export BUILD_JOBS
export CMAKE_FLAGS

bash -x ${PATHTR}/build.sh

mv ${BUILD_DIR}/ufs_model ${PATHTR}/tests/${BUILD_NAME}.exe
if [[ $MACHINE_ID == linux ]]; then
  cp ${PATHTR}/modulefiles/ufs_${MACHINE_ID}.${RT_COMPILER}       ${PATHTR}/tests/modules.${BUILD_NAME}
else
  cp ${PATHTR}/modulefiles/ufs_${MACHINE_ID}.${RT_COMPILER}.lua       ${PATHTR}/tests/modules.${BUILD_NAME}.lua
fi

if [ $clean_after = YES ] ; then
  rm -rf ${BUILD_DIR}
fi

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${CMAKE_FLAGS} finished"
echo "Compile ${COMPILE_NR/#_} elapsed time $elapsed seconds. ${CMAKE_FLAGS}" > compile${COMPILE_NR}_time.log
