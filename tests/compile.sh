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
  echo "Usage: $0 MACHINE_ID [ MAKE_OPT [ BUILD_NR ] [ clean_before ] [ clean_after ] ]"
  echo Valid MACHINE_IDs:
  echo $( ls -1 ../cmake/configure_* | sed s:.*configure_::g | sed s:\.cmake:: ) | fold -sw72
  exit 1
else
  MACHINE_ID=$1
  MAKE_OPT=${2:-}
  BUILD_NAME=fv3${3:+_$3}
  clean_before=${4:-YES}
  clean_after=${5:-YES}
fi

PATHTR=${PATHTR:-$( cd ${MYDIR}/.. && pwd )}
BUILD_DIR=$(pwd)/build_${BUILD_NAME}

# ----------------------------------------------------------------------
# Make sure we have reasonable number of threads.

if [[ $MACHINE_ID == cheyenne.* ]] ; then
    BUILD_JOBS=${BUILD_JOBS:-3}
elif [[ $MACHINE_ID == wcoss_dell_p3 ]] ; then
    BUILD_JOBS=${BUILD_JOBS:-1}
fi

BUILD_JOBS=${BUILD_JOBS:-8}

hostname

set +x
if [[ $MACHINE_ID == macosx.* ]] || [[ $MACHINE_ID == linux.* ]]; then
  source $PATHTR/modulefiles/${MACHINE_ID}/fv3
else
  if [[ $MACHINE_ID == wcoss2 ]]; then
    source /apps/prod/lmodules/startLmod
  fi
  # Activate lua environment for gaea
  if [[ $MACHINE_ID == gaea.* ]] ; then
    source /lustre/f2/pdata/esrl/gsd/contrib/lua-5.1.4.9/init/init_lmod.sh
  fi
  # Load fv3 module
  module use $PATHTR/modulefiles/${MACHINE_ID}
  modulefile="fv3"
  if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
    [[ -f $PATHTR/modulefiles/${MACHINE_ID}/fv3_debug ]] && modulefile="fv3_debug"
  fi
  module load $modulefile
  module list
fi
set -x

echo "Compiling ${MAKE_OPT} into $BUILD_NAME.exe on $MACHINE_ID"

# set CMAKE_FLAGS based on $MAKE_OPT

CMAKE_FLAGS=''

if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
  CMAKE_FLAGS="${CMAKE_FLAGS} -DDEBUG=Y"
elif [[ "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
  CMAKE_FLAGS="${CMAKE_FLAGS} -DREPRO=Y"
fi

if [[ "${MAKE_OPT}" == *"32BIT=Y"* ]]; then
  CMAKE_FLAGS="${CMAKE_FLAGS} -D32BIT=Y"
fi

if [[ "${MAKE_OPT}" == *"OPENMP=N"* ]]; then
  CMAKE_FLAGS="${CMAKE_FLAGS} -DOPENMP=OFF"
fi

if [[ "${MAKE_OPT}" == *"MULTI_GASES=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DMULTI_GASES=ON"
else
    CMAKE_FLAGS="${CMAKE_FLAGS} -DMULTI_GASES=OFF"
fi

if [[ "${MAKE_OPT}" == *"CCPP=Y"* ]]; then

  # FIXME - create CCPP include directory before building FMS to avoid
  # gfortran warnings of non-existent include directory (adding
  # -Wno-missing-include-dirs) to the GNU compiler flags does not work,
  # see also https://gcc.gnu.org/bugzilla/show_bug.cgi?id=55534);
  # this line can be removed once FMS becomes a pre-installed library
  mkdir -p $PATHTR/FV3/ccpp/include
  # Similar for this directory, which apparently never gets populated
  mkdir -p $PATHTR/FMS/fms2_io/include

  CMAKE_FLAGS="${CMAKE_FLAGS} -DCCPP=ON -DMPI=ON"

  if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
  elif [[ "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Bitforbit"
  else
    CMAKE_FLAGS="${CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Release"
    if [[ "${MACHINE_ID}" == "jet.intel" ]]; then
      CMAKE_FLAGS="${CMAKE_FLAGS} -DSIMDMULTIARCH=ON"
    fi
  fi

  if [[ "${MAKE_OPT}" == *"32BIT=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DDYN32=ON"
  else
    CMAKE_FLAGS="${CMAKE_FLAGS} -DDYN32=OFF"
  fi

    # Check if suites argument is provided or not
  set +ex
  TEST=$( echo $MAKE_OPT | grep -e "SUITES=" )
  if [[ $? -eq 0 ]]; then
    CCPP_SUITES=$( echo $MAKE_OPT | sed 's/.* SUITES=//' | sed 's/ .*//' )
    echo "Compiling suites ${CCPP_SUITES}"
  fi
  set -ex

fi

if [[ "${MAKE_OPT}" == *"WW3=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DWW3=Y"
fi

if [[ "${MAKE_OPT}" == *"S2S=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DS2S=Y"
fi

if [[ "${MAKE_OPT}" == *"DATM=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DDATM=Y"
fi

if [[ "${MAKE_OPT}" == *"S2S=Y"* ]] || [[ ${MAKE_OPT} == *"DATM=Y"* ]]; then
    CMAKE_FLAGS="${CMAKE_FLAGS} -DMOM6SOLO=ON"
fi

CMAKE_FLAGS=$(trim "${CMAKE_FLAGS}")

if [ $clean_before = YES ] ; then
  rm -rf ${BUILD_DIR}
fi

export BUILD_VERBOSE=1
export BUILD_DIR
export BUILD_JOBS
export CCPP_SUITES
export CMAKE_FLAGS

bash -x ${PATHTR}/build.sh

mv ${BUILD_DIR}/ufs_model ${PATHTR}/tests/${BUILD_NAME}.exe
if [[ "${MAKE_OPT}" == "DEBUG=Y" ]]; then
  cp ${PATHTR}/modulefiles/${MACHINE_ID}/fv3_debug ${PATHTR}/tests/modules.${BUILD_NAME}
else
  cp ${PATHTR}/modulefiles/${MACHINE_ID}/fv3 ${PATHTR}/tests/modules.${BUILD_NAME}
fi

if [ $clean_after = YES ] ; then
  rm -rf ${BUILD_DIR}
fi

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${MAKE_OPT} finished"
