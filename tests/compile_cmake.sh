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

readonly MYDIR=$( dirname $(readlink -f $0) )

# ----------------------------------------------------------------------
# Parse arguments.

readonly ARGC=$#

cd ${MYDIR}

if [[ $ARGC -eq 0 ]]; then
  COMPILER=intel
  . detect_machine.sh
  PATHTR=$(cd -P ${MYDIR}/.. && pwd)
  MAKE_OPT=''
  BUILD_NAME=fv3
  clean_before=YES
  clean_after=YES
elif [[ $ARGC -lt 2 ]]; then
  echo "Usage: $0 PATHTR MACHINE_ID [ MAKE_OPT [ BUILD_NR ] [ clean_before ] [ clean_after ]  ]"
  echo Valid MACHINE_IDs:
  echo $( ls -1 ../conf/configure.fv3.* | sed s,.*fv3\.,,g ) | fold -sw72
  exit 1
else
  PATHTR=$1
  MACHINE_ID=$2
  MAKE_OPT=${3:-}
  BUILD_NAME=fv3${4:+_$4}
  clean_before=${5:-YES}
  clean_after=${6:-YES}
fi
BUILD_DIR=build_${BUILD_NAME}

# ----------------------------------------------------------------------
# Make sure we have reasonable number of threads.

if [[ $MACHINE_ID == cheyenne.* ]] ; then
    MAKE_THREADS=${MAKE_THREADS:-3}
elif [[ $MACHINE_ID == wcoss_dell_p3 ]] ; then
    MAKE_THREADS=${MAKE_THREADS:-4}
fi

MAKE_THREADS=${MAKE_THREADS:-8}

hostname

cd ${PATHTR}/tests

# ----------------------------------------------------------------------

echo "Compiling ${MAKE_OPT} into $BUILD_NAME.exe on $MACHINE_ID"

if [ $clean_before = YES ] ; then
  rm -rf ${BUILD_DIR}
fi

mkdir -p ${BUILD_DIR}

# set CCPP_CMAKE_FLAGS based on $MAKE_OPT

CCPP_CMAKE_FLAGS=""

if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DDEBUG=Y"
elif [[ "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DREPRO=Y"
fi

if [[ "${MAKE_OPT}" == *"32BIT=Y"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -D32BIT=Y"
fi

if [[ "${MAKE_OPT}" == *"OPENMP=N"* ]]; then
  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DOPENMP=OFF"
fi

if [[ "${MAKE_OPT}" == *"CCPP=Y"* ]]; then

  CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCCPP=ON -DMPI=ON"

  if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Debug"
  elif [[ "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Bitforbit"
  else
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DCMAKE_BUILD_TYPE=Release"
    if [[ "${MACHINE_ID}" == "jet.intel" ]]; then
      CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DSIMDMULTIARCH=ON"
    fi
  fi

  if [[ "${MAKE_OPT}" == *"32BIT=Y"* ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DDYN32=ON"
  else
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DDYN32=OFF"
  fi

  if [[ "${MAKE_OPT}" == *"STATIC=Y"* ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DSTATIC=ON"
  else
    echo "Error, dynamic CCPP build not supported"
    exit 1
  fi

  if [[ "${MAKE_OPT}" == *"MULTI_GASES=Y"* ]]; then
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DMULTI_GASES=ON"
  else
    CCPP_CMAKE_FLAGS="${CCPP_CMAKE_FLAGS} -DMULTI_GASES=OFF"
  fi

  (
    SUITES=$( echo $MAKE_OPT | sed 's/.* SUITES=//' | sed 's/ .*//' )
    cd ${PATHTR}
    ./FV3/ccpp/framework/scripts/ccpp_prebuild.py --config=FV3/ccpp/config/ccpp_prebuild_config.py --static --suites=${SUITES} --builddir=tests/${BUILD_DIR}/FV3
  )

  # Read list of schemes, caps, and static API
  source ${BUILD_DIR}/FV3/ccpp/physics/CCPP_SCHEMES.sh
  source ${BUILD_DIR}/FV3/ccpp/physics/CCPP_CAPS.sh
  source ${BUILD_DIR}/FV3/ccpp/physics/CCPP_STATIC_API.sh

 fi

CCPP_CMAKE_FLAGS=$(trim "${CCPP_CMAKE_FLAGS}")

(
  source $PATHTR/NEMS/src/conf/module-setup.sh.inc
  if [[ $MACHINE_ID == macosx.* ]] || [[ $MACHINE_ID == linux.* ]]; then
    source $PATHTR/modulefiles/${MACHINE_ID}/fv3
  else
    module use $PATHTR/modulefiles/${MACHINE_ID}
    module load fv3
    module list
  fi

  cd ${BUILD_DIR}

  cmake ${PATHTR} ${CCPP_CMAKE_FLAGS}
  make -j ${MAKE_THREADS}
  mv NEMS.exe ../${BUILD_NAME}.exe
  cp ${PATHTR}/modulefiles/${MACHINE_ID}/fv3 ../modules.${BUILD_NAME}
  cd ..
)

if [ $clean_after = YES ] ; then
  rm -rf ${BUILD_DIR}
fi

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${MAKE_OPT} finished"

