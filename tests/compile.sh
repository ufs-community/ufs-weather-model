#!/bin/bash
set -eux

SECONDS=0

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 PATHTR BUILD_TARGET [ MAKE_OPT [ BUILD_NR ] [ clean_before ] [ clean_after ]  ]"
  echo Valid BUILD_TARGETs:
  echo $( ls -1 ../conf/configure.fv3.* | sed s,.*fv3\.,,g ) | fold -sw72
  exit 1
fi

# ----------------------------------------------------------------------
# Parse arguments.

readonly PATHTR=$1
readonly BUILD_TARGET=$2
MAKE_OPT=${3:-}
readonly BUILD_NAME=fv3${4:+_$4}

readonly clean_before=${5:-YES}
readonly clean_after=${6:-YES}

hostname

# ----------------------------------------------------------------------

echo "Compiling ${MAKE_OPT} into $BUILD_NAME.exe on $BUILD_TARGET"

# ----------------------------------------------------------------------
# Make sure we have a "make" and reasonable threads.

gnu_make=gmake
if ( ! which $gnu_make ) ; then
    echo WARNING: Cannot find gmake in \$PATH.  I will use \"make\" instead.
    gnu_make=make
    if ( ! $gnu_make --version 2>&1 | grep -i gnu > /dev/null 2>&1 ) ; then
       echo WARNING: The build system requires GNU Make. Things may break.
    fi
fi

if [[ $BUILD_TARGET == cheyenne.* || $BUILD_TARGET == stampede.* ]] ; then
    MAKE_THREADS=${MAKE_THREADS:-3}
fi

MAKE_THREADS=${MAKE_THREADS:-8}

if [[ "$MAKE_THREADS" -gt 1 ]] ; then
    echo Using \$MAKE_THREADS=$MAKE_THREADS threads to build FV3 and FMS.
    echo Consider reducing \$MAKE_THREADS if you hit memory or process limits.
    gnu_make="$gnu_make -j $MAKE_THREADS"
fi

# ----------------------------------------------------------------------
# Configure NEMS and components

# Configure NEMS
cd "$PATHTR/../NEMS"

COMPONENTS=""

if [[ "${MAKE_OPT}" == *"FV3=Y"* ]]; then
  COMPONENTS="FMS,FV3,$COMPONENTS"
fi

if [[ "${MAKE_OPT}" == *"CCPP=Y"* ]]; then
  COMPONENTS="CCPP,$COMPONENTS"
fi

if [[ "${MAKE_OPT}" == *"WW3=Y"* ]]; then
  COMPONENTS="WW3,$COMPONENTS"
fi

if [[ "${MAKE_OPT}" == *"HYCOM=Y"* ]]; then
  COMPONENTS="HYCOM,$COMPONENTS"
fi

if [[ "${MAKE_OPT}" == *"CMEPS=Y"* ]]; then
  COMPONENTS="CMEPS,$COMPONENTS"
fi

if [[ "${MAKE_OPT}" == *"CDEPS=Y"* ]]; then
  COMPONENTS="CDEPS,$COMPONENTS"
fi

if [[ "${MAKE_OPT}" == *"CDEPS_DATM=Y"* ]]; then
  COMPONENTS="CDEPS_DATM,$COMPONENTS"
fi

if [[ "${MAKE_OPT}" == *"CDEPS_DOCN=Y"* ]]; then
  COMPONENTS="CDEPS_DOCN,$COMPONENTS"
fi

# Make variables:
#   COMPONENTS = list of components to build
#   BUILD_ENV = theia.intel, wcoss_dell_p3, etc.
#   FV3_MAKEOPT = build options to send to FV3, CCPP, and FMS
#   TEST_BUILD_NAME = requests copying of modules.nems and
#      NEMS.x into the tests/ directory using the given build name.

# FIXME: add -j $MAKE_THREADS once FV3 bug is fixed

# Pass DEBUG or REPRO flags to NEMS
if [[ "${MAKE_OPT}" == *"DEBUG=Y"* && "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
  echo "ERROR in compile.sh: options DEBUG=Y and REPRO=Y are mutually exclusive"
  exit 1
elif [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
  NEMS_BUILDOPT="DEBUG=Y"
elif [[ "${MAKE_OPT}" == *"REPRO=Y"* ]]; then
  NEMS_BUILDOPT="REPRO=Y"
else
  NEMS_BUILDOPT=""
fi

if [ $clean_before = YES ] ; then
  $gnu_make -k COMPONENTS="$COMPONENTS" TEST_BUILD_NAME="$BUILD_NAME" \
           BUILD_ENV="$BUILD_TARGET" FV3_MAKEOPT="$MAKE_OPT" \
           NEMS_BUILDOPT="$NEMS_BUILDOPT" distclean
fi

# FIXME - create CCPP include directory before building FMS to avoid
# gfortran warnings of non-existent include directory (adding
# -Wno-missing-include-dirs) to the GNU compiler flags does not work,
# see also https://gcc.gnu.org/bugzilla/show_bug.cgi?id=55534);
# this line can be removed once FMS becomes a pre-installed library
if [[ "${MAKE_OPT}" == *"CCPP=Y"* ]]; then
  mkdir -p $PATHTR/ccpp/include
fi

  $gnu_make -k COMPONENTS="$COMPONENTS" TEST_BUILD_NAME="$BUILD_NAME" \
           BUILD_ENV="$BUILD_TARGET" FV3_MAKEOPT="$MAKE_OPT" \
           NEMS_BUILDOPT="$NEMS_BUILDOPT" build

if [ $clean_after = YES ] ; then
  $gnu_make -k COMPONENTS="$COMPONENTS" TEST_BUILD_NAME="$BUILD_NAME" \
           BUILD_ENV="$BUILD_TARGET" FV3_MAKEOPT="$MAKE_OPT" \
           NEMS_BUILDOPT="$NEMS_BUILDOPT" clean
fi

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${MAKE_OPT} finished"

