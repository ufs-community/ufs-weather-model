#!/bin/bash

SECONDS=0

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 PATHTR BUILD_TARGET [ MAKE_OPT [ BUILD_NR ] [ clean_before ] [ clean_after ]  ]"
  echo Valid BUILD_TARGETs:
  echo $( ls -1 ../conf/configure.fv3.* | sed s,.*fv3\.,,g ) | fold -sw72
  exit 1
fi

set -xeu

readonly PATHTR=$1
readonly BUILD_TARGET=$2
readonly MAKE_OPT=${3:-}
#readonly BUILD_NAME=${4:-}
readonly BUILD_NAME=fv3${4:+_$4}

clean_before=${5:-YES}
clean_after=${6:-YES}

readonly PATHNEMS="$PATHTR/../NEMS"

#set +x
hostname

echo "Compiling ${MAKE_OPT} into $BUILD_NAME.exe on $BUILD_TARGET"

# Configure NEMS
cd $PATHNEMS
# Current Intel 15 compiler is incompatible with latest features needed
# for CCPP (dynamic loading, ISO_C bindings etc.) and causes segfaults
# when trying to load the shared library libccpp.so when launching fv3.exe
if [[ "${MAKE_OPT}" == *"CCPP=Y"* && "${BUILD_TARGET}" == "theia.intel" ]]; then
  src/configure configure.fv3.$BUILD_TARGET $BUILD_TARGET/fv3.intel-18.0.1.163
else
  src/configure configure.fv3.$BUILD_TARGET $BUILD_TARGET/fv3
fi

# Load the "module" command, purge modules, and load modules
source $PATHNEMS/src/conf/modules.nems.sh
if [[ $BUILD_TARGET != "macosx.gnu" ]]; then
  module list
fi

# Build CCPP
if [[ "${MAKE_OPT}" == *"CCPP=Y"* ]]; then
  export readonly PATH_CCPP="$PATHTR/../ccpp"
  export readonly NEMS_CCPP_CPPFLAGS="-DCCPP"
  export readonly NEMS_CCPP_LDFLAGS="-L${PATH_CCPP}/lib -lccpp"
  # Run ccpp_prebuild.py from the top-level directory before building the CCPP framework and physics
  cd ${PATHTR}/..
  if [[ "${MAKE_OPT}" == *"DEBUG=Y"* ]]; then
    ./ccpp/framework/scripts/ccpp_prebuild.py --model=FV3v1 --debug
  else
    ./ccpp/framework/scripts/ccpp_prebuild.py --model=FV3v1
  fi
  # Build CCPP framework and physics
  cd ${PATH_CCPP}
  ./build_ccpp.sh ${BUILD_TARGET} ${PATH_CCPP} "${MAKE_OPT}" ${clean_before} ${clean_after}
fi

# Copy configuration to FV3:
cp -fp $PATHNEMS/src/conf/configure.nems $PATHTR/conf/configure.fv3
cp -fp $PATHNEMS/src/conf/modules.nems $PATHTR/conf/modules.fv3

# Build FV3
cd $PATHTR
if [[ $BUILD_TARGET == "macosx.gnu" ]]; then
  if [ $clean_before = YES ] ; then make clean ; fi
  make ${MAKE_OPT} -j 8 nemsinstall
else
  if [ $clean_before = YES ] ; then gmake clean ; fi
  # CISL kills your shell if using too much resources on the login node
  if [[ $BUILD_TARGET == cheyenne.* ]]; then
    gmake ${MAKE_OPT} -j 3 nemsinstall
  else
    gmake ${MAKE_OPT} -j 8 nemsinstall
  fi
fi

#Build NEMS
cd $PATHNEMS/src
export COMP=FV3
export COMP_SRCDIR=$PATHTR
export COMP_BINDIR=$PATHTR/FV3_INSTALL
if [[ $BUILD_TARGET == "macosx.gnu" ]]; then
  make nems COMP=,fv3, FV3_DIR=$PATHTR ${MAKE_OPT}
else
  gmake nems COMP=,fv3, FV3_DIR=$PATHTR ${MAKE_OPT}
fi

cp ../exe/NEMS.x ${PATHTR}/../tests/$BUILD_NAME.exe
cp $PATHNEMS/src/conf/modules.nems ${PATHTR}/../tests/modules.$BUILD_NAME

if [ $clean_after = YES ] ; then
  if [[ $BUILD_TARGET == "macosx.gnu" ]]; then
    make clean
    cd $PATHTR
    make cleanall
  else
    gmake clean
    cd $PATHTR
    gmake cleanall
  fi
# A few things that "cleanall" doesn't clean:
 rm -rf FV3_INSTALL
 rm -rf nems_dir
 find $COMP_SRCDIR/fms/ -name '*.o' -o -name '*.mod' | xargs rm -f
fi

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${MAKE_OPT} finished"

