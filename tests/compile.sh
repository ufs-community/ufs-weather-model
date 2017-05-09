#!/bin/bash
set -xeu

SECONDS=0

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 PATHTR MACHINE_ID [ MAKE_OPT [ BUILD_NR ] ]"
  exit 1
fi

readonly PATHTR=$1
readonly MACHINE_ID=$2
readonly MAKE_OPT=${3:- }
readonly BUILD_NAME=fv3${4:+_$4}

readonly PATHNEMS="$PATHTR/../NEMS"

#set +x
hostname

echo "Compiling ${MAKE_OPT} into $BUILD_NAME.exe on $MACHINE_ID"

# Configure NEMS
cd $PATHNEMS
src/configure configure.fv3.$MACHINE_ID \
              $MACHINE_ID/fv3

# Load the "module" command, purge modules, and load modules
source $PATHNEMS/src/conf/modules.nems.sh
module list

# Copy configuration to FV3:
cp -fp $PATHNEMS/src/conf/configure.nems $PATHTR/conf/configure.fv3
cp -fp $PATHNEMS/src/conf/modules.nems $PATHTR/conf/modules.fv3

# Build FV3
cd $PATHTR
gmake clean
gmake ${MAKE_OPT} -j 8 nemsinstall

# Build NEMS
cd $PATHNEMS/src
export COMP=FV3                     
export COMP_SRCDIR=$PATHTR           
export COMP_BINDIR=$PATHTR/FV3_INSTALL
gmake nems COMP=,fv3, FV3_DIR=$PATHTR

cp ../exe/NEMS.x ${PATHTR}/../tests/$BUILD_NAME.exe
cp $PATHNEMS/src/conf/modules.nems ${PATHTR}/../tests/modules.$BUILD_NAME

gmake clean
cd $PATHTR
gmake cleanall

# A few things that "cleanall" doesn't clean:
rm -rf FV3_INSTALL
rm -rf nems_dir
find $COMP_SRCDIR/fms/ -name '*.o' -o -name '*.mod' | xargs rm -f

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${MAKE_OPT} finished"
