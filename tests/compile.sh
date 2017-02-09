#!/bin/bash
set -eu

SECONDS=0

if [[ $# != 4 ]]; then
  echo "Usage: $0 PATHTR MACHINE_ID MAKE_OPT BUILD_NR"
  exit 1
fi

readonly PATHTR=$1
readonly MACHINE_ID=$2
readonly MAKE_OPT=$3
readonly BUILD_NAME=fv3${4:+_$4}

readonly PATHNEMS="$PATHTR/../NEMS"

set +x
hostname

echo "Compiling ${MAKE_OPT}"

# Configure FV3
cd ${PATHTR}
./configure ${MACHINE_ID}

# Configure NEMS
cd $PATHNEMS
src/configure conf/configure.fv3.$MACHINE_ID \
              modulefiles/$MACHINE_ID/fv3

# Load the "module" command, purge modules, and load modules
source $PATHNEMS/src/conf/modules.nems.sh
module list

# Build FV3
cd $PATHTR
gmake clean
gmake ${MAKE_OPT} -j 8 nemsinstall

# Build NEMS
cd $PATHNEMS/src
export COMP=FV3                     
export COMP_SRCDIR=$PATHTR           
export COMP_BINDIR=$PATHTR/FV3_INSTALL
gmake nems COMP=,fv3, FV3_DIR=$COMP_BINDIR

cp ../exe/NEMS.x ${PATHTR}/../tests/$BUILD_NAME.exe
cp $PATHNEMS/conf/modules.nems ${PATHTR}/../tests/modules.$BUILD_NAME

gmake clean
cd $PATHTR
gmake cleanall

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${MAKE_OPT} finished"
