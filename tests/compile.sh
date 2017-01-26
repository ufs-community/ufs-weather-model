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
readonly BUILD_NR=$4

set +x
hostname

echo "Compiling ${MAKE_OPT}"
cd ${PATHTR}

./configure ${MACHINE_ID}
cp conf/configure.fv3 $PATHTR/../NEMS/src/conf
cat $PATHTR/../NEMS/src/conf/configure.nems.NUOPC conf/configure.fv3 >$PATHTR/../NEMS/src/conf/configure.nems
cat << EOF >> $PATHTR/../NEMS/src/conf/configure.nems
  FFLAGS += \$(ESMF_INC)
  CPPFLAGS += -traditional
  EXTLIBS = \$(NCEPLIBS) \$(ESMF_LIB) \$(LDFLAGS)
EOF

set +x
if [[ $MACHINE_ID = wcoss ]]; then
  source /usrx/local/Modules/default/init/sh
elif [[ $MACHINE_ID = wcoss_cray ]]; then
  source /opt/modules/default/init/sh
elif [[ $MACHINE_ID = theia ]]; then
  source /apps/lmod/lmod/init/sh
fi
source conf/modules.fv3
module list
set -x

gmake clean

gmake ${MAKE_OPT} -j 8 nemsinstall

cd $PATHTR/../NEMS/src
touch conf/externals.nems
cp ESMFVersionDefine_ESMF_NUOPC.h ESMFVersionDefine.h
export COMP=FV3                     
export COMP_SRCDIR=$PATHTR           
export COMP_BINDIR=$PATHTR/FV3_INSTALL
gmake nems COMP=,fv3, FV3_DIR=$COMP_BINDIR

cp ../exe/NEMS.x ${PATHTR}/../tests/fv3_${BUILD_NR}.exe
cp ${PATHTR}/conf/modules.fv3 ${PATHTR}/../tests/modules.fv3_${BUILD_NR}

gmake clean
cd $PATHTR
gmake cleanall

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling ${MAKE_OPT} finished"
