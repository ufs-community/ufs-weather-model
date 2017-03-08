#!/bin/bash

set -xeu

SECONDS=0

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 PATHTR MACHINE_ID [ MAKE_OPT [ BUILD_NR ] ]"
  exit 1
fi

readonly PATHTR=$1
readonly APP=$2
readonly BUILD_NAME=fv3${3:+_$3}

hostname

echo "Compiling app $APP into $BUILD_NAME.exe"

cd ${PATHTR}/..

rm -f $PATHTR/../NEMS/exe/NEMS.x
rm -f $PATHTR/../NEMS/src/conf/modules.nems

./NEMS/NEMSAppBuilder app="$APP"

cp $PATHTR/../NEMS/exe/NEMS.x ${PATHTR}/../tests/$BUILD_NAME.exe
cp $PATHTR/../NEMS/src/conf/modules.nems ${PATHTR}/../tests/modules.$BUILD_NAME

cd $PATHTR/../NEMS/src
gmake clean
cd $PATHTR
gmake cleanall

# A few things that "cleanall" doesn't clean:
rm -rf FV3_INSTALL
rm -rf nems_dir

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compiling app ${APP} finished"
