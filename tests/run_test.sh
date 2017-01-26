#!/bin/bash
set -eux

SECONDS=0

if [[ $# != 5 ]]; then
  echo "Usage: $0 PATHRT RUNDIR_ROOT TEST_NAME TEST_NR COMPILE_NR"
  exit 1
fi

export PATHRT=$1
export RUNDIR_ROOT=$2
export TEST_NAME=$3
export TEST_NR=$4
export COMPILE_NR=$5

cd ${PATHRT}

[[ -e run_test.env ]] && source run_test.env
source default_vars.sh
source tests/$TEST_NAME

export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}
JBNME=$(basename $RUNDIR_ROOT)_${TEST_NR}
export JBNME

export FV3X=fv3_${COMPILE_NR}.exe
export REGRESSIONTEST_LOG=${LOG_DIR}/rt_${TEST_NR}_${TEST_NAME}.log

# Submit the actual test run script
echo "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR}"
trap 'echo "run_test.sh: Test ${TEST_NAME} killed"; kill $(jobs -p); wait; trap 0; exit' 1 2 3 4 5 6 7 8 10 12 13 15
./${RUN_SCRIPT} > ${RUNDIR_ROOT}/${TEST_NAME}.log 2>&1 &
wait

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Test ${TEST_NAME}"
