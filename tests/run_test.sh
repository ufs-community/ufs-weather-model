#!/bin/bash
set -eux

write_fail_test() {
  if [[ ${UNIT_TEST} == true ]]; then
    echo ${TEST_NR} $TEST_NAME >> $PATHRT/fail_unit_test
  fi
}

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

[[ -e ${RUNDIR_ROOT}/run_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/run_test_${TEST_NR}.env
source default_vars.sh
source tests/$TEST_NAME
[[ -e ${RUNDIR_ROOT}/unit_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/unit_test_${TEST_NR}.env

# Save original CNTL_DIR name as INPUT_DIR for regression
# tests that try to copy input data from CNTL_DIR
export INPUT_DIR=${CNTL_DIR}
# Append RT_SUFFIX to RUNDIR, and BL_SUFFIX to CNTL_DIR
export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}${RT_SUFFIX}
export CNTL_DIR=${CNTL_DIR}${BL_SUFFIX}

JBNME=$(basename $RUNDIR_ROOT)_${TEST_NR}
export JBNME

export FV3X=fv3_${COMPILE_NR}.exe

UNIT_TEST=${UNIT_TEST:-false}
if [[ ${UNIT_TEST} == false ]]; then
  REGRESSIONTEST_LOG=${LOG_DIR}/rt_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log
else
  REGRESSIONTEST_LOG=${LOG_DIR}/ut_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log
fi
export REGRESSIONTEST_LOG

# Submit the actual test run script
echo "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR}"
trap 'echo "run_test.sh: Test ${TEST_NAME} killed"; kill $(jobs -p); wait; trap 0; exit' 1 2 3 4 5 6 7 8 10 12 13 15
trap '[ "$?" -eq 0 ] || write_fail_test' EXIT
#trap 'echo "run_test.sh: Test ${TEST_NAME} error"; echo "${TEST_NAME}" >> ${PATHRT}/fail_test; trap 0; exit' ERR

if [[ $CI_TEST = true ]]; then
  ./${RUN_SCRIPT} >${RUNDIR_ROOT}/${TEST_NAME}${RT_SUFFIX}.log 2> >(tee -a ${RUNDIR_ROOT}/${TEST_NAME}${RT_SUFFIX}.log >&3)
else
  ./${RUN_SCRIPT} > ${RUNDIR_ROOT}/${TEST_NAME}${RT_SUFFIX}.log 2>&1
fi

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Test ${TEST_NAME}"
