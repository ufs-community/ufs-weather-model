#!/bin/bash
set -eux

echo "PID=$$"
SECONDS=0

trap '[ "$?" -eq 0 ] || write_fail_test' EXIT
trap 'echo "run_compile.sh interrupted PID=$$"; cleanup' INT
trap 'echo "run_compile.sh terminated PID=$$";  cleanup' TERM

cleanup() {
  [[ $ROCOTO = 'false' ]] && interrupt_job
  trap 0
  exit
}

write_fail_test() {
  if [[ ${UNIT_TEST} == true ]]; then
    echo compile ${COMPILE_NR} >> $PATHRT/fail_unit_test
  else
    echo "compile_${COMPILE_NR} failed in run_compile" >> $PATHRT/fail_test
  fi
  exit 1
}

if [[ $# != 4 ]]; then
  echo "Usage: $0 PATHRT RUNDIR_ROOT MAKE_OPT COMPILE_NR"
  exit 1
fi

export PATHRT=$1
export RUNDIR_ROOT=$2
export MAKE_OPT=$3
export COMPILE_NR=$4

cd ${PATHRT}

[[ -e ${RUNDIR_ROOT}/compile_${COMPILE_NR}.env ]] && source ${RUNDIR_ROOT}/compile_${COMPILE_NR}.env
source default_vars.sh


export TEST_NAME=compile
export TEST_NR=${COMPILE_NR}
export JBNME="compile_${COMPILE_NR}"
export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}_${TEST_NR}

echo -n "${JBNME}, $( date +%s )," > ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

source rt_utils.sh
source atparse.bash

mkdir -p ${RUNDIR}
cd $RUNDIR

if [[ $SCHEDULER = 'slurm' ]]; then
  atparse < $PATHRT/fv3_conf/compile_slurm.IN > job_card
elif [[ $SCHEDULER = 'lsf' ]]; then
  atparse < $PATHRT/fv3_conf/compile_bsub.IN > job_card
elif [[ $SCHEDULER = 'pbs' ]]; then
  atparse < $PATHRT/fv3_conf/compile_qsub.IN > job_card
fi

################################################################################
# Submit compile job
################################################################################

if [[ $ROCOTO = 'false' ]]; then
  submit_and_wait job_card
else
  chmod u+x job_card
  ./job_card
fi

ls -l ${PATHTR}/tests/fv3_${COMPILE_NR}.exe

cp ${RUNDIR}/compile_*_time.log ${LOG_DIR}
cat ${RUNDIR}/job_timestamp.txt >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt
################################################################################
# End compile job
################################################################################

echo " $( date +%s )" >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compile ${COMPILE_NR}"
