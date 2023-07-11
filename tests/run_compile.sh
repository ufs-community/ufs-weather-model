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
  echo "compile_${COMPILE_NR} failed in run_compile" >> $PATHRT/fail_compile_${COMPILE_NR}
  exit 1
}

remove_fail_test() {
    echo "Removing test failure flag file for compile_${COMPILE_NR}"
    rm -f $PATHRT/fail_compile_${COMPILE_NR}
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
remove_fail_test

[[ -e ${RUNDIR_ROOT}/compile_${COMPILE_NR}.env ]] && source ${RUNDIR_ROOT}/compile_${COMPILE_NR}.env
source default_vars.sh
[[ -e ${RUNDIR_ROOT}/compile_${COMPILE_NR}.env ]] && source ${RUNDIR_ROOT}/compile_${COMPILE_NR}.env



export TEST_NAME=compile
export TEST_NR=${COMPILE_NR}
export JBNME="compile_${COMPILE_NR}"
export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}_${TEST_NR}

echo -n "${JBNME}, $( date +%s )," > ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

export RT_LOG=${LOG_DIR}/compile_${TEST_NR}.log

source rt_utils.sh
source atparse.bash

rm -rf ${RUNDIR}
mkdir -p ${RUNDIR}
cd $RUNDIR

if [[ $SCHEDULER = 'pbs' ]]; then
  if [[ -e $PATHRT/fv3_conf/compile_qsub.IN_${MACHINE_ID} ]]; then 
    atparse < $PATHRT/fv3_conf/compile_qsub.IN_${MACHINE_ID} > job_card
  else
    echo "Looking for fv3_conf/compile_qsub.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
elif [[ $SCHEDULER = 'slurm' ]]; then
  if [[ -e $PATHRT/fv3_conf/compile_slurm.IN_${MACHINE_ID} ]]; then
    atparse < $PATHRT/fv3_conf/compile_slurm.IN_${MACHINE_ID} > job_card
  else
    echo "Looking for fv3_conf/compile_slurm.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
elif [[ $SCHEDULER = 'lsf' ]]; then
  if [[ -e $PATHRT/fv3_conf/compile_bsub.IN_${MACHINE_ID} ]]; then
    atparse < $PATHRT/fv3_conf/compile_bsub.IN_${MACHINE_ID} > job_card
  else
    echo "Looking for fv3_conf/compile_bsub.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
fi

################################################################################
# Submit compile job
################################################################################

if [[ $ROCOTO = 'false' ]]; then
  submit_and_wait job_card
else
  chmod u+x job_card
  ( ./job_card 2>&1 1>&3 3>&- | tee err ) 3>&1 1>&2 | tee out
  # The above shell redirection copies stdout to "out" and stderr to "err"
  # while still sending them to stdout and stderr. It does this without
  # relying on bash-specific extensions or non-standard OS features.
fi

ls -l ${PATHTR}/tests/fv3_${COMPILE_NR}.exe

cp ${RUNDIR}/compile_*_time.log ${LOG_DIR}
cat ${RUNDIR}/job_timestamp.txt >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

remove_fail_test

################################################################################
# End compile job
################################################################################

echo " $( date +%s ), 1" >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Compile ${COMPILE_NR}"
