#!/bin/bash
set -eux

echo kimmmmmmmmm2

echo "PID=$$"
SECONDS=0

trap '[ "$?" -eq 0 ] || write_fail_test' EXIT
trap 'echo "run_test.sh interrupted PID=$$"; cleanup' INT
trap 'echo "run_test.sh terminated PID=$$";  cleanup' TERM

cleanup() {
  [[ $ROCOTO = 'false' ]] && interrupt_job
  trap 0
  exit
}

write_fail_test() {
  if [[ ${OPNREQ_TEST} == true ]]; then
    echo "${TEST_NAME} ${TEST_NR} failed in run_test" >> $PATHRT/fail_opnreq_test_${TEST_NR}
  else
    echo "${TEST_NAME} ${TEST_NR} failed in run_test" >> $PATHRT/fail_test_${TEST_NR}
  fi
  exit 1
}

remove_fail_test() {
    echo "Removing test failure flag file for ${TEST_NAME} ${TEST_NR}"
    if [[ ${OPNREQ_TEST} == true ]] ; then
        rm -f $PATHRT/fail_opnreq_test_${TEST_NR}
    else
        rm -f $PATHRT/fail_test_${TEST_NR}
    fi
}

function compute_petbounds() {

  # each test MUST define ${COMPONENT}_tasks variable for all components it is using
  # and MUST NOT define those that it's not using or set the value to 0.

  # ATM is a special case since it is running on the sum of compute and io tasks.
  # CHM component and mediator are running on ATM compute tasks only.

  local n=0
  unset atm_petlist_bounds ocn_petlist_bounds ice_petlist_bounds wav_petlist_bounds chm_petlist_bounds med_petlist_bounds aqm_petlist_bounds

  # ATM
  ATM_io_tasks=${ATM_io_tasks:-0}
  if [[ $((ATM_compute_tasks + ATM_io_tasks)) -gt 0 ]]; then
     atm_petlist_bounds="${n} $((n + ATM_compute_tasks + ATM_io_tasks -1))"
     n=$((n + ATM_compute_tasks + ATM_io_tasks))
  fi

  # OCN
  if [[ ${OCN_tasks:-0} -gt 0 ]]; then
     ocn_petlist_bounds="${n} $((n + OCN_tasks - 1))"
     n=$((n + OCN_tasks))
  fi

  # ICE
  if [[ ${ICE_tasks:-0} -gt 0 ]]; then
     ice_petlist_bounds="${n} $((n + ICE_tasks - 1))"
     n=$((n + ICE_tasks))
  fi

  # WAV
  if [[ ${WAV_tasks:-0} -gt 0 ]]; then
     wav_petlist_bounds="${n} $((n + WAV_tasks - 1))"
     n=$((n + WAV_tasks))
  fi

  # CHM
  chm_petlist_bounds="0 $((ATM_compute_tasks - 1))"

  # MED
  med_petlist_bounds="0 $((ATM_compute_tasks - 1))"

  # AQM
  aqm_petlist_bounds="0 $((ATM_compute_tasks - 1))"

  UFS_tasks=${n}

  echo "ATM_petlist_bounds: ${atm_petlist_bounds:-}"
  echo "OCN_petlist_bounds: ${ocn_petlist_bounds:-}"
  echo "ICE_petlist_bounds: ${ice_petlist_bounds:-}"
  echo "WAV_petlist_bounds: ${wav_petlist_bounds:-}"
  echo "CHM_petlist_bounds: ${chm_petlist_bounds:-}"
  echo "MED_petlist_bounds: ${med_petlist_bounds:-}"
  echo "AQM_petlist_bounds: ${aqm_petlist_bounds:-}"
  echo "UFS_tasks         : ${UFS_tasks:-}"

}

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
OPNREQ_TEST=${OPNREQ_TEST:-false}
remove_fail_test

[[ -e ${RUNDIR_ROOT}/run_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/run_test_${TEST_NR}.env
source default_vars.sh
source tests/$TEST_NAME
[[ -e ${RUNDIR_ROOT}/opnreq_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/opnreq_test_${TEST_NR}.env

#----
# Save original CNTL_DIR name as INPUT_DIR for regression
# tests that try to copy input data from CNTL_DIR
export INPUT_DIR=${CNTL_DIR}
# Append RT_SUFFIX to RUNDIR, and BL_SUFFIX to CNTL_DIR
export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}${RT_SUFFIX}
export CNTL_DIR=${CNTL_DIR}${BL_SUFFIX}

export JBNME=$(basename $RUNDIR_ROOT)_${TEST_NR}

echo -n "${TEST_NAME}, $( date +%s )," > ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

if [[ ${OPNREQ_TEST} == false ]]; then
  REGRESSIONTEST_LOG=${LOG_DIR}/rt_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log
else
  REGRESSIONTEST_LOG=${LOG_DIR}/opnReqTest_${TEST_NAME}${RT_SUFFIX}.log
fi
export REGRESSIONTEST_LOG

rm -f ${REGRESSIONTEST_LOG}

echo "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR}"

source rt_utils.sh
source atparse.bash

rm -rf ${RUNDIR}
mkdir -p ${RUNDIR}
cd $RUNDIR

#------------
###############################################################################
# Make configure and run files
###############################################################################

# FV3 executable:
cp ${PATHRT}/fv3_${COMPILE_NR}.exe                 fv3.exe

# modulefile for FV3 prerequisites:
cp ${PATHRT}/modules.fv3_${COMPILE_NR}             modules.fv3
cp ${PATHTR}/modulefiles/ufs_common*               .

# Get the shell file that loads the "module" command and purges modules:
cp ${PATHRT}/module-setup.sh                       module-setup.sh

SRCD="${PATHTR}"
RUND="${RUNDIR}"

# FV3_RUN could have multiple entry seperated by space
for i in ${FV3_RUN:-fv3_run.IN}
do
  atparse < ${PATHRT}/fv3_conf/${i} >> fv3_run
done

if [[ $DATM_CDEPS = 'true' ]] || [[ $FV3 = 'true' ]] || [[ $S2S = 'true' ]]; then
  if [[ $HAFS = 'false' ]] || [[ $FV3 = 'true' && $HAFS = 'true' ]]; then
    atparse < ${PATHRT}/parm/${INPUT_NML:-input.nml.IN} > input.nml
  fi
fi

atparse < ${PATHRT}/parm/${MODEL_CONFIGURE:-model_configure.IN} > model_configure

if [[ $DATM_CDEPS = 'false' ]]; then
  if [[ ${ATM_compute_tasks:-0} -eq 0 ]]; then
    ATM_compute_tasks=$((INPES * JNPES * NTILES))
  fi
  if [[ $QUILTING = '.true.' ]]; then
    ATM_io_tasks=$((WRITE_GROUP * WRTTASK_PER_GROUP))
  fi
fi

compute_petbounds

atparse < ${PATHRT}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure

# remove after all tests pass
if [[ $TASKS -ne $UFS_tasks ]]; then
   echo "$TASKS -ne $UFS_tasks "
  exit 1
fi

if [[ "Q${INPUT_NEST02_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST02; JNPES_NEST=$JNPES_NEST02
    NPX_NEST=$NPX_NEST02; NPY_NEST=$NPY_NEST02
    K_SPLIT_NEST=$K_SPLIT_NEST02; N_SPLIT_NEST=$N_SPLIT_NEST02
    atparse < ${PATHRT}/parm/${INPUT_NEST02_NML} > input_nest02.nml
else
    sed -i -e "/<output_grid_02>/,/<\/output_grid_02>/d" model_configure
fi

if [[ "Q${INPUT_NEST03_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST03; JNPES_NEST=$JNPES_NEST03
    NPX_NEST=$NPX_NEST03; NPY_NEST=$NPY_NEST03
    K_SPLIT_NEST=$K_SPLIT_NEST03; N_SPLIT_NEST=$N_SPLIT_NEST03
    atparse < ${PATHRT}/parm/${INPUT_NEST03_NML} > input_nest03.nml
else
    sed -i -e "/<output_grid_03>/,/<\/output_grid_03>/d" model_configure
fi

if [[ "Q${INPUT_NEST04_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST04; JNPES_NEST=$JNPES_NEST04
    NPX_NEST=$NPX_NEST04; NPY_NEST=$NPY_NEST04
    K_SPLIT_NEST=$K_SPLIT_NEST04; N_SPLIT_NEST=$N_SPLIT_NEST04
    atparse < ${PATHRT}/parm/${INPUT_NEST04_NML} > input_nest04.nml
else
    sed -i -e "/<output_grid_04>/,/<\/output_grid_04>/d" model_configure
fi

if [[ "Q${INPUT_NEST05_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST05; JNPES_NEST=$JNPES_NEST05
    NPX_NEST=$NPX_NEST05; NPY_NEST=$NPY_NEST05
    K_SPLIT_NEST=$K_SPLIT_NEST05; N_SPLIT_NEST=$N_SPLIT_NEST05
    atparse < ${PATHRT}/parm/${INPUT_NEST05_NML} > input_nest05.nml
else
    sed -i -e "/<output_grid_05>/,/<\/output_grid_05>/d" model_configure
fi

if [[ "Q${INPUT_NEST06_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST06; JNPES_NEST=$JNPES_NEST06
    NPX_NEST=$NPX_NEST06; NPY_NEST=$NPY_NEST06
    K_SPLIT_NEST=$K_SPLIT_NEST06; N_SPLIT_NEST=$N_SPLIT_NEST06
    atparse < ${PATHRT}/parm/${INPUT_NEST06_NML} > input_nest06.nml
else
    sed -i -e "/<output_grid_06>/,/<\/output_grid_06>/d" model_configure
fi





# Set up the run directory
#source ./fv3_run

#--------------------------jk
