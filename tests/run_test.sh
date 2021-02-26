#!/bin/bash
set -eux

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
  if [[ ${UNIT_TEST} == true ]]; then
    echo ${TEST_NR} $TEST_NAME >> $PATHRT/fail_unit_test
  else
    echo "${TEST_NAME} ${TEST_NR} failed in run_test" >> $PATHRT/fail_test
  fi
  exit 1
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

export JBNME=$(basename $RUNDIR_ROOT)_${TEST_NR}

UNIT_TEST=${UNIT_TEST:-false}
if [[ ${UNIT_TEST} == false ]]; then
  REGRESSIONTEST_LOG=${LOG_DIR}/rt_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log
else
  REGRESSIONTEST_LOG=${LOG_DIR}/ut_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log
fi
export REGRESSIONTEST_LOG

echo "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR}"

source rt_utils.sh
source atparse.bash
source edit_inputs.sh

mkdir -p ${RUNDIR}
cd $RUNDIR

###############################################################################
# Make configure and run files
###############################################################################

# FV3 executable:
cp ${PATHRT}/fv3_${COMPILE_NR}.exe                 fv3.exe

# modulefile for FV3 prerequisites:
cp ${PATHRT}/modules.fv3_${COMPILE_NR}             modules.fv3

# Get the shell file that loads the "module" command and purges modules:
cp ${PATHRT}/../NEMS/src/conf/module-setup.sh.inc  module-setup.sh
if [[ $FV3 = 'true' ]]; then
  cp ${PATHRT}/parm/post_itag itag
  cp ${PATHRT}/parm/postxconfig-NT.txt postxconfig-NT.txt
  cp ${PATHRT}/parm/postxconfig-NT_FH00.txt postxconfig-NT_FH00.txt
  cp ${PATHRT}/parm/params_grib2_tbl_new params_grib2_tbl_new
fi

SRCD="${PATHTR}"
RUND="${RUNDIR}"

# FV3_RUN could have multiple entry seperated by space
for i in ${FV3_RUN:-fv3_run.IN}
do
  atparse < ${PATHRT}/fv3_conf/${i} >> fv3_run
done

if [[ $FV3 = 'true' ]] || [[ $DATM = 'true' ]]; then
atparse < ${PATHRT}/parm/${INPUT_NML:-input.nml.IN} > input.nml
fi

atparse < ${PATHRT}/parm/${MODEL_CONFIGURE:-model_configure.IN} > model_configure

atparse < ${PATHRT}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure

if [[ "Q${INPUT_NEST02_NML:-}" != Q ]] ; then
    atparse < ${PATHRT}/parm/${INPUT_NEST02_NML} > input_nest02.nml
fi

if [[ $CDEPS_DATM = 'true' ]]; then
  atparse < ${PATHRT}/parm/${DATM_CONFIGURE_A:-datm_in} > datm_in
  atparse < ${PATHRT}/parm/${DATM_CONFIGURE_B:-datm.streams.xml} > datm.streams.xml
fi

if [[ $CDEPS_DOCN = 'true' ]]; then
  atparse < ${PATHRT}/parm/${DOCN_CONFIGURE_A:-docn_in} > docn_in
  atparse < ${PATHRT}/parm/${DOCN_CONFIGURE_B:-docn.streams.xml} > docn.streams.xml
fi

# Set up the run directory
source ./fv3_run

# CMEPS
if [[ $DATM = 'true' ]] || [[ $CDEPS_DATM = 'true' ]] || [[ $CDEPS_DOCN = 'true' ]] || [[ $S2S = 'true' ]]; then
  cp ${PATHRT}/parm/fd_nems.yaml fd_nems.yaml
fi

if [[ $DATM = 'true' ]] || [[ $S2S = 'true' ]]; then
  edit_ice_in     < ${PATHRT}/parm/ice_in_template > ice_in
  edit_mom_input  < ${PATHRT}/parm/${MOM_INPUT:-MOM_input_template_$OCNRES} > INPUT/MOM_input
  edit_diag_table < ${PATHRT}/parm/diag_table_template > diag_table
  edit_data_table < ${PATHRT}/parm/data_table_template > data_table
fi

if [[ $DATM = 'true' ]]; then
  cp ${PATHRT}/parm/datm_data_table.IN datm_data_table
fi

if [[ $SCHEDULER = 'pbs' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_qsub.IN > job_card
elif [[ $SCHEDULER = 'slurm' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_slurm.IN > job_card
elif [[ $SCHEDULER = 'lsf' ]]; then
  if (( TASKS < TPN )); then
    TPN=${TASKS}
  fi
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_bsub.IN > job_card
fi

atparse < ${PATHRT}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure

################################################################################
# Submit test job
################################################################################

if [[ $SCHEDULER = 'none' ]]; then

  ulimit -s unlimited
  mpiexec -n ${TASKS} ./fv3.exe >out 2> >(tee err >&3)

else

  if [[ $ROCOTO = 'false' ]]; then
    submit_and_wait job_card
  else
    chmod u+x job_card
    ./job_card
  fi

fi

check_results

################################################################################
# End test
################################################################################

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Test ${TEST_NAME}"
