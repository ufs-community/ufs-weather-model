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
  if [[ ${OPNREQ_TEST} == true ]]; then
    echo ${TEST_NR} $TEST_NAME >> $PATHRT/fail_opnreq_test
  else
    echo "${TEST_NAME} ${TEST_NR} failed in run_test" >> $PATHRT/fail_test_${TEST_NR}
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
rm -f fail_test_${TEST_NR}

[[ -e ${RUNDIR_ROOT}/run_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/run_test_${TEST_NR}.env
source default_vars.sh
source tests/$TEST_NAME
[[ -e ${RUNDIR_ROOT}/opnreq_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/opnreq_test_${TEST_NR}.env

# Save original CNTL_DIR name as INPUT_DIR for regression
# tests that try to copy input data from CNTL_DIR
export INPUT_DIR=${CNTL_DIR}
# Append RT_SUFFIX to RUNDIR, and BL_SUFFIX to CNTL_DIR
export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}${RT_SUFFIX}
export CNTL_DIR=${CNTL_DIR}${BL_SUFFIX}

export JBNME=$(basename $RUNDIR_ROOT)_${TEST_NR}

echo -n "${TEST_NAME}, $( date +%s )," > ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

OPNREQ_TEST=${OPNREQ_TEST:-false}
if [[ ${OPNREQ_TEST} == false ]]; then
  REGRESSIONTEST_LOG=${LOG_DIR}/rt_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log
else
  REGRESSIONTEST_LOG=${LOG_DIR}/opnReqTest_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log
fi
export REGRESSIONTEST_LOG

echo "Test ${TEST_NR} ${TEST_NAME} ${TEST_DESCR}"

source rt_utils.sh
source atparse.bash

rm -rf ${RUNDIR}
mkdir -p ${RUNDIR}
cd $RUNDIR

###############################################################################
# Make configure and run files
###############################################################################

# FV3 executable:
cp ${PATHRT}/fv3_${COMPILE_NR}.exe                 fv3.exe

# modulefile for FV3 prerequisites:
cp ${PATHRT}/modules.fv3_${COMPILE_NR}             modules.fv3
cp ${PATHTR}/modulefiles/ufs_common*               .

# Get the shell file that loads the "module" command and purges modules:
cp ${PATHRT}/../NEMS/src/conf/module-setup.sh.inc  module-setup.sh

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

atparse < ${PATHRT}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure

if [[ "Q${INPUT_NEST02_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST02; JNPES_NEST=$JNPES_NEST02
    NPX_NEST=$NPX_NEST02; NPY_NEST=$NPY_NEST02
    K_SPLIT_NEST=$K_SPLIT_NEST02; N_SPLIT_NEST=$N_SPLIT_NEST02
    atparse < ${PATHRT}/parm/${INPUT_NEST02_NML} > input_nest02.nml
fi

if [[ "Q${INPUT_NEST03_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST03; JNPES_NEST=$JNPES_NEST03
    NPX_NEST=$NPX_NEST03; NPY_NEST=$NPY_NEST03
    K_SPLIT_NEST=$K_SPLIT_NEST03; N_SPLIT_NEST=$N_SPLIT_NEST03
    atparse < ${PATHRT}/parm/${INPUT_NEST03_NML} > input_nest03.nml
fi

if [[ "Q${INPUT_NEST04_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST04; JNPES_NEST=$JNPES_NEST04
    NPX_NEST=$NPX_NEST04; NPY_NEST=$NPY_NEST04
    K_SPLIT_NEST=$K_SPLIT_NEST04; N_SPLIT_NEST=$N_SPLIT_NEST04
    atparse < ${PATHRT}/parm/${INPUT_NEST04_NML} > input_nest04.nml
fi

if [[ "Q${INPUT_NEST05_NML:-}" != Q ]] ; then
    INPES_NEST=$INPES_NEST05; JNPES_NEST=$JNPES_NEST05
    NPX_NEST=$NPX_NEST05; NPY_NEST=$NPY_NEST05
    K_SPLIT_NEST=$K_SPLIT_NEST05; N_SPLIT_NEST=$N_SPLIT_NEST05
    atparse < ${PATHRT}/parm/${INPUT_NEST05_NML} > input_nest05.nml
fi

# diag table
if [[ "Q${DIAG_TABLE:-}" != Q ]] ; then
  cp ${PATHRT}/parm/diag_table/${DIAG_TABLE} diag_table
fi
# Field table
if [[ "Q${FIELD_TABLE:-}" != Q ]] ; then
  cp ${PATHRT}/parm/field_table/${FIELD_TABLE} field_table
fi

# fix files
if [[ $FV3 == true ]]; then
  cp ${INPUTDATA_ROOT}/FV3_fix/*.txt .
  cp ${INPUTDATA_ROOT}/FV3_fix/*.f77 .
  cp ${INPUTDATA_ROOT}/FV3_fix/*.dat .
  cp ${INPUTDATA_ROOT}/FV3_fix/fix_co2_proj/* .
  if [[ $TILEDFIX != .true. ]]; then
    cp ${INPUTDATA_ROOT}/FV3_fix/*.grb .
  fi
fi

# Field Dictionary
cp ${PATHRT}/parm/fd_nems.yaml fd_nems.yaml

# Set up the run directory
source ./fv3_run

if [[ $CPLWAV == .true. ]]; then
  atparse < ${PATHRT}/parm/ww3_multi.inp.IN > ww3_multi.inp
fi

if [[ $DATM_CDEPS = 'true' ]] || [[ $S2S = 'true' ]]; then
  if [[ $HAFS = 'false' ]]; then
    atparse < ${PATHRT}/parm/ice_in_template > ice_in
    atparse < ${PATHRT}/parm/${MOM_INPUT:-MOM_input_template_$OCNRES} > INPUT/MOM_input
    atparse < ${PATHRT}/parm/diag_table/${DIAG_TABLE:-diag_table_template} > diag_table
    atparse < ${PATHRT}/parm/data_table_template > data_table
  fi
fi

if [[ $HAFS = 'true' ]] && [[ $DATM_CDEPS = 'false' ]]; then
  atparse < ${PATHRT}/parm/diag_table/${DIAG_TABLE:-diag_table_template} > diag_table
fi

if [[ "${DIAG_TABLE_ADDITIONAL:-}Q" != Q ]] ; then
  # Append diagnostic outputs, to support tests that vary from others
  # only by adding diagnostics.
  atparse < "${PATHRT}/parm/diag_table/${DIAG_TABLE_ADDITIONAL:-}" >> diag_table
fi

if [[ $DATM_CDEPS = 'true' ]]; then
  atparse < ${PATHRT}/parm/${DATM_IN_CONFIGURE:-datm_in} > datm_in
  atparse < ${PATHRT}/parm/${DATM_STREAM_CONFIGURE:-datm.streams.IN} > datm.streams
fi

if [[ $DOCN_CDEPS = 'true' ]]; then
  atparse < ${PATHRT}/parm/${DOCN_IN_CONFIGURE:-docn_in} > docn_in
  atparse < ${PATHRT}/parm/${DOCN_STREAM_CONFIGURE:-docn.streams.IN} > docn.streams
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

################################################################################
# Submit test job
################################################################################

if [[ $SCHEDULER = 'none' ]]; then

  ulimit -s unlimited
  if [[ $CI_TEST = 'true' ]]; then
    eval ${OMP_ENV} mpiexec -n ${TASKS} ${MPI_PROC_BIND} ./fv3.exe >out 2> >(tee err >&3)
  else
    mpiexec -n ${TASKS} ./fv3.exe >out 2> >(tee err >&3)
  fi

else

  if [[ $ROCOTO = 'false' ]]; then
    submit_and_wait job_card
  else
    chmod u+x job_card
    ./job_card
  fi

fi

check_results

if [[ $SCHEDULER != 'none' ]]; then
  cat ${RUNDIR}/job_timestamp.txt >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt
fi
################################################################################
# End test
################################################################################

echo " $( date +%s )" >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Test ${TEST_NAME}"
