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

# TASKS is now set to UFS_TASKS
export TASKS=$UFS_tasks

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

# diag table
if [[ "Q${DIAG_TABLE:-}" != Q ]] ; then
  atparse < ${PATHRT}/parm/diag_table/${DIAG_TABLE} > diag_table
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

# AQM
if [[ $AQM == .true. ]]; then
  cp ${PATHRT}/parm/aqm/aqm.rc .
fi

# Field Dictionary
cp ${PATHRT}/parm/fd_nems.yaml fd_nems.yaml

# Set up the run directory
source ./fv3_run

if [[ $CPLWAV == .true. ]]; then
  if [[ $MULTIGRID = 'true' ]]; then
    atparse < ${PATHRT}/parm/ww3_multi.inp.IN > ww3_multi.inp
  else
    atparse < ${PATHRT}/parm/ww3_shel.inp.IN > ww3_shel.inp
  fi
fi

if [[ $CPLCHM == .true. ]]; then
  cp ${PATHRT}/parm/gocart/*.rc .
  atparse < ${PATHRT}/parm/gocart/AERO_HISTORY.rc.IN > AERO_HISTORY.rc
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

# ATMAERO
if [[ $CPLCHM == .true. ]] && [[ $S2S = 'false' ]]; then
  atparse < ${PATHRT}/parm/diag_table/${DIAG_TABLE:-diag_table_template} > diag_table
fi

if [[ $DATM_CDEPS = 'true' ]]; then
  atparse < ${PATHRT}/parm/${DATM_IN_CONFIGURE:-datm_in} > datm_in
  atparse < ${PATHRT}/parm/${DATM_STREAM_CONFIGURE:-datm.streams.IN} > datm.streams
fi

if [[ $DOCN_CDEPS = 'true' ]]; then
  atparse < ${PATHRT}/parm/${DOCN_IN_CONFIGURE:-docn_in} > docn_in
  atparse < ${PATHRT}/parm/${DOCN_STREAM_CONFIGURE:-docn.streams.IN} > docn.streams
fi

TPN=$(( TPN / THRD ))
if (( TASKS < TPN )); then
  TPN=${TASKS}
fi
NODES=$(( TASKS / TPN ))
if (( NODES * TPN < TASKS )); then
  NODES=$(( NODES + 1 ))
fi
TASKS=$(( NODES * TPN ))

if [[ $SCHEDULER = 'pbs' ]]; then
  atparse < $PATHRT/fv3_conf/fv3_qsub.IN > job_card
elif [[ $SCHEDULER = 'slurm' ]]; then
  atparse < $PATHRT/fv3_conf/fv3_slurm.IN > job_card
elif [[ $SCHEDULER = 'lsf' ]]; then
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
    ( ./job_card 2>&1 1>&3 3>&- | tee err ) 3>&1 1>&2 | tee out
    # The above shell redirection copies stdout to "out" and stderr to "err"
    # while still sending them to stdout and stderr. It does this without
    # relying on bash-specific extensions or non-standard OS features.
  fi

fi

if [[ $skip_check_results = false ]]; then
  check_results
else
  echo                                               >> ${REGRESSIONTEST_LOG}
  grep "The total amount of wall time" ${RUNDIR}/out >> ${REGRESSIONTEST_LOG}
  grep "The maximum resident set size" ${RUNDIR}/out >> ${REGRESSIONTEST_LOG}
  echo                                               >> ${REGRESSIONTEST_LOG}
  echo "Test ${TEST_NR} ${TEST_NAME} RUN_SUCCESS"    >> ${REGRESSIONTEST_LOG}
  echo;echo;echo                                     >> ${REGRESSIONTEST_LOG}
fi

if [[ $SCHEDULER != 'none' ]]; then
  cat ${RUNDIR}/job_timestamp.txt >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt
fi

remove_fail_test

################################################################################
# End test
################################################################################

echo " $( date +%s ), ${NODES}" >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

################################################################################
# Remove RUN_DIRs if they are no longer needed by other tests
################################################################################
if [[ ${delete_rundir} = true ]]; then
  keep_run_dir=false
  while  read -r line; do
    keep_test=$(echo $line| sed -e 's/^ *//' -e 's/ *$//')
    if [[ $TEST_NAME == ${keep_test} ]]; then
      keep_run_dir=true
    fi
  done < ${PATHRT}/keep_tests.tmp

  if [[ ${keep_run_dir} == false ]]; then
    rm -rf ${RUNDIR}
  fi
fi

elapsed=$SECONDS
echo "Elapsed time $elapsed seconds. Test ${TEST_NAME}"
