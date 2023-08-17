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
  echo "${TEST_NAME}_${RT_COMPILER} ${TEST_NR} failed in run_test" >> $PATHRT/fail_test_${TEST_NR}
  exit 1
}

remove_fail_test() {
    echo "Removing test failure flag file for ${TEST_NAME}_${RT_COMPILER} ${TEST_NR}"
    rm -f $PATHRT/fail_test_${TEST_NR}
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

echo "PATHRT: ${PATHRT}"
echo "RUNDIR_ROOT: ${RUNDIR_ROOT}"
echo "TEST_NAME: ${TEST_NAME}"
echo "TEST_NR: ${TEST_NR}"
echo "COMPILE_NR: ${COMPILE_NR}"

cd ${PATHRT}


unset MODEL_CONFIGURE
unset NEMS_CONFIGURE

[[ -e ${RUNDIR_ROOT}/run_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/run_test_${TEST_NR}.env
source default_vars.sh
[[ -e ${RUNDIR_ROOT}/run_test_${TEST_NR}.env ]] && source ${RUNDIR_ROOT}/run_test_${TEST_NR}.env
source tests/$TEST_NAME

remove_fail_test

# Save original CNTL_DIR name as INPUT_DIR for regression
# tests that try to copy input data from CNTL_DIR

export INPUT_DIR=${CNTL_DIR}

# Append RT_SUFFIX to RUNDIR, and BL_SUFFIX to CNTL_DIR
export RUNDIR=${RUNDIR_ROOT}/${TEST_NAME}_${RT_COMPILER}${RT_SUFFIX}
export CNTL_DIR=${CNTL_DIR}${BL_SUFFIX}

export JBNME=$(basename $RUNDIR_ROOT)_${TEST_NR}

echo -n "${TEST_NAME}_${RT_COMPILER}, $( date +%s )," > ${LOG_DIR}/job_${JOB_NR}_timestamp.txt

export RT_LOG=${LOG_DIR}/rt_${TEST_NR}_${TEST_NAME}_${RT_COMPILER}${RT_SUFFIX}.log
echo "Test ${TEST_NR} ${TEST_NAME}_${RT_COMPILER} ${TEST_DESCR}"

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
mkdir -p modulefiles
if [[ $MACHINE_ID == linux ]]; then
  cp ${PATHRT}/modules.fv3_${COMPILE_NR}             ./modulefiles/modules.fv3
else
  cp ${PATHRT}/modules.fv3_${COMPILE_NR}.lua             ./modulefiles/modules.fv3.lua
fi
cp ${PATHTR}/modulefiles/ufs_common*               ./modulefiles/.

# Get the shell file that loads the "module" command and purges modules:
cp ${PATHRT}/module-setup.sh                       module-setup.sh

# load nccmp module
if [[ " hera orion gaea jet cheyenne acorn wcoss2 " =~ " $MACHINE_ID " ]]; then
  if [[ " wcoss2 acorn " =~ " ${MACHINE_ID} " ]] ; then
    module load intel/19.1.3.304 netcdf/4.7.4
  fi
  module load nccmp
fi

SRCD="${PATHTR}"
RUND="${RUNDIR}"

# FV3_RUN could have multiple entry seperated by space
if [ ! -z "$FV3_RUN" ]; then
  for i in ${FV3_RUN}
  do
    atparse < ${PATHRT}/fv3_conf/${i} >> fv3_run
  done
else
  echo "No FV3_RUN set in test file"
  exit 1
fi

if [[ $DATM_CDEPS = 'true' ]] || [[ $FV3 = 'true' ]] || [[ $S2S = 'true' ]]; then
  if [[ $HAFS = 'false' ]] || [[ $FV3 = 'true' && $HAFS = 'true' ]]; then
    atparse < ${PATHRT}/parm/${INPUT_NML:-input.nml.IN} > input.nml
  fi
fi

if [[ -f ${PATHRT}/parm/${MODEL_CONFIGURE} ]]; then
  atparse < ${PATHRT}/parm/${MODEL_CONFIGURE} > model_configure
else
  echo "Cannot find file ${MODEL_CONFIGURE} set by variable MODEL_CONFIGURE"
  exit 1
fi

compute_petbounds_and_tasks

if [[ -f ${PATHRT}/parm/${NEMS_CONFIGURE} ]]; then
  atparse < ${PATHRT}/parm/${NEMS_CONFIGURE} > nems.configure
else
  echo "Cannot find file ${NEMS_CONFIGURE} set by variable NEMS_CONFIGURE"
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

# NoahMP table file
  cp ${PATHRT}/parm/noahmptable.tbl .


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
    atparse < ${PATHRT}/parm/ww3_shel.nml.IN > ww3_shel.nml
    cp ${PATHRT}/parm/ww3_points.list .
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
export TPN

NODES=$(( TASKS / TPN ))
if (( NODES * TPN < TASKS )); then
  NODES=$(( NODES + 1 ))
fi
export NODES

TASKS=$(( NODES * TPN ))
export TASKS

if [[ $SCHEDULER = 'pbs' ]]; then
  if [[ -e $PATHRT/fv3_conf/fv3_qsub.IN_${MACHINE_ID} ]]; then 
    atparse < $PATHRT/fv3_conf/fv3_qsub.IN_${MACHINE_ID} > job_card
  else
    echo "Looking for fv3_conf/fv3_qsub.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
elif [[ $SCHEDULER = 'slurm' ]]; then
  if [[ -e $PATHRT/fv3_conf/fv3_slurm.IN_${MACHINE_ID} ]]; then
    atparse < $PATHRT/fv3_conf/fv3_slurm.IN_${MACHINE_ID} > job_card
  else
    echo "Looking for fv3_conf/fv3_slurm.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
elif [[ $SCHEDULER = 'lsf' ]]; then
  if [[ -e $PATHRT/fv3_conf/fv3_bsub.IN_${MACHINE_ID} ]]; then
    atparse < $PATHRT/fv3_conf/fv3_bsub.IN_${MACHINE_ID} > job_card
  else
    echo "Looking for fv3_conf/fv3_bsub.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
fi

################################################################################
# Submit test job
################################################################################
export OMP_ENV=${OMP_ENV:-""}
if [[ $SCHEDULER = 'none' ]]; then

  ulimit -s unlimited
  if [[ $CI_TEST = 'true' ]]; then
    eval ${OMP_ENV} mpiexec -n ${TASKS} ./fv3.exe >out 2> >(tee err >&3)
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
  echo                                               >> ${RT_LOG}
  grep "The total amount of wall time" ${RUNDIR}/out >> ${RT_LOG}
  grep "The maximum resident set size" ${RUNDIR}/out >> ${RT_LOG}
  echo                                               >> ${RT_LOG}
  echo "Test ${TEST_NR} ${TEST_NAME}_${RT_COMPILER} RUN_SUCCESS"    >> ${RT_LOG}
  echo;echo;echo                                     >> ${RT_LOG}
fi

if [[ $SCHEDULER != 'none' ]]; then
  cat ${RUNDIR}/job_timestamp.txt >> ${LOG_DIR}/job_${JOB_NR}_timestamp.txt
fi

if [[ $ROCOTO = true ]]; then
  remove_fail_test
fi

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
echo "Elapsed time $elapsed seconds. Test ${TEST_NAME}_${RT_COMPILER}"
