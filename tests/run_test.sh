#!/bin/bash
set -eux
set -o pipefail

echo "PID=$$"
SECONDS=0

trap '[ "$?" -eq 0 ] || write_fail_test' EXIT
trap 'echo "run_test.sh interrupted PID=$$"; cleanup' INT
trap 'echo "run_test.sh terminated PID=$$";  cleanup' TERM

cleanup() {
  [[ ${ROCOTO} = 'false' ]] && interrupt_job
  trap 0
  exit
}

write_fail_test() {
  echo "${TEST_ID} failed in run_test" >> "${PATHRT}/fail_test_${TEST_ID}"
  if [[ ${ROCOTO:-false} == true ]] || [[ ${ECFLOW:-false} == true ]]; then
    # if this script has been submitted by a workflow return non-zero exit status
    # so that workflow can resubmit it
    exit 1
  else
    # if this script has been executed interactively, return zero exit status
    # so that rt.sh can continue running, and hope that rt.sh's generate_log
    # will catch failed tests
    exit 0
  fi
}

if [[ $# != 5 ]]; then
  echo "Usage: $0 PATHRT RUNDIR_ROOT TEST_NAME TEST_ID COMPILE_ID"
  exit 1
fi

export PATHRT=$1
export RUNDIR_ROOT=$2
export TEST_NAME=$3
export TEST_ID=$4
export COMPILE_ID=$5

echo "PATHRT: ${PATHRT}"
echo "RUNDIR_ROOT: ${RUNDIR_ROOT}"
echo "TEST_NAME: ${TEST_NAME}"
echo "TEST_ID: ${TEST_ID}"
echo "COMPILE_ID: ${COMPILE_ID}"

cd "${PATHRT}"


unset MODEL_CONFIGURE
unset UFS_CONFIGURE

[[ -e ${RUNDIR_ROOT}/run_test_${TEST_ID}.env ]] && source "${RUNDIR_ROOT}/run_test_${TEST_ID}.env"
source default_vars.sh
[[ -e ${RUNDIR_ROOT}/run_test_${TEST_ID}.env ]] && source "${RUNDIR_ROOT}/run_test_${TEST_ID}.env"
source "tests/${TEST_NAME}"

rm -f "${PATHRT}/fail_test_${TEST_ID}"

# Save original CNTL_DIR name as INPUT_DIR for regression
# tests that try to copy input data from CNTL_DIR

export INPUT_DIR=${CNTL_DIR}

# Append RT_SUFFIX to RUNDIR, and BL_SUFFIX to CNTL_DIR
export RUNDIR=${RUNDIR_ROOT}/${TEST_ID}${RT_SUFFIX}
export CNTL_DIR=${CNTL_DIR}${BL_SUFFIX}

JBNME="run_${TEST_ID}"
export JBNME
date_s=$( date +%s )
echo -n "${TEST_ID}, ${date_s}," > "${LOG_DIR}/${JBNME}_timestamp.txt"

export RT_LOG=${LOG_DIR}/rt_${TEST_ID}${RT_SUFFIX}.log
echo "Test ${TEST_ID} ${TEST_DESCR}"

source rt_utils.sh
source atparse.bash

rm -rf "${RUNDIR}"
mkdir -p "${RUNDIR}"
cd "${RUNDIR}"

###############################################################################
# Make configure and run files
###############################################################################

# FV3 executable:
cp "${PATHRT}/fv3_${COMPILE_ID}.exe" "fv3.exe"

# modulefile for FV3 prerequisites:
mkdir -p modulefiles
if [[ ${MACHINE_ID} == linux ]]; then
  cp "${PATHRT}/modules.fv3_${COMPILE_ID}" "./modulefiles/modules.fv3"
else
  cp "${PATHRT}/modules.fv3_${COMPILE_ID}.lua" "./modulefiles/modules.fv3.lua"
fi
cp "${PATHTR}/modulefiles/ufs_common.lua" "./modulefiles/."

# Get the shell file that loads the "module" command and purges modules:
cp "${PATHRT}/module-setup.sh" "module-setup.sh"

case ${MACHINE_ID} in
  wcoss2|acorn)
    module load intel/19.1.3.304 netcdf/4.7.4
    module load nccmp
    ;;
  s4)
    module use /data/prod/jedi/spack-stack/spack-stack-1.4.1/envs/ufs-pio-2.5.10/install/modulefiles/Core
    module load stack-intel/2021.5.0 stack-intel-oneapi-mpi/2021.5.0
    module load miniconda/3.9.12
    module load nccmp/1.9.0.1
    ;;
  stampede|expanse|noaacloud)
    echo "No special nccmp load necessary"
    ;;
  gaea)
    module use /ncrc/proj/epic/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core
    module load stack-intel/2023.1.0 stack-cray-mpich/8.1.25
    module load nccmp/1.9.0.1
    ;;
  derecho)
    module load nccmp
    ;;
  *)
    module use modulefiles
    module load modules.fv3
    ;;
esac

# FV3_RUN could have multiple entry seperated by space
if [[ -n "${FV3_RUN}" ]]; then
  for i in ${FV3_RUN}
  do
    atparse < "${PATHRT}/fv3_conf/${i}" >> fv3_run
  done
else
  echo "No FV3_RUN set in test file"
  exit 1
fi

# Magic to handle namelist versions of &cires_ugwp_nml
if [[ ${DO_UGWP_V1:-.false.} == .true. ]] ; then
  export HIDE_UGWPV0='!'
  export HIDE_UGWPV1=' '
else
  export HIDE_UGWPV0=' '
  export HIDE_UGWPV1='!'
fi

if [[ ${DATM_CDEPS} = 'true' ]] || [[ ${FV3} = 'true' ]] || [[ ${S2S} = 'true' ]]; then
  if [[ ${HAFS} = 'false' ]] || [[ ${FV3} = 'true' && ${HAFS} = 'true' ]]; then
    atparse < "${PATHRT}/parm/${INPUT_NML:-input.nml.IN}" > input.nml
  fi
fi

if [[ -f ${PATHRT}/parm/${MODEL_CONFIGURE} ]]; then
  atparse < "${PATHRT}/parm/${MODEL_CONFIGURE}" > model_configure
else
  echo "Cannot find file ${MODEL_CONFIGURE} set by variable MODEL_CONFIGURE"
  exit 1
fi

compute_petbounds_and_tasks

if [[ -f ${PATHRT}/parm/${UFS_CONFIGURE} ]]; then
  atparse < "${PATHRT}/parm/${UFS_CONFIGURE}" > ufs.configure
else
  echo "Cannot find file ${UFS_CONFIGURE} set by variable UFS_CONFIGURE"
  exit 1
fi

if [[ "Q${INPUT_NEST02_NML:-}" != Q ]]; then
    export INPES_NEST=${INPES_NEST02:-}
    export JNPES_NEST=${JNPES_NEST02:-}
    export NPX_NEST=${NPX_NEST02:-}
    export NPY_NEST=${NPY_NEST02:-}
    export K_SPLIT_NEST=${K_SPLIT_NEST02:-}
    export N_SPLIT_NEST=${N_SPLIT_NEST02:-}
    atparse < "${PATHRT}/parm/${INPUT_NEST02_NML}" > input_nest02.nml
else
    sed -i -e "/<output_grid_02>/,/<\/output_grid_02>/d" model_configure
fi

if [[ "Q${INPUT_NEST03_NML:-}" != Q ]]; then
    export INPES_NEST=${INPES_NEST03:-}
    export JNPES_NEST=${JNPES_NEST03:-}
    export NPX_NEST=${NPX_NEST03:-}
    export NPY_NEST=${NPY_NEST03:-}
    export K_SPLIT_NEST=${K_SPLIT_NEST03:-}
    export N_SPLIT_NEST=${N_SPLIT_NEST03:-}
    atparse < "${PATHRT}/parm/${INPUT_NEST03_NML}" > input_nest03.nml
else
    sed -i -e "/<output_grid_03>/,/<\/output_grid_03>/d" model_configure
fi

if [[ "Q${INPUT_NEST04_NML:-}" != Q ]]; then
    export INPES_NEST=${INPES_NEST04:-}
    export JNPES_NEST=${JNPES_NEST04:-}
    export NPX_NEST=${NPX_NEST04:-}
    export NPY_NEST=${NPY_NEST04:-}
    export K_SPLIT_NEST=${K_SPLIT_NEST04:-}
    export N_SPLIT_NEST=${N_SPLIT_NEST04:-}
    atparse < "${PATHRT}/parm/${INPUT_NEST04_NML}" > input_nest04.nml
else
    sed -i -e "/<output_grid_04>/,/<\/output_grid_04>/d" model_configure
fi

if [[ "Q${INPUT_NEST05_NML:-}" != Q ]]; then
    export INPES_NEST=${INPES_NEST05:-}
    export JNPES_NEST=${JNPES_NEST05:-}
    export NPX_NEST=${NPX_NEST05:-}
    export NPY_NEST=${NPY_NEST05:-}
    export K_SPLIT_NEST=${K_SPLIT_NEST05:-}
    export N_SPLIT_NEST=${N_SPLIT_NEST05:-}
    atparse < "${PATHRT}/parm/${INPUT_NEST05_NML}" > input_nest05.nml
else
    sed -i -e "/<output_grid_05>/,/<\/output_grid_05>/d" model_configure
fi

if [[ "Q${INPUT_NEST06_NML:-}" != Q ]]; then
    export INPES_NEST=${INPES_NEST06:-}
    export JNPES_NEST=${JNPES_NEST06:-}
    export NPX_NEST=${NPX_NEST06:-}
    export NPY_NEST=${NPY_NEST06:-}
    export K_SPLIT_NEST=${K_SPLIT_NEST06:-}
    export N_SPLIT_NEST=${N_SPLIT_NEST06:-}
    atparse < "${PATHRT}/parm/${INPUT_NEST06_NML}" > input_nest06.nml
else
    sed -i -e "/<output_grid_06>/,/<\/output_grid_06>/d" model_configure
fi

# diag table
if [[ "Q${DIAG_TABLE:-}" != Q ]]; then
  atparse < "${PATHRT}/parm/diag_table/${DIAG_TABLE}" > diag_table
fi
# Field table
if [[ "Q${FIELD_TABLE:-}" != Q ]]; then
  cp "${PATHRT}/parm/field_table/${FIELD_TABLE}" field_table
fi
    
# fix files
if [[ ${FV3} == true ]]; then
  cp "${INPUTDATA_ROOT}"/FV3_fix/*.txt .
  cp "${INPUTDATA_ROOT}"/FV3_fix/*.f77 .
  cp "${INPUTDATA_ROOT}"/FV3_fix/*.dat .
  cp "${INPUTDATA_ROOT}"/FV3_fix/fix_co2_proj/* .
  if [[ ${TILEDFIX} != .true. ]]; then
    cp "${INPUTDATA_ROOT}"/FV3_fix/*.grb .
  fi
fi

# NoahMP table file
  cp "${PATHRT}/parm/noahmptable.tbl" .


# AQM
if [[ ${AQM} == .true. ]]; then
  cp "${PATHRT}/parm/aqm/aqm.rc" .
fi

# Field Dictionary
cp "${PATHRT}/parm/fd_ufs.yaml" fd_ufs.yaml

# Set up the run directory
source ./fv3_run

if [[ ${CPLWAV} == .true. ]]; then
  if [[ ${WW3_MULTIGRID} = 'true' ]]; then
    atparse < "${PATHRT}/parm/ww3_multi.inp.IN" > ww3_multi.inp
  else
    atparse < "${PATHRT}/parm/ww3_shel.nml.IN" > ww3_shel.nml
    cp "${PATHRT}/parm/ww3_points.list" .
  fi
fi

if [[ ${CPLCHM} == .true. ]]; then
  cp "${PATHRT}"/parm/gocart/*.rc .
  atparse < "${PATHRT}/parm/gocart/AERO_HISTORY.rc.IN" > AERO_HISTORY.rc
fi

#TODO: this logic needs to be cleaned up for datm applications w/o
#ocean or ice
if [[ ${DATM_CDEPS} = 'true' ]] || [[ ${S2S} = 'true' ]]; then
  if [[ ${HAFS} = 'false' ]]; then
    atparse < "${PATHRT}/parm/ice_in.IN" > ice_in
    atparse < "${PATHRT}/parm/${MOM6_INPUT:-MOM_input_${OCNRES}.IN}" > INPUT/MOM_input
    atparse < "${PATHRT}/parm/diag_table/${DIAG_TABLE:-diag_table_template}" > diag_table
    atparse < "${PATHRT}/parm/MOM6_data_table.IN" > data_table
  fi
fi

if [[ ${HAFS} = 'true' ]] && [[ ${DATM_CDEPS} = 'false' ]]; then
  atparse < "${PATHRT}/parm/diag_table/${DIAG_TABLE:-diag_table_template}" > diag_table
fi

if [[ "${DIAG_TABLE_ADDITIONAL:-}Q" != Q ]]; then
  # Append diagnostic outputs, to support tests that vary from others
  # only by adding diagnostics.
  atparse < "${PATHRT}/parm/diag_table/${DIAG_TABLE_ADDITIONAL:-}" >> diag_table
fi

if [[ "${FIELD_TABLE_ADDITIONAL:-}Q" != Q ]] ; then
    # Append field table
    atparse < "${PATHRT}/parm/field_table/${FIELD_TABLE_ADDITIONAL:-}" >> field_table
fi

# ATMAERO
if [[ ${CPLCHM} == .true. ]] && [[ ${S2S} = 'false' ]]; then
  atparse < "${PATHRT}/parm/diag_table/${DIAG_TABLE:-diag_table_template}" > diag_table
fi

if [[ ${DATM_CDEPS} = 'true' ]]; then
  atparse < "${PATHRT}/parm/${DATM_IN_CONFIGURE:-datm_in.IN}" > datm_in
  atparse < "${PATHRT}/parm/${DATM_STREAM_CONFIGURE:-datm.streams.IN}" > datm.streams
fi

if [[ ${DOCN_CDEPS} = 'true' ]]; then
  atparse < "${PATHRT}/parm/${DOCN_IN_CONFIGURE:-docn_in.IN}" > docn_in
  atparse < "${PATHRT}/parm/${DOCN_STREAM_CONFIGURE:-docn.streams.IN}" > docn.streams
fi

if [[ ${DICE_CDEPS} = 'true' ]]; then
  atparse < "${PATHRT}/parm/${DICE_IN_CONFIGURE:-dice_in.IN}" > dice_in
  atparse < "${PATHRT}/parm/${DICE_STREAM_CONFIGURE:-dice.streams.IN}" > dice.streams
fi

if [[ ${CICE_PRESCRIBED} = 'true' ]]; then
  atparse < "${PATHRT}"/parm/ice_in.IN > ice_in
fi

if [[ ${CDEPS_INLINE} = 'true' ]]; then
  atparse < "${PATHRT}/parm/${CDEPS_INLINE_CONFIGURE:-stream.config.IN}" > stream.config
fi

if [[ ${FIRE_BEHAVIOR} = 'true' ]]; then
  atparse < "${PATHRT}/parm/${FIRE_NML:-namelist.fire.IN}" > namelist.fire
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

UFS_TASKS=${TASKS}
TASKS=$(( NODES * TPN ))
export TASKS

PPN=$(( UFS_TASKS / NODES ))
if (( UFS_TASKS - ( PPN * NODES ) > 0 )); then
  PPN=$((PPN + 1))
fi
export PPN
export UFS_TASKS

if [[ ${SCHEDULER} = 'pbs' ]]; then
  if [[ -e ${PATHRT}/fv3_conf/fv3_qsub.IN_${MACHINE_ID} ]]; then
    atparse < "${PATHRT}/fv3_conf/fv3_qsub.IN_${MACHINE_ID}" > job_card
  else
    echo "Looking for fv3_conf/fv3_qsub.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
elif [[ ${SCHEDULER} = 'slurm' ]]; then
  if [[ -e ${PATHRT}/fv3_conf/fv3_slurm.IN_${MACHINE_ID} ]]; then
    atparse < "${PATHRT}/fv3_conf/fv3_slurm.IN_${MACHINE_ID}" > job_card
  else
    echo "Looking for fv3_conf/fv3_slurm.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
fi

# This "if" block is part of the rt.sh self-tests in error-test.conf.
# It emulates run_test.sh not being able to populate the work directory.
if [[ "${JOB_SHOULD_FAIL:-NO}" == WHEN_COPYING ]] ; then
    echo "The job should abort now, with exit status 1." 1>&2
    echo "If error checking is working, the metascheduler should mark the job as failed." 1>&2
    false
fi

################################################################################
# Submit test job
################################################################################
export OMP_ENV=${OMP_ENV:-""}
if [[ ${SCHEDULER} = 'none' ]]; then
  ulimit -s unlimited
  if [[ ${CI_TEST} = 'true' ]]; then
    eval "${OMP_ENV}" redirect_out_err mpiexec -n "${TASKS}" ./fv3.exe
  else
    redirect_out_err mpiexec -n "${TASKS}" ./fv3.exe
  fi

else

  if [[ ${ROCOTO} = 'false' ]]; then
    submit_and_wait job_card
  else
    chmod u+x job_card
    redirect_out_err ./job_card
  fi

fi
skip_check_results=${skip_check_results:-false}
if [[ ${skip_check_results} == false ]]; then

  test_status='PASS'

  {
  echo
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}_${RT_COMPILER}"
  echo "working dir  = ${RUNDIR}"
  echo "Checking test ${TEST_ID} results ...."
  } > "${RT_LOG}"
  echo
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}_${RT_COMPILER}"
  echo "working dir  = ${RUNDIR}"
  echo "Checking test ${TEST_ID} results ...."

  if [[ ${CREATE_BASELINE} = false ]]; then
    #
    # --- regression test comparison
    #
    for i in ${LIST_FILES} ; do
      printf %s " Comparing ${i} ....." >> "${RT_LOG}"
      printf %s " Comparing ${i} ....."

      if [[ ! -f ${RUNDIR}/${i} ]] ; then

        echo ".......MISSING file" >> "${RT_LOG}"
        echo ".......MISSING file"
        test_status='FAIL'

      elif [[ ! -f ${RTPWD}/${CNTL_DIR}_${RT_COMPILER}/${i} ]] ; then

        echo ".......MISSING baseline" >> "${RT_LOG}"
        echo ".......MISSING baseline"
        test_status='FAIL'

      else
        if [[ ${i##*.} == nc* ]] ; then
          if [[ " orion hercules hera wcoss2 acorn derecho gaea jet s4 noaacloud " =~ ${MACHINE_ID} ]]; then
            printf "USING NCCMP.." >> "${RT_LOG}"
            printf "USING NCCMP.."
              if [[ ${CMP_DATAONLY} == false ]]; then
                nccmp -d -S -q -f -g -B --Attribute=checksum --warn=format "${RTPWD}/${CNTL_DIR}_${RT_COMPILER}/${i}" "${RUNDIR}/${i}" > "${i}_nccmp.log" 2>&1 && d=$? || d=$?
              else
                nccmp -d -S -q -f -B --Attribute=checksum --warn=format "${RTPWD}/${CNTL_DIR}_${RT_COMPILER}/${i}" "${RUNDIR}/${i}" > "${i}_nccmp.log" 2>&1 && d=$? || d=$?
              fi
              if [[ ${d} -ne 0 && ${d} -ne 1 ]]; then
                printf "....ERROR" >> "${RT_LOG}"
                printf "....ERROR"
                test_status='FAIL'
              fi
          fi
        else
          printf "USING CMP.." >> "${RT_LOG}"
          printf "USING CMP.."
          cmp "${RTPWD}/${CNTL_DIR}_${RT_COMPILER}/${i}" "${RUNDIR}/${i}" >/dev/null 2>&1 && d=$? || d=$?
          if [[ ${d} -eq 2 ]]; then
            printf "....ERROR" >> "${RT_LOG}"
            printf "....ERROR"
            test_status='FAIL'
          fi

        fi

        if [[ ${d} -ne 0 ]]; then
          echo "....NOT IDENTICAL" >> "${RT_LOG}"
          echo "....NOT IDENTICAL"
          test_status='FAIL'
        else
          echo "....OK" >> "${RT_LOG}"
          echo "....OK"
        fi

      fi

    done

  else
    #
    # --- create baselines
    #
    echo;echo "Moving baseline ${TEST_ID} files ...."
    echo;echo "Moving baseline ${TEST_ID} files ...." >> "${RT_LOG}"

    for i in ${LIST_FILES} ; do
      printf %s " Moving ${i} ....."
      printf %s " Moving ${i} ....."   >> "${RT_LOG}"
      if [[ -f ${RUNDIR}/${i} ]] ; then
        mkdir -p "${NEW_BASELINE}/${CNTL_DIR}_${RT_COMPILER}/$(dirname "${i}")"
        cp "${RUNDIR}/${i}" "${NEW_BASELINE}/${CNTL_DIR}_${RT_COMPILER}/${i}"
        echo "....OK" >> "${RT_LOG}"
        echo "....OK"
      else
        echo "....NOT OK. Missing ${RUNDIR}/${i}" >> "${RT_LOG}"
        echo "....NOT OK. Missing ${RUNDIR}/${i}"
        test_status='FAIL'
      fi
    done

  fi

  {
  echo
  grep "The total amount of wall time" "${RUNDIR}/out"
  grep "The maximum resident set size" "${RUNDIR}/out"
  echo
  echo "Test ${TEST_ID} ${test_status}"
  echo
  } >> "${RT_LOG}"

  echo "Test ${TEST_ID} ${test_status}"
  echo

  if [[ ${test_status} = 'FAIL' ]]; then
    echo "${TEST_ID} failed in check_result" >> "${PATHRT}/fail_test_${TEST_ID}"
    write_fail_test
  fi

else
  {
  echo
  grep "The total amount of wall time" "${RUNDIR}/out"
  grep "The maximum resident set size" "${RUNDIR}/out"
  echo
  echo "Test ${TEST_ID} RUN_SUCCESS"
  echo;echo;echo
  } >> "${RT_LOG}"
fi

if [[ ${SCHEDULER} != 'none' ]]; then
  cat "${RUNDIR}/job_timestamp.txt" >> "${LOG_DIR}/${JBNME}_timestamp.txt"
fi

################################################################################
# End test
################################################################################

date_s=$( date +%s )
echo " ${date_s}, ${NODES}" >> "${LOG_DIR}/${JBNME}_timestamp.txt"

################################################################################
# Remove RUN_DIRs if they are no longer needed by other tests
################################################################################
delete_rundir=${delete_rundir:-false}
if [[ ${delete_rundir} = true ]]; then
  keep_run_dir=false
  while  read -r line; do
    keep_test=$(echo "${line}" | sed -e 's/^ *//' -e 's/ *$//')
    if [[ ${TEST_NAME} == "${keep_test}" ]]; then
      keep_run_dir=true
    fi
  done < "${PATHRT}/keep_tests.tmp"

  if [[ ${keep_run_dir} == false ]]; then
    rm -rf "${RUNDIR}"
  fi
fi

elapsed=${SECONDS}
echo "run_test.sh: Test ${TEST_ID} Completed."
echo "run_test.sh: Test ${TEST_ID} Elapsed time ${elapsed} seconds."
