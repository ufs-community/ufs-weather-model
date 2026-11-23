#!/bin/bash
set -eu
set -o pipefail
[[ "${RTVERBOSE:-false}" == true ]] && set -x

echo "PID=$$"
SECONDS=0

trap '[ "$?" -eq 0 ] || write_fail_test' EXIT
trap 'echo "run_compile.sh interrupted PID=$$"; cleanup' INT
trap 'echo "run_compile.sh terminated PID=$$";  cleanup' TERM

cleanup() {
  [[ ${ROCOTO} = 'false' ]] && interrupt_job
  trap 0
  exit
}

write_fail_test() {
  echo "${JBNME} failed in run_compile" >> "${PATHRT}/fail_${JBNME}"
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

remove_fail_test() {
    echo "Removing test failure flag file for ${JBNME}"
    rm -f "${PATHRT}/fail_${JBNME}"
}

if [[ $# != 4 ]]; then
  echo "Usage: $0 PATHRT RUNDIR_ROOT MAKE_OPT COMPILE_ID"
  exit 1
fi

export PATHRT=$1
export RUNDIR_ROOT=$2
export MAKE_OPT=$3
export COMPILE_ID=$4

export JBNME="compile_${COMPILE_ID}"

cd "${PATHRT}"
remove_fail_test

# shellcheck disable=SC1090
[[ -e ${RUNDIR_ROOT}/${JBNME}.env ]] && source "${RUNDIR_ROOT}/${JBNME}.env"
source default_vars.sh
# shellcheck disable=SC1090
[[ -e ${RUNDIR_ROOT}/${JBNME}.env ]] && source "${RUNDIR_ROOT}/${JBNME}.env"

export RUNDIR=${RUNDIR_ROOT}/${JBNME}
date_s=$( date +%s )
echo -n "${JBNME}, ${date_s}," > "${LOG_DIR}/${JBNME}_timestamp.txt"

export RT_LOG=${LOG_DIR}/${JBNME}.log

source rt_utils.sh
source atparse.bash

if [[ ${USE_CONTAINER:-false} == true && -d ${RUNDIR} ]]; then
  # A community platform (-P) reuses a fixed RUNDIR_ROOT across runs; keep
  # a prior compile directory instead of clobbering it (rt.sh only invokes
  # this compile at all when the executable wasn't already found reusable).
  mv "${RUNDIR}" "${RUNDIR}_old_$(date +%Y%m%d%H%M)"
else
  rm -rf "${RUNDIR}"
fi
mkdir -p "${RUNDIR}"
cd "${RUNDIR}"

if [[ ${MACHINE_ID} = container ]]; then
    mkdir -p modulefiles
    if [[ -f "${PATHTR}/modulefiles/ufs_container.${RT_COMPILER}.lua" ]]; then
        cp "${PATHTR}/modulefiles/ufs_container.${RT_COMPILER}.lua" "modulefiles/modules.fv3.lua"
    else
        echo "ERROR: no inside-container build module found" >&2
        echo "       Provide ${PATHTR}/modulefiles/ufs_container.${RT_COMPILER}.lua" >&2
        exit 1
    fi
    # Host-side runtime modulefile — loaded on the host before the container launches.
    cp "${PATHTR}/modulefiles/ufs_container.runtime.lua" "modulefiles/modules.fv3.runtime.lua"
    [[ -f "${PATHTR}/modulefiles/ufs_common.lua" ]] && \
        cp "${PATHTR}/modulefiles/ufs_common.lua" "modulefiles/ufs_common.lua"
    cp "${PATHRT}/module-setup.sh" "module-setup.sh"
fi


if [[ ${SCHEDULER} = 'pbs' ]]; then
  if [[ -e ${PATHRT}/fv3_conf/compile_qsub.IN_${MACHINE_ID} ]]; then
    atparse < "${PATHRT}/fv3_conf/compile_qsub.IN_${MACHINE_ID}" > job_card
  else
    echo "Looking for fv3_conf/compile_qsub.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
elif [[ ${SCHEDULER} = 'slurm' ]]; then
  if [[ -e ${PATHRT}/fv3_conf/compile_slurm.IN_${MACHINE_ID} ]]; then
    atparse < "${PATHRT}/fv3_conf/compile_slurm.IN_${MACHINE_ID}" > job_card
  else
    echo "Looking for fv3_conf/compile_slurm.IN_${MACHINE_ID} but it is not found. Exiting"
    exit 1
  fi
elif [[ ${SCHEDULER:-none} = 'none' && ${MACHINE_ID} = 'container' ]]; then
  # No job scheduler: write a compile script with no scheduler headers.
  # container_exec runs it inside the container.
  # Unquoted heredoc: ${...} variables expand at write time (COMPILE_ID etc. are baked in).
  cat > job_card << COMPILE_EOF
#!/bin/bash
set -e
MACHINE_ID=container
source ${PWD}/module-setup.sh
module purge
module use ${PWD}/modulefiles
module load modules.fv3
module list
"${PATHRT}/compile.sh" "${MACHINE_ID}" "${MAKE_OPT}" "${COMPILE_ID}" "${RT_COMPILER}"
COMPILE_EOF
  chmod u+x job_card
fi

################################################################################
# Submit compile job
################################################################################

if [[ ${SCHEDULER:-none} = 'none' && ${MACHINE_ID} = 'container' ]]; then
  # No scheduler: compile on this host, inside the container.
  container_exec job_card
elif [[ ${ROCOTO} = 'false' ]]; then
  submit_and_wait job_card
else
  chmod u+x job_card
  redirect_out_err ./job_card
fi
#ls -l "${PATHTR}/tests/fv3_${COMPILE_ID}.exe"

cp "${RUNDIR}/${JBNME}_time.log" "${LOG_DIR}"
cat "${RUNDIR}/job_timestamp.txt" >> "${LOG_DIR}/${JBNME}_timestamp.txt"

remove_fail_test

################################################################################
# End compile job
################################################################################
date_s=$( date +%s )
echo " ${date_s}, 1" >> "${LOG_DIR}/${JBNME}_timestamp.txt"

elapsed=${SECONDS}
echo "run_compile.sh: Compile ${COMPILE_ID} Completed."
echo "run_compile.sh: Compile ${COMPILE_ID} Elapsed time ${elapsed} seconds."
