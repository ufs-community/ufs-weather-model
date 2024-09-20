#!/bin/bash
set -eu

if [[ "$0" = "${BASH_SOURCE[0]}" ]]; then
  echo "$0 must be sourced"
  exit 1
fi

ECFLOW_RUNNING=false

# Note: this file must only contain subroutines, and variables that
# are not dependent on the caller. Most regression test variables
# (such as ACCNR) are not set until after rt.sh sources this file.

jobid=0

redirect_out_err() {
    ( set -e -o pipefail ; ( "$@" 2>&1 1>&3 3>&- | tee err ) 3>&1 1>&2 | tee out )
    # The above shell redirection copies stdout to "out" and stderr to "err"
    # while still sending them to stdout and stderr. It ensures the entire
    # redirect_out_err command will return non-zero if "$@" or tee return non-zero.
}

function compute_petbounds_and_tasks() {

  # each test MUST define ${COMPONENT}_tasks variable for all components it is using
  # and MUST NOT define those that it's not using or set the value to 0.

  # ATM is a special case since it is running on the sum of compute and io tasks.
  # CHM component and mediator are running on ATM compute tasks only.

  if [[ ${DATM_CDEPS} = 'false' ]]; then
    if [[ ${ATM_compute_tasks:-0} -eq 0 ]]; then
      ATM_compute_tasks=$((INPES * JNPES * NTILES))
    fi
    if [[ ${QUILTING} = '.true.' ]]; then
      ATM_io_tasks=$((WRITE_GROUP * WRTTASK_PER_GROUP))
    fi
  fi

  local n=0
  unset atm_petlist_bounds ocn_petlist_bounds ice_petlist_bounds wav_petlist_bounds chm_petlist_bounds med_petlist_bounds aqm_petlist_bounds fbh_petlist_bounds

  # ATM
  ATM_io_tasks=${ATM_io_tasks:-0}
  if [[ $((ATM_compute_tasks + ATM_io_tasks)) -gt 0 ]]; then
     atm_petlist_bounds="${n} $((n + ATM_compute_tasks*atm_omp_num_threads + ATM_io_tasks*atm_omp_num_threads - 1))"
     n=$((n + ATM_compute_tasks*atm_omp_num_threads + ATM_io_tasks*atm_omp_num_threads))
  fi

  # OCN
  if [[ ${OCN_tasks:-0} -gt 0 ]]; then
     OCN_tasks=$((OCN_tasks * ocn_omp_num_threads))
     ocn_petlist_bounds="${n} $((n + OCN_tasks - 1))"
     n=$((n + OCN_tasks))
  fi

  # ICE
  if [[ ${ICE_tasks:-0} -gt 0 ]]; then
     ICE_tasks=$((ICE_tasks * ice_omp_num_threads))
     ice_petlist_bounds="${n} $((n + ICE_tasks - 1))"
     n=$((n + ICE_tasks))
  fi

  # WAV
  if [[ ${WAV_tasks:-0} -gt 0 ]]; then
     WAV_tasks=$((WAV_tasks * wav_omp_num_threads))
     wav_petlist_bounds="${n} $((n + WAV_tasks - 1))"
     n=$((n + WAV_tasks))
  fi

  # CHM
  chm_petlist_bounds="0 $((ATM_compute_tasks * atm_omp_num_threads - 1))"

  # MED
  med_petlist_bounds="0 $((ATM_compute_tasks * atm_omp_num_threads - 1))"

  # AQM
  aqm_petlist_bounds="0 $((ATM_compute_tasks * atm_omp_num_threads - 1))"

  # LND
  if [[ ${LND_tasks:-0} -gt 0 ]]; then
     LND_tasks=$((LND_tasks * lnd_omp_num_threads))
     lnd_petlist_bounds="${n} $((n + LND_tasks - 1))"
     n=$((n + LND_tasks))
  fi

  # FBH
  if [[ ${FBH_tasks:-0} -gt 0 ]]; then
     FBH_tasks=$((FBH_tasks * fbh_omp_num_threads))
     fbh_petlist_bounds="${n} $((n + FBH_tasks - 1))"
     n=$((n + FBH_tasks))
  fi

  UFS_tasks=${n}

  if [[ ${RTVERBOSE} == true ]]; then
    echo "ATM_petlist_bounds: ${atm_petlist_bounds:-}"
    echo "OCN_petlist_bounds: ${ocn_petlist_bounds:-}"
    echo "ICE_petlist_bounds: ${ice_petlist_bounds:-}"
    echo "WAV_petlist_bounds: ${wav_petlist_bounds:-}"
    echo "CHM_petlist_bounds: ${chm_petlist_bounds:-}"
    echo "MED_petlist_bounds: ${med_petlist_bounds:-}"
    echo "AQM_petlist_bounds: ${aqm_petlist_bounds:-}"
    echo "LND_petlist_bounds: ${lnd_petlist_bounds:-}"
    echo "FBH_petlist_bounds: ${fbh_petlist_bounds:-}"
    echo "UFS_tasks         : ${UFS_tasks:-}"
  fi

  # TASKS is now set to UFS_TASKS
  export TASKS=${UFS_tasks}
}

interrupt_job() {
  echo "rt_utils.sh: Job ${jobid} interrupted"
  case ${SCHEDULER} in
    pbs)
      qdel "${jobid}"
      ;;
    slurm)
      scancel "${jobid}"
      ;;
    *)
      echo "Unsupported scheduler, job may have not terminated properly."
      ;;
  esac
}

submit_and_wait() {
  echo "rt_utils.sh: Submitting job on scheduler: ${SCHEDULER}"
  [[ -z $1 ]] && exit 1

  local -r job_card=$1

  case ${SCHEDULER} in
    pbs)
      qsubout=$( qsub "${job_card}" )
      re='^([0-9]+)(\.[a-zA-Z0-9\.-]+)$'
      [[ "${qsubout}" =~ ${re} ]] && jobid=${BASH_REMATCH[1]}
      ;;
    slurm)
      slurmout=$( sbatch "${job_card}" )
      re='Submitted batch job ([0-9]+)'
      [[ "${slurmout}" =~ ${re} ]] && jobid=${BASH_REMATCH[1]}
      ;;
    *)
      echo "Unsupported scheduler: ${SCHEDULER}"
      exit 1
      ;;
  esac

  echo "rt_utils.sh: Submitted Job. ID is ${jobid}."
  sleep 10
  # wait for the job to enter the queue
  local count=0
  local job_running=''
  echo "rt_utils.sh: Job is waiting to enter the queue..."
  until [[ ${job_running} == 'true' ]]
  do
    case ${SCHEDULER} in
      pbs)
        set +e
        job_info=$( qstat "${jobid}" )
        set -e
        ;;
      slurm)
        job_info=$( squeue -u "${USER}" -j "${jobid}" )
        ;;
      *)
        ;;
    esac
    if grep -q "${jobid}" <<< "${job_info}"; then
      job_running=true
      continue
    else
      job_running=false
    fi

    sleep 5
    (( count=count+1 ))
    if [[ ${count} -eq 13 ]]; then echo "No job in queue after one minute, exiting..."; exit 2; fi
  done
  echo "rt_utils.sh Job (${jobid}) is now in the queue."

  # wait for the job to finish and compare results
  local n=1
  until [[ ${job_running} == 'false' ]]
  do
    case ${SCHEDULER} in
      pbs)
        set +e
        job_info=$( qstat "${jobid}" )
        set -e
        if grep -q "${jobid}" <<< "${job_info}"; then
          job_running=true
          # Getting the status letter from scheduler info
          status=$( grep "${jobid}" <<< "${job_info}" )
          status=$( awk '{print $5}' <<< "${status}" )
        else
          job_running=false
          status='COMPLETED'
          set +e
          exit_status=$( qstat "${jobid}" -x -f | grep Exit_status | awk '{print $3}')
          set -e
          if [[ ${exit_status} != 0 ]]; then
            status='FAILED'
          fi
        fi
        ;;
      slurm)
        job_info=$( squeue -u "${USER}" -j "${jobid}" -o '%i %T' )
        if grep -q "${jobid}" <<< "${job_info}"; then
          job_running=true
        else
          job_running=false
          job_info=$( sacct -n -j "${jobid}" --format=JobID,state%20,Jobname%128 | grep "^${jobid}" | grep "${JBNME}" )
        fi
        # Getting the status letter from scheduler info
        status=$( grep "${jobid}" <<< "${job_info}" )
        status=$( awk '{print $2}' <<< "${status}" )
        ;;
      *)
        ;;
    esac

    case ${status} in
      #waiting cases
      #pbs: Q
      #Slurm: (old: PD, new: PENDING)
      Q|PD|PENDING)
        status_label='Job waiting to start'
        ;;
      #running cases
      #pbs: R
      #slurm: (old: R, new: RUNNING)
      R|RUNNING|COMPLETING)
        status_label='Job running'
        ;;
      #held cases
      #pbs only: H
      H)
        status_label='Job being held'
        echo "rt_utils.sh: *** WARNING ***: Job in a HELD state. Might want to stop manually."
        ;;
      #fail/completed cases
      #slurm: F/FAILED TO/TIMEOUT CA/CANCELLED
      F|TO|CA|FAILED|TIMEOUT|CANCELLED)
        echo "rt_utils.sh: !!!!!!!!!!JOB TERMINATED!!!!!!!!!! status=${status}"
        job_running=false #Trip the loop to end with these status flags
        interrupt_job
        exit 1
        ;;
      #completed
      #pbs: C-Complete E-Exiting
      #slurm: CD/COMPLETED
      C|E|CD|COMPLETED)
        status_label='Completed'
        ;;
      *)
        status_label="Unknown"
        echo "rt_utils.sh: *** WARNING ***: Job status unsupported: ${status}"
        echo "rt_utils.sh: *** WARNING ***: Status might be non-terminating, please manually stop if needed"
        ;;
    esac

    echo "${n} min. ${SCHEDULER^} Job ${jobid} Status: ${status_label} (${status})"

    (( n=n+1 ))
    sleep 60 & wait $!
  done
}

kill_job() {
  echo "rt_utils.sh: Killing job: ${jobid} on ${SCHEDULER}..."
  [[ -z $1 ]] && exit 1

  local -r jobid=$1

  if [[ ${SCHEDULER} = 'pbs' ]]; then
    qdel "${jobid}"
  elif [[ ${SCHEDULER} = 'slurm' ]]; then
    scancel "${jobid}"
  fi
}

rocoto_create_compile_task() {
  echo "rt_utils.sh: ${COMPILE_ID}: Creating ROCOTO compile task."
  new_compile=true
  if [[ ${in_metatask} == true ]]; then
    in_metatask=false
    echo "  </metatask>" >> "${ROCOTO_XML}"
  fi

  NATIVE=""
  BUILD_CORES=8
  BUILD_WALLTIME="00:30:00"
  if [[ ${MACHINE_ID} == jet ]]; then
    BUILD_WALLTIME="02:00:00"
  fi
  if [[ ${MACHINE_ID} == hera ]]; then
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == orion ]]; then
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == hercules ]]; then
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == s4 ]]; then
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == gaea ]]; then
    BUILD_WALLTIME="01:00:00"
  fi


  cat << EOF >> "${ROCOTO_XML}"
  <task name="compile_${COMPILE_ID}" maxtries="${ROCOTO_COMPILE_MAXTRIES:-3}">
    <command>bash -c 'set -xe -o pipefail ; &PATHRT;/run_compile.sh &PATHRT; &RUNDIR_ROOT; "${MAKE_OPT}" ${COMPILE_ID} 2>&amp;1 | tee &LOG;/compile_${COMPILE_ID}.log'</command>
    <jobname>compile_${COMPILE_ID}</jobname>
    <account>${ACCNR}</account>
    <queue>${COMPILE_QUEUE}</queue>
EOF

  if [[ "${MACHINE_ID}" == gaea ]] ; then
    cat << EOF >> "${ROCOTO_XML}"
    <native>--clusters=es</native>
    <partition>eslogin_c5</partition>
EOF
  elif [[ -n "${PARTITION}" || ${MACHINE_ID} != hera ]] ; then
    cat << EOF >> "${ROCOTO_XML}"
    <partition>${PARTITION}</partition>
EOF
  fi

  cat << EOF >> "${ROCOTO_XML}"
    <nodes>1:ppn=${BUILD_CORES}</nodes>
    <walltime>${BUILD_WALLTIME}</walltime>
    <join>&RUNDIR_ROOT;/compile_${COMPILE_ID}.log</join>
    ${NATIVE}
  </task>
EOF
}

rocoto_create_run_task() {
  echo "rt_utils.sh: ${TEST_ID}: Creating ROCOTO run task."
  if [[ ${DEP_RUN} != '' ]]; then
    DEP_STRING="<and> <taskdep task=\"compile_${COMPILE_ID}\"/> <taskdep task=\"${DEP_RUN}\"/> </and>"
  else
    DEP_STRING="<taskdep task=\"compile_${COMPILE_ID}\"/>"
  fi

  CORES=$((TASKS*THRD))
  if (( TPN > CORES )); then
    TPN=${CORES}
  fi

  NATIVE=""

  cat << EOF >> "${ROCOTO_XML}"
    <task name="${TEST_ID}${RT_SUFFIX}" maxtries="${ROCOTO_TEST_MAXTRIES:-3}">
      <dependency> ${DEP_STRING} </dependency>
      <command>bash -c 'set -xe -o pipefail ; &PATHRT;/run_test.sh &PATHRT; &RUNDIR_ROOT; ${TEST_NAME} ${TEST_ID} ${COMPILE_ID} 2>&amp;1 | tee &LOG;/run_${TEST_ID}${RT_SUFFIX}.log' </command>
      <jobname>${TEST_ID}${RT_SUFFIX}</jobname>
      <account>${ACCNR}</account>
      ${ROCOTO_NODESIZE:+<nodesize>${ROCOTO_NODESIZE}</nodesize>}
EOF

  if [[ "${MACHINE_ID}" == gaea ]] ; then
    cat << EOF >> "${ROCOTO_XML}"
      <native>--clusters=${PARTITION}</native>
      <native>--partition=batch</native>
EOF

  elif [[ -n "${PARTITION}" || ${MACHINE_ID} != hera ]] ; then
    cat << EOF >> "${ROCOTO_XML}"
      <queue>${QUEUE}</queue>
      <partition>${PARTITION}</partition>
EOF
  fi

  cat << EOF >> "${ROCOTO_XML}"
      <nodes>${NODES}:ppn=${TPN}</nodes>
      <walltime>00:${WLCLK}:00</walltime>
      <join>&RUNDIR_ROOT;/${TEST_ID}${RT_SUFFIX}.log</join>
      ${NATIVE}
    </task>
EOF

}

rocoto_kill() {
  echo "rt_utils.sh: Killing ROCOTO Workflow..."
  job_id_in=$( "${ROCOTOSTAT}" -w "${ROCOTO_XML}" -d "${ROCOTO_DB}" )
  job_id_in=$(grep 197001010000 <<< "${job_id_in}" )
  job_id_in=$(grep -E 'QUEUED|RUNNING' <<< "${job_id_in}" )
  job_id_in=$(awk -F" " '{print $3}' <<< "${job_id_in}" )
   for jobid in ${job_id_in}; do
      kill_job "${jobid}"
   done
}

rocoto_step() {
    echo "rt_utils.sh: Running one iteration of rocotorun and rocotostat..."
    echo "Unknown" > rocoto_workflow.state
    # Run one iteration of rocotorun and rocotostat.
    ${ROCOTORUN} -v 10 -w "${ROCOTO_XML}" -d "${ROCOTO_DB}"
    sleep 1
    #   Is it done?
    state=$( "${ROCOTOSTAT}" -w "${ROCOTO_XML}" -d "${ROCOTO_DB}" -s )
    state=$( grep 197001010000 <<< "${state}" )
    state=$( awk -F" " '{print $2}' <<< "${state}" )
    echo "${state}" > "${ROCOTO_STATE}"
}

rocoto_run() {
  echo "rt_utils.sh: Running ROCOTO workflow"
  # Run the rocoto workflow until it is complete
  local naptime=60
  local step_attempts=0
  local max_step_attempts=100 # infinite loop safeguard; should never reach this
  local start_time=0
  local now_time=0
  local max_time=3600 # seconds to wait for Rocoto to start working again
  local result=0
  state="Active"
  while [[ ${state} != "Done" ]]; do
      # Run one iteration of rocotorun and rocotostat.  Use an
      # exponential backoff algorithm to handle temporary system
      # failures breaking Rocoto.
      start_time=$( env TZ=UTC date +%s )
      for step_attempts in $( seq 1 "${max_step_attempts}" ) ; do
          now_time=$( env TZ=UTC date +%s )

          set +e
          ( rocoto_step )
          result=$?
          state=$( cat "${ROCOTO_STATE}" )
          set -e

          if [[ "${state:-Unknown}" == Done ]] ; then
              echo "Rocoto workflow has completed."
              return 0
          elif [[ ${result} == 0 ]] ; then
              break # rocoto_step succeeded
          elif (( now_time-start_time > max_time || step_attempts >= max_step_attempts )) ; then
              hostnamein=$(hostname)
              echo "Rocoto commands have failed ${step_attempts} times, for $(( (now_time-start_time+30)/60 )) minutes."
              echo "There may be something wrong with the ${hostnamein} node or the batch system."
              echo "I'm giving up. Sorry."
              return 2
          fi
          sleep $(( naptime * 2**((step_attempts-1)%4) * RANDOM/32767 ))
      done
      sleep "${naptime}"
  done
}


ecflow_create_compile_task() {
  echo "rt_utils.sh: ${COMPILE_ID}: Creating ECFLOW compile task"
  export new_compile=true

  cat << EOF > "${ECFLOW_RUN}/${ECFLOW_SUITE}/compile_${COMPILE_ID}.ecf"
%include <head.h>
(
cd "${LOG_DIR}"
ln -sf "compile_${COMPILE_ID}.log.\${ECF_TRYNO}" "compile_${COMPILE_ID}.log"
)
${PATHRT}/run_compile.sh "${PATHRT}" "${RUNDIR_ROOT}" "${MAKE_OPT}" "${COMPILE_ID}" > "${LOG_DIR}/compile_${COMPILE_ID}.log.\${ECF_TRYNO}" 2>&1 &
%include <tail.h>
EOF
  {
  echo "  task compile_${COMPILE_ID}"
  echo "      label build_options '${MAKE_OPT}'"
  echo "      inlimit max_builds"
  } >> "${ECFLOW_RUN}/${ECFLOW_SUITE}.def"
}

ecflow_create_run_task() {
  echo "rt_utils.sh: ${TEST_ID}: Creating ECFLOW run task"
  cat << EOF > "${ECFLOW_RUN}/${ECFLOW_SUITE}/${TEST_ID}${RT_SUFFIX}.ecf"
%include <head.h>
(
cd "${LOG_DIR}"
ln -sf "run_${TEST_ID}${RT_SUFFIX}.log.\${ECF_TRYNO}" "${LOG_DIR}/run_${TEST_ID}${RT_SUFFIX}.log"
)
${PATHRT}/run_test.sh "${PATHRT}" "${RUNDIR_ROOT}" "${TEST_NAME}" "${TEST_ID}" "${COMPILE_ID}" > "${LOG_DIR}/run_${TEST_ID}${RT_SUFFIX}.log.\${ECF_TRYNO}" 2>&1 &
%include <tail.h>
EOF
  {
  echo "    task ${TEST_ID}${RT_SUFFIX}"
  echo "      inlimit max_jobs"
  } >> "${ECFLOW_RUN}/${ECFLOW_SUITE}.def"
  if [[ ${DEP_RUN} != '' ]]; then
    echo "      trigger compile_${COMPILE_ID} == complete and ${DEP_RUN} == complete" >> "${ECFLOW_RUN}/${ECFLOW_SUITE}.def"
  else
    echo "      trigger compile_${COMPILE_ID} == complete" >> "${ECFLOW_RUN}/${ECFLOW_SUITE}.def"
  fi

}

ecflow_run() {
  echo "rt_utils.sh: Starting ECFLOW run"
  # NOTE: ECFLOW IS NOT SAFE TO RUN WITH set -e, PLEASE AVOID
  #ECF_HOST="${ECF_HOST:-${HOSTNAME}}"


  # Make sure ECF_HOST and ECF_PORT are set/ready on systems that have an
  # explicit ecflow node
  if [[ ${MACHINE_ID} == wcoss2 || ${MACHINE_ID} == acorn ]]; then
    if [[ "${HOST::1}" == "a" ]]; then
      ECF_HOST=aecflow01
    elif [[ "${HOST::1}" == "c" ]]; then
      ECF_HOST=cdecflow01
    elif [[ "${HOST::1}" == "d" ]]; then
      ECF_HOST=ddecflow01
    fi
  elif [[ ${MACHINE_ID} == hera || ${MACHINE_ID} == jet ]]; then
    module load ecflow
  fi
  if [[ -z ${ECF_HOST} || -z ${ECF_PORT} ]]; then
    echo "ERROR: ECF_HOST or ECF_PORT are not set, and rt.sh cannot continue with ECFLOW"
    exit 1
  else
    echo "ECF_HOST: ${ECF_HOST}, ECF_PORT: ${ECF_PORT}"
    export ECF_HOST
    export ECF_PORT
  fi

  # Start the ecflow_server
  echo "rt_utils.sh: Checking status of the ecflow_server..."
  set +e
  ecflow_client --ping --host="${ECF_HOST}" --port="${ECF_PORT}"
  not_running=$?
  set -e

  if [[ ${not_running} -eq 1 ]]; then
    echo "rt_utils.sh: ecflow_server is not running on ${ECF_HOST}:${ECF_PORT}"
    echo "rt_utils.sh: attempting to start ecflow_server..."

    save_traps=$(trap)
    trap "" SIGINT  # Ignore INT signal during ecflow startup
    case ${MACHINE_ID} in
      wcoss2|acorn|hera|jet)
        #shellcheck disable=SC2029
        ssh "${ECF_HOST}" "bash -l -c \"module load ecflow && ${ECFLOW_START} -p ${ECF_PORT}\""
        ;;
      *)
        ${ECFLOW_START} -p "${ECF_PORT}" -d "${RUNDIR_ROOT}/ecflow_server"
        ;;
    esac

    ECFLOW_RUNNING=true
    eval "${save_traps}"
    # Try pinging ecflow server now, and erroring out if not there.
    set +e
    ecflow_client --ping --host="${ECF_HOST}" --port="${ECF_PORT}"
    not_running=$?
    set -e

    if [[ ${not_running} -eq 1 ]]; then
      echo "rt_utils.sh: ERROR -- Failure to start ecflow. Exiting..."
      exit 1
    fi
  else
    echo "rt_utils.sh: Confirmed: ecflow_server is running on ${ECF_HOST}:${ECF_PORT}"
    ECFLOW_RUNNING=true
  fi

  echo "rt_utils.sh: Starting ECFLOW tasks..."
  set +e
  ecflow_client --load="${ECFLOW_RUN}/${ECFLOW_SUITE}.def" --host="${ECF_HOST}" --port="${ECF_PORT}"
  ecflow_client --begin="${ECFLOW_SUITE}" --host="${ECF_HOST}" --port="${ECF_PORT}"
  ecflow_client --restart --host="${ECF_HOST}" --port="${ECF_PORT}"
  set -e
  sleep 10

  active_tasks=1
  max_active_tasks=$( ecflow_client --get_state "/${ECFLOW_SUITE}" )
  max_active_tasks=$( grep "task " <<< "${max_active_tasks}" )
  max_active_tasks=$( grep -cP 'state:active|state:submitted|state:queued' <<< "${max_active_tasks}" )
  echo "rt_utils.sh: Total number of tasks processed -- ${max_active_tasks}"
  prev_active_tasks=${active_tasks}
  while [[ "${active_tasks}" -ne 0 ]]
  do
    sleep 10 & wait $!
    set +e
    active_tasks=$( ecflow_client --get_state "/${ECFLOW_SUITE}" )
    active_tasks=$( grep "task " <<< "${active_tasks}" )
    active_tasks=$( grep -cP 'state:active|state:submitted|state:queued' <<< "${active_tasks}" )
    set -e
    if [[ ${active_tasks} -ne ${prev_active_tasks} ]]; then
      echo
      echo -n "ECFLOW Tasks Remaining: ${active_tasks}/${max_active_tasks} "
      prev_active_tasks=${active_tasks}
    else
      echo -n "."
    fi
    "${PATHRT}/abort_dep_tasks.py"
  done
  echo

  sleep 65 # wait one ECF_INTERVAL plus 5 seconds
  echo "rt_utils.sh: ECFLOW tasks completed, cleaning up suite"
  set +e
  ecflow_client --delete=force yes "/${ECFLOW_SUITE}"
  set -e
  sleep 5
}

ecflow_kill() {
  [[ ${ECFLOW_RUNNING:-false} == true ]] || return
  echo "rt_utils.sh: Deleting ECFLOW suite: ${ECFLOW_SUITE}"
  set +e
  ecflow_client --suspend "/${ECFLOW_SUITE}"
  ecflow_client --kill "/${ECFLOW_SUITE}"
  sleep 20
  ecflow_client --delete=force yes "/${ECFLOW_SUITE}"
  set -e
}

ecflow_stop() {
  [[ ${ECFLOW_RUNNING:-false} == true ]] || return
  echo "rt_utils.sh: Checking whether to stop ecflow_server..."
  set +e
  SUITES=$( ecflow_client --get )
  SUITES=$( grep "^suite" <<< "${SUITES}" )
  if [[ -z "${SUITES}" ]]; then
    echo "rt_utils.sh: No other suites running, stopping ecflow_server"
    ecflow_client --halt=yes
    ecflow_client --check_pt
    ecflow_client --terminate=yes
  else
    echo "rt_utils.sh: Active suites running, NOT stopping ecflow_server..."
    echo "SUITES are: ${SUITES}"
  fi
  set -e
}
