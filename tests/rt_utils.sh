set -eu

if [[ "$0" = "${BASH_SOURCE[0]}" ]]; then
  echo "$0 must be sourced"
  exit 1
fi

OPNREQ_TEST=${OPNREQ_TEST:-false}

qsub_id=0
slurm_id=0
bsub_id=0

interrupt_job() {
  set -x
  if [[ $SCHEDULER = 'pbs' ]]; then
    echo "run_util.sh: interrupt_job qsub_id = ${qsub_id}"
    qdel ${qsub_id}
  elif [[ $SCHEDULER = 'slurm' ]]; then
    echo "run_util.sh: interrupt_job slurm_id = ${slurm_id}"
    scancel ${slurm_id}
  elif [[ $SCHEDULER = 'lsf' ]]; then
    echo "run_util.sh: interrupt_job bsub_id = ${bsub_id}"
    bkill ${bsub_id}
  else
    echo "run_util.sh: interrupt_job unknown SCHEDULER $SCHEDULER"
  fi
}

submit_and_wait() {

  [[ -z $1 ]] && exit 1

  [ -o xtrace ] && set_x='set -x' || set_x='set +x'
  set +x

  local -r job_card=$1

  ROCOTO=${ROCOTO:-false}
  ECFLOW=${ECFLOW:-false}

  local test_status='PASS'

  if [[ $SCHEDULER = 'pbs' ]]; then
    qsubout=$( qsub $job_card )
    re='^([0-9]+)(\.[a-zA-Z0-9\.-]+)$'
    [[ "${qsubout}" =~ $re ]] && qsub_id=${BASH_REMATCH[1]}
    echo "Job id ${qsub_id}"
  elif [[ $SCHEDULER = 'slurm' ]]; then
    slurmout=$( sbatch $job_card )
    re='Submitted batch job ([0-9]+)'
    [[ "${slurmout}" =~ $re ]] && slurm_id=${BASH_REMATCH[1]}
    echo "Job id ${slurm_id}"
  elif [[ $SCHEDULER = 'lsf' ]]; then
    bsubout=$( bsub < $job_card )
    re='Job <([0-9]+)> is submitted to queue <(.+)>.'
    [[ "${bsubout}" =~ $re ]] && bsub_id=${BASH_REMATCH[1]}
    echo "Job id ${bsub_id}"
  else
    echo "Unknown SCHEDULER $SCHEDULER"
    exit 1
  fi

  # wait for the job to enter the queue
  local count=0
  local job_running=0
  until [[ $job_running -eq 1 ]]
  do
    echo "TEST ${TEST_NR} ${TEST_NAME} is waiting to enter the queue"
    [[ ${ECFLOW:-false} == true ]] && ecflow_client --label=job_status "waiting to enter the queue"
    if [[ $SCHEDULER = 'pbs' ]]; then
      job_running=$( qstat ${qsub_id} | grep ${qsub_id} | wc -l )
    elif [[ $SCHEDULER = 'slurm' ]]; then
      job_running=$( squeue -u ${USER} -j ${slurm_id} | grep ${slurm_id} | wc -l)
    elif [[ $SCHEDULER = 'lsf' ]]; then
      job_running=$( bjobs ${bsub_id} | grep ${bsub_id} | wc -l)
    else
      echo "Unknown SCHEDULER $SCHEDULER"
      exit 1
    fi
    sleep 5
    (( count=count+1 ))
    if [[ $count -eq 13 ]]; then echo "No job in queue after one minute, exiting..."; exit 2; fi
  done

  # find jobid
  if [[ $SCHEDULER = 'pbs' ]]; then
    jobid=${qsub_id}
  elif [[ $SCHEDULER = 'slurm' ]]; then
    jobid=${slurm_id}
  elif [[ $SCHEDULER = 'lsf' ]]; then
    jobid=${bsub_id}
  else
    echo "Unknown SCHEDULER $SCHEDULER"
    exit 1
  fi
  echo "TEST ${TEST_NR} ${TEST_NAME} is submitted "
  if [[ ${ECFLOW:-false} == true ]]; then
    ecflow_client --label=job_id "${jobid}"
    ecflow_client --label=job_status "submitted"
  fi

  # wait for the job to finish and compare results
  job_running=1
  local n=1
  until [[ $job_running -eq 0 ]]
  do

    if [[ $SCHEDULER = 'pbs' ]]; then
      job_running=$( qstat ${qsub_id} | grep ${qsub_id} | wc -l )
    elif [[ $SCHEDULER = 'slurm' ]]; then
      job_running=$( squeue -u ${USER} -j ${slurm_id} | grep ${slurm_id} | wc -l)
    elif [[ $SCHEDULER = 'lsf' ]]; then
      job_running=$( bjobs ${bsub_id} | grep ${bsub_id} | wc -l)
    else
      echo "Unknown SCHEDULER $SCHEDULER"
      exit 1
    fi

    if [[ $SCHEDULER = 'pbs' ]]; then

      status=$( qstat ${qsub_id} | grep ${qsub_id} | awk '{print $5}' ); status=${status:--}
      if   [[ $status = 'Q' ]];  then
        status_label='waiting in a queue'
      elif [[ $status = 'H' ]];  then
        status_label='held in a queue'
      elif [[ $status = 'R' ]];  then
        status_label='running'
      elif [[ $status = 'E' ]] || [[ $status = 'C' ]] || [[ $status = '-' ]];  then
        status_label='finished'
        test_status='DONE'
        exit_status=$( qstat ${jobid} -x -f | grep Exit_status | awk '{print $3}')
        if [[ $exit_status != 0 ]]; then
          test_status='FAIL'
        fi
      else
        status_label='finished'
      fi

    elif [[ $SCHEDULER = 'slurm' ]]; then

      status=$( squeue -u ${USER} -j ${slurm_id} 2>/dev/null | grep ${slurm_id} | awk '{print $5}' ); status=${status:--}
      if   [[ $status = 'R'  ]];  then
        status_label='running'
      elif [[ $status = 'PD' ]];  then
        status_label='pending'
      elif [[ $status = 'F'  ]];  then
        status_label='failed'
        test_status='FAIL'
      elif [[ $status = 'C'  ]];  then
        status_label='finished'
        test_status='DONE'
      else
        echo "Slurm unknown status ${status}. Check sacct ..."
        sacct -n -j ${slurm_id} --format=JobID,state%20,Jobname%20
        status_label=$( sacct -n -j ${slurm_id} --format=JobID,state%20,Jobname%20 | grep "^${slurm_id}" | grep ${JBNME} | awk '{print $2}' )
        if [[ $status_label = 'FAILED' ]] || [[ $status_label = 'TIMEOUT' ]] || [[ $status_label = 'CANCELLED' ]] ; then
            test_status='FAIL'
        fi
      fi

    elif [[ $SCHEDULER = 'lsf' ]]; then

      status=$( bjobs ${bsub_id} 2>/dev/null | grep ${bsub_id} | awk '{print $3}' ); status=${status:--}
      if   [[ $status = 'PEND' ]];  then
        status_label='pending'
      elif [[ $status = 'RUN'  ]];  then
        status_label='running'
      elif [[ $status = 'DONE'  ]];  then
        status_label='finished'
        test_status='DONE'
      elif [[ $status = 'EXIT' ]];  then
        status_label='failed'
        test_status='FAIL'
      else
        echo "bsub unknown status ${status}"
        status_label='finished'
        test_status='DONE'
        exit_status=$( bjobs ${bsub_id} 2>/dev/null | grep ${bsub_id} | awk '{print $3}' ); status=${status:--}
        if [[ $exit_status = 'EXIT' ]];  then
          status_label='failed'
          test_status='FAIL'
        fi
      fi

    else
      echo "Unknown SCHEDULER $SCHEDULER"
      exit 1

    fi

    echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is ${status_label},  status: $status jobid ${jobid}"
    [[ ${ECFLOW:-false} == true ]] && ecflow_client --label=job_status "$status_label"

    if [[ $test_status = 'FAIL' || $test_status = 'DONE' ]]; then
      break
    fi

    (( n=n+1 ))
    sleep 60 & wait $!
  done

  if [[ $test_status = 'FAIL' ]]; then
    if [[ ${OPNREQ_TEST} == false ]]; then
      echo "Test ${TEST_NR} ${TEST_NAME} FAIL" >> ${REGRESSIONTEST_LOG}
      echo;echo;echo                           >> ${REGRESSIONTEST_LOG}
      echo "Test ${TEST_NR} ${TEST_NAME} FAIL"
    else
      echo ${TEST_NR} $TEST_NAME >> $PATHRT/fail_opnreq_test
    fi

    if [[ $ROCOTO == true || $ECFLOW == true ]]; then
      exit 1
    fi
  fi

  eval "$set_x"
}

check_results() {

  [ -o xtrace ] && set_x='set -x' || set_x='set +x'
  set +x

  ROCOTO=${ROCOTO:-false}
  ECFLOW=${ECFLOW:-false}

  local test_status='PASS'

  # Give one minute for data to show up on file system
  #sleep 60

  echo                                                       >  ${REGRESSIONTEST_LOG}
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}"                 >> ${REGRESSIONTEST_LOG}
  echo "working dir  = ${RUNDIR}"                            >> ${REGRESSIONTEST_LOG}
  echo "Checking test ${TEST_NR} ${TEST_NAME} results ...."  >> ${REGRESSIONTEST_LOG}
  echo
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}"
  echo "working dir  = ${RUNDIR}"
  echo "Checking test ${TEST_NR} ${TEST_NAME} results ...."

  if [[ ${CREATE_BASELINE} = false ]]; then
    #
    # --- regression test comparison
    #
    for i in ${LIST_FILES} ; do
      printf %s " Comparing " $i " ....." >> ${REGRESSIONTEST_LOG}
      printf %s " Comparing " $i " ....."

      if [[ ! -f ${RUNDIR}/$i ]] ; then

        echo ".......MISSING file" >> ${REGRESSIONTEST_LOG}
        echo ".......MISSING file"
        test_status='FAIL'

      elif [[ ! -f ${RTPWD}/${CNTL_DIR}/$i ]] ; then

        echo ".......MISSING baseline" >> ${REGRESSIONTEST_LOG}
        echo ".......MISSING baseline"
        test_status='FAIL'

      elif [[ $RT_COMPILER == "gnu" && $i == "RESTART/fv_core.res.nc" ]] ; then

        # Although identical in ncdiff, RESTART/fv_core.res.nc differs in byte 469, line 3,
        # for the fv3_control_32bit test between each run (without changing the source code)
        # for GNU compilers - skip comparison.
        echo ".......SKIP for gnu compilers" >> ${REGRESSIONTEST_LOG}
        echo ".......SKIP for gnu compilers"

      else

        cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i >/dev/null 2>&1 && d=$? || d=$?
        if [[ $d -eq 2 ]]; then
          echo "....CMP ERROR" >> ${REGRESSIONTEST_LOG}
          echo "....CMP ERROR"
          exit 1
        fi

        if [[ $d -eq 1 && ${i##*.} == 'nc' ]] ; then
          if [[ ${MACHINE_ID} =~ orion || ${MACHINE_ID} =~ hera || ${MACHINE_ID} =~ wcoss_dell_p3 || ${MACHINE_ID} =~ wcoss_cray || ${MACHINE_ID} =~ cheyenne || ${MACHINE_ID} =~ gaea || ${MACHINE_ID} =~ jet || ${MACHINE_ID} =~ s4 ]] ; then
            printf ".......ALT CHECK.." >> ${REGRESSIONTEST_LOG}
            printf ".......ALT CHECK.."
            ${PATHRT}/compare_ncfile.py ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i >/dev/null 2>&1 && d=$? || d=$?
            if [[ $d -eq 1 ]]; then
              echo "....ERROR" >> ${REGRESSIONTEST_LOG}
              echo "....ERROR"
              exit 1
            fi
          fi
        fi

        if [[ $d -ne 0 ]]; then
          echo "....NOT OK" >> ${REGRESSIONTEST_LOG}
          echo "....NOT OK"
          test_status='FAIL'
        else
          echo "....OK" >> ${REGRESSIONTEST_LOG}
          echo "....OK"
        fi

      fi

    done

  else
    #
    # --- create baselines
    #
    echo;echo "Moving baseline ${TEST_NR} ${TEST_NAME} files ...."
    echo;echo "Moving baseline ${TEST_NR} ${TEST_NAME} files ...." >> ${REGRESSIONTEST_LOG}

    for i in ${LIST_FILES} ; do
      printf %s " Moving " $i " ....."
      printf %s " Moving " $i " ....."   >> ${REGRESSIONTEST_LOG}
      if [[ -f ${RUNDIR}/$i ]] ; then
        mkdir -p ${NEW_BASELINE}/${CNTL_DIR}/$(dirname ${i})
        cp ${RUNDIR}/${i} ${NEW_BASELINE}/${CNTL_DIR}/${i}
        echo "....OK" >>${REGRESSIONTEST_LOG}
        echo "....OK"
      else
        echo "....NOT OK. Missing " ${RUNDIR}/$i >>${REGRESSIONTEST_LOG}
        echo "....NOT OK. Missing " ${RUNDIR}/$i
        test_status='FAIL'
      fi
    done

  fi

  echo                                               >> ${REGRESSIONTEST_LOG}
  grep "The total amount of wall time" ${RUNDIR}/out >> ${REGRESSIONTEST_LOG}
  echo                                               >> ${REGRESSIONTEST_LOG}

  TRIES=''
  if [[ $ECFLOW == true ]]; then
    if [[ $ECF_TRYNO -gt 1 ]]; then
      TRIES=" Tries: $ECF_TRYNO"
    fi
  fi
  echo "Test ${TEST_NR} ${TEST_NAME} ${test_status}${TRIES}" >> ${REGRESSIONTEST_LOG}
  echo                                                       >> ${REGRESSIONTEST_LOG}
  echo "Test ${TEST_NR} ${TEST_NAME} ${test_status}${TRIES}"
  echo

  if [[ $test_status = 'FAIL' ]]; then
    if [[ ${OPNREQ_TEST} == false ]]; then
      echo "${TEST_NAME} ${TEST_NR} failed in check_result" >> $PATHRT/fail_test_${TEST_NR}
    else
      echo ${TEST_NR} $TEST_NAME >> $PATHRT/fail_opnreq_test
    fi

    if [[ $ROCOTO = true || $ECFLOW == true ]]; then
      exit 1
    fi
  fi

  eval "$set_x"
}


kill_job() {

  [[ -z $1 ]] && exit 1

  local -r jobid=$1

  if [[ $SCHEDULER = 'pbs' ]]; then
    qdel ${jobid}
  elif [[ $SCHEDULER = 'slurm' ]]; then
    scancel ${jobid}
  elif [[ $SCHEDULER = 'lsf' ]]; then
    bkill ${jobid}
  fi
}

rocoto_create_compile_task() {

  new_compile=true
  if [[ $in_metatask == true ]]; then
    in_metatask=false
    echo "  </metatask>" >> $ROCOTO_XML
  fi

  # serialize WW3 builds. FIXME
  DEP_STRING=""
  if [[ ${MAKE_OPT^^} =~ "WW3=Y" && ${COMPILE_PREV_WW3_NR} != '' ]]; then
    DEP_STRING="<dependency><taskdep task=\"compile_${COMPILE_PREV_WW3_NR}\"/></dependency>"
  fi

  NATIVE=""
  BUILD_CORES=8
  BUILD_WALLTIME="00:30:00"
  if [[ ${MACHINE_ID} == wcoss_dell_p3 ]]; then
    BUILD_CORES=1
    NATIVE="<memory>8G</memory> <native>-R 'affinity[core(1)]'</native>"
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == wcoss_cray ]]; then
    BUILD_CORES=24
    NATIVE="<exclusive></exclusive> <envar><name>PATHTR</name><value>&PATHTR;</value></envar>"
  fi
  if [[ ${MACHINE_ID} == jet ]]; then
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == orion.* ]]; then
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == s4.* ]]; then
    BUILD_WALLTIME="01:00:00"
  fi

  cat << EOF >> $ROCOTO_XML
  <task name="compile_${COMPILE_NR}" maxtries="3">
    $DEP_STRING
    <command>&PATHRT;/run_compile.sh &PATHRT; &RUNDIR_ROOT; "${MAKE_OPT}" ${COMPILE_NR}</command>
    <jobname>compile_${COMPILE_NR}</jobname>
    <account>${ACCNR}</account>
    <queue>${COMPILE_QUEUE}</queue>
    <partition>${PARTITION}</partition>
    <cores>${BUILD_CORES}</cores>
    <walltime>${BUILD_WALLTIME}</walltime>
    <stdout>&RUNDIR_ROOT;/compile_${COMPILE_NR}/out</stdout>
    <stderr>&RUNDIR_ROOT;/compile_${COMPILE_NR}/err</stderr>
    ${NATIVE}
  </task>
EOF
}

rocoto_create_run_task() {

  if [[ $DEP_RUN != '' ]]; then
    DEP_STRING="<and> <taskdep task=\"compile_${COMPILE_NR}\"/> <taskdep task=\"${DEP_RUN}${RT_SUFFIX}\"/> </and>"
  else
    DEP_STRING="<taskdep task=\"compile_${COMPILE_NR}\"/>"
  fi

  CORES=$(( ${TASKS} * ${THRD} ))
  if (( TPN > CORES )); then
    TPN=$CORES
  fi

  NATIVE=""
  if [[ ${MACHINE_ID} == wcoss_dell_p3 ]]; then
    NATIVE="<nodesize>28</nodesize><native>-R 'affinity[core(${THRD})]'</native>"
  fi
  if [[ ${MACHINE_ID} == wcoss_cray ]]; then
    NATIVE="<exclusive></exclusive>"
  fi

  cat << EOF >> $ROCOTO_XML
    <task name="${TEST_NAME}${RT_SUFFIX}" maxtries="1">
      <dependency> $DEP_STRING </dependency>
      <command>&PATHRT;/run_test.sh &PATHRT; &RUNDIR_ROOT; ${TEST_NAME} ${TEST_NR} ${COMPILE_NR} </command>
      <jobname>${TEST_NAME}${RT_SUFFIX}</jobname>
      <account>${ACCNR}</account>
      <queue>${QUEUE}</queue>
      <partition>${PARTITION}</partition>
      <nodes>${NODES}:ppn=${TPN}</nodes>
      <walltime>00:${WLCLK}:00</walltime>
      <stdout>&RUNDIR_ROOT;/${TEST_NAME}${RT_SUFFIX}/out</stdout>
      <stderr>&RUNDIR_ROOT;/${TEST_NAME}${RT_SUFFIX}/err</stderr>
      ${NATIVE}
    </task>
EOF

}

rocoto_kill() {
   for jobid in $( $ROCOTOSTAT -w $ROCOTO_XML -d $ROCOTO_DB | grep 197001010000 | grep -E 'QUEUED|RUNNING' | awk -F" " '{print $3}' ); do
      kill_job ${jobid}
   done
}

rocoto_run() {

  state="Active"
  while [[ $state != "Done" ]]
  do
    $ROCOTORUN -v 10 -w $ROCOTO_XML -d $ROCOTO_DB
    sleep 10
    state=$($ROCOTOSTAT -w $ROCOTO_XML -d $ROCOTO_DB -s | grep 197001010000 | awk -F" " '{print $2}')
    dead_compile=$($ROCOTOSTAT -w $ROCOTO_XML -d $ROCOTO_DB | grep compile_ | grep DEAD | head -1 | awk -F" " '{print $2}')
    if [[ ! -z ${dead_compile} ]]; then
      echo "y" | ${ROCOTOCOMPLETE} -w $ROCOTO_XML -d $ROCOTO_DB -m ${dead_compile}_tasks
      ${ROCOTOCOMPLETE} -w $ROCOTO_XML -d $ROCOTO_DB -t ${dead_compile}
    fi
    sleep 20
  done

}


ecflow_create_compile_task() {

  new_compile=true


  cat << EOF > ${ECFLOW_RUN}/${ECFLOW_SUITE}/compile_${COMPILE_NR}.ecf
%include <head.h>
$PATHRT/run_compile.sh ${PATHRT} ${RUNDIR_ROOT} "${MAKE_OPT}" $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1 &
%include <tail.h>
EOF

  echo "  task compile_${COMPILE_NR}" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      label build_options '${MAKE_OPT}'" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      label job_id ''" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      label job_status ''" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      inlimit max_builds" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  # serialize WW3 builds. FIXME
  if [[ ${MAKE_OPT^^} =~ "-DAPP=ATMW" ]] || [[ ${MAKE_OPT^^} =~ "-DAPP=S2SW" ]] || [[ ${MAKE_OPT^^} =~ "-DAPP=HAFSW" ]] || [[ ${MAKE_OPT^^} =~ "-DAPP=HAFS-ALL" ]] ; then
    if [[ ${COMPILE_PREV_WW3_NR} != '' ]]; then
      echo "    trigger compile_${COMPILE_PREV_WW3_NR} == complete"  >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
    fi
  fi
}

ecflow_create_run_task() {

  cat << EOF > ${ECFLOW_RUN}/${ECFLOW_SUITE}/${TEST_NAME}${RT_SUFFIX}.ecf
%include <head.h>
$PATHRT/run_test.sh ${PATHRT} ${RUNDIR_ROOT} ${TEST_NAME} ${TEST_NR} ${COMPILE_NR} > ${LOG_DIR}/run_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log 2>&1 &
%include <tail.h>
EOF

  echo "    task ${TEST_NAME}${RT_SUFFIX}" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      label job_id ''" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      label job_status ''" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      inlimit max_jobs" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  if [[ $DEP_RUN != '' ]]; then
    if [[ ${OPNREQ_TEST} == false ]]; then
      echo "      trigger compile_${COMPILE_NR} == complete and ${DEP_RUN}${RT_SUFFIX} == complete" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
    else
      echo "      trigger compile_${COMPILE_NR} == complete and ${DEP_RUN} == complete" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
    fi
  else
    echo "      trigger compile_${COMPILE_NR} == complete" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  fi
}

ecflow_run() {

  # in rare instances when UID is greater then 58500 (like Ratko's UID on theia)
  [[ $ECF_PORT -gt 49151 ]] && ECF_PORT=12179

  ECF_HOST=$( hostname )

  set +e
  ecflow_client --ping --host=${ECF_HOST} --port=${ECF_PORT}
  not_running=$?
  if [[ $not_running -eq 1 ]]; then
    echo "ecflow_server is NOT running on ${ECF_HOST}:${ECF_PORT}"
    if [[ ${MACHINE_ID} == wcoss2 ]]; then
      # Annoying "Has NCO assigned port $ECF_PORT for use by this account? (yes/no) ".
      echo yes | ${ECFLOW_START} -p ${ECF_PORT}
    else
      ${ECFLOW_START} -p ${ECF_PORT}
    fi
  else
    echo "ecflow_server is already running on ${ECF_HOST}:${ECF_PORT}"
  fi
  set -e

  ECFLOW_RUNNING=true

  export ECF_PORT
  export ECF_HOST

  ecflow_client --load=${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  ecflow_client --begin=${ECFLOW_SUITE}
  ecflow_client --restart

  active_tasks=1
  while [[ $active_tasks -ne 0 ]]
  do
    sleep 10 & wait $!
    active_tasks=$( ecflow_client --get_state /${ECFLOW_SUITE} | grep "task " | grep -E 'state:active|state:submitted|state:queued' | wc -l )
    echo "ecflow tasks remaining: ${active_tasks}"
    ${PATHRT}/abort_dep_tasks.py
  done
  sleep 65 # wait one ECF_INTERVAL plus 5 seconds
  ecflow_client --delete=yes /${ECFLOW_SUITE}
  sleep 5
}

ecflow_kill() {
   [[ ${ECFLOW_RUNNING:-false} == true ]] || return
   set +e
   ecflow_client --suspend /${ECFLOW_SUITE}
   ecflow_client --kill /${ECFLOW_SUITE}
   sleep 20
   ecflow_client --delete=force yes /${ECFLOW_SUITE}
}

ecflow_stop() {
   [[ ${ECFLOW_RUNNING:-false} == true ]] || return
   set +e
   SUITES=$( ecflow_client --get | grep "^suite" )
   echo "SUITES=${SUITES}"
   if [ -z "${SUITES}" ]; then
     ecflow_client --halt=yes
     ecflow_client --check_pt
     ecflow_client --terminate=yes
   fi
}
