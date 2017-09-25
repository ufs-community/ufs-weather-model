set -ex

if [[ "$0" = "${BASH_SOURCE[0]}" ]]; then
  echo "$0 must be sourced"
  exit 1
fi

submit_and_wait() {

  [[ -z $1 ]] && exit 1

  local -r job_card=$1

  ROCOTO=${ROCOTO:-false}

  local test_status='PASS'

  if [[ $SCHEDULER = 'moab' ]]; then
    msub $job_card
  elif [[ $SCHEDULER = 'pbs' ]]; then
    qsubout=$( qsub $job_card )
    re='^([0-9]+\.[a-zA-Z0-9]+)$'
    qsub_id=0
    [[ "${qsubout}" =~ $re ]] && qsub_id=${BASH_REMATCH[1]}
  elif [[ $SCHEDULER = 'lsf' ]]; then
    bsubout=$( bsub < $job_card )
    re='Job <([0-9]+)> is submitted to queue <(.+)>.'
    bsub_id=0
    [[ "${bsubout}" =~ $re ]] && bsub_id=${BASH_REMATCH[1]}
  fi

  # wait for the job to enter the queue
  local count=0
  local job_running=0
  until [[ $job_running -eq 1 ]]
  do
    echo "TEST is waiting to enter the queue"
    if [[ $SCHEDULER = 'moab' ]]; then
      job_running=$( showq -u ${USER} -n | grep ${JBNME} | wc -l); sleep 5
    elif [[ $SCHEDULER = 'pbs' ]]; then
      job_running=$( qstat -u ${USER} -n | grep ${JBNME} | wc -l); sleep 5
    elif [[ $SCHEDULER = 'lsf' ]]; then
      job_running=$( bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | wc -l); sleep 5
    fi
    (( count=count+1 ))
    if [[ $count -eq 13 ]]; then echo "No job in queue after one minute, exiting..."; exit 2; fi
  done

  # find jobid
  if [[ $SCHEDULER = 'moab' ]]; then
     :
  elif [[ $SCHEDULER = 'pbs' ]]; then
    jobid=$( qstat -u ${USER} | grep ${JBNME} | awk '{print $1}' )
    trap 'echo "Job ${jobid} killed"; qdel ${jobid}; trap 0; exit' 1 2 3 4 5 6 7 8 10 12 13 15
    if [[ ${qsub_id} != ${jobid} ]]; then
      echo "Warning: qsub_id is not equal to jobid"
    fi
  elif [[ $SCHEDULER = 'lsf' ]]; then
    jobid=$( bjobs -u ${USER} -J ${JBNME} -noheader -o "jobid" )
    trap 'echo "Job ${jobid} killed"; bkill ${jobid}; trap 0; exit' 1 2 3 4 5 6 7 8 10 12 13 15
    if [[ ${bsub_id} -ne ${jobid} ]]; then
      echo "Warning: bsub_id is not equal to jobid"
    fi
  fi

  # wait for the job to finish and compare results
  job_running=1
  local n=1
  until [[ $job_running -eq 0 ]]
  do

    sleep 60 & wait $!
    if [[ $SCHEDULER = 'moab' ]]; then
      job_running=$( showq -u ${USER} -n | grep ${JBNME} | wc -l)
    elif [[ $SCHEDULER = 'pbs' ]]; then
      job_running=$( qstat -u ${USER} -n | grep ${JBNME} | wc -l)
    elif [[ $SCHEDULER = 'lsf' ]]; then
      job_running=$( bjobs -u ${USER} -J ${JBNME} 2>/dev/null | wc -l)
    fi

    if [[ $SCHEDULER = 'moab' ]]; then

      status=$( showq -u ${USER} -n | grep ${JBNME} | awk '{print $3}'); status=${status:--}
      if [[ -f ${RUNDIR}/err ]] ; then FnshHrs=$( grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $9 }'); fi
      FnshHrs=${FnshHrs:-0}
      if   [[ $status = 'Idle' ]];       then echo "$n min. TEST ${TEST_NR} is waiting in a queue, Status: $status"
      elif [[ $status = 'Running' ]];    then echo "$n min. TEST ${TEST_NR} is running,            Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'Starting' ]];   then echo "$n min. TEST ${TEST_NR} is ready to run,       Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'Completed' ]];  then echo "$n min. TEST ${TEST_NR} is finished,           Status: $status" ; job_running=0
      else                                    echo "$n min. TEST ${TEST_NR} is finished,           Status: $status , Finished $FnshHrs hours"
      fi

    elif [[ $SCHEDULER = 'pbs' ]]; then

      #status=$( qstat -u ${USER} -n | grep ${JBNME} | awk '{print $"10"}' ); status=${status:--}  PJP comment out to speed up regression test
      status=$( qstat -u ${USER} -n | grep ${JBNME} | awk '{print $10}' ); status=${status:--}
      if [[ -f ${RUNDIR}/err ]] ; then FnshHrs=$( tail -100 ${RUNDIR}/err | grep Finished | tail -1 | awk '{ print $9 }' ); fi
      FnshHrs=${FnshHrs:-0}
      if   [[ $status = 'Q' ]];  then echo "$n min. TEST ${TEST_NR} is waiting in a queue, Status: $status jobid ${jobid}"
      elif [[ $status = 'H' ]];  then echo "$n min. TEST ${TEST_NR} is held in a queue,    Status: $status"
      elif [[ $status = 'R' ]];  then echo "$n min. TEST ${TEST_NR} is running,            Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'E' ]] || [[ $status = 'C' ]];  then
        jobid=$( qstat -u ${USER} | grep ${JBNME} | awk '{print $1}')
        exit_status=$( qstat ${jobid} -f | grep exit_status | awk '{print $3}')
        if [[ $exit_status != 0 ]]; then
          echo "Test ${TEST_NR} FAIL" >> ${REGRESSIONTEST_LOG}
          echo                        >> ${REGRESSIONTEST_LOG}
          echo "Test ${TEST_NR} FAIL"
          echo
          test_status='FAIL'
          break
        fi
        echo "$n min. TEST ${TEST_NR} is finished,           Status: $status"
        job_running=0
      elif [[ $status = 'C' ]];  then echo "$n min. TEST ${TEST_NR} is finished,           Status: $status" ; job_running=0
      else                            echo "$n min. TEST ${TEST_NR} is finished,           Status: $status , Finished $FnshHrs hours"
      fi

    elif [[ $SCHEDULER = 'lsf' ]]; then

      status=$( bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | awk '{print $3}' ); status=${status:--}
      if [[ -f ${RUNDIR}/err ]] ; then FnshHrs=$( grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $9 }' ) ; fi
      FnshHrs=${FnshHrs:-0}
      if   [[ $status = 'PEND' ]];  then echo "$n min. TEST ${TEST_NR} is waiting in a queue, Status: $status"
      elif [[ $status = 'RUN'  ]];  then echo "$n min. TEST ${TEST_NR} is running,            Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'EXIT' ]];  then
        echo "Test ${TEST_NR} FAIL" >> ${REGRESSIONTEST_LOG}
        echo;echo;echo              >> ${REGRESSIONTEST_LOG}
        echo "Test ${TEST_NR} FAIL"
        echo;echo;echo
        test_status='FAIL'
        break
      else                               echo "$n min. TEST ${TEST_NR} is finished,           Status: $status , Finished $FnshHrs hours"
        exit_status=$( bjobs -u ${USER} -J ${JBNME} -a 2>/dev/null | grep $QUEUE | awk '{print $3}' )
        if [[ $exit_status = 'EXIT' ]];  then
          echo "Test ${TEST_NR} FAIL" >> ${REGRESSIONTEST_LOG}
          echo;echo;echo              >> ${REGRESSIONTEST_LOG}
          echo "Test ${TEST_NR} FAIL"
          echo;echo;echo
          test_status='FAIL'
          break
        fi
      fi

    fi
    (( n=n+1 ))
  done

  if [[ $test_status = 'FAIL' ]]; then
    echo $TEST_NAME >> $PATHRT/fail_test
    if [[ $ROCOTO == true ]]; then
      exit 2
    fi
  fi

}

check_results() {

  ROCOTO=${ROCOTO:-false}

  local test_status='PASS'

  # Give one minute for data to show up on file system
  #sleep 60

  echo                                          >> ${REGRESSIONTEST_LOG}
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}"    >> ${REGRESSIONTEST_LOG}
  echo "working dir  = ${RUNDIR}"               >> ${REGRESSIONTEST_LOG}
  echo "Checking test ${TEST_NR} results ...."  >> ${REGRESSIONTEST_LOG}
  echo
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}"
  echo "working dir  = ${RUNDIR}"
  echo "Checking test ${TEST_NR} results ...."

  if [[ ${CREATE_BASELINE} = false ]]; then
    #
    # --- regression test comparison ----
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

      else

        d=$( cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i | wc -l )

        if [[ $d -ne 0 ]] ; then
          echo ".......NOT OK" >> ${REGRESSIONTEST_LOG}
          echo ".......NOT OK"
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
    echo;echo;echo "Moving set ${TEST_NR} files ...."
    if [[ ! -d ${NEW_BASELINE}/${CNTL_DIR} ]]; then
      echo " mkdir -p ${NEW_BASELINE}/${CNTL_DIR}" >> ${REGRESSIONTEST_LOG}
      mkdir -p ${NEW_BASELINE}/${CNTL_DIR}/RESTART
    fi

    for i in ${LIST_FILES} ; do
      printf %s " Moving " $i " ....."   >> ${REGRESSIONTEST_LOG}
      if [[ -f ${RUNDIR}/$i ]] ; then
        cp ${RUNDIR}/${i} ${NEW_BASELINE}/${CNTL_DIR}/${i}  
      else
        echo "Missing " ${RUNDIR}/$i " output file"
        echo;echo " Set ${TEST_NR} failed"
        test_status='FAIL'
      fi
    done

  fi

  echo "Test ${TEST_NR} ${test_status}" >> ${REGRESSIONTEST_LOG}
  echo                                  >> ${REGRESSIONTEST_LOG}
  echo "Test ${TEST_NR} ${test_status}"
  echo

  if [[ $test_status = 'FAIL' ]]; then
    echo $TEST_NAME >> $PATHRT/fail_test
    if [[ $ROCOTO = true ]]; then
      exit 2
    fi
  fi

}

kill_job() {

  [[ -z $1 ]] && exit 1

  local -r jobid=$1

  if [[ $SCHEDULER = 'moab' ]]; then
     :
  elif [[ $SCHEDULER = 'pbs' ]]; then
    qdel ${jobid}
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

  if [[ "Q$APP" != Q ]] ; then
      rocoto_cmd="&PATHRT;/appbuild.sh &PATHTR; $APP $COMPILE_NR"
  else
      rocoto_cmd='&PATHRT;/compile.sh &PATHTR; $MACHINE_ID "${NEMS_VER}" $COMPILE_NR'
  fi

  if [[ ${COMPILE_NR_DEP} -gt 0 ]]; then
    cat << EOF >> $ROCOTO_XML
  <task name="compile_${COMPILE_NR}" maxtries="1">
    <dependency> <taskdep task="compile_${COMPILE_NR_DEP}"/> </dependency>
    <command>$rocoto_cmd</command>
    <jobname>compile_${COMPILE_NR}</jobname>
    <account>${ACCNR}</account>
    <queue>${COMPILE_QUEUE}</queue>
    <cores>1</cores>
    <memory>10G</memory>
    <walltime>01:00:00</walltime>
    <join>&LOG;/compile_${COMPILE_NR}.log</join>
    <envar> <name>SITE</name> <value>${MACHINE_ID}</value> </envar>
  </task>
EOF
  else
    cat << EOF >> $ROCOTO_XML
  <task name="compile_${COMPILE_NR}" maxtries="1">
    <command>$rocoto_cmd</command>
    <jobname>compile_${COMPILE_NR}</jobname>
    <account>${ACCNR}</account>
    <queue>${COMPILE_QUEUE}</queue>
    <memory>10G</memory>
    <cores>1</cores>
    <walltime>01:00:00</walltime>
    <join>&LOG;/compile_${COMPILE_NR}.log</join>
    <envar> <name>SITE</name> <value>${MACHINE_ID}</value> </envar>
  </task>
EOF
  fi
}


rocoto_create_run_task() {

  if [[ $DEP_RUN != '' ]]; then
    DEP_STRING="<and> <taskdep task=\"compile_${COMPILE_NR}\"/> <taskdep task=\"${DEP_RUN}\"/> </and>"
  else
    DEP_STRING="<taskdep task=\"compile_${COMPILE_NR}\"/>"
  fi

  NATIVE=""
  if [[ ${MACHINE_ID} == wcoss ]]; then
    NATIVE="<native>-a poe -R span[ptile=${TPN}]</native>"
  fi

  cat << EOF >> $ROCOTO_XML
    <task name="${TEST_NAME}" maxtries="1">
      <dependency> $DEP_STRING </dependency>
      <command>&PATHRT;/run_test.sh &PATHRT; &RUNDIR_ROOT; ${TEST_NAME} ${TEST_NR} ${COMPILE_NR} </command>
      <jobname>${TEST_NAME}</jobname>
      <account>${ACCNR}</account>
      <queue>${QUEUE}</queue>
      <memory></memory>
      <cores>${TASKS}</cores>
      <walltime>00:${WLCLK}:00</walltime>
      <join>&LOG;/run_${TEST_NAME}.log</join>
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
    $ROCOTORUN -w $ROCOTO_XML -d $ROCOTO_DB
    sleep 10 & wait $!
    state=$($ROCOTOSTAT -w $ROCOTO_XML -d $ROCOTO_DB -s | grep 197001010000 | awk -F" " '{print $2}')
    sleep 20 & wait $!
  done

}

ecflow_create_compile_task() {

  new_compile=true

  if [[ "Q$APP" != Q ]] ; then
      ecflow_cmd="$PATHRT/appbuild.sh $PATHTR $APP $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1"
  else
      ecflow_cmd="$PATHRT/compile.sh $PATHTR $MACHINE_ID \"${NEMS_VER}\" $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1"
  fi

  cat << EOF > ${ECFLOW_RUN}/regtest/compile_${COMPILE_NR}.ecf
%include <head.h>
$PATHRT/compile.sh $PATHTR $MACHINE_ID "${NEMS_VER}" $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1
%include <tail.h>
EOF

  echo "  task compile_${COMPILE_NR}" >> ${ECFLOW_RUN}/regtest.def
  if [[ ${COMPILE_NR_DEP} -gt 0 ]]; then
    echo "    trigger compile_${COMPILE_NR_DEP} == complete"  >> ${ECFLOW_RUN}/regtest.def
  fi

}

ecflow_create_run_task() {

  cat << EOF > ${ECFLOW_RUN}/regtest/${TEST_NAME}.ecf
%include <head.h>
$PATHRT/run_test.sh ${PATHRT} ${RUNDIR_ROOT} ${TEST_NAME} ${TEST_NR} ${COMPILE_NR} > ${LOG_DIR}/run_${TEST_NAME}.log 2>&1 &
%include <tail.h>
EOF

  echo "    task ${TEST_NAME}" >> ${ECFLOW_RUN}/regtest.def
  echo "      inlimit max_jobs" >> ${ECFLOW_RUN}/regtest.def
  if [[ $DEP_RUN != '' ]]; then
    echo "      trigger compile_${COMPILE_NR} == complete and ${DEP_RUN} == complete" >> ${ECFLOW_RUN}/regtest.def
  else
    echo "      trigger compile_${COMPILE_NR} == complete" >> ${ECFLOW_RUN}/regtest.def
  fi
}

ecflow_run() {

  ECF_PORT=$(( $(id -u) + 1500 ))
  # in rare instances when UID is greater then 58500 (like Ratko's UID on theia)
  [[ $ECF_PORT -ge 60000 ]] && ECF_PORT=12179

  ECF_NODE=$( hostname )

  set +e
  i=0
  max_atempts=5
  while [[ $i -lt $max_atempts ]]; do
    ecflow_client --ping --host=${ECF_NODE} --port=${ECF_PORT}
    not_running=$?
    if [[ $not_running -eq 1 ]]; then
      echo "ecflow_server is NOT running on ${ECF_NODE}:${ECF_PORT}"
      sh ${ECFLOW_START} -d ${ECFLOW_RUN} -p ${ECF_PORT} >> ${ECFLOW_RUN}/ecflow.log 2>&1
      break
    else
      echo "ecflow_server is already running on ${ECF_NODE}:${ECF_PORT}"
      ECF_PORT=$(( ECF_PORT + 1 ))
    fi
    i=$(( i + 1 ))
    if [[ $i -eq $max_atempts ]]; then
       echo "You already have $max_atempts ecFlow servers running on this node"
       exit 1
    fi
  done
  set -e

  ECFLOW_RUNNING=true
 
  export ECF_PORT
  export ECF_NODE

  ecflow_client --load=${ECFLOW_RUN}/regtest.def            >> ${ECFLOW_RUN}/ecflow.log 2>&1
  ecflow_client --begin=regtest                             >> ${ECFLOW_RUN}/ecflow.log 2>&1

  active_tasks=1
  while [[ $active_tasks -ne 0 ]]
  do
    sleep 10 & wait $!
    active_tasks=$( ecflow_client --get_state /regtest | grep "task " | grep -E 'state:active|state:submitted|state:queued' | wc -l )
    ${PATHRT}/abort_dep_tasks.py
  done
}

ecflow_kill() {
   [[ ${ECFLOW_RUNNING:-false} == true ]] || return
   set +e
   wait
   ecflow_client --kill /regtest
   sleep 10
}

ecflow_stop() {
   [[ ${ECFLOW_RUNNING:-false} == true ]] || return
   set +e
   wait
   ecflow_client --halt=yes
   ecflow_client --check_pt
   ecflow_client --terminate=yes
}
