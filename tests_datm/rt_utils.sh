#
# DH* TODO - COMBINE SBATCH AND SLURM?
#
set -eu

if [[ "$0" = "${BASH_SOURCE[0]}" ]]; then
  echo "$0 must be sourced"
  exit 1
fi

UNIT_TEST=${UNIT_TEST:-false}

submit_and_wait() {

  [[ -z $1 ]] && exit 1

  [ -o xtrace ] && set_x='set -x' || set_x='set +x'
  set +x

  local -r job_card=$1

  ROCOTO=${ROCOTO:-false}

  local test_status='PASS'

  if [[ $SCHEDULER = 'moab' ]]; then
    msub $job_card
  elif [[ $SCHEDULER = 'pbs' ]]; then
    qsubout=$( qsub $job_card )
    if [[ ${MACHINE_ID} = cheyenne.* ]]; then
      re='^([0-9]+\.[a-zA-Z0-9\.]+)$'
    else
      re='^([0-9]+\.[a-zA-Z0-9]+)$'
    fi
    qsub_id=0
    [[ "${qsubout}" =~ $re ]] && qsub_id=${BASH_REMATCH[1]}
    if [[ ${MACHINE_ID} = cheyenne.* ]]; then
      qsub_id="${qsub_id%.chadm*}"
    fi
  elif [[ $SCHEDULER = 'sbatch' ]]; then
    qsubout=$( sbatch $job_card )
    re='^([0-9]+\.[a-zA-Z0-9]+)$'
    qsub_id=0
    [[ "${qsubout}" =~ $re ]] && qsub_id=${BASH_REMATCH[1]}
    if [[ ${MACHINE_ID} = stampede.* ]]; then
      qsub_id="${qsub_id}"
    fi
  elif [[ $SCHEDULER = 'slurm' ]]; then
    slurmout=$( sbatch $job_card )
    re='Submitted batch job ([0-9]+)'
    slurm_id=0
    [[ "${slurmout}" =~ $re ]] && slurm_id=${BASH_REMATCH[1]}
  elif [[ $SCHEDULER = 'lsf' ]]; then
    bsubout=$( bsub < $job_card )
    re='Job <([0-9]+)> is submitted to queue <(.+)>.'
    bsub_id=0
    [[ "${bsubout}" =~ $re ]] && bsub_id=${BASH_REMATCH[1]}
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
    if [[ $SCHEDULER = 'moab' ]]; then
      job_running=$( showq -u ${USER} -n | grep ${JBNME} | wc -l); sleep 5
    elif [[ $SCHEDULER = 'pbs' ]]; then
      if [[ ${MACHINE_ID} = cheyenne.* ]]; then
        job_running=$( qstat ${qsub_id} | grep ${qsub_id} | wc -l); sleep 5
      else
        job_running=$( qstat -u ${USER} -n | grep ${JBNME} | wc -l); sleep 5
      fi
    elif [[ $SCHEDULER = 'sbatch' ]]; then
      if [[ ${MACHINE_ID} = stampede.* ]]; then
        job_running=$( squeue ${qsub_id} | grep ${qsub_id} | wc -l); sleep 5
      else
        job_running=$( squeue -u ${USER} -n | grep ${JBNME} | wc -l); sleep 5
      fi
    elif [[ $SCHEDULER = 'slurm' ]]; then
      job_running=$( squeue -u ${USER} -j ${slurm_id} | grep ${slurm_id} | wc -l); sleep 5
    elif [[ $SCHEDULER = 'lsf' ]]; then
      job_running=$( bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | wc -l); sleep 5
    else
      echo "Unknown SCHEDULER $SCHEDULER"
      exit 1
    fi
    (( count=count+1 ))
    if [[ $count -eq 13 ]]; then echo "No job in queue after one minute, exiting..."; exit 2; fi
  done

  # find jobid
  if [[ $SCHEDULER = 'moab' ]]; then
     :
  elif [[ $SCHEDULER = 'pbs' ]]; then
    if [[ ${MACHINE_ID} = cheyenne.* ]]; then
      jobid=$( qstat ${qsub_id} | grep ${qsub_id} | awk '{print $1}' )
      jobid="${jobid%.chadm*}"
    else
      jobid=$( qstat -u ${USER} | grep ${JBNME} | awk '{print $1}' )
    fi
    trap 'echo "Job ${jobid} killed"; qdel ${jobid}; trap 0; exit' 1 2 3 4 5 6 7 8 10 12 13 15
    if [[ ${qsub_id} != ${jobid} ]]; then
      echo "Warning: qsub_id is not equal to jobid"
    fi
  elif [[ $SCHEDULER = 'sbatch' ]]; then
    if [[ ${MACHINE_ID} = stampede.* ]]; then
      jobid=$( squeue ${qsub_id} | grep ${qsub_id} | awk '{print $1}' )
      jobid="${jobid}"
    else
      jobid=$( squeue -u ${USER} | grep ${JBNME} | awk '{print $1}' )
    fi
    trap 'echo "Job ${jobid} killed"; qdel ${jobid}; trap 0; exit' 1 2 3 4 5 6 7 8 10 12 13 15
    if [[ ${qsub_id} != ${jobid} ]]; then
      echo "Warning: qsub_id is not equal to jobid"
    fi
  elif [[ $SCHEDULER = 'slurm' ]]; then
    jobid=${slurm_id}
  elif [[ $SCHEDULER = 'lsf' ]]; then
    jobid=$( bjobs -u ${USER} -J ${JBNME} -noheader -o "jobid" )
    trap 'echo "Job ${jobid} killed"; bkill ${jobid}; trap 0; exit' 1 2 3 4 5 6 7 8 10 12 13 15
    if [[ ${bsub_id} -ne ${jobid} ]]; then
      echo "Warning: bsub_id is not equal to jobid"
    fi
  else
    echo "Unknown SCHEDULER $SCHEDULER"
    exit 1
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
      if [[ ${MACHINE_ID} = cheyenne.* ]]; then
        job_running=$( qstat ${qsub_id} | grep ${qsub_id} | wc -l); sleep 5
      else
        job_running=$( qstat -u ${USER} -n | grep ${JBNME} | wc -l)
      fi
    elif [[ $SCHEDULER = 'sbatch' ]]; then
      if [[ ${MACHINE_ID} = stampede.* ]]; then
        job_running=$( squeue ${qsub_id} | grep ${qsub_id} | wc -l); sleep 5
      else
        job_running=$( squeue -u ${USER} -n | grep ${JBNME} | wc -l)
      fi
    elif [[ $SCHEDULER = 'slurm' ]]; then
      job_running=$( squeue -u ${USER} -j ${slurm_id} | grep ${slurm_id} | wc -l)
    elif [[ $SCHEDULER = 'lsf' ]]; then
      job_running=$( bjobs -u ${USER} -J ${JBNME} 2>/dev/null | wc -l)
    else
      echo "Unknown SCHEDULER $SCHEDULER"
      exit 1
    fi

    if [[ $SCHEDULER = 'moab' ]]; then

      status=$( showq -u ${USER} -n | grep ${JBNME} | awk '{print $3}'); status=${status:--}
      if [[ -f ${RUNDIR}/err ]] ; then FnshHrs=$( grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $9 }'); fi
      FnshHrs=${FnshHrs:-0}
      if   [[ $status = 'Idle' ]];       then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is waiting in a queue, Status: $status"
      elif [[ $status = 'Running' ]];    then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is running,            Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'Starting' ]];   then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is ready to run,       Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'Completed' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is finished,           Status: $status" ; job_running=0
      else                                    echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is finished,           Status: $status , Finished $FnshHrs hours"
      fi

    elif [[ $SCHEDULER = 'pbs' ]]; then

      #status=$( qstat -u ${USER} -n | grep ${JBNME} | awk '{print $"10"}' ); status=${status:--}  PJP comment out to speed up regression test
      if [[ ${MACHINE_ID} = cheyenne.* ]]; then
        status=$( qstat ${qsub_id} | grep ${qsub_id} | awk '{print $5}' ); status=${status:--}
      else
        status=$( qstat -u ${USER} -n | grep ${JBNME} | awk '{print $10}' ); status=${status:--}
      fi
      if [[ -f ${RUNDIR}/err ]] ; then FnshHrs=$( tail -100 ${RUNDIR}/err | grep Finished | tail -1 | awk '{ print $9 }' ); fi
      FnshHrs=${FnshHrs:-0}
      if   [[ $status = 'Q' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is waiting in a queue, Status: $status jobid ${jobid}"
      elif [[ $status = 'H' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is held in a queue,    Status: $status"
      elif [[ $status = 'R' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is running,            Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'E' ]] || [[ $status = 'C' ]];  then
        if [[ ${MACHINE_ID} = cheyenne.* ]]; then
          exit_status=$( qstat ${jobid} -x -f | grep Exit_status | awk '{print $3}')
        else
          jobid=$( qstat -u ${USER} | grep ${JBNME} | awk '{print $1}')
          exit_status=$( qstat ${jobid} -f | grep exit_status | awk '{print $3}')
        fi
        if [[ $exit_status != 0 ]]; then
          echo "Test ${TEST_NR} ${TEST_NAME} FAIL" >> ${REGRESSIONTEST_LOG}
          echo                                     >> ${REGRESSIONTEST_LOG}
          echo "Test ${TEST_NR} ${TEST_NAME} FAIL"
          echo
          test_status='FAIL'
          break
        fi
        echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is finished,           Status: $status"
        job_running=0
      elif [[ $status = 'C' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is finished,           Status: $status" ; job_running=0
      else                            echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is finished,           Status: $status , Finished $FnshHrs hours"
      fi

    elif [[ $SCHEDULER = 'slurm' ]]; then

      status=$( squeue -u ${USER} -j ${slurm_id} 2>/dev/null | grep ${slurm_id} | awk '{print $5}' ); status=${status:--}
      if [[ $status = 'R' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is running,            Status: $status"
      elif [[ $status = 'F' ]];  then
        echo "Test ${TEST_NR} ${TEST_NAME} FAIL" >> ${REGRESSIONTEST_LOG}
        echo                                     >> ${REGRESSIONTEST_LOG}
        echo "Test ${TEST_NR} ${TEST_NAME} FAIL"
        echo
        echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is failed,             Status: $status"
        job_running=0
      elif [[ $status = 'C' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is finished,           Status: $status" ; job_running=0
      else
        state=$( sacct -n -j ${slurm_id}.batch --format=JobID,state,Jobname | grep ${slurm_id} | awk '{print $2}' )
        echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is ${state}"
      fi

    elif [[ $SCHEDULER = 'sbatch' ]]; then

      #status=$( qstat -u ${USER} -n | grep ${JBNME} | awk '{print $"10"}' ); status=${status:--}  PJP comment out to speed up regression test
      if [[ ${MACHINE_ID} = stampede.* ]]; then
        status=$( squeue ${qsub_id} | grep ${qsub_id} | awk '{print $5}' ); status=${status:--}
      else
        status=$( squeue -u ${USER} -n | grep ${JBNME} | awk '{print $10}' ); status=${status:--}
      fi
      if [[ -f ${RUNDIR}/err ]] ; then FnshHrs=$( tail -100 ${RUNDIR}/err | grep Finished | tail -1 | awk '{ print $9 }' ); fi
      FnshHrs=${FnshHrs:-0}
      if   [[ $status = 'Q' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is waiting in a queue, Status: $status jobid ${jobid}"
      elif [[ $status = 'H' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is held in a queue,    Status: $status"
      elif [[ $status = 'R' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is running,            Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'E' ]] || [[ $status = 'C' ]];  then
        if [[ ${MACHINE_ID} = stampede.* ]]; then
          exit_status=$( squeue ${jobid} -x -f | grep Exit_status | awk '{print $3}')
        else
          jobid=$( squeue -u ${USER} | grep ${JBNME} | awk '{print $1}')
          exit_status=$( qstat ${jobid} -f | grep exit_status | awk '{print $3}')
        fi
        if [[ $exit_status != 0 ]]; then
          echo "Test ${TEST_NR}  ${TEST_NAME} FAIL" >> ${REGRESSIONTEST_LOG}
          echo                                      >> ${REGRESSIONTEST_LOG}
          echo "Test ${TEST_NR}  ${TEST_NAME} FAIL"
          echo
          test_status='FAIL'
          break
        fi
        echo "$n min. TEST ${TEST_NR}  ${TEST_NAME} is finished,           Status: $status"
        job_running=0
      elif [[ $status = 'C' ]];  then echo "$n min. TEST ${TEST_NR}  ${TEST_NAME} is finished,           Status: $status" ; job_running=0
      else                            echo "$n min. TEST ${TEST_NR}  ${TEST_NAME} is finished,           Status: $status , Finished $FnshHrs hours"
      fi

    elif [[ $SCHEDULER = 'lsf' ]]; then

      status=$( bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | awk '{print $3}' ); status=${status:--}
      if [[ -f ${RUNDIR}/err ]] ; then FnshHrs=$( grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $9 }' ) ; fi
      FnshHrs=${FnshHrs:-0}
      if   [[ $status = 'PEND' ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is waiting in a queue, Status: $status"
      elif [[ $status = 'RUN'  ]];  then echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is running,            Status: $status , Finished $FnshHrs hours"
      elif [[ $status = 'EXIT' ]];  then
        echo "Test ${TEST_NR} ${TEST_NAME} FAIL" >> ${REGRESSIONTEST_LOG}
        echo;echo;echo                           >> ${REGRESSIONTEST_LOG}
        echo "Test ${TEST_NR} ${TEST_NAME} FAIL"
        echo;echo;echo
        test_status='FAIL'
        break
      else                               echo "$n min. TEST ${TEST_NR} ${TEST_NAME} is finished,           Status: $status , Finished $FnshHrs hours"
        exit_status=$( bjobs -u ${USER} -J ${JBNME} -a 2>/dev/null | grep $QUEUE | awk '{print $3}' )
        if [[ $exit_status = 'EXIT' ]];  then
          echo "Test ${TEST_NR} ${TEST_NAME} FAIL" >> ${REGRESSIONTEST_LOG}
          echo;echo;echo                           >> ${REGRESSIONTEST_LOG}
          echo "Test ${TEST_NR} ${TEST_NAME} FAIL"
          echo;echo;echo
          test_status='FAIL'
          break
        fi
      fi

    else
      echo "Unknown SCHEDULER $SCHEDULER"
      exit 1

    fi
    (( n=n+1 ))
  done

  if [[ $test_status = 'FAIL' ]]; then
    if [[ ${UNIT_TEST} == false ]]; then
      echo $TEST_NAME >> $PATHRT/fail_test
    else
      echo ${TEST_NR} $TEST_NAME >> $PATHRT/fail_unit_test
    fi

    if [[ $ROCOTO == true ]]; then
      exit 2
    fi
  fi

  eval "$set_x"
}

check_results() {

  [ -o xtrace ] && set_x='set -x' || set_x='set +x'
  set +x

  ROCOTO=${ROCOTO:-false}

  # Default compiler "intel"
  export COMPILER=${NEMS_COMPILER:-intel}

  local test_status='PASS'

  # Give one minute for data to show up on file system
  #sleep 60

  echo                                                       >  ${REGRESSIONTEST_LOG}
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}"                 >> ${REGRESSIONTEST_LOG}
  echo "mediator baseline dir = ${RTPWD}/${CNTLMED_DIR}"     >> ${REGRESSIONTEST_LOG}
  echo "working dir  = ${RUNDIR}"                            >> ${REGRESSIONTEST_LOG}
  echo "Checking test ${TEST_NR} ${TEST_NAME} results ...."  >> ${REGRESSIONTEST_LOG}
  echo
  echo "baseline dir = ${RTPWD}/${CNTL_DIR}"
  echo "mediator baseline dir = ${RTPWD}/${CNTLMED_DIR}"
  echo "working dir  = ${RUNDIR}"
  echo "Checking test ${TEST_NR} ${TEST_NAME} results ...."

  if [[ ${CREATE_BASELINE} = false ]]; then
    #
    # --- regression test comparison ----
    #
    for i in ${LIST_FILES} ; do
      printf %s " Comparing " $i " ....." >> ${REGRESSIONTEST_LOG}
      printf %s " Comparing " $i " ....."

      crst=''
      if [[ $i =~ RESTART/ ]]; then
        crst=RESTART/$(basename $i)
      fi

      if [[ ! -f ${RUNDIR}/$i ]] ; then

        echo ".......MISSING file" >> ${REGRESSIONTEST_LOG}
        echo ".......MISSING file"
        test_status='FAIL'

      elif [[ ! -f ${RTPWD}/${CNTL_DIR}/$i && ! -f ${RTPWD}/${CNTLMED_DIR}/$i && ! -f ${RTPWD}/${CNTL_DIR}/$crst ]] ; then

        echo ".......MISSING baseline" >> ${REGRESSIONTEST_LOG}
        echo ".......MISSING baseline"
        test_status='FAIL'

      elif [[ ( $COMPILER == "gnu" || $COMPILER == "pgi" ) ]] ; then
        echo ".......SKIP for gnu/pgi compilers" >> ${REGRESSIONTEST_LOG}
        echo ".......SKIP for gnu/pgi compilers"

      else

        if [[ $i =~ mediator ]]; then
          d=$( cmp ${RTPWD}/${CNTLMED_DIR}/$i ${RUNDIR}/$i | wc -l )
        elif [[ $i =~ RESTART/ ]]; then
          d=$( cmp ${RTPWD}/${CNTL_DIR}/$crst ${RUNDIR}/$i | wc -l )
        else
          d=$( cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i | wc -l )
        fi

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
    echo;echo;echo "Moving set ${TEST_NR} ${TEST_NAME} files ...."
    if [[ ! -d ${NEW_BASELINE}/${CNTL_DIR}/RESTART ]] ; then
      echo " mkdir -p ${NEW_BASELINE}/${CNTL_DIR}/RESTART" >> ${REGRESSIONTEST_LOG}
      mkdir -p ${NEW_BASELINE}/${CNTL_DIR}/RESTART
    fi
    if [[ ${CNTLMED_DIR} =~ MEDIATOR && ! -d ${NEW_BASELINE}/${CNTLMED_DIR} ]]; then
      echo " mkdir -p ${NEW_BASELINE}/${CNTLMED_DIR}" >> ${REGRESSIONTEST_LOG}
      mkdir -p ${NEW_BASELINE}/${CNTLMED_DIR}
    fi

    for i in ${LIST_FILES} ; do
      printf %s " Moving " $i " ....."   >> ${REGRESSIONTEST_LOG}
      if [[ -f ${RUNDIR}/$i ]] ; then
        if [[ $i =~ RESTART/ ]]; then
          cp ${RUNDIR}/$i ${NEW_BASELINE}/${CNTL_DIR}/RESTART/$(basename $i)
        elif [[ $i =~ mediator ]]; then
          cp ${RUNDIR}/$i ${NEW_BASELINE}/${CNTLMED_DIR}
        else
          cp ${RUNDIR}/${i} ${NEW_BASELINE}/${CNTL_DIR}/${i}
        fi
      else
        echo "Missing " ${RUNDIR}/$i " output file"
        echo;echo " Set ${TEST_NR} ${TEST_NAME} failed"
        test_status='FAIL'
      fi
    done

  fi

  echo "Test ${TEST_NR} ${TEST_NAME} ${test_status}" >> ${REGRESSIONTEST_LOG}
  echo                                               >> ${REGRESSIONTEST_LOG}
  echo "Test ${TEST_NR} ${TEST_NAME} ${test_status}"
  echo

  if [[ $test_status = 'FAIL' ]]; then
    if [[ ${UNIT_TEST} == false ]]; then
      echo $TEST_NAME >> $PATHRT/fail_test
    else
      echo ${TEST_NR} $TEST_NAME >> $PATHRT/fail_unit_test
    fi

    if [[ $ROCOTO = true ]]; then
      exit 2
    fi
  fi

  eval "$set_x"
}

kill_job() {

  [[ -z $1 ]] && exit 1

  local -r jobid=$1

  if [[ $SCHEDULER = 'moab' ]]; then
     :
  elif [[ $SCHEDULER = 'pbs' ]]; then
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

  if [[ "Q$APP" != Q ]] ; then
      rocoto_cmd="&PATHRT;/appbuild.sh &PATHTR;/FV3 $APP $COMPILE_NR"
  else
      #rocoto_cmd="&PATHRT;/compile_cmake.sh &PATHTR; $MACHINE_ID \"${NEMS_VER}\" $COMPILE_NR"
      rocoto_cmd="&PATHRT;/compile.sh ${NEMS_VER} $COMPILE_NR"
  fi

  NATIVE=""
  BUILD_CORES=8
  if [[ ${MACHINE_ID} == wcoss_dell_p3 ]]; then
    BUILD_CORES=1
    NATIVE="<memory>8G</memory> <native>-R 'affinity[core(1)]'</native>"
  fi
  if [[ ${MACHINE_ID} == wcoss_cray ]]; then
    BUILD_CORES=24
    rocoto_cmd="aprun -n 1 -j 1 -N 1 -d $BUILD_CORES $rocoto_cmd"
    NATIVE="<exclusive></exclusive>"
  fi
  BUILD_WALLTIME="00:30:00"
  if [[ ${MACHINE_ID} == jet ]]; then
    BUILD_WALLTIME="01:00:00"
  fi
  if [[ ${MACHINE_ID} == orion.* ]]; then
    BUILD_WALLTIME="01:00:00"
  fi

  if [[ ${COMPILE_NR_DEP} -gt 0 ]]; then
  cat << EOF >> $ROCOTO_XML
  <task name="compile_${COMPILE_NR}" maxtries="3">
    <dependency> <taskdep task="compile_${COMPILE_NR_DEP}"/></dependency>
    <command>$rocoto_cmd</command>
    <jobname>compile_${COMPILE_NR}</jobname>
    <account>${ACCNR}</account>
    <queue>${COMPILE_QUEUE}</queue>
    <partition>${PARTITION}</partition>
    <cores>${BUILD_CORES}</cores>
    <walltime>${BUILD_WALLTIME}</walltime>
    <join>&LOG;/compile_${COMPILE_NR}.log</join>
    ${NATIVE}
  </task>
EOF
  else
  cat << EOF >> $ROCOTO_XML
  <task name="compile_${COMPILE_NR}" maxtries="3">
    <command>$rocoto_cmd</command>
    <jobname>compile_${COMPILE_NR}</jobname>
    <account>${ACCNR}</account>
    <queue>${COMPILE_QUEUE}</queue>
    <partition>${PARTITION}</partition>
    <cores>${BUILD_CORES}</cores>
    <walltime>${BUILD_WALLTIME}</walltime>
    <join>&LOG;/compile_${COMPILE_NR}.log</join>
    ${NATIVE}
  </task>
EOF
  fi
}


rocoto_create_run_task() {

  if [[ $CREATE_BASELINE == true && $DEP_RUN != '' ]] || [[ $DEP_RUN != '' ]]; then
    DEP_STRING="<and> <taskdep task=\"compile_${COMPILE_NR}\"/> <taskdep task=\"${DEP_RUN}${RT_SUFFIX}\"/> </and>"
  else
    DEP_STRING="<taskdep task=\"compile_${COMPILE_NR}\"/>"
  fi

  CORES=$(( ${TASKS} * ${THRD} ))
  if (( TPN > CORES )); then
    TPN=$CORES
  fi

  NATIVE=""
  if [[ ${MACHINE_ID} == wcoss ]]; then
    NATIVE="<native>-a poe -R span[ptile=${TPN}]</native>"
  fi
  if [[ ${MACHINE_ID} == wcoss_dell_p3 ]]; then
    NATIVE="<native>-R span[ptile=${TPN}]</native>"
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
      <cores>${CORES}</cores>
      <walltime>00:${WLCLK}:00</walltime>
      <join>&LOG;/run_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log</join>
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
    sleep 10 & wait $!
    state=$($ROCOTOSTAT -w $ROCOTO_XML -d $ROCOTO_DB -s | grep 197001010000 | awk -F" " '{print $2}')
    dead_compile=$($ROCOTOSTAT -w $ROCOTO_XML -d $ROCOTO_DB | grep compile_ | grep DEAD | head -1 | awk -F" " '{print $2}')
    if [[ ! -z ${dead_compile} ]]; then
      echo "y" | ${ROCOTOCOMPLETE} -w $ROCOTO_XML -d $ROCOTO_DB -m ${dead_compile}_tasks
      ${ROCOTOCOMPLETE} -w $ROCOTO_XML -d $ROCOTO_DB -t ${dead_compile}
    fi
    sleep 20 & wait $!
  done

}

ecflow_create_compile_task() {

  new_compile=true

  if [[ "Q$APP" != Q ]] ; then
      ecflow_cmd="$PATHRT/appbuild.sh ${PATHTR}/FV3 $APP $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1"
  else
      #ecflow_cmd="$PATHRT/compile_cmake.sh ${PATHTR} $MACHINE_ID \"${NEMS_VER}\" $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1"
      ecflow_cmd="$PATHRT/compile.sh ${NEMS_VER} $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1"
  fi

  cat << EOF > ${ECFLOW_RUN}/${ECFLOW_SUITE}/compile_${COMPILE_NR}.ecf
%include <head.h>
$ecflow_cmd
%include <tail.h>
EOF

  echo "  task compile_${COMPILE_NR}" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      inlimit max_builds" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  # serialize WW3 builds. FIXME
  if [[ ${NEMS_VER^^} =~ "WW3=Y" && ${COMPILE_PREV_WW3_NR} != '' ]]; then
    echo "    trigger compile_${COMPILE_PREV_WW3_NR} == complete"  >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  fi
}

ecflow_create_run_task() {

  cat << EOF > ${ECFLOW_RUN}/${ECFLOW_SUITE}/${TEST_NAME}${RT_SUFFIX}.ecf
%include <head.h>
$PATHRT/run_test.sh ${PATHRT} ${RUNDIR_ROOT} ${TEST_NAME} ${TEST_NR} ${COMPILE_NR} > ${LOG_DIR}/run_${TEST_NR}_${TEST_NAME}${RT_SUFFIX}.log 2>&1
%include <tail.h>
EOF

  echo "    task ${TEST_NAME}${RT_SUFFIX}" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  echo "      inlimit max_jobs" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  if [[ ${UNIT_TEST} == true && $DEP_RUN != '' ]]; then
    echo "      trigger compile_${COMPILE_NR} == complete and ${DEP_RUN} == complete" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  elif [[ $CREATE_BASELINE == true && $DEP_RUN != '' ]] || [[ $DEP_RUN != '' ]]; then
    echo "      trigger compile_${COMPILE_NR} == complete and ${DEP_RUN}${RT_SUFFIX} == complete" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  else
    echo "      trigger compile_${COMPILE_NR} == complete" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  fi
}

ecflow_run() {

  # in rare instances when UID is greater then 58500 (like Ratko's UID on theia)
  [[ $ECF_PORT -gt 49151 ]] && ECF_PORT=12179

  ECF_HOST=$( hostname )

  sh ${ECFLOW_START} -d ${ECFLOW_RUN} -p ${ECF_PORT} >> ${ECFLOW_RUN}/ecflow.log 2>&1

#  set +e
#  i=0
#  max_atempts=5
#  while [[ $i -lt $max_atempts ]]; do
#    ecflow_client --ping --host=${ECF_HOST} --port=${ECF_PORT}
#    not_running=$?
#    if [[ $not_running -eq 1 ]]; then
#      echo "ecflow_server is NOT running on ${ECF_HOST}:${ECF_PORT}"
#      sh ${ECFLOW_START} -d ${ECFLOW_RUN} -p ${ECF_PORT} >> ${ECFLOW_RUN}/ecflow.log 2>&1
#      break
#    else
#      echo "ecflow_server is already running on ${ECF_HOST}:${ECF_PORT}"
#      ECF_PORT=$(( ECF_PORT + 1 ))
#    fi
#    i=$(( i + 1 ))
#    if [[ $i -eq $max_atempts ]]; then
#       echo "You already have $max_atempts ecFlow servers running on this node"
#       exit 1
#    fi
#  done
#  set -e

  ECFLOW_RUNNING=true

  export ECF_PORT
  export ECF_HOST

  ecflow_client --load=${ECFLOW_RUN}/${ECFLOW_SUITE}.def            >> ${ECFLOW_RUN}/ecflow.log 2>&1
  ecflow_client --begin=${ECFLOW_SUITE}                             >> ${ECFLOW_RUN}/ecflow.log 2>&1

  active_tasks=1
  while [[ $active_tasks -ne 0 ]]
  do
    sleep 10 & wait $!
    active_tasks=$( ecflow_client --get_state /${ECFLOW_SUITE} | grep "task " | grep -E 'state:active|state:submitted|state:queued' | wc -l )
    echo "ecflow tasks remaining: ${active_tasks}"
    ${PATHRT}/abort_dep_tasks.py
  done
  sleep 65 # wait one ECF_INTERVAL plus 5 seconds
}

ecflow_kill() {
   [[ ${ECFLOW_RUNNING:-false} == true ]] || return
   set +e
   wait
   ecflow_client --kill /${ECFLOW_SUITE}
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
