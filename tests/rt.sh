#!/bin/bash
set -eux

SECONDS=0

hostname

die() { echo "$@" >&2; exit 1; }
usage() {
  set +x
  echo
  echo "Usage: $0 -a <account> | -b <file> | -c | -d | -e | -h | -k | -l <file> | -m | -n <name> | -o | -r | -w"
  echo
  echo "  -a  <account> to use on for HPC queue"
  echo "  -b  create new baselines only for tests listed in <file>"
  echo "  -c  create new baseline results"
  echo "  -d  delete run direcotries that are not used by other tests"
  echo "  -e  use ecFlow workflow manager"
  echo "  -h  display this help"
  echo "  -k  keep run directory after rt.sh is completed"
  echo "  -l  runs test specified in <file>"
  echo "  -m  compare against new baseline results"
  echo "  -n  run single test <name>"
  echo "  -o  compile only, skip tests"
  echo "  -r  use Rocoto workflow manager"
  echo "  -w  for weekly_test, skip comparing baseline results"
  echo
  set -x
  exit 1
}

[[ $# -eq 0 ]] && usage

update_rtconf() {

  find_match() {
    # This function finds if a test in $TESTS_FILE matches one 
    # in our list of tests to be run.
    THIS_TEST_WITH_COMPILER=$1
    shift
    TWC=("$@")
    FOUND=false
    for i in "${!TWC[@]}"; do
      if [[ "${TWC[$i]}" == "${THIS_TEST_WITH_COMPILER}" ]]; then
        FOUND=true
        echo "${i}"
        return
      fi
    done
    if [[ $FOUND == false ]]; then
      echo "-1"
    fi
  }

  # This script will update the rt.conf ($TESTS_FILE) if needed by the
  # -b or -n options being called/used.

  # THE USER CHOSE THE -b OPTION
  if [[ $NEW_BASELINES_FILE != '' ]]; then
    [[ -s "$NEW_BASELINES_FILE" ]] || die "${NEW_BASELINES_FILE} is empty, exiting..."
    TEST_WITH_COMPILE=()
    readarray -t TEST_WITH_COMPILE < "$NEW_BASELINES_FILE"
  # else USER CHOSE THE -l OPTION
  elif [[ $DEFINE_CONF_FILE == true ]]; then
    echo "No update needed to TESTS_FILE"
    return
  # else USER CHOSE THE -n OPTION
  elif [[ $RUN_SINGLE_TEST == true ]]; then
    TEST_WITH_COMPILE=("${SRT_NAME} ${SRT_COMPILER}")
  else
    echo "No update needed to rt.conf"
    return
  fi

  RT_TEMP_CONF="rt_temp.conf"
  rm -f $RT_TEMP_CONF && touch $RT_TEMP_CONF
  local compile_line=''
  
  while read -r line || [ "$line" ]; do
    line="${line#"${line%%[![:space:]]*}"}"
    [[ -n $line ]] || continue
    [[ ${#line} == 0 ]] && continue
    [[ $line == \#* ]] && continue
    
    if [[ $line == COMPILE* ]] ; then
      COMPILE_LINE_USED=false
      MACHINES=$(echo $line | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
      RT_COMPILER_IN=$(echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')

      if [[ ${MACHINES} == '' ]]; then
        compile_line=$line
      elif [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || compile_line=$line
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && compile_line=$line
      fi

    fi

    if [[ $line =~ RUN ]]; then
      tmp_test=$(echo "$line" | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
      TEST_IDX=$(find_match "$tmp_test $RT_COMPILER_IN" "${TEST_WITH_COMPILE[@]}")
      
      if [[ $TEST_IDX != -1 ]]; then
        if [[ $COMPILE_LINE_USED == false ]]; then
          echo "$compile_line" >> $RT_TEMP_CONF
          COMPILE_LINE_USED=true
        fi
        dep_test=$(echo "$line" | grep -w "$tmp_test" | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
        
        if [[ $dep_test != '' ]]; then
          if [[ $(find_match "$dep_test $RT_COMPILER_IN" "${TEST_WITH_COMPILE[@]}") == -1 ]]; then

            dep_line=$(grep -w "$dep_test" rt.conf | grep -v "$tmp_test")
            dep_line="${dep_line#"${dep_line%%[![:space:]]*}"}"
            dep_line=$(echo "${dep_line}" | tr -d '\n')
            CORRECT_LINE[1]=$(awk -F'RUN|RUN' '{print $2}' <<< "$dep_line")
            CORRECT_LINE[2]=$(awk -F'RUN|RUN' '{print $3}' <<< "$dep_line")
            
            if [[ $RT_COMPILER_IN == "intel" ]]; then
              echo "RUN ${CORRECT_LINE[1]}" >> $RT_TEMP_CONF
            elif [[ $RT_COMPILER_IN == "gnu" ]]; then
              echo "RUN ${CORRECT_LINE[2]}" >> $RT_TEMP_CONF
            fi
          fi
        fi
        echo "$line" >> $RT_TEMP_CONF
      fi
    fi
  done < "$TESTS_FILE"

  if [[ ! -s $RT_TEMP_CONF ]]; then
    echo "The tests listed/chosen do not exist or cannot be run on $MACHINE_ID"
    exit 1
  else
    TESTS_FILE=$RT_TEMP_CONF
  fi
}

generate_log() {

  COMPILE_COUNTER=0
  FAILED_COMPILES=()
  TEST_COUNTER=0
  FAILED_TESTS=()
  FAILED_TEST_ID=()
  FAILED_COMPILE_LOGS=()
  FAILED_TEST_LOGS=()
  TEST_CHANGES_LOG="test_changes.list"
  TEST_END_TIME="$(date '+%Y%m%d %T')"
  cat << EOF > "${REGRESSIONTEST_LOG}"
====START OF ${MACHINE_ID^^} REGRESSION TESTING LOG====

UFSWM hash used in testing:
$(git rev-parse HEAD)

Submodule hashes used in testing:
EOF
  cd ..
  git submodule status >> "${REGRESSIONTEST_LOG}"
  echo; echo >> "${REGRESSIONTEST_LOG}"
  cd tests
  
  cat << EOF >> "${REGRESSIONTEST_LOG}"

NOTES:
[Times](Memory) are at the end of each compile/test in format [MM:SS](Size).
The first time is for the full script (prep+run+finalize).
The second time is specifically for the run phase.
Times/Memory will be empty for failed tests.

BASELINE DIRECTORY: ${RTPWD}
COMPARISON DIRECTORY: ${RUNDIR_ROOT}

RT.SH OPTIONS USED:
EOF

  [[ -n $ACCNR ]] && echo "* (-a) - HPC PROJECT ACCOUNT: ${ACCNR}" >> "${REGRESSIONTEST_LOG}"
  [[ -n $NEW_BASELINES_FILE ]] && echo "* (-b) - NEW BASELINES FROM FILE: ${NEW_BASELINES_FILE}" >> "${REGRESSIONTEST_LOG}"
  [[ $CREATE_BASELINE == true ]] && echo "* (-c) - CREATE NEW BASELINES" >> "${REGRESSIONTEST_LOG}"
  [[ $DEFINE_CONF_FILE == true ]] && echo "* (-l) - USE CONFIG FILE: ${TESTS_FILE}" >> "${REGRESSIONTEST_LOG}"
  [[ $RTPWD_NEW_BASELINE == true ]] && echo "* (-m) - COMPARE AGAINST CREATED BASELINES" >> "${REGRESSIONTEST_LOG}"
  [[ $RUN_SINGLE_TEST == true ]] && echo "* (-n) - RUN SINGLE TEST: ${SINGLE_OPTS}" >> "${REGRESSIONTEST_LOG}"
  [[ $COMPILE_ONLY == true ]]&& echo "* 9 (-o) COMPILE ONLY, SKIP TESTS" >> "${REGRESSIONTEST_LOG}"
  [[ $delete_rundir == true ]] && echo "* (-d) - DELETE RUN DIRECTORY" >> "${REGRESSIONTEST_LOG}"
  [[ $skip_check_results == true ]] && echo "* (-w) - SKIP RESULTS CHECK" >> "${REGRESSIONTEST_LOG}"
  [[ $KEEP_RUNDIR == true ]] && echo "* (-k) - KEEP RUN DIRECTORY" >> "${REGRESSIONTEST_LOG}"
  [[ $ROCOTO == true ]] && echo "* (-r) - USE ROCOTO" >> "${REGRESSIONTEST_LOG}"
  [[ $ECFLOW == true ]] && echo "* (-e) - USE ECFLOW" >> "${REGRESSIONTEST_LOG}"


  [[ -f "${TEST_CHANGES_LOG}" ]] && rm ${TEST_CHANGES_LOG}
  touch ${TEST_CHANGES_LOG}
  while read -r line || [ "$line" ]; do
    line="${line#"${line%%[![:space:]]*}"}"
    [[ -n "$line" ]] || continue
    [[ ${#line} == 0 ]] && continue
    [[ $line == \#* ]] && continue
    local valid_compile=false
    local valid_test=false

    if [[ $line == COMPILE* ]] ; then
      
      CMACHINES=$(echo "$line" | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
      COMPILER=$(echo "$line" | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
      COMPILE_NAME=$(echo "$line" | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
      COMPILE_ID=${COMPILE_NAME}_${COMPILER}
      
      if [[ ${CMACHINES} == '' ]]; then
        valid_compile=true
      elif [[ ${CMACHINES} == -* ]]; then
        [[ ${CMACHINES} =~ ${MACHINE_ID} ]] || valid_compile=true
      elif [[ ${CMACHINES} == +* ]]; then
        [[ ${CMACHINES} =~ ${MACHINE_ID} ]] && valid_compile=true
      fi

      if [[ $valid_compile == true ]]; then
        COMPILE_COUNTER=$((COMPILE_COUNTER+1))
        FAIL_LOG=""
        COMPILE_RESULT=""
        TIME_FILE=""
        COMPILE_TIME=""
        RT_COMPILE_TIME=""
        if [[ ! -f "${LOG_DIR}/compile_${COMPILE_ID}.log" ]]; then
          COMPILE_RESULT="MISSING"
          FAIL_LOG="N/A"
        elif [[ -f fail_compile_${COMPILE_ID} ]]; then
          COMPILE_RESULT="FAIL TO RUN"
          FAIL_LOG="${LOG_DIR}/compile_${COMPILE_ID}.log"
        else 
          if grep -q "quota" "${LOG_DIR}/compile_${COMPILE_ID}.log"; then
            COMPILE_RESULT="FAIL FROM DISK QUOTA"
            FAIL_LOG="${LOG_DIR}/compile_${COMPILE_ID}.log"
          elif grep -q "timeout" "${LOG_DIR}/compile_${COMPILE_ID}.log"; then
            COMPILE_RESULT="FAIL FROM TIMEOUT"
            FAIL_LOG="${LOG_DIR}/compile_${COMPILE_ID}.log"
          else
            COMPILE_RESULT="PASS"
            TIME_FILE="${LOG_DIR}/compile_${COMPILE_ID}_timestamp.txt"
            if [[ -f "${TIME_FILE}" ]]; then
              while read -r times || [ "$times" ]; do
                  times="${times#"${times%%[![:space:]]*}"}"

                  DATE1=$(echo "$times" | cut -d ',' -f2 )
                  DATE2=$(echo "$times" | cut -d ',' -f3 )
                  DATE3=$(echo "$times" | cut -d ',' -f4 )
                  DATE4=$(echo "$times" | cut -d ',' -f5 )

                  COMPILE_TIME=$(date --date=@$((DATE3 - DATE2)) +'%M:%S')
                  RT_COMPILE_TIME=$(date --date=@$((DATE4 - DATE1)) +'%M:%S')

              done < "$TIME_FILE"
            fi
          fi
        fi
        echo >> "${REGRESSIONTEST_LOG}"
        echo "${COMPILE_RESULT} -- COMPILE '${COMPILE_ID}' [${RT_COMPILE_TIME}, ${COMPILE_TIME}]" >> "${REGRESSIONTEST_LOG}"
        [[ -n $FAIL_LOG ]] && FAILED_COMPILES+=("COMPILE ${COMPILE_ID}: ${COMPILE_RESULT}")
        [[ -n $FAIL_LOG ]] && FAILED_COMPILE_LOGS+=("${FAIL_LOG}")
      fi

    elif [[ $line =~ RUN ]]; then

      if [[ $COMPILE_ONLY == true ]]; then
        continue
      fi

      RMACHINES=$(echo "$line" | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
      TEST_NAME=$(echo "$line" | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
      GEN_BASELINE=$(echo "$line" | cut -d'|' -f4 | sed -e 's/^ *//' -e 's/ *$//')
      
      if [[ ${RMACHINES} == '' ]]; then
        valid_test=true
      elif [[ ${RMACHINES} == -* ]]; then
        [[ ${RMACHINES} =~ ${MACHINE_ID} ]] || valid_test=true
      elif [[ ${RMACHINES} == +* ]]; then
        [[ ${RMACHINES} =~ ${MACHINE_ID} ]] && valid_test=true
      fi

      if [[ $valid_test == true ]]; then
        TEST_COUNTER=$((TEST_COUNTER+1))
        GETMEMFROMLOG=""
        FAIL_LOG=""
        TEST_RESULT=""
        TIME_FILE=""
        TEST_TIME=""
        RT_TEST_TIME=""
        RT_TEST_MEM=""
        if [[ $CREATE_BASELINE == true && $GEN_BASELINE != "baseline" ]]; then
          TEST_RESULT="SKIPPED (TEST DOES NOT GENERATE BASELINE)"
        elif [[ ! -f "${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log" ]]; then
          TEST_RESULT="MISSING"
          FAIL_LOG="N/A"
        elif [[ -f fail_test_${TEST_NAME}_${COMPILER} ]]; then
          if [[ -f "${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log" ]]; then
            TEST_RESULT="FAIL TO COMPARE"
            FAIL_LOG="${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log"
          else
            TEST_RESULT="FAIL TO RUN"
            FAIL_LOG="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"
          fi
        else
          if grep -q "quota" "${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"; then
            TEST_RESULT="FAIL FROM DISK QUOTA"
            FAIL_LOG="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"
          elif grep -q "timeout" "${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"; then
            TEST_RESULT="FAIL FROM TIMEOUT"
            FAIL_LOG="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"
          else
            
            TEST_RESULT="PASS"
            TIME_FILE="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}_timestamp.txt"
            GETMEMFROMLOG=$(grep "The maximum resident set size" "${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log")
            RT_TEST_MEM=$(echo "${GETMEMFROMLOG:9:${#GETMEMFROMLOG}-1}" | tr -dc '0-9')
            RT_TEST_MEM=$((RT_TEST_MEM/1000))
            if [[ -f "${TIME_FILE}" ]]; then
              while read -r times || [ "$times" ]; do
                  times="${times#"${times%%[![:space:]]*}"}"

                  DATE1=$(echo "$times" | cut -d ',' -f2 )
                  DATE2=$(echo "$times" | cut -d ',' -f3 )
                  DATE3=$(echo "$times" | cut -d ',' -f4 )
                  DATE4=$(echo "$times" | cut -d ',' -f5 )

                  TEST_TIME=$(date --date=@$((DATE3 - DATE2)) +'%M:%S')
                  RT_TEST_TIME=$(date --date=@$((DATE4 - DATE1)) +'%M:%S')

              done < "$TIME_FILE"
            fi
          fi
        fi

        echo "${TEST_RESULT} -- TEST '${TEST_NAME}_${COMPILER}' [${RT_TEST_TIME}, ${TEST_TIME}](${RT_TEST_MEM} MB)" >> "${REGRESSIONTEST_LOG}"
        [[ -n $FAIL_LOG ]] && FAILED_TESTS+=("TEST ${TEST_NAME}_${COMPILER}: ${TEST_RESULT}")
        [[ -n $FAIL_LOG ]] && FAILED_TEST_LOGS+=("${FAIL_LOG}")
        [[ -n $FAIL_LOG ]] && FAILED_TEST_ID+=("${TEST_NAME} ${COMPILER}")
      fi
    fi
  done < "$TESTS_FILE"
  
  elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )
  
  cat << EOF >> "${REGRESSIONTEST_LOG}"

SYNOPSIS:
Starting Date/Time: ${TEST_START_TIME}
Ending Date/Time: ${TEST_END_TIME}
Total Time: ${elapsed_time}
Compiles Completed: $((COMPILE_COUNTER-${#FAILED_COMPILES[@]}))/${COMPILE_COUNTER}
Tests Completed: $((TEST_COUNTER-${#FAILED_TESTS[@]}))/${TEST_COUNTER}
EOF
  # PRINT FAILED COMPILES
  if [[ "${#FAILED_COMPILES[@]}" -ne "0" ]]; then
    echo "Failed Compiles:" >> "${REGRESSIONTEST_LOG}"
    for i in "${!FAILED_COMPILES[@]}"; do
      echo "* ${FAILED_COMPILES[$i]}" >> "${REGRESSIONTEST_LOG}"
      echo "-- LOG: ${FAILED_COMPILE_LOGS[$i]}" >> "${REGRESSIONTEST_LOG}"
    done
  fi
  
  # PRINT FAILED TESTS
  if [[ "${#FAILED_TESTS[@]}" -ne "0" ]]; then
    
    echo "Failed Tests:" >> ${REGRESSIONTEST_LOG}
    for j in "${!FAILED_TESTS[@]}"; do
      echo "* ${FAILED_TESTS[$j]}" >> "${REGRESSIONTEST_LOG}"
      echo "-- LOG: ${FAILED_TEST_LOGS[$j]}" >> "${REGRESSIONTEST_LOG}"
    done
    
  fi
  
  # WRITE FAILED_TEST_ID LIST TO TEST_CHANGES_LOG
  if [[ "${#FAILED_TESTS[@]}" -ne "0" ]]; then
    for item in "${FAILED_TEST_ID[@]}"; do
      echo "$item" >> "${TEST_CHANGES_LOG}"
    done
  fi

  if [[ "${#FAILED_COMPILES[@]}" -eq "0" && "${#FAILED_TESTS[@]}" -eq "0" ]]; then
    cat << EOF >> "${REGRESSIONTEST_LOG}"

NOTES:
A file '${TEST_CHANGES_LOG}' was generated but is empty.
If you are using this log as a pull request verification, please commit '${TEST_CHANGES_LOG}'.

Result: SUCCESS

====END OF ${MACHINE_ID^^} REGRESSION TESTING LOG====
EOF
    echo "Performing Cleanup..."
    rm -f fv3_*.x fv3_*.exe modules.fv3_* modulefiles/modules.fv3_* keep_tests.tmp
    [[ ${KEEP_RUNDIR} == false ]] && rm -rf "${RUNDIR_ROOT}" && rm "${PATHRT}/run_dir"
    [[ ${ROCOTO} == true ]] && rm -f "${ROCOTO_XML}" "${ROCOTO_DB}" "${ROCOTO_STATE}" *_lock.db
    [[ ${TEST_35D} == true ]] && rm -f tests/cpld_bmark*_20*
    echo "REGRESSION TEST RESULT: SUCCESS"
  else
    cat << EOF >> "${REGRESSIONTEST_LOG}"

NOTES:
A file '${TEST_CHANGES_LOG}' was generated with list of all failed tests.
You can use './rt.sh -c -b test_changes.list' to create baselines for the failed tests.
If you are using this log as a pull request verification, please commit '${TEST_CHANGES_LOG}'.

Result: FAILURE

====END OF ${MACHINE_ID^^} REGRESSION TESTING LOG====
EOF
    echo "REGRESSION TEST RESULT: FAILURE"
  fi

}

create_or_run_compile_task() {

  cat << EOF > ${RUNDIR_ROOT}/compile_${COMPILE_ID}.env
export JOB_NR=${JOB_NR}
export COMPILE_ID=${COMPILE_ID}
export MACHINE_ID=${MACHINE_ID}
export RT_COMPILER=${RT_COMPILER}
export PATHRT=${PATHRT}
export PATHTR=${PATHTR}
export SCHEDULER=${SCHEDULER}
export ACCNR=${ACCNR}
export QUEUE=${COMPILE_QUEUE}
export PARTITION=${PARTITION}
export ROCOTO=${ROCOTO}
export ECFLOW=${ECFLOW}
export REGRESSIONTEST_LOG=${REGRESSIONTEST_LOG}
export LOG_DIR=${LOG_DIR}
EOF

  if [[ $ROCOTO == true ]]; then
    rocoto_create_compile_task
  elif [[ $ECFLOW == true ]]; then
    ecflow_create_compile_task
  else
    ./run_compile.sh ${PATHRT} ${RUNDIR_ROOT} "${MAKE_OPT}" ${COMPILE_ID} > ${LOG_DIR}/compile_${COMPILE_ID}.log 2>&1
  fi

  RT_SUFFIX=""
  BL_SUFFIX=""
}

rt_35d() {
if [[ $TEST_NAME =~ '35d' ]] ; then
  local sy=$(echo ${DATE_35D} | cut -c 1-4)
  local sm=$(echo ${DATE_35D} | cut -c 5-6)
  local new_test_name="tests/${TEST_NAME}_${DATE_35D}"
  rm -f $new_test_name
  cp tests/$TEST_NAME $new_test_name

  sed -i -e "s/\(export SYEAR\)/\1=\"$sy\"/" $new_test_name
  sed -i -e "s/\(export SMONTH\)/\1=\"$sm\"/" $new_test_name

  TEST_NAME=${new_test_name#tests/}
fi
}

rt_trap() {
  [[ ${ROCOTO:-false} == true ]] && rocoto_kill
  [[ ${ECFLOW:-false} == true ]] && ecflow_kill
  cleanup
}

cleanup() {
  [[ $(awk '{print $2}' < "${LOCKDIR}/PID") == $$ ]] && rm -rf ${LOCKDIR}
  [[ ${ECFLOW:-false} == true ]] && ecflow_stop
  trap 0
  exit
}

trap '{ echo "rt.sh interrupted"; rt_trap ; }' INT
trap '{ echo "rt.sh quit"; rt_trap ; }' QUIT
trap '{ echo "rt.sh terminated"; rt_trap ; }' TERM
trap '{ echo "rt.sh error on line $LINENO"; cleanup ; }' ERR
trap '{ echo "rt.sh finished"; cleanup ; }' EXIT

# PATHRT - Path to regression tests directory
readonly PATHRT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
cd ${PATHRT}

# PATHTR - Path to nmmb trunk directory
readonly PATHTR=$( cd ${PATHRT}/.. && pwd )

# make sure only one instance of rt.sh is running
readonly LOCKDIR="${PATHRT}"/lock
if mkdir "${LOCKDIR}" ; then
  echo $(hostname) $$ > "${LOCKDIR}/PID"
else
  echo "Only one instance of rt.sh can be running at a time"
  exit 1
fi

source detect_machine.sh # Note: this does not set ACCNR. The "if" block below does.
source rt_utils.sh
source module-setup.sh

CREATE_BASELINE=false
ROCOTO=false
ECFLOW=false
KEEP_RUNDIR=false
TEST_35D=false
export skip_check_results=false
export delete_rundir=false
SKIP_ORDER=false
COMPILE_ONLY=false
RTPWD_NEW_BASELINE=false
TESTS_FILE='rt.conf'
NEW_BASELINES_FILE=''
DEFINE_CONF_FILE=false
RUN_SINGLE_TEST=false
ACCNR=${ACCNR:-""}

while getopts ":a:b:cl:mn:dwkreoh" opt; do
  case $opt in
    a)
      ACCNR=$OPTARG
      ;;
    b)
      NEW_BASELINES_FILE=$OPTARG
      ;;
    c)
      CREATE_BASELINE=true
      ;;
    l)
      DEFINE_CONF_FILE=true
      TESTS_FILE=$OPTARG
      grep -q '[^[:space:]]' < "$TESTS_FILE" ||  die "${TESTS_FILE} empty, exiting..."
      ;;
    o)
      COMPILE_ONLY=true
      ;;
    m)
      # redefine RTPWD to point to newly created baseline outputs
      RTPWD_NEW_BASELINE=true
      ;;
    n)
      RUN_SINGLE_TEST=true
      IFS=' ' read -r -a SINGLE_OPTS <<< $OPTARG

      if [[ ${#SINGLE_OPTS[@]} != 2 ]]; then
        die 'The -n option needs <testname> AND <compiler> in quotes, i.e. -n "control_p8 intel"'
      fi

      SRT_NAME="${SINGLE_OPTS[0]}"
      SRT_COMPILER="${SINGLE_OPTS[1]}"

      if [[ "${SRT_COMPILER}" != "intel" ]] && [[ "${SRT_COMPILER}" != "gnu" ]]; then
        die "COMPILER MUST BE 'intel' OR 'gnu'"
      fi
      ;;
    d)
      export delete_rundir=true
      awk -F "|" '{print $5}' rt.conf | grep "\S" > keep_tests.tmp
      ;;
    w)
      export skip_check_results=true
      ;;
    k)
      KEEP_RUNDIR=true
      ;;
    r)
      ROCOTO=true
      ECFLOW=false
      ;;
    e)
      ECFLOW=true
      ROCOTO=false
      ;;
    h)
      usage
      ;;
    \?)
      usage
      die "Invalid option: -$OPTARG"
      ;;
    :)
      usage
      die "Option -$OPTARG requires an argument."
      ;;
  esac
done

#Check to error out if incompatible options are chosen together
[[ $KEEP_RUNDIR == true && $delete_rundir == true ]] && die "-k and -d options cannot be used at the same time"
[[ $ECFLOW == true && $ROCOTO == true ]] && die "-r and -e options cannot be used at the same time"
[[ $CREATE_BASELINE == true && $RTPWD_NEW_BASELINE == true ]] && die "-c and -m options cannot be used at the same time"

if [[ -z "$ACCNR" ]]; then
  echo "Please use -a <account> to set group account to use on HPC"
  exit 1
fi

# Display the machine and account using the format detect_machine.sh used:
echo "Machine: " $MACHINE_ID "    Account: " $ACCNR

if [[ $MACHINE_ID = wcoss2 ]]; then

  module load ecflow/5.6.0.13
  module load intel/19.1.3.304 python/3.8.6
  ECFLOW_START=${ECF_ROOT}/scripts/server_check.sh
  export ECF_OUTPUTDIR=${PATHRT}/ecf_outputdir
  export ECF_COMDIR=${PATHRT}/ecf_comdir
  rm -rf ${ECF_OUTPUTDIR} ${ECF_COMDIR}
  mkdir -p ${ECF_OUTPUTDIR}
  mkdir -p ${ECF_COMDIR}
  export colonifnco=":output"  # hack

  DISKNM=/lfs/h2/emc/nems/noscrub/emc.nems/RT
  QUEUE=dev
  COMPILE_QUEUE=dev
  ROCOTO_SCHEDULER=pbs
  PARTITION=
  STMP=/lfs/h2/emc/ptmp
  PTMP=/lfs/h2/emc/ptmp
  SCHEDULER=pbs

elif [[ $MACHINE_ID = acorn ]]; then

  module load ecflow/5.6.0.13
  module load intel/19.1.3.304 python/3.8.6
  ECFLOW_START=${ECF_ROOT}/scripts/server_check.sh
  export ECF_OUTPUTDIR=${PATHRT}/ecf_outputdir
  export ECF_COMDIR=${PATHRT}/ecf_comdir
  rm -rf ${ECF_OUTPUTDIR} ${ECF_COMDIR}
  mkdir -p ${ECF_OUTPUTDIR}
  mkdir -p ${ECF_COMDIR}
  export colonifnco=":output"  # hack

  DISKNM=/lfs/h1/emc/nems/noscrub/emc.nems/RT
  QUEUE=dev
  COMPILE_QUEUE=dev
  ROCOTO_SCHEDULER=pbs
  PARTITION=
  STMP=/lfs/h2/emc/ptmp
  PTMP=/lfs/h2/emc/ptmp
  SCHEDULER=pbs

elif [[ $MACHINE_ID = gaea ]]; then

  module use /ncrc/proj/epic/rocoto/modulefiles
  module load rocoto
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  
  module load PrgEnv-intel/8.3.3
  module load intel-classic/2023.1.0
  module load cray-mpich/8.1.25
  module load python/3.9.12
  module use /ncrc/proj/epic/spack-stack/modulefiles
  module load gcc/12.2.0
  module load ecflow/5.8.4
  ECFLOW_START=/ncrc/proj/epic/spack-stack/ecflow-5.8.4/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  DISKNM=/gpfs/f5/epic/world-shared/UFS-WM_RT
  QUEUE=normal
  COMPILE_QUEUE=normal
  PARTITION=c5
  STMP=/gpfs/f5/epic/scratch
  PTMP=/gpfs/f5/epic/scratch

  SCHEDULER=slurm

elif [[ $MACHINE_ID = hera ]]; then

  module load rocoto
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  module load ecflow/5.5.3
  ECFLOW_START=ecflow_start.sh

  QUEUE=batch
  COMPILE_QUEUE=batch

  PARTITION=
  dprefix=/scratch1/NCEPDEV
  DISKNM=/scratch2/NAGAPE/epic/UFS-WM_RT
  STMP=$dprefix/stmp4
  PTMP=$dprefix/stmp2

  SCHEDULER=slurm

elif [[ $MACHINE_ID = orion ]]; then

  module load git/2.28.0
  module load gcc/10.2.0
  module load python/3.9.2

  module load contrib rocoto
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  module use /work/noaa/epic/role-epic/spack-stack/orion/modulefiles
  module load ecflow/5.8.4
  ECFLOW_START=/work/noaa/epic/role-epic/spack-stack/orion/ecflow-5.8.4/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  QUEUE=batch
  COMPILE_QUEUE=batch
  PARTITION=orion
  dprefix=/work/noaa/stmp/${USER}
  DISKNM=/work/noaa/epic/UFS-WM_RT
  STMP=$dprefix/stmp
  PTMP=$dprefix/stmp

  SCHEDULER=slurm

elif [[ $MACHINE_ID = hercules ]]; then

  module load contrib rocoto
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  module use /work/noaa/epic/role-epic/spack-stack/hercules/modulefiles
  module load ecflow/5.8.4
  ECFLOW_START=/work/noaa/epic/role-epic/spack-stack/hercules/ecflow-5.8.4/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  QUEUE=windfall
  COMPILE_QUEUE=windfall
  PARTITION=hercules
  dprefix=/work2/noaa/stmp/${USER}
  DISKNM=/work/noaa/epic/hercules/UFS-WM_RT
  STMP=$dprefix/stmp
  PTMP=$dprefix/stmp

  SCHEDULER=slurm
  cp fv3_conf/fv3_slurm.IN_hercules fv3_conf/fv3_slurm.IN
  cp fv3_conf/compile_slurm.IN_hercules fv3_conf/compile_slurm.IN

elif [[ $MACHINE_ID = jet ]]; then

  module load rocoto
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  module load ecflow/5.5.3
  ECFLOW_START=/apps/ecflow/5.5.3/bin/ecflow_start.sh

  QUEUE=batch
  COMPILE_QUEUE=batch
  PARTITION=xjet
  DISKNM=/mnt/lfs4/HFIP/hfv3gfs/role.epic/RT
  dprefix=${dprefix:-/lfs4/HFIP/$ACCNR/$USER}
  STMP=${STMP:-$dprefix/RT_BASELINE}
  PTMP=${PTMP:-$dprefix/RT_RUNDIRS}

  SCHEDULER=slurm

elif [[ $MACHINE_ID = s4 ]]; then

  module load rocoto/1.3.2
  module load ecflow/5.6.0
  module load miniconda/3.8-s4
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  module use /data/prod/jedi/spack-stack/modulefiles
  module load ecflow/5.8.4
  ECFLOW_START=/data/prod/jedi/spack-stack/ecflow-5.8.4/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  QUEUE=s4
  COMPILE_QUEUE=s4

  PARTITION=s4
  dprefix=/data/prod
  DISKNM=$dprefix/emc.nemspara/RT
  STMP=/scratch/short/users
  PTMP=/scratch/users

  SCHEDULER=slurm

elif [[ $MACHINE_ID = derecho ]]; then

  module use /glade/work/epicufsrt/contrib/derecho/rocoto/modulefiles
  module load rocoto
  module use /glade/work/epicufsrt/contrib/spack-stack/derecho/modulefiles
  module load ecflow/5.8.4
  module unload ncarcompilers
  module use /glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core
  module load stack-intel/2021.10.0
  module load stack-python/3.10.8
#  export PYTHONPATH=/glade/p/ral/jntp/tools/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/glade/p/ral/jntp/tools/miniconda3/4.8.3/lib/python3.8/site-packages
  ECFLOW_START=/glade/work/epicufsrt/contrib/spack-stack/derecho/ecflow-5.8.4/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  QUEUE=main
  COMPILE_QUEUE=main
  PARTITION=
  dprefix=/glade/derecho/scratch
  DISKNM=/glade/derecho/scratch/epicufsrt/ufs-weather-model/RT/
  STMP=$dprefix
  PTMP=$dprefix
  SCHEDULER=pbs
  cp fv3_conf/fv3_qsub.IN_derecho fv3_conf/fv3_qsub.IN
  cp fv3_conf/compile_qsub.IN_derecho fv3_conf/compile_qsub.IN

  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=pbspro

elif [[ $MACHINE_ID = stampede ]]; then

  export PYTHONPATH=
  ECFLOW_START=
  QUEUE=skx-normal
  COMPILE_QUEUE=skx-dev
  PARTITION=
  dprefix=$SCRATCH/ufs-weather-model/run
  DISKNM=/work2/07736/minsukji/stampede2/ufs-weather-model/RT
  STMP=$dprefix
  PTMP=$dprefix
  SCHEDULER=slurm
  MPIEXEC=ibrun
  MPIEXECOPTS=

elif [[ $MACHINE_ID = expanse ]]; then

  export PYTHONPATH=
  ECFLOW_START=
  QUEUE=compute
  COMPILE_QUEUE=shared
  PARTITION=
  dprefix=/expanse/lustre/scratch/$USER/temp_project/run
  DISKNM=/expanse/lustre/scratch/domh/temp_project/RT
  STMP=$dprefix
  PTMP=$dprefix
  SCHEDULER=slurm

 elif [[ $MACHINE_ID = noaacloud ]]; then

  export PATH=/contrib/EPIC/bin:$PATH
  module use /apps/modules/modulefiles
  module load rocoto/1.3.3

  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  QUEUE=batch
  COMPILE_QUEUE=batch
  PARTITION=
  dprefix=/lustre/
  DISKNM=/contrib/ufs-weather-model/RT
  STMP=$dprefix/stmp4
  PTMP=$dprefix/stmp2
  SCHEDULER=slurm


else
  die "Unknown machine ID, please edit detect_machine.sh file"
fi

mkdir -p ${STMP}/${USER}

NEW_BASELINE=${STMP}/${USER}/FV3_RT/REGRESSION_TEST

# Overwrite default RUNDIR_ROOT if environment variable RUNDIR_ROOT is set
RUNDIR_ROOT=${RUNDIR_ROOT:-${PTMP}/${USER}/FV3_RT}/rt_$$
mkdir -p ${RUNDIR_ROOT}
if [[ -e ${PATHRT}/run_dir ]]; then
  rm ${PATHRT}/run_dir
fi
echo "Linking ${RUNDIR_ROOT} to ${PATHRT}/run_dir"
ln -s ${RUNDIR_ROOT} ${PATHRT}/run_dir
echo "Run regression test in: ${RUNDIR_ROOT}"

update_rtconf

if [[ $TESTS_FILE =~ '35d' ]] || [[ $TESTS_FILE =~ 'weekly' ]]; then
  TEST_35D=true
fi

source bl_date.conf

if [[ "$RTPWD_NEW_BASELINE" == true ]] ; then
  RTPWD=${NEW_BASELINE}
else
  RTPWD=${RTPWD:-$DISKNM/NEMSfv3gfs/develop-${BL_DATE}}
fi

if [[ "$CREATE_BASELINE" == false ]] ; then
  if [[ ! -d "$RTPWD" ]] ; then
    echo "Baseline directory does not exist:"
    echo "   $RTPWD"
    exit 1
  elif [[ $( ls -1 "$RTPWD/" | wc -l ) -lt 1 ]] ; then
    echo "Baseline directory is empty:"
    echo "   $RTPWD"
    exit 1
  fi
fi

INPUTDATA_ROOT=${INPUTDATA_ROOT:-$DISKNM/NEMSfv3gfs/input-data-20221101}
INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT}/WW3_input_data_20220624
INPUTDATA_ROOT_BMIC=${INPUTDATA_ROOT_BMIC:-$DISKNM/NEMSfv3gfs/BM_IC-20220207}

shift $((OPTIND-1))
[[ $# -gt 1 ]] && usage

if [[ $CREATE_BASELINE == true ]]; then
  # PREPARE NEW REGRESSION TEST DIRECTORY
  rm -rf "${NEW_BASELINE}"
  mkdir -p "${NEW_BASELINE}"

fi

if [[ $skip_check_results == true ]]; then
  REGRESSIONTEST_LOG=${PATHRT}/logs/RegressionTests_weekly_$MACHINE_ID.log
else
  REGRESSIONTEST_LOG=${PATHRT}/logs/RegressionTests_$MACHINE_ID.log
fi

export TEST_START_TIME="$(date '+%Y%m%d %T')"

source default_vars.sh

JOB_NR=0
COMPILE_COUNTER=0
rm -f fail_test* fail_compile*

export LOG_DIR=${PATHRT}/logs/log_$MACHINE_ID
rm -rf ${LOG_DIR}
mkdir -p ${LOG_DIR}

if [[ $ROCOTO == true ]]; then

  ROCOTO_XML=${PATHRT}/rocoto_workflow.xml
  ROCOTO_STATE=${PATHRT}/rocoto_workflow.state
  ROCOTO_DB=${PATHRT}/rocoto_workflow.db

  rm -f $ROCOTO_XML $ROCOTO_DB $ROCOTO_STATE *_lock.db

  if [[ $MACHINE_ID = stampede || $MACHINE_ID = expanse ]]; then
    die "Rocoto is not supported on this machine: $MACHINE_ID"
  fi

  cat << EOF > $ROCOTO_XML
<?xml version="1.0"?>
<!DOCTYPE workflow
[
  <!ENTITY PATHRT         "${PATHRT}">
  <!ENTITY LOG            "${LOG_DIR}">
  <!ENTITY PATHTR         "${PATHTR}">
  <!ENTITY RTPWD          "${RTPWD}">
  <!ENTITY INPUTDATA_ROOT "${INPUTDATA_ROOT}">
  <!ENTITY INPUTDATA_ROOT_WW3 "${INPUTDATA_ROOT_WW3}">
  <!ENTITY INPUTDATA_ROOT_BMIC "${INPUTDATA_ROOT_BMIC}">
  <!ENTITY RUNDIR_ROOT    "${RUNDIR_ROOT}">
  <!ENTITY NEW_BASELINE   "${NEW_BASELINE}">
]>
<workflow realtime="F" scheduler="${ROCOTO_SCHEDULER}" taskthrottle="20">
  <cycledef>197001010000 197001010000 01:00:00</cycledef>
  <log>&LOG;/workflow.log</log>
EOF

fi

if [[ $ECFLOW == true ]]; then

  # Default maximum number of compile and run jobs
  MAX_BUILDS=10
  MAX_JOBS=30

  # Default number of tries to run jobs
  ECF_TRIES=2

  # Reduce maximum number of compile jobs on jet and s4 because of licensing issues
  if [[ $MACHINE_ID = jet ]]; then
    MAX_BUILDS=5
  elif [[ $MACHINE_ID = s4 ]]; then
    MAX_BUILDS=1
  fi

  ECFLOW_RUN=${PATHRT}/ecflow_run
  ECFLOW_SUITE=regtest_$$
  rm -rf ${ECFLOW_RUN}
  mkdir -p ${ECFLOW_RUN}/${ECFLOW_SUITE}
  cp head.h tail.h ${ECFLOW_RUN}
  cat << EOF > ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
suite ${ECFLOW_SUITE}
    edit ECF_HOME '${ECFLOW_RUN}'
    edit ECF_INCLUDE '${ECFLOW_RUN}'
    edit ECF_KILL_CMD kill -15 %ECF_RID% > %ECF_JOB%.kill 2>&1
    edit ECF_TRIES ${ECF_TRIES}
    label src_dir '${PATHTR}'
    label run_dir '${RUNDIR_ROOT}'
    limit max_builds ${MAX_BUILDS}
    limit max_jobs ${MAX_JOBS}
EOF

  if [[ $MACHINE_ID = stampede || $MACHINE_ID = expanse ]]; then
    die "ecFlow is not supported on this machine: $MACHINE_ID"
  fi

fi

##
## read rt.conf and then either execute the test script directly or create
## workflow description file
##

new_compile=false
in_metatask=false

[[ -f $TESTS_FILE ]] || die "$TESTS_FILE does not exist"

LAST_COMPILER_NR=-9999
COMPILE_PREV=''

declare -A compiles

while read -r line || [ "$line" ]; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ $line == \#* ]] && continue

  JOB_NR=$( printf '%03d' $(( 10#$JOB_NR + 1 )) )

  if [[ $line == COMPILE* ]]; then

    COMPILE_NAME=$( echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
    RT_COMPILER=$(echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
    MAKE_OPT=$(   echo $line | cut -d'|' -f4 | sed -e 's/^ *//' -e 's/ *$//')
    MACHINES=$(   echo $line | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
    CB=$(         echo $line | cut -d'|' -f6)
    COMPILE_ID=${COMPILE_NAME}_${RT_COMPILER}
    COMPILE_PREV=${COMPILE_ID}

    set +u
    if [[ ! -z ${compiles[$COMPILE_ID]} ]] ; then
        echo "Error! Duplicated compilation $COMPILE_NAME for compiler $RT_COMPILER!"
        exit 1
    fi
    set -u
    compiles[$COMPILE_ID]=$COMPILE_ID
    echo "COMPILING ${compiles[${COMPILE_ID}]}"

    [[ $CREATE_BASELINE == true && $CB != *fv3* ]] && continue

    if [[ ${MACHINES} != '' ]]; then
      if [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && continue
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || continue
      else
        echo "MACHINES=|${MACHINES}|"
        die "MACHINES spec must be either an empty string or start with either '+' or '-'"
      fi
    fi

    create_or_run_compile_task

    continue

  elif [[ $line == RUN* ]] ; then

    if [[ $COMPILE_ONLY == true ]]; then
      continue
    fi

    TEST_NAME=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
    MACHINES=$( echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
    CB=$(       echo $line | cut -d'|' -f4)
    DEP_RUN=$(  echo $line | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
    DATE_35D=$( echo $line | cut -d'|' -f6 | sed -e 's/^ *//' -e 's/ *$//')

    if [[ $DEP_RUN != '' ]]; then
      DEP_RUN=${DEP_RUN}_${RT_COMPILER}
    fi

    export TEST_ID=${TEST_NAME}_${RT_COMPILER}

    [[ -e "tests/$TEST_NAME" ]] || die "run test file tests/$TEST_NAME does not exist"
    [[ $CREATE_BASELINE == true && $CB != *baseline* ]] && continue

    if [[ ${MACHINES} != '' ]]; then
      if [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && continue
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || continue
      else
        echo "MACHINES=|${MACHINES}|"
        die "MACHINES spec must be either an empty string or start with either '+' or '-'"
      fi
    fi

    COMPILE_METATASK_NAME=${COMPILE_ID}

    # 35 day tests
    [[ $TEST_35D == true ]] && rt_35d

    # Avoid uninitialized RT_SUFFIX/BL_SUFFIX (see definition above)
    RT_SUFFIX=${RT_SUFFIX:-""}
    BL_SUFFIX=${BL_SUFFIX:-""}

    if [[ $ROCOTO == true && $new_compile == true ]]; then
      new_compile=false
      in_metatask=true
      cat << EOF >> $ROCOTO_XML
  <metatask name="compile_${COMPILE_METATASK_NAME}_tasks"><var name="zero">0</var>
EOF
    fi

    (
      source ${PATHRT}/tests/$TEST_NAME

      compute_petbounds_and_tasks

      TPN=$(( TPN / THRD ))
      NODES=$(( TASKS / TPN ))
      if (( NODES * TPN < TASKS )); then
        NODES=$(( NODES + 1 ))
      fi

      PPN=$(( TASKS / NODES ))
      if (( TASKS - ( PPN * NODES ) > 0 )); then
          PPN=$((PPN + 1))
      fi

      cat << EOF > ${RUNDIR_ROOT}/run_test_${TEST_ID}.env
export JOB_NR=${JOB_NR}
export TEST_ID=${TEST_ID}
export MACHINE_ID=${MACHINE_ID}
export RT_COMPILER=${RT_COMPILER}
export RTPWD=${RTPWD}
export INPUTDATA_ROOT=${INPUTDATA_ROOT}
export INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT_WW3}
export INPUTDATA_ROOT_BMIC=${INPUTDATA_ROOT_BMIC}
export PATHRT=${PATHRT}
export PATHTR=${PATHTR}
export NEW_BASELINE=${NEW_BASELINE}
export CREATE_BASELINE=${CREATE_BASELINE}
export RT_SUFFIX=${RT_SUFFIX}
export BL_SUFFIX=${BL_SUFFIX}
export SCHEDULER=${SCHEDULER}
export ACCNR=${ACCNR}
export QUEUE=${QUEUE}
export PARTITION=${PARTITION}
export ROCOTO=${ROCOTO}
export ECFLOW=${ECFLOW}
export REGRESSIONTEST_LOG=${REGRESSIONTEST_LOG}
export LOG_DIR=${LOG_DIR}
export DEP_RUN=${DEP_RUN}
export skip_check_results=${skip_check_results}
export delete_rundir=${delete_rundir}
export WLCLK=${WLCLK}
EOF
      if [[ $MACHINE_ID = jet ]]; then
        cat << EOF >> ${RUNDIR_ROOT}/run_test_${TEST_ID}.env
export PATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/bin:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/bin:$PATH
export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/lib/python3.8/site-packages
EOF
      fi

      if [[ $ROCOTO == true ]]; then
        rocoto_create_run_task
      elif [[ $ECFLOW == true ]]; then
        ecflow_create_run_task
      else
        ./run_test.sh ${PATHRT} ${RUNDIR_ROOT} ${TEST_NAME} ${TEST_ID} ${COMPILE_ID} > ${LOG_DIR}/run_${TEST_ID}${RT_SUFFIX}.log 2>&1
      fi
    )

    continue
  else
    die "Unknown command $line"
  fi
done < $TESTS_FILE

##
## run regression test workflow (currently Rocoto or ecFlow are supported)
##

if [[ $ROCOTO == true ]]; then
  if [[ $in_metatask == true ]]; then
    echo "  </metatask>" >> $ROCOTO_XML
  fi
  echo "</workflow>" >> $ROCOTO_XML
  # run rocoto workflow until done
  rocoto_run
fi

if [[ $ECFLOW == true ]]; then
  echo "endsuite" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  # run ecflow workflow until done
  ecflow_run
fi

# IF -c AND -b; LINK VERIFIED BASELINES TO NEW_BASELINE
if [[ $CREATE_BASELINE == true && $NEW_BASELINES_FILE != '' ]]; then
  for dir in "${RTPWD}"/*/; do
    dir=${dir%*/}
    [[ -d "${NEW_BASELINE}/${dir##*/}" ]] && continue
    ln -s "${dir%*/}" "${NEW_BASELINE}/"
  done
fi

## Lets verify all tests were run and that they passed
generate_log
