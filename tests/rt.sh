#!/bin/bash
set -eu
set -o errexit #Lets trap exit info as error for logging
echo "******Regression Testing Script Started******"
SECONDS=0

hostname

die() { echo "$@" >&2; exit 1; }

usage() {
  set +x #No reason to print out a bunch of echo statements here
  echo
  echo "Usage: $0 -a <account> | -b <file> | -c | -d | -e | -h | -k | -l <file> | -m | -n <name> | -o | -p | -P <file> | -r | -v | -w | -x"
  echo
  echo "  -a  <account> to use on for HPC queue"
  echo "  -b  create new baselines only for tests listed in <file>"
  echo "  -c  create new baseline results"
  echo "  -d  delete run directories that are not used by other tests"
  echo "  -e  use ecFlow workflow manager"
  echo "  -h  display this help"
  echo "  -k  keep run directory after rt.sh is completed"
  echo "  -l  runs test specified in <file>"
  echo "  -m  compare against new baseline results"
  echo "  -n  run single test <name>"
  echo "  -o  compile only, skip tests"
  echo "  -p  build and run inside the GNU/Intel container staged on this Tier 1"
  echo "      platform (compiler taken from each COMPILE line); container"
  echo "      baselines and logs are kept separate from the native-stack ones"
  echo "  -P  <file> build and run on a community (non-Tier-1) platform defined"
  echo "      in <file>; a portability check, not a regression test -- always"
  echo "      sequential (no -r/-e) and always skips baseline comparison"
  echo "  -r  use Rocoto workflow manager"
  echo "  -v  verbose output"
  echo "  -w  for weekly_test, skip comparing baseline results"
  echo "  -x  dry-run"
  echo
}

[[ $# -eq 0 ]] && usage

update_rtconf() {
  echo "rt.sh: Checking & Updating test configuration..."
  find_match() {
    # This function finds if a test in $TESTS_FILE matches one
    # in our list of tests to be run.
    THIS_TEST_WITH_COMPILER=$1
    shift
    TWC=("$@")
    FOUND=false
    for i in "${!TWC[@]}"; do
      if [[ "${TWC[${i}]}" == "${THIS_TEST_WITH_COMPILER}" ]]; then
        FOUND=true
        echo "${i}"
        return
      fi
    done
    if [[ ${FOUND} == false ]]; then
      echo "-1"
    fi
  }

  # This script will update the rt.conf ($TESTS_FILE) if needed by the
  # -b or -n options being called/used.

  # THE USER CHOSE THE -b OPTION
  if [[ ${NEW_BASELINES_FILE} != '' ]]; then
    [[ -s "${NEW_BASELINES_FILE}" ]] || die "${NEW_BASELINES_FILE} is empty, exiting..."
    TEST_WITH_COMPILE=()
    readarray -t TEST_WITH_COMPILE < "${NEW_BASELINES_FILE}"
  # else USER CHOSE THE -n OPTION
  elif [[ ${RUN_SINGLE_TEST} == true ]]; then
    TEST_WITH_COMPILE=("${SRT_NAME} ${SRT_COMPILER}")
  else
    echo "No update needed to rt.conf"
    return
  fi

  RT_TEMP_CONF="rt_temp.conf"
  rm -f "${RT_TEMP_CONF}" && touch "${RT_TEMP_CONF}"
  local compile_line=''
  while read -r line || [[ -n "${line}" ]]; do
    line="${line#"${line%%[![:space:]]*}"}"
    [[ -n "${line}" ]] || continue
    [[ ${#line} == 0 ]] && continue
    [[ ${line} == \#* ]] && continue

    if [[ ${line} =~ COMPILE ]] ; then
      MACHINES=$(cut -d'|' -f5 <<< "${line}")
      MACHINES=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${MACHINES}")
      RT_COMPILER_IN=$(cut -d'|' -f3 <<< "${line}")
      RT_COMPILER_IN=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${RT_COMPILER_IN}")
      machines_allow_run "${MACHINES}" && compile_line=${line}
      COMPILE_LINE_USED=false

    fi

    if [[ ${line} =~ RUN ]]; then
      to_run_test=false
      tmp_test=$(cut -d'|' -f2 <<< "${line}")
      tmp_test=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${tmp_test}")
      MACHINES=$(cut -d'|' -f3 <<< "${line}")
      MACHINES=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${MACHINES}")
      machines_allow_run "${MACHINES}" && to_run_test=true
      if [[ ${to_run_test} == true ]]; then
        TEST_IDX=$(set -e; find_match "${tmp_test} ${RT_COMPILER_IN}" "${TEST_WITH_COMPILE[@]}")

        if [[ ${TEST_IDX} != -1 ]]; then
          if [[ ${COMPILE_LINE_USED} == false ]]; then
            echo -en '\n' >> "${RT_TEMP_CONF}"
            echo "${compile_line}" >> "${RT_TEMP_CONF}"

            COMPILE_LINE_USED=true
          fi
          dep_test=$(grep -w "${tmp_test}" <<< "${line}")
          dep_test=$(cut -d'|' -f5 <<< "${dep_test}")
          dep_test=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${dep_test}")

          if [[ ${dep_test} != '' ]]; then
            find_match_result=$(set -e; find_match "${dep_test} ${RT_COMPILER_IN}" "${TEST_WITH_COMPILE[@]}")
            if [[ ${find_match_result} == -1 ]]; then
              dep_line=$(grep -w "${dep_test}" rt.conf)
              dep_line=$(grep -v "${tmp_test}" <<< "${dep_line}")
              dep_line="${dep_line#"${dep_line%%[![:space:]]*}"}"
              dep_line=$(tr -d '\n' <<< "${dep_line}")
              CORRECT_LINE[1]=$(awk -F'RUN|RUN' '{print $2}' <<< "${dep_line}")
              CORRECT_LINE[2]=$(awk -F'RUN|RUN' '{print $3}' <<< "${dep_line}")

              if [[ ${RT_COMPILER_IN} == "intel" ]]; then
                echo "RUN ${CORRECT_LINE[1]}" >> "${RT_TEMP_CONF}"
              elif [[ ${RT_COMPILER_IN} == "gnu" ]]; then
                echo "RUN ${CORRECT_LINE[2]}" >> "${RT_TEMP_CONF}"
              fi
            fi
          fi
          echo "${line}" >> "${RT_TEMP_CONF}"
        fi
      fi
    fi
  done < "${TESTS_FILE}"

  if [[ ! -s ${RT_TEMP_CONF} ]]; then
    echo "The tests listed/chosen do not exist or cannot be run on ${MACHINE_ID}"
    exit 1
  else
    TESTS_FILE=${RT_TEMP_CONF}
  fi
}

generate_log() {
  echo "rt.sh: Generating Regression Testing Log..."
  COMPILE_COUNTER=0
  FAILED_COMPILES=()
  TEST_COUNTER=0
  FAILED_TESTS=()
  SKIPPED_TESTS=()
  FAILED_TEST_ID=()
  FAILED_COMPILE_LOGS=()
  FAILED_TEST_LOGS=()
  TEST_CHANGES_LOG="test_changes.list"
  TEST_END_TIME="$(date '+%Y%m%d %T')"
  GIT_HASHES=$(git rev-parse HEAD)
  cat << EOF > "${REGRESSIONTEST_LOG}"
====START OF ${MACHINE_ID^^} REGRESSION TESTING LOG====

UFSWM hash used in testing:
${GIT_HASHES}

Submodule hashes used in testing:
EOF
  cd ..
  if  [[ ${MACHINE_ID} != hera  ]]; then
    git submodule status --recursive >> "${REGRESSIONTEST_LOG}"
  else
    git submodule status >> "${REGRESSIONTEST_LOG}"
  fi
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

  [[ -n ${ACCNR} ]] && echo "* (-a) - HPC PROJECT ACCOUNT: ${ACCNR}" >> "${REGRESSIONTEST_LOG}"
  [[ -n ${NEW_BASELINES_FILE} ]] && echo "* (-b) - NEW BASELINES FROM FILE: ${NEW_BASELINES_FILE}" >> "${REGRESSIONTEST_LOG}"
  [[ ${CREATE_BASELINE} == true ]] && echo "* (-c) - CREATE NEW BASELINES" >> "${REGRESSIONTEST_LOG}"
  [[ ${DEFINE_CONF_FILE} == true ]] && echo "* (-l) - USE CONFIG FILE: ${TESTS_FILE}" >> "${REGRESSIONTEST_LOG}"
  [[ ${RTPWD_NEW_BASELINE} == true ]] && echo "* (-m) - COMPARE AGAINST CREATED BASELINES" >> "${REGRESSIONTEST_LOG}"
  [[ ${RUN_SINGLE_TEST} == true ]] && echo "* (-n) - RUN SINGLE TEST: ${SINGLE_OPTS}" >> "${REGRESSIONTEST_LOG}"
  [[ ${COMPILE_ONLY} == true ]]&& echo "* (-o) - COMPILE ONLY, SKIP TESTS" >> "${REGRESSIONTEST_LOG}"
  [[ ${CONTAINER_USE} == true ]] && echo "* (-p) - CONTAINER MODE ON ${MACHINE_ID} (baselines: ${RTPWD})" >> "${REGRESSIONTEST_LOG}"
  [[ ${COMMUNITY_PLATFORM_USE} == true ]] && echo "* (-P) - COMMUNITY PLATFORM: ${MACHINE_ID} (${COMMUNITY_PLATFORM_FILE}); no baseline comparison" >> "${REGRESSIONTEST_LOG}"
  [[ ${delete_rundir} == true ]] && echo "* (-d) - DELETE RUN DIRECTORY" >> "${REGRESSIONTEST_LOG}"
  [[ ${skip_check_results} == true ]] && echo "* (-w) - SKIP RESULTS CHECK" >> "${REGRESSIONTEST_LOG}"
  [[ ${KEEP_RUNDIR} == true ]] && echo "* (-k) - KEEP RUN DIRECTORY" >> "${REGRESSIONTEST_LOG}"
  [[ ${ROCOTO} == true ]] && echo "* (-r) - USE ROCOTO" >> "${REGRESSIONTEST_LOG}"
  [[ ${ECFLOW} == true ]] && echo "* (-e) - USE ECFLOW" >> "${REGRESSIONTEST_LOG}"
  [[ ${RTVERBOSE} == true ]] && echo "* (-v) - VERBOSE OUTPUT" >> "${REGRESSIONTEST_LOG}"


  [[ -f "${TEST_CHANGES_LOG}" ]] && rm "${TEST_CHANGES_LOG}"
  touch "${TEST_CHANGES_LOG}"
  while read -r line; do
    line="${line#"${line%%[![:space:]]*}"}"
    [[ -n "${line}" ]] || continue
    [[ ${#line} == 0 ]] && continue
    [[ ${line} == \#* ]] && continue
    local valid_compile=false
    local valid_test=false

    if [[ ${line} == COMPILE* ]] ; then

      CMACHINES=$(cut -d'|' -f5 <<< "${line}")
      CMACHINES=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${CMACHINES}")

      COMPILER=$(cut -d'|' -f3 <<< "${line}")
      COMPILER=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${COMPILER}")

      COMPILE_NAME=$(cut -d'|' -f2 <<< "${line}")
      COMPILE_NAME=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${COMPILE_NAME}")

      COMPILE_ID=${COMPILE_NAME}_${COMPILER}

      machines_allow_run "${CMACHINES}" && valid_compile=true

      if [[ (${CONTAINER_USE} == true || ${COMMUNITY_PLATFORM_USE} == true) && ${valid_compile} == true ]]; then
        RT_COMPILER=${COMPILER}
        resolve_container_image || valid_compile=false
      fi

      if [[ ${valid_compile} == true ]]; then
        COMPILE_COUNTER=$((COMPILE_COUNTER+1))
        FAIL_LOG=""
        COMPILE_RESULT=""
        TIME_FILE=""
        COMPILE_TIME=""
        RT_COMPILE_TIME=""
        COMPILE_WARNINGS=""
        if [[ ! -f "${LOG_DIR}/compile_${COMPILE_ID}.log" ]]; then
          COMPILE_RESULT="FAILED: UNABLE TO START COMPILE"
          FAIL_LOG="N/A"
        elif [[ -f fail_compile_${COMPILE_ID} ]]; then
          COMPILE_RESULT="FAILED: UNABLE TO FINISH COMPILE"
          FAIL_LOG="${LOG_DIR}/compile_${COMPILE_ID}.log"
          if grep -q "quota" "${LOG_DIR}/compile_${COMPILE_ID}.log"; then
            COMPILE_RESULT="FAILED: DISK QUOTA ISSUE"
            FAIL_LOG="${LOG_DIR}/compile_${COMPILE_ID}.log"
          elif grep -q "TIME LIMIT" "${RUNDIR_ROOT}/compile_${COMPILE_ID}/err"; then
            COMPILE_RESULT="FAILED: COMPILE TIMED OUT"
            FAIL_LOG="${RUNDIR_ROOT}/compile_${COMPILE_ID}/err"
          fi
        else
          COMPILE_RESULT="PASS"
          if [[ ${COMPILER} == "intel" ]]; then
            COMPILE_NUM_WARNINGS=$(grep -c ": warning #" "${RUNDIR_ROOT}/compile_${COMPILE_ID}/err" || true)
            COMPILE_NUM_REMARKS=$(grep -c ": remark #" "${RUNDIR_ROOT}/compile_${COMPILE_ID}/err" || true)
            if [[ ${COMPILE_NUM_WARNINGS} -gt 0 || ${COMPILE_NUM_REMARKS} -gt 0 ]]; then
               COMPILE_WARNINGS+=" ("
               [[ ${COMPILE_NUM_WARNINGS} -gt 0 ]] && COMPILE_WARNINGS+=" ${COMPILE_NUM_WARNINGS} warnings"
               [[ ${COMPILE_NUM_REMARKS}  -gt 0 ]] && COMPILE_WARNINGS+=" ${COMPILE_NUM_REMARKS} remarks"
               COMPILE_WARNINGS+=" )"
            fi
          elif [[ ${COMPILER} == "gnu" ]]; then
            COMPILE_NUM_WARNINGS=$(grep -c "^Warning: " "${RUNDIR_ROOT}/compile_${COMPILE_ID}/err" || true)
            if [[ ${COMPILE_NUM_WARNINGS} -gt 0 ]]; then
               COMPILE_WARNINGS+=" ("
               [[ ${COMPILE_NUM_WARNINGS} -gt 0 ]] && COMPILE_WARNINGS+=" ${COMPILE_NUM_WARNINGS} warnings"
               COMPILE_WARNINGS+=" )"
            fi
          fi
          TIME_FILE="${LOG_DIR}/compile_${COMPILE_ID}_timestamp.txt"
          if [[ -f "${TIME_FILE}" ]]; then
            while read -r times || [[ -n "${times}" ]]; do
                times="${times#"${times%%[![:space:]]*}"}"

                DATE1=$(cut -d ',' -f2 <<< "${times}")
                DATE2=$(cut -d ',' -f3 <<< "${times}")
                DATE3=$(cut -d ',' -f4 <<< "${times}")
                DATE4=$(cut -d ',' -f5 <<< "${times}")

                COMPILE_TIME=$(date --date=@$((DATE3 - DATE2)) +'%M:%S')
                RT_COMPILE_TIME=$(date --date=@$((DATE4 - DATE1)) +'%M:%S')

            done < "${TIME_FILE}"

          fi
        fi
        echo >> "${REGRESSIONTEST_LOG}"
        echo "${COMPILE_RESULT} -- COMPILE '${COMPILE_ID}' [${RT_COMPILE_TIME}, ${COMPILE_TIME}]${COMPILE_WARNINGS}" >> "${REGRESSIONTEST_LOG}"
        [[ -n ${FAIL_LOG} ]] && FAILED_COMPILES+=("COMPILE ${COMPILE_ID}: ${COMPILE_RESULT}")
        [[ -n ${FAIL_LOG} ]] && FAILED_COMPILE_LOGS+=("${FAIL_LOG}")
      fi

    elif [[ ${line} =~ RUN ]]; then

      if [[ ${COMPILE_ONLY} == true ]]; then
        continue
      fi

      RMACHINES=$(cut -d '|' -f3 <<< "${line}")
      RMACHINES=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${RMACHINES}")
      TEST_NAME=$(cut -d '|' -f2 <<< "${line}")
      TEST_NAME=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${TEST_NAME}")
      GEN_BASELINE=$(cut -d '|' -f4 <<< "${line}")
      GEN_BASELINE=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${GEN_BASELINE}")

      machines_allow_run "${RMACHINES}" && valid_test=true

      if [[ (${CONTAINER_USE} == true || ${COMMUNITY_PLATFORM_USE} == true) && ${valid_test} == true ]]; then
        RT_COMPILER=${COMPILER}
        resolve_container_image || valid_test=false
      fi

      if [[ ${valid_test} == true ]]; then
        TEST_COUNTER=$((TEST_COUNTER+1))
        GETMEMFROMLOG=""
        FAIL_LOG=""
        TEST_RESULT=""
        TIME_FILE=""
        TEST_TIME=""
        RT_TEST_TIME=""
        RT_TEST_MEM=""
        if [[ ${CREATE_BASELINE} == true && ${GEN_BASELINE} != "baseline" ]]; then
          TEST_RESULT="SKIPPED: TEST DOES NOT GENERATE BASELINE"
          SKIPPED_TESTS+=("TEST ${TEST_NAME}_${COMPILER}: ${TEST_RESULT}")
        elif [[ ${COMPILE_RESULT} =~ FAILED ]]; then
          TEST_RESULT="SKIPPED: ASSOCIATED COMPILE FAILED"
          SKIPPED_TESTS+=("TEST ${TEST_NAME}_${COMPILER}: ${TEST_RESULT}")
        elif [[ ! -f "${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log" ]]; then
          TEST_RESULT="FAILED: UNABLE TO START TEST"
          FAIL_LOG="N/A"
        elif [[ -f fail_test_${TEST_NAME}_${COMPILER} ]]; then
          if [[ -f "${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log" ]]; then
            if grep -q "FAIL" "${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log"; then
              TEST_RESULT="FAILED: UNABLE TO COMPLETE COMPARISON"
              FAIL_LOG="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"
            # We need to catch a "PASS" in rt_*.log even if a fail_test_* files exists
            # I am not sure why this can happen.
            elif grep -q "PASS" "${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log"; then
              TEST_RESULT="PASS"
            else
              TEST_RESULT="FAILED: UNSUCCESSFUL BASELINE COMPARISON"
              FAIL_LOG="${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log"
            fi
          else
            TEST_RESULT="FAILED: RUN DID NOT COMPLETE"
            FAIL_LOG="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"
          fi
          if grep -q "quota" "${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"; then
            TEST_RESULT="FAILED: DISK QUOTA ISSUE"
            FAIL_LOG="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}.log"
          elif grep -q "TIME LIMIT" "${RUNDIR_ROOT}/${TEST_NAME}_${COMPILER}/err"; then
            TEST_RESULT="FAILED: TEST TIMED OUT"
            FAIL_LOG="${RUNDIR_ROOT}/${TEST_NAME}_${COMPILER}/err"
          fi
        else
          TEST_RESULT="PASS"
        fi
        if [[ ${TEST_RESULT} == "PASS" ]]; then
          TIME_FILE="${LOG_DIR}/run_${TEST_NAME}_${COMPILER}_timestamp.txt"
          GETMEMFROMLOG=$(grep "The maximum resident set size" "${LOG_DIR}/rt_${TEST_NAME}_${COMPILER}.log")
          RT_TEST_MEM=$(echo "${GETMEMFROMLOG:9:${#GETMEMFROMLOG}-1}" | tr -dc '0-9')
          RT_TEST_MEM=$((RT_TEST_MEM/1000))
          if [[ -f "${TIME_FILE}" ]]; then
            while read -r times || [[ -n "${times}" ]]; do
                times="${times#"${times%%[![:space:]]*}"}"

                DATE1=$(cut -d ',' -f2 <<< "${times}")
                DATE2=$(cut -d ',' -f3 <<< "${times}")
                DATE3=$(cut -d ',' -f4 <<< "${times}")
                DATE4=$(cut -d ',' -f5 <<< "${times}")

                TEST_TIME=$(date --date=@$((DATE3 - DATE2)) +'%M:%S')
                RT_TEST_TIME=$(date --date=@$((DATE4 - DATE1)) +'%M:%S')

            done < "${TIME_FILE}"
          fi
        fi

        echo "${TEST_RESULT} -- TEST '${TEST_NAME}_${COMPILER}' [${RT_TEST_TIME}, ${TEST_TIME}](${RT_TEST_MEM} MB)" >> "${REGRESSIONTEST_LOG}"
        [[ -n ${FAIL_LOG} ]] && FAILED_TESTS+=("TEST ${TEST_NAME}_${COMPILER}: ${TEST_RESULT}")
        [[ -n ${FAIL_LOG} ]] && FAILED_TEST_LOGS+=("${FAIL_LOG}")
        [[ -n ${FAIL_LOG} ]] && FAILED_TEST_ID+=("${TEST_NAME} ${COMPILER}")
      fi
    fi
  done < "${TESTS_FILE}"

  elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )

  cat << EOF >> "${REGRESSIONTEST_LOG}"

SYNOPSIS:
Starting Date/Time: ${TEST_START_TIME}
Ending Date/Time: ${TEST_END_TIME}
Total Time: ${elapsed_time}
Compiles Completed: $((COMPILE_COUNTER-${#FAILED_COMPILES[@]}))/${COMPILE_COUNTER}
Tests Completed: $((TEST_COUNTER-${#FAILED_TESTS[@]}-${#SKIPPED_TESTS[@]}))/${TEST_COUNTER}
EOF
  # PRINT FAILED COMPILES
  if [[ "${#FAILED_COMPILES[@]}" -ne "0" ]]; then
    echo "Failed Compiles:" >> "${REGRESSIONTEST_LOG}"
    for i in "${!FAILED_COMPILES[@]}"; do
      echo "* ${FAILED_COMPILES[${i}]}" >> "${REGRESSIONTEST_LOG}"
      echo "-- LOG: ${FAILED_COMPILE_LOGS[${i}]}" >> "${REGRESSIONTEST_LOG}"
    done
  fi

  # PRINT FAILED TESTS
  if [[ "${#FAILED_TESTS[@]}" -ne "0" ]]; then

    echo "Failed Tests:" >> "${REGRESSIONTEST_LOG}"
    for j in "${!FAILED_TESTS[@]}"; do
      echo "* ${FAILED_TESTS[${j}]}" >> "${REGRESSIONTEST_LOG}"
      echo "-- LOG: ${FAILED_TEST_LOGS[${j}]}" >> "${REGRESSIONTEST_LOG}"
    done

  fi

  # WRITE FAILED_TEST_ID LIST TO TEST_CHANGES_LOG
  if [[ "${#FAILED_TESTS[@]}" -ne "0" ]]; then
    for item in "${FAILED_TEST_ID[@]}"; do
      echo "${item}" >> "${TEST_CHANGES_LOG}"
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
    [[ ${ROCOTO} == true ]] && rm -f "${ROCOTO_XML}" "${ROCOTO_DB}" "${ROCOTO_STATE}" ./*_lock.db
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

# Whether a COMPILE/RUN line's MACHINES field allows it under the real
# MACHINE_ID and the current -p/-P state. "+<tag>"/"-<tag>" (PLATFORM_TAG --
# 'container' for -p, or the -P file's declared platform name) explicitly
# mark a line to run (or not run) on the current container/community
# platform.
machines_allow_run() {
  local machines=$1
  local has_tag=false
  local has_no_tag=false
  [[ ${machines} == *"+${PLATFORM_TAG}"* ]] && has_tag=true
  [[ ${machines} == *"-${PLATFORM_TAG}"* ]] && has_no_tag=true

  if [[ ${CONTAINER_USE} == true || ${COMMUNITY_PLATFORM_USE} == true ]]; then
    [[ ${has_tag} == true && ${has_no_tag} == false ]] && return 0 || return 1
  fi

  local native_machines=${machines//"+${PLATFORM_TAG}"/}
  native_machines=${native_machines//"-${PLATFORM_TAG}"/}
  native_machines=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${native_machines}")
  [[ ${native_machines} == '+' || ${native_machines} == '-' ]] && native_machines=''

  [[ ${native_machines} == '' ]] && return 0
  if [[ ${native_machines} == -* ]]; then
    [[ ${native_machines} =~ ${MACHINE_ID} ]] && return 1 || return 0
  elif [[ ${native_machines} == +* ]]; then
    [[ ${native_machines} =~ ${MACHINE_ID} ]] && return 0 || return 1
  else
    echo "MACHINES=|${machines}|" >&2
    die "MACHINES spec must be either an empty string or start with either '+' or '-'"
  fi
}

# Resolves RT_CONTAINER_IMG for the current RT_COMPILER; returns 1 if no
# image is staged for this machine/compiler (caller skips the line), or -- for
# -P -- if the line's compiler doesn't match the community platform's one
# declared compiler. A community platform may have no container at all (a
# native stack), in which case RT_CONTAINER_IMG is left empty.
resolve_container_image() {
  RT_CONTAINER_IMG=''

  if [[ ${COMMUNITY_PLATFORM_USE} == true ]]; then
    [[ ${RT_COMPILER} == "${COMMUNITY_PLATFORM_COMPILER}" ]] || return 1
    RT_CONTAINER_IMG=${COMMUNITY_PLATFORM_CONTAINER_IMG}
  else
    local img_name=''
    case ${RT_COMPILER} in
      intel) img_name=${CONTAINER_IMG_INTEL} ;;
      gnu)   img_name=${CONTAINER_IMG_GNU} ;;
      *)     die "resolve_container_image: unexpected RT_COMPILER='${RT_COMPILER}'" ;;
    esac
    [[ -n ${CONTAINER_PATH} && -n ${img_name} ]] || return 1
    RT_CONTAINER_IMG=${CONTAINER_PATH}/${img_name}
  fi

  [[ -z ${RT_CONTAINER_IMG} ]] && return 0

  [[ ${DRY_RUN} == true ]] && return 0

  [[ -f ${RT_CONTAINER_IMG} ]] || return 1
  [[ -f ${PATHTR}/modulefiles/ufs_container.${RT_COMPILER}.lua ]] \
    || die "modulefiles/ufs_container.${RT_COMPILER}.lua not found under ${PATHTR}"
  [[ -f ${PATHTR}/modulefiles/ufs_container.runtime.lua ]] \
    || die "modulefiles/ufs_container.runtime.lua not found under ${PATHTR}"
  return 0
}

# Parses a -P community-platform definition file (a 4-line pipe-delimited
# header; rt.conf remains the only test source, so there is no compile/test
# list here). Sets MACHINE_ID, RT_COMPILER (the platform's one and only
# compiler), and the platform's paths/scheduler info.
parse_platform_def() {
  local file=$1
  local line header_lines_read=0
  local f1 f2 f3 f4 f5 f6 _rest
  while IFS= read -r line || [[ -n "${line}" ]]; do
    line=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${line}")
    [[ -z ${line} ]] && continue
    [[ ${line} == \#* ]] && continue

    case ${header_lines_read} in
      0)
        IFS='|' read -r f1 f2 f3 f4 _rest <<< "${line}"
        MACHINE_ID=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f1:-}")
        COMMUNITY_PLATFORM_COMPILER=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f2:-}")
        COMMUNITY_PLATFORM_CONTAINER_IMG=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f3:-}")
        CONTAINER_BIND_DIRS=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f4:-}")
        ;;
      1)
        IFS='|' read -r f1 f2 f3 f4 f5 f6 _rest <<< "${line}"
        CONTAINER_TPN=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f1:-}")
        SCHEDULER=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f2:-}")
        ACCNR=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f3:-${ACCNR}}")
        PARTITION=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f4:-}")
        QUEUE=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f5:-}")
        MPI_LAUNCH=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f6:-mpirun}")
        ;;
      2)
        RUNDIR_ROOT=${line}
        ;;
      3)
        IFS='|' read -r f1 f2 f3 f4 _rest <<< "${line}"
        INPUTDATA_ROOT=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f1:-}")
        INPUTDATA_ROOT_WW3=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f2:-}")
        INPUTDATA_LM4=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f3:-}")
        INPUTDATA_GFSv17opn=$(sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' <<< "${f4:-}")
        ;;
    esac
    header_lines_read=$((header_lines_read + 1))
    [[ ${header_lines_read} -ge 4 ]] && break
  done < "${file}"

  [[ ${header_lines_read} -ge 4 ]] || die "${file}: expected 4 header lines, found ${header_lines_read}"
  [[ -n ${MACHINE_ID} ]] || die "${file}: platform name (header line 1, field 1) is required"
  [[ -n ${COMMUNITY_PLATFORM_COMPILER} ]] || die "${file}: compiler (header line 1, field 2) is required"
  [[ -n ${SCHEDULER} ]] || die "${file}: scheduler (header line 2, field 2) is required"
  [[ -n ${RUNDIR_ROOT} ]] || die "${file}: RUNDIR_ROOT (header line 3) is required"
  [[ -n ${INPUTDATA_ROOT} ]] || die "${file}: INPUTDATA_ROOT (header line 4, field 1) is required"
}

create_or_run_compile_task() {
  cat << EOF > "${RUNDIR_ROOT}/compile_${COMPILE_ID}.env"
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
export RTVERBOSE=${RTVERBOSE}
export CONTAINER_IMG=${RT_CONTAINER_IMG}
export CONTAINER_BIND_FLAGS="${CONTAINER_BIND_FLAGS}"
export COMMUNITY_PLATFORM=${COMMUNITY_PLATFORM_USE}
EOF

  if [[ ${CONTAINER_USE} == true || ${COMMUNITY_PLATFORM_USE} == true ]]; then
    cat << EOF >> "${RUNDIR_ROOT}/compile_${COMPILE_ID}.env"
export TPN=${CONTAINER_TPN:-}
EOF
  fi

  if [[ ${ROCOTO} == true ]]; then
    rocoto_create_compile_task
  elif [[ ${ECFLOW} == true ]]; then
    ecflow_create_compile_task
  else
    echo "rt.sh: Running compile ${COMPILE_ID}"
    ./run_compile.sh "${PATHRT}" "${RUNDIR_ROOT}" "${MAKE_OPT}" "${COMPILE_ID}" > "${LOG_DIR}/compile_${COMPILE_ID}.log" 2>&1
    echo "rt.sh: Compile ${COMPILE_ID} completed."
  fi

  RT_SUFFIX=""
  BL_SUFFIX=""
}

rt_35d() {
  echo "rt.sh: Running 35day Regression Test..."
  local sy
  local sm
if [[ ${TEST_NAME} =~ '35d' ]] ; then
  sy=$(cut -c 1-4 <<< "${DATE_35D}")
  sm=$(cut -c 5-6 <<< "${DATE_35D}")
  local new_test_name="tests/${TEST_NAME}_${DATE_35D}"
  rm -f "${new_test_name}"
  cp tests/"${TEST_NAME}" "${new_test_name}"

  sed -i -e "s/\(export SYEAR\)/\1=\"${sy}\"/" "${new_test_name}"
  sed -i -e "s/\(export SMONTH\)/\1=\"${sm}\"/" "${new_test_name}"

  TEST_NAME=${new_test_name#tests/}
fi
}

handle_error() {
  echo "rt.sh: Getting error information..."
  local exit_code=$1
  local exit_line=$2
  echo "Exited at line ${exit_line} having code ${exit_code}"
  rt_trap
}

rt_trap() {
  echo "rt.sh: Exited abnormally, killing workflow and cleaning up"
  trap "" SIGINT
  [[ ${ROCOTO:-false} == true ]] && rocoto_kill
  [[ ${ECFLOW:-false} == true ]] && ecflow_kill
  cleanup
}

cleanup() {
  echo "rt.sh: Cleaning up..."
  awk_info=$(awk '{print $2}' < "${LOCKDIR}/PID")
  [[ ${awk_info} == "$$" ]] && rm -rf "${LOCKDIR}"
  [[ ${ECFLOW:-false} == true ]] && ecflow_stop
  trap 0
  echo "rt.sh: Exiting."
  exit
}

trap '{ echo "rt.sh interrupted"; rt_trap ; }' INT
trap '{ echo "rt.sh quit"; rt_trap ; }' QUIT
trap '{ echo "rt.sh terminated"; rt_trap ; }' TERM
trap '{ handle_error $? $LINENO ; }' ERR
trap '{ echo "rt.sh finished"; cleanup ; }' EXIT


# PATHRT - Path to regression tests directory
PATHRT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
readonly PATHRT
cd "${PATHRT}"

# PATHTR - Path to nmmb trunk directory
PATHTR=$( cd "${PATHRT}/.." && pwd )
readonly PATHTR

# make sure only one instance of rt.sh is running
readonly LOCKDIR="${PATHRT}"/lock
HOSTNAME_IN=$(hostname)
if mkdir "${LOCKDIR}" ; then
  echo "${HOSTNAME_IN}" $$ > "${LOCKDIR}/PID"
else
  echo "Only one instance of rt.sh can be running at a time"
  exit 1
fi

ls -l detect_machine.sh rt_utils.sh
source rt_utils.sh

CREATE_BASELINE=false
ROCOTO=false
ECFLOW=false
KEEP_RUNDIR=false
TEST_35D=false
export skip_check_results=false
export delete_rundir=false

COMPILE_ONLY=false
RTPWD_NEW_BASELINE=false
TESTS_FILE='rt.conf'
NEW_BASELINES_FILE=''
DEFINE_CONF_FILE=false
RUN_SINGLE_TEST=false
RTVERBOSE=false
export RTVERBOSE
export STOP_ECFLOW_AT_END=false
export DRY_RUN=false
ACCNR=${ACCNR:-""}

# -p: build/run inside the GNU/Intel container staged on this Tier 1 host.
CONTAINER_USE=false
CONTAINER_SUFFIX=''

# Default image filenames; overridden per host below where needed.
CONTAINER_IMG_INTEL='rocky9-oneapi2024.2-ss192.sif'
CONTAINER_IMG_GNU='rocky9-gcc13-ss192-ompi416.sif'
# Set per host in the "case ${MACHINE_ID}" block below.
CONTAINER_PATH=''
CONTAINER_BIND_DIRS=''
CONTAINER_TPN=''
RT_CONTAINER_IMG=''

# -P <file>: build/run on a community (non-Tier-1) platform defined in <file>.
# PLATFORM_TAG is the "+<tag>"/"-<tag>" name machines_allow_run() checks in
# rt.conf; it stays 'container' for -p, and becomes the platform's own
# declared name for -P.
COMMUNITY_PLATFORM_USE=false
COMMUNITY_PLATFORM_FILE=''
COMMUNITY_PLATFORM_COMPILER=''
COMMUNITY_PLATFORM_CONTAINER_IMG=''
PLATFORM_TAG='container'

while getopts ":a:b:cl:mn:dwkpP:reovhx" opt; do
  case ${opt} in
    a)
      ACCNR=${OPTARG}
      ;;
    b)
      NEW_BASELINES_FILE=${OPTARG}
      ;;
    c)
      CREATE_BASELINE=true
      ;;
    l)
      DEFINE_CONF_FILE=true
      TESTS_FILE=${OPTARG}
      grep -q '[^[:space:]]' < "${TESTS_FILE}" ||  die "${TESTS_FILE} empty, exiting..."
      ;;
    o)
      COMPILE_ONLY=true
      ;;
    p)
      CONTAINER_USE=true
      CONTAINER_SUFFIX='_container'
      ;;
    P)
      COMMUNITY_PLATFORM_USE=true
      COMMUNITY_PLATFORM_FILE=${OPTARG}
      [[ -s ${COMMUNITY_PLATFORM_FILE} ]] || die "${COMMUNITY_PLATFORM_FILE} empty or not found, exiting..."
      ;;
    m)
      # redefine RTPWD to point to newly created baseline outputs
      RTPWD_NEW_BASELINE=true
      ;;
    n)
      RUN_SINGLE_TEST=true
      IFS=' ' read -r -a SINGLE_OPTS <<< "${OPTARG}"

      if [[ ${#SINGLE_OPTS[@]} != 2 ]]; then
        die 'The -n option needs [testname] AND [compiler] in quotes, i.e. -n "control_p8 intel"'
      fi

      SRT_NAME="${SINGLE_OPTS[0]}"
      SRT_COMPILER="${SINGLE_OPTS[1]}"

      if [[ "${SRT_COMPILER}" != "intel" ]] && [[ "${SRT_COMPILER}" != "intelllvm" ]] && [[ "${SRT_COMPILER}" != "gnu" ]]; then
        die "COMPILER MUST BE 'intel' OR 'intelllvm' OR 'gnu'"
      fi
      ;;
    d)
      export delete_rundir=true
      AWK_OUT=$(awk -F "|" '{print $5}' rt.conf)
      grep "\S" <<< "${AWK_OUT}" > keep_tests.tmp
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
    v)
      RTVERBOSE=true
      ;;
    x)
      DRY_RUN=true
      ;;
    h)
      usage
      die ""
      ;;
    \?)
      usage
      die "Invalid option: -${OPTARG}"
      ;;
    :)
      usage
      die "Option -${OPTARG} requires an argument."
      ;;
    *)
      usage
      die "Arguments are required."
      ;;
  esac
done

#Check to error out if incompatible options are chosen together
[[ ${KEEP_RUNDIR} == true && ${delete_rundir} == true ]] && die "-k and -d options cannot be used at the same time"
[[ ${ECFLOW} == true && ${ROCOTO} == true ]] && die "-r and -e options cannot be used at the same time"
[[ ${CREATE_BASELINE} == true && ${RTPWD_NEW_BASELINE} == true ]] && die "-c and -m options cannot be used at the same time"
#B&N not run together
[[ ${NEW_BASELINES_FILE} != '' && ${RUN_SINGLE_TEST} == true ]] && die "-b and -n options cannot be used at the same time"
#P&p not run together; a community platform run is always sequential and never compares baselines
[[ ${COMMUNITY_PLATFORM_USE} == true && ${CONTAINER_USE} == true ]] && die "-p and -P options cannot be used at the same time"
if [[ ${COMMUNITY_PLATFORM_USE} == true ]]; then
  [[ ${ROCOTO} == false ]] || die "-P should not be used with -r"
  [[ ${ECFLOW} == false ]] || die "-P should not be used with -e"
  [[ ${CREATE_BASELINE} == false ]] || die "-P should not be used with -c"
  [[ ${RTPWD_NEW_BASELINE} == false ]] || die "-P should not be used with -m"
fi

if [[ ${DRY_RUN} == true ]]; then
   [[ ${NEW_BASELINES_FILE} == '' ]] || die "-x should not be used with -b"
   [[ ${CREATE_BASELINE} == false ]] || die "-x should not be used with -c"
   [[ ${delete_rundir} == false ]] || die "-x should not be used with -d"
   [[ ${ECFLOW} == false ]] || die "-x should not be used with -e"
   [[ ${RTPWD_NEW_BASELINE} == false ]] || die "-x should not be used with -m"
   [[ ${COMPILE_ONLY} == false ]] || die "-x should not be used with -o"
   [[ ${ROCOTO} == false ]] || die "-x should not be used with -r"
   [[ ${skip_check_results} == false ]] || die "-x should not be used with -w"
fi

if [[ ${RTVERBOSE} == true ]]; then
  set -x
fi

if [[ -z "${ACCNR}" ]]; then
  echo "Please use -a <account> to set group account to use on HPC"
  exit 1
fi

if [[ ${COMMUNITY_PLATFORM_USE} == true ]]; then
  parse_platform_def "${COMMUNITY_PLATFORM_FILE}"
  PLATFORM_TAG=${MACHINE_ID}
else
  source detect_machine.sh
fi
# shellcheck disable=SC1091
source module-setup.sh

# Display the machine and account using the format detect_machine.sh used:
echo "Machine: ${MACHINE_ID}"
echo "Account: ${ACCNR}"

if [[ ${COMMUNITY_PLATFORM_USE} == true ]]; then
  # Community platform: no per-host case block -- everything came from the
  # -P file. Fill in what unrelated downstream code still references.
  DISKNM=''
  STMP=${RUNDIR_ROOT}
  PTMP=${RUNDIR_ROOT}
  COMPILE_QUEUE=${QUEUE}
  ROCOTO=false
  ECFLOW=false
  export skip_check_results=true
else
case ${MACHINE_ID} in
  wcoss2|acorn)
    echo "rt.sh: Setting up WCOSS2/Acorn"
    if [[ "${ECFLOW:-false}" == true ]] ; then
      module load ecflow/5.6.0.13
    fi
    module load intel/19.1.3.304 python/3.8.6

    DISKNM="/lfs/h2/emc/nems/noscrub/emc.nems/RT"
    QUEUE="dev"
    COMPILE_QUEUE="dev"
    if [[ "${ROCOTO:-false}" == true ]] ; then
      ROCOTO_SCHEDULER="pbs"
    fi
    PARTITION=
    STMP="/lfs/h2/emc/ptmp"
    PTMP="/lfs/h2/emc/ptmp"
    SCHEDULER="pbs"

    # no container image staged yet
    # CONTAINER_PATH=                 # directory holding the *.sif images
    # CONTAINER_BIND_DIRS=         # comma-separated host dirs to bind
    # CONTAINER_TPN=128
    ;;
  gaeac5)
    echo "rt.sh: Setting up gaea c5..."
    if [[ "${ROCOTO:-false}" == true ]] ; then
      module use /ncrc/proj/epic/rocoto/modulefiles
      module load rocoto
      ROCOTO_SCHEDULER="slurm"
    fi

    export LD_PRELOAD=/usr/lib64/libstdc++.so.6
    module use /ncrc/proj/epic/spack-stack/c5/spack-stack-1.9.1/envs/ue-intel-2023.2.0/install/modulefiles/Core
    module load stack-intel/2023.2.0
    module load cray-mpich/8.1.30
    module load python/3.11
    module use /ncrc/proj/epic/spack-stack/modulefiles
    if [[ "${ECFLOW:-false}" == true ]] ; then
      module load ecflow/5.8.4
      ECF_HOST=$(hostname)
      ECF_PORT=$(( $(id -u) + 1500 ))
      export ECF_PORT ECF_HOST
    fi

    DISKNM=/gpfs/f5/epic/world-shared/UFS-WM_RT
    QUEUE=normal
    COMPILE_QUEUE=normal
    PARTITION=c5
    dprefix=${dprefix:-/gpfs/f5/${ACCNR}/scratch/${USER}}
    STMP=${STMP:-${dprefix}/RT_BASELINE}
    PTMP=${PTMP:-${dprefix}/RT_RUNDIRS}

    SCHEDULER="slurm"

    # gaea c5 no longer supported
    # CONTAINER_PATH=                 # directory holding the *.sif images
    # CONTAINER_BIND_DIRS=         # comma-separated host dirs to bind
    # CONTAINER_TPN=128
    ;;
  gaeac6)
    echo "rt.sh: Setting up gaea c6..."
    if [[ "${ROCOTO:-false}" == true ]] ; then
      module use /ncrc/proj/epic/c6/modulefiles
      module load rocoto/1.3.7
      ROCOTO_SCHEDULER="slurm"
    fi

    export LD_PRELOAD=/usr/lib64/libstdc++.so.6
    module use /ncrc/proj/epic/spack-stack/c6/spack-stack-1.9.2/envs/ue-intel-2023.2.0/install/modulefiles/Core
    module load stack-intel/2023.2.0
    module load cray-mpich/8.1.30
    module load python/3.11
    if [[ "${ECFLOW:-false}" == true ]] ; then
      module use /ncrc/proj/epic/spack-stack/c6/spack-stack-1.9.2/envs/ue-intel-2023.2.0/install/modulefiles/gcc/12.3.0
      module load ecflow/5.11.4
      ECF_HOST=$(hostname)
      ECF_PORT=$(( $(id -u) + 1500 ))
      export ECF_PORT ECF_HOST
    fi

    DISKNM=/gpfs/f6/bil-fire8/world-shared/role.epic/UFS-WM_RT
    QUEUE=normal
    COMPILE_QUEUE=normal
    PARTITION=c6
    dprefix=${dprefix:-/gpfs/f6/${ACCNR}/proj-shared/${USER}}
    STMP=${STMP:-${dprefix}/RT_BASELINE}
    PTMP=${PTMP:-${dprefix}/RT_RUNDIRS}

    SCHEDULER="slurm"

    CONTAINER_PATH=/gpfs/f6/bil-fire8/world-shared/containers   # directory holding the *.sif images
    CONTAINER_BIND_DIRS="/gpfs,/ncrc/home2"                  # comma-separated host dirs to bind
    CONTAINER_TPN=192
    ;;
  hera)
    echo "rt.sh: Setting up hera..."
    if [[ "${ROCOTO:-false}" == true ]] ; then
      module load rocoto
      ROCOTO_SCHEDULER=slurm
    fi

    if [[ "${ECFLOW:-false}" == true ]] ; then
      module load ecflow/5.11.4
    fi

    QUEUE="batch"
    COMPILE_QUEUE="batch"

    PARTITION=
    dprefix=${dprefix:-"/scratch3/NCEPDEV/stmp/${USER}"}
    DISKNM="/scratch3/NAGAPE/epic/role.epic/UFS-WM_RT"
    STMP="${dprefix}/RT_BASELINE"
    PTMP="${dprefix}/RT_RUNDIRS"

    SCHEDULER=slurm

    # hera no longer supported
    # CONTAINER_PATH=                 # directory holding the *.sif images
    # CONTAINER_BIND_DIRS=         # comma-separated host dirs to bind
    # CONTAINER_TPN=40
    ;;
  ursa)
    echo "rt.sh: Setting up ursa..."
    if [[ "${ROCOTO:-false}" == true ]] ; then
      module load rocoto
      ROCOTO_SCHEDULER=slurm
    fi

    if [[ "${ECFLOW:-false}" == true ]] ; then
      module load ecflow/5.11.4
      ECF_HOST="uecflow01"
      ECF_PORT="$(( $(id -u) + 1500 ))"
      export ECF_HOST ECF_PORT
    fi

    QUEUE="batch"
    COMPILE_QUEUE="batch"

    PARTITION="u1-compute"
    dprefix="/scratch4/NCEPDEV/stmp/${USER}"
    if [[ "${ACCNR}" == 'epic' ]] ; then
      dprefix="/scratch4/NAGAPE/epic/${USER}/stmp"
    fi
    DISKNM="/scratch4/NAGAPE/epic/role-epic/UFS-WM_RT"
    STMP="${STMP:-${dprefix}/RT_BASELINE}"
    PTMP="${PTMP:-${dprefix}/RT_RUNDIRS}"

    SCHEDULER=slurm

    CONTAINER_PATH=/scratch3/NCEPDEV/nems/role.epic/containers     # directory holding the *.sif images
    CONTAINER_BIND_DIRS="/scratch3,/scratch4"                   # comma-separated host dirs to bind
    CONTAINER_TPN=192

    ;;
  orion)
    echo "rt.sh: Setting up orion..."

    if [[ "${ROCOTO:-false}" == true ]] ; then
      module load contrib ruby/3.2.3 rocoto/1.3.7
      ROCOTO_SCHEDULER="slurm"
    fi

    module use /work/noaa/epic/role-epic/spack-stack/orion/modulefiles
    if [[ "${ECFLOW:-false}" == true ]] ; then
      module load ecflow/5.8.4
      ECF_HOST=$(hostname)
      ECF_PORT="$(( $(id -u) + 1500 ))"
      export ECF_PORT ECF_HOST
    fi

    QUEUE="batch"
    COMPILE_QUEUE="batch"
    PARTITION="orion"
    dprefix=${dprefix:-"/work/noaa/stmp/${USER}"}
    DISKNM="/work2/noaa/epic/UFS-WM_RT"
    STMP="${dprefix}/stmp"
    PTMP="${dprefix}/stmp"

    SCHEDULER="slurm"

    cp fv3_conf/fv3_slurm.IN_orion fv3_conf/fv3_slurm.IN
    cp fv3_conf/compile_slurm.IN_orion fv3_conf/compile_slurm.IN

    CONTAINER_PATH=/work/noaa/epic/role-epic/contrib/containers   # directory holding the *.sif images
    CONTAINER_BIND_DIRS="/work,/work2,/local"                  # comma-separated host dirs to bind
    CONTAINER_TPN=40
    ;;
  hercules)
    echo "rt.sh: Setting up hercules..."
    if [[ "${ROCOTO:-false}" == true ]] ; then
      module load contrib rocoto
      ROCOTO_SCHEDULER="slurm"
    fi

    module use /work/noaa/epic/role-epic/spack-stack/hercules/modulefiles
    if [[ "${ECFLOW:-false}" == true ]] ; then
      module load ecflow/5.8.4
      ECF_HOST=$(hostname)
      ECF_PORT="$(( $(id -u) + 1500 ))"
      export ECF_PORT ECF_HOST
    fi

    QUEUE="batch"
    COMPILE_QUEUE="batch"
    PARTITION="hercules"
    dprefix=${dprefix:-"/work2/noaa/stmp/${USER}"}
    DISKNM="/work2/noaa/epic/hercules/UFS-WM_RT"
    STMP="${dprefix}/stmp"
    PTMP="${dprefix}/stmp"

    SCHEDULER="slurm"
    cp fv3_conf/fv3_slurm.IN_hercules fv3_conf/fv3_slurm.IN
    cp fv3_conf/compile_slurm.IN_hercules fv3_conf/compile_slurm.IN

    CONTAINER_PATH=/work/noaa/epic/role-epic/contrib/containers
    CONTAINER_BIND_DIRS="/work,/work2,/local"
    CONTAINER_TPN=80
    ;;
  derecho)
    echo "rt.sh: Setting up derecho..."
    if [[ "${ROCOTO:-false}" == true ]] ; then
      module use /glade/work/epicufsrt/contrib/derecho/modulefiles
      module load rocoto/1.3.7
    fi
    if [[ "${ECFLOW:-false}" == true ]] ; then
      module use /glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/oneapi/2024.2.1
      module load stack-python/3.11.7
      module load ecflow/5.11.4
      ECF_HOST=$(hostname)
      ECF_PORT=$(( $(id -u) + 1500 ))
      export ECF_PORT ECF_HOST
    fi

    QUEUE="main"
    COMPILE_QUEUE="main"
    PARTITION=
    dprefix=${dprefix:-"/glade/derecho/scratch"}
    DISKNM="/glade/derecho/scratch/epicufsrt/ufs-weather-model/RT/"
    STMP="${dprefix}"
    PTMP="${dprefix}"
    SCHEDULER="pbs"
    cp fv3_conf/fv3_qsub.IN_derecho fv3_conf/fv3_qsub.IN
    cp fv3_conf/compile_qsub.IN_derecho fv3_conf/compile_qsub.IN


    if [[ "${ROCOTO:-false}" == true ]] ; then
      ROCOTO_SCHEDULER="pbspro"
    fi

    CONTAINER_IMG_GNU="rocky9-gcc13-ss192-ompi507.sif"
    CONTAINER_PATH=/glade/work/epicufsrt/contrib/containers  # directory holding the *.sif images
    CONTAINER_BIND_DIRS="/glade"                                # comma-separated host dirs to bind
    CONTAINER_TPN=128
    ;;
  noaacloud)
    echo "rt.sh: Setting up noaacloud..."
    export PATH="/contrib/EPIC/bin:${PATH}"
    module use /apps/modules/modulefiles

    if [[ "${ROCOTO:-false}" == true ]] ; then
      module load rocoto/1.3.7
      ROCOTO_SCHEDULER=slurm
    fi

    QUEUE="batch"
    COMPILE_QUEUE="batch"
    PARTITION=
    dprefix=${dprefix:-"/lustre/"}
    DISKNM="/contrib/ufs-weather-model/RT"
    STMP="${dprefix}/stmp4"
    PTMP="${dprefix}/stmp2"
    SCHEDULER="slurm"

    CONTAINER_PATH=/contrib/EPIC/containers     # directory holding the *.sif images
    CONTAINER_BIND_DIRS="/contrib,/lustre"   # comma-separated host dirs to bind
    CONTAINER_TPN=36                            # may need to be specified for different cloud platforms
    ;;
  *)
    die "Unknown machine ID, please edit detect_machine.sh file"
    ;;
esac
fi

if [[ ${CONTAINER_USE} == true && -z ${CONTAINER_PATH} ]]; then
  echo "rt.sh: WARNING -- no container images are staged for ${MACHINE_ID};"
  echo "                 all container tests will be skipped."
fi

# Resolve CONTAINER_BIND_DIRS into apptainer/singularity "-B dir" flags once.
CONTAINER_BIND_FLAGS=''
if [[ -n ${CONTAINER_BIND_DIRS} ]]; then
  IFS=',' read -r -a _container_bind_dirs <<< "${CONTAINER_BIND_DIRS}"
  for _dir in "${_container_bind_dirs[@]}"; do
    CONTAINER_BIND_FLAGS="${CONTAINER_BIND_FLAGS} -B ${_dir}"
  done
fi

mkdir -p "${STMP}/${USER}"

NEW_BASELINE=${STMP}/${USER}/FV3_RT/REGRESSION_TEST${CONTAINER_SUFFIX:-}

# Overwrite default RUNDIR_ROOT if environment variable RUNDIR_ROOT is set
RUNDIR_ROOT=${RUNDIR_ROOT:-${PTMP}/${USER}/FV3_RT}/rt_$$
mkdir -p "${RUNDIR_ROOT}"
rm -rf "${PATHRT}/run_dir"
echo "Linking ${RUNDIR_ROOT} to ${PATHRT}/run_dir"
ln -s "${RUNDIR_ROOT}" "${PATHRT}/run_dir"
echo "Run regression test in: ${RUNDIR_ROOT}"

# BEFORE MOVING ANY FURTHER LETS CHECK THAT DISKNM/STMP/PTMP ALL EXIST
# (a community platform, -P, has no DISKNM/baseline area at all -- skip)
if [[ ${COMMUNITY_PLATFORM_USE} == false ]]; then
  [[ -d ${DISKNM} ]] || die "ERROR: DISKNM: ${DISKNM} -- DOES NOT EXIST"
  [[ -d ${STMP} ]] || die "ERROR: STMP: ${STMP} -- DOES NOT EXIST"
  [[ -d ${PTMP} ]] || die "ERROR: PTMP: ${PTMP} -- DOES NOT EXIST"
fi

update_rtconf

if [[ ${TESTS_FILE} =~ '35d' ]] || [[ ${TESTS_FILE} =~ 'weekly' ]]; then
  TEST_35D=true
fi

source bl_date.conf

if [[ "${RTPWD_NEW_BASELINE}" == true ]] ; then
  RTPWD=${NEW_BASELINE}
else
  RTPWD=${RTPWD:-${DISKNM}/NEMSfv3gfs/develop-${BL_DATE}${CONTAINER_SUFFIX:-}}
fi

# A community platform (-P) always skips baseline comparison -- no baseline
# directory to check.
if [[ "${CREATE_BASELINE}" == false && ${COMMUNITY_PLATFORM_USE} == false ]] ; then
  EMPTY_CHECK=$(find "${RTPWD}/" -type d -prune -empty)
  if [[ ! -d "${RTPWD}" ]] ; then
    echo "Baseline directory does not exist:"
    echo "   ${RTPWD}"
    exit 1
  elif [[ -n ${EMPTY_CHECK} ]] ; then
    echo "Baseline directory is empty:"
    echo "   ${RTPWD}"
    exit 1
  fi
fi

INPUTDATA_ROOT=${INPUTDATA_ROOT:-${DISKNM}/NEMSfv3gfs/input-data-20260617}
INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT_WW3:-${INPUTDATA_ROOT}/WW3_input_data_20260811}
INPUTDATA_LM4=${INPUTDATA_LM4:-${INPUTDATA_ROOT}/LM4_input_data}
INPUTDATA_GFSv17opn=${INPUTDATA_GFSv17opn:-${DISKNM}/NEMSfv3gfs/GFSv17opn_20251014}

shift $((OPTIND-1))
if [[ $# -gt 1 ]]; then
  usage
  die ""
fi

if [[ ${CREATE_BASELINE} == true ]]; then
  # PREPARE NEW REGRESSION TEST DIRECTORY
  echo "rt.sh: Preparing RT Directory..."
  rm -rf "${NEW_BASELINE}"
  mkdir -p "${NEW_BASELINE}"

fi

if [[ ${skip_check_results} == true ]]; then
  REGRESSIONTEST_LOG=${PATHRT}/logs/RegressionTests_weekly_${MACHINE_ID}${CONTAINER_SUFFIX:-}.log
else
  REGRESSIONTEST_LOG=${PATHRT}/logs/RegressionTests_${MACHINE_ID}${CONTAINER_SUFFIX:-}.log
fi

[ -f "${REGRESSIONTEST_LOG}" ] && cp "${REGRESSIONTEST_LOG}" "${REGRESSIONTEST_LOG}.bak"
rm -f "${REGRESSIONTEST_LOG}"

TEST_START_TIME="$(date '+%Y%m%d %T')"
export TEST_START_TIME

source default_vars.sh

COMPILE_COUNTER=0
rm -f fail_test* fail_compile*

LOG_DIR=${PATHRT}/logs/log_${MACHINE_ID}${CONTAINER_SUFFIX:-}
export LOG_DIR

rm -rf "${LOG_DIR}"
mkdir -p "${LOG_DIR}"

if [[ ${ROCOTO} == true ]]; then

  echo "rt.sh: Verifying ROCOTO support..."

  case ${MACHINE_ID} in
    wcoss2|acorn)
      die "Rocoto not supported on this machine, please do not use '-r'."
      ;;
    *)
      ;;
  esac

  ROCOTORUN="$(command -v rocotorun)"
  ROCOTOSTAT="$(command -v rocotostat)"
  ROCOTOCOMPLETE="$(command -v rocotocomplete)"
  export ROCOTOCOMPLETE ROCOTOSTAT ROCOTORUN

  ROCOTO_XML=${PATHRT}/rocoto_workflow.xml
  ROCOTO_STATE=${PATHRT}/rocoto_workflow.state
  ROCOTO_DB=${PATHRT}/rocoto_workflow.db

  rm -f "${ROCOTO_XML}" "${ROCOTO_DB}" "${ROCOTO_STATE}" ./*_lock.db

  cat << EOF > "${ROCOTO_XML}"
<?xml version="1.0"?>
<!DOCTYPE workflow
[
  <!ENTITY PATHRT         "${PATHRT}">
  <!ENTITY LOG            "${LOG_DIR}">
  <!ENTITY PATHTR         "${PATHTR}">
  <!ENTITY RTPWD          "${RTPWD}">
  <!ENTITY INPUTDATA_ROOT "${INPUTDATA_ROOT}">
  <!ENTITY INPUTDATA_ROOT_WW3 "${INPUTDATA_ROOT_WW3}">
  <!ENTITY RUNDIR_ROOT    "${RUNDIR_ROOT}">
  <!ENTITY NEW_BASELINE   "${NEW_BASELINE}">
]>
<workflow realtime="F" scheduler="${ROCOTO_SCHEDULER}" taskthrottle="10">
  <cycledef>197001010000 197001010000 01:00:00</cycledef>
  <log>&LOG;/workflow.log</log>
EOF

fi

if [[ ${ECFLOW} == true ]]; then
  echo "Verifying ECFLOW support..."
  case ${MACHINE_ID} in
    noaacloud)
      die "ECFLOW not supported on this machine, please do not use '-e'."
      ;;
    *)
      ECFLOW_START="$(command -v ecflow_start.sh)"
      ;;
  esac
  export ECFLOW_START

  #export ECF_OUTPUTDIR="${PATHRT}/ecf_outputdir"
  #export ECF_COMDIR="${PATHRT}/ecf_comdir"
  #rm -rf "${ECF_OUTPUTDIR}" "${ECF_COMDIR}"
  #mkdir -p "${ECF_OUTPUTDIR}"
  #mkdir -p "${ECF_COMDIR}"
  # Default maximum number of compile and run jobs
  MAX_BUILDS=10 #Max build jobs
  MAX_JOBS=30   #Max test/run jobs
  ECF_TRIES=2   #Tries before failure

  ECFLOW_RUN=${PATHRT}/ecflow_run
  ECFLOW_SUITE=regtest_$$
  rm -rf "${ECFLOW_RUN}"
  mkdir -p "${ECFLOW_RUN}/${ECFLOW_SUITE}"
  cp head.h tail.h "${ECFLOW_RUN}"
  cat << EOF > "${ECFLOW_RUN}/${ECFLOW_SUITE}.def"
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

fi

##
## read rt.conf and then either execute the test script directly or create
## workflow description file
##

new_compile=false
in_metatask=false

[[ -f ${TESTS_FILE} ]] || die "${TESTS_FILE} does not exist"

declare -A compiles

while read -r line || [[ -n "${line}" ]]; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ ${line} == \#* ]] && continue

  if [[ ${line} == COMPILE* ]]; then

    COMPILE_NAME=$(cut -d '|' -f2 <<< "${line}")
    COMPILE_NAME=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${COMPILE_NAME}")

    RT_COMPILER=$(cut -d '|' -f3  <<< "${line}")
    RT_COMPILER=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${RT_COMPILER}")

    MAKE_OPT=$(cut -d '|' -f4  <<< "${line}")
    MAKE_OPT=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${MAKE_OPT}")

    MACHINES=$(cut -d '|' -f5  <<< "${line}")
    MACHINES=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${MACHINES}")

    CB=$(cut -d '|' -f6  <<< "${line}")
    COMPILE_ID=${COMPILE_NAME}_${RT_COMPILER}

    set +u
    if [[ -n ${compiles[${COMPILE_ID}]} ]] ; then
        echo "Error! Duplicated compilation ${COMPILE_NAME} for compiler ${RT_COMPILER}!"
        exit 1
    fi
    set -u
    compiles[${COMPILE_ID}]=${COMPILE_ID}

    [[ ${CREATE_BASELINE} == true && ${CB} != *fv3* ]] && continue

    machines_allow_run "${MACHINES}" || continue

    if [[ ${CONTAINER_USE} == true || ${COMMUNITY_PLATFORM_USE} == true ]] && ! resolve_container_image; then
      [[ ${RUN_SINGLE_TEST} == true ]] && die "No ${RT_COMPILER} container/platform match on ${MACHINE_ID} for -n test"
      echo "rt.sh: SKIP compile ${COMPILE_ID} -- compiler ${RT_COMPILER} not available on ${MACHINE_ID}"
      continue
    fi

    [[ ${DRY_RUN} == true ]] && continue

    create_or_run_compile_task
    continue

  elif [[ ${line} == RUN* ]]; then

    [[ ${COMPILE_ONLY} == true ]] && continue

    TEST_NAME=$(cut -d'|' -f2 <<< "${line}")
    TEST_NAME=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${TEST_NAME}")

    MACHINES=$(cut -d'|' -f3 <<< "${line}")
    MACHINES=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${MACHINES}")

    CB=$(cut -d'|' -f4 <<< "${line}")

    DEP_RUN=$(cut -d'|' -f5 <<< "${line}")
    DEP_RUN=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${DEP_RUN}")

    DATE_35D=$(cut -d'|' -f6 <<< "${line}")
    DATE_35D=$(sed -e 's/^ *//' -e 's/ *$//' <<< "${DATE_35D}")

    if [[ ${DEP_RUN} != '' ]]; then
      DEP_RUN=${DEP_RUN}_${RT_COMPILER}
    fi

    export TEST_ID=${TEST_NAME}_${RT_COMPILER}

    [[ -e "tests/${TEST_NAME}" ]] || die "run test file tests/${TEST_NAME} does not exist"
    [[ ${CREATE_BASELINE} == true && ${CB} != *baseline* ]] && continue

    machines_allow_run "${MACHINES}" || continue

    if [[ ${CONTAINER_USE} == true || ${COMMUNITY_PLATFORM_USE} == true ]] && ! resolve_container_image; then
      echo "rt.sh: SKIP test ${TEST_ID} -- compiler ${RT_COMPILER} not available on ${MACHINE_ID}"
      continue
    fi

    COMPILE_METATASK_NAME=${COMPILE_ID}

    # 35 day tests
    [[ ${TEST_35D} == true ]] && rt_35d

    # Avoid uninitialized RT_SUFFIX/BL_SUFFIX (see definition above)
    RT_SUFFIX=${RT_SUFFIX:-""}
    BL_SUFFIX=${BL_SUFFIX:-""}

    if [[ ${ROCOTO} == true && ${new_compile} == true ]]; then
      new_compile=false
      in_metatask=true
      cat << EOF >> "${ROCOTO_XML}"
  <metatask name="compile_${COMPILE_METATASK_NAME}_tasks"><var name="zero">0</var>
EOF
    fi

    (
      # shellcheck source=/github/workspace/tests/tests/control_c48
      source "${PATHRT}/tests/${TEST_NAME}"

      if [[ ${ESMF_THREADING} == true ]]; then
        compute_petbounds_and_tasks_esmf_threading
      else
        compute_petbounds_and_tasks_traditional_threading
      fi

      TPN=$(( TPN / THRD ))
      NODES=$(( TASKS / TPN ))
      if (( NODES * TPN < TASKS )); then
        NODES=$(( NODES + 1 ))
      fi

      PPN=$(( TASKS / NODES ))
      if (( TASKS - ( PPN * NODES ) > 0 )); then
          PPN=$((PPN + 1))
      fi

      cat << EOF > "${RUNDIR_ROOT}/run_test_${TEST_ID}.env"
export TEST_ID=${TEST_ID}
export MACHINE_ID=${MACHINE_ID}
export RT_COMPILER=${RT_COMPILER}
export RTPWD=${RTPWD}
export INPUTDATA_ROOT=${INPUTDATA_ROOT}
export INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT_WW3}
export INPUTDATA_LM4=${INPUTDATA_LM4}
export INPUTDATA_GFSv17opn=${INPUTDATA_GFSv17opn}
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
export RTVERBOSE=${RTVERBOSE}
export delete_rundir=${delete_rundir}
export WLCLK=${WLCLK}
export DRY_RUN=${DRY_RUN}
export CONTAINER_IMG=${RT_CONTAINER_IMG}
export CONTAINER_BIND_FLAGS="${CONTAINER_BIND_FLAGS}"
export COMMUNITY_PLATFORM=${COMMUNITY_PLATFORM_USE}
EOF

      if [[ ${CONTAINER_USE} == true || ${COMMUNITY_PLATFORM_USE} == true ]]; then
        cat << EOF >> "${RUNDIR_ROOT}/run_test_${TEST_ID}.env"
export TPN=${CONTAINER_TPN:-}
EOF
      fi

      if [[ ${ROCOTO} == true ]]; then
        rocoto_create_run_task
      elif [[ ${ECFLOW} == true ]]; then
        ecflow_create_run_task
      else
        echo "rt.sh: Running test ${TEST_ID} using compile ${COMPILE_ID}"
        ./run_test.sh "${PATHRT}" "${RUNDIR_ROOT}" "${TEST_NAME}" "${TEST_ID}" "${COMPILE_ID}" > "${LOG_DIR}/run_${TEST_ID}${RT_SUFFIX}.log" 2>&1
        echo "rt.sh: Run with test ${TEST_ID} completed."
      fi
    )
    continue
  else
    die "Unknown command ${line}"
  fi
done < "${TESTS_FILE}"

##
## run regression test workflow (currently Rocoto or ecFlow are supported)
##

if [[ ${ROCOTO} == true ]]; then
  if [[ ${in_metatask} == true ]]; then
    echo "  </metatask>" >> "${ROCOTO_XML}"
  fi
  echo "</workflow>" >> "${ROCOTO_XML}"
  # run rocoto workflow until done
  rocoto_run
fi

if [[ ${ECFLOW} == true ]]; then
  echo "endsuite" >> "${ECFLOW_RUN}/${ECFLOW_SUITE}.def"
  # run ecflow workflow until done
  ecflow_run
fi

# IF -c AND -b; LINK VERIFIED BASELINES TO NEW_BASELINE
if [[ ${CREATE_BASELINE} == true && ${NEW_BASELINES_FILE} != '' ]]; then
  for dir in "${RTPWD}"/*/; do
    dir=${dir%*/}
    [[ -d "${NEW_BASELINE}/${dir##*/}" ]] && continue
    ln -s "${dir%*/}" "${NEW_BASELINE}/"
  done
fi

if [[ ${DRY_RUN} == true ]]; then
  echo "Successful dry run"
  exit 0
fi

## Lets verify all tests were run and that they passed
generate_log
echo "******Regression Testing Script Completed******"
