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
    
    if [[ $line =~ COMPILE ]] ; then
      MACHINES=$(echo $line | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
      RT_COMPILER_IN=$(echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
      if [[ ${MACHINES} == '' ]]; then
        compile_line=$line
        COMPILE_LINE_USED=false
      elif [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || compile_line=$line; COMPILE_LINE_USED=false
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && compile_line=$line; COMPILE_LINE_USED=false
      fi

    fi

    if [[ $line =~ RUN ]]; then
      to_run_test=false
      tmp_test=$(echo "$line" | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
      MACHINES=$(echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
      if [[ ${MACHINES} == '' ]]; then
        to_run_test=true
      elif [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || to_run_test=true
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && to_run_test=true
      fi
      if [[ $to_run_test == true ]]; then
        TEST_IDX=$(find_match "$tmp_test $RT_COMPILER_IN" "${TEST_WITH_COMPILE[@]}")
        
        if [[ $TEST_IDX != -1 ]]; then
          if [[ $COMPILE_LINE_USED == false ]]; then
              echo -en '\n' >> $RT_TEMP_CONF
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
    fi
  done < "$TESTS_FILE"

  if [[ ! -s $RT_TEMP_CONF ]]; then
    echo "The tests listed/chosen do not exist or cannot be run on $MACHINE_ID"
    exit 1
  else
    TESTS_FILE=$RT_TEMP_CONF
  fi
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

    echo 'WCOSS2'

elif [[ $MACHINE_ID = acorn ]]; then

    echo 'ACORN'

elif [[ $MACHINE_ID = hera ]]; then

  if [[ "${ROCOTO:-false}" == true ]] ; then
    module load rocoto
    ROCOTORUN=$(which rocotorun)
    ROCOTOSTAT=$(which rocotostat)
    ROCOTOCOMPLETE=$(which rocotocomplete)
    ROCOTO_SCHEDULER=slurm
  fi

  if [[ "${ECFLOW:-false}" == true ]] ; then
    module load ecflow/5.11.4
    ECFLOW_START=ecflow_start.sh
  fi
  
else
  die "Unknown machine ID, please edit detect_machine.sh file"
fi

update_rtconf

source bl_date.conf

shift $((OPTIND-1))
[[ $# -gt 1 ]] && usage

export TEST_START_TIME="$(date '+%Y%m%d %T')"

source default_vars.sh

JOB_NR=0
COMPILE_COUNTER=0
rm -f fail_test* fail_compile*

if [[ $ROCOTO == true ]]; then
  ROCOTO_XML=${PATHRT}/rocoto_workflow.xml
  ROCOTO_STATE=${PATHRT}/rocoto_workflow.state
  ROCOTO_DB=${PATHRT}/rocoto_workflow.db
  rm -f $ROCOTO_XML $ROCOTO_DB $ROCOTO_STATE *_lock.db
fi

[[ -f $TESTS_FILE ]] || die "$TESTS_FILE does not exist"

LAST_COMPILER_NR=-9999
COMPILE_PREV=''

declare -A compiles

    #create_or_run_compile_task
    export ROCOTO_SCHEDULER
    export ACCNR
    export ROCOTO_XML
    export PATHRT
    export ROCOTO
    export ECFLOW
    export MACHINE_ID
    export RTPWD_NEW_BASELINE
    export NEW_BASELINE
    export CREATE_BASELINE
    #RTVERBOSE = false
    export RTVERBOSE
    
    python -c "import create_xml; create_xml.main_loop()"

##
## run regression test workflow (currently Rocoto or ecFlow are supported)
##

if [[ $ROCOTO == true ]]; then
  rocoto_run
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
#generate_log
