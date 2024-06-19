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
  echo "  -e  use ecFlow workflow manager (this option is not fully functional yet)"
  echo "  -h  display this help"
  echo "  -k  keep run directory after ufs_test.sh is completed"
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

rt_trap() {
  [[ ${ROCOTO:-false} == true ]] && rocoto_kill
  [[ ${ECFLOW:-false} == true ]] && ecflow_kill
  cleanup
}

cleanup() {
  PID_LOCK=$(awk '{print $2}' < "${LOCKDIR}/PID")
  [[ ${PID_LOCK} == "$$" ]] && rm -rf "${LOCKDIR}"
  [[ ${ECFLOW:-false} == true ]] && ecflow_stop
  trap 0
  exit
}

trap '{ echo "ufs_test.sh interrupted"; rt_trap ; }' INT
trap '{ echo "ufs_test.sh quit"; rt_trap ; }' QUIT
trap '{ echo "ufs_test.sh terminated"; rt_trap ; }' TERM
trap '{ echo "ufs_test.sh error on line $LINENO"; cleanup ; }' ERR
trap '{ echo "ufs_test.sh finished"; cleanup ; }' EXIT

# PATHRT - Path to regression tests directory
readonly PATHRT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
cd "${PATHRT}"

# make sure only one instance of ufs_test.sh is running
readonly LOCKDIR="${PATHRT}"/lock
if mkdir "${LOCKDIR}" ; then
  HOSTNAME=$(hostname)
  echo "${HOSTNAME} $$" > "${LOCKDIR}/PID"
else
  echo "Only one instance of ufs_test.sh can be running at a time"
  exit 1
fi

source detect_machine.sh # Note: this does not set ACCNR. The "if" block below does.
source rt_utils.sh
source module-setup.sh

CREATE_BASELINE=false
ROCOTO=false
ECFLOW=false
KEEP_RUNDIR=false
export skip_check_results=false
export delete_rundir=false
COMPILE_ONLY=false
RTPWD_NEW_BASELINE=false
TESTS_FILE='ufs_test.yaml'
NEW_BASELINES_FILE=''
RUN_SINGLE_TEST=false
ACCNR=${ACCNR:-""}
UFS_TEST_YAML="ufs_test.yaml"
export UFS_TEST_YAML

while getopts ":a:b:cl:mn:dwkreoh" opt; do
  case ${opt} in
    a)
	ACCNR=${OPTARG}
	;;
    b)
	NEW_BASELINES_FILE=${OPTARG}
	export NEW_BASELINES_FILE
	python -c "import create_yml; create_yml.update_testyaml_b()"
	UFS_TEST_YAML="ufs_test_temp.yaml"
	export UFS_TEST_YAML
	;;
    c)
	CREATE_BASELINE=true
	;;
    l)
	TESTS_FILE=${OPTARG}
	grep -q '[^[:space:]]' < "${TESTS_FILE}" ||  die "${TESTS_FILE} empty, exiting..."
	UFS_TEST_YAML=${TESTS_FILE}
	export UFS_TEST_YAML
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
	IFS=' ' read -r -a SINGLE_OPTS <<< "${OPTARG}"

	if [[ ${#SINGLE_OPTS[@]} != 2 ]]; then
            die 'The -n option needs <testname> AND <compiler> in quotes, i.e. -n "control_p8 intel"'
	fi

	SRT_NAME="${SINGLE_OPTS[0]}"
	SRT_COMPILER="${SINGLE_OPTS[1]}"

	if [[ "${SRT_COMPILER}" != "intel" ]] && [[ "${SRT_COMPILER}" != "gnu" ]]; then
            die "COMPILER MUST BE 'intel' OR 'gnu'"
	fi

	export SRT_NAME
	export SRT_COMPILER
	python -c "import create_yml; create_yml.update_testyaml_n()"
	UFS_TEST_YAML="ufs_test_temp.yaml"
	export UFS_TEST_YAML
	;;
    d)
	export delete_rundir=true
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
	die "Work-in-progress to support for ECFLOW. Please, use the ROCOTO workflow manamegment option (-r)"
	;;
    h)
	usage
	;;
    :)
	die "Option -${OPTARG} requires an argument."
	;;
    ?)
	die "Invalid option: -${OPTARG}"
	;;
    *)
        die "Invalid runtime options: no parameter included with argument-${OPTARG}"
	;;
  esac
done

#Check to error out if incompatible options are chosen together
[[ ${KEEP_RUNDIR} == true && ${delete_rundir} == true ]] && die "-k and -d options cannot be used at the same time"
[[ ${ECFLOW} == true && ${ROCOTO} == true ]] && die "-r and -e options cannot be used at the same time"
[[ ${CREATE_BASELINE} == true && ${RTPWD_NEW_BASELINE} == true ]] && die "-c and -m options cannot be used at the same time"

if [[ -z "${ACCNR}" ]]; then
  echo "Please use -a <account> to set group account to use on HPC"
  exit 1
fi

# Display the machine and account using the format detect_machine.sh used:
echo "Machine:  "${MACHINE_ID}"    Account: "${ACCNR}" "

check_machine=false
platforms=( hera orion hercules gaea jet derecho noaacloud s4 )
for name in "${platforms[@]}"
do
  if [[ ${MACHINE_ID} == "${name}" ]]; then
    check_machine=true
    break
  fi
done

if [[ ${check_machine} == true ]]; then
    source "${PATHRT}"/machine_config/machine_"${MACHINE_ID}".config
else
    die "*** Current support of ufs_test.sh only for hera orion hercules gaea jet derecho noaacloud s4 ! ***"
fi

shift $((OPTIND-1))
[[ $# -gt 1 ]] && usage

TEST_START_TIME="$(date '+%Y%m%d %T')"
export TEST_START_TIME

rm -f fail_test* fail_compile*

if [[ ${ROCOTO} == true ]]; then
  ROCOTO_XML="${PATHRT}"/rocoto_workflow.xml
  ROCOTO_STATE="${PATHRT}"/rocoto_workflow.state
  ROCOTO_DB="${PATHRT}"/rocoto_workflow.db
  rm -f "${ROCOTO_XML}" "${ROCOTO_DB}" "${ROCOTO_STATE}" *_lock.db
fi

[[ -f ${TESTS_FILE} ]] || die "${TESTS_FILE} does not exist"

export ROCOTO_SCHEDULER
export ACCNR
export ROCOTO_XML
export PATHRT
export ROCOTO
export ECFLOW
export MACHINE_ID
export RTPWD_NEW_BASELINE
export CREATE_BASELINE
export RTVERBOSE

export TESTS_FILE
export NEW_BASELINES_FILE
export RUN_SINGLE_TEST
export COMPILE_ONLY
export delete_rundir
export skip_check_results
export KEEP_RUNDIR  

if ! python -c "import create_xml; create_xml.xml_loop()"
then
  echo "*** experiment setup didn't run successfully! ***"
fi

##
## run regression test workflow (currently Rocoto or ecFlow are supported)
##
if [[ ${ROCOTO} == true ]]; then
    rocoto_run
fi

# IF -c AND -b; LINK VERIFIED BASELINES TO NEW_BASELINE
if [[ ${CREATE_BASELINE} == true && ${NEW_BASELINES_FILE} != '' ]]; then
    python -c "import ufs_test_utils; ufs_test_utils.link_new_baselines()"
fi

TEST_END_TIME="$(date '+%Y%m%d %T')"
export TEST_END_TIME

## Lets verify all tests were run and that they passed
python -c "import create_log; create_log.finish_log()"
