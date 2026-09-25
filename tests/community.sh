#!/usr/bin/env bash
# community.sh -
#    UFS Weather-model test driver for container-based runs (default) or
#    runs with natively installed software stack on community platforms ("-p")
#
# Reads a configuration file (e.g., community.conf)
#   that include header lines with container or community platform configuration
#   options, followed by one or more compile configurations/test blocks.  
#   All the compile jobs and tests run sequentially
#   NOTE: - no Rocoto or ECFlow workflow manager is involved
#         - no baseline test created or comparison made
# Test reported as failed if the compile or test fails, and a fail flag file is created.
# 
#   For each compile configuration:
#     1. Compile the model inside the software container or on the community platform
#     2. Run every listed test case sequentially, waiting for each test to finish
#
# Users many need to adapt modulefiles before running:
#
# Modulefiles in the ./modulefiles directory:
#   container option (MACHINE_ID=container):
#     ufs_container.runtime.lua   — host-side runtime module 
#     NOTE: no need to change a modulefile for inside-container build
#        and run environment, already included in the container image: 
#        modulefiles/ufs_container.<compiler>.lua
#
#   community platform with a native software stack ("-p" flag, MACHINE_ID=<COMMUNITY>):
#     ufs_<COMMUNITY>.<compiler>.lua     — community platform module
#
# Job scheduler job templates for Slurm or PBS (if job scheduler is used) in
#    ./tests/fv3_conf/ directory:
#
#   container option (MACHINE_ID=container):
#     compile_slurm.IN_container  or  compile_qsub.IN_container
#     fv3_slurm.IN_container      or  fv3_qsub.IN_container
#
#   community platform option ("-p" flag, MACHINE_ID=<COMMUNITY>):
#     compile_slurm.IN_<COMMUNITY>  or  compile_qsub.IN_<COMMUNITY>
#     fv3_slurm.IN_<COMMUNITY>      or  fv3_qsub.IN_<COMMUNITY>
#

set -uo pipefail

###############################################################################
# Locate this script; set PATHRT and PATHTR
###############################################################################

SCRIPT_REALPATH=$(realpath "${BASH_SOURCE[0]}")
PATHRT=$(dirname "${SCRIPT_REALPATH}")
PATHTR=$(cd "${PATHRT}/.." && pwd)

readonly PATHRT PATHTR

###############################################################################
# Helpers
###############################################################################

trim() {
    local s="$1"
    s="${s#"${s%%[![:space:]]*}"}"
    s="${s%"${s##*[![:space:]]}"}"
    echo "${s}"
}

usage() {
    echo ""
    echo "Usage: $(basename "$0") [options] <conf_file>"
    echo ""
    echo "  <conf_file>   Required. Path to the configuration file (e.g. community.conf)."
    echo "                Defines the platform, scheduler, data paths, and the list of"
    echo "                compile configurations and test cases to run."
    echo ""
    echo "Options:"
    echo "  -p            platform mode: run without a container with a natively installed software stack"
    echo "  -d            delete run directories after each test completes"
    echo "  -h            display this help"
    echo "  -n <name>     run only the single test named <name>;"
    echo "                <name> must match a test case listed in <conf_file>"
    echo "  -o            compile only, skip all tests"
    echo "  -v            verbose output (set -x)"
    echo ""
}

###############################################################################
# Flag defaults
###############################################################################

export delete_rundir=false
COMMUNITY_PLATFORM=false
COMPILE_ONLY=false
RUN_SINGLE_TEST=false
SRT_NAME=''
export skip_check_results=true
export RTVERBOSE=false
TESTS_FILE=''

###############################################################################
# Parse command-line options
###############################################################################

# getopts handles only single-character flags; scan for --help before it runs.
for _arg in "$@"; do
    [[ "${_arg}" == "--help" ]] && { usage; exit 0; }
    [[ "${_arg}" == "--"    ]] && break
done
unset _arg

while getopts ":pdhn:ov" opt; do
    case ${opt} in
        p) COMMUNITY_PLATFORM=true ;;
        d) export delete_rundir=true ;;
        h) usage; exit 0 ;;
        n) RUN_SINGLE_TEST=true; SRT_NAME=${OPTARG} ;;
        o) COMPILE_ONLY=true ;;
        v) export RTVERBOSE=true ;;
        \?) usage; echo "ERROR: invalid option -${OPTARG}" >&2; exit 1 ;;
        :)  usage; echo "ERROR: option -${OPTARG} requires an argument" >&2; exit 1 ;;
        *)  usage; echo "ERROR: unknown arguments" >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

# Conf file is a required positional argument.
if [[ -n "${1:-}" ]]; then
    TESTS_FILE="$1"
else
    usage
    echo "ERROR: configuration file is required." >&2
    echo "       Provide the path to a community.conf-style file as the last argument." >&2
    exit 1
fi

if [[ "${RTVERBOSE}" == true ]]; then
    set -x
fi

###############################################################################
# Configuration file
###############################################################################

input_file="${TESTS_FILE}"

if [[ ! -f "${input_file}" ]]; then
    echo "ERROR: configuration file not found: ${input_file}" >&2
    exit 1
fi

###############################################################################
# Parse community.conf
#
# Header line 1: MACHINE_ID | RT_COMPILER | CONTAINER_IMG | BIND_DIRS
# Header line 2: TPN | SCHEDULER | ACCNR | PARTITION | QUEUE | MPI_LAUNCH
# Header line 3: RUNDIR_ROOT
# Header line 4: INPUTDATA_ROOT | INPUTDATA_ROOT_WW3 | INPUTDATA_LM4 | INPUTDATA_GFSv17opn
#
# Configuration line (after headers): compile_id | MAKE_OPT   (contains '|')
# Case line         (after headers):  test_case_name           (no '|')
###############################################################################

declare -a compile_ids=()
declare -a make_opts=()
declare -A tests_by_compile_id

header_lines_read=0
current_compile_id=""

while IFS= read -r line || [[ -n "${line}" ]]; do

    line="${line%$'\r'}"

    # Skip blank lines; a blank line also closes the current configuration block.
    if [[ -z "${line//[[:space:]]/}" ]]; then
        current_compile_id=""
        continue
    fi

    # Skip comment lines.
    if [[ "${line}" =~ ^[[:space:]]*# ]]; then
        continue
    fi

    # ------------------------------------------------------------------
    # Header line 1: MACHINE_ID | RT_COMPILER | CONTAINER_IMG | BIND_DIRS
    # ------------------------------------------------------------------
    if [[ ${header_lines_read} -eq 0 ]]; then
        IFS='|' read -r f1 f2 f3 f4 _rest <<< "${line}"
        MACHINE_ID=$(trim "${f1:-container}")
        RT_COMPILER=$(trim "${f2:-}")
        CONTAINER_IMG=$(trim "${f3:-}")
        CONTAINER_BIND=$(trim "${f4:-}")
        header_lines_read=1
        if [[ "${RTVERBOSE}" == true ]]; then
            echo "MACHINE_ID=${MACHINE_ID}"
            echo "RT_COMPILER=${RT_COMPILER}"
            echo "CONTAINER_IMG=${CONTAINER_IMG}"
            echo "BIND_DIRS=${CONTAINER_BIND}"
        fi
        continue
    fi

    # ------------------------------------------------------------------
    # Header line 2: TPN | SCHEDULER | ACCNR | PARTITION | QUEUE | MPI_LAUNCH
    # ------------------------------------------------------------------
    if [[ ${header_lines_read} -eq 1 ]]; then
        IFS='|' read -r f1 f2 f3 f4 f5 f6 _rest <<< "${line}"
        TPN=$(trim "${f1:-40}")
        SCHEDULER=$(trim "${f2:-}")
        ACCNR=$(trim "${f3:-}")
        PARTITION=$(trim "${f4:-}")
        QUEUE=$(trim "${f5:-}")
        MPI_LAUNCH=$(trim "${f6:-mpirun}")
        ROCOTO=false
        ECFLOW=false
        header_lines_read=2
        if [[ "${RTVERBOSE}" == true ]]; then
            echo "TPN=${TPN}  SCHEDULER=${SCHEDULER}"
            echo "ACCNR=${ACCNR}  PARTITION=${PARTITION}  QUEUE=${QUEUE}"
            [[ "${SCHEDULER:-none}" == none ]] && echo "MPI_LAUNCH=${MPI_LAUNCH}"
        fi
        continue
    fi

    # ------------------------------------------------------------------
    # Header line 3: RUNDIR_ROOT
    # ------------------------------------------------------------------
    if [[ ${header_lines_read} -eq 2 ]]; then
        RUNDIR_ROOT=$(trim "${line}")
        header_lines_read=3
        [[ "${RTVERBOSE}" == true ]] && echo "RUNDIR_ROOT=${RUNDIR_ROOT}"
        continue
    fi

    # ------------------------------------------------------------------
    # Header line 4: INPUTDATA_ROOT | INPUTDATA_ROOT_WW3 | INPUTDATA_LM4 | INPUTDATA_GFSv17opn
    # ------------------------------------------------------------------
    if [[ ${header_lines_read} -eq 3 ]]; then
        IFS='|' read -r f1 f2 f3 f4 _rest <<< "${line}"
        INPUTDATA_ROOT=$(trim "${f1:-}")
        INPUTDATA_ROOT_WW3=$(trim "${f2:-}")
        INPUTDATA_LM4=$(trim "${f3:-}")
        INPUTDATA_GFSv17opn=$(trim "${f4:-}")
        header_lines_read=4
        [[ "${RTVERBOSE}" == true ]] && echo "INPUTDATA_ROOT=${INPUTDATA_ROOT}"
        [[ "${RTVERBOSE}" == true ]] && echo "INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT_WW3}"
        [[ "${RTVERBOSE}" == true ]] && echo "INPUTDATA_LM4=${INPUTDATA_LM4}"
        [[ "${RTVERBOSE}" == true ]] && echo "INPUTDATA_GFSv17opn=${INPUTDATA_GFSv17opn}"
        continue
    fi

    # ------------------------------------------------------------------
    # Configuration line: compile_id | MAKE_OPT   (contains '|')
    # Case line:          test_case_name           (no '|')
    # ------------------------------------------------------------------
    if [[ "${line}" == *"|"* ]]; then
        IFS='|' read -r f1 f2 _rest <<< "${line}"
        current_compile_id=$(trim "${f1:-}")
        make_opt="${f2:-}"
        make_opt="${make_opt#"${make_opt%%[![:space:]]*}"}"  # ltrim only
        compile_ids+=( "${current_compile_id}" )
        make_opts+=( "${make_opt}" )
        tests_by_compile_id["${current_compile_id}"]=""
    else
        test_case=$(trim "${line}")
        if [[ -z "${current_compile_id}" ]]; then
            echo "WARNING: case '${test_case}' has no active configuration — skipping" >&2
            continue
        fi
        tests_by_compile_id["${current_compile_id}"]+="${test_case}"$'\n'
    fi

done < "${input_file}"

if [[ ${header_lines_read} -lt 4 ]]; then
    echo "ERROR: fewer than 4 header lines found in ${input_file}" >&2
    exit 1
fi

if [[ ${#compile_ids[@]} -eq 0 ]]; then
    echo "ERROR: no compile configuration sections found in ${input_file}" >&2
    exit 1
fi

echo "====================================================================="
echo "Configuration summary:"
echo "  MACHINE_ID:     ${MACHINE_ID}"
echo "  RT_COMPILER:    ${RT_COMPILER}"
[[ -n "${CONTAINER_IMG}" ]]       && echo "  CONTAINER_IMG:  ${CONTAINER_IMG}"
[[ -n "${CONTAINER_BIND}" ]]      && echo "  BIND_DIRS:      ${CONTAINER_BIND}"
echo "  TPN:            ${TPN}"
[[ -n "${SCHEDULER}" ]]           && echo "  SCHEDULER:      ${SCHEDULER}"
[[ -n "${ACCNR}" ]]               && echo "  ACCNR:          ${ACCNR}"
[[ -n "${PARTITION}" ]]           && echo "  PARTITION:      ${PARTITION}"
[[ -n "${QUEUE}" ]]               && echo "  QUEUE:          ${QUEUE}"
[[ -n "${MPI_LAUNCH}" ]]          && echo "  MPI_LAUNCH:     ${MPI_LAUNCH}"
echo "  RUNDIR_ROOT:    ${RUNDIR_ROOT}"
echo "  INPUTDATA_ROOT: ${INPUTDATA_ROOT}"
[[ -n "${INPUTDATA_ROOT_WW3}" ]]  && echo "  INPUTDATA_ROOT_WW3:  ${INPUTDATA_ROOT_WW3}"
[[ -n "${INPUTDATA_LM4}" ]]       && echo "  INPUTDATA_LM4:       ${INPUTDATA_LM4}"
[[ -n "${INPUTDATA_GFSv17opn}" ]] && echo "  INPUTDATA_GFSv17opn: ${INPUTDATA_GFSv17opn}"
echo "====================================================================="
echo ""
echo "Found ${#compile_ids[@]} compile configuration(s):"
for i in "${!compile_ids[@]}"; do
    echo "  $((i+1)): ${compile_ids[$i]}  [${make_opts[$i]}]"
    while IFS= read -r _t; do
        [[ -z "${_t}" ]] && continue
        echo "       - ${_t}"
    done <<< "${tests_by_compile_id[${compile_ids[$i]}]}"
done
echo ""

###############################################################################
# Derived paths (after conf is parsed)
###############################################################################

# Derive INPUTDATA_ROOT sub-paths (same convention as rt.sh).
# All three can be overridden by setting the variable before invoking this script.
INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT_WW3:-${INPUTDATA_ROOT}/WW3_input_data_20250807}
INPUTDATA_LM4=${INPUTDATA_LM4:-${INPUTDATA_ROOT}/LM4_input_data}
INPUTDATA_GFSv17opn=${INPUTDATA_GFSv17opn:-${INPUTDATA_ROOT}/GFSv17opn_input_data}

###############################################################################
# Validate prerequisites
###############################################################################

if [[ "${MACHINE_ID}" == container && "${COMMUNITY_PLATFORM}" == true ]]; then
    echo "ERROR: MACHINE_ID=container and -p (COMMUNITY_PLATFORM) are mutually exclusive." >&2
    echo "       Use MACHINE_ID=container without -p for Singularity/Apptainer container runs." >&2
    echo "       Use a custom MACHINE_ID with -p for native-stack community platform runs." >&2
    exit 1
fi

if [[ "${MACHINE_ID}" == container ]]; then
    if [[ ! -f "${PATHTR}/modulefiles/ufs_container.${RT_COMPILER}.lua" ]]; then
        echo "ERROR: modulefiles/ufs_container.${RT_COMPILER}.lua not found under ${PATHTR}" >&2
        echo "       Provide this modulefile for the inside-container build/run environment." >&2
        exit 1
    fi
    if [[ ! -f "${PATHTR}/modulefiles/ufs_container.runtime.lua" ]]; then
        echo "ERROR: modulefiles/ufs_container.runtime.lua not found under ${PATHTR}" >&2
        echo "       Provide this host-side runtime modulefile before running." >&2
        exit 1
    fi
fi

if [[ "${COMMUNITY_PLATFORM}" = true ]]; then
    if [[ ! -f "${PATHTR}/modulefiles/ufs_${MACHINE_ID}.${RT_COMPILER}.lua" ]]; then
        echo "ERROR: modulefiles/ufs_${MACHINE_ID}.${RT_COMPILER}.lua not found under ${PATHTR}" >&2
        echo "       Provide this modulefile for the community platform build/run environment." >&2
        exit 1
    fi
fi

if [[ "${MACHINE_ID}" == container ]]; then
    if [[ ! -f "${CONTAINER_IMG}" ]]; then
        echo "ERROR: container image not found: ${CONTAINER_IMG}" >&2
        exit 1
    fi
elif [[ "${COMMUNITY_PLATFORM}" != true ]]; then
    echo "ERROR: MACHINE_ID='${MACHINE_ID}' is not 'container' but -p was not specified." >&2
    echo "       For native-stack community platform runs with a custom MACHINE_ID," >&2
    echo "       re-run with the -p option:  $(basename "$0") -p ... ${input_file}" >&2
    exit 1
fi

###############################################################################
# Set up run directory and logging
###############################################################################

LOG_DIR="${RUNDIR_ROOT}/logs"
mkdir -p "${LOG_DIR}"

# Create a convenience symlink $PATHRT/run_dir -> RUNDIR_ROOT, unless the user
# already set RUNDIR_ROOT to that exact path.
if [[ "${RUNDIR_ROOT}" != "${PATHRT}/run_dir" ]]; then
    ln -sfn "${RUNDIR_ROOT}" "${PATHRT}/run_dir"
fi

[[ "${RTVERBOSE}" == true ]] && echo "LOG_DIR=${LOG_DIR}"

# When -d is set, run_test.sh reads keep_tests.tmp to decide which dirs to keep.
# An empty file means every test dir gets removed.
if [[ "${delete_rundir}" == true ]]; then
    : > "${PATHRT}/keep_tests.tmp"
fi

###############################################################################
# Export variables used by run_compile.sh, run_test.sh, and atparse
###############################################################################

export MACHINE_ID
export COMMUNITY_PLATFORM
export RT_COMPILER TPN CONTAINER_IMG CONTAINER_BIND
export SCHEDULER ROCOTO ECFLOW ACCNR PARTITION QUEUE MPI_LAUNCH
export INPUTDATA_ROOT INPUTDATA_ROOT_WW3 INPUTDATA_LM4 INPUTDATA_GFSv17opn
export PATHRT PATHTR RUNDIR_ROOT LOG_DIR

###############################################################################
# Source atparse (used by run_compile.sh/run_test.sh; sourced here so the
# script can verify it is present before starting any work)
###############################################################################

source "${PATHRT}/atparse.bash"

###############################################################################
# Print active flags and confirm before starting
###############################################################################

if [[ "${RTVERBOSE}" == true ]]; then
    echo "Active flags:"
    [[ "${delete_rundir}" == true ]]   && echo "  -d  delete run directories after each test"
    [[ "${RUN_SINGLE_TEST}" == true ]] && echo "  -n  single test: ${SRT_NAME}"
    [[ "${COMPILE_ONLY}" == true ]]    && echo "  -o  compile only, skip tests"
    echo "  -v  verbose output"
    echo ""
    read -r -n 1 -s -p "Configuration above — press any key to start, or Ctrl-C to abort... "
    echo ""
    echo ""
fi

###############################################################################
# Tracking
###############################################################################

declare -a failed_compiles=()
declare -a failed_tests=()
declare -a skipped_tests=()

###############################################################################
# Main loop: compile then run all tests for each configuration, sequentially
###############################################################################

for i in "${!compile_ids[@]}"; do

    export COMPILE_ID="${compile_ids[$i]}"
    export MAKE_OPT="${make_opts[$i]}"
    export BUILD_NAME="fv3_${COMPILE_ID}"
    export JBNME="compile_${COMPILE_ID}"

    # When running a single test, skip compiles that don't include it.
    if [[ "${RUN_SINGLE_TEST}" == true ]]; then
        if [[ $'\n'"${tests_by_compile_id[$COMPILE_ID]}" != *$'\n'"${SRT_NAME}"$'\n'* ]]; then
            echo "Skipping compile ${COMPILE_ID} (no matching test for -n ${SRT_NAME})"
            continue
        fi
    fi

    echo "====================================================================="
    echo "Compile $((i+1)) of ${#compile_ids[@]}: COMPILE_ID=${COMPILE_ID}"
    echo "MAKE_OPT=${MAKE_OPT}"
    echo "====================================================================="

    # ------------------------------------------------------------------
    # Skip compile if the executable already exists in PATHRT.
    # ------------------------------------------------------------------
    if [[ -f "${PATHRT}/fv3_${COMPILE_ID}.exe" ]]; then
        echo "Found existing ${PATHRT}/fv3_${COMPILE_ID}.exe — skipping compile"
        echo ""
    else
        # ------------------------------------------------------------------
        # Write the compile environment file (sourced by run_compile.sh).
        # env file is sourced twice in run_compile.sh (before and after
        # default_vars.sh), so variables here survive the defaults reset.
        # ------------------------------------------------------------------
        mkdir -p "${RUNDIR_ROOT}"
        cat > "${RUNDIR_ROOT}/${JBNME}.env" << ENV_EOF
export MACHINE_ID=${MACHINE_ID}
export PATHTR=${PATHTR}
export PATHRT=${PATHRT}
export RUNDIR_ROOT=${RUNDIR_ROOT}
export LOG_DIR=${LOG_DIR}
export RT_COMPILER=${RT_COMPILER}
export SCHEDULER=${SCHEDULER}
export ROCOTO=false
export ECFLOW=false
export ACCNR=${ACCNR}
export PARTITION=${PARTITION}
export QUEUE=${QUEUE}
export TPN=${TPN}
export CONTAINER_IMG=${CONTAINER_IMG}
export CONTAINER_BIND=${CONTAINER_BIND}
export COMMUNITY_PLATFORM=${COMMUNITY_PLATFORM}
export MPI_LAUNCH=${MPI_LAUNCH}
export RTVERBOSE=${RTVERBOSE}
export WLCLK=120
ENV_EOF

        # ------------------------------------------------------------------
        # Submit compile job and wait (submit_and_wait inside run_compile.sh
        # blocks until the scheduler job completes).
        # run_compile.sh exits 0 even on failure in non-rocoto mode; check
        # the fail file instead.
        # ------------------------------------------------------------------
        "${PATHRT}/run_compile.sh" \
            "${PATHRT}" "${RUNDIR_ROOT}" "${MAKE_OPT}" "${COMPILE_ID}" || true

        if [[ -f "${PATHRT}/fail_${JBNME}" ]]; then
            echo "ERROR: compile failed for COMPILE_ID=${COMPILE_ID}" >&2
            echo "       See ${PATHRT}/fail_${JBNME} and ${LOG_DIR}/${JBNME}.log" >&2
            failed_compiles+=("${COMPILE_ID}")
            continue
        fi

        if [[ ! -f "${PATHRT}/fv3_${COMPILE_ID}.exe" ]]; then
            echo "WARNING: fv3_${COMPILE_ID}.exe not found after compile — skipping tests" >&2
            skipped_tests+=("${COMPILE_ID}:ALL:missing_executable")
            continue
        fi

        echo "Compile ${COMPILE_ID} succeeded — ${PATHRT}/fv3_${COMPILE_ID}.exe"
        echo ""
    fi

    # -o: compile only, do not run tests.
    if [[ "${COMPILE_ONLY}" == true ]]; then
        echo "(-o) Skipping tests for COMPILE_ID=${COMPILE_ID}"
        continue
    fi

    # ------------------------------------------------------------------
    # Run test cases for this configuration, one at a time
    # ------------------------------------------------------------------
    mapfile -t test_cases <<< "${tests_by_compile_id[$COMPILE_ID]}"

    if [[ ${#test_cases[@]} -eq 0 ]]; then
        echo "WARNING: no test cases listed for COMPILE_ID=${COMPILE_ID}" >&2
        skipped_tests+=("${COMPILE_ID}:ALL:no_tests")
        continue
    fi

    for TEST_NAME in "${test_cases[@]}"; do
        [[ -z "${TEST_NAME}" ]] && continue

        # -n: run only the named test.
        if [[ "${RUN_SINGLE_TEST}" == true && "${TEST_NAME}" != "${SRT_NAME}" ]]; then
            continue
        fi

        export TEST_NAME
        export TEST_ID="${TEST_NAME}_${RT_COMPILER}"

        echo "---------------------------------------------------------------------"
        echo "Test: ${TEST_NAME}  (COMPILE_ID=${COMPILE_ID})"
        echo "---------------------------------------------------------------------"

        # Write the test environment file (sourced by run_test.sh).
        cat > "${RUNDIR_ROOT}/run_test_${TEST_ID}.env" << ENV_EOF
export MACHINE_ID=${MACHINE_ID}
export PATHTR=${PATHTR}
export PATHRT=${PATHRT}
export RUNDIR_ROOT=${RUNDIR_ROOT}
export LOG_DIR=${LOG_DIR}
export RT_COMPILER=${RT_COMPILER}
export SCHEDULER=${SCHEDULER}
export ROCOTO=false
export ECFLOW=false
export ACCNR=${ACCNR}
export PARTITION=${PARTITION}
export QUEUE=${QUEUE}
export TPN=${TPN}
export CONTAINER_IMG=${CONTAINER_IMG}
export CONTAINER_BIND=${CONTAINER_BIND}
export COMMUNITY_PLATFORM=${COMMUNITY_PLATFORM}
export MPI_LAUNCH=${MPI_LAUNCH}
export INPUTDATA_ROOT=${INPUTDATA_ROOT}
export INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT_WW3}
export INPUTDATA_LM4=${INPUTDATA_LM4}
export INPUTDATA_GFSv17opn=${INPUTDATA_GFSv17opn}
export CNTL_DIR=${TEST_NAME}
export RT_SUFFIX=
export BL_SUFFIX=
export CREATE_BASELINE=false
export skip_check_results=true
export delete_rundir=${delete_rundir}
export RTVERBOSE=${RTVERBOSE}
export WLCLK=60
ENV_EOF

        rm -f "${PATHRT}/fail_test_${TEST_ID}"

        "${PATHRT}/run_test.sh" \
            "${PATHRT}" "${RUNDIR_ROOT}" "${TEST_NAME}" "${TEST_ID}" "${COMPILE_ID}" || true

        if [[ -f "${PATHRT}/fail_test_${TEST_ID}" ]]; then
            echo "FAIL: ${TEST_NAME}" >&2
            failed_tests+=("${COMPILE_ID}:${TEST_NAME}")
        else
            echo "PASS: ${TEST_NAME}"
        fi

    done

    echo ""

done

###############################################################################
# Cleanup
###############################################################################

[[ "${delete_rundir}" == true ]] && rm -f "${PATHRT}/keep_tests.tmp"

###############################################################################
# Summary
###############################################################################

echo ""
echo "====================================================================="
echo "community.sh finished"
echo "====================================================================="

if [[ ${#skipped_tests[@]} -gt 0 ]]; then
    echo ""
    echo "Skipped:"
    for s in "${skipped_tests[@]}"; do
        IFS=':' read -r sc st sr <<< "${s}"
        echo "  COMPILE_ID=${sc}  TEST=${st}  reason=${sr}"
    done
fi

if [[ ${#failed_compiles[@]} -gt 0 ]]; then
    echo ""
    echo "Failed compiles:"
    for fc in "${failed_compiles[@]}"; do
        echo "  ${fc}"
    done
fi

if [[ ${#failed_tests[@]} -gt 0 ]]; then
    echo ""
    echo "Failed tests:"
    for ft in "${failed_tests[@]}"; do
        IFS=':' read -r fc ft_name <<< "${ft}"
        echo "  COMPILE_ID=${fc}  TEST=${ft_name}"
    done
    exit 1
fi

if [[ ${#failed_compiles[@]} -gt 0 ]]; then
    exit 1
fi

echo ""
echo "All tests completed successfully."
exit 0
