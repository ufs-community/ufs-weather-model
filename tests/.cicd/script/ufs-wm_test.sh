#!/usr/bin/env bash
#
# A unified test script for the WM application. This script is expected to
# test the WM application for all supported platforms. NOTE: At this time,
# this script is a placeholder for a more robust test framework.
#
set -e -u -x

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"

# Get repository root from Jenkins WORKSPACE variable if set, otherwise, set
# relative to script directory.
declare workspace
if [[ -n "${WORKSPACE}" ]]; then
    workspace="${WORKSPACE}"
else
    workspace="$(cd -- "${script_dir}/../.." && pwd)"
fi

# Normalize Parallel Works cluster platform value.
declare platform
platform="orion"
WM_COMPILER="intel"

# Test directories
we2e_experiment_base_dir="${workspace}/expt_dirs"
we2e_test_dir="${workspace}/tests/WE2E"
nco_dir="${workspace}/nco_dirs"

# Run the end-to-end tests.
if "${WM_WE2E_COMPREHENSIVE_TESTS}"; then
    test_type="comprehensive"
else
    test_type="fundamental"
fi

cd ${we2e_test_dir}
./setup_WE2E_tests.sh ${platform} ${WM_PROJECT} ${WM_COMPILER} ${test_type} ${we2e_experiment_base_dir} ${nco_dir}

# Allow the tests to start before checking for status.
# TODO: Create a parameter that sets the initial start delay.
sleep 300

# Progress file
progress_file="${workspace}/we2e_test_results-${platform}-${WM_COMPILER}.txt"

# Wait for all tests to complete.
while true; do

    # Check status of all experiments
    ./get_expts_status.sh expts_basedir="${we2e_experiment_base_dir}" \
         verbose="FALSE" | tee ${progress_file}

    # Exit loop only if there are not tests in progress
    set +e
    grep -q "Workflow status:  IN PROGRESS" ${progress_file}
    exit_code=$?
    set -e

    if [[ $exit_code -ne 0 ]]; then
       break
    fi

    # TODO: Create a paremeter that sets the poll frequency.
    sleep 60
done

# Allow we2e cron jobs time to complete and clean up themselves
# TODO: Create parameter that sets the interval for the we2e cron jobs; this
# value should be some factor of that interval to ensure the cron jobs execute
# before the workspace is cleaned up.
sleep 600

# Set exit code to number of failures
set +e
failures=$(grep "Workflow status:  FAILURE" ${progress_file} | wc -l)
set -e
exit ${failures}