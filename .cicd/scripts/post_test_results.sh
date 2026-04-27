#!/bin/bash
# RT - RegressionTest label
# BL - Baseline label

set -x
export machine=${1:-${NODE_NAME}}
label=${2:-${WM_TEST_LABEL}}
[[ -n "${label}" ]] || exit 1

export PATH=${PATH}:~/bin
echo "USER=${USER}"
echo "WORKSPACE=${WORKSPACE}"
export ACCNR=epic

export account="-a ${ACCNR}"

which jq

set -eu

SCRIPT_REALPATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPTS_DIR=$(dirname "${SCRIPT_REALPATH}")
UFS_MODEL_DIR=$(realpath "${SCRIPTS_DIR}/../..")
readonly UFS_MODEL_DIR
echo "UFS MODEL DIR: ${UFS_MODEL_DIR}"

export CC=${CC:-mpicc}
export CXX=${CXX:-mpicxx}
export FC=${FC:-mpif90}

BUILD_DIR=${BUILD_DIR:-${UFS_MODEL_DIR}/build}
TESTS_DIR=${TESTS_DIR:-${UFS_MODEL_DIR}/tests}

function post_test() {
	local machine=${1:-${NODE_NAME}}
	#local machine_id=${machine,,} # tolower
	#local machine_name_logs=$(echo "${machine}" | awk '{ print tolower($1) }')
	local label=${2:-"undef"}
	local WORKSPACE
	WORKSPACE="$(pwd)"
	GIT_URL=${GIT_URL:-"ufs-weather-model"}
	CHANGE_ID=${CHANGE_ID:-"develop"}

	git config user.email "ecc.platform@noaa.gov"
	git config user.name "epic-cicd-jenkins"

	git add tests/logs/RegressionTests_${machine,,}.log
	git status
        
        ##Check regression test logs results 
        if grep -q "Result: SUCCESS" ${UFS_MODEL_DIR}/tests/logs/RegressionTests_${machine,,}.log && grep -q "status=0" ${UFS_MODEL_DIR}/${machine,,}-status; then
           git commit -m "[AutoRT] ${machine} Job Completed Successfully.\n\n\n on-behalf-of @ufs-community <ecc.platform@noaa.gov>"
        else
           git commit --allow-empty -m "[AutoRT] ${machine} Job Failed! \n\n\n on-behalf-of @ufs-community <ecc.platform@noaa.gov>"
        fi

	SSH_ORIGIN=$(curl --silent "https://api.github.com/repos/ufs-community/ufs-weather-model/pulls/${CHANGE_ID}" | jq -r '.head.repo.ssh_url')
	git remote -v | grep -w sshorigin > /dev/null 2>&1 && git remote remove sshorigin > /dev/null 2>&1
	git remote add sshorigin ${SSH_ORIGIN} > /dev/null 2>&1 || return 0

	FORK_BRANCH=$(curl --silent "https://api.github.com/repos/ufs-community/ufs-weather-model/pulls/${CHANGE_ID}" | jq -r '.head.ref')
	git pull sshorigin ${FORK_BRANCH} || return 0
	git status
	git push sshorigin HEAD:${FORK_BRANCH} || return 0
}

pwd
post_test "${machine}" "${label}"
