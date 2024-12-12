#!/bin/bash
set -eu
export UFS_PLATFORM=${UFS_PLATFORM:-${NODE_NAME,,}}
export UFS_COMPILER=${UFS_COMPILER:-intel}

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

cd "${UFS_MODEL_DIR}"
echo "UFS_PLATFORM=<${UFS_PLATFORM}>"
echo "UFS_COMPILER=<${UFS_COMPILER}>"

pwd
echo "NODE_NAME=${NODE_NAME}"
echo "UFS_PLATFORM=${UFS_PLATFORM}"
echo "UFS_COMPILER=${UFS_COMPILER}"
workspace=$(pwd)
export workspace
machine=${NODE_NAME}
echo "machine=<${machine}>"
machine_id=${UFS_PLATFORM,,}
if [[ ${UFS_PLATFORM} =~ clusternoaa ]] ; then
	machine_id="noaacloud"
fi
echo "machine_id=<${machine_id}>"

/usr/bin/time -p \
	-o ${workspace}/${UFS_PLATFORM}-${UFS_COMPILER}-time-wm_init.json \
	-f '{\n  "cpu": "%P"\n, "memMax": "%M"\n, "mem": {"text": "%X", "data": "%D", "swaps": "%W", "context": "%c", "waits": "%w"}\n, "pagefaults": {"major": "%F", "minor": "%R"}\n, "filesystem": {"inputs": "%I", "outputs": "%O"}\n, "time": {"real": "%e", "user": "%U", "sys": "%S"}\n}' \
	pwd
