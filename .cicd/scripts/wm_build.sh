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

cd "${UFS_MODEL_DIR}"
pwd
ls -l ./build.sh

BUILD_DIR=${BUILD_DIR:-${UFS_MODEL_DIR}/build}
TESTS_DIR=${TESTS_DIR:-${UFS_MODEL_DIR}/tests}
mkdir -p "${BUILD_DIR}"

(
	cd "${BUILD_DIR}"
	pwd
)

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
	sed -e "s|EPIC/spack-stack/spack-stack-1.5.0|spack-stack/spack-stack-1.5.1|g" -i modulefiles/ufs_noaacloud.intel.lua
fi
echo "machine_id=<${machine_id}>"

if [[ ${UFS_PLATFORM} = derecho ]] ; then
	export ACCNR=nral0032
else
	export ACCNR=epic
fi
echo "ACCNR=${ACCNR}"

export LMOD_SH_DBG_ON=0
echo "LMOD_VERSION=${LMOD_VERSION}"
if [[ ${UFS_PLATFORM} = gaea ]] ; then
	source /gpfs/f5/epic/scratch/role.epic/contrib/Lmod_init_C5.sh
	echo "LMOD_VERSION=${LMOD_VERSION}"
fi
set +x
module use ${PWD}/modulefiles >/dev/null 2>&1
module load ufs_${machine_id}.${UFS_COMPILER} || true
[[ ${UFS_PLATFORM} = gaea ]] && module load cmake/3.23.1 || true
module list

echo "Pipeline Building WM on ${UFS_PLATFORM} ${UFS_COMPILER} with Account=${ACCNR}."
export CMAKE_FLAGS="-DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16"
/usr/bin/time -p \
	-o ${workspace}/${UFS_PLATFORM}-${UFS_COMPILER}-time-wm_build.json \
	-f '{\n  "cpu": "%P"\n, "memMax": "%M"\n, "mem": {"text": "%X", "data": "%D", "swaps": "%W", "context": "%c", "waits": "%w"}\n, "pagefaults": {"major": "%F", "minor": "%R"}\n, "filesystem": {"inputs": "%I", "outputs": "%O"}\n, "time": {"real": "%e", "user": "%U", "sys": "%S"}\n}' \
	./build.sh | tee ${workspace}/${UFS_PLATFORM}-${UFS_COMPILER}-wm_build-log.txt
status=${PIPESTATUS[0]}
echo "Pipeline Completed WM build on ${UFS_PLATFORM} ${UFS_COMPILER}. status=${status}"

ls -l build/ufs_model
