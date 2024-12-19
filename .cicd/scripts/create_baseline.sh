#!/bin/bash -x

export PATH=${PATH}:~/bin
echo "USER=${USER}"
echo "WORKSPACE=${WORKSPACE}"
export ACCNR=epic

export account="-a ${ACCNR}"
export workflow="-e"
	#[[ ${UFS_PLATFORM} =  jet         ]] && workflow="-r"
	#[[ ${UFS_PLATFORM} =  hera        ]] && workflow="-r"
	#[[ ${UFS_PLATFORM} =~ clusternoaa ]] && workflow=""

export opt="-l"
export suite="rt.conf"
	[[ -n ${WM_OPERATIONAL_TESTS}                 ]] && opt="-n" && suite="${WM_OPERATIONAL_TESTS} ${UFS_COMPILER}" || return 0
	[[    ${WM_OPERATIONAL_TESTS} = default       ]] && opt="-n" && suite="control_p8 ${UFS_COMPILER}"
	[[    ${WM_OPERATIONAL_TESTS} = comprehensive ]] && opt="-l" && suite="rt.conf"
	[[    ${WM_OPERATIONAL_TESTS} = rt.conf       ]] && opt="-l" && suite="rt.conf"
	[[   "${suite}"               = rt.conf       ]] && opt="-l"

set -eu

export machine=${NODE_NAME}

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

pwd
ls -al .cicd/*
ls -al ${TESTS_DIR}/rt.sh

function create_baseline() {
	local machine=${1:-${NODE_NAME}}
	local machine_id=${machine,,} # tolower
	local WORKSPACE
	WORKSPACE="$(pwd)"
	local status=0

	git submodule update --init --recursive
	pwd
	ls -al .cicd/*
	cd tests
		pwd

		[[ ${UFS_PLATFORM} =~ clusternoaa ]] && echo "export BL_DATE=20240426" > bl_date.conf || cat bl_date.conf

		mkdir -p logs/
		BL_DATE=$(cut -d '=' -f2 bl_date.conf)
		export BL_DATE

		if [[ ${machine} =~ "Jet" ]]
		then
		    echo "Creating baselines on ${machine}"
		    export dprefix=/lfs5/NAGAPE/${ACCNR}/${USER}
		    sed 's|/lfs4/HFIP/${ACCNR}/${USER}|/lfs4/HFIP/hfv3gfs/${USER}|g' -i rt.sh
		    sed 's|/lfs5/HFIP/${ACCNR}/${USER}|/lfs5/NAGAPE/${ACCNR}/${USER}|g' -i rt.sh
		    local workflow="-r"
		    ./rt.sh -a "${ACCNR}" -c "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		elif [[ ${machine} =~ "Hercules" ]]
		then
		    echo "Creating baselines on ${machine}"
		    export dprefix=/work2/noaa/${ACCNR}/${USER}
		    sed "s|/noaa/stmp/|/noaa/${ACCNR}/stmp/|g" -i rt.sh
		    export ACCNR=epic
		    ./rt.sh -a "${ACCNR}" -c "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    export DISKNM=/work/noaa/epic/hercules/UFS-WM_RT
		    cd ${DISKNM}/NEMSfv3gfs/
		    mkdir -p develop-${BL_DATE}
		    cd /work2/noaa/epic/stmp/role-epic/stmp/role-epic/FV3_RT
		    ls -l REGRESSION_TEST/.
		    rsync -a --no-t REGRESSION_TEST/ "${DISKNM}/NEMSfv3gfs/develop-${BL_DATE}" || echo "#### Warning! rsync $(pwd)/REGRESSION_TEST/ incomplete."
		    cd ${DISKNM}/NEMSfv3gfs/
		    ./adjust_permissions.sh hercules "develop-${BL_DATE}" || :
		    chgrp noaa-hpc "develop-${BL_DATE}" || :
		    cd ${WORKSPACE}/tests
		    ./rt.sh -a "${ACCNR}" "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    cd logs/
		    cp "RegressionTests_${machine_id}.log" "$(dirname "${WORKSPACE}")" #/work/noaa/epic/role-epic/jenkins/workspace
		    git remote -v
		    git fetch --no-recurse-submodules origin
		    git reset FETCH_HEAD --hard
		    cd .. && cd .. && cd ..
		    pwd
		    cp "$(dirname "${WORKSPACE}")/RegressionTests_${machine_id}.log" "${WORKSPACE}/tests/logs/"
		    cd ${WORKSPACE}/tests/
		elif [[ ${machine} =~ "Orion" ]]
		then
		    cd ..
		    #module load git/2.28.0
			git --version
		    git submodule update --init --recursive
		    cd tests
		    echo "Creating baselines on ${machine}"
		    export dprefix=/work2/noaa/${ACCNR}/${USER}
		    sed "s|/noaa/stmp/|/noaa/${ACCNR}/stmp/|g" -i rt.sh
		    export ACCNR=epic
		    ./rt.sh -a "${ACCNR}" -c "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    export DISKNM=/work/noaa/epic/UFS-WM_RT
		    cd ${DISKNM}/NEMSfv3gfs/
		    mkdir -p develop-${BL_DATE}
		    cd  /work/noaa/epic/stmp/role-epic/stmp/role-epic/FV3_RT/
		    ls -l REGRESSION_TEST/.
		    rsync -a --no-t REGRESSION_TEST/ "${DISKNM}/NEMSfv3gfs/develop-${BL_DATE}" || echo "#### Warning! rsync $(pwd)/REGRESSION_TEST/ incomplete."
		    cd ${DISKNM}/NEMSfv3gfs/
		    ./adjust_permissions.sh orion "develop-${BL_DATE}" || :
		    chgrp noaa-hpc "develop-${BL_DATE}" || :
		    cd ${WORKSPACE}/tests
		    ./rt.sh -a "${ACCNR}" "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    cd logs/
		    cp "RegressionTests_${machine_id}.log" "$(dirname "${WORKSPACE}")" #/work/noaa/epic/role-epic/jenkins/workspace
		    git remote -v
		    git fetch --no-recurse-submodules origin
		    git reset FETCH_HEAD --hard
		    cd .. && cd .. && cd ..
		    pwd
		    cp "$(dirname "${WORKSPACE}")/RegressionTests_${machine_id}.log" "${WORKSPACE}/tests/logs/"
		    cd ${WORKSPACE}/tests/
		elif [[ ${machine} =~ "Gaea" ]]
		then
		    echo "Creating baselines on ${machine}"
		    ./rt.sh -a "${ACCNR}" -c "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    unset LD_LIBRARY_PATH
		    export DISKNM=/gpfs/f5/epic/world-shared/UFS-WM_RT
		    cd ${DISKNM}/NEMSfv3gfs/
		    mkdir -p develop-${BL_DATE}
		    cd /gpfs/f5/epic/scratch/role.epic/FV3_RT
		    ls -l REGRESSION_TEST/.
		    rsync -a --no-t REGRESSION_TEST/ "${DISKNM}/NEMSfv3gfs/develop-${BL_DATE}" || echo "#### Warning! rsync $(pwd)/REGRESSION_TEST/ incomplete."
		    cd ${DISKNM}/NEMSfv3gfs/
		    chgrp ncep "develop-${BL_DATE}" || :
		    cd ${WORKSPACE}/tests
		    ./rt.sh -a "${ACCNR}" "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    cd logs/
		    cp "RegressionTests_${machine_id}.log" "$(dirname "${WORKSPACE}")" #/gpfs/f5/epic/scratch/role.epic/jenkins/workspace
		    git remote -v
		    git fetch --no-recurse-submodules origin
		    git reset FETCH_HEAD --hard
		    cd .. && cd .. && cd ..
		    pwd
		    cp "$(dirname "${WORKSPACE}")/RegressionTests_${machine_id}.log" "${WORKSPACE}/tests/logs/"
		    cd ${WORKSPACE}/tests/
		elif [[ ${machine} =~ "Hera" ]]
		then
		    echo "Creating baselines on ${machine}"
		    export ACCNR=epic
		    sed "s|QUEUE=batch|QUEUE=windfall|g" -i rt.sh
		    local workflow="-r"
		    ./rt.sh -a "${ACCNR}" -c "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    export DISKNM=/scratch2/NAGAPE/epic/UFS-WM_RT
		    cd ${DISKNM}/NEMSfv3gfs/
		    mkdir -p develop-${BL_DATE}
		    cd  /scratch1/NCEPDEV/stmp4/role.epic/FV3_RT
		    ls -l REGRESSION_TEST/.
		    rsync -a --no-t REGRESSION_TEST/ "${DISKNM}/NEMSfv3gfs/develop-${BL_DATE}" || echo "#### Warning! rsync $(pwd)/REGRESSION_TEST/ incomplete."
		    cd ${WORKSPACE}/tests
		    ./rt.sh -a "${ACCNR}" "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    cd logs/
		    cp "RegressionTests_${machine_id}.log" "$(dirname "${WORKSPACE}")" #/scratch2/NAGAPE/epic/role.epic/jenkins/workspace
		    git remote -v
		    git fetch --no-recurse-submodules origin
		    git reset FETCH_HEAD --hard
		    cd .. && cd .. && cd ..
		    pwd
		    cp "$(dirname "${WORKSPACE}")/RegressionTests_${machine_id}.log" "${WORKSPACE}/tests/logs/"
		    cd ${WORKSPACE}/tests/
		elif [[ ${machine} =~ "Derecho" ]]
		then
		    echo "Creating baselines on ${machine}"
		    export ACCNR=nral0032
		    ./rt.sh -a "${ACCNR}" -c "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    export DISKNM=/glade/derecho/scratch/epicufsrt/ufs-weather-model/RT/
		    cd ${DISKNM}/NEMSfv3gfs/
		    mkdir -p develop-${BL_DATE}
		    cd /glade/derecho/scratch/epicufsrt/FV3_RT
		    ls -l REGRESSION_TEST/.
		    rsync -a --no-t REGRESSION_TEST/ "${DISKNM}/NEMSfv3gfs/develop-${BL_DATE}" || echo "#### Warning! rsync $(pwd)/REGRESSION_TEST/ incomplete."
		    cd ${WORKSPACE}/tests
		    ./rt.sh -a "${ACCNR}" "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		    cd logs/
		    cp "RegressionTests_${machine_id}.log" "$(dirname "${WORKSPACE}")" #/glade/derecho/scratch/epicufsrt/jenkins/workspace
		    git remote -v
		    git fetch --no-recurse-submodules origin
		    git reset FETCH_HEAD --hard
		    cd .. && cd .. && cd ..
		    pwd
		    cp "$(dirname "${WORKSPACE}")/RegressionTests_${machine_id}.log" "${WORKSPACE}/tests/logs/"
		    cd ${WORKSPACE}/tests/
		else
		    echo "Creating baselines on ${machine}"
		    local workflow="-r"
		    ./rt.sh -a "${ACCNR}" -c "${workflow}" "${opt}" "${suite}" | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
		    status=${PIPESTATUS[0]}
		fi

	cd ${WORKSPACE}

	echo "Testing concluded for ${machine}. status=${status}"
	return ${status}
}

create_baseline "${machine}"
