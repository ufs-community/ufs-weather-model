#!/bin/bash -x

set_up_env() {

   export PATH=${PATH}:~/bin
   export ACCNR=epic
   export machine=${NODE_NAME}
   export opt="-l"
   export suite="rt.conf"
      [[ -n ${WM_OPERATIONAL_TESTS}                 ]] && opt="-n" && suite="${WM_OPERATIONAL_TESTS} ${UFS_COMPILER}" || return 0
      [[    ${WM_OPERATIONAL_TESTS} = default       ]] && opt="-n" && suite="control_p8 ${UFS_COMPILER}"
      [[    ${WM_OPERATIONAL_TESTS} = comprehensive ]] && opt="-l" && suite="rt.conf"
      [[    ${WM_OPERATIONAL_TESTS} = rt.conf       ]] && opt="-l" && suite="rt.conf"
      [[   "${suite}"               = rt.conf       ]] && opt="-l"

   set -eu

   SCRIPT_REALPATH=$(realpath "${BASH_SOURCE[0]}")
   SCRIPTS_DIR=$(dirname "${SCRIPT_REALPATH}")
   UFS_MODEL_DIR=$(realpath "${SCRIPTS_DIR}/../..")
   readonly UFS_MODEL_DIR
   echo "UFS MODEL DIR: ${UFS_MODEL_DIR}"

   BUILD_DIR=${BUILD_DIR:-${UFS_MODEL_DIR}/build}
   TESTS_DIR=${TESTS_DIR:-${UFS_MODEL_DIR}/tests}

   pwd
   echo "USER=${USER}"
   echo "ACCOUNT=${ACCNR}"
   echo "MACHINE=${machine}"
   echo "WORKSPACE=${WORKSPACE}"

}

function update_baselines() {
   echo "Creating baselines on ${machine}"
   # shellcheck disable=SC2155
   export PREVIOUS_BL_DATE=$(ls -d develop-* | tail -n 1 | cut -d '-' -f2)
   ./rt.sh -a "${ACCNR}" -ecs test_changes.list | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
   status=${PIPESTATUS[0]}
   if [[ "${machine}" = "GaeaC6" ]]; then
      unset LD_LIBRARY_PATH
   fi
   cd ${DISKNM}/NEMSfv3gfs/
   if [[ ! -d "develop-${BL_DATE}" ]]; then
      mkdir -p "develop-${BL_DATE}"
   fi
   cd "develop-${BL_DATE}"
   rsync -rvh "../develop-${PREVIOUS_BL_DATE}/" .
   cd "${STMP}"
   rsync -a --no-t REGRESSION_TEST/ "${DISKNM}/NEMSfv3gfs/develop-${BL_DATE}" || echo "#### Warning! rsync $(pwd)/REGRESSION_TEST/ incomplete."
   cd ${DISKNM}/NEMSfv3gfs/"develop-${BL_DATE}"
   chmod -R 755 .
   cd "${WORKSPACE}/tests"
   ./rt.sh -a "${ACCNR}" -e | tee "${WORKSPACE}/tests/logs/RT-run-${machine}.log"
   status=${PIPESTATUS[0]}
   cd logs/
   cp "RegressionTests_${machine_id}.log" "$(dirname "${WORKSPACE}")" #/work2/noaa/epic/role-epic/jenkins/workspace
   git fetch --no-recurse-submodules origin
   git reset FETCH_HEAD --hard
   cd ../../..
   pwd
   cp "$(dirname "${WORKSPACE}")/RegressionTests_${machine_id}.log" "${WORKSPACE}/tests/logs/"
   cd "${WORKSPACE}/tests/"

}

function create_baseline() {
   local machine=${1:-${NODE_NAME}}
   local machine_id=${machine,,} # tolower
   local WORKSPACE
   WORKSPACE="$(pwd)"
   local status=0

   git submodule update --init --recursive
   cd tests
   pwd
   # shellcheck disable=SC2015
   [[ ${UFS_PLATFORM} =~ clusternoaa ]] && echo "export BL_DATE=20240426" > bl_date.conf || cat bl_date.conf

   mkdir -p logs/
   BL_DATE=$(cut -d '=' -f2 bl_date.conf)
   export BL_DATE

   if [[ ${machine} =~ "Hercules" ]]
   then
       export dprefix=/work2/noaa/${ACCNR}/${USER}
       sed "s|/noaa/stmp/|/noaa/${ACCNR}/stmp/|g" -i rt.sh
       export DISKNM=/work2/noaa/epic/hercules/UFS-WM_RT
       export STMP="${dprefix}/stmp/role-epic/FV3_RT"
       update_baselines
   elif [[ ${machine} =~ "Orion" ]]
   then
      cd ..
      git submodule update --init --recursive
      cd tests
      export dprefix=/work/noaa/${ACCNR}/${USER}
      sed "s|/noaa/stmp/|/noaa/${ACCNR}/stmp/|g" -i rt.sh
      export DISKNM=/work2/noaa/epic/UFS-WM_RT
      export STMP="${dprefix}/stmp/role-epic/FV3_RT"
      update_baselines
   elif [[ ${machine} =~ "Gaea" ]]
   then
      export DISKNM=/gpfs/f6/bil-fire8/world-shared/role.epic/UFS-WM_RT
      export STMP="/gpfs/f6/epic/proj-shared/role.epic/RT_BASELINE/role.epic/FV3_RT"
      update_baselines
   elif [[ ${machine} =~ "Ursa" ]]
   then
      sed "s|QUEUE=batch|QUEUE=windfall|g" -i rt.sh
      export DISKNM=/scratch4/NAGAPE/epic/role-epic/UFS-WM_RT
      export STMP="/scratch4/NAGAPE/epic/role.epic/stmp/RT_BASELINE/role.epic/FV3_RT"
      update_baselines
   else
      echo "${machine} not recognized."
      exit 1
      status=${PIPESTATUS[0]}
   fi

   cd "${WORKSPACE}"

   echo "Testing concluded for ${machine}. status=${status}"
   return "${status}"
}

set_up_env
create_baseline "${machine}"
