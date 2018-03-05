#!/bin/bash
set -eux

hostname

die() { echo "$@" >&2; exit 1; }
usage() {
  echo
  echo "Usage: $0 -c <model> | -f | -s | -l <file> | -m | -r | -e | -h"
  echo
  echo "  -c  create new baseline results for <model>"
  echo "  -f  run full suite of regression tests"
  echo "  -s  run standard suite of regression tests"
  echo "  -l  runs test specified in <file>"
  echo "  -m  compare against new baseline results"
  echo "  -r  use Rocoto workflow manager"
  echo "  -e  use ecFlow workflow manager"
  echo "  -h  display this help"
  echo
  exit 1
}

[[ $# -eq 0 ]] && usage

rt_trap() {
  [[ ${ROCOTO:-false} == true ]] && rocoto_kill
  [[ ${ECFLOW:-false} == true ]] && { ecflow_kill; ecflow_stop; }
  cleanup
}

cleanup() {
  rm -rf ${LOCKDIR}
  trap 0
  exit
}

trap '{ echo "rt.sh interrupted"; rt_trap ; }' INT
trap '{ echo "rt.sh quit"; rt_trap ; }' QUIT
trap '{ echo "rt.sh terminated"; rt_trap ; }' TERM
trap '{ echo "rt.sh error on line $LINENO"; rt_trap ; }' ERR
trap '{ echo "rt.sh finished"; cleanup ; }' EXIT

# PATHRT - Path to regression tests directory
readonly PATHRT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P )"
cd ${PATHRT}

# PATHTR - Path to nmmb trunk directory
readonly PATHTR=$( cd ${PATHRT}/.. && pwd )

# make sure only one instance of rt.sh is running
readonly LOCKDIR="${PATHRT}"/lock
if mkdir "${LOCKDIR}" ; then
  echo $(hostname) $$ > "${LOCKDIR}/PID"
else
  echo "Only one instance of rt.sh can be running at a time"
  exit 1
fi

source detect_machine.sh
source rt_utils.sh

if [[ $MACHINE_ID = wcoss ]]; then

  source $PATHTR/NEMS/src/conf/module-setup.sh.inc

  set +u
  source /usrx/local/ecflow/setup.sh
  ECFLOW_START=/usrx/local/ecflow/bin/ecflow_start.sh
  set -u
  ROCOTORUN="/u/Christopher.W.Harrop/rocoto/bin/rocotorun"
  ROCOTOSTAT="/u/Christopher.W.Harrop/rocoto/bin/rocotostat"
  DISKNM=/nems/noscrub/emc.nemspara/RT
  MDISK=/global/noscrub
  QUEUE=debug
  ACCNR=GFS-T2O
  STMP=/ptmpp$pex
  PTMP=/ptmpp$pex
  SCHEDULER=lsf
  MPIEXEC=mpirun.lsf
  MPIEXECOPTS=""
# cp fv3_conf/fv3_bsub.IN_wcoss fv3_conf/fv3_bsub.IN

elif [[ $MACHINE_ID = wcoss_cray ]]; then

  source $PATHTR/NEMS/src/conf/module-setup.sh.inc
  module load xt-lsfhpc

  export PATH=/gpfs/hps/nco/ops/ecf/ecfdir/ecflow.v4.1.0.intel/bin:$PATH
  export PYTHONPATH=/gpfs/hps/nco/ops/ecf/ecfdir/ecflow.v4.1.0.intel/lib/python2.6/site-packages
  ECFLOW_START=/gpfs/hps/nco/ops/ecf/ecfdir/ecflow.v4.1.0.intel/bin/ecflow_start.sh
  DISKNM=/gpfs/hps3/emc/nems/noscrub/emc.nemspara/RT
  MDISK=/gpfs/hps3/emc/global/noscrub
  QUEUE=debug
  ACCNR=dev
  if [[ -d /gpfs/hps3/ptmp ]] ; then
      STMP=/gpfs/hps3/stmp
      PTMP=/gpfs/hps3/stmp
  else
      STMP=/gpfs/hps3/stmp
      PTMP=/gpfs/hps3/ptmp
  fi
  SCHEDULER=lsf
  MPIEXEC=aprun
  MPIEXECOPTS="\"-j 1 -n @[TASKS] -N @[TPN] -d 1\""
  cp fv3_conf/fv3_bsub.IN_wcoss_cray fv3_conf/fv3_bsub.IN

elif [[ $MACHINE_ID = theia ]]; then

  source $PATHTR/NEMS/src/conf/module-setup.sh.inc

  module load rocoto
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  export PATH=/scratch4/NCEPDEV/meso/save/Dusan.Jovic/ecflow/bin:$PATH
  export PYTHONPATH=/scratch4/NCEPDEV/meso/save/Dusan.Jovic/ecflow/lib/python2.6/site-packages
  ECFLOW_START=/scratch4/NCEPDEV/meso/save/Dusan.Jovic/ecflow/bin/ecflow_start.sh
  QUEUE=debug
  dprefix=/scratch4/NCEPDEV
  DISKNM=$dprefix/nems/noscrub/emc.nemspara/RT
  STMP=$dprefix/stmp4
  PTMP=$dprefix/stmp3
  SCHEDULER=pbs
  MPIEXEC=mpirun
  MPIEXECOPTS=
  cp fv3_conf/fv3_qsub.IN_theia fv3_conf/fv3_qsub.IN

else
  die "Unknown machine ID, please edit detect_machine.sh file"
fi

mkdir -p ${STMP}/${USER}
mkdir -p ${PTMP}/${USER}

 NEW_BASELINE=${STMP}/${USER}/FV3_RT/REGRESSION_TEST
#NEW_BASELINE=/gpfs/hps3/emc/global/noscrub/Shrinivas.Moorthi/REGRESSION_TEST

RUNDIR_ROOT=${PTMP}/${USER}/FV3_RT/rt_$$
mkdir -p ${RUNDIR_ROOT}

CREATE_BASELINE=false
ROCOTO=false
ECFLOW=false

TESTS_FILE='rt.conf'
SET_ID='standard'
while getopts ":cfsl:mreh" opt; do
  case $opt in
    c)
      CREATE_BASELINE=true
      SET_ID=' '
      ;;
    f)
      SET_ID=' '
      ;;
    s)
      SET_ID='standard'
      ;;
    l)
      TESTS_FILE=$OPTARG
      SET_ID=' '
      ;;
    m)
      # redefine RTPWD to point to newly created baseline outputs
      RTPWD=${NEW_BASELINE}
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

RTPWD=${RTPWD:-$DISKNM/NEMSfv3gfs/trunk-20180226}

shift $((OPTIND-1))
[[ $# -gt 1 ]] && usage

if [[ $CREATE_BASELINE == true ]]; then
  #
  # prepare new regression test directory
  #
  rm -rf "${NEW_BASELINE}"
  mkdir -p "${NEW_BASELINE}"
  echo "copy REGRESSION_TEST_baselines"
  cp -rf "$RTPWD"/* "$NEW_BASELINE"/.

fi

COMPILE_LOG=${PATHRT}/Compile_$MACHINE_ID.log
REGRESSIONTEST_LOG=${PATHRT}/RegressionTests_$MACHINE_ID.log

date > ${REGRESSIONTEST_LOG}
echo "Start Regression test" >> ${REGRESSIONTEST_LOG}
echo                         >> ${REGRESSIONTEST_LOG}

source default_vars.sh

TEST_NR=0
COMPILE_NR=0
rm -f fail_test

LOG_DIR=${PATHRT}/log_$MACHINE_ID
rm -rf ${LOG_DIR}
mkdir ${LOG_DIR}

rm -f ../fv3.exe

if [[ $ROCOTO == true ]]; then

  ROCOTO_XML=${PATHRT}/rocoto_workflow.xml
  ROCOTO_DB=${PATHRT}/rocoto_workflow.db

  rm -f $ROCOTO_XML $ROCOTO_DB

  if [[ $MACHINE_ID = wcoss ]]; then
    QUEUE=dev
    COMPILE_QUEUE=dev
    ROCOTO_SCHEDULER=lsf
  elif [[ $MACHINE_ID = theia ]]; then
    QUEUE=batch
    COMPILE_QUEUE=service
    ROCOTO_SCHEDULER=moabtorque
  else
    die "Rocoto is not supported on this machine $MACHINE_ID"
  fi

  cat << EOF > $ROCOTO_XML
<?xml version="1.0"?>
<!DOCTYPE workflow
[
  <!ENTITY PATHRT       "${PATHRT}">
  <!ENTITY LOG          "${LOG_DIR}">
  <!ENTITY PATHTR       "${PATHTR}">
  <!ENTITY RTPWD        "${RTPWD}">
  <!ENTITY RUNDIR_ROOT  "${RUNDIR_ROOT}">
  <!ENTITY NEW_BASELINE "${NEW_BASELINE}">
]>
<workflow realtime="F" scheduler="${ROCOTO_SCHEDULER}" taskthrottle="20">
  <cycledef>197001010000 197001010000 01:00:00</cycledef>
  <log>&LOG;/workflow.log</log>
EOF

fi

if [[ $ECFLOW == true ]]; then

  ECFLOW_RUN=${PATHRT}/ecflow_run
  rm -rf ${ECFLOW_RUN}
  mkdir -p ${ECFLOW_RUN}/regtest
  cp head.h tail.h ${ECFLOW_RUN}
  > ${ECFLOW_RUN}/regtest.def
  cat << EOF >> ${ECFLOW_RUN}/regtest.def
suite regtest
    edit ECF_HOME '${ECFLOW_RUN}'
    edit ECF_INCLUDE '${ECFLOW_RUN}'
    edit ECF_KILL_CMD kill -15 %ECF_RID% > %ECF_JOB%.kill 2>&1
    edit ECF_TRIES 1
    limit max_jobs 10
EOF

  if [[ $MACHINE_ID = wcoss ]]; then
    QUEUE=dev
  elif [[ $MACHINE_ID = wcoss_cray ]]; then
    QUEUE=dev
  elif [[ $MACHINE_ID = theia ]]; then
    QUEUE=batch
  else
    die "ecFlow is not supported on this machine $MACHINE_ID"
  fi

fi

##
## read rt.conf and then either execute the test script directly or create
## workflow description file
##

new_compile=false
in_metatask=false

while read -r line; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ $line == \#* ]] && continue

  if [[ $line == COMPILE* ]] ; then

      unset APP
      NEMS_VER=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
      SET=$(     echo $line | cut -d'|' -f3)
      MACHINES=$(echo $line | cut -d'|' -f4 | sed -e 's/^ *//' -e 's/ *$//')

      [[ $SET_ID != ' ' && $SET != *${SET_ID}* ]] && continue
      [[ $MACHINES != ' ' && $MACHINES != "${MACHINE_ID}" ]] && continue

      COMPILE_NR_DEP=${COMPILE_NR}
      (( COMPILE_NR += 1 ))

      if [[ $ROCOTO == true ]]; then
        rocoto_create_compile_task
      elif [[ $ECFLOW == true ]]; then
        ecflow_create_compile_task
      else
        ./compile.sh $PATHTR/FV3 $MACHINE_ID "${NEMS_VER}" $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1
#       ./compile.sh $PATHTR/FV3 $MACHINE_ID DEBUG=Y  $COMPILE_NR > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1
        echo " bash Compile is done"
      fi

    continue

  elif [[ $line == APPBUILD* ]] ; then

      unset NEMS_VER
      APP=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
      SET=$(     echo $line | cut -d'|' -f3)
      MACHINES=$(echo $line | cut -d'|' -f4 | sed -e 's/^ *//' -e 's/ *$//')

      [[ $SET_ID != ' ' && $SET != *${SET_ID}* ]] && continue
      [[ $MACHINES != ' ' && $MACHINES != "${MACHINE_ID}" ]] && continue

      COMPILE_NR_DEP=${COMPILE_NR}
      (( COMPILE_NR += 1 ))

      if [[ $ROCOTO == true ]]; then
        rocoto_create_compile_task
      elif [[ $ECFLOW == true ]]; then
        ecflow_create_compile_task
      else
          echo test  > "${LOG_DIR}/compile_${COMPILE_NR}.log" 2>&1
          test -s ./appbuild.sh
          test -x ./appbuild.sh
        ./appbuild.sh "$PATHTR/FV3" "$APP" "$COMPILE_NR" 2>&1 | tee "${LOG_DIR}/compile_${COMPILE_NR}.log"
        echo " bash NEMSAppBuilder is done"
      fi

      unset APP

  elif [[ $line == RUN* ]] ; then

    TEST_NAME=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
    SET=$(      echo $line | cut -d'|' -f3)
    MACHINES=$( echo $line | cut -d'|' -f4)
    CB=$(       echo $line | cut -d'|' -f5)
    DEP_RUN=$(  echo $line | cut -d'|' -f6 | sed -e 's/^ *//' -e 's/ *$//')
    [[ -e "tests/$TEST_NAME" ]] || die "run test file tests/$TEST_NAME does not exist"
    [[ $SET_ID != ' ' && $SET != *${SET_ID}* ]] && continue
    [[ $MACHINES != ' ' && $MACHINES != *${MACHINE_ID}* ]] && continue
    [[ $CREATE_BASELINE == true && $CB != *fv3* ]] && continue

    if [[ $ROCOTO == true && $new_compile == true ]]; then
      new_compile=false
      in_metatask=true
      cat << EOF >> $ROCOTO_XML
  <metatask name="${NEMS_VER}"><var name="zero">0</var>
EOF
    fi

    TEST_NR=$( printf '%02d' $(( 10#$TEST_NR + 1 )) )

    (
      source ${PATHRT}/tests/$TEST_NAME

      cat << EOF > run_test.env
      export MACHINE_ID=${MACHINE_ID}
      export RTPWD=${RTPWD}
      export PATHRT=${PATHRT}
      export PATHTR=${PATHTR}
      export NEW_BASELINE=${NEW_BASELINE}
      export MPIEXEC=${MPIEXEC}
      export MPIEXECOPTS=${MPIEXECOPTS}
      export CREATE_BASELINE=${CREATE_BASELINE}
      export SCHEDULER=${SCHEDULER}
      export ACCNR=${ACCNR}
      export QUEUE=${QUEUE}
      export ROCOTO=${ROCOTO}
      export LOG_DIR=${LOG_DIR}
EOF

      if [[ $ROCOTO == true ]]; then
        rocoto_create_run_task
      elif [[ $ECFLOW == true ]]; then
        ecflow_create_run_task
      else
        ./run_test.sh ${PATHRT} ${RUNDIR_ROOT} ${TEST_NAME} ${TEST_NR} ${COMPILE_NR} > ${LOG_DIR}/run_${TEST_NAME}.log 2>&1
      fi
    )

    continue
  else
    die "Unknown command $line"
  fi
done < $TESTS_FILE

##
## run regression test workflow (currently Rocoto or ecFlow are supported)
##

if [[ $ROCOTO == true ]]; then
  if [[ $in_metatask == true ]]; then
    echo "  </metatask>" >> $ROCOTO_XML
  fi
  echo "</workflow>" >> $ROCOTO_XML
  # run rocoto workflow until done
  rocoto_run
fi

if [[ $ECFLOW == true ]]; then
  echo "endsuite" >> ${ECFLOW_RUN}/regtest.def
  # run ecflow workflow until done
  ecflow_run
  ecflow_stop
fi

##
## regression test is either failed or successful
##
set +e
cat ${LOG_DIR}/compile_*.log                   >  ${COMPILE_LOG}
cat ${LOG_DIR}/rt_*.log                        >> ${REGRESSIONTEST_LOG}
if [[ -e fail_test ]]; then
  echo "FAILED TESTS: "
  echo "FAILED TESTS: "                        >> ${REGRESSIONTEST_LOG}
  while read -r failed_test_name
  do
    echo "Test ${failed_test_name} failed "
    echo "Test ${failed_test_name} failed "    >> ${REGRESSIONTEST_LOG}
  done < fail_test
   echo ; echo REGRESSION TEST FAILED
  (echo ; echo REGRESSION TEST FAILED)         >> ${REGRESSIONTEST_LOG}
else
   echo ; echo REGRESSION TEST WAS SUCCESSFUL
  (echo ; echo REGRESSION TEST WAS SUCCESSFUL) >> ${REGRESSIONTEST_LOG}

  rm -f fv3_*.x fv3_*.exe modules.fv3_* run_test.env
  [[ ${ROCOTO:-false} == true ]] && rm -f ${ROCOTO_XML} ${ROCOTO_DB}
  [[ ${ECFLOW:-false} == true ]] && rm -rf ${ECFLOW_RUN}
fi

date >> ${REGRESSIONTEST_LOG}

exit
