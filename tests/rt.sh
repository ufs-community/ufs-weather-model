#!/bin/bash
set -eux

SECONDS=0

hostname

die() { echo "$@" >&2; exit 1; }
usage() {
  set +x
  echo
  echo "Usage: $0 -c | -e | -h | -k | -w | -d | -l <file> | -m | -n <name> | -r "
  echo
  echo "  -c  create new baseline results"
  echo "  -e  use ecFlow workflow manager"
  echo "  -h  display this help"
  echo "  -k  keep run directory after rt.sh is completed"
  echo "  -l  runs test specified in <file>"
  echo "  -m  compare against new baseline results"
  echo "  -n  run single test <name>"
  echo "  -r  use Rocoto workflow manager"
  echo "  -w  for weekly_test, skip comparing baseline results"
  echo "  -d  delete run direcotries that are not used by other tests"
  echo
  set -x
  exit 1
}

[[ $# -eq 0 ]] && usage

rt_single() {
  rm -f $RT_SINGLE_CONF
  local compile_line=''
  local run_line=''
  while read -r line || [ "$line" ]; do
    line="${line#"${line%%[![:space:]]*}"}"
    [[ ${#line} == 0 ]] && continue
    [[ $line == \#* ]] && continue

    if [[ $line == COMPILE* ]] ; then
      MACHINES=$(echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
      if [[ ${MACHINES} == '' ]]; then
        compile_line=$line
      elif [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || compile_line=$line
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && compile_line=$line
      fi
    fi

    if [[ $line =~ RUN ]]; then
      tmp_test=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
      if [[ $SINGLE_NAME == $tmp_test && $compile_line != '' ]]; then
        echo $compile_line > $RT_SINGLE_CONF
        dep_test=$(echo $line | grep -w $tmp_test | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
        if [[ $dep_test != '' ]]; then
          dep_line=$(cat rt.conf | grep -w "$dep_test" | grep -v "$tmp_test")
          dep_line="${dep_line#"${dep_line%%[![:space:]]*}"}"
          echo $dep_line >> $RT_SINGLE_CONF
        fi
        echo $line >> $RT_SINGLE_CONF
        break
      fi
    fi
  done < $TESTS_FILE

  if [[ ! -f $RT_SINGLE_CONF ]]; then
    echo "$SINGLE_NAME does not exist or cannot be run on $MACHINE_ID"
    exit 1
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
  rm -rf ${LOCKDIR}
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

readonly RT_SINGLE_CONF='rt_single.conf'

# Default compiler "intel"
export RT_COMPILER=${RT_COMPILER:-intel}

source detect_machine.sh # Note: this does not set ACCNR. The "if" block below does.
source rt_utils.sh

source module-setup.sh

if [[ $MACHINE_ID = wcoss2.* ]]; then

  #module use /usrx/local/dev/emc_rocoto/modulefiles
  #module load ruby/2.5.1 rocoto/1.3.0rc2
  #ROCOTORUN=$(which rocotorun)
  #ROCOTOSTAT=$(which rocotostat)
  #ROCOTOCOMPLETE=$(which rocotocomplete)
  #ROCOTO_SCHEDULER=lsf

  module load ecflow/5.6.0.13
  module load gcc/10.3.0 python/3.8.6
  ECFLOW_START=${ECF_ROOT}/scripts/server_check.sh
  export ECF_OUTPUTDIR=${PATHRT}/ecf_outputdir
  export ECF_COMDIR=${PATHRT}/ecf_comdir
  rm -rf ${ECF_OUTPUTDIR} ${ECF_COMDIR}
  mkdir -p ${ECF_OUTPUTDIR}
  mkdir -p ${ECF_COMDIR}
  export colonifnco=":output"  # hack

  DISKNM=/lfs/h2/emc/nems/noscrub/emc.nems/RT
  QUEUE=dev
  COMPILE_QUEUE=dev
  PARTITION=
  ACCNR="${ACCNR:-GFS-DEV}"
  STMP=/lfs/h2/emc/ptmp
  PTMP=/lfs/h2/emc/ptmp
  SCHEDULER=pbs
  cp fv3_conf/fv3_qsub.IN_wcoss2 fv3_conf/fv3_qsub.IN
  cp fv3_conf/compile_qsub.IN_wcoss2 fv3_conf/compile_qsub.IN

elif [[ $MACHINE_ID = acorn.* ]]; then

  module load ecflow/5.6.0.13
  module load gcc/10.3.0 python/3.8.6
  ECFLOW_START=${ECF_ROOT}/scripts/server_check.sh
  export ECF_OUTPUTDIR=${PATHRT}/ecf_outputdir
  export ECF_COMDIR=${PATHRT}/ecf_comdir
  rm -rf ${ECF_OUTPUTDIR} ${ECF_COMDIR}
  mkdir -p ${ECF_OUTPUTDIR}
  mkdir -p ${ECF_COMDIR}
  export colonifnco=":output"  # hack

  DISKNM=/lfs/h1/emc/nems/noscrub/emc.nems/RT
  QUEUE=dev
  COMPILE_QUEUE=dev
  PARTITION=
  ACCNR="${ACCNR:-GFS-DEV}"
  STMP=/lfs/h2/emc/ptmp
  PTMP=/lfs/h2/emc/ptmp
  SCHEDULER=pbs
  cp fv3_conf/fv3_qsub.IN_wcoss2 fv3_conf/fv3_qsub.IN
  cp fv3_conf/compile_qsub.IN_wcoss2 fv3_conf/compile_qsub.IN

elif [[ $MACHINE_ID = gaea.* ]]; then

  export PATH=/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/envs/ufs-weather-model/bin:/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/bin:$PATH
  export PYTHONPATH=/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/lib/python3.8/site-packages
  ECFLOW_START=/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/envs/ufs-weather-model/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  DISKNM=/lustre/f2/pdata/ncep_shared/emc.nemspara/RT
  QUEUE=normal
  COMPILE_QUEUE=normal
#  ACCNR="${ACCNR:-cmp}"
  PARTITION=c4
  STMP=/lustre/f2/scratch
  PTMP=/lustre/f2/scratch

  SCHEDULER=slurm
  cp fv3_conf/fv3_slurm.IN_gaea fv3_conf/fv3_slurm.IN
  cp fv3_conf/compile_slurm.IN_gaea fv3_conf/compile_slurm.IN

elif [[ $MACHINE_ID = hera.* ]]; then

  module load rocoto
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  PYTHONHOME=/scratch1/NCEPDEV/nems/emc.nemspara/soft/miniconda3_new_20210629
  export PATH=$PYTHONHOME/bin:$PATH
  export PYTHONPATH=$PYTHONHOME/lib/python3.7/site-packages

  module load ecflow
  ECFLOW_START=ecflow_start.sh

  QUEUE=batch
  COMPILE_QUEUE=batch

  #ACCNR="${ACCNR:-fv3-cpu}
  PARTITION=
  dprefix=/scratch1/NCEPDEV
  DISKNM=$dprefix/nems/emc.nemspara/RT
  STMP=$dprefix/stmp4
  PTMP=$dprefix/stmp2

  SCHEDULER=slurm
  cp fv3_conf/fv3_slurm.IN_hera fv3_conf/fv3_slurm.IN
  cp fv3_conf/compile_slurm.IN_hera fv3_conf/compile_slurm.IN

elif [[ $MACHINE_ID = orion.* ]]; then

  module load gcc/8.3.0

  module load contrib rocoto/1.3.1
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  export PATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/bin:$PATH
  export PYTHONPATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/lib/python3.8/site-packages
  ECFLOW_START=/work/noaa/nems/emc.nemspara/soft/miniconda3/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  QUEUE=batch
  COMPILE_QUEUE=batch
  PARTITION=orion
  dprefix=/work/noaa/stmp/${USER}
  DISKNM=/work/noaa/nems/emc.nemspara/RT
  STMP=$dprefix/stmp
  PTMP=$dprefix/stmp

  SCHEDULER=slurm
  cp fv3_conf/fv3_slurm.IN_orion fv3_conf/fv3_slurm.IN
  cp fv3_conf/compile_slurm.IN_orion fv3_conf/compile_slurm.IN

elif [[ $MACHINE_ID = jet.* ]]; then

  module load rocoto/1.3.2
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  export PATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/bin:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/bin:$PATH
  export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/lib/python3.8/site-packages
  module load ecflow
  ECFLOW_START=/apps/ecflow/5.5.3/bin/ecflow_start.sh

  QUEUE=batch
  COMPILE_QUEUE=batch
  ACCNR="${ACCNR:-h-nems}"
  PARTITION=xjet
  DISKNM=/mnt/lfs4/HFIP/hfv3gfs/role.epic/RT
  dprefix=${dprefix:-/lfs4/HFIP/$ACCNR/$USER}
  STMP=${STMP:-$dprefix/RT_BASELINE}
  PTMP=${PTMP:-$dprefix/RT_RUNDIRS}

  SCHEDULER=slurm
  cp fv3_conf/fv3_slurm.IN_jet fv3_conf/fv3_slurm.IN
  cp fv3_conf/compile_slurm.IN_jet fv3_conf/compile_slurm.IN

elif [[ $MACHINE_ID = s4.* ]]; then

  module load rocoto/1.3.2
  module load ecflow/5.6.0
  module load miniconda/3.8-s4
  ROCOTORUN=$(which rocotorun)
  ROCOTOSTAT=$(which rocotostat)
  ROCOTOCOMPLETE=$(which rocotocomplete)
  ROCOTO_SCHEDULER=slurm

  ECFLOW_START=/opt/ecflow/5.6.0/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  QUEUE=s4
  COMPILE_QUEUE=s4

  ACCNR="${ACCNR:-star}"
  PARTITION=s4
  dprefix=/data/prod
  DISKNM=$dprefix/emc.nemspara/RT
  STMP=/scratch/short/users
  PTMP=/scratch/users

  SCHEDULER=slurm
  cp fv3_conf/fv3_slurm.IN_s4 fv3_conf/fv3_slurm.IN
  cp fv3_conf/compile_slurm.IN_s4 fv3_conf/compile_slurm.IN

elif [[ $MACHINE_ID = cheyenne.* ]]; then

  export PATH=/glade/p/ral/jntp/tools/miniconda3/4.8.3/envs/ufs-weather-model/bin:/glade/p/ral/jntp/tools/miniconda3/4.8.3/bin:$PATH
  export PYTHONPATH=/glade/p/ral/jntp/tools/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/glade/p/ral/jntp/tools/miniconda3/4.8.3/lib/python3.8/site-packages
  ECFLOW_START=/glade/p/ral/jntp/tools/miniconda3/4.8.3/envs/ufs-weather-model/bin/ecflow_start.sh
  ECF_PORT=$(( $(id -u) + 1500 ))

  QUEUE=regular
  COMPILE_QUEUE=regular
  PARTITION=
  dprefix=/glade/scratch
  DISKNM=/glade/scratch/epicufsrt/GMTB/ufs-weather-model/RT
  STMP=$dprefix
  PTMP=$dprefix
  SCHEDULER=pbs
  cp fv3_conf/fv3_qsub.IN_cheyenne fv3_conf/fv3_qsub.IN
  cp fv3_conf/compile_qsub.IN_cheyenne fv3_conf/compile_qsub.IN

elif [[ $MACHINE_ID = stampede.* ]]; then

  export PYTHONPATH=
  ECFLOW_START=
  QUEUE=skx-normal
  COMPILE_QUEUE=skx-dev
  PARTITION=
  ACCNR="${ACCNR:-TG-EES200015}"
  dprefix=$SCRATCH/ufs-weather-model/run
  DISKNM=/work2/07736/minsukji/stampede2/ufs-weather-model/RT
  STMP=$dprefix
  PTMP=$dprefix
  SCHEDULER=slurm
  MPIEXEC=ibrun
  MPIEXECOPTS=
  cp fv3_conf/fv3_slurm.IN_stampede fv3_conf/fv3_slurm.IN

elif [[ $MACHINE_ID = expanse.* ]]; then

  export PYTHONPATH=
  ECFLOW_START=
  QUEUE=compute
  COMPILE_QUEUE=shared
  PARTITION=
  ACCNR="${ACCNR:-TG-EES200015}"
  dprefix=/expanse/lustre/scratch/$USER/temp_project/run
  DISKNM=/expanse/lustre/scratch/domh/temp_project/RT
  STMP=$dprefix
  PTMP=$dprefix
  SCHEDULER=slurm
  cp fv3_conf/fv3_slurm.IN_expanse fv3_conf/fv3_slurm.IN

else
  die "Unknown machine ID, please edit detect_machine.sh file"
fi

# If account is unspecified, assume the machine has a "nems"
# accounting code.
export ACCNR="${ACCNR:-nems}"

# Display the machine and account using the format detect_machine.sh used:
echo "Machine: " $MACHINE_ID "    Account: " $ACCNR

mkdir -p ${STMP}/${USER}

NEW_BASELINE=${STMP}/${USER}/FV3_RT/REGRESSION_TEST_${RT_COMPILER^^}

# Overwrite default RUNDIR_ROOT if environment variable RUNDIR_ROOT is set
RUNDIR_ROOT=${RUNDIR_ROOT:-${PTMP}/${USER}/FV3_RT}/rt_$$
mkdir -p ${RUNDIR_ROOT}

CREATE_BASELINE=false
ROCOTO=false
ECFLOW=false
KEEP_RUNDIR=false
SINGLE_NAME=''
TEST_35D=false
export skip_check_results=false
export delete_rundir=false

TESTS_FILE='rt.conf'

while getopts ":cl:mn:dwkreh" opt; do
  case $opt in
    c)
      CREATE_BASELINE=true
      ;;
    l)
      TESTS_FILE=$OPTARG
      ;;
    m)
      # redefine RTPWD to point to newly created baseline outputs
      RTPWD=${NEW_BASELINE}
      ;;
    n)
      SINGLE_NAME=$OPTARG
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

if [[ $SINGLE_NAME != '' ]]; then
  rt_single
  TESTS_FILE=$RT_SINGLE_CONF
fi

if [[ $TESTS_FILE =~ '35d' ]] || [[ $TESTS_FILE =~ 'weekly' ]]; then
  TEST_35D=true
fi


BL_DATE=20230324

RTPWD=${RTPWD:-$DISKNM/NEMSfv3gfs/develop-${BL_DATE}/${RT_COMPILER^^}}

INPUTDATA_ROOT=${INPUTDATA_ROOT:-$DISKNM/NEMSfv3gfs/input-data-20221101}
INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT}/WW3_input_data_20220624
INPUTDATA_ROOT_BMIC=${INPUTDATA_ROOT_BMIC:-$DISKNM/NEMSfv3gfs/BM_IC-20220207}

shift $((OPTIND-1))
[[ $# -gt 1 ]] && usage

if [[ $CREATE_BASELINE == true ]]; then
  #
  # prepare new regression test directory
  #
  rm -rf "${NEW_BASELINE}"
  mkdir -p "${NEW_BASELINE}"
fi

if [[ $skip_check_results == true ]]; then
  REGRESSIONTEST_LOG=${PATHRT}/RegressionTests_weekly_$MACHINE_ID.log
else
  REGRESSIONTEST_LOG=${PATHRT}/RegressionTests_$MACHINE_ID.log
fi

date > ${REGRESSIONTEST_LOG}
echo "Start Regression test" >> ${REGRESSIONTEST_LOG}
echo                         >> ${REGRESSIONTEST_LOG}

source default_vars.sh

JOB_NR=0
TEST_NR=0
COMPILE_NR=0
rm -f fail_test* fail_compile*

export LOG_DIR=${PATHRT}/log_$MACHINE_ID
rm -rf ${LOG_DIR}
mkdir ${LOG_DIR}

if [[ $ROCOTO == true ]]; then

  ROCOTO_XML=${PATHRT}/rocoto_workflow.xml
  ROCOTO_STATE=${PATHRT}/rocoto_workflow.state
  ROCOTO_DB=${PATHRT}/rocoto_workflow.db

  rm -f $ROCOTO_XML $ROCOTO_DB $ROCOTO_STATE *_lock.db

  if [[ $MACHINE_ID = wcoss2.* || $MACHINE_ID = acorn.* ]]; then
    QUEUE=dev
    COMPILE_QUEUE=dev
    ROCOTO_SCHEDULER=pbs
  elif [[ $MACHINE_ID = hera.* ]]; then
    QUEUE=batch
    COMPILE_QUEUE=batch
    ROCOTO_SCHEDULER=slurm
  elif [[ $MACHINE_ID = orion.* ]]; then
    QUEUE=batch
    COMPILE_QUEUE=batch
    ROCOTO_SCHEDULER=slurm
  elif [[ $MACHINE_ID = s4.* ]]; then
    QUEUE=s4
    COMPILE_QUEUE=s4
    ROCOTO_SCHEDULER=slurm
  elif [[ $MACHINE_ID = jet.* ]]; then
    QUEUE=batch
    COMPILE_QUEUE=batch
    ROCOTO_SCHEDULER=slurm
  else
    die "Rocoto is not supported on this machine $MACHINE_ID"
  fi

  cat << EOF > $ROCOTO_XML
<?xml version="1.0"?>
<!DOCTYPE workflow
[
  <!ENTITY PATHRT         "${PATHRT}">
  <!ENTITY LOG            "${LOG_DIR}">
  <!ENTITY PATHTR         "${PATHTR}">
  <!ENTITY RTPWD          "${RTPWD}">
  <!ENTITY INPUTDATA_ROOT "${INPUTDATA_ROOT}">
  <!ENTITY INPUTDATA_ROOT_WW3 "${INPUTDATA_ROOT_WW3}">
  <!ENTITY INPUTDATA_ROOT_BMIC "${INPUTDATA_ROOT_BMIC}">
  <!ENTITY RUNDIR_ROOT    "${RUNDIR_ROOT}">
  <!ENTITY NEW_BASELINE   "${NEW_BASELINE}">
]>
<workflow realtime="F" scheduler="${ROCOTO_SCHEDULER}" taskthrottle="20">
  <cycledef>197001010000 197001010000 01:00:00</cycledef>
  <log>&LOG;/workflow.log</log>
EOF

fi

if [[ $ECFLOW == true ]]; then

  # Default maximum number of compile and run jobs
  MAX_BUILDS=10
  MAX_JOBS=30

  # Default number of tries to run jobs - on wcoss, no error tolerance
  ECF_TRIES=2
  if [[ $MACHINE_ID = wcoss* ]]; then
    ECF_TRIES=1
  fi

  # Reduce maximum number of compile jobs on jet.intel and s4.intel because of licensing issues
  if [[ $MACHINE_ID = jet.intel ]]; then
    MAX_BUILDS=5
  elif [[ $MACHINE_ID = s4.intel ]]; then
    MAX_BUILDS=1
  fi

  if [[ $MACHINE_ID = hera.* ]] && [[ ! $HOSTNAME = hecflow* ]]; then
    echo "ERROR: To use ECFlow on Hera we must be logged into 'hecflow01' login node."
    exit 1
  fi

  ECFLOW_RUN=${PATHRT}/ecflow_run
  ECFLOW_SUITE=regtest_$$
  rm -rf ${ECFLOW_RUN}
  mkdir -p ${ECFLOW_RUN}/${ECFLOW_SUITE}
  cp head.h tail.h ${ECFLOW_RUN}
  cat << EOF > ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
suite ${ECFLOW_SUITE}
    edit ECF_HOME '${ECFLOW_RUN}'
    edit ECF_INCLUDE '${ECFLOW_RUN}'
    edit ECF_KILL_CMD kill -15 %ECF_RID% > %ECF_JOB%.kill 2>&1
    edit ECF_TRIES ${ECF_TRIES}
    label src_dir '${PATHTR}'
    label run_dir '${RUNDIR_ROOT}'
    limit max_builds ${MAX_BUILDS}
    limit max_jobs ${MAX_JOBS}
EOF

  if [[ $MACHINE_ID = wcoss2.* || $MACHINE_ID = acorn.* ]]; then
    QUEUE=dev
  elif [[ $MACHINE_ID = hera.* ]]; then
    QUEUE=batch
  elif [[ $MACHINE_ID = orion.* ]]; then
    QUEUE=batch
  elif [[ $MACHINE_ID = jet.* ]]; then
    QUEUE=batch
  elif [[ $MACHINE_ID = s4.* ]]; then
    QUEUE=s4
  elif [[ $MACHINE_ID = gaea.* ]]; then
    QUEUE=normal
  elif [[ $MACHINE_ID = cheyenne.* ]]; then
    QUEUE=regular
  else
    die "ecFlow is not supported on this machine $MACHINE_ID"
  fi

else

  if [[ $MACHINE_ID = hera.* ]] && [[ $HOSTNAME = hecflow* ]]; then
    echo "ERROR: To run without using ECFlow on Hera, please do not use ecflow node."
    exit 1
  fi

fi

##
## read rt.conf and then either execute the test script directly or create
## workflow description file
##

new_compile=false
in_metatask=false

[[ -f $TESTS_FILE ]] || die "$TESTS_FILE does not exist"

while read -r line || [ "$line" ]; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ $line == \#* ]] && continue

  JOB_NR=$( printf '%03d' $(( 10#$JOB_NR + 1 )) )

  if [[ $line == COMPILE* ]] ; then

    MAKE_OPT=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
    MACHINES=$(echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
    CB=$(      echo $line | cut -d'|' -f4)

    [[ $CREATE_BASELINE == true && $CB != *fv3* ]] && continue

    if [[ ${MACHINES} != '' ]]; then
      if [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && continue
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || continue
      else
        echo "MACHINES=|${MACHINES}|"
        die "MACHINES spec must be either an empty string or start with either '+' or '-'"
      fi
    fi

    export COMPILE_NR=$( printf '%03d' $(( 10#$COMPILE_NR + 1 )) )

    cat << EOF > ${RUNDIR_ROOT}/compile_${COMPILE_NR}.env
    export JOB_NR=${JOB_NR}
    export COMPILE_NR=${COMPILE_NR}
    export MACHINE_ID=${MACHINE_ID}
    export RT_COMPILER=${RT_COMPILER}
    export PATHRT=${PATHRT}
    export PATHTR=${PATHTR}
    export SCHEDULER=${SCHEDULER}
    export ACCNR=${ACCNR}
    export QUEUE=${COMPILE_QUEUE}
    export PARTITION=${PARTITION}
    export ROCOTO=${ROCOTO}
    export ECFLOW=${ECFLOW}
    export REGRESSIONTEST_LOG=${REGRESSIONTEST_LOG}
    export LOG_DIR=${LOG_DIR}
EOF

    if [[ $ROCOTO == true ]]; then
      rocoto_create_compile_task
    elif [[ $ECFLOW == true ]]; then
      ecflow_create_compile_task
    else
      ./run_compile.sh ${PATHRT} ${RUNDIR_ROOT} "${MAKE_OPT}" ${COMPILE_NR} > ${LOG_DIR}/compile_${COMPILE_NR}.log 2>&1
    fi

    RT_SUFFIX=""
    BL_SUFFIX=""

    continue

  elif [[ $line == RUN* ]] ; then

    TEST_NAME=$(echo $line | cut -d'|' -f2 | sed -e 's/^ *//' -e 's/ *$//')
    MACHINES=$( echo $line | cut -d'|' -f3 | sed -e 's/^ *//' -e 's/ *$//')
    CB=$(       echo $line | cut -d'|' -f4)
    DEP_RUN=$(  echo $line | cut -d'|' -f5 | sed -e 's/^ *//' -e 's/ *$//')
    DATE_35D=$( echo $line | cut -d'|' -f6 | sed -e 's/^ *//' -e 's/ *$//')

    [[ -e "tests/$TEST_NAME" ]] || die "run test file tests/$TEST_NAME does not exist"
    [[ $CREATE_BASELINE == true && $CB != *fv3* ]] && continue

    if [[ ${MACHINES} != '' ]]; then
      if [[ ${MACHINES} == -* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] && continue
      elif [[ ${MACHINES} == +* ]]; then
        [[ ${MACHINES} =~ ${MACHINE_ID} ]] || continue
      else
        echo "MACHINES=|${MACHINES}|"
        die "MACHINES spec must be either an empty string or start with either '+' or '-'"
      fi
    fi

    # 35 day tests
    [[ $TEST_35D == true ]] && rt_35d

    # Avoid uninitialized RT_SUFFIX/BL_SUFFIX (see definition above)
    RT_SUFFIX=${RT_SUFFIX:-""}
    BL_SUFFIX=${BL_SUFFIX:-""}

    if [[ $ROCOTO == true && $new_compile == true ]]; then
      new_compile=false
      in_metatask=true
      cat << EOF >> $ROCOTO_XML
  <metatask name="compile_${COMPILE_NR}_tasks"><var name="zero">0</var>
EOF
    fi

    TEST_NR=$( printf '%03d' $(( 10#$TEST_NR + 1 )) )

    (
      source ${PATHRT}/tests/$TEST_NAME

      compute_petbounds_and_tasks

      TPN=$(( TPN / THRD ))
      NODES=$(( TASKS / TPN ))
      if (( NODES * TPN < TASKS )); then
        NODES=$(( NODES + 1 ))
      fi

      cat << EOF > ${RUNDIR_ROOT}/run_test_${TEST_NR}.env
      export JOB_NR=${JOB_NR}
      export MACHINE_ID=${MACHINE_ID}
      export RT_COMPILER=${RT_COMPILER}
      export RTPWD=${RTPWD}
      export INPUTDATA_ROOT=${INPUTDATA_ROOT}
      export INPUTDATA_ROOT_WW3=${INPUTDATA_ROOT_WW3}
      export INPUTDATA_ROOT_BMIC=${INPUTDATA_ROOT_BMIC}
      export PATHRT=${PATHRT}
      export PATHTR=${PATHTR}
      export NEW_BASELINE=${NEW_BASELINE}
      export CREATE_BASELINE=${CREATE_BASELINE}
      export RT_SUFFIX=${RT_SUFFIX}
      export BL_SUFFIX=${BL_SUFFIX}
      export SCHEDULER=${SCHEDULER}
      export ACCNR=${ACCNR}
      export QUEUE=${QUEUE}
      export PARTITION=${PARTITION}
      export ROCOTO=${ROCOTO}
      export ECFLOW=${ECFLOW}
      export REGRESSIONTEST_LOG=${REGRESSIONTEST_LOG}
      export LOG_DIR=${LOG_DIR}
      export DEP_RUN=${DEP_RUN}
      export skip_check_results=${skip_check_results}
      export delete_rundir=${delete_rundir}
EOF
      if [[ $MACHINE_ID = jet.* ]]; then
        cat << EOF >> ${RUNDIR_ROOT}/run_test_${TEST_NR}.env
      export PATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/bin:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/bin:$PATH
      export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/lib/python3.8/site-packages
EOF
      fi

      if [[ $ROCOTO == true ]]; then
        rocoto_create_run_task
      elif [[ $ECFLOW == true ]]; then
        ecflow_create_run_task
      else
        ./run_test.sh ${PATHRT} ${RUNDIR_ROOT} ${TEST_NAME} ${TEST_NR} ${COMPILE_NR} > ${LOG_DIR}/run_${TEST_NAME}${RT_SUFFIX}.log 2>&1
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
  echo "endsuite" >> ${ECFLOW_RUN}/${ECFLOW_SUITE}.def
  # run ecflow workflow until done
  ecflow_run
fi

##
## regression test is either failed or successful
##
set +e
cat ${LOG_DIR}/compile_*_time.log              >> ${REGRESSIONTEST_LOG}
cat ${LOG_DIR}/rt_*.log                        >> ${REGRESSIONTEST_LOG}

FILES="fail_test_* fail_compile_*"
for f in $FILES; do
  if [[ -f "$f" ]]; then
    cat "$f" >> fail_test
  fi
done

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

  rm -f fv3_*.x fv3_*.exe modules.fv3_* modulefiles/modules.fv3_* keep_tests.tmp
  [[ ${KEEP_RUNDIR} == false ]] && rm -rf ${RUNDIR_ROOT}
  [[ ${ROCOTO} == true ]] && rm -f ${ROCOTO_XML} ${ROCOTO_DB} ${ROCOTO_STATE} *_lock.db
  [[ ${TEST_35D} == true ]] && rm -f tests/cpld_bmark*_20*
  [[ ${SINGLE_NAME} != '' ]] && rm -f $RT_SINGLE_CONF
fi

date >> ${REGRESSIONTEST_LOG}

elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )
echo "Elapsed time: ${elapsed_time}. Have a nice day!" >> ${REGRESSIONTEST_LOG}
echo "Elapsed time: ${elapsed_time}. Have a nice day!"
