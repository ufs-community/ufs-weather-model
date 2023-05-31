set -eu

if [[ $application == 'global' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=3
    JNPES=2
    WRITE_GROUP=1
    WRTTASK_PER_GROUP=12
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
    NODES=$(((TASKS+TPN-1)/TPN))
  fi
  RESTART_N=$(( FHMAX/2 ))
  RESTART_INTERVAL="${RESTART_N} -1"
elif [[ $application == 'regional' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=10
    JNPES=11
    NTILES=1
    WRTTASK_PER_GROUP=10
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
    NODES=$(((TASKS+TPN-1)/TPN))
  fi
elif [[ $application == 'cpld' ]]; then
  if [ $CI_TEST == 'true' ]; then
    INPES=2
    JNPES=2
    OCN_tasks=10
    ICE_tasks=6
    NPROC_ICE=$ICE_tasks
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP + OCN_tasks + ICE_tasks))
    NODES=$(((TASKS+TPN-1)/TPN))
  fi
  RESTART_N=$(( FHMAX/2 ))
  RESTART_INTERVAL="${RESTART_N} -1"
elif [[ $application == 'atmw' ]]; then
  if [ $CI_TEST == 'true' ]; then
    INPES=3
    JNPES=8
    WRTTASK_PER_GROUP=6
    TASKS=$((INPES*JNPES + WRITE_GROUP*WRTTASK_PER_GROUP))
    NODES=$(((TASKS+TPN-1)/TPN))
  fi
  RESTART_N=$(( FHMAX/2 ))
  RESTART_INTERVAL="${RESTART_N} -1"
  WW3RSTDTHR=3
  DT_2_RST="$(printf "%02d" $(( ${WW3RSTDTHR}*3600 )))"
fi


#outdated (test $CI_TEST == 'true') && source $PATHRT/opnReqTests/cmp_proc_bind.sh
if [[ $RT_SUFFIX =~ std ]]; then
  source $PATHRT/opnReqTests/wrt_env.sh
fi

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export DT_2_RST=${DT_2_RST:-}
EOF
