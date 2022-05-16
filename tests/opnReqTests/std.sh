set -eu

if [[ $application == 'global' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=3
    JNPES=2
    WRITE_GROUP=1
    WRTTASK_PER_GROUP=12
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
  fi
  RESTART_N=$(( FHMAX/2 ))
  RESTART_INTERVAL="${RESTART_N} -1"
elif [[ $application == 'regional' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=4
    JNPES=6
    NTILES=1
    WRTTASK_PER_GROUP=8
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
  fi
elif [[ $application == 'cpld' ]]; then
  if [ $CI_TEST == 'true' ]; then
    INPES=2
    JNPES=2
    OCN_tasks=10
    ICE_tasks=6
    NPROC_ICE=$ICE_tasks
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP + OCN_tasks + ICE_tasks))
  fi
  RESTART_N=$(( FHMAX/2 ))
  RESTART_INTERVAL="${RESTART_N} -1"
fi

NODES=$(((TASKS+TPN-1)/TPN))
(test $CI_TEST == 'true') && source $PATHRT/opnReqTests/cmp_proc_bind.sh
if [[ $RT_SUFFIX =~ std ]]; then
  source $PATHRT/opnReqTests/wrt_env.sh
fi
