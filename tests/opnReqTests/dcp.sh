set -eu
source $PATHRT/opnReqTests/std.sh

if [[ $application == 'global' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    temp=$INPES
    INPES=$JNPES
    JNPES=$temp
  else
    INPES=6
    JNPES=4
  fi
elif [[ $application == 'regional' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=10
    JNPES=3
    NTILES=1
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
    NODES=$(((TASKS+TPN-1)/TPN))
  else
    if [[ $TEST_NAME == 'regional_control' ]]; then
      INPES=5
      JNPES=12
    elif [[ $TEST_NAME =~ 'regional_3km' ]]; then
      temp=$INPES
      INPES=$JNPES
      JNPES=$temp
    fi
  fi
elif [[ $application == 'cpld' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=3
    JNPES=1
    OCN_tasks=10
    ICE_tasks=6
    NPROC_ICE=$ICE_tasks
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP + OCN_tasks + ICE_tasks))
  else
    temp=$INPES
    INPES=$JNPES
    JNPES=$temp
  fi
fi

(test $CI_TEST == 'true') && source $PATHRT/opnReqTests/cmp_proc_bind.sh
source $PATHRT/opnReqTests/wrt_env.sh
