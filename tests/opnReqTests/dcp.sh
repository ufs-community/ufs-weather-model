set -eu
source $PATHRT/opnReqTests/std.sh

if [[ $application == 'global' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    temp=$INPES
    INPES=$JNPES
    JNPES=$temp
  else
    temp=$INPES
    INPES=$JNPES
    JNPES=$temp
  fi
elif [[ $application == 'regional' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=5
    JNPES=12
    NTILES=1
    WRTTASK_PER_GROUP=10
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
    NODES=$(((TASKS+TPN-1)/TPN))
  else
    INPES=5
    JNPES=12
    NTILES=1
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


source $PATHRT/opnReqTests/wrt_env.sh
