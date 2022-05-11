set -eu
source $PATHRT/opnReqTests/std.sh

THRD=2
TPN=$((TPN/THRD))
if [[ $application == 'global' ]]; then
  JNPES=$((JNPES/THRD))
  TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'regional' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=4
    JNPES=4
    NTILES=1
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
  fi
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'cpld' ]]; then
  if [[ $CI_TEST != 'true' ]]; then
    if [[ $TEST_NAME =~ 'cpld_control_c96_p8' ]]; then
      INPES=3
      JNPES=4
      OCN_tasks=30
      ICE_tasks=12
      NPROC_ICE=$ICE_tasks
      TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP + OCN_tasks + ICE_tasks))
      NODES=$(((TASKS+TPN-1)/TPN))
    elif [[ $TEST_NAME =~ 'cpld_control_c96_noaero_p8' ]]; then
      INPES=3
      JNPES=4
      OCN_tasks=30
      ICE_tasks=12
      NPROC_ICE=$ICE_tasks
      TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP + OCN_tasks + ICE_tasks))
      NODES=$(((TASKS+TPN-1)/TPN))
    elif [[ $TEST_NAME =~ 'cpld_control_p8' ]]; then
      INPES=3
      JNPES=4
      OCN_tasks=20
      ICE_tasks=10
      WAV_tasks=12
      NPROC_ICE=$ICE_tasks
      TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP + OCN_tasks + ICE_tasks + WAV_tasks))
      NODES=$(((TASKS+TPN-1)/TPN))
    elif [[ $TEST_NAME == 'cpld_bmark_p8' ]]; then
      NODES=$(((TASKS+TPN-1)/TPN))
    else
      echo "This test is not yet set up for the thread test"
      exit 1
    fi
  fi
fi

(test $CI_TEST == 'true') && source $PATHRT/opnReqTests/cmp_proc_bind.sh
source $PATHRT/opnReqTests/wrt_env.sh
