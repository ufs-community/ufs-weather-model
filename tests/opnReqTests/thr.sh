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
    INPES=5
    JNPES=11
    NTILES=1
    WRTTASK_PER_GROUP=10
    TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP))
    NODES=$(((TASKS+TPN-1)/TPN))
  fi
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
    elif [[ $TEST_NAME =~ 'cpld_control_p8' ]] || [[ $TEST_NAME =~ 'cpld_control_ciceC_p8' ]] || [[ $TEST_NAME =~ 'cpld_control_gfsv17' ]]; then
      INPES=3
      JNPES=4
      OCN_tasks=20
      ICE_tasks=10
      WAV_tasks=12
      NPROC_ICE=$ICE_tasks
      TASKS=$((INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP + OCN_tasks + ICE_tasks + WAV_tasks))
      NODES=$(((TASKS+TPN-1)/TPN))
    elif [[ $TEST_NAME == 'cpld_bmark_p8' ]]; then
      #need to overhaul NODES=$(((TASKS+TPN-1)/TPN))
      echo $TEST_NAME
    else
      echo "This test is not yet set up for the thread test, skipping..."
      continue 1
    fi
  fi
elif [[ $application == 'atmw' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=5
    JNPES=11
    WRTTASK_PER_GROUP=10
    TASKS=$((INPES*JNPES + WRITE_GROUP*WRTTASK_PER_GROUP))
    NODES=$(((TASKS+TPN-1)/TPN))
  fi
fi


source $PATHRT/opnReqTests/wrt_env.sh
