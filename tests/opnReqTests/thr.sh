set -eu
source $PATHRT/opnReqTests/std.sh

THRD=2
TPN=$((TPN/THRD))
if [[ $application == 'global' ]]; then
  JNPES=$((JNPES/THRD))
  TASKS=$((INPES*JNPES*6 + WRITE_GROUP*WRTTASK_PER_GROUP))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'regional' ]]; then
  if [[ $CI_TEST == 'true' ]]; then
    INPES=4
    JNPES=4
    TASKS=$((INPES*JNPES + WRITE_GROUP*WRTTASK_PER_GROUP))
  fi
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'cpld' ]]; then
  if [[ $CI_TEST != 'true' ]]; then
    if [[ $TEST_NAME == 'cpld_control_c96_p8' ]]; then
      INPES=3
      JNPES=4
      med_petlist_bounds="0 71"
      chm_petlist_bounds="0 71"
      atm_petlist_bounds="0 77"
      ocn_petlist_bounds="78 107"
      ice_petlist_bounds="108 119"
      TASKS=$((INPES*JNPES*6 + WRITE_GROUP*WRTTASK_PER_GROUP + 30 + 12))
      NODES=$(((TASKS+TPN-1)/TPN))
    elif [[ $TEST_NAME == 'cpld_control_c96_noaero_p8' ]]; then
      INPES=3
      JNPES=4
      med_petlist_bounds="0 71"
      atm_petlist_bounds="0 77"
      ocn_petlist_bounds="78 107"
      ice_petlist_bounds="108 119"
      TASKS=$((INPES*JNPES*6 + WRITE_GROUP*WRTTASK_PER_GROUP + 30 + 12))
      NODES=$(((TASKS+TPN-1)/TPN))
    elif [[ $TEST_NAME == 'cpld_control_p8' ]]; then
      INPES=3
      JNPES=4
      med_petlist_bounds="0 71"
      chm_petlist_bounds="0 71"
      atm_petlist_bounds="0 77"
      ocn_petlist_bounds="78 97"
      ice_petlist_bounds="98 107"
      wav_petlist_bounds="108 119"
      TASKS=$((INPES*JNPES*6 + WRITE_GROUP*WRTTASK_PER_GROUP + 20 + 10 + 12))
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
