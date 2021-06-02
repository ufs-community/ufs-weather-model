set -eu
source $PATHRT/utests/std.sh

THRD=2
if [[ $application == 'global' ]]; then
  JNPES=$((JNPES/THRD))
  TASKS=$((INPES*JNPES*6 + WRITE_GROUP*WRTTASK_PER_GROUP))
  TPN=$((TPN/THRD))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'regional' ]]; then
  TPN=$((TPN/THRD))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'cpld' ]]; then
  if [[ $CI_TEST != 'true' ]]; then
    INPES=3
    JNPES=4
    med_petlist_bounds="0 71"
    atm_petlist_bounds="0 77"
    ocn_petlist_bounds="78 107"
    ice_petlist_bounds="108 119"
    TASKS=$((INPES*JNPES*6 + WRTIE_GROUP*WRTTASK_PER_GROUP + 30 + 12))
    TPN=$((TPN/THRD))
    NODES=$(((TASKS+TPN-1)/TPN))
  fi
fi

(test $CI_TEST == 'true') && source $PATHRT/utests/cmp_proc_bind.sh
source $PATHRT/utests/wrt_env.sh
