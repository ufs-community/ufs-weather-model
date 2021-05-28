set -eu
source $PATHRT/utests/std.sh

THRD=2
if [[ $application == 'global' ]]; then
  JNPES=$((JNPES/THRD))
  TASKS=$((INPES*JNPES*6 + WRITE_GROUP*WRTTASK_PER_GROUP))
  TPN=$((TPN/THRD))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'regional' ]]; then
  #JNPES=$((JNPES/THRD))
  #TASKS=$((INPES*JNPES))
  TPN=$((TPN/THRD))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'cpld' ]]; then
  JNPES=$((JNPES/THRD))
  TASKS=$((INPES*JNPES*6 + WRTIE_GROUP*WRTTASK_PER_GROUP + 10 + 6))
  TPN=$((TPN/THRD))
  NODES=$(((TASKS+TPN-1)/TPN))
fi

(test $CI_TEST == 'true') && source $PATHRT/utests/cmp_proc_bind.sh
source $PATHRT/utests/wrt_env.sh
