source $PATHRT/utests/std.sh

if [[ $application == 'global' ]]; then
  temp=$INPES
  INPES=$JNPES
  JNPES=$temp
elif [[ $application == 'regional' ]]; then
  JNPES=$((JNPES/2))
  TASKS=$((INPES*JNPES))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'cpld' ]]; then
  echo "Coupled application is not implemented in utest yet"
  exit 1
fi

(test $CI_TEST == 'true') && source $PATHRT/utests/cmp_proc_bind.sh
#source $PATHRT/utests/cmp_proc_bind.sh
source $PATHRT/utests/wrt_env.sh
