set -eu
source $PATHRT/opnReqTests/std.sh

if [[ $application == 'global' ]]; then
  if [ $CI_TEST == 'true' ]; then
    INPES=2
    JNPES=2
  else
    JNPES=$((JNPES/2))
  fi
  WRITE_GROUP=2
  WRTTASK_PER_GROUP=12
  TASKS=$(( INPES*JNPES*NTILES + WRITE_GROUP*WRTTASK_PER_GROUP ))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'regional' ]]; then
  echo "Regional application not yet implemented for mpi"
  exit 1
elif [[ $application == 'cpld' ]]; then
  echo "Coupled application not yet implemented for mpi"
  exit 1
fi

(test $CI_TEST == 'true') && source $PATHRT/opnReqTests/cmp_proc_bind.sh
source $PATHRT/opnReqTests/wrt_env.sh
