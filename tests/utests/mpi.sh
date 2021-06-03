set -eu
source $PATHRT/utests/std.sh

if [[ $application == 'global' ]]; then
  if [ $CI_TEST == 'true' ]; then
    INPES=2
    JNPES=2
  else
    JNPES=$((JNPES/2))
  fi
  WRITE_GROUP=2
  WRTTASK_PER_GROUP=12
  TASKS=$(( INPES*JNPES*6 + WRITE_GROUP*WRTTASK_PER_GROUP ))
  NODES=$(((TASKS+TPN-1)/TPN))
elif [[ $application == 'regional' ]]; then
  echo "Regional application not yet implemented for mpi"
  exit 1
elif [[ $application == 'cpld' ]]; then
  echo "Coupled application not yet implemented for mpi"
  exit 1
fi

(test $CI_TEST == 'true') && source $PATHRT/utests/cmp_proc_bind.sh
source $PATHRT/utests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/unit_test${RT_SUFFIX}.env
export WRITE_GROUP=${WRITE_GROUP}
export WRTTASK_PER_GROUP=${WRTTASK_PER_GROUP}
EOF
