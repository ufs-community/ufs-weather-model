set -eu

MPI_PROC_BIND="-bind-to user:"
export OMP_ENV=${OMP_ENV:-""}
if [[ $rc != 'thr' ]]; then
  for i in $(seq 0 $((TASKS-1))); do
    MPI_PROC_BIND="$MPI_PROC_BIND$i,"
  done
elif [[ $rc == 'thr' ]]; then
  OMP_ENV="OMP_PROC_BIND=true OMP_NUM_THREADS=$THRD"
  if [[ $application != 'regional' ]]; then
    atm_total=$((INPES*JNPES*6))
  elif [[ $application == 'regional' ]]; then
    atm_total=$((INPES*JNPES))
  fi

  for i in $(seq 0 $((atm_total-1))); do
    MPI_PROC_BIND="$MPI_PROC_BIND$((i*2))+$((i*2+1)),"
  done
  for i in $(seq $((atm_total*2)) $((TASKS+atm_total-1))); do
    MPI_PROC_BIND="$MPI_PROC_BIND$i,"
  done
fi

MPI_PROC_BIND=$(echo -n $MPI_PROC_BIND | sed -e "s/,$//")
