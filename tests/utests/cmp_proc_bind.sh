MPI_PROC_BIND="-bind-to user:"
for i in $(seq 0 $((TASKS-1))); do
  MPI_PROC_BIND="$MPI_PROC_BIND$i,"
done
MPI_PROC_BIND=$(echo -n $MPI_PROC_BIND | sed -e "s/,$//")
