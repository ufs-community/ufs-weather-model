
### tail.h start
wait                      # wait for background process to stop

ecflow_client --ping --host=${ECF_HOST} --port=${ECF_PORT}
not_running=$?
if [[ $not_running -eq 0 ]]; then
  ecflow_client --complete  # Notify ecFlow of a normal end
fi

trap 0                    # Remove all traps
exit 0                    # End the shell
### tail.h end
