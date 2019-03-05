### head.h start

set -e # stop the shell on first error
set -u # fail when using an undefined variable
set -x # echo script lines as they are executed


# Defines the variables that are needed for any communication with ECF
export ECF_PORT=%ECF_PORT%    # The server port number
export ECF_HOST=%ECF_HOST%    # The name of ecf host that issued this task
export ECF_NAME=%ECF_NAME%    # The name of this current task
export ECF_PASS=%ECF_PASS%    # A unique password
export ECF_TRYNO=%ECF_TRYNO%  # Current try number of the task
export ECF_RID=$$             # record the process id. Also used for zombie detection

# Define the path where to find ecflow_client
# make sure client and server use the *same* version.
# Important when there are multiple versions of ecFlow
#export PATH=....:$PATH

# Tell ecFlow we have started
ecflow_client --init=$$


# Define a error handler
ERROR() {
   set +e                      # Clear -e flag, so we don't fail
   kill $(jobs -p)
   wait                        # wait for background process to stop

   ecflow_client --ping --host=${ECF_HOST} --port=${ECF_PORT}
   not_running=$?
   if [[ $not_running -eq 0 ]]; then
     export ECF_TIMEOUT=5
     ecflow_client --abort=trap  # Notify ecFlow that something went wrong, using 'trap' as the reason
   fi
   sleep 5
   trap 0                      # Remove the trap
   exit 0                      # End the script
}


# Trap any calls to exit and errors caught by the -e flag
trap ERROR 0

# Trap any signal that may cause the script to fail
trap '{ echo "$0 Killed by a signal"; ERROR ; }' 1 2 3 4 5 6 7 8   10    12 13    15

### head.h end
