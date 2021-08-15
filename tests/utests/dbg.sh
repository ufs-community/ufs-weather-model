set -eu
source $PATHRT/utests/std.sh

if [[ $application == 'global' ]]; then
  LIST_FILES="sfcf000.nc sfcf001.nc atmf000.nc atmf001.nc"
  FHMAX=1
elif [[ $application == 'regional' ]]; then
  echo "Regional application not yet implemented for debug"
  exit 1
elif [[ $application == 'cpld' ]]; then
  FHMAX=6
  DAYS=0.25
  NFHOUT_HF=1
  RESTART_INTERVAL=${FHMAX}
  RESTART_N=${FHMAX}
  LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/sfcf024/sfcf006/g" \
                                   | sed -E "s/atmf024/atmf006/g" \
                                   | sed -E "s/2016-10-04-00000/2016-10-03-21600/g" \
                                   | sed -e "s/^ *//" -e "s/ *$//")
fi

source $PATHRT/utests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/unit_test${RT_SUFFIX}.env
export WLCLK=${WLCLK}
EOF
