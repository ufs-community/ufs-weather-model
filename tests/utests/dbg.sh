set -eu
source $PATHRT/utests/std.sh

if [[ $application == 'global' ]]; then
  LIST_FILES="sfcf000.nc sfcf001.nc atmf000.nc atmf001.nc"
  FHMAX=1
  WLCLK=60
elif [[ $application == 'regional' ]]; then
  echo "Regional application not yet implemented for debug"
  exit 1
elif [[ $application == 'cpld' ]]; then
  echo "Coupled application not yet implemented for debug"
  exit 1
fi

source $PATHRT/utests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/unit_test${RT_SUFFIX}.env
export WLCLK=${WLCLK}
EOF
