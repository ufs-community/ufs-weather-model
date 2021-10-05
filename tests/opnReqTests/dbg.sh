set -eu
source $PATHRT/opnReqTests/std.sh

if [[ $application == 'global' ]]; then
  FHMAX=1
  OUTPUT_FH="0 1"
  LIST_FILES="sfcf000.nc sfcf001.nc atmf000.nc atmf001.nc"
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
                                   | sed -E "s/2021-03-23-21600/2021-03-22-43200/g" \
                                   | sed -e "s/^ *//" -e "s/ *$//")
fi

source $PATHRT/opnReqTests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export WLCLK=${WLCLK}
EOF
