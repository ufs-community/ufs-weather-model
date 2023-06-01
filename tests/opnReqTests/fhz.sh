set -eu
source $PATHRT/opnReqTests/std.sh

if [[ $application == 'global' ]]; then
  FHZERO=3
  LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/phyf024.nc ?//g" \
                                   | sed -E "s/dynf024.nc ?//g" \
                                   | sed -E "s/sfcf024.nc ?//g" \
                                   | sed -E "s/atmf024.nc ?//g" \
                                   | sed -E "s/GFSFLX.GrbF24 ?//g" \
                                   | sed -E "s/GFSPRS.GrbF24 ?//g" \
                                   | sed -e "s/^ *//" -e "s/ *$//")
elif [[ $application == 'cpld' ]]; then
  if [[ $TEST_NAME =~ 'cpld_control_c96_p8' ]] || [[ $TEST_NAME =~ 'cpld_control_p8' ]] || [[ $TEST_NAME =~ 'cpld_control_c96_noaero_p8' ]] || [[ $TEST_NAME =~ 'cpld_control_ciceC_p8' ]]; then
    FHZERO=3
    LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/sfcf024.tile[1-6].nc ?//g" \
                                     | sed -E "s/atmf024.tile[1-6].nc ?//g" \
                                     | sed -e "s/^ *//" -e "s/ *$//")
  else
    echo "This test is not yet set up for the fhzero test, skipping..."
    continue 1
  fi
else
  echo "This test is not yet set up for the fhzero test, skipping..."
  continue 1
fi

source $PATHRT/opnReqTests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export FHZERO=${FHZERO}
EOF
