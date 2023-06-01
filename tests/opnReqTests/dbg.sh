set -eu
source $PATHRT/opnReqTests/std.sh

if [[ $application == 'global' ]]; then
  FHMAX=1
  OUTPUT_FH="0 1"

  FHMAX_2D=$(printf "%02d" $FHMAX)
  if [[ $TEST_NAME == 'fv3_gsd' ]]; then
    LIST_FILES="sfcf0$FHMAX_2D.tile1.nc sfcf0$FHMAX_2D.tile2.nc sfcf0$FHMAX_2D.tile3.nc \
                sfcf0$FHMAX_2D.tile4.nc sfcf0$FHMAX_2D.tile5.nc sfcf0$FHMAX_2D.tile6.nc \
                atmf0$FHMAX_2D.tile1.nc atmf0$FHMAX_2D.tile2.nc atmf0$FHMAX_2D.tile3.nc \
                atmf0$FHMAX_2D.tile4.nc atmf0$FHMAX_2D.tile5.nc atmf0$FHMAX_2D.tile6.nc"
  else
    LIST_FILES="sfcf0$FHMAX_2D.nc sfcf0$FHMAX_2D.nc"
  fi
  WRITE_DOPOST=.false.
elif [[ $application == 'regional' ]]; then
  echo "Regional application not yet implemented for debug, skipping..."
  continue 1
elif [[ $application == 'cpld' ]]; then
  FHMAX=3
  DAYS=0.125
  NFHOUT_HF=1
  RESTART_INTERVAL=${FHMAX}
  RESTART_N=${FHMAX}
  OUTPUT_FH="0 ${FHMAX}"
  AOD_FRQ=030000
  LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/sfcf024/sfcf003/g" \
                                   | sed -E "s/atmf024/atmf003/g" \
                                   | sed -E "s/2021-03-23-21600/2021-03-22-32400/g" \
                                   | sed -E "s/20210323\.060000/20210322\.090000/g" \
                                   | sed -E "s/sfcf021\.tile[1-6]\.nc ?//g" \
                                   | sed -E "s/atmf021\.tile[1-6]\.nc ?//g" \
                                   | sed -E "s/(gocart\.inst_aod\.202103)23_0600z\.nc4/\122_0900z\.nc4/g" \
                                   | sed -E "s/20210323\.060000\.out_pnt\.ww3/20210322\.090000\.out_pnt\.ww3/g" \
                                   | sed -E "s/20210323\.060000\.out_grd\.ww3/20210322\.090000\.out_grd\.ww3/g" \
                                   | sed -e "s/^ *//" -e "s/ *$//")
elif [[ $application == 'atmw' ]]; then
  FHMAX=3
  WW3RSTDTHR=3
  DT_2_RST="$(printf "%02d" $(( ${WW3RSTDTHR}*3600 )))"
  DAYS=0.125
  NFHOUT_HF=1
  RESTART_INTERVAL=${FHMAX}
  RESTART_N=${FHMAX}
  OUTPUT_FH="0 ${FHMAX}"
  LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/sfcf012/sfcf003/g" \
                                   | sed -E "s/atmf012/atmf003/g" \
                                   | sed -E "s/2021-03-22-64800/2021-03-22-32400/g" \
                                   | sed -E "s/20210322\.180000/20210322\.090000/g" \
                                   | sed -E "s/20210322\.180000\.out_pnt\.ww3/20210322\.090000\.out_pnt\.ww3/g" \
                                   | sed -E "s/20210322\.180000\.out_grd\.ww3/20210322\.090000\.out_grd\.ww3/g" \
                                   | sed -e "s/^ *//" -e "s/ *$//")
                                   
fi

source $PATHRT/opnReqTests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export WLCLK=${WLCLK}
export DT_2_RST=${DT_2_RST:-}
EOF
