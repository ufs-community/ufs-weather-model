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
elif [[ $application == 'regional' ]]; then
  FHMAX=1
  OUTPUT_FH="0 1"

  FHMAX_2D=$(printf "%02d" $FHMAX)
  LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/phyf0[0-9][1-9]/phyf0$FHMAX_2D/g" \
                                   | sed -E "s/dynf0[0-9][1-9]/dynf0$FHMAX_2D/g" \
                                   | sed -E "s/PRSLEV.GrbF[0-9][0-9]//g" \
                                   | sed -E "s/NATLEV.GrbF[0-9][0-9]//g" \
                                   | sed -e "s/^ *//" -e "s/ *$//")

  WLCLK=30
  WRITE_DOPOST=.false.
elif [[ $application == 'cpld' ]]; then
  FHMAX=6
  DAYS=0.25
  NFHOUT_HF=1
  RESTART_INTERVAL=${FHMAX}
  RESTART_N=${FHMAX}
  OUTPUT_FH="0 ${FHMAX}"
  LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/sfcf024/sfcf006/g" \
                                   | sed -E "s/atmf024/atmf006/g" \
                                   | sed -E "s/2021-03-23-21600/2021-03-22-43200/g" \
                                   | sed -E "s/sfcf021.tile[1-6].nc ?//g" \
                                   | sed -E "s/atmf021.tile[1-6].nc ?//g" \
                                   | sed -E "s/(gocart.inst_aod.202103)23_0600z.nc4/\122_1200z.nc4/g" \
                                   | sed -e "s/^ *//" -e "s/ *$//")
fi

source $PATHRT/opnReqTests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export WLCLK=${WLCLK}
EOF
