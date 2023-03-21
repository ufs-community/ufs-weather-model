set -eu
source $PATHRT/opnReqTests/std.sh

DEP_RUN=${TEST_NAME}

if [[ ! -z $NSTF_NAME ]]; then
  second_value=$(echo $NSTF_NAME | cut -d ',' -f 2)
  if [[ $second_value -eq 1 ]]; then
    NSTF_NAME=$(echo $NSTF_NAME | awk 'BEGIN{FS=OFS=","} { $2-=1; print}')
  fi
fi

if [[ $application == 'global' ]]; then
  FHROT=$(( FHMAX/2 ))
  OUTPUT_FH="3 -1"
  if [[ $(( SHOUR + FHROT )) -lt 24 ]]; then
    RESTART_FILE_PREFIX="${SYEAR}${SMONTH}$(printf "%02d" ${SDAY}).$(printf "%02d" $(( SHOUR + FHROT  )))0000"
  else
    RESTART_FILE_PREFIX="${SYEAR}${SMONTH}$(printf "%02d" $((SDAY+1))).$(printf "%02d" $(( SHOUR + FHROT - 24 )))0000"
  fi

elif [[ $application == 'regional' ]]; then
  echo "Regional application not yet implemented for restart"
  exit 1
elif [[ $application == 'cpld' ]]; then
  FHROT=$(( FHMAX/2 ))

  CICERUNTYPE='continue'
  RUNTYPE='continue'
  USE_RESTART_TIME='.true.'
  MOM6_RESTART_SETTING="r"
  RESTART_N=$(( FHMAX - FHROT ))
  RESTART_FILE_PREFIX="${SYEAR}${SMONTH}${SDAY}.$(printf "%02d" $(( SHOUR + FHROT  )))0000"
  RESTART_FILE_SUFFIX_SECS="${SYEAR}-${SMONTH}-${SDAY}-$(printf "%05d" $(( (SHOUR + FHROT)* 3600 )))"
  RUN_BEG="${SYEAR}${SMONTH}${SDAY} $(printf "%02d" $(( ${FHROT}+${SHOUR} )))0000"
fi

WARM_START=.T.
NGGPS_IC=.F.
EXTERNAL_IC=.F.
MAKE_NH=.F.
MOUNTAIN=.T.
NA_INIT=0

FHMAX_2D=$(printf "%02d" $FHMAX)
LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/phyf0[0-9][0-9]/phyf0$FHMAX_2D/g" \
                                 | sed -E "s/dynf0[0-9][0-9]/dynf0$FHMAX_2D/g" \
                                 | sed -E "s/sfcf0[0-9][0-9]/sfcf0$FHMAX_2D/g" \
                                 | sed -E "s/atmf0[0-9][0-9]/atmf0$FHMAX_2D/g" \
                                 | sed -E "s/GFSFLX.GrbF[0-9][0-9]/GFSFLX.GrbF$FHMAX_2D/g" \
                                 | sed -E "s/GFSPRS.GrbF[0-9][0-9]/GFSPRS.GrbF$FHMAX_2D/g" \
                                 | sed -E "s/atmos_4xdaily\.tile[1-6]\.nc ?//g" \
                                 | sed -e "s/^ *//" -e "s/ *$//")
LIST_FILES=$(echo $LIST_FILES | xargs -n1 | sort -u | xargs)



source $PATHRT/opnReqTests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export FHROT=${FHROT}
export RESTART_FILE_PREFIX=${RESTART_FILE_PREFIX}
export NSTF_NAME=${NSTF_NAME}
export CICERUNTYPE=${CICERUNTYPE:-}
export RUNTYPE=${RUNTYPE:-}
export USE_RESTART_TIME=${USE_RESTART_TIME:-}
export MOM6_RESTART_SETTING=${MOM6_RESTART_SETTING:-}
export RESTART_N=${RESTART_N:-}
export RESTART_FILE_SUFFIX_SECS=${RESTART_FILE_SUFFIX_SECS:-}
export RUN_BEG="${RUN_BEG:-}"
export OUT_BEG="${RUN_BEG:-}"
export RST_BEG="${RUN_BEG:-}"
export RST_2_BEG="${RUN_BEG:-}"
export DEP_RUN=${DEP_RUN:-}
export WARM_START=${WARM_START}
export NGGPS_IC=${NGGPS_IC}
export EXTERNAL_IC=${EXTERNAL_IC}
export MAKE_NH=${MAKE_NH}
export MOUNTAIN=${MOUNTAIN}
export NA_INIT=${NA_INIT}
EOF
