set -eu
source $PATHRT/utests/std.sh

# Set up date and time of restart files for restart run
FHROT=12
RESTART_FILE_PREFIX="${SYEAR}${SMONTH}${SDAY}.$(printf "%02d" $(( SHOUR + FHROT  )))0000"

if [[ $application == 'global' ]]; then
  OUTPUT_FH="3 -1"
elif [[ $application == 'regional' ]]; then
  echo "Regional application not yet implemented for restart"
  exit 1
elif [[ $application == 'cpld' ]]; then
  CICERUNTYPE='continue'
  RUNTYPE='continue'
  USE_RESTART_TIME='.true.'
  MOM6_RESTART_SETTING="r"
  RESTART_FILE_SUFFIX_HRS="${SYEAR}-${SMONTH}-${SDAY}-$(printf "%02d" $(( ${FHROT} )))"
  RESTART_FILE_SUFFIX_SECS="${SYEAR}-${SMONTH}-${SDAY}-$(printf "%02d" $(( ${FHROT}*3600 )))"
  RESTART_N=$((FHMAX-$FHROT))
  DEP_RUN=${TEST_NAME}
fi

WARM_START=.T.
NGGPS_IC=.F.
EXTERNAL_IC=.F.
MAKE_NH=.F.
MOUNTAIN=.T.
NA_INIT=0

LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/phyf00(00|21)\.(tile.\.nc|nemsio|nc) ?//g" \
                                 | sed -E "s/dynf00(00|21)\.(tile.\.nc|nemsio|nc) ?//g" \
                                 | sed -E "s/sfcf0(00|21).nc ?//g" | sed -E "s/atmf0(00|21).nc ?//g" \
                                 | sed -E "s/GFSFLX.GrbF(00|21) ?//g" | sed -E "s/GFSPRS.GrbF(00|21) ?//g" \
                                 | sed -E "s/atmos_4xdaily\.tile[1-6]\.nc ?//g" | sed -e "s/^ *//" -e "s/ *$//")

(test $CI_TEST == 'true') && source $PATHRT/utests/cmp_proc_bind.sh
source $PATHRT/utests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/unit_test${RT_SUFFIX}.env
export FHROT=${FHROT}
export RESTART_FILE_PREFIX=${RESTART_FILE_PREFIX}
export NSTF_NAME=${NSTF_NAME}
export CICERUNTYPE=${CICERUNTYPE:-}
export RUNTYPE=${RUNTYPE:-}
export USE_RESTART_TIME=${USE_RESTART_TIME:-}
export MOM6_RESTART_SETTING=${MOM6_RESTART_SETTING:-}
export RESTART_FILE_SUFFIX_HRS=${RESTART_FILE_SUFFIX_HRS:-}
export RESTART_FILE_SUFFIX_SECS=${RESTART_FILE_SUFFIX_SECS:-}
export DEP_RUN=${DEP_RUN:-}
export WARM_START=${WARM_START}
export NGGPS_IC=${NGGPS_IC}
export EXTERNAL_IC=${EXTERNAL_IC}
export MAKE_NH=${MAKE_NH}
export MOUNTAIN=${MOUNTAIN}
export NA_INIT=${NA_INIT}
export OUTPUT_FH="${OUTPUT_FH:-}"
EOF
