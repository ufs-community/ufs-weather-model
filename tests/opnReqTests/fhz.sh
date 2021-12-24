set -eu
source $PATHRT/opnReqTests/std.sh

FHZERO=3

FHMAX_2D=$(printf "%02d" $FHMAX)
LIST_FILES=$(echo -n $LIST_FILES | sed -E "s/phyf024.nc ?//g" \
                                 | sed -E "s/dynf024.nc ?//g" \
                                 | sed -E "s/sfcf024.nc ?//g" \
                                 | sed -E "s/atmf024.nc ?//g" \
                                 | sed -E "s/GFSFLX.GrbF24 ?//g" \
                                 | sed -E "s/GFSPRS.GrbF24 ?//g" \
                                 | sed -e "s/^ *//" -e "s/ *$//")

source $PATHRT/opnReqTests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export FHZERO=${FHZERO}
EOF
