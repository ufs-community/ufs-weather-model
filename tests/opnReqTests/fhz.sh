set -eu
source $PATHRT/opnReqTests/std.sh

FHZERO=3

source $PATHRT/opnReqTests/wrt_env.sh

cat <<EOF >>${RUNDIR_ROOT}/opnreq_test${RT_SUFFIX}.env
export FHZERO=${FHZERO}
EOF
