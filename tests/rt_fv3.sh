#!/bin/bash
set -eux

source rt_utils.sh
source atparse.bash

mkdir -p ${RUNDIR}
cd $RUNDIR

###############################################################################
# Make configure and run files
###############################################################################

# FV3 executable:
cp ${PATHRT}/$FV3X                                 fv3.exe

# modulefile for FV3 prerequisites:
cp ${PATHRT}/modules.fv3_${COMPILE_NR}             modules.fv3

# Get the shell file that loads the "module" command and purges modules:
cp ${PATHRT}/../NEMS/src/conf/module-setup.sh.inc  module-setup.sh

SRCD="${PATHTR}"
RUND="${RUNDIR}"

atparse < ${PATHRT}/fv3_conf/${FV3_RUN:-fv3_run.IN} > fv3_run

atparse < ${PATHTR}/parm/${INPUT_NML:-input.nml.IN} > input.nml

atparse < ${PATHTR}/parm/${MODEL_CONFIGURE:-model_configure.IN} > model_configure

if [[ "Q${INPUT_NEST02_NML:-}" != Q ]] ; then
    atparse < ${PATHTR}/parm/${INPUT_NEST02_NML} > input_nest02.nml
fi

source ./fv3_run

if [[ $SCHEDULER = 'moab' ]]; then
  atparse < $PATHRT/fv3_conf/fv3_msub.IN > job_card
elif [[ $SCHEDULER = 'pbs' ]]; then
  NODES=$(( TASKS / TPN ))
  atparse < $PATHRT/fv3_conf/fv3_qsub.IN > job_card
elif [[ $SCHEDULER = 'lsf' ]]; then
  atparse < $PATHRT/fv3_conf/fv3_bsub.IN > job_card
fi

################################################################################
# Submit test
################################################################################

if [[ $ROCOTO = 'false' ]]; then
  submit_and_wait job_card
else
  chmod u+x job_card
  ./job_card
fi

check_results

################################################################################
# End test
################################################################################

exit 0
