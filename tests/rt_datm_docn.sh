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

atparse < ${PATHRT}/cdeps_conf/${DATM_RUN:-datm_run.IN} > datm_run
atparse < ${PATHRT}/cdeps_conf/${DOCN_RUN:-docn_run.IN} > docn_run

# Set up the run directory
source ./datm_run
source ./docn_run

if [[ $SCHEDULER = 'pbs' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_qsub.IN > job_card
elif [[ $SCHEDULER = 'slurm' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_slurm.IN > job_card
elif [[ $SCHEDULER = 'lsf' ]]; then
  if (( TASKS < TPN )); then
    TPN=${TASKS}
  fi
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_bsub.IN > job_card
fi

# Driver
atparse < ${PATHTR}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure
atparse < ${PATHTR}/parm/${MODEL_CONFIGURE:-model_configure.cdeps.IN} > model_configure

# CMEPS
cp ${PATHTR}/parm/fd_nems.yaml fd_nems.yaml
cp ${PATHTR}/parm/pio_in pio_in
cp ${PATHTR}/parm/med_modelio.nml med_modelio.nml

# CDEPS, DATM and DOCN
cp ${PATHTR}/parm/atm_modelio.nml atm_modelio.nml
cp ${PATHTR}/parm/ocn_modelio.nml ocn_modelio.nml
atparse < ${PATHTR}/parm/${DATM_CONFIGURE_A:-datm_in} > datm_in
atparse < ${PATHTR}/parm/${DATM_CONFIGURE_B:-datm.streams.xml} > datm.streams.xml
atparse < ${PATHTR}/parm/${DOCN_CONFIGURE_A:-docn_in} > docn_in
atparse < ${PATHTR}/parm/${DOCN_CONFIGURE_B:-docn.streams.xml} > docn.streams.xml

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
