#!/bin/bash
set -eux

source rt_utils.sh
source atparse.bash
source edit_inputs.sh

mkdir -p ${RUNDIR}
cd $RUNDIR

###############################################################################
# Make configure and run files
###############################################################################

# executable:
cp ${PATHRT}/$FV3X datm_mom6_cice.exe

# modulefile for prerequisites:
#cp ${PATHRT}/modules.datm_mom6_cice_${COMPILE_NR} modules.datm
cp ${PATHRT}/modules.fcst modules.datm

# Get the shell file that loads the "module" command and purges modules:
cp ${PATHRT}/../NEMS/src/conf/module-setup.sh.inc  module-setup.sh

SRCD="${PATHTR}"
RUND="${RUNDIR}"

# Set up the run directory
atparse < ${PATHRT}/datm_conf/${FV3_RUN} > datm_run
source ./datm_run
atparse < ${PATHRT}/parm/input.mom6.nml.IN > input.nml
atparse < ${PATHRT}/parm/model_configure.IN > model_configure
atparse < ${PATHRT}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure

edit_ice_in < ${PATHRT}/parm/ice_in_template > ice_in
edit_diag_table < ${PATHRT}/parm/diag_table_template > diag_table
cp ${PATHRT}/parm/data_table data_table
cp ${PATHRT}/parm/datm_data_table.IN datm_data_table

if [[ $MEDCOMP != '' ]]; then
cp ${PATHRT}/parm/fd_nems.yaml fd_nems.yaml
cp ${PATHRT}/parm/pio_in pio_in
cp ${PATHRT}/parm/med_modelio.nml med_modelio.nml
fi

if [[ $OCNRES = '025' ]]; then
edit_mom_input < ${PATHRT}/parm/MOM_input_template > INPUT/MOM_input
elif [[ $OCNRES = '100' ]]; then
edit_mom_input < ${PATHRT}/parm/MOM_input_template100 > INPUT/MOM_input
fi

if [[ $SCHEDULER = 'moab' ]]; then
  atparse < $PATHRT/datm_conf/datm_msub.IN > job_card
elif [[ $SCHEDULER = 'pbs' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/datm_conf/datm_qsub.IN > job_card
elif [[ $SCHEDULER = 'sbatch' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/datm_conf/datm_qsub.IN > job_card
elif [[ $SCHEDULER = 'slurm' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/datm_conf/datm_slurm.IN > job_card
elif [[ $SCHEDULER = 'lsf' ]]; then
  if (( TASKS < TPN )); then
    TPN=${TASKS}
  fi
  atparse < $PATHRT/datm_conf/datm_bsub.IN > job_card
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
