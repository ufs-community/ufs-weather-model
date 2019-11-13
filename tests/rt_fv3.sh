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
cp ${PATHTR}/parm/post_itag itag
cp ${PATHTR}/parm/postxconfig-NT.txt postxconfig-NT.txt
cp ${PATHTR}/parm/params_grib2_tbl_new params_grib2_tbl_new

SRCD="${PATHTR}"
RUND="${RUNDIR}"

atparse < ${PATHRT}/fv3_conf/${FV3_RUN:-fv3_run.IN} > fv3_run

atparse < ${PATHTR}/parm/${INPUT_NML:-input.nml.IN} > input.nml

atparse < ${PATHTR}/parm/${MODEL_CONFIGURE:-model_configure.IN} > model_configure

atparse < ${PATHTR}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure

if [[ "Q${INPUT_NEST02_NML:-}" != Q ]] ; then
    atparse < ${PATHTR}/parm/${INPUT_NEST02_NML} > input_nest02.nml
fi

# Allow fv3_run to proceed despite errors when setting up the run
# directory. With "set -e", a failure in 'source ./fv3_run; (e.g.
# if an input file to copy is not found) aborts rt_fv3.sh and the
# calling run_test.sh, and consequently rt.sh reports that the
# test ran successfully. Allowing fv3_run to proceed will in such
# cases lead to errors when the model is run (crash or different
# results), because a necessary input file would be missing.
set +e
source ./fv3_run
set -e

if [[ $SCHEDULER = 'moab' ]]; then
  atparse < $PATHRT/fv3_conf/fv3_msub.IN > job_card
elif [[ $SCHEDULER = 'pbs' ]]; then
  NODES=$(( TASKS / TPN ))
  if (( NODES * TPN < TASKS )); then
    NODES=$(( NODES + 1 ))
  fi
  atparse < $PATHRT/fv3_conf/fv3_qsub.IN > job_card
elif [[ $SCHEDULER = 'sbatch' ]]; then
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
  atparse < $PATHRT/fv3_conf/fv3_bsub.IN > job_card
fi

atparse < ${PATHTR}/parm/${NEMS_CONFIGURE:-nems.configure} > nems.configure

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
