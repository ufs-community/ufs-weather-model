cat <<EOF >${RUNDIR_ROOT}/unit_test${RT_SUFFIX}.env
export UNIT_TEST=${UNIT_TEST}
export CI_TEST=${CI_TEST}
export RT_COMPILER=${RT_COMPILER}
export FHMAX=${FHMAX}
export DAYS=${DAYS}
export RESTART_INTERVAL="${RESTART_INTERVAL:-}"
export RESTART_N=${RESTART_N:-}
export INPES=${INPES}
export JNPES=${JNPES}
export WRITE_GROUP=${WRITE_GROUP}
export WRTTASK_PER_GROUP=${WRTTASK_PER_GROUP}
export NPROC_ICE=${NPROC_ICE:-}
export med_petlist_bounds="${med_petlist_bounds:-}"
export atm_petlist_bounds="${atm_petlist_bounds:-}"
export ocn_petlist_bounds="${ocn_petlist_bounds:-}"
export ice_petlist_bounds="${ice_petlist_bounds:-}"
export THRD=${THRD}
export TASKS=${TASKS}
export TPN=${TPN}
export NODES=${NODES}
export OMP_ENV="${OMP_ENV:-}"
export MPI_PROC_BIND="${MPI_PROC_BIND:-}"
export NFHOUT=${NFHOUT}
export NFHMAX_HF=${NFHMAX_HF}
export NFHOUT_HF=${NFHOUT_HF}
export LIST_FILES="${LIST_FILES}"
EOF
