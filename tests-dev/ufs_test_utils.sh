#!/bin/bash
set -eux

function set_run_task() {
    source default_vars.sh
    source rt_utils.sh
    source "${PATHRT}"/tests/"${TEST_NAME}"
    #jkim compute_petbounds_and_tasks
    if [[ ${ESMF_THREADING} == true ]]; then
        compute_petbounds_and_tasks_esmf_threading
    else
        compute_petbounds_and_tasks_traditional_threading
    fi
    
    TPN=$(( TPN / THRD ))
    NODES=$(( TASKS / TPN ))
    if (( NODES * TPN < TASKS )); then
	NODES=$(( NODES + 1 ))
    fi

    PPN=$(( TASKS / NODES ))
    if (( TASKS - ( PPN * NODES ) > 0 )); then
       PPN=$((PPN + 1))
    fi

    export WLCLK
     
    python -c "import create_xml; create_xml.write_runtest_env()"
    rocoto_create_run_task
   
}

function link_new_baselines() {
    for dir in "${RTPWD}"/*/; do
	dir=${dir%*/}
	[[ -d "${NEW_BASELINE}/${dir##*/}" ]] && continue
	ln -s "${dir%*/}" "${NEW_BASELINE}/"
    done
}
