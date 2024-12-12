#!/usr/bin/env bash

# Output a CSV report of disk usage on subdirs of some path
# Usage: 
#    [JOB_NAME=<ci_job>] [BUILD_NUMBER=<n>] [UFS_COMPILER=<intel>] [UFS_PLATFORM=<machine_id>] disk_usage path depth size outfile.csv
#
# args:
#    directory=$1
#    depth=$2
#    size=$3
#    outfile=$4

export UFS_PLATFORM=${UFS_PLATFORM:-${NODE_NAME,,}}
export UFS_COMPILER=${UFS_COMPILER:-intel}
[[ -n ${WORKSPACE} ]] || WORKSPACE="$(pwd)"
[[ -n ${UFS_PLATFORM} ]] || UFS_PLATFORM="$(hostname -s 2>/dev/null)" || UFS_PLATFORM="$(hostname 2>/dev/null)"
[[ -n ${UFS_COMPILER} ]] || UFS_COMPILER="compiler"
echo "STAGE_NAME=${STAGE_NAME%% *}" # from pipeline

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" > /dev/null 2>&1 && pwd)"
echo "script_dir=${script_dir}"

# Get repository root from Jenkins WORKSPACE variable if set, otherwise, set
# relative to script directory.
declare workspace
if [[ -d "${WORKSPACE}/${UFS_PLATFORM}" ]]; then
    workspace="${WORKSPACE}/${UFS_PLATFORM}"
    outfile="${4:-${workspace}-${UFS_PLATFORM}-${UFS_COMPILER}-disk-usage${STAGE_NAME%% *}.csv}"
else
    workspace="$(cd -- "${script_dir}/../.." && pwd)"
    outfile="${4:-${workspace}/${UFS_PLATFORM}-${UFS_COMPILER}-disk-usage${STAGE_NAME%% *}.csv}"
fi
echo "workspace=${workspace}"
echo "outfile=${outfile}"

function disk_usage() {
    local directory="${1:-${PWD}}"
    local depth="${2:-1}"
    local size="${3:-k}"
    echo "Disk usage: ${JOB_NAME:-ci}/${UFS_PLATFORM}/$(basename ${directory})"
    (
    cd ${directory} || exit 1
    echo "Platform,Build,Owner,Group,Inodes,${size:-k}bytes,Access Time,Filename"
    du -Px -d ${depth:-1} --inode --exclude='./workspace' | \
        while read -r line ; do
            read -ra arr<<<"${line}"; inode="${arr[0]}"; filename="${arr[1]}";
            echo "${UFS_PLATFORM}-${UFS_COMPILER:-compiler},${JOB_NAME:-ci}/${BUILD_NUMBER:-0},$(stat -c '%U,%G' "${filename:-.}" || true),${inode:-0},$(du -Px -s -${size:-k} --time "${filename:-null}" 2>/dev/null || true)" | tr '\t' ',' || true;
        done | sort -t, -k5 -n #-r
    )
    echo ""
}

disk_usage "${1}" "${2}" "${3}" | tee ${outfile}
