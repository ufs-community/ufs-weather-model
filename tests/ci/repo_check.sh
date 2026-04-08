#!/bin/bash
set -eu

get_shas () {
    cwd=$(pwd)
    # Get sha-1's of the top of develop and feature branches
    app="Accept: application/vnd.github.v3+json"
    url=$1
    gitapi=$2
    branch=$3
    base_sha=$(curl -sS -H "$app" $gitapi | jq -r '.commit.sha')
    workspace=$4
    cd $workspace
    git remote add upstream $url
    git fetch -q upstream $branch
    common=$(git merge-base $base_sha @)
    echo $common $base_sha $workspace
    if [[ "$common" != "$base_sha" ]]; then
        printf "%s\n\n" "** $workspace **NOT** up to date"
        flag_sync=false
    fi
    cd $cwd
}

flag_sync=true

declare -A urls branches paths
# UPP, ccpp-framework, and gocart are intentionally excluded because they update at a different cadence 
# and periodically bring in changes. 
submodules="base ufsatm mom6 cice ww3 stoch cmeps cdeps cmake ccpp_physics aqm noahmp cubed_sphere lm4"

urls[base]='https://github.com/ufs-community/ufs-weather-model'
branches[base]='develop'
paths[base]=''

urls[ufsatm]='https://github.com/NOAA-EMC/ufsatm'
branches[ufsatm]='develop'
paths[ufsatm]='UFSATM'

urls[mom6]='https://github.com/NOAA-EMC/MOM6'
branches[mom6]='dev/emc'
paths[mom6]='MOM6-interface/MOM6'

urls[cice]='https://github.com/NOAA-EMC/CICE'
branches[cice]='develop'
paths[cice]='CICE-interface/CICE'

urls[ww3]='https://github.com/NOAA-EMC/WW3'
branches[ww3]='dev/ufs-weather-model'
paths[ww3]='WW3'

urls[stoch]='https://github.com/noaa-psl/stochastic_physics'
branches[stoch]='master'
paths[stoch]='stochastic_physics'

urls[cmeps]='https://github.com/NOAA-EMC/CMEPS'
branches[cmeps]='emc/develop'
paths[cmeps]='CMEPS-interface/CMEPS'

urls[cdeps]='https://github.com/NOAA-EMC/CDEPS'
branches[cdeps]='develop'
paths[cdeps]='CDEPS-interface/CDEPS'

urls[cmake]='https://github.com/NOAA-EMC/CMakeModules'
branches[cmake]='develop'
paths[cmake]='CMakeModules'

urls[ccpp_physics]='https://github.com/ufs-community/ccpp-physics'
branches[ccpp_physics]='ufs/dev'
paths[ccpp_physics]='UFSATM/ccpp/physics'

urls[aqm]='https://github.com/NOAA-EMC/AQM'
branches[aqm]='develop'
paths[aqm]='AQM'

urls[noahmp]='https://github.com/NOAA-EMC/noahmp'
branches[noahmp]='develop'
paths[noahmp]='NOAHMP-interface/noahmp'

urls[cubed_sphere]='https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere'
branches[cubed_sphere]='dev/emc'
paths[cubed_sphere]='UFSATM/fv3/atmos_cubed_sphere'

urls[lm4]='https://github.com/NOAA-GFDL/LM4-NUOPC-driver'
branches[lm4]='develop'
paths[lm4]='LM4-driver'


for submodule in $submodules; do
    url=${urls[$submodule]}
    branch=${branches[$submodule]}
    workspace=${GITHUB_WORKSPACE}'/'${paths[$submodule]}
    gitapi=$(echo "$url" | sed 's/github.com/api.github.com\/repos/g')'/branches/'$branch
    get_shas $url $gitapi $branch $workspace

    if [[ "$flag_sync" == "false" ]]; then
       echo "** ${GITHUB_WORKSPACE} **NOT** up to date"
       exit 1
    fi
done


echo "** ${GITHUB_WORKSPACE} up to date **"

exit 0
