#!/bin/bash
#set -eu

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
    if [[ $common != $base_sha ]]; then
        printf "%s\n\n" "** $workspace **NOT** up to date"
        flag_sync=false
    fi
    cd $cwd
}

flag_sync=true

declare -A urls branches pathes
submodules="base fv3 mom6 cice ww3 stoch cmeps cdeps hycom ccpp_physics aqm noahmp cubed_sphere"

urls[base]='https://github.com/ufs-community/ufs-weather-model'
branches[base]='develop'
pathes[base]=''

urls[fv3]='https://github.com/NOAA-EMC/fv3atm'
branches[fv3]='develop'
pathes[fv3]='FV3'

urls[mom6]='https://github.com/NOAA-EMC/MOM6'
branches[mom6]='dev/emc'
pathes[mom6]='MOM6-interface/MOM6'

urls[cice]='https://github.com/NOAA-EMC/CICE'
branches[cice]='emc/develop'
pathes[cice]='CICE-interface/CICE'

urls[ww3]='https://github.com/NOAA-EMC/WW3'
branches[ww3]='dev/ufs-weather-model'
pathes[ww3]='WW3'

urls[stoch]='https://github.com/noaa-psl/stochastic_physics'
branches[stoch]='master'
pathes[stoch]='stochastic_physics'

urls[gocart]='https://github.com/GEOS-ESM/GOCART'
branches[gocart]='develop'
pathes[gocart]='GOCART'

urls[cmeps]='https://github.com/NOAA-EMC/CMEPS'
branches[cmeps]='emc/develop'
pathes[cmeps]='CMEPS-interface/CMEPS'

urls[cdeps]='https://github.com/NOAA-EMC/CDEPS'
branches[cdeps]='develop'
pathes[cdeps]='CDEPS-interface/CDEPS'

urls[hycom]='https://github.com/NOAA-EMC/HYCOM-src'
branches[hycom]='emc/develop'
pathes[hycom]='HYCOM-interface/HYCOM'

urls[cmake]='https://github.com/NOAA-EMC/CMakeModules'
branches[cmake]='develop'
pathes[cmake]='CMakeModules'

urls[ccpp_physics]='https://github.com/ufs-community/ccpp-physics'
branches[ccpp_physics]='ufs/dev'
pathes[ccpp_physics]='FV3/ccpp/physics'

urls[ccpp_framework]='https://github.com/NCAR/ccpp-framework'
branches[ccpp_framework]='main'
pathes[ccpp_framework]='FV3/ccpp/framework'

urls[aqm]='https://github.com/NOAA-EMC/AQM'
branches[aqm]='develop'
pathes[aqm]='AQM'

urls[noahmp]='https://github.com/NOAA-EMC/noahmp'
branches[noahmp]='develop'
pathes[noahmp]='NOAHMP-interface/noahmp'

#urls[upp]='https://github.com/NOAA-EMC/UPP'
#branches[upp]='develop'
#pathes[upp]='upp'

urls[cubed_sphere]='https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere'
branches[cubed_sphere]='dev/emc'
pathes[cubed_sphere]='FV3/atmos_cubed_sphere'

for submodule in $submodules; do
    url=${urls[$submodule]}
    branch=${branches[$submodule]}
    workspace=${GITHUB_WORKSPACE}'/'${pathes[$submodule]}
    gitapi=$(echo "$url" | sed 's/github.com/api.github.com\/repos/g')'/branches/'$branch
    get_shas $url $gitapi $branch $workspace
done

if [[ ! $flag_sync ]]; then
    echo "** ${GITHUB_WORKSPACE} **NOT** up to date"
    exit 1
fi

echo "** ${GITHUB_WORKSPACE} up to date **"

exit 0
