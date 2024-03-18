#!/bin/bash
set -eu

# This script checks if head repo of PR is up to date with ufs-weather-model develop
# Checks for top level (ufs-weather-model) and next level components (submodules)
result() {
  if [[ -n $comment ]]; then
    logID=$1
    comment="@$logID please bring these up to date with respective authoritative repositories\n"$comment
    printf %s "$comment"
    #exit 1
  fi
}

# Declare variables
declare -A base fv3 mom6 cice ww3 stoch gocart cmeps cdeps hycom cmake ccpp-framework ccpp-physics upp atmos_cubed_sphere
submodules="fv3 mom6 cice ww3 stoch gocart cmeps cdeps hycom cmake ccpp-framework ccpp-physics upp atmos_cubed_sphere"
comment=''
ownerID=$1

# Base branch: this is the top of develop of ufs-weather-model
base[repo]='https://github.com/ufs-community/ufs-weather-model'
base[branch]='develop'

# Submodules to check
fv3[repo]='https://github.com/NOAA-EMC/fv3atm'
fv3[branch]='develop'
fv3[dir]='FV3'

mom6[repo]='https://github.com/NOAA-EMC/MOM6'
mom6[branch]='dev/emc'
mom6[dir]='MOM6-interface/MOM6'

cice[repo]='https://github.com/NOAA-EMC/CICE'
cice[branch]='emc/develop'
cice[dir]='CICE-interface/CICE'

ww3[repo]='https://github.com/NOAA-EMC/WW3'
ww3[branch]='dev/ufs-weather-model'
ww3[dir]='WW3'

stoch[repo]='https://github.com/noaa-psl/stochastic_physics'
stoch[branch]='master'
stoch[dir]='stochastic_physics'

gocart[repo]='https://github.com/GEOS-ESM/GOCART'
gocart[branch]='main'
gocart[dir]='GOCART'

cmeps[repo]='https://github.com/NOAA-EMC/CMEPS'
cmeps[branch]='emc/develop'
cmeps[dir]='CMEPS-interface/CMEPS'

cdeps[repo]='https://github.com/NOAA-EMC/CDEPS'
cdeps[branch]='develop'
cdeps[dir]='CDEPS-interface/CDEPS'

hycom[repo]='https://github.com/NOAA-EMC/HYCOM-src'
hycom[branch]='emc/develop'
hycom[dir]='HYCOM-interface/HYCOM'

cmake[repo]='https://github.com/NOAA-EMC/CMakeModules'
cmake[branch]='develop'
cmake[dir]='CMakeModules'

ccpp-framework[repo]='https://github.com/NCAR/ccpp-framework'
ccpp-framework[branch]='main'
ccpp-framework[dir]='ccpp/framework'

ccpp-physics[repo]='https://github.com/ufs-community/ccpp-physics'
ccpp-physics[branch]='ufs/dev'
ccpp-physics[dir]='ccpp/physics'

upp[repo]='https://github.com/NOAA-EMC/UPP'
upp[branch]='develop'
upp[dir]='upp'

atmos_cubed_sphere[repo]='https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere'
atmos_cubed_sphere[branch]='main'
atmos_cubed_sphere[dir]='atmos_cubed_sphere'

# Get sha-1's of the top of develop of ufs-weather-model
app="Accept: application/vnd.github.v3+json"
url="https://api.github.com/repos/ufs-community/ufs-weather-model/branches/develop"
base[sha]=$(curl -sS -H "$app" $url | jq -r '.commit.sha')
for submodule in $submodules; do
  eval url=https://api.github.com/repos/ufs-community/ufs-weather-model/contents/'${'$submodule'[dir]}'
  eval $submodule'[sha]=$(curl -sS -H "$app" $url | jq -r '.sha')'
done

# Check if the head branch is up to date with the base branch
cd ${GITHUB_WORKSPACE}
git remote add upstream ${base[repo]}
git fetch -q upstream ${base[branch]}
common=$(git merge-base ${base[sha]} @)
if [[ $common != ${base[sha]} ]]; then
  comment="* ufs-weather-model **NOT** up to date\n"
fi

for submodule in $submodules; do
  eval cd ${GITHUB_WORKSPACE}/'${'$submodule'[dir]}'
  eval git remote add upstream '${'$submodule'[repo]}'
  eval git fetch -q upstream '${'$submodule'[branch]}'
  common=$(eval git merge-base '${'$submodule'[sha]}' @)
  if (eval test $common != '${'$submodule'[sha]}'); then
    comment+="* $submodule **NOT** up to date\n"
  fi
done

result $ownerID
exit 0
