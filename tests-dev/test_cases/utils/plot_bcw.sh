#!/bin/bash

###############################################################
# User input:
###############################################################

exp_name=C192-BCW
level=850
# set forecast hour at which to generate plots
fcst_time=150

###############################################################
# initialize module
. "${LMOD_ROOT}"/lmod/init/bash

# update path with current directory
export PATH=.:${PATH}

# if grads-to-control app is not present, get it from web:
[[ -f g2ctl ]] || wget -q https://ftp.cpc.ncep.noaa.gov/wd51we/g2ctl/g2ctl
chmod 755 g2ctl

# get colorbar scripts from github
rm -rf gscript
git clone https://github.com/kodamail/gscript.git
cp gscript/color.gs .
cp gscript/xcbar.gs .
rm -rf gscript

# load modules grads and wgrib2
HOSTNAME=$(hostname)
if [[ ${HOSTNAME} == gaea6[1-9] ]]; then module load Core/24.11 ; fi
module load grads wgrib2

# check if model output file exists:
nfiles=$(find GFSPRS.GrbF* | wc -l || true)
if ls GFSPRS.GrbF* >/dev/null 2>&1
then
  echo Using files: GFSPRS.GrbF\*
else
  echo "No model output (GFSPRS.GrbF*) ... exiting"
  exit
fi
# check if plotting fcst time is <= existing forecast time
last_fcst=$(((nfiles-1)*6))
 if (( fcst_time > last_fcst )); then
  echo "Plot time ""${fcst_time}" "is larger than existing fcst time ""${last_fcst}"
  echo "Exiting... "
  exit
fi

# Create grads control and index files
g2ctl -0 GFSPRS.GrbF%f2 GFSPRS.idx > GFSPRS.ctl
sed -i 's/tdef 1/tdef '"${nfiles}"'/g' GFSPRS.ctl
sed -i 's/ 1mo/ 6hr/g' GFSPRS.ctl
gribmap -i GFSPRS.ctl

###############################################################
# Plot baroclinic case
###############################################################
cat << EOF > bcw.gs
exp="${exp_name}"
var="hcurl(ugrdprs,vgrdprs)"
lev=${level}
tt=${fcst_time}
t=tt/6+1
'reinit'
'open GFSPRS.ctl'
'set t 't
'set xlopts 1 6 0.16'
'set ylopts 1 6 0.16'
'set gxout shaded'
'set display color white'
'c'
'set grads off'
'set lat 0 90'
'set lon 90 270'
'color.gs -5 5 1 -kind blue->white->brown'
'set lev ' lev
'd hcurl(ugrdprs,vgrdprs)*1.0e5'
'set gxout contour'
'set cint 3'
'd vgrdprs'
'xcbar.gs 1 10 0.5 0.8'
'draw title 'exp' vort(1.e5) T='tt'hr'
'printim 'exp'_'tt'.png  x1200 y1000 white '
'c'
'quit'
EOF
grads -blc "run bcw.gs" > /dev/null 2>&1
###############################################################
