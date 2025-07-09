#!/bin/bash
###############################################################
step=3
###############################################################
# define borders depending on DOMAIN
  lat1=5
  lat2=30
  lon1=160
  lon2=200
###############################################################
# initialize module
. "${LMOD_ROOT}"/lmod/init/bash

# update path with current directory
export PATH=.:${PATH}

# if grads to control app is not present, get it from web
[[ -f g2ctl ]] || wget -q https://ftp.cpc.ncep.noaa.gov/wd51we/g2ctl/g2ctl
chmod 755 g2ctl

# load modules grads and wgrib2
module load grads wgrib2

# check existance of model output file:
if [[ ! -f HURPRS.GrbF00 ]] ; then echo "No model output (HURPRS.GrbF00)" ; exit ; fi
nfiles=$(find . -maxdepth 1 -type f -name 'HURPRS.GrbF*' ! -name '*idx*' ! -name '*ctl*' | wc -l) || true
echo "=== Using model file: HURPRS.GrbF\*\*"
echo "=== Number of files: ${nfiles} Step: ${step} hours"

# Create control and index files
echo "=== G2CTL:"
g2ctl HURPRS.GrbF00 HURPRS.idx > HURPRS.ctl

sed -i "s/options pascals/options pascals template/" HURPRS.ctl
sed -i "s/dset ^HURPRS.GrbF00/dset ^HURPRS.GrbF%f2/" HURPRS.ctl
sed -i "s/1mn/${step}hr/" HURPRS.ctl
sed -i "s/tdef 2/tdef ${nfiles}/" HURPRS.ctl

echo "=== GRIBMAP:"
gribmap -i HURPRS.ctl > /dev/null 2>&1

###############################################################
# Plot 10m wind
###############################################################
echo "=== Plot 10m wind "
ntime=0
while [[ "${ntime}" -lt "${nfiles}" ]]; do
  nhour=$((ntime * step))
  sett=$((ntime + 1))
  printf -v nhour "%03d" "${nhour}"
  echo "=== Plotting hour ${nhour} file: w10-${nhour}.png"

if [[ -f w10-"${nhour}".png ]] ; then rm w10-"${nhour}".png ; fi
cat << EOF > plot.gs
'reinit'
'open HURPRS.ctl'
'set gxout shaded'
'set display color white'
'c'
'set t ${sett}'
'set grads off'
'set lat ${lat1} ${lat2}'
'set lon ${lon1} ${lon2}'
'set clevs 1 2 3 4 6 10 20 30 40'
'd SQRT(UGRD10m*UGRD10m+VGRD10m*VGRD10m)'
'run cbar.gs'
'draw title Wind 10m ${nhour} HR'
'printim w10-${nhour}.png x1200 y1000 '
'c'
'quit'
EOF
grads -blc "run plot.gs" > /dev/null 2>&1
rm -f plot.gs

  ntime=$((ntime + 1))
done
echo "=== Done!"
###############################################################
