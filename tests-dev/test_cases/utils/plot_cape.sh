#!/bin/bash

###############################################################
DATE=2020072400
DOMAIN=CONUS # CONUS or global
ANL_RES=1p00 # Analysis resolution (1p00 - 1deg, 0p50 - 0.5 deg, 0p25 - .25 deg)
###############################################################
if [[ "${ANL_RES}" = 1p00 ]]; then
  res_txt="1x1 deg"
elif [[ "${ANL_RES}" = 0p50 ]]; then
  res_txt="0.5x0.5 deg"
elif [[ "${ANL_RES}" = 0p25 ]]; then
  res_txt="0.25x0.25 deg"
else
  echo "Wrong analysis resolution choice, exiting"
  exit
fi
###############################################################
YY=$(echo "${DATE}" | cut -c1-4)
MM=$(echo "${DATE}" | cut -c5-6)
DD=$(echo "${DATE}" | cut -c7-8)
HH=$(echo "${DATE}" | cut -c9-10)
echo "       YYYY MM DD HH"
echo "Date: " "${YY}" "${MM}" "${DD}" "${HH}"
###############################################################
# define borders depending on DOMAIN
if [[ "${DOMAIN}" = CONUS ]]; then
  lat1=20
  lat2=55
  lon1=-130
  lon2=-65
elif [[ "${DOMAIN}" = global ]]; then
  lat1=-90
  lat2=90
  lon1=0
  lon2=360
else
  echo "Wrong DOMAIN, exiting"
  exit
fi
echo
###############################################################
# initialize module
. "${LMOD_ROOT}"/lmod/init/bash

# update path with current directory
export PATH=.:${PATH}

# if grads-to-control app is not present, get it from web:
[[ -f g2ctl ]] || wget -q https://ftp.cpc.ncep.noaa.gov/wd51we/g2ctl/g2ctl
chmod 755 g2ctl

# if color bar script is not present, get it from web:
[[ -f cbar.gs ]] || wget -q http://cola.gmu.edu/grads/scripts/cbar.gs

# load modules grads and wgrib2
HOSTNAME=$(hostname)
if [[ ${HOSTNAME} == gaea6[1-9] ]]; then module load Core/24.11 ; fi
module load grads wgrib2

# check if model output file exxists:
if [[ ! -f GFSPRS.GrbF24 ]] ; then echo "No model output (GFSPRS.GrbF24)" ; exit ; fi
echo Using model file: GFSPRS.GrbF24

# detect machine we are on
case $(hostname -f) in
  gaea5*)       HSD_path=/gpfs/f5/epic/world-shared/HSD_INPUT_DATA        ;; ## gaea5
  gaea6*)       HSD_path=/gpfs/f6/bil-fire8/world-shared/HSD_INPUT_DATA   ;; ## gaea6
  hfe*)         HSD_path=/scratch1/NCEPDEV/nems/role.epic/HSD_INPUT_DATA  ;; ## hera
  hecflow*)     HSD_path=/scratch1/NCEPDEV/nems/role.epic/HSD_INPUT_DATA  ;; ## hera
  fe*)          HSD_path=/mnt/lfs5/HFIP/hfv3gfs/role.epic/HSD_INPUT_DATA  ;; ## jet
  tfe*)         HSD_path=/mnt/lfs5/HFIP/hfv3gfs/role.epic/HSD_INPUT_DATA  ;; ## tjet
  [Oo]rion*)    HSD_path=/work/noaa/epic/role-epic/contrib/HSD_INPUT_DATA ;; ## orion
  [Hh]ercules*) HSD_path=/work/noaa/epic/role-epic/contrib/HSD_INPUT_DATA ;; ## hercules
  derecho*)     HSD_path=/glade/work/epicufsrt/contrib/HSD_INPUT_DATA     ;; ## derecho
  *)            echo "Unknown or unsupported machine " ; exit             ;; ## Unknown platform
esac

# check if analysis file exists:
#   analysis files retrieved from hpss:
#   htar -xvf \
#   /NCEPPROD/hpssprod/runhistory/rh2020/202007/20200724/com_gfs_prod_gfs.20200724_00.gfs_pgrb2.tar \
#   ./gfs.20200724/00/gfs.t00z.pgrb2.0p25.anl
if [[ ! -f gfs.t00z.pgrb2."${ANL_RES}".anl ]] ; then
 echo "No analysis file, copy file from ${HSD_path} directory:"
 cp "${HSD_path}"/"${DATE}"/gfs.t00z.pgrb2."${ANL_RES}".anl .
 [[ -f gfs.t00z.pgrb2."${ANL_RES}".anl ]] || exit
fi
echo Using analysis file: gfs.t00z.pgrb2."${ANL_RES}".anl

# Create grads control and index files
g2ctl gfs.t00z.pgrb2."${ANL_RES}".anl gfs.t00z.pgrb2."${ANL_RES}".idx > gfs.t00z.pgrb2."${ANL_RES}".ctl
gribmap -i gfs.t00z.pgrb2."${ANL_RES}".ctl > /dev/null 2>&1
g2ctl GFSPRS.GrbF24 GFSPRS.idx > GFSPRS.ctl
gribmap -i GFSPRS.ctl > /dev/null 2>&1

###############################################################
# Plot CAPE analysis
###############################################################
echo "Plot CAPE analysis " "${DOMAIN}"
if [[ -f CAPE-analysis.png ]] ; then rm CAPE-analysis.png ; fi
cat << EOF > plot.gs
'reinit'
'open gfs.t00z.pgrb2.${ANL_RES}.ctl'
'set gxout shaded'
'set display color white'
'c'
'set grads off'
'set lat ${lat1} ${lat2}'
'set lon ${lon1} ${lon2}'
'set clevs 0 500 1000 1500 2000 2500 3000 3500'
'd CAPEsfc'
'run cbar.gs'
'draw title CAPE Analysis ${res_txt} valid ${YY}/${MM}/${DD}/${HH}Z'
'printim CAPE-analysis_${DOMAIN}.png x1200 y1000 '
'c'
'quit'
EOF
grads -blc "run plot.gs" > /dev/null 2>&1
###############################################################
# Plot CAPE model output
###############################################################
echo "Plot CAPE model output " "${DOMAIN}"
if [[ -f CAPE-model.png ]] ; then rm CAPE-model.png ; fi
if [[ -f plot.gs ]] ; then rm plot.gs ; fi
cat << EOF > plot.gs
'reinit'
'open GFSPRS.ctl'
'set gxout shaded'
'set display color white'
'c'
'set grads off'
'set lat ${lat1} ${lat2}'
'set lon ${lon1} ${lon2}'
'set clevs 0 500 1000 1500 2000 2500 3000 3500'
'd CAPEsfc'
'run cbar.gs'
'draw title CAPE Model fcst (+24h) valid ${YY}/${MM}/${DD}/${HH}Z'
'printim CAPE-model_${DOMAIN}.png x1200 y1000 '
'c'
'quit'
EOF
grads -blc "run plot.gs" > /dev/null 2>&1
###############################################################
# Plot MSL pressure analysis
###############################################################
echo "Plot PRSMSL analysis " "${DOMAIN}"
if [[ -f MSL-analysis.png ]] ; then rm MSL-analysis.png ; fi
cat << EOF > plot.gs
'reinit'
'open gfs.t00z.pgrb2.${ANL_RES}.ctl'
'set gxout shaded'
'set display color white'
'c'
'set grads off'
'set lat ${lat1} ${lat2}'
'set lon ${lon1} ${lon2}'
'set clevs 1000 1002 1004 1006 1008 1010 1012 1014 1016 1018 1020 1022'
'd PRMSLmsl/100'
'run cbar.gs'
'draw title PRSMSL [mb] Analysis ${res_txt} valid ${YY}/${MM}/${DD}/${HH}Z'
 'printim MSL-analysis_${DOMAIN}.png x1200 y1000 '
'c'
'quit'
EOF
grads -blc "run plot.gs" > /dev/null 2>&1
###############################################################
# Plot MSL pressure model output
###############################################################
echo "Plot PRSMSL model output " "${DOMAIN}"
if [[ -f MSL-model.png ]] ; then rm MSL-model.png ; fi
cat << EOF > plot.gs
'reinit'
'open GFSPRS.ctl'
'set gxout shaded'
'set display color white'
'c'
'set grads off'
'set lat ${lat1} ${lat2}'
'set lon ${lon1} ${lon2}'
'set clevs 1000 1002 1004 1006 1008 1010 1012 1014 1016 1018 1020 1022'
'd PRMSLmsl/100'
'run cbar.gs'
'draw title PRSMSL [mb] Model fcst (+24) valid ${YY}/${MM}/${DD}/${HH}Z'
'printim MSL-model_${DOMAIN}.png x1200 y1000 '
'c'
'quit'
EOF
grads -blc "run plot.gs" > /dev/null 2>&1
###############################################################
# Plot Vorticity analysis
###############################################################
echo "Plot vorticity analysis " "${DOMAIN}"
if [[ -f VORT-analysis.png ]] ; then rm VORT-analysis.png ; fi
cat << EOF > plot.gs
'reinit'
'open gfs.t00z.pgrb2.${ANL_RES}.ctl'
'set gxout shaded'
'set display color white'
'c'
'set grads off'
'set lev 500'
'set lat ${lat1} ${lat2}'
'set lon ${lon1} ${lon2}'
'set rgb 20 255 255 255'
'set rgb 21 255 250 170'
'set rgb 22 255 232 100'
'set rgb 23 255 192   0'
'set rgb 24 255 160   0'
'set rgb 25 255  96   0'
'set rgb 26 255  50   0'
'set rgb 27 225  20   0'
'set rgb 28 192   0   0'
'set ccols 20 21 22 23 24 25 26 27 28'
 'set clevs 2 4 6 8 10 15 20 25'
'define h=hcurl(ugrdprs,vgrdprs)*1.0e5'
'd sqrt(h*h)'
'run cbar.gs'
'set gxout contour'
'd hgtprs'
'draw title Analysis 500mb Height and Abs Vorticity ${res_txt} \ Valid ${YY}/${MM}/${DD}/${HH}Z'
'printim VORT-analysis_${DOMAIN}.png x1200 y1000 '
'c'
'quit'
EOF
grads -blc "run plot.gs" > /dev/null 2>&1
###############################################################
# Plot Vorticity model output
###############################################################
echo "Plot vorticity model output " "${DOMAIN}"
if [[ -f VORT-model.png ]] ; then rm VORT-model.png ; fi
cat << EOF > plot.gs
'reinit'
'open GFSPRS.ctl'
'set gxout shaded'
'set display color white'
'c'
'set grads off'
'set lev 500'
'set lat ${lat1} ${lat2}'
'set lon ${lon1} ${lon2}'
'set rgb 20 255 255 255'
'set rgb 21 255 250 170'
'set rgb 22 255 232 100'
'set rgb 23 255 192   0'
'set rgb 24 255 160   0'
'set rgb 25 255  96   0'
'set rgb 26 255  50   0'
'set rgb 27 225  20   0'
'set rgb 28 192   0   0'
'set ccols 20 21 22 23 24 25 26 27 28'
 'set clevs 2 4 6 8 10 15 20 25'
'define h=hcurl(ugrdprs,vgrdprs)*1.0e5'
'd sqrt(h*h)'
'run cbar.gs'
'set gxout contour'
'd hgtprs'
'draw title Model fcst (+24) 500mb Height and Abs Vorticity \ Valid ${YY}/${MM}/${DD}/${HH}Z'
'printim VORT-model_${DOMAIN}.png x1200 y1000 '
'c'
'quit'
EOF
grads -blc "run plot.gs" > /dev/null 2>&1
###############################################################
rm -f plot.gs
