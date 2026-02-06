#!/bin/bash

#========================================================
# User input:
#========================================================
#--- start hour of the season:
winter_start=2184
spring_start=3984
summer_start=5784
fall_start=7584

#--- duration of season in days (usually 90):
season_len=75

#--- path to output (usually ./) if you copy this script to run directory
#--- staged data on HPSS: /5year/NCEPDEV/emc-meso/Ratko.Vasic/AQUAPLANET/1yr-results.tar
out_pth=/scratch3/NAGAPE/epic/Ratko.Vasic/GFS_v17/control_c48_intel/tmp-002
#out_pth=./


#========================================================
# initialize module
. "${LMOD_ROOT}"/lmod/init/bash

# load grads module
HOSTNAME=$(hostname)
if [[ ${HOSTNAME} == gaea6[1-9] ]]; then module load Core/24.11 ; fi
module load grads

# if color bar script is not present, get it from github:
[[ -f cbar.gs ]] || wget -q //raw.githubusercontent.com/RatkoVasic-NOAA/Aquaplanet/refs/heads/main/utils/cbar.gs

#========================================================
# Plot Jet stream, four seasons
#========================================================
for season in Winter Spring Summer Fall
do

echo "Jet " $season

case ${season} in
  Winter) hour=$winter_start ;;
  Spring) hour=$spring_start ;;
  Summer) hour=$summer_start ;;
  Fall)   hour=$fall_start   ;;
esac

echo reinit                                  > plot.j
echo set gxout shaded                       >> plot.j
echo set display color white                >> plot.j
echo c                                      >> plot.j
echo set grads off                          >> plot.j

i=1
while (( i <= $season_len ))
do

echo sdfopen $out_pth/atmf$hour.nc          >> plot.j
echo set z 73                               >> plot.j
if (( i == 1 )); then
echo define utot=ugrd.1\(t=1\)              >> plot.j
else
echo define utot=utot+ugrd.$i\(t=1\)        >> plot.j
fi

((i++))
hour=$((hour + 24))
done

echo define ut=utot/$season_len             >> plot.j
echo set clevs -5 0 5 10 15 20 30 40        >> plot.j
echo d ut                                   >> plot.j
echo run cbar.gs                            >> plot.j
echo draw title $season Jet stream          >> plot.j
echo printim $season-jet.png x1200 y1000    >> plot.j
echo c                                      >> plot.j

echo \'exec plot.j\'       > plot.gs
echo \'quit\'             >> plot.gs

grads -blc "run plot.gs" > /dev/null 2>&1

done
#========================================================

#========================================================
# Plot Temp 500hPa, four seasons
#========================================================
for season in Winter Spring Summer Fall
do

echo "Temp " $season

case ${season} in
  Winter) hour=$winter_start ;;
  Spring) hour=$spring_start ;;
  Summer) hour=$summer_start ;;
  Fall)   hour=$fall_start   ;;
esac

echo reinit                                  > plot.j
echo set gxout shaded                       >> plot.j
echo set display color white                >> plot.j
echo c                                      >> plot.j
echo set grads off                          >> plot.j

i=1
while (( i <= $season_len ))
do

echo sdfopen $out_pth/atmf$hour.nc          >> plot.j
echo set z 49                               >> plot.j
if (( i == 1 )); then
echo define ttot=tmp.1\(t=1\)               >> plot.j
else
echo define ttot=ttot+tmp.$i\(t=1\)         >> plot.j
fi

((i++))
hour=$((hour + 24))
done

echo define tt=ttot/$season_len             >> plot.j
echo set clevs -24 -20 -16 -12 -8 -4 0      >> plot.j
echo d tt-273.15                            >> plot.j
echo run cbar.gs                            >> plot.j
echo draw title $season Temp 500hPa         >> plot.j
echo printim $season-t500.png x1200 y1000   >> plot.j
echo c                                      >> plot.j

echo \'exec plot.j\'       > plot.gs
echo \'quit\'             >> plot.gs

grads -blc "run plot.gs" > /dev/null 2>&1

done
#========================================================

#========================================================
# Plot Precipitation, four seasons
#========================================================
for season in Winter Spring Summer Fall
do

echo "Prec " $season

case ${season} in
  Winter) hour=$winter_start ;;
  Spring) hour=$spring_start ;;
  Summer) hour=$summer_start ;;
  Fall)   hour=$fall_start   ;;
esac

echo reinit                                  > plot.j
echo set gxout shaded                       >> plot.j
echo set display color white                >> plot.j
echo c                                      >> plot.j
echo set grads off                          >> plot.j
echo set rgb 40 128 0 160                   >> plot.j
echo set rgb 42 128 0 208                   >> plot.j
echo set rgb 44 128 0 255                   >> plot.j
echo set rgb 46 96 0 224                    >> plot.j
echo set rgb 48 0 0 192                     >> plot.j
echo set rgb 50 0 88 208                    >> plot.j
echo set rgb 52 0 144 224                   >> plot.j
echo set rgb 54 0 200 240                   >> plot.j
echo set rgb 56 0 255 255                   >> plot.j
echo set rgb 58 128 255 64                  >> plot.j
echo set rgb 60 192 255 0                   >> plot.j

i=1
while (( i <= $season_len ))
do

echo sdfopen $out_pth/sfcf$hour.nc          >> plot.j
if (( i == 1 )); then
echo define ptot=prate_ave.1\(t=1\)         >> plot.j
else
echo define ptot=ptot+prate_ave.$i\(t=1\)   >> plot.j
fi

((i++))
hour=$((hour + 24))
done

echo define pt=ptot/$season_len             >> plot.j
echo set clevs 0 1 2 3 4 5 6 7 8 9          >> plot.j
echo set ccols 60 58 56 54 52 50 48 46 44 42 40 >> plot.j
echo d pt*86400                             >> plot.j
echo run cbar.gs                            >> plot.j
echo draw title $season Precip mm/day       >> plot.j
echo printim $season-prec.png x1200 y1000   >> plot.j
echo c                                      >> plot.j

echo \'exec plot.j\'       > plot.gs
echo \'quit\'             >> plot.gs

grads -blc "run plot.gs" > /dev/null 2>&1

done
#========================================================
#--- clean:
rm -f plot.j plot.gs
#========================================================
