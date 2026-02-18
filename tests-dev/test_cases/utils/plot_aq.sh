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

#--- path to model output (usually ./) if you copy this script to run directory
#--- staged model output data on HPSS: /5year/NCEPDEV/emc-meso/Ratko.Vasic/AQUAPLANET/1yr-results.tar
data_path=/scratch3/NAGAPE/epic/role.epic/Aquaplanet
#data_path=./


#========================================================
# initialize module
. "${LMOD_ROOT}"/lmod/init/bash

# load grads module
HOSTNAME=$(hostname)
if [[ ${HOSTNAME} == gaea6[1-9] ]]; then module load Core/24.11 ; fi
module load grads

# if color bar script is not present, get it from github:
[[ -f cbar.gs ]] || wget -q https://raw.githubusercontent.com/NOAA-EPIC/Aquaplanet/refs/heads/main/utils/cbar.gs

#========================================================
# Plot Jet stream, four seasons
#========================================================
for season in Winter Spring Summer Fall
do

echo "Jet " "${season}"

case ${season} in
  Winter) hour=${winter_start} ;;
  Spring) hour=${spring_start} ;;
  Summer) hour=${summer_start} ;;
  Fall)   hour=${fall_start}   ;;
  *)      echo "There must be at least one season" ; exit   ;;
esac

cat << EOF > plot.j
reinit
set gxout shaded
set display color white
c
set grads off
EOF

i=1
 while (( i <= season_len ))
do

echo sdfopen "${data_path}"/atmf"${hour}".nc      >> plot.j
echo set z 73                                   >> plot.j
if (( i == 1 )); then
echo define utot=ugrd.1\(t=1\)                  >> plot.j
else
echo define utot=utot+ugrd."${i}"\(t=1\)        >> plot.j
fi

((i++))
hour=$((hour + 24))
done

cat << EOF >> plot.j
define ut=utot/"${season_len}"
set clevs -5 0 5 10 15 20 30 40
d ut
run cbar.gs
draw title ${season} Jet stream
printim ${season}-jet.png x1200 y1000
c
EOF

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

echo "Temp " "${season}"

case ${season} in
  Winter) hour=${winter_start} ;;
  Spring) hour=${spring_start} ;;
  Summer) hour=${summer_start} ;;
  Fall)   hour=${fall_start}   ;;
  *)      echo "There must be at least one season" ; exit   ;;
esac

cat << EOF > plot.j
reinit
set gxout shaded
set display color white
c
set grads off
EOF

i=1
while (( i <= season_len ))
do

echo sdfopen "${data_path}"/atmf"${hour}".nc      >> plot.j
echo set z 49                                   >> plot.j
if (( i == 1 )); then
echo define ttot=tmp.1\(t=1\)                   >> plot.j
else
echo define ttot=ttot+tmp."${i}"\(t=1\)         >> plot.j
fi

((i++))
hour=$((hour + 24))
done

cat << EOF >> plot.j
define tt=ttot/"${season_len}"
set clevs -24 -20 -16 -12 -8 -4 0
d tt-273.15
run cbar.gs
draw title ${season} Temp 500hPa
printim ${season}-t500.png x1200 y1000
c
EOF

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

echo "Prec " "${season}"

case ${season} in
  Winter) hour=${winter_start} ;;
  Spring) hour=${spring_start} ;;
  Summer) hour=${summer_start} ;;
  Fall)   hour=${fall_start}   ;;
  *)      echo "There must be at least one season" ; exit   ;;
esac

cat << EOF > plot.j
reinit
set gxout shaded
set display color white
c
set grads off
set rgb 40 128 0 160
set rgb 42 128 0 208
set rgb 44 128 0 255
set rgb 46 96 0 224
set rgb 48 0 0 192
set rgb 50 0 88 208
set rgb 52 0 144 224
set rgb 54 0 200 240
set rgb 56 0 255 255
set rgb 58 128 255 64
set rgb 60 192 255 0
EOF

i=1
while (( i <= season_len ))
do

echo sdfopen "${data_path}"/sfcf"${hour}".nc      >> plot.j
if (( i == 1 )); then
echo define ptot=prate_ave.1\(t=1\)             >> plot.j
else
echo define ptot=ptot+prate_ave."${i}"\(t=1\)   >> plot.j
fi

((i++))
hour=$((hour + 24))
done

cat << EOF >> plot.j
define pt=ptot/"${season_len}"
set clevs 0 1 2 3 4 5 6 7 8 9
set ccols 60 58 56 54 52 50 48 46 44 42 40
d pt*86400
run cbar.gs
draw title ${season} Precip mm/day
printim ${season}-prec.png x1200 y1000
c
EOF

echo \'exec plot.j\'       > plot.gs
echo \'quit\'             >> plot.gs

grads -blc "run plot.gs" > /dev/null 2>&1

done
#========================================================
#--- clean:
rm -f plot.j plot.gs
#========================================================
