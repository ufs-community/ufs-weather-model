#!/bin/ksh
set -x

machine_name=${1:-wcoss_cray}
homedir=${2:-`pwd`/../../..}
cd $homedir/tests

compile_option='HYDRO=N 32BIT=Y'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 1
cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_nh.prod.32bit.x
rm $homedir/tests/fv3_1.exe

compile_option='HYDRO=N 32BIT=N'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 1
cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_nh.prod.64bit.x
rm $homedir/tests/fv3_1.exe

compile_option='HYDRO=Y 32BIT=Y'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 1
cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_hydro.prod.32bit.x
rm $homedir/tests/fv3_1.exe

compile_option='HYDRO=Y 32BIT=N'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 1
cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_hydro.prod.64bit.x
rm $homedir/tests/fv3_1.exe
