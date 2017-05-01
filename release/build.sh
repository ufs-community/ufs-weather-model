#!/bin/ksh
set -x

homedir=${1:-`pwd`/trunk}
machine_name=${2:-wcoss_cray}
cd $homedir/tests

compile_option='HYDRO=N 32BIT=Y'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 1
cp $homedir/tests/fv3_1.exe ../NEMS/exe/fv3_gfs_nh.prod.32bit.x

compile_option='HYDRO=N 32BIT=N'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 2
cp $homedir/tests/fv3_2.exe ../NEMS/exe/fv3_gfs_nh.prod.64bit.x

compile_option='HYDRO=Y 32BIT=Y'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 3
cp $homedir/tests/fv3_3.exe ../NEMS/exe/fv3_gfs_hydro.prod.32bit.x

compile_option='HYDRO=Y 32BIT=N'
./compile.sh $homedir/FV3 $machine_name "$compile_option" 4
cp $homedir/tests/fv3_4.exe ../NEMS/exe/fv3_gfs_hydro.prod.64bit.x


#NEMS.x in nemsfv3_baseline/NEMS/exe, and it is also copied to nemsfv3_baseline/fv3_1.exe for regression test.
#/gpfs/hps/emc/global/noscrub/Jun.Wang/nems/20161110/nems4fv3/FV3/jobs
#RUN_fv3_gfs.sh  submit_fv3.sh
