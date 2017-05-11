#!/bin/sh
#
# compare test results wih baseline
# 20170505  Jun Wang
#

machine=${1:-wcoss_cray}
CASE=${2:-C96}
if [ $machine = "wcoss_cray" ]; then
  dir1=/gpfs/hps/ptmp/$LOGNAME/fv3gfs/C96/gfs.20160929/00
  dir2=/gpfs/hps/emc/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/baseline/fv3gfs_nh_32bit/${CASE}/gfs.20160929/00
elif [ $machine = "theia" ]; then
  dir1=/scratch4/NCEPDEV/stmp3/${LOGNAME}/fv3gfs/${CASE}/gfs.20160929/00
  dir2=/scratch4/NCEPDEV/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE/baseline/fv3gfs_nh_32bit/${CASE}/gfs.20160929/00
else
   echo "Platform are not supported!"
   exit 0
fi
logfile=`pwd`/log_$CASE

cd $dir2
filelist=`ls`

echo "compare files in $dir1 and $dir2" >$logfile
for fl in $filelist
do

  if [ -d $fl ]; then
    filelistd=`ls $fl`
    echo "compare files in subdirectory $fl"  >> $logfile
    for fld in $filelistd
    do
      echo " cmp $fl/$fld $dir1/$fl/$fld "    >>$logfile
      cmp $fl/$fld $dir1/$fl/$fld             >>$logfile
    done
  else
    echo " cmp $fl $dir1/$fl " >>$logfile
    cmp $fl $dir1/$fl          >>$logfile
  fi

done

ndiff=0
ndiff=`cat $logfile |grep -i diff | wc -l`
if [ $ndiff -gt 0 ]; then
  echo "WARNING: test results are different from baseline!!! Please see the log file log_$CASE"
  exit 8
else
  echo "Succeed! Test results are identical to baseline!"
  rm $logfile
  exit 0
fi

