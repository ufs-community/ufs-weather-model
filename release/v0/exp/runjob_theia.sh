#!/bin/sh -l
#PBS -o out
#PBS -e err
#PBS -N fv3gfs
#PBS -A nems
#PBS -q dev
#PBS -d .
#PBS -l nodes=16:ppn=24
#PBS -l walltime=00:10:00

set -x

#--------------------------------------------
#--------------------------------------------
# Running NEMS FV3GFS on WCOSS CRAY
# Fanglin.Yang@noaa.gov, April 2017
#--------------------------------------------
#--------------------------------------------

export machine=theia              ;#WCOSS_C, theia, etc
export PSLOT=fv3gfs               ;#user-defined experiment name
export CASE=C384                  ;#resolution, C96, C384 or C768
export CDATE=2016092900           ;#initial condition dates  2016092900 2016011812 2016081200               

export BASE_GSM=/gpfs/hps/emc/global/noscrub/$LOGNAME/svn/FV3GFS_V0_RELEASE   ;# source directory
export FIX_FV3=$BASE_GSM/fix/fix_fv3                  ;#model fixed fields
export IC_DIR=$BASE_GSM/ICs                           ;#forecast initial conditions 

# temporary running directory
export DATA=/gpfs/hps/stmp/$LOGNAME/${CASE}${PSLOT}${CDATE}     
if [ -d $DATA ]; then rm -rf $DATA ; fi

# directory to save output
export ROTDIR=/gpfs/hps/ptmp/$LOGNAME/$PSLOT/$CASE                    

# NEMS FV3GFS forecast executable directory
FV3DIR=${1:-/gpfs/hps/emc/global/noscrub/Fanglin.Yang/svn/fv3gfs/NEMSfv3gfs/trunk}
export FCSTEXECDIR=$FV3DIR/NEMS/exe

export FHMAX=240                                      ;#maximum forecast hours
export FHOUT=3                                        ;#forecast output frequency in hours
#---------------------------------------------------------
#---------------------------------------------------------
case $CASE in
  C96)  export DELTIM=1800; export layout_x=3; export layout_y=8;  export NODES=24  ;;
  C384) export DELTIM=450 ; export layout_x=3; export layout_y=8;  export NODES=96 ;;
  C768) export DELTIM=225 ; export layout_x=6; export layout_y=16; export NODES=64;;
  *)    echo "grid $CASE not supported, exit"
        exit ;;
esac

export PARM_FV3DIAG=$BASE_GSM/parm/parm_fv3diag
export FORECASTSH=$FV3DIR/release/v0/exglobal_fcst_nemsfv3gfs.sh         

#---determine task configuration
export nth_f=2                             # number of threads 
export npe_node_f=24                       # number of pes per node 
export task_per_node=$((npe_node_f/nth_f))
export tasks=$((NODES*task_per_node))      # number of tasks 


export MODE=32bit           # dycore precision:   32bit, 64bit
export TYPE=nh              # hydrostatic option: nh, hydro
export HYPT=off             # hyperthreading:     on, off  
export COMP="prod"          # compiling mode:     debug, repro, prod
if [ ${HYPT} = on ]; then
   export hyperthread=".true."
   export j_opt="-j 2"
else
   export hyperthread=".false."
   export j_opt="-j 1"
fi

export FCSTEXEC=fv3_gfs_${TYPE}.${COMP}.${MODE}.x
export APRUN="aprun -n $tasks -N $task_per_node -d $nth_f $j_opt -cc depth" 
#--------------------------------------------------------------------------

#--execute the forecast
$FV3dir/release/exglobal_fcst_nemsfv3gfs.sh
if [ $? != 0 ]; then echo "forecast failed, exit"; exit; fi


#-------------------------------------------------------------------------
#--convert 6-tile output to global arrary in netCDF format
ymd=`echo $CDATE |cut -c 1-8`
cyc=`echo $CDATE |cut -c 9-10`
export DATA=$ROTDIR/gfs.$ymd/$cyc
export IPD4=YES
export REMAPSH=$FV3DIR/release/v0/fv3gfs_remap.sh            #remap 6-tile output to global array in netcdf
export REMAPEXE=$FV3DIR/release/exec/fregrid_parallel
export master_grid=0p25deg                              #1deg 0p5deg 0p25deg 0p125deg etc
if [ $nth_f -eq 1 ] ; then
  export FCST_LAUNCHER="env LD_LIBRARY_PATH=$LD_LIBRARY_PATH $mpiexec -n $PE1"
else
  export FCST_LAUNCHER="env LD_LIBRARY_PATH=$LD_LIBRARY_PATH $mpiexec -np $PE1"
fi

#export APRUN_REMAP="aprun -n 48 -N 12 -j 1 -d 2 -cc depth"

$FV3_DIR/release/v0/util/fv3gfs_remap.sh


exit

