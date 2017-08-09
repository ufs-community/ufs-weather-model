#!/bin/sh --login
#BSUB -L /bin/sh
#BSUB -P FV3GFS-T2O
#BSUB -e err_cray                                        
#BSUB -o out_cray
#BSUB -J fv3gfs
#BSUB -q debug
#BSUB -M 256
#BSUB -extsched 'CRAYLINUX[]'
#BSUB -W 00:30
set -eux

#--------------------------------------------
#--------------------------------------------
# Running NEMS FV3GFS on WCOSS CRAY
# Fanglin.Yang@noaa.gov, April 2017
#--------------------------------------------
#--------------------------------------------

export machine=WCOSS_C            #WCOSS_C, THEIA, etc
export PSLOT=fv3gfs               #user-defined experiment name
export CASE=C96                   #resolution, C96, C384 or C768
export CDATE=2016092900           #initial condition dates  2016092900 2016011812 2016081200               

export BASE_DATA=/gpfs/hps/emc/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE   # source directory
export FIX_FV3=$BASE_DATA/fix/fix_fv3                  #model fixed fields
export IC_DIR=$BASE_DATA/ICs                           #forecast initial conditions 

# temporary running directory
export DATA=/gpfs/hps/stmp/$LOGNAME/${CASE}${PSLOT}${CDATE}     
if [ -d $DATA ]; then rm -rf $DATA ; fi

# directory to save output
export ROTDIR=/gpfs/hps/ptmp/$LOGNAME/$PSLOT/$CASE                    

# NEMS FV3 directory, fv3 release directory  and forecast excutable directory
export FV3DIR=`pwd`/../../..
export FV3DIR_RELEASE=`pwd`/..
export FCSTEXECDIR=$FV3DIR/NEMS/exe

export FHMAX=48                                       #maximum forecast hours
export FHOUT=3                                        #forecast output frequency in hours
#---------------------------------------------------------
#---------------------------------------------------------
case $CASE in
  C96)  export DELTIM=1800; export layout_x=4; export layout_y=8;  export NODES=16;
        export master_grid=1deg;   export REMAP_TASKS=48 ;;
  C384) export DELTIM=450 ; export layout_x=4; export layout_y=8;  export NODES=16;
        export master_grid=0p5deg; export REMAP_TASKS=96 ;;
  C768) export DELTIM=225 ; export layout_x=8; export layout_y=16; export NODES=64;
        export master_grid=0p5deg; export REMAP_TASKS=192 ;;
  *)    echo "grid $CASE not supported, exit"
        exit ;;
esac

export PARM_FV3DIAG=$FV3DIR_RELEASE/parm/parm_fv3diag
export FORECASTSH=$FV3DIR_RELEASE/scripts/exglobal_fcst_nemsfv3gfs.sh         

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

#--------------------------------------------------------------------------
. $MODULESHOME/init/sh 2>/dev/null
cp $FV3DIR/NEMS/src/conf/module-setup.sh.inc module-setup.sh
cp $FV3DIR/NEMS/src/conf/modules.nems modules.fv3
source ./module-setup.sh
module use $( pwd -P )
module load modules.fv3
export IOBUF_PARAMS=${IOBUF_PARAMS:-'*:size=8M:verbose'}
export MPICH_GNI_COLL_OPT_OFF=${MPICH_GNI_COLL_OPT_OFF:-MPI_Alltoallv}
export MKL_CBWR=AVX2
module list

export FCST_LAUNCHER="aprun -n $tasks -N $task_per_node -d $nth_f $j_opt -cc depth" 

#--NSST optins
export nstf_name="2,0,1,0,5"

#--execute the forecast
$FORECASTSH
if [ $? != 0 ]; then echo "forecast failed, exit"; exit; fi


#-------------------------------------------------------------------------
#--convert 6-tile output to global arrary in netCDF format
ymd=`echo $CDATE |cut -c 1-8`
cyc=`echo $CDATE |cut -c 9-10`
export DATA=$ROTDIR/gfs.$ymd/$cyc
export IPD4=YES
export REMAPSH=$FV3DIR_RELEASE/ush/fv3gfs_remap.sh            #remap 6-tile output to global array in netcdf
export REMAPEXE=$FV3DIR_RELEASE/exec/fregrid_parallel
#export master_grid=0p25deg                              #1deg 0p5deg 0p25deg 0p125deg etc
export REMAP_LAUNCHER="aprun -n 48 -N 12 -j 1 -d 2 -cc depth"

$REMAPSH


exit

