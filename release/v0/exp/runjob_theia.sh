#!/bin/sh -l
#PBS -o out_theia
#PBS -e err_theia
#PBS -N fv3gfs
#PBS -A nems
###PBS -q dev
#PBS -d .
#PBS -l nodes=24:ppn=12
#PBS -l walltime=04:00:00

set -x

#------------------------------------------------------------------
#------------------------------------------------------------------
# Running NEMS FV3GFS on Theia
#------------------------------------------------------------------
#notes:
# this job card is for C96 case. If you are running
# C384 or C768 cases, please make the following change:
#
#   for C384, change line 8 and line 33 to:
#      #PBS -l nodes=96:ppn=12      
#      export CASE=C384            
#
#   for C768, change line 8 and line 33 to:
#      #PBS -l nodes=192:ppn=12
#      export CASE=C768            
#
#------------------------------------------------------------------

export machine=theia              #WCOSS_C, theia, etc
export PSLOT=fv3gfs               #user-defined experiment name
export CASE=C96                   #resolution, C96, C384 or C768
export CDATE=2016092900           #initial condition dates  2016092900 2016011812 2016081200               

export BASE_DATA=/scratch4/NCEPDEV/nems/noscrub/emc.nemspara/FV3GFS_V0_RELEASE ;# data directory
export FIX_FV3=$BASE_DATA/fix/fix_fv3                 #model fixed fields
export IC_DIR=$BASE_DATA/ICs                          #forecast initial conditions 

# temporary running directory
export DATA=/scratch4/NCEPDEV/stmp3/$LOGNAME/${CASE}${PSLOT}${CDATE}     
if [ -d $DATA ]; then rm -rf $DATA ; fi

# directory to save output
export ROTDIR=/scratch4/NCEPDEV/stmp3/$LOGNAME/$PSLOT/$CASE                    

# NEMS FV3GFS forecast executable directory
FV3DIR=${1:-`pwd`/../../..}
FV3DIR_RELEASE=${1:-`pwd`/..}
export FCSTEXECDIR=$FV3DIR/NEMS/exe

export FHMAX=48                                       #maximum forecast hours
export FHOUT=3                                        #forecast output frequency in hours
#---------------------------------------------------------
#---------------------------------------------------------
case $CASE in
  C96)  export DELTIM=1800; export layout_x=6; export layout_y=8;  export NODES=24; 
        export master_grid=1deg;   export REMAP_TASKS=48 ;;
  C384) export DELTIM=450 ; export layout_x=12; export layout_y=16;  export NODES=96;
        export master_grid=0p5deg; export REMAP_TASKS=96 ;;
  C768) export DELTIM=225 ; export layout_x=16; export layout_y=24; export NODES=192;
        export master_grid=0p25deg; export REMAP_TASKS=384 ;;
  *)    echo "grid $CASE not supported, exit"
        exit ;;
esac

export PARM_FV3DIAG=$FV3DIR_RELEASE/parm/parm_fv3diag
export FORECASTSH=$FV3DIR_RELEASE/scripts/exglobal_fcst_nemsfv3gfs.sh         

#---determine task configuration
export nth_f=2                                         # number of threads 
export npe_node_f=24                                   # number of pes per node 
export task_per_node=$((npe_node_f/nth_f))
export tasks=$((NODES*task_per_node))                  # number of tasks 
export NTHREADS_REMAP=$nth_f


export MODE=32bit           			       # dycore precision:   32bit, 64bit
export TYPE=nh         				       # hydrostatic option: nh, hydro
export HYPT=off           			       # hyperthreading:     on, off  
export COMP="prod"        			       # compiling mode:     debug, repro, prod
if [ ${HYPT} = on ]; then
   export hyperthread=".true."
   export j_opt="-j 2"
else
   export hyperthread=".false."
   export j_opt="-j 1"
fi
export FCSTEXEC=fv3_gfs_${TYPE}.${COMP}.${MODE}.x

cp $FV3DIR/NEMS/src/conf/module-setup.sh.inc module-setup.sh
cp $FV3DIR/NEMS/src/conf/modules.nems modules.fv3
source ./module-setup.sh
module use $( pwd -P )
module load modules.fv3
module list

export mpiexec=`which mpirun`
export FCST_LAUNCHER="$mpiexec -prepend-rank -np $PBS_NP"

echo "Model started:  " `date`
export MPI_TYPE_DEPTH=20
export OMP_STACKSIZE=512M
export OMP_NUM_THREADS=$nth_f
export ESMF_RUNTIME_COMPLIANCECHECK=OFF:depth=4
#--------------------------------------------------------------------------

#--execute the forecast
$FORECASTSH
if [ $? != 0 ]; then echo "forecast failed, exit"; exit; fi

echo "fcst job is done"

#-------------------------------------------------------------------------
#--convert 6-tile output to global arrary in netCDF format
ymd=`echo $CDATE |cut -c 1-8`
cyc=`echo $CDATE |cut -c 9-10`
export DATA=$ROTDIR/gfs.$ymd/$cyc
export IPD4=YES
export REMAPSH=$FV3DIR_RELEASE/ush/fv3gfs_remap.sh            #remap 6-tile output to global array in netcdf
export REMAPEXE=$FV3DIR_RELEASE/exec/fregrid_parallel
export REMAP_LAUNCHER="mpirun -prepend-rank -np $REMAP_TASKS"

cp $FV3DIR_RELEASE/modulefiles/fv3gfs/fre-nctools.${machine} module.fre-nctools
module load module.fre-nctools
module list

$REMAPSH

echo "Remap job is done!"

exit

