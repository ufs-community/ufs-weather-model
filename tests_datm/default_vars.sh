
###############################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
###############################################################################

if [ $MACHINE_ID = wcoss_cray ]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4

elif [ $MACHINE_ID = wcoss_dell_p3 ]; then

  TASKS_dflt=150 ; TPN_dflt=28 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=14 ; INPES_thrd=3 ; JNPES_thrd=4

elif [[ $MACHINE_ID = orion.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4

  TASKS_cpl_dflt=280; TPN_cpl_dflt=40; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="16 111"; APB_cpl_dflt="0 15"
  OPB_cpl_dflt="112 231"; IPB_cpl_dflt="232 279"

  TASKS_cpl_thrd=246; TPN_cpl_thrd=40; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 77";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 197";  IPB_cpl_thrd="198 245"

  TASKS_cpl_wwav=520; TPN_cpl_wwav=40; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=1; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 519"

elif [[ $MACHINE_ID = hera.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=20 ; INPES_thrd=3 ; JNPES_thrd=4

  TASKS_cpl_dflt=280; TPN_cpl_dflt=40; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="16 111"; APB_cpl_dflt="0 15"
  OPB_cpl_dflt="112 231"; IPB_cpl_dflt="232 279"

  TASKS_cpl_thrd=246; TPN_cpl_thrd=40; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 77";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 197";  IPB_cpl_thrd="198 245"

  TASKS_cpl_wwav=520; TPN_cpl_wwav=40; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=1; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 519"

elif [[ $MACHINE_ID = jet.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4

elif [[ $MACHINE_ID = gaea.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=36 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=18 ; INPES_thrd=3 ; JNPES_thrd=4

elif [[ $MACHINE_ID = cheyenne.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=36 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=18 ; INPES_thrd=3 ; JNPES_thrd=4

elif [[ $MACHINE_ID = stampede.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=48 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=24 ; INPES_thrd=3 ; JNPES_thrd=4

else

  echo "Unknown MACHINE_ID ${MACHINE_ID}"
  exit 1

fi

# Re-instantiate COMPILER in case it gets deleted by module purge
COMPILER=${NEMS_COMPILER:-intel}

# Longer default walltime for GNU and PGI
if [[ $COMPILER = gnu ]] || [[ $COMPILER = pgi ]]; then
    WLCLK_dflt=30
else
    WLCLK_dflt=30
fi

export_datm ()
{
export THRD=1
export WLCLK=$WLCLK_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt
export RESTART_INTERVAL=840
export CPL=.F.
export CPLFLX=.F.
export CPLWAV=.F.
export CPLWAV2ATM=.F.
export IATM=1536
export JATM=768
export WARM_START=.F.
export ENS_NUM=1
export SYEAR='2011'
export SMONTH='10'
export SDAY='01'
export SHOUR='00'
export CDATE=${SYEAR}${SMONTH}${SDAY}${SHOUR}
export NFHOUT=6
export DAYS=0.041666666
export FHMAX=1
export FHMAX=${FHMAX:-`expr $DAYS \* 24`}
export DT_ATMOS=900
export DATM_SRC="GEFS"
export FILENAME_BASE='gefs.'
}

export_cpl ()
{ 
export TASKS=$TASKS_cpl_dflt
export TPN=$TPN_cpl_dflt
export THRD=$THRD_cpl_dflt

export med_petlist_bounds=$MPB_cpl_dflt
export atm_petlist_bounds=$APB_cpl_dflt
export ocn_petlist_bounds=$OPB_cpl_dflt
export ice_petlist_bounds=$IPB_cpl_dflt

export OCNRES='025'
export INPUT_NML="input.mom6.nml.IN"
export FIELD_TABLE="field_table"
export NSOUT='-1'
export CPLFLX='.T.'
export CPL='.true.'
export MOM6_RESTART_SETTING='r'
export MOM6_RIVER_RUNOFF='True'
export DT_DYNAM_MOM6='900'
export DT_THERM_MOM6='1800'

export med_model="nems"
export atm_model="datm"
export ocn_model="mom6"
export ice_model="cice6"
export cap_dbug_flag="0"
export use_coldstart=".true."

export NPROC_ICE='48'
# defaults for CICE runtype and restart writing
export RUNTYPE='startup' 
export RUNID='unknown'
export DUMPFREQ='d' 
export DUMPFREQ_N='1' 
export NX_GLB='1440'
export NY_GLB='1080'
export BLCKX='60'
export BLCKY='540'
export grid_cice_NEMS_mx='grid_cice_NEMS_mx025.nc' 
export kmtu_cice_NEMS_mx='kmtu_cice_NEMS_mx025.nc' 
export USE_RESTART_TIME='.false.'
export MESHICE="mesh.mx025.nc"
# setting to true will allow Frazil FW and Salt to be
# included in fluxes sent to ocean
export FRAZIL_FWSALT='.true.'
# default to write CICE average history files
export CICE_HIST_AVG='.true.'
export MEDCOMP=''
}
