
###############################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
###############################################################################

if [ $MACHINE_ID = wcoss ]; then

  TASKS_dflt=144 ; TPN_dflt=16 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

elif [ $MACHINE_ID = wcoss_cray ]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=24 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=24 ; INPES_strnest=2 ; JNPES_strnest=4

elif [ $MACHINE_ID = wcoss_dell_p3 ]; then

  TASKS_dflt=150 ; TPN_dflt=28 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=14 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=28 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=28 ; INPES_strnest=2 ; JNPES_strnest=4

elif [[ $MACHINE_ID = hera.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=20 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

elif [[ $MACHINE_ID = jet.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

elif [[ $MACHINE_ID = gaea.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=36 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=18 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=18 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=18 ; INPES_strnest=2 ; JNPES_strnest=4

elif [[ $MACHINE_ID = cheyenne.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=36 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=18 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=18 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=18 ; INPES_strnest=2 ; JNPES_strnest=4

elif [[ $MACHINE_ID = stampede.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=48 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=24 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4

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
    WLCLK_dflt=15
fi

export_fv3 ()
{
export THRD=1
export WLCLK=$WLCLK_dflt
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt
export QUILTING=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'dyn' 'phy'"
export OUTPUT_GRID="'cubed_sphere_grid'"
export OUTPUT_FILE="'netcdf'"
export WRITE_NEMSIOFLIP=.false.
export WRITE_FSYNCFLAG=.false.
export IMO=384
export JMO=190

# Coldstart/warmstart
export WARM_START=.false.
export READ_INCREMENT=.false.
export NGGPS_IC=.true.
export EXTERNAL_IC=.true.
export MAKE_NH=.true.
export MOUNTAIN=.false.
export NA_INIT=1

# Microphysics
export IMP_PHYSICS=11
# GFDL MP
export DNATS=1
export DO_SAT_ADJ=.true.
export LHEATSTRG=.false.
export LGFDLMPRAD=.false.
export EFFR_IN=.false.
# Thompson MP
export LRADAR=.true.
export LTAEROSOL=.true.

# GWD
export LDIAG_UGWP=.false.
export DO_UGWP=.false.
export DO_TOFD=.false.
export GWD_OPT=1

# PBL
export SATMEDMF=.false.
export ISATMEDMF=0
export HYBEDMF=.true.
export SHINHONG=.false.
export DO_YSU=.false.
export DO_MYNNEDMF=.false.
export DO_MYJPBL=.false.

# Shallow/deep convection
export IMFSHALCNV=2
export IMFDEEPCNV=2

# SFC
export DO_MYJSFC=.false.

# LSM
export LSM=1
export LSOIL_LSM=4
export LANDICE=.true.

# Ozone / stratospheric H2O
export OZ_PHYS_OLD=.true.
export OZ_PHYS_NEW=.false.
export H2O_PHYS=.false.

export CPL=.false.
export CPLFLX=.false.
export CPLWAV=.false.
export DAYS=1
export NPX=97
export NPY=97
export NPZ=64
export NPZP=65
export NSTF_NAME=2,1,1,0,5
export FDIAG=0,1,2,3,4,5,6,7,8,9,10,11,12,15,18,21,24
export NFHOUT=3
export NFHMAX_HF=12
export NFHOUT_HF=1
export FNALBC="'global_snowfree_albedo.bosu.t126.384.190.rg.grb',"
export FNVETC="'global_vegtype.igbp.t126.384.190.rg.grb',"
export FNSOTC="'global_soiltype.statsgo.t126.384.190.rg.grb',"
export FNSMCC="'global_soilmgldas.t126.384.190.grb',"
export FNABSC="'global_mxsnoalb.uariz.t126.384.190.rg.grb',"


export ENS_NUM=1
export SYEAR=2016
export SMONTH=10
export SDAY=03
export SHOUR=00
export FHMAX=${FHMAX:-`expr $DAYS \* 24`}
export DT_ATMOS=1800
export FHCYC=24

# Stochastic physics
export DO_SPPT=.false.
export DO_SHUM=.false.
export DO_SKEB=.false.
export DO_SFCPERTS=.false.
export SKEB=-999.
export SPPT=-999.
export SHUM=-999.
export IAU_INC_FILES="''"
}
