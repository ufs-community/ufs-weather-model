
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

elif [[ $MACHINE_ID = theia.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

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
export WARM_START=.F.
export READ_INCREMENT=.F.
export NGGPS_IC=.T.
export EXTERNAL_IC=.T.
export MAKE_NH=.T.
export MOUNTAIN=.F.
export NA_INIT=1

# Radiation
export ICLOUD=0
export IOVR_LW=1
export IOVR_SW=1

# Microphysics
export IMP_PHYSICS=11
# GFDL MP
export DNATS=1
export DO_SAT_ADJ=.T.
export LHEATSTRG=.F.
export LGFDLMPRAD=.F.
export EFFR_IN=.F.
# Thompson MP
export LRADAR=.T.
export LTAEROSOL=.T.

# GWD
export LDIAG_UGWP=.F.
export DO_UGWP=.F.
export DO_TOFD=.F.
export GWD_OPT=1

# PBL
export SATMEDMF=.F.
export ISATMEDMF=0
export HYBEDMF=.T.
export SHINHONG=.F.
export DO_YSU=.F.
export DO_MYNNEDMF=.F.
export DO_MYJPBL=.F.
export HURR_PBL=.F.
export MONINQ_FAC=1.0

# Shallow/deep convection
export IMFSHALCNV=2
export HWRF_SAMFSHAL=.F.
export IMFDEEPCNV=2
export HWRF_SAMFDEEP=.F.

# SFC
export DO_MYJSFC=.F.

# LSM
export LSM=1
export LSOIL_LSM=4
export LANDICE=.T.

# Ozone / stratospheric H2O
export OZ_PHYS_OLD=.T.
export OZ_PHYS_NEW=.F.
export H2O_PHYS=.F.

export CPL=.F.
export CPLFLX=.F.
export CPLWAV=.F.
export CPLWAV2ATM=.F.
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
export DO_SPPT=.F.
export DO_SHUM=.F.
export DO_SKEB=.F.
export DO_SFCPERTS=.F.
export SKEB=-999.
export SPPT=-999.
export SHUM=-999.

#IAU
export IAU_INC_FILES="''"
export IAU_DRYMASSFIXER=.false.
}
