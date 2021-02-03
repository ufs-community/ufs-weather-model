
###############################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
###############################################################################

if [[ $MACHINE_ID = wcoss_cray ]]; then

  TASKS_dflt=150 ; TPN_dflt=24 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=12 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=24 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=24 ; INPES_strnest=2 ; JNPES_strnest=4

elif [[ $MACHINE_ID = wcoss_dell_p3 || $MACHINE_ID = wcoss2 ]]; then

  TASKS_dflt=150 ; TPN_dflt=28 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=14 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=28 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=28 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_dflt=192; TPN_cpl_dflt=28; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 179"; IPB_cpl_dflt="180 191"

  TASKS_cpl_thrd=120; TPN_cpl_thrd=14; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 77";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 107";  IPB_cpl_thrd="108 119"

  TASKS_cpl_bmrk=480; TPN_cpl_bmrk=28; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"

  TASKS_cpl_wwav=520; TPN_cpl_wwav=28; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=1; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 519"

  TASKS_cpl_c192=288; TPN_cpl_c192=28; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"

  TASKS_cpl_c384=318; TPN_cpl_c384=28; INPES_cpl_c384=3; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=6;  MPB_cpl_c384="0 143"; APB_cpl_c384="0 149"
  OPB_cpl_c384="150 269"; IPB_cpl_c384="270 317"

  TASKS_datm_100=120; TPN_datm_100=28
  MPB_datm_100="16 77"; APB_datm_100="0 15"
  OPB_datm_100="78 107"; IPB_datm_100="108 119"

  TASKS_datm_025=208; TPN_datm_025=28
  MPB_datm_025="0 39"; APB_datm_025="0 39"
  OPB_datm_025="40 159"; IPB_datm_025="160 207"

elif [[ $MACHINE_ID = orion.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=20 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_dflt=192; TPN_cpl_dflt=40; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 179"; IPB_cpl_dflt="180 191"

  TASKS_cpl_thrd=120; TPN_cpl_thrd=40; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 77";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 107";  IPB_cpl_thrd="108 119"

  TASKS_cpl_bmrk=480; TPN_cpl_bmrk=40; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"

  TASKS_cpl_wwav=520; TPN_cpl_wwav=40; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=1; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 519"

  TASKS_cpl_c192=288; TPN_cpl_c192=40; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"

  TASKS_cpl_c384=318; TPN_cpl_c384=40; INPES_cpl_c384=3; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=6;  MPB_cpl_c384="0 143"; APB_cpl_c384="0 149"
  OPB_cpl_c384="150 269"; IPB_cpl_c384="270 317"

  TASKS_datm_100=120; TPN_datm_100=40
  MPB_datm_100="16 77"; APB_datm_100="0 15"
  OPB_datm_100="78 107"; IPB_datm_100="108 119"

  TASKS_datm_025=208; TPN_datm_025=40
  MPB_datm_025="0 39"; APB_datm_025="0 39"
  OPB_datm_025="40 159"; IPB_datm_025="160 207"

elif [[ $MACHINE_ID = hera.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=20 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4
  TASKS_strnest=96 ; TPN_strnest=12 ; INPES_strnest=2 ; JNPES_strnest=4

  TASKS_cpl_dflt=192; TPN_cpl_dflt=40; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 179"; IPB_cpl_dflt="180 191"

  TASKS_cpl_thrd=120; TPN_cpl_thrd=40; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 77";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 107";  IPB_cpl_thrd="108 119"

  TASKS_cpl_bmrk=480; TPN_cpl_bmrk=40; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"

  TASKS_cpl_wwav=520; TPN_cpl_wwav=40; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=1; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 519"

  TASKS_cpl_c192=288; TPN_cpl_c192=40; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"

  TASKS_cpl_c384=318; TPN_cpl_c384=40; INPES_cpl_c384=3; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=6;  MPB_cpl_c384="0 143"; APB_cpl_c384="0 149"
  OPB_cpl_c384="150 269"; IPB_cpl_c384="270 317"

  TASKS_datm_100=120; TPN_datm_100=40
  MPB_datm_100="16 77"; APB_datm_100="0 15"
  OPB_datm_100="78 107"; IPB_datm_100="108 119"

  TASKS_datm_025=208; TPN_datm_025=40
  MPB_datm_025="0 39"; APB_datm_025="0 39"
  OPB_datm_025="40 159"; IPB_datm_025="160 207"

elif [[ $MACHINE_ID = linux.* ]]; then

  if [[ $CI_TEST = true ]]; then
  TASKS_dflt=12 ; TPN_dflt=16 ; INPES_dflt=1 ; JNPES_dflt=1
  else
  TASKS_dflt=150 ; TPN_dflt=40 ; INPES_dflt=3 ; JNPES_dflt=8
  fi
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

  TASKS_cpl_dflt=192; TPN_cpl_dflt=36; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 179"; IPB_cpl_dflt="180 191"

  TASKS_cpl_thrd=120; TPN_cpl_thrd=36; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 77";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 107";  IPB_cpl_thrd="108 119"

  TASKS_cpl_bmrk=480; TPN_cpl_bmrk=36; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"

  TASKS_cpl_wwav=520; TPN_cpl_wwav=36; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=1; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 519"

  TASKS_cpl_c192=288; TPN_cpl_c192=36; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"

  TASKS_cpl_c384=318; TPN_cpl_c384=36; INPES_cpl_c384=3; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=6;  MPB_cpl_c384="0 143"; APB_cpl_c384="0 149"
  OPB_cpl_c384="150 269"; IPB_cpl_c384="270 317"

  TASKS_datm_100=120; TPN_datm_100=36
  MPB_datm_100="16 77"; APB_datm_100="0 15"
  OPB_datm_100="78 107"; IPB_datm_100="108 119"

  TASKS_datm_025=208; TPN_datm_025=36
  MPB_datm_025="0 39"; APB_datm_025="0 39"
  OPB_datm_025="40 159"; IPB_datm_025="160 207"

elif [[ $MACHINE_ID = stampede.* ]]; then

  TASKS_dflt=150 ; TPN_dflt=48 ; INPES_dflt=3 ; JNPES_dflt=8
  TASKS_thrd=84  ; TPN_thrd=24 ; INPES_thrd=3 ; JNPES_thrd=4
  TASKS_stretch=48 ; TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4

  TASKS_cpl_dflt=192; TPN_cpl_dflt=48; INPES_cpl_dflt=3; JNPES_cpl_dflt=8
  THRD_cpl_dflt=1; WPG_cpl_dflt=6;  MPB_cpl_dflt="0 143"; APB_cpl_dflt="0 149"
  OPB_cpl_dflt="150 179"; IPB_cpl_dflt="180 191"

  TASKS_cpl_thrd=120; TPN_cpl_thrd=48; INPES_cpl_thrd=3; JNPES_cpl_thrd=4
  THRD_cpl_thrd=2; WPG_cpl_thrd=6;  MPB_cpl_thrd="0 77";  APB_cpl_thrd="0 77"
  OPB_cpl_thrd="78 107";  IPB_cpl_thrd="108 119"

  TASKS_cpl_bmrk=480; TPN_cpl_bmrk=48; INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8
  THRD_cpl_bmrk=1; WPG_cpl_bmrk=24; MPB_cpl_bmrk="0 287"; APB_cpl_bmrk="0 311"
  OPB_cpl_bmrk="312 431"; IPB_cpl_bmrk="432 479"

  TASKS_cpl_wwav=520; TPN_cpl_wwav=48; INPES_cpl_wwav=6; JNPES_cpl_wwav=8
  THRD_cpl_wwav=1; WPG_cpl_wwav=24; MPB_cpl_wwav="0 287"; APB_cpl_wwav="0 311"
  OPB_cpl_wwav="312 431"; IPB_cpl_wwav="432 479"; WPB_cpl_wwav="480 519"

  TASKS_cpl_c192=288; TPN_cpl_c192=40; INPES_cpl_c192=4; JNPES_cpl_c192=8
  THRD_cpl_c192=1; WPG_cpl_c192=12;  MPB_cpl_c192="0 191"; APB_cpl_c192="0 203"
  OPB_cpl_c192="204 263"; IPB_cpl_c192="264 287"

  TASKS_cpl_c384=318; TPN_cpl_c384=48; INPES_cpl_c384=3; JNPES_cpl_c384=8
  THRD_cpl_c384=1; WPG_cpl_c384=6;  MPB_cpl_c384="0 143"; APB_cpl_c384="0 149"
  OPB_cpl_c384="150 269"; IPB_cpl_c384="270 317"

  TASKS_datm_100=120; TPN_datm_100=48
  MPB_datm_100="16 77"; APB_datm_100="0 15"
  OPB_datm_100="78 107"; IPB_datm_100="108 119"

  TASKS_datm_025=208; TPN_datm_025=48
  MPB_datm_025="0 39"; APB_datm_025="0 39"
  OPB_datm_025="40 159"; IPB_datm_025="160 207"

else

  echo "Unknown MACHINE_ID ${MACHINE_ID}"
  exit 1

fi

# Longer default walltime for GNU
if [[ ${RT_COMPILER:-} = gnu ]]; then
  WLCLK_dflt=30
else
  WLCLK_dflt=15
fi

export_fv3 ()
{
export FV3=true
export S2S=false
export DATM=false
export CDEPS_DATM=false
export CDEPS_DOCN=false
export THRD=1
export WLCLK=$WLCLK_dflt
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt
export RESTART_INTERVAL=0
export QUILTING=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'dyn' 'phy'"
export OUTPUT_GRID="'cubed_sphere_grid'"
export OUTPUT_FILE="'netcdf'"
export IDEFLATE=0
export NBITS=0
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
export DO_RRTMGP=.F.
export ICLOUD=0
export IOVR=1

# Microphysics
export IMP_PHYSICS=11
export NWAT=6
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
export DO_UGWP_V0=.F.
export DO_UGWP_V0_OROG_ONLY=.F.
export DO_GSL_DRAG_LS_BL=.F.
export DO_GSL_DRAG_SS=.F.
export DO_GSL_DRAG_TOFD=.F.
export DO_UGWP_V1=.F.
export DO_UGWP_V1_OROG_ONLY=.F.

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
export DO_MYNNSFCLAY=.F.

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
export FDIAG=0,1,2,3,4,5,6,7,8,9,10,11,12,15,18,21,24,27,30,33,36,39,42,45,48
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
export FHROT=0
export LDIAG3D=.F.
export QDIAG3D=.F.
export MAX_OUTPUT_FIELDS=300

# Stochastic physics
export DO_SPPT=.F.
export DO_SHUM=.F.
export DO_SKEB=.F.
export LNDP_TYPE=0
export N_VAR_LNDP=0
export SKEB=-999.
export SPPT=-999.
export SHUM=-999.

#IAU
export IAU_INC_FILES="''"

#Cellular automata
export DO_CA=.F.
export CA_SGS=.F.
export CA_GLOBAL=.F.

export IAU_DRYMASSFIXER=.false.

# Regional
export WRITE_RESTART_WITH_BCS=.false.

export coupling_interval_fast_sec=0
}

export_cpl ()
{
export FV3=true
export S2S=true
export DATM=false

export DAYS="1"
export FHMAX="24"
export FDIAG="6"
export WLCLK=30

# default atm/ocn/ice resolution
export ATMRES='C96'
export OCNRES='100'
export ICERES='1.00'
export NX_GLB=360
export NY_GLB=320

# default resources
export TASKS=$TASKS_cpl_dflt
export TPN=$TPN_cpl_dflt
export INPES=$INPES_cpl_dflt
export JNPES=$JNPES_cpl_dflt
export THRD=$THRD_cpl_dflt
export WRTTASK_PER_GROUP=$WPG_cpl_dflt

export med_petlist_bounds=$MPB_cpl_dflt
export atm_petlist_bounds=$APB_cpl_dflt
export ocn_petlist_bounds=$OPB_cpl_dflt
export ice_petlist_bounds=$IPB_cpl_dflt

# component and coupling timesteps
export DT_ATMOS='900'
export DT_CICE=${DT_ATMOS}
export DT_DYNAM_MOM6='1800'
export DT_THERM_MOM6='3600'
export CPL_SLOW=${DT_THERM_MOM6}
export CPL_FAST=${DT_ATMOS}

# nems.configure defaults
export NEMS_CONFIGURE="nems.configure.cpld.IN"
export med_model="nems"
export atm_model="fv3"
export ocn_model="mom6"
export ice_model="cice6"
export wav_model="ww3"

export coupling_interval_slow_sec=${CPL_SLOW}
export coupling_interval_fast_sec=${CPL_FAST}

export RESTART_N=${FHMAX}
export CPLMODE='nems_orig'
export cap_dbug_flag="0"
export use_coldstart="false"
export RUNTYPE='startup'

# FV3 defaults
# to use new oro and ics created from 1deg ocean mask on c96 tiles
# set frac_grid=.F. but FRAC_GRID_INPUT=.T.
# to repro existing tests set both frac_grid and FRAC_GRID_INPUT to .F.
# to run frac_grid, set both frac_grid and FRAC_GRID_INPUTs to .T.
export FRAC_GRID='.F.'
export FRAC_GRID_INPUT='.T.'
export SUITE_NAME="FV3_GFS_2017_coupled"
export INPUT_NML=input.mom6_ccpp.nml.IN
export FIELD_TABLE="field_table"

export FHROT='0'
export NSOUT='-1'
export FDIAG='6'
export NFHOUT='6'
#no high freq fv3 output
export NFHMAX_HF='-1'
export NFHOUT_HF='-1'

export CPLFLX='.T.'
export CPL='.true.'
export NSTF_NAME='0,0,0,0,0'

# for FV3: default values will be changed if doing a warm-warm restart
export WARM_START='.F.'
export MAKE_NH='.T.'
export NA_INIT='1'
export EXTERNAL_IC='.T.'
export NGGPS_IC='.T.'
export MOUNTAIN='.F.'

# MOM6 defaults; 1 degree
export MOM_INPUT=MOM_input_template_100
export MOM6_RESTART_SETTING='n'
export MOM6_RIVER_RUNOFF='False'
export FRUNOFF=""
export CHLCLIM="seawifs_1998-2006_smoothed_2X.nc"
# this must be set False for restart repro
export MOM6_REPRO_LA='False'
# since CPL_SLOW is set to DT_THERM, this should be always be false
export MOM6_THERMO_SPAN='False'
# no WW3
export MOM6_USE_WAVES='False'

# CICE6 defaults; 1 degree
export NPROC_ICE='12'
export MESHICE="mesh.mx${OCNRES}.nc"
export CICEGRID="grid_cice_NEMS_mx${OCNRES}.nc"
export CICEMASK="kmtu_cice_NEMS_mx${OCNRES}.nc"
export RUNID='unknown'
# set large; restart frequency now controlled by restart_n in nems.configure
export DUMPFREQ='d'
export DUMPFREQ_N=1000
export USE_RESTART_TIME='.false.'
export RESTART_EXT='.false.'
# setting to true will allow Frazil FW and Salt to be
# included in fluxes sent to ocean
export FRAZIL_FWSALT='.true.'
# default to write CICE average history files
export CICE_HIST_AVG='.true.'

# checkpoint restarts
export RESTART_FILE_PREFIX=''
export RESTART_FILE_SUFFIX_HRS=''
export RESTART_FILE_SUFFIX_SECS=''
export RT35D=''
}
export_35d_run ()
{
export CNTL_DIR=""
export LIST_FILES=""
}
export_datm ()
{
export FV3=false
export S2S=false
export DATM=true
export CDEPS_DATM=false
export CDEPS_DOCN=false
export DAYS=1
export FHMAX=24
export WLCLK=30
export THRD=1
export FHROT='0'
export WARM_START=.F.

# atm/ocn/ice resolution
# GEFS
export DATM_SRC="GEFS"
export FILENAME_BASE='gefs.'
export IATM=1536
export JATM=768
export ATMRES='C96'
export OCNRES='100'
export ICERES='1.00'
export NX_GLB=360
export NY_GLB=320

# nems.configure
export NEMS_CONFIGURE="nems.configure.datm.IN"
export med_model="nems"
export atm_model="datm"
export ocn_model="mom6"
export ice_model="cice6"
export atm_petlist_bounds=$APB_datm_100
export med_petlist_bounds=$MPB_datm_100
export ocn_petlist_bounds=$OPB_datm_100
export ice_petlist_bounds=$IPB_datm_100
export TASKS=$TASKS_datm_100
export TPN=$TPN_datm_100
export NPROC_ICE='12'

export ENS_NUM=1
export SYEAR='2011'
export SMONTH='10'
export SDAY='01'
export SHOUR='00'
export CDATE=${SYEAR}${SMONTH}${SDAY}${SHOUR}

export NFHOUT=6
export FDIAG=6
export DT_ATMOS='900'
export DT_CICE=${DT_ATMOS}
export DT_DYNAM_MOM6='1800'
export DT_THERM_MOM6='3600'
export CPL_SLOW=${DT_THERM_MOM6}
export CPL_FAST=${DT_ATMOS}
export coupling_interval_slow_sec=${CPL_SLOW}
export coupling_interval_fast_sec=${CPL_FAST}

export RESTART_N=${FHMAX}
export CPLMODE='nems_orig_data'
export cap_dbug_flag="0"
export use_coldstart=".false."
export RUNTYPE='startup'
export flux_convergence='0.0'
export flux_iteration='2'
export flux_scheme='0'

export INPUT_NML=input.mom6.nml.IN
export MODEL_CONFIGURE=datm_configure.IN
export FIELD_TABLE="field_table"

# MOM6 defaults; 1 degree
export MOM_INPUT=MOM_input_template_100
export MOM6_RESTART_SETTING='n'
export MOM6_RIVER_RUNOFF='False'
export FRUNOFF=""
export CHLCLIM='"seawifs_1998-2006_smoothed_2X.nc"'
# this must be set False for restart repro
export MOM6_REPRO_LA='False'
# since CPL_SLOW is set to DT_THERM, this should be always be false
export MOM6_THERMO_SPAN='False'
# no WW3
export MOM6_USE_WAVES='False'

# CICE6 defaults; 1 degree
export MESHICE="mesh.mx${OCNRES}.nc"
export CICEGRID="grid_cice_NEMS_mx${OCNRES}.nc"
export CICEMASK="kmtu_cice_NEMS_mx${OCNRES}.nc"
export RUNID='unknown'
# set large; restart frequency now controlled by restart_n in nems.configure
export DUMPFREQ='d'
export DUMPFREQ_N=1000
export USE_RESTART_TIME='.false.'
export RESTART_EXT='.false.'
# setting to true will allow Frazil FW and Salt to be
# included in fluxes sent to ocean
export FRAZIL_FWSALT='.true.'
# default to write CICE average history files
export CICE_HIST_AVG='.true.'
export BL_SUFFIX=""
export RT_SUFFIX=""
}
export_cdeps_datm ()
{
export FV3=false
export CDEPS_DATM=true
export DATM=false
export THRD=1
export WLCLK=$WLCLK_dflt
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt

export atm_model="datm"
}
export_cdeps_docn ()
{
export CDEPS_DOCN=true
export THRD=1
export WLCLK=$WLCLK_dflt
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export TASKS=$TASKS_dflt
export TPN=$TPN_dflt

export ocn_model="docn"
export ocn_datamode="sstdata"
}
export_cpl_regional ()
{
export S2S=false

# model_configure
export SYEAR='2019'
export SMONTH='08'
export SDAY='29'
export SHOUR='00'
export FHMAX=24
export ENS_NUM=1
export DT_ATMOS='900'
export CPL='.true.'
export RESTART_INTERVAL=0
export FHROT=0
export coupling_interval_fast_sec=0
export QUILTING=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'dyn' 'phy'"
export OUTPUT_GRID="'regional_latlon'"
export OUTPUT_FILE="'netcdf'"
export IDEFLATE=0
export NBITS=0
export WRITE_NEMSIOFLIP=.false.
export WRITE_FSYNCFLAG=.false.
export NFHOUT=3
export NFHMAX_HF=-1
export NFHOUT_HF=3

# nems.configure
export med_model="nems"
export RESTART_N=${FHMAX} 
export CPLMODE="hafs"
export RUNTYPE="startup"
export USE_COLDSTART="false"
}
