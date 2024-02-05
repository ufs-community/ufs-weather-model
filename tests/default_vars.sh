
###############################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
###############################################################################

  THRD=1

  INPES_atmaero=4; JNPES_atmaero=8; WPG_atmaero=6

  THRD_cpl_atmw=1
  INPES_cpl_atmw=3; JNPES_cpl_atmw=8; WPG_cpl_atmw=6
  WAV_tasks_cpl_atmw=30
  WAV_thrds_cpl_atmw=1

  THRD_cpl_c48=1
  INPES_cpl_c48=1; JNPES_cpl_c48=1; WPG_cpl_c48=6
  OCN_tasks_cpl_c48=4
  ICE_tasks_cpl_c48=4

  THRD_cpl_dflt=1
  INPES_cpl_dflt=3; JNPES_cpl_dflt=8; WPG_cpl_dflt=6
  OCN_tasks_cpl_dflt=20
  ICE_tasks_cpl_dflt=10
  WAV_tasks_cpl_dflt=20

  THRD_cpl_thrd=2
  INPES_cpl_thrd=3; JNPES_cpl_thrd=4; WPG_cpl_thrd=6
  OCN_tasks_cpl_thrd=20
  OCN_thrds_cpl_thrd=1
  ICE_tasks_cpl_thrd=10
  ICE_thrds_cpl_thrd=1
  WAV_tasks_cpl_thrd=12
  WAV_thrds_cpl_thrd=2

  THRD_cpl_dcmp=1
  INPES_cpl_dcmp=4; JNPES_cpl_dcmp=6; WPG_cpl_dcmp=6
  OCN_tasks_cpl_dcmp=20
  ICE_tasks_cpl_dcmp=10
  WAV_tasks_cpl_dcmp=20

  THRD_cpl_mpi=1
  INPES_cpl_mpi=4; JNPES_cpl_mpi=8; WPG_cpl_mpi=6
  OCN_tasks_cpl_mpi=34
  ICE_tasks_cpl_mpi=20
  WAV_tasks_cpl_mpi=28

  THRD_cpl_bmrk=2
  INPES_cpl_bmrk=8; JNPES_cpl_bmrk=8; WPG_cpl_bmrk=48
  OCN_tasks_cpl_bmrk=120
  OCN_thrds_cpl_bmrk=1
  ICE_tasks_cpl_bmrk=48
  ICE_thrds_cpl_bmrk=1
  WAV_tasks_cpl_bmrk=80
  WAV_thrds_cpl_bmrk=2

  THRD_cpl_c192=2
  INPES_cpl_c192=6; JNPES_cpl_c192=8; WPG_cpl_c192=12
  OCN_tasks_cpl_c192=60
  ICE_tasks_cpl_c192=24
  WAV_tasks_cpl_c192=80

  ATM_compute_tasks_cdeps_100=12
  OCN_tasks_cdeps_100=16
  ICE_tasks_cdeps_100=12

  ATM_compute_tasks_cdeps_025=40
  OCN_tasks_cdeps_025=120
  ICE_tasks_cdeps_025=48

  INPES_aqm=33; JNPES_aqm=8

  THRD_cpl_unstr=1
  INPES_cpl_unstr=3; JNPES_cpl_unstr=8; WPG_cpl_unstr=6
  OCN_tasks_cpl_unstr=20
  ICE_tasks_cpl_unstr=10
  WAV_tasks_cpl_unstr=60

  THRD_cpl_unstr_mpi=1
  INPES_cpl_unstr_mpi=4; JNPES_cpl_unstr_mpi=8; WPG_cpl_unstr_mpi=6
  OCN_tasks_cpl_unstr_mpi=34
  ICE_tasks_cpl_unstr_mpi=20
  WAV_tasks_cpl_unstr_mpi=50

  aqm_omp_num_threads=1
  atm_omp_num_threads=1
  chm_omp_num_threads=1
  ice_omp_num_threads=1
  lnd_omp_num_threads=1
  med_omp_num_threads=1
  ocn_omp_num_threads=1
  wav_omp_num_threads=1

if [[ $MACHINE_ID = wcoss2 || $MACHINE_ID = acorn ]]; then

  TPN=128

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=8 ; JNPES_c384=6  ; THRD_c384=2
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ $MACHINE_ID = orion ]]; then

  TPN=40

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=8 ; JNPES_c384=6  ; THRD_c384=2
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ $MACHINE_ID = hercules ]]; then

  TPN=80

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=8 ; JNPES_c384=6  ; THRD_c384=2
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248


elif [[ $MACHINE_ID = hera ]]; then

  TPN=40

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=6 ; JNPES_c384=8  ; THRD_c384=2
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=4

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ $MACHINE_ID = linux ]]; then

  TPN=40

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4

  THRD_cpl_dflt=1
  INPES_cpl_dflt=3; JNPES_cpl_dflt=8; WPG_cpl_dflt=6
  OCN_tasks_cpl_dflt=20
  ICE_tasks_cpl_dflt=10
  WAV_tasks_cpl_dflt=20

  THRD_cpl_thrd=2
  INPES_cpl_thrd=3; JNPES_cpl_thrd=4; WPG_cpl_thrd=6
  OCN_tasks_cpl_thrd=20
  ICE_tasks_cpl_thrd=10
  WAV_tasks_cpl_thrd=12

elif [[ $MACHINE_ID = jet ]]; then

  TPN=24

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=6 ; JNPES_c384=12 ; THRD_c384=1
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2
  WRTTASK_PER_GROUP_c384=84
  WRTTASK_PER_GROUP_c384gdas=88

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=240

  # run only in weekly test
  THRD_cpl_bmrk=2
  INPES_cpl_bmrk=16; JNPES_cpl_bmrk=16; WPG_cpl_bmrk=48
  OCN_tasks_cpl_bmrk=100
  ICE_tasks_cpl_bmrk=48
  WAV_tasks_cpl_bmrk=100
  WLCLK_cpl_bmrk=120

  # run only in weekly test
  THRD_cpl_c192=2
  INPES_cpl_c192=12; JNPES_cpl_c192=16; WPG_cpl_c192=24
  OCN_tasks_cpl_c192=100
  ICE_tasks_cpl_c192=48
  WAV_tasks_cpl_c192=80
  WLCLK_cpl_c192=120

elif [[ $MACHINE_ID = s4 ]]; then

  TPN=32

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=6 ; JNPES_c384=8 ; THRD_c384=2
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=1

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

  THRD_cpl_bmrk=2
  INPES_cpl_bmrk=6; JNPES_cpl_bmrk=8; WPG_cpl_bmrk=24
  OCN_tasks_cpl_bmrk=120
  ICE_tasks_cpl_bmrk=48
  WAV_tasks_cpl_bmrk=80

elif [[ $MACHINE_ID = gaea ]]; then

  TPN=128

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=6 ; JNPES_c384=8  ; THRD_c384=1
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

  THRD_cpl_atmw_gdas=3
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=264

elif [[ $MACHINE_ID = derecho ]]; then

  TPN=128
  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=8 ; JNPES_c384=6  ; THRD_c384=2
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ $MACHINE_ID = stampede ]]; then

  echo "Unknown MACHINE_ID ${MACHINE_ID}. Please update tasks configurations in default_vars.sh"
  exit 1

  TPN_dflt=48 ; INPES_dflt=3 ; JNPES_dflt=8
  TPN_thrd=24 ; INPES_thrd=3 ; JNPES_thrd=4
  TPN_c384=20 ; INPES_c384=8 ; JNPES_c384=6
  TPN_c768=20 ; INPES_c768=8 ; JNPES_c768=16
  TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4

  TPN_cpl_atmw_gdas=12; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=4; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

elif [[ ${MACHINE_ID} = noaacloud ]] ; then

    if [[ $PW_CSP == aws ]]; then
    TPN=36
    elif [[ $PW_CSP == azure ]]; then
    TPN=44
    elif [[ $PW_CSP == google ]]; then
    TPN=30
    fi

    INPES_dflt=3 ; JNPES_dflt=8
    INPES_thrd=3 ; JNPES_thrd=4

    INPES_c384=8 ; JNPES_c384=6  ; THRD_c384=2
    INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

    THRD_cpl_dflt=1
    INPES_cpl_dflt=3; JNPES_cpl_dflt=8; WPG_cpl_dflt=6
    OCN_tasks_cpl_dflt=20
    ICE_tasks_cpl_dflt=10
    WAV_tasks_cpl_dflt=20

    THRD_cpl_thrd=2
    INPES_cpl_thrd=3; JNPES_cpl_thrd=4; WPG_cpl_thrd=6
    OCN_tasks_cpl_thrd=20
    ICE_tasks_cpl_thrd=10
    WAV_tasks_cpl_thrd=12

elif [[ $MACHINE_ID = expanse ]]; then

  echo "Unknown MACHINE_ID ${MACHINE_ID}. Please update tasks configurations in default_vars.sh"
  exit 1

  TPN_dflt=64 ; INPES_dflt=3 ; JNPES_dflt=8
  TPN_thrd=64 ; INPES_thrd=3 ; JNPES_thrd=4
  TPN_stretch=12 ; INPES_stretch=2 ; JNPES_stretch=4

  TPN_cpl_atmw_gdas=12; INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8
  THRD_cpl_atmw_gdas=2; WPG_cpl_atmw_gdas=24; APB_cpl_atmw_gdas="0 311"; WPB_cpl_atmw_gdas="312 559"

else

  echo "Unknown MACHINE_ID ${MACHINE_ID}"
  exit 1

fi

WLCLK_dflt=30

export WLCLK=$WLCLK_dflt
export CMP_DATAONLY=false

# Defaults for ufs.configure
export esmf_logkind="ESMF_LOGKIND_MULTI"
export DumpFields="false"

export_fv3 ()
{
# ufs.configure defaults
export UFS_CONFIGURE=ufs.configure.atm.IN
export MODEL_CONFIGURE=model_configure.IN
export atm_model=fv3

export FV3=true
export S2S=false
export HAFS=false
export AQM=false
export DATM_CDEPS=false
export DOCN_CDEPS=false
export CDEPS_INLINE=false
export POSTAPP='global'
export USE_MERRA2=.false.

export NTILES=6
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export RESTART_INTERVAL=0
export QUILTING=.true.
export QUILTING_RESTART=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export ITASKS=1
export OUTPUT_HISTORY=.true.
export HISTORY_FILE_ON_NATIVE_GRID=.false.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'atm' 'sfc'"
export OUTPUT_GRID="'cubed_sphere_grid'"
export OUTPUT_FILE="'netcdf'"
export ZSTANDARD_LEVEL=0
export IDEFLATE=0
export QUANTIZE_NSD=0
export ICHUNK2D=0
export JCHUNK2D=0
export ICHUNK3D=0
export JCHUNK3D=0
export KCHUNK3D=0
export IMO=384
export JMO=190
export WRITE_NSFLIP=.false.

#input file
export DIAG_TABLE=diag_table_gfsv16
export FIELD_TABLE=field_table_gfsv16

export DOMAINS_STACK_SIZE=3000000

# Coldstart/warmstart
#rt script for ICs
export MODEL_INITIALIZATION=false
#namelist variable
export WARM_START=.false.
export READ_INCREMENT=.false.
export RES_LATLON_DYNAMICS="''"
export NGGPS_IC=.true.
export EXTERNAL_IC=.true.
export MAKE_NH=.true.
export MOUNTAIN=.false.
export NA_INIT=1

# Radiation
export DO_RRTMGP=.false.
export DOGP_CLDOPTICS_LUT=.false.
export DOGP_LWSCAT=.false.
export USE_LW_JACOBIAN=.false.
export DAMP_LW_FLUXADJ=.false.
export RRTMGP_LW_PHYS_BLKSZ=2
export ICLOUD=0
export IAER=111
export ICLIQ_SW=1
export IOVR=1
export LFNC_K=-999
export LFNC_P0=-999

# Microphysics
export IMP_PHYSICS=11
export NWAT=6
# GFDL MP
export DNATS=1
export DO_SAT_ADJ=.true.
export LHEATSTRG=.true.
export LSEASPRAY=.false.
export LGFDLMPRAD=.false.
export EFFR_IN=.false.
# Thompson MP
export LRADAR=.true.
export LTAEROSOL=.true.
export EXT_DIAG_THOMPSON=.false.
export SEDI_SEMI=.true.
export DECFL=10
# NSSL MP
export NSSL_CCCN=0.6e9
export NSSL_ALPHAH=0.0
export NSSL_ALPHAHL=1.0
export NSSL_HAIL_ON=.false.
export NSSL_CCN_ON=.true.
export NSSL_INVERTCCN=.true.

# Smoke
export RRFS_SMOKE=.false.
export SMOKE_FORECAST=0
export RRFS_RESTART=NO
export SEAS_OPT=2

# GWD
export LDIAG_UGWP=.false.
export DO_UGWP=.false.
export DO_TOFD=.false.
export GWD_OPT=1
export DO_UGWP_V0=.false.
export DO_UGWP_V1_W_GSLDRAG=.false.
export DO_UGWP_V0_OROG_ONLY=.false.
export DO_GSL_DRAG_LS_BL=.false.
export DO_GSL_DRAG_SS=.false.
export DO_GSL_DRAG_TOFD=.false.
export DO_UGWP_V1=.false.
export DO_UGWP_V1_OROG_ONLY=.false.
export KNOB_UGWP_DOKDIS=1
export KNOB_UGWP_NDX4LH=1
export KNOB_UGWP_VERSION=0
export KNOB_UGWP_PALAUNCH=500.e2
export KNOB_UGWP_NSLOPE=1

# resolution dependent settings
export CDMBWD_c48='0.071,2.1,1.0,1.0'
export CDMBWD_c96='0.14,1.8,1.0,1.0'
export CDMBWD_c192='0.23,1.5,1.0,1.0'
export CDMBWD_c384='1.1,0.72,1.0,1.0'
export CDMBWD_c768='4.0,0.15,1.0,1.0'

#DT_INNER=(Time step)/2
export DT_INNER_c96=360
export DT_INNER_c192=300
export DT_INNER_c384=150
export DT_INNER_c768=75

# set default
export CDMBWD=${CDMBWD_c96}
export DT_INNER=${DT_INNER_c96}

# PBL
export SATMEDMF=.false.
export ISATMEDMF=0
export HYBEDMF=.true.
export SHINHONG=.false.
export DO_YSU=.false.
export DO_MYNNEDMF=.false.
export HURR_PBL=.false.
export MONINQ_FAC=1.0
export SFCLAY_COMPUTE_FLUX=.false.

# Shallow/deep convection
export DO_DEEP=.true.
export SHAL_CNV=.true.
export IMFSHALCNV=2
export HWRF_SAMFSHAL=.false.
export IMFDEEPCNV=2
export HWRF_SAMFDEEP=.false.
export RAS=.false.
export RANDOM_CLDS=.false.
export CNVCLD=.true.
export PROGSIGMA=.false.
export BETASCU=8.0
export BETAMCU=1.0
export BETADCU=2.0

# Aerosol convective scavenging
export FSCAV_AERO='"*:0.3","so2:0.0","msa:0.0","dms:0.0","nh3:0.4","nh4:0.6","bc1:0.6","bc2:0.6","oc1:0.4","oc2:0.4","dust1:0.6","dust2:0.6","dust3:0.6","dust4:0.6","dust5:0.6","seas1:0.5","seas2:0.5","seas3:0.5","seas4:0.5","seas5:0.5"'

# SFC
export DO_MYJSFC=.false.
export DO_MYNNSFCLAY=.false.
export BL_MYNN_TKEADVECT=.false.

# LSM
export LSM=1
export LSOIL_LSM=4
export LANDICE=.true.
export KICE=2
export IALB=1
export IEMS=1

# Ozone / stratospheric H2O
export OZ_PHYS_OLD=.true.
export OZ_PHYS_NEW=.false.
export H2O_PHYS=.false.

# Lake models
export LKM=0 # 0=no lake, 1=run lake model, 2=run both lake and nsst on lake points
export IOPT_LAKE=2 # 1=flake, 2=clm lake
export LAKEFRAC_THRESHOLD=0.0 # lake fraction must be higher for lake model to run it
export LAKEDEPTH_THRESHOLD=1.0 # lake must be deeper (in meters) for a lake model to run it
export FRAC_ICE=.true. # should be false for flake, true for clm_lake

export CPL=.false.
export CPLCHM=.false.
export CPLFLX=.false.
export CPLICE=.false.
export CPLWAV=.false.
export CPLWAV2ATM=.false.
export CPLLND=.false.
export CPLLND2ATM=.false.
export USE_MED_FLUX=.false.
export DAYS=1
export NPX=97
export NPY=97
export NPZ=64
export NPZP=65
export NSTF_NAME=2,1,1,0,5
export OUTPUT_FH="12 -1"
export FHZERO=6
export FNALBC="'global_snowfree_albedo.bosu.t126.384.190.rg.grb'"
export FNVETC="'global_vegtype.igbp.t126.384.190.rg.grb'"
export FNSOTC="'global_soiltype.statsgo.t126.384.190.rg.grb'"
export FNSOCC="''"
export FNSMCC="'global_soilmgldas.t126.384.190.grb'"
export FNSMCC_control="'global_soilmgldas.statsgo.t1534.3072.1536.grb'"
export FNMSKH_control="'global_slmask.t1534.3072.1536.grb'"
export FNABSC="'global_mxsnoalb.uariz.t126.384.190.rg.grb'"

# Dynamical core
export FV_CORE_TAU=0.
export RF_CUTOFF=30.0
export FAST_TAU_W_SEC=0.0

# Tiled Fix files
export ATMRES=C96
export TILEDFIX=.false.

export ENS_NUM=1
export SYEAR=2016
export SMONTH=10
export SDAY=03
export SHOUR=00
export SECS=`expr $SHOUR \* 3600`
export FHMAX=$(( DAYS*24 ))
export DT_ATMOS=1800
export FHCYC=24
export FHROT=0
export LDIAG3D=.false.
export QDIAG3D=.false.
export PRINT_DIFF_PGR=.false.
export MAX_OUTPUT_FIELDS=310

# Stochastic physics
export STOCHINI=.false.
export DO_SPPT=.false.
export DO_SHUM=.false.
export DO_SKEB=.false.
export LNDP_TYPE=0
export N_VAR_LNDP=0
export SKEB=-999.
export SPPT=-999.
export SHUM=-999.
export LNDP_VAR_LIST="'XXX'"
export LNDP_PRT_LIST=-999
export LNDP_MODEL_TYPE=0

#IAU
export IAU_INC_FILES="''"
export IAU_DELTHRS=0
export IAUFHRS=-1
export IAU_OFFSET=0

export FH_DFI_RADAR='-2e10'

#Cellular automata
export DO_CA=.false.
export CA_SGS=.false.
export CA_GLOBAL=.false.

#waves
export WW3_RSTDTHR=12
export WW3_DT_2_RST="$(printf "%02d" $(( ${WW3_RSTDTHR}*3600 )))"
export WW3_OUTDTHR=1
export WW3_DTFLD="$(printf "%02d" $(( ${WW3_OUTDTHR}*3600 )))"
export WW3_DTPNT="$(printf "%02d" $(( ${WW3_OUTDTHR}*3600 )))"
export DTRST=0
export RSTTYPE=T
export GOFILETYPE=1
export POFILETYPE=1
export WW3_OUTPARS="WND HS FP DP PHS PTP PDIR"
export CPLILINE='$'
export ICELINE='$'
export WINDLINE='$'
export CURRLINE='$'
export NFGRIDS=0
export NMGRIDS=1
export WW3GRIDLINE="'glo_1deg'  'no' 'no' 'CPL:native' 'no' 'no' 'no' 'no' 'no' 'no'   1  1  0.00 1.00  F"
export FUNIPNT=T
export IOSRV=1
export FPNTPROC=T
export FGRDPROC=T
export UNIPOINTS='points'
export FLAGMASKCOMP=' F'
export FLAGMASKOUT=' F'
export RUN_BEG="${SYEAR}${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
export RUN_END="2100${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
export OUT_BEG=$RUN_BEG
export OUT_END=$RUN_END
export RST_BEG=$RUN_BEG
export RST_2_BEG=$RUN_BEG
export RST_END=$RUN_END
export RST_2_END=$RUN_END
export WW3_CUR='F'
export WW3_ICE='F'
export WW3_IC1='F'
export WW3_IC5='F'
# ATMW
export WW3_MULTIGRID=true
export WW3_MODDEF=mod_def.glo_1deg
export MESH_WAV=mesh.glo_1deg.nc

# ATMA
export AOD_FRQ=060000

# Regional
export WRITE_RESTART_WITH_BCS=.false.

# Diagnostics
export PRINT_DIFF_PGR=.false.

# Coupling
export coupling_interval_fast_sec=0
export CHOUR=06
export MOM6_OUTPUT_DIR=./MOM6_OUTPUT
}

# Defaults for the CICE6 model namelist, mx100
export_cice6() {
export SECS=`expr $SHOUR \* 3600`
export DT_CICE=${DT_ATMOS}
export CICE_NPT=999
export CICE_RUNTYPE=initial
export CICE_RUNID=unknown
export CICE_USE_RESTART_TIME=.false.
export CICE_RESTART_DIR=./RESTART/
export CICE_RESTART_FILE=iced
export CICE_DUMPFREQ=d
export CICE_DUMPFREQ_N=1000
export CICE_DIAGFREQ=`expr $FHMAX \* 3600 / $DT_CICE`
export CICE_HISTFREQ_N="0, 0, 6, 1, 1"
export CICE_HIST_AVG=.true.
export CICE_HISTORY_DIR=./history/
export CICE_INCOND_DIR=./history/
export CICE_GRID=grid_cice_NEMS_mx${OCNRES}.nc
export CICE_MASK=kmtu_cice_NEMS_mx${OCNRES}.nc
export CICE_GRIDATM=A
export CICE_GRIDOCN=A
export CICE_GRIDICE=B
export CICE_TR_POND_LVL=.true.
export CICE_RESTART_POND_LVL=.false.
# setting to true will allow Frazil FW and Salt to be included in fluxes sent to ocean
export CICE_FRAZIL_FWSALT=.true.
export CICE_KTHERM=2
export CICE_TFREEZE_OPTION=mushy
# SlenderX2
export CICE_NPROC=$ICE_tasks
export np2=`expr $CICE_NPROC / 2`
export CICE_BLCKX=`expr $NX_GLB / $np2`
export CICE_BLCKY=`expr $NY_GLB / 2`
export CICE_DECOMP=slenderX2
}

# Defaults for the MOM6 model namelist, mx100
export_mom6() {
export DT_DYNAM_MOM6=1800
export DT_THERM_MOM6=3600
export MOM6_INPUT=MOM_input_100.IN
export MOM6_OUTPUT_DIR=./MOM6_OUTPUT
export MOM6_RESTART_DIR=./RESTART/
export MOM6_RESTART_SETTING=n
export MOM6_RIVER_RUNOFF=False
export MOM6_FRUNOFF=''
export MOM6_CHLCLIM=seawifs_1998-2006_smoothed_2X.nc
export MOM6_USE_LI2016=True
export MOM6_TOPOEDITS=''
# since CPL_SLOW is set to DT_THERM, this should be always be false
export MOM6_THERMO_SPAN=False
export MOM6_USE_WAVES=True
export MOM6_ALLOW_LANDMASK_CHANGES=False
# MOM6 diag
export MOM6_DIAG_COORD_DEF_Z_FILE=interpolate_zgrid_40L.nc
export MOM6_DIAG_MISVAL='-1e34'
# MOM6 IAU
export ODA_INCUPD=False
export ODA_INCUPD_NHOURS=6
export ODA_TEMPINC_VAR="'pt_inc'"
export ODA_SALTINC_VAR="'s_inc'"
export ODA_THK_VAR="'h_fg'"
export ODA_INCUPD_UV=False
export ODA_UINC_VAR="'u_inc'"
export ODA_VINC_VAR="'v_inc'"
# MOM6 stochastics
export DO_OCN_SPPT=False
export PERT_EPBL=False
export OCN_SPPT=-999.
export EPBL=-999.
}

# Defaults for the WW3 global model
export_ww3() {
export WW3_DOMAIN=mx${OCNRES}
export WW3_MODDEF=mod_def.mx${OCNRES}
export WW3_RSTDTHR=3
export WW3_DT_2_RST="$(printf "%02d" $(( ${WW3_RSTDTHR}*3600 )))"
export WW3_OUTDTHR=3
export WW3_DTFLD="$(printf "%02d" $(( ${WW3_OUTDTHR}*3600 )))"
export WW3_DTPNT="$(printf "%02d" $(( ${WW3_OUTDTHR}*3600 )))"
export WW3_CUR='C'
export WW3_ICE='C'
export WW3_IC1='F'
export WW3_IC5='F'
export WW3_user_sets_restname="true"
}

# Defaults for the coupled 5-component
export_cmeps() {
export UFS_CONFIGURE=ufs.configure.s2swa_fast_esmf.IN
export med_model=cmeps
export atm_model=fv3
export chm_model=gocart
export ocn_model=mom6
export ice_model=cice6
export wav_model=ww3
export lnd_model=noahmp
export coupling_interval_slow_sec=${DT_THERM_MOM6}
export coupling_interval_fast_sec=${DT_ATMOS}
export MESH_OCN=mesh.mx${OCNRES}.nc
export MESH_ICE=mesh.mx${OCNRES}.nc
export MESH_WAV=mesh.${WW3_DOMAIN}.nc
export CPLMODE=ufs.frac
export pio_rearranger=box
export RUNTYPE=startup
export RESTART_N=${FHMAX}
export CMEPS_RESTART_DIR=./RESTART/
export cap_dbug_flag=0
# MOM6 attributes
export use_coldstart=false
export use_mommesh=true
# CICE attributes
export eps_imesh=1.0e-1
# mediator AO flux
export flux_convergence=0.0
export flux_iteration=2
export flux_scheme=0
# mediator ocean albedo
export ocean_albedo_limit=0.06
export use_mean_albedos=.false.
# WW3 (used in run_test only)
export WW3_MULTIGRID=false
}

export_cpl ()
{
export FV3=true
export S2S=true
export HAFS=false
export AQM=false
export DATM_CDEPS=false
export DOCN_CDEPS=false
export CDEPS_INLINE=false
export FV3BMIC='p8c'
export BMIC=.false.
export DAYS=1

#model configure
export MODEL_CONFIGURE=model_configure.IN
export SYEAR=2021
export SMONTH=03
export SDAY=22
export SHOUR=06
export CHOUR=06
export FHMAX=24
export FHROT=0
export DT_ATMOS=720
export QUILTING_RESTART=.false.
export WRTTASK_PER_GROUP=$WPG_cpl_dflt
export WRITE_NSFLIP=.true.
export OUTPUT_FH='6 -1'

# default atm/ocn/ice resolution
export ATMRES=C96
export OCNRES=100
export ICERES=1.00
export NX_GLB=360
export NY_GLB=320
export NPZ=127
export NPZP=128

# default resources
export DOMAINS_STACK_SIZE=8000000
export INPES=$INPES_cpl_dflt
export JNPES=$JNPES_cpl_dflt
export THRD=$THRD_cpl_dflt
OCN_tasks=$OCN_tasks_cpl_dflt
ICE_tasks=$ICE_tasks_cpl_dflt
WAV_tasks=$WAV_tasks_cpl_dflt

# Set CICE6 component defaults
export_cice6

# Set MOM6 component defaults
export_mom6

# Set WW3 component defaults
export_ww3

# Set CMEPS component defauls
export_cmeps
export ATMTILESIZE=`expr $NPX - 1`

# FV3 defaults
export FRAC_GRID=.true.
export CCPP_SUITE=FV3_GFS_v17_coupled_p8
export INPUT_NML=cpld_control.nml.IN
export FIELD_TABLE=field_table_thompson_noaero_tke_GOCART
export DIAG_TABLE=diag_table_cpld.IN
export DIAG_TABLE_ADDITIONAL=''

export FHZERO=6
export DT_INNER=${DT_ATMOS}

# P7 default
export IALB=2
export IEMS=2
export LSM=2
export IOPT_DVEG=4
export IOPT_CRS=2
export IOPT_RAD=3
export IOPT_ALB=1
export IOPT_STC=3
# P8
export IOPT_SFC=3
export IOPT_TRS=2
export IOPT_DIAG=2

# FV3 P7 settings
export D2_BG_K1=0.20
export D2_BG_K2=0.04
#export DZ_MIN=2
export PSM_BC=1
export DDDMP=0.1

#P8
export DZ_MIN=6

# P7 Merra2 Aerosols & NSST
export USE_MERRA2=.true.
export IAER=1011
export NSTF_NAME=2,0,0,0,0

export LHEATSTRG=.false.
export LSEASPRAY=.true.

# P7 UGWP1
export GWD_OPT=2
export KNOB_UGWP_NSLOPE=1
export DO_GSL_DRAG_LS_BL=.true.
export DO_GSL_DRAG_SS=.true.
export DO_UGWP_V1_OROG_ONLY=.false.
export DO_UGWP_V0_NST_ONLY=.false.
export LDIAG_UGWP=.false.
#P8
export DO_GSL_DRAG_TOFD=.false.
export CDMBWD=${CDMBWD_c96}

# P8 RRTMGP
export DO_RRTMGP=.false.
export DOGP_CLDOPTICS_LUT=.true.
export DOGP_LWSCAT=.true.
export DOGP_SGS_CNV=.true.

#P8 UGWD
export DO_UGWP_V0=.true.
export DO_UGWP_V1=.false.
export DO_GSL_DRAG_LS_BL=.false.
export KNOB_UGWP_VERSION=0

# P7 CA
export DO_CA=.true.
export CA_SGS=.true.
export CA_GLOBAL=.false.
export NCA=1
export NCELLS=5
export NLIVES=12
export NTHRESH=18
export NSEED=1
export NFRACSEED=0.5
export CA_TRIGGER=.true.
export NSPINUP=1
export ISEED_CA=12345

# P7 settings
export FNALBC="'C96.snowfree_albedo.tileX.nc'"
export FNALBC2="'C96.facsf.tileX.nc'"
export FNTG3C="'C96.substrate_temperature.tileX.nc'"
export FNVEGC="'C96.vegetation_greenness.tileX.nc'"
export FNVETC="'C96.vegetation_type.tileX.nc'"
export FNSOTC="'C96.soil_type.tileX.nc'"
export FNSOCC="'C96.soil_color.tileX.nc'"
export FNSMCC=${FNSMCC_control}
export FNMSKH=${FNMSKH_control}
export FNVMNC="'C96.vegetation_greenness.tileX.nc'"
export FNVMXC="'C96.vegetation_greenness.tileX.nc'"
export FNSLPC="'C96.slope_type.tileX.nc'"
export FNABSC="'C96.maximum_snow_albedo.tileX.nc'"
export LANDICE=".false."
#P8
export FSICL=0
export FSICS=0

# P8
export USE_CICE_ALB=.true.
export MIN_SEAICE=1.0e-6
export DNATS=2
export IMP_PHYSICS=8
export LGFDLMPRAD=.false.
export DO_SAT_ADJ=.false.
export SATMEDMF=.true.

# P7 default
export CPLFLX=.true.
export CPLICE=.true.
export CPL=.true.
export CPLWAV=.true.
export CPLWAV2ATM=.true.
export USE_MED_FLUX=.false.
export CPLCHM=.true.
export CPLLND=.false.

# for FV3: default values will be changed if doing a warm-warm restart
export WARM_START=.false.
export MAKE_NH=.true.
export NA_INIT=1
export EXTERNAL_IC=.true.
export NGGPS_IC=.true.
export MOUNTAIN=.false.
# gocart inst_aod output; uses AERO_HIST.rc.IN from parm/gocart directory
export AOD_FRQ=060000

# checkpoint restarts
export RESTART_FILE_PREFIX=''
export RESTART_FILE_SUFFIX_SECS=''
export RT35D=''
}
export_35d_run ()
{
export CNTL_DIR=""
export LIST_FILES=""
}
export_datm_cdeps ()
{
export FV3=false
export S2S=false
export HAFS=false
export AQM=false
export DATM_CDEPS=true
export DOCN_CDEPS=false
export CDEPS_INLINE=false
export DAYS=1

# model configure
export MODEL_CONFIGURE=datm_cdeps_configure.IN
export SYEAR=2011
export SMONTH=10
export SDAY=01
export SHOUR=00
export FHMAX=24
export DT_ATMOS=900
export FHROT=0

# required but unused
export WARM_START=.false.
export CPLWAV=.false.
export CPLCHM=.false.

# atm/ocn/ice resolution
export IATM=1760
export JATM=880
export ATM_NX_GLB=$IATM
export ATM_NY_GLB=$JATM
export ATMRES=${IATM}x${JATM}
export OCNRES=100
export ICERES=1.00
export NX_GLB=360
export NY_GLB=320

# default resources
export ATM_compute_tasks=$ATM_compute_tasks_cdeps_100
export OCN_tasks=$OCN_tasks_cdeps_100
export ICE_tasks=$ICE_tasks_cdeps_100

# Set CICE6 component defaults
export_cice6
# default non-mushy thermo for CICE
export CICE_KTHERM=1
export CICE_TFREEZE_OPTION=linear_salt

# Set MOM6 component defaults
export_mom6
# default no waves
export MOM6_USE_LI2016=False
export MOM6_USE_WAVES=False
export WW3_DOMAIN=''

# Set CMEPS component defauls
export_cmeps
# default configure
export UFS_CONFIGURE=ufs.configure.datm_cdeps.IN
export atm_model=datm
export CPLMODE=ufs.nfrac.aoflux

# datm defaults
export INPUT_NML=input.mom6.nml.IN
export DIAG_TABLE=diag_table_template
export DATM_SRC=CFSR
export FILENAME_BASE=cfsr.
export MESH_ATM=${FILENAME_BASE//.}_mesh.nc
export atm_datamode=${DATM_SRC}
export stream_files=INPUT/${FILENAME_BASE}201110.nc
export EXPORT_ALL=.false.
export STREAM_OFFSET=0

export BL_SUFFIX=""
export RT_SUFFIX=""
}
export_hafs_datm_cdeps ()
{
export FV3=false
export S2S=false
export HAFS=true
export AQM=false
export DATM_CDEPS=true
export DOCN_CDEPS=false
export CDEPS_INLINE=false
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export NTILES=1

export atm_model=datm
export DATM_IN_CONFIGURE=datm_in.IN
export DATM_STREAM_CONFIGURE=hafs_datm.streams.era5.IN
export EXPORT_ALL=.false.
}
export_hafs_docn_cdeps ()
{
export FV3=true
export S2S=false
export HAFS=true
export AQM=false
export DOCN_CDEPS=true
export CDEPS_INLINE=false
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export NTILES=1

export ocn_model=docn
export ocn_datamode=sstdata
export pio_rearranger=box
export DOCN_IN_CONFIGURE=docn_in.IN
export DOCN_STREAM_CONFIGURE=hafs_docn.streams.IN
}
export_hafs_regional ()
{
export FV3=true
export S2S=false
export HAFS=true
export AQM=false
export DATM_CDEPS=false
export DOCN_CDEPS=false
export CDEPS_INLINE=false
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export NTILES=1

# model_configure
export SYEAR=2019
export SMONTH=08
export SDAY=29
export SHOUR=00
export SECS=`expr $SHOUR \* 3600`
export FHMAX=6
export ENS_NUM=1
export DT_ATMOS=900
export CPL=.true.
export RESTART_INTERVAL=0
export FHROT=0
export coupling_interval_fast_sec=0
export QUILTING=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'atm' 'sfc'"
export OUTPUT_GRID="'regional_latlon'"
export OUTPUT_FILE="'netcdf'"
export IDEFLATE=0
export QUANTIZE_NSD=0
export CEN_LON=-62.0
export CEN_LAT=25.0
export LON1=-114.5
export LAT1=-5.0
export LON2=-9.5
export LAT2=55.0
export DLON=0.03
export DLAT=0.03

# shel.inp
# input.nml
export CPL_IMP_MRG=.true.
export DIAG_TABLE=diag_table_hafs
export FIELD_TABLE=field_table_hafs

export OCNRES=''
export ICERES=''
export DT_THERM_MOM6=''

# Set WW3 component defaults
export_ww3
# default hafs with no ice
export WW3_DOMAIN=natl_6m
export WW3_MODDEF=mod_def.${WW3_DOMAIN}
export WW3_ICE='F'
export WW3_OUTPARS="WND HS T01 T02 DIR FP DP PHS PTP PDIR UST CHA USP"

# Set CMEPS component defaults
export_cmeps
# default hafs
export ocn_model=hycom
export CPLMODE=hafs
export MESH_WAV=mesh.hafs.nc
}

export_hafs ()
{
export FV3=true
export S2S=false
export HAFS=true
export AQM=false
export DATM_CDEPS=false
export DOCN_CDEPS=false
export CDEPS_INLINE=false
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export NTILES=1
export IMFSHALCNV=2
export IMFDEEPCNV=2
export HYBEDMF=.false.
export SATMEDMF=.true.
export MONINQ_FAC=-1.0
export HURR_PBL=.true.
export ISATMEDMF=1
export IOPT_SFC=1
export IOPT_DVEG=2
export IOPT_CRS=1
export IOPT_RAD=1
export IOPT_ALB=2
export IOPT_STC=1
export LSM=1
export DO_GSL_DRAG_LS_BL=.true.
export DO_GSL_DRAG_SS=.true.
export DO_GSL_DRAG_TOFD=.true.
export IMP_PHYSICS=11
export IAER=111
export CNVGWD=.false.
export LTAEROSOL=.false.
export CDMBWD=1.0,1.0,1.0,1.0
export LHEATSTRG=.false.
export LRADAR=.true.

export FV_CORE_TAU=5.
export RF_CUTOFF=30.e2
export RF_CUTOFF_NEST=50.e2

export IS_MOVING_NEST=".false."
export VORTEX_TRACKER=0
export NTRACK=0
export MOVE_CD_X=0
export MOVE_CD_Y=0
export CPL_IMP_MRG=.true.

export OUTPUT_GRID=''
export IMO=''
export JMO=''
export CEN_LON=''
export CEN_LAT=''
export LON1=''
export LAT1=''
export LON2=''
export LAT2=''
export DLON=''
export DLAT=''
export STDLAT1=''
export STDLAT2=''
export NX=''
export NY=''
export DX=''
export DY=''

export OUTPUT_GRID_2=''
export IMO_2=''
export JMO_2=''
export CEN_LON_2=''
export CEN_LAT_2=''
export LON1_2=''
export LAT1_2=''
export LON2_2=''
export LAT2_2=''
export DLON_2=''
export DLAT_2=''
export STDLAT1_2=''
export STDLAT2_2=''
export NX_2=''
export NY_2=''
export DX_2=''
export DY_2=''

export OUTPUT_GRID_3=''
export IMO_3=''
export JMO_3=''
export CEN_LON_3=''
export CEN_LAT_3=''
export LON1_3=''
export LAT1_3=''
export LON2_3=''
export LAT2_3=''
export DLON_3=''
export DLAT_3=''
export STDLAT1_3=''
export STDLAT2_3=''
export NX_3=''
export NY_3=''
export DX_3=''
export DY_3=''

export OUTPUT_GRID_4=''
export IMO_4=''
export JMO_4=''
export CEN_LON_4=''
export CEN_LAT_4=''
export LON1_4=''
export LAT1_4=''
export LON2_4=''
export LAT2_4=''
export DLON_4=''
export DLAT_4=''
export STDLAT1_4=''
export STDLAT2_4=''
export NX_4=''
export NY_4=''
export DX_4=''
export DY_4=''

export OUTPUT_GRID_5=''
export IMO_5=''
export JMO_5=''
export CEN_LON_5=''
export CEN_LAT_5=''
export LON1_5=''
export LAT1_5=''
export LON2_5=''
export LAT2_5=''
export DLON_5=''
export DLAT_5=''
export STDLAT1_5=''
export STDLAT2_5=''
export NX_5=''
export NY_5=''
export DX_5=''
export DY_5=''

export OUTPUT_GRID_6=''
export IMO_6=''
export JMO_6=''
export CEN_LON_6=''
export CEN_LAT_6=''
export LON1_6=''
export LAT1_6=''
export LON2_6=''
export LAT2_6=''
export DLON_6=''
export DLAT_6=''
export STDLAT1_6=''
export STDLAT2_6=''
export NX_6=''
export NY_6=''
export DX_6=''
export DY_6=''

export OUTPUT_FH='3 -1'
}
export_hrrr() {
export_fv3
export NPZ=127
export NPZP=128
export DT_ATMOS=300
export SYEAR=2021
export SMONTH=03
export SDAY=22
export SHOUR=06
export OUTPUT_GRID='gaussian_grid'
export NSTF_NAME='2,0,0,0,0'
export WRITE_DOPOST=.true.
export IAER=5111
export FHMAX=12

export FRAC_GRID=.false.
export FRAC_ICE=.true.

export FV_CORE_TAU=10.
export RF_CUTOFF=7.5e2

export FV3_RUN=lake_control_run.IN
export CCPP_SUITE=FV3_HRRR
export INPUT_NML=rap.nml.IN
export FIELD_TABLE=field_table_thompson_aero_tke
export NEW_DIAGTABLE=diag_table_rap

export SFCLAY_COMPUTE_FLUX=.true.

export LKM=1
export IOPT_LAKE=2
export IMP_PHYSICS=8
export DNATS=0
export DO_SAT_ADJ=.false.
export LRADAR=.true.
export LTAEROSOL=.true.
export IALB=2
export IEMS=2
export HYBEDMF=.false.
export DO_MYNNEDMF=.true.
export DO_MYNNSFCLAY=.true.
export DO_DEEP=.false.
export SHAL_CNV=.false.
export IMFSHALCNV=-1
export IMFDEEPCNV=-1
export LHEATSTRG=.false.
export LSM=3
export LSOIL_LSM=9
export KICE=9

export GWD_OPT=3
export DO_UGWP_V0=.false.
export DO_UGWP_V0_OROG_ONLY=.false.
export DO_GSL_DRAG_LS_BL=.true.
export DO_GSL_DRAG_SS=.true.
export DO_GSL_DRAG_TOFD=.true.
export DO_UGWP_V1=.false.
export DO_UGWP_V1_OROG_ONLY=.false.
}
export_hrrr_conus13km()
{
export_fv3
export SYEAR=2021
export SMONTH=05
export SDAY=12
export SHOUR=16
export FHMAX=2
export DT_ATMOS=120
export RESTART_INTERVAL=1
export QUILTING=.true.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export NTILES=1
export WRITE_DOPOST=.false.
export OUTPUT_HISTORY=.true.
export OUTPUT_GRID=lambert_conformal
export OUTPUT_FILE="'netcdf'"

# Revert these two to GFS_typedefs defaults to avoid a crash:
export SEDI_SEMI=.false.
export DECFL=8

export RRFS_SMOKE=.true.
export SEAS_OPT=0

export LKM=1
export SFCLAY_COMPUTE_FLUX=.true.
export IALB=2
export ICLIQ_SW=2
export IEMS=2
export IOVR=3
export KICE=9
export LSM=3
export LSOIL_LSM=9
export DO_MYNNSFCLAY=.true.
export DO_MYNNEDMF=.true.
export HYBEDMF=.false.
export SHAL_CNV=.false.
export DO_SAT_ADJ=.false.
export DO_DEEP=.false.
export CCPP_SUITE='FV3_HRRR'
export INPES=12
export JNPES=12
export NPX=397
export NPY=233
export NPZ=65
export MAKE_NH=.false.
export NA_INIT=0
export DNATS=0
export EXTERNAL_IC=.false.
export NGGPS_IC=.false.
export MOUNTAIN=.true.
export WARM_START=.true.
export READ_INCREMENT=.false.
export RES_LATLON_DYNAMICS="'fv3_increment.nc'"
export NPZP=66
export FHZERO=1.0
export IMP_PHYSICS=8
export LDIAG3D=.false.
export QDIAG3D=.false.
export PRINT_DIFF_PGR=.true.
export FHCYC=0.0
export IAER=1011
export LHEATSTRG=.false.
export RANDOM_CLDS=.false.
export CNVCLD=.false.
export IMFSHALCNV=-1
export IMFDEEPCNV=-1
export CDMBWD='3.5,1.0'
export DO_SPPT=.false.
export DO_SHUM=.false.
export DO_SKEB=.false.
export LNDP_TYPE=0
export N_VAR_LNDP=0

export GWD_OPT=3
export DO_UGWP_V0=.false.
export DO_UGWP_V0_OROG_ONLY=.false.
export DO_GSL_DRAG_LS_BL=.true.
export DO_GSL_DRAG_SS=.true.
export DO_GSL_DRAG_TOFD=.true.
export DO_UGWP_V1=.false.
export DO_UGWP_V1_OROG_ONLY=.false.

export FV3_RUN=rrfs_warm_run.IN
export INPUT_NML=rrfs_conus13km_hrrr.nml.IN
export FIELD_TABLE=field_table_thompson_aero_tke_smoke
export DIAG_TABLE=diag_table_hrrr
export MODEL_CONFIGURE=model_configure_rrfs_conus13km.IN
export DIAG_TABLE_ADDITIONAL=diag_additional_rrfs_smoke
export FRAC_ICE=.true.
}
export_rap_common()
{
export_fv3
export NPZ=127
export NPZP=128
export DT_ATMOS=300
export SYEAR=2021
export SMONTH=03
export SDAY=22
export SHOUR=06
export OUTPUT_GRID='gaussian_grid'
export NSTF_NAME='2,0,0,0,0'
export WRITE_DOPOST=.true.
export IAER=5111

export FV_CORE_TAU=10.
export RF_CUTOFF=7.5e2

export FV3_RUN=control_run.IN
export INPUT_NML=rap.nml.IN
export FIELD_TABLE=field_table_thompson_aero_tke

export LHEATSTRG=.false.
export IMP_PHYSICS=8
export DNATS=0
export DO_SAT_ADJ=.false.
export LRADAR=.true.
export LTAEROSOL=.true.
export IALB=2
export IEMS=2
export HYBEDMF=.false.
export DO_MYNNEDMF=.true.
export DO_MYNNSFCLAY=.true.
}
export_rap()
{
export_rap_common

export DIAG_TABLE=diag_table_rap
export CCPP_SUITE=FV3_RAP

export IMFSHALCNV=3
export IMFDEEPCNV=3
export LSM=3
export LSOIL_LSM=9
export KICE=9

export GWD_OPT=3
export DO_UGWP_V0=.false.
export DO_UGWP_V0_OROG_ONLY=.false.
export DO_GSL_DRAG_LS_BL=.true.
export DO_GSL_DRAG_SS=.true.
export DO_GSL_DRAG_TOFD=.true.
export DO_UGWP_V1=.false.
export DO_UGWP_V1_OROG_ONLY=.false.
}
export_rrfs_v1()
{
export_rap_common

export CCPP_SUITE=FV3_RRFS_v1beta
export DIAG_TABLE=diag_table_rap_noah

export DO_DEEP=.false.
export SHAL_CNV=.false.
export IMFSHALCNV=-1
export IMFDEEPCNV=-1
export LHEATSTRG=.false.
export LSM=2
export LSOIL_LSM=4
}
