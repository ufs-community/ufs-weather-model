
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

  TPN=24

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=6 ; JNPES_c384=8  ; THRD_c384=1
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

  THRD_cpl_atmw_gdas=3
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=264

elif [[ $MACHINE_ID = cheyenne ]]; then

  TPN=36

  INPES_dflt=3 ; JNPES_dflt=8
  INPES_thrd=3 ; JNPES_thrd=4
  INPES_c384=8 ; JNPES_c384=6  ; THRD_c384=2
  INPES_c768=8 ; JNPES_c768=16 ; THRD_c768=2

  THRD_cpl_atmw_gdas=3
  INPES_cpl_atmw_gdas=6; JNPES_cpl_atmw_gdas=8; WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=264

  THRD_cpl_dflt=1
  INPES_cpl_dflt=6; JNPES_cpl_dflt=8; WPG_cpl_dflt=6
  OCN_tasks_cpl_dflt=20
  ICE_tasks_cpl_dflt=10
  WAV_tasks_cpl_dflt=20

  THRD_cpl_thrd=2
  INPES_cpl_thrd=6; JNPES_cpl_thrd=4; WPG_cpl_thrd=6
  OCN_tasks_cpl_thrd=20
  ICE_tasks_cpl_thrd=10
  WAV_tasks_cpl_thrd=20

  THRD_cpl_dcmp=1
  INPES_cpl_dcmp=4; JNPES_cpl_dcmp=12; WPG_cpl_dcmp=6
  OCN_tasks_cpl_dcmp=20
  ICE_tasks_cpl_dcmp=10
  WAV_tasks_cpl_dcmp=20

  THRD_cpl_mpi=1
  INPES_cpl_mpi=6; JNPES_cpl_mpi=12; WPG_cpl_mpi=6
  OCN_tasks_cpl_mpi=24
  ICE_tasks_cpl_mpi=12
  WAV_tasks_cpl_mpi=24

  THRD_cpl_bmrk=2
  INPES_cpl_bmrk=12; JNPES_cpl_bmrk=8; WPG_cpl_bmrk=48
  OCN_tasks_cpl_bmrk=120
  ICE_tasks_cpl_bmrk=48
  WAV_tasks_cpl_bmrk=100

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

export_fv3 ()
{
# nems.configure defaults
export NEMS_CONFIGURE=nems.configure.atm.IN
export MODEL_CONFIGURE=model_configure.IN
export atm_model=fv3

export FV3=true
export S2S=false
export HAFS=false
export AQM=false
export DATM_CDEPS=false
export DOCN_CDEPS=false
export POSTAPP='global'
export USE_MERRA2=.false.

export NTILES=6
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export RESTART_INTERVAL=0
export QUILTING=.true.
export QUILTING_RESTART=.false.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export ITASKS=1
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'atm' 'sfc'"
export OUTPUT_GRID="'cubed_sphere_grid'"
export OUTPUT_FILE="'netcdf'"
export IDEFLATE=0
export NBITS=0
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
export DO_UGWP_V0_OROG_ONLY=.false.
export DO_GSL_DRAG_LS_BL=.false.
export DO_GSL_DRAG_SS=.false.
export DO_GSL_DRAG_TOFD=.false.
export DO_UGWP_V1=.false.
export DO_UGWP_V1_OROG_ONLY=.false.

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
export DO_MYJPBL=.false.
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

# Aerosol convective scavenging
export FSCAV_AERO='"*:0.3","so2:0.0","msa:0.0","dms:0.0","nh3:0.4","nh4:0.6","bc1:0.6","bc2:0.6","oc1:0.4","oc2:0.4","dust1:0.6","dust2:0.6","dust3:0.6","dust4:0.6","dust5:0.6","seas1:0.5","seas2:0.5","seas3:0.5","seas4:0.5","seas5:0.5"'

# SFC
export DO_MYJSFC=.false.
export DO_MYNNSFCLAY=.false.

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
export FRAC_ICE=.false.

export CPL=.false.
export CPLCHM=.false.
export CPLFLX=.false.
export CPLICE=.false.
export CPLWAV=.false.
export CPLWAV2ATM=.false.
export CPLLND=.false.
export USE_MED_FLUX=.false.
export DAYS=1
export NPX=97
export NPY=97
export NPZ=64
export NPZP=65
export NSTF_NAME=2,1,1,0,5
export OUTPUT_FH="12 -1"
export NFHOUT=12
export NFHMAX_HF=12
export NFHOUT_HF=6
export IAU_OFFSET=0
export FHZERO=6
export FNALBC="'global_snowfree_albedo.bosu.t126.384.190.rg.grb'"
export FNVETC="'global_vegtype.igbp.t126.384.190.rg.grb'"
export FNSOTC="'global_soiltype.statsgo.t126.384.190.rg.grb'"
export FNSMCC="'global_soilmgldas.t126.384.190.grb'"
export FNSMCC_control="'global_soilmgldas.statsgo.t1534.3072.1536.grb'"
export FNMSKH_control="'global_slmask.t1534.3072.1536.grb'"
export FNABSC="'global_mxsnoalb.uariz.t126.384.190.rg.grb'"

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

export FH_DFI_RADAR='-2e10'

#Cellular automata
export DO_CA=.false.
export CA_SGS=.false.
export CA_GLOBAL=.false.

export IAU_DRYMASSFIXER=.false.

#waves
export WW3RSTDTHR=12
export DT_2_RST="$(printf "%02d" $(( ${WW3RSTDTHR}*3600 )))"
export WW3OUTDTHR=1
export DTFLD="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
export DTPNT="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
export DTRST=0
export RSTTYPE=T
export GOFILETYPE=1
export POFILETYPE=1
export OUTPARS_WAV="WND HS FP DP PHS PTP PDIR"
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
export WAV_CUR='F'
export WAV_ICE='F'
export WAV_IC1='F'
export WAV_IC5='F'
# ATMW
export MULTIGRID=true
export MODDEF_WAV=mod_def.glo_1deg
export MESH_WAV=mesh.glo_1deg.nc

# ATMA
export AOD_FRQ=060000

# Regional
export WRITE_RESTART_WITH_BCS=.false.

# Diagnostics
export PRINT_DIFF_PGR=.false.

# Coupling
export coupling_interval_fast_sec=0
}

export_cpl ()
{
export FV3=true
export S2S=true
export HAFS=false
export AQM=false
export DATM_CDEPS=false
export DOCN_CDEPS=false

export FV3BMIC='p8c'

export SYEAR=2021
export SMONTH=03
export SDAY=22
export SHOUR=06
export SECS=`expr $SHOUR \* 3600`
export BMIC=.false.

export DAYS=1
export FHMAX=24
export FDIAG=6
export FHZERO=6

export WRITE_NSFLIP=.true.

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
export WRTTASK_PER_GROUP=$WPG_cpl_dflt

OCN_tasks=$OCN_tasks_cpl_dflt
ICE_tasks=$ICE_tasks_cpl_dflt
WAV_tasks=$WAV_tasks_cpl_dflt

# component and coupling timesteps
export DT_ATMOS=720
export DT_INNER=${DT_ATMOS}
export DT_CICE=${DT_ATMOS}
export DT_DYNAM_MOM6=1800
export DT_THERM_MOM6=3600

# nems.configure defaults
export NEMS_CONFIGURE=nems.configure.cpld.IN
export med_model=cmeps
export atm_model=fv3
export chm_model=gocart
export ocn_model=mom6
export ice_model=cice6
export wav_model=ww3

export coupling_interval_slow_sec=${DT_THERM_MOM6}
export coupling_interval_fast_sec=${DT_ATMOS}

export RESTART_N=${FHMAX}
export CPLMODE=nems_frac
export cap_dbug_flag=0
export use_coldstart=false
export use_mommesh=true
export RUNTYPE=startup
export CICERUNTYPE=initial
export eps_imesh=1.0e-1
export ATMTILESIZE=`expr $NPX - 1`
export MULTIGRID=false

# FV3 defaults
export FRAC_GRID=.true.
export CCPP_SUITE=FV3_GFS_v17_coupled_p8
export INPUT_NML=cpld_control.nml.IN
export FIELD_TABLE=field_table_thompson_noaero_tke_GOCART
export DIAG_TABLE=diag_table_p8_template

export DIAG_TABLE_ADDITIONAL=''

export FHROT=0
export NSOUT=-1
export OUTPUT_FH='6 -1'

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

# P7 default mushy thermo
export KTHERM=2
export TFREEZE_OPTION=mushy
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

# MOM6 defaults; 1 degree
export MOM_INPUT=MOM_input_template_100
export MOM6_RESTART_SETTING=n
export MOM6_RIVER_RUNOFF=False
export FRUNOFF=''
export CHLCLIM=seawifs_1998-2006_smoothed_2X.nc
export MOM6_USE_LI2016=True
# since CPL_SLOW is set to DT_THERM, this should be always be false
export MOM6_THERMO_SPAN=False
export MOM6_USE_WAVES=True
export MOM6_ALLOW_LANDMASK_CHANGES=False
# MOM6 IAU
export MOM_IAU=False
export MOM_IAU_HRS=6
# MOM6 stochastics
export DO_OCN_SPPT=False
export PERT_EPBL=False
export OCN_SPPT=-999.
export EPBL=-999.

# CICE6 defaults; 1 degree
export CICE_DECOMP=slenderX2
export NPROC_ICE=$ICE_tasks
# SlenderX2
export CICE_DECOMP=slenderX2
export np2=`expr $NPROC_ICE / 2`
export BLCKX=`expr $NX_GLB / $np2`
export BLCKY=`expr $NY_GLB / 2`
export MESHOCN_ICE=mesh.mx${OCNRES}.nc
export WAVDOMAIN=mx${OCNRES}
export MODDEF_WAV=mod_def.mx${OCNRES}
export MESH_WAV=mesh.${WAVDOMAIN}.nc
export CICEGRID=grid_cice_NEMS_mx${OCNRES}.nc
export CICEMASK=kmtu_cice_NEMS_mx${OCNRES}.nc
export RUNID=unknown
# set large; restart frequency now controlled by restart_n in nems.configure
export DUMPFREQ=d
export DUMPFREQ_N=1000
export USE_RESTART_TIME=.false.
export RESTART_EXT=.false.
# setting to true will allow Frazil FW and Salt to be
# included in fluxes sent to ocean
export FRAZIL_FWSALT=.true.
# default to write CICE average history files
export CICE_HIST_AVG=.true.
# default CICE B-grid, ATM and OCN are provided by cap on A grid
export GRIDATM=A
export GRIDOCN=A
export GRIDICE=B

#wave
export WW3RSTDTHR=3
export DT_2_RST="$(printf "%02d" $(( ${WW3RSTDTHR}*3600 )))"
export WW3OUTDTHR=3
export DTFLD="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
export DTPNT="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
# waves when using shel.nml.IN
export WAV_CUR='C'
export WAV_ICE='C'
export WAV_IC1='F'
export WAV_IC5='F'
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
export CPLWAV=.false.
export CPLCHM=.false.
export DAYS=1
export FHMAX=24
export FHROT=0
export WARM_START=.false.
export WRITE_NSFLIP=.true.

# atm/ocn/ice resolution
export IATM=1760
export JATM=880
export ATM_NX_GLB=$IATM
export ATM_NY_GLB=$JATM
export ATMRES=1760x880
export OCNRES=100
export ICERES=1.00
export NX_GLB=360
export NY_GLB=320

# nems.configure
export NEMS_CONFIGURE=nems.configure.datm_cdeps.IN
export med_model=cmeps
export atm_model=datm
export ocn_model=mom6
export ice_model=cice6
export ATM_compute_tasks=$ATM_compute_tasks_cdeps_100
export OCN_tasks=$OCN_tasks_cdeps_100
export ICE_tasks=$ICE_tasks_cdeps_100
# SlenderX2
export CICE_DECOMP=slenderX2
export NPROC_ICE=12
export np2=`expr $NPROC_ICE / 2`
export BLCKX=`expr $NX_GLB / $np2`
export BLCKY=`expr $NY_GLB / 2`

export ENS_NUM=1
export SYEAR=2011
export SMONTH=10
export SDAY=01
export SHOUR=00
export SECS=`expr $SHOUR \* 3600`
export CDATE=${SYEAR}${SMONTH}${SDAY}${SHOUR}

export NFHOUT=6
export DT_ATMOS=900
export DT_CICE=${DT_ATMOS}
export DT_DYNAM_MOM6=1800
export DT_THERM_MOM6=3600
export coupling_interval_slow_sec=${DT_THERM_MOM6}
export coupling_interval_fast_sec=${DT_ATMOS}

export RESTART_N=${FHMAX}
export CPLMODE=nems_orig_data
export cap_dbug_flag=0
export use_coldstart=false
export use_mommesh=true
export RUNTYPE=startup
export CICERUNTYPE=initial
export eps_imesh=1.0e-1
export flux_convergence=0.0
export flux_iteration=2
export flux_scheme=0

export INPUT_NML=input.mom6.nml.IN
export MODEL_CONFIGURE=datm_cdeps_configure.IN
export DIAG_TABLE=diag_table_template

# atm defaults
export DATM_SRC=CFSR
export FILENAME_BASE=cfsr.
export mesh_file=cfsr_mesh.nc
export MESH_ATM=DATM_INPUT/${mesh_file}
export atm_datamode=${DATM_SRC}
export stream_files=DATM_INPUT/${FILENAME_BASE}201110.nc

# MOM6 defaults; 1 degree
export MOM_INPUT=MOM_input_template_100
export MOM6_RESTART_SETTING=n
export MOM6_RIVER_RUNOFF=False
export FRUNOFF=''
export CHLCLIM=seawifs_1998-2006_smoothed_2X.nc
# MOM6 IAU
export MOM_IAU=False
export MOM_IAU_HRS=6
export MOM6_USE_LI2016=False
# MOM6 stochastics
export DO_OCN_SPPT=False
export PERT_EPBL=False
export OCN_SPPT=-999.
export EPBL=-999.
# since coupling_interval_slow is set to DT_THERM, this should be always be false
export MOM6_THERMO_SPAN=False
export MOM6_USE_WAVES=False
export MOM6_ALLOW_LANDMASK_CHANGES=False

# CICE6 defaults; 1 degree
export MESHOCN_ICE=mesh.mx${OCNRES}.nc
export CICEGRID=grid_cice_NEMS_mx${OCNRES}.nc
export CICEMASK=kmtu_cice_NEMS_mx${OCNRES}.nc
export RUNID=unknown
# set large; restart frequency now controlled by restart_n in nems.configure
export DUMPFREQ=d
export DUMPFREQ_N=1000
export USE_RESTART_TIME=.false.
export RESTART_EXT=.false.
# setting to true will allow Frazil FW and Salt to be
# included in fluxes sent to ocean
export FRAZIL_FWSALT=.true.
# default to write CICE average history files
export CICE_HIST_AVG=.true.
export GRIDATM=A
export GRIDOCN=A
export GRIDICE=B
# default non-mushy thermo
export KTHERM=1
export TFREEZE_OPTION=linear_salt
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
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export NTILES=1

export atm_model=datm

export DATM_IN_CONFIGURE=datm_in
export DATM_STREAM_CONFIGURE=hafs_datm.streams.era5.IN
}
export_hafs_docn_cdeps ()
{
export FV3=true
export S2S=false
export HAFS=true
export AQM=false
export DOCN_CDEPS=true
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export NTILES=1

export ocn_model=docn
export ocn_datamode=sstdata

export DOCN_IN_CONFIGURE=docn_in
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
export QUILTING_RESTART=.false.
export WRITE_GROUP=1
export WRTTASK_PER_GROUP=6
export OUTPUT_HISTORY=.true.
export WRITE_DOPOST=.false.
export NUM_FILES=2
export FILENAME_BASE="'atm' 'sfc'"
export OUTPUT_GRID="'regional_latlon'"
export OUTPUT_FILE="'netcdf'"
export IDEFLATE=0
export NBITS=0
export NFHOUT=3
export NFHMAX_HF=-1
export NFHOUT_HF=3
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
export WW3RSTDTHR=${FHMAX}
export DT_2_RST="$(printf "%02d" $(( ${WW3RSTDTHR}*3600 )))"
export WW3OUTDTHR=3
export DTFLD="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
export DTPNT="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
export OUTPARS_WAV="WND HS T01 T02 DIR FP DP PHS PTP PDIR UST CHA USP"
export WAV_CUR='C'

# nems.configure
export med_model=cmeps
export CAP_DBUG_FLAG=0
export RESTART_N=${FHMAX}
export CPLMODE=hafs
export RUNTYPE=startup
export USE_COLDSTART=false
export MESH_WAV=mesh.hafs.nc
export MODDEF_WAV=mod_def.natl_6m
export MULTIGRID=false
}

export_hafs ()
{
export FV3=true
export S2S=false
export HAFS=true
export AQM=false
export DATM_CDEPS=false
export DOCN_CDEPS=false
export INPES=$INPES_dflt
export JNPES=$JNPES_dflt
export NTILES=1

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

export NFHOUT=3
export NFHMAX_HF=-1
export NFHOUT_HF=3
export NSOUT=-1
export OUTPUT_FH=-1
}
