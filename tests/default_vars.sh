#!/bin/bash
# shellcheck disable=SC2034
###############################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
###############################################################################

THRD=1

INPES_atmaero=4
JNPES_atmaero=8
WPG_atmaero=6

THRD_cpl_atmw=1
INPES_cpl_atmw=3
JNPES_cpl_atmw=8
WPG_cpl_atmw=6
WAV_tasks_cpl_atmw=30
WAV_thrds_cpl_atmw=1

THRD_cpl_c48=1
INPES_cpl_c48=1
JNPES_cpl_c48=1
WPG_cpl_c48=6
OCN_tasks_cpl_c48=4
ICE_tasks_cpl_c48=4
WAV_tasks_cpl_c48=4

THRD_cpl_dflt=1
INPES_cpl_dflt=3
JNPES_cpl_dflt=8
WPG_cpl_dflt=6
OCN_tasks_cpl_dflt=20
ICE_tasks_cpl_dflt=10
WAV_tasks_cpl_dflt=20

THRD_cpl_thrd=2
INPES_cpl_thrd=3
JNPES_cpl_thrd=4
WPG_cpl_thrd=6
OCN_tasks_cpl_thrd=20
OCN_thrds_cpl_thrd=1
ICE_tasks_cpl_thrd=10
ICE_thrds_cpl_thrd=1
WAV_tasks_cpl_thrd=12
WAV_thrds_cpl_thrd=2

THRD_cpl_dcmp=1
INPES_cpl_dcmp=4
JNPES_cpl_dcmp=6
WPG_cpl_dcmp=6
OCN_tasks_cpl_dcmp=20
ICE_tasks_cpl_dcmp=10
WAV_tasks_cpl_dcmp=20

THRD_cpl_mpi=1
INPES_cpl_mpi=4
JNPES_cpl_mpi=8
WPG_cpl_mpi=6
OCN_tasks_cpl_mpi=34
ICE_tasks_cpl_mpi=20
WAV_tasks_cpl_mpi=28

THRD_cpl_bmrk=2
INPES_cpl_bmrk=8
JNPES_cpl_bmrk=8
WPG_cpl_bmrk=48
OCN_tasks_cpl_bmrk=120
OCN_thrds_cpl_bmrk=1
ICE_tasks_cpl_bmrk=48
ICE_thrds_cpl_bmrk=1
WAV_tasks_cpl_bmrk=120
WAV_thrds_cpl_bmrk=2

THRD_cpl_c192=2
INPES_cpl_c192=6
JNPES_cpl_c192=8
WPG_cpl_c192=12
OCN_tasks_cpl_c192=60
ICE_tasks_cpl_c192=24
WAV_tasks_cpl_c192=80

ATM_compute_tasks_cdeps_100=12
OCN_tasks_cdeps_100=16
ICE_tasks_cdeps_100=12

ATM_compute_tasks_cdeps_025=40
OCN_tasks_cdeps_025=120
ICE_tasks_cdeps_025=48

INPES_aqm=33
JNPES_aqm=8

INPES_sfs=4
JNPES_sfs=6
THRD_sfs=1
WPG_sfs=24
OCN_tasks_sfs=168
ICE_tasks_sfs=48

THRD_cpl_unstr=1
INPES_cpl_unstr=3
JNPES_cpl_unstr=8
WPG_cpl_unstr=6
OCN_tasks_cpl_unstr=20
ICE_tasks_cpl_unstr=10
WAV_tasks_cpl_unstr=60

THRD_cpl_unstr_mpi=1
INPES_cpl_unstr_mpi=4
JNPES_cpl_unstr_mpi=8
WPG_cpl_unstr_mpi=6
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
fbh_omp_num_threads=1

histaux_enabled=.false.
BMIC=.false.

GFSv17opn=.false.
SFS=.false.

EXCLUSIVE_NODES=.false.
MEM_PER_NODE_OPT=""

# Enable CATChem coupling when true. The CPLCAT collection is a superset of CPLCHM. CPLCHM should
# always be true if CPLCAT is true.
CPLCAT=.false.

if [[ ${MACHINE_ID} = wcoss2 || ${MACHINE_ID} = acorn ]]; then

  TPN=128
  EXCLUSIVE_NODES=.true.

  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=8
  JNPES_c384=6
  THRD_c384=2
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ ${MACHINE_ID} = orion ]]; then

  TPN=40
  EXCLUSIVE_NODES=.true.

  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=8
  JNPES_c384=6
  THRD_c384=2
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ ${MACHINE_ID} = hercules ]]; then

  TPN=80
  EXCLUSIVE_NODES=.true.
  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=8
  JNPES_c384=6
  THRD_c384=2
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ ${MACHINE_ID} = hera ]]; then

  TPN=40

  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=6
  JNPES_c384=8
  THRD_c384=2
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=4

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ ${MACHINE_ID} = ursa ]]; then

  TPN=192
  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=6
  JNPES_c384=8
  THRD_c384=2
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=4

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ ${MACHINE_ID} = linux ]]; then

  TPN=40

  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4

  THRD_cpl_dflt=1
  INPES_cpl_dflt=3
  JNPES_cpl_dflt=8
  WPG_cpl_dflt=6
  OCN_tasks_cpl_dflt=20
  ICE_tasks_cpl_dflt=10
  WAV_tasks_cpl_dflt=20

  THRD_cpl_thrd=2
  INPES_cpl_thrd=3
  JNPES_cpl_thrd=4
  WPG_cpl_thrd=6
  OCN_tasks_cpl_thrd=20
  ICE_tasks_cpl_thrd=10
  WAV_tasks_cpl_thrd=12

elif [[ ${MACHINE_ID} = gaeac5 ]]; then

  TPN=128

  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=6
  JNPES_c384=8
  THRD_c384=1
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=2

  THRD_cpl_atmw_gdas=3
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=264

elif [[ ${MACHINE_ID} = gaeac6 ]]; then

  TPN=192

  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=6
  JNPES_c384=8
  THRD_c384=1
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=2

  THRD_cpl_atmw_gdas=3
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=264
elif [[ ${MACHINE_ID} = derecho ]]; then

  TPN=128
  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4
  INPES_c384=8
  JNPES_c384=6
  THRD_c384=2
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=2

  THRD_cpl_atmw_gdas=2
  INPES_cpl_atmw_gdas=6
  JNPES_cpl_atmw_gdas=8
  WPG_cpl_atmw_gdas=24
  WAV_tasks_atmw_gdas=248

elif [[ ${MACHINE_ID} = noaacloud ]]; then

  if [[ ${PW_CSP} == aws ]]; then
    TPN=36
  elif [[ ${PW_CSP} == azure ]]; then
    TPN=44
  elif [[ ${PW_CSP} == google ]]; then
    TPN=30
  fi

  EXCLUSIVE_NODES=.true.
  INPES_dflt=3
  JNPES_dflt=8
  INPES_thrd=3
  JNPES_thrd=4

  INPES_c384=8
  JNPES_c384=6
  THRD_c384=2
  INPES_c768=8
  JNPES_c768=16
  THRD_c768=2

  THRD_cpl_dflt=1
  INPES_cpl_dflt=3
  JNPES_cpl_dflt=8
  WPG_cpl_dflt=6
  OCN_tasks_cpl_dflt=20
  ICE_tasks_cpl_dflt=10
  WAV_tasks_cpl_dflt=20

  THRD_cpl_thrd=2
  INPES_cpl_thrd=3
  JNPES_cpl_thrd=4
  WPG_cpl_thrd=6
  OCN_tasks_cpl_thrd=20
  ICE_tasks_cpl_thrd=10
  WAV_tasks_cpl_thrd=12

else

  echo "Unknown MACHINE_ID ${MACHINE_ID}"
  exit 1

fi

WLCLK_dflt=30

WLCLK=${WLCLK_dflt}
CMP_DATAONLY=false
nccmp_exclude=""
nccmp_exclude_attr=""

# Defaults for ufs.configure
esmf_logkind="ESMF_LOGKIND_MULTI"
ESMF_THREADING=true
DumpFields="false"
MED_history_n=1000000
RESTART_FH=" "

function set_restart_file_prefix() {
  local restart_file_prefix
  restart_file_prefix=$(date -u +"%Y%m%d.%H0000" -d "${SYEAR}${SMONTH}${SDAY} ${SHOUR} ${FHROT} hours")
  echo "${restart_file_prefix}"
}

function set_restart_file_suffix_secs() {
  local restart_valid_time
  local restart_secs
  local restart_file_suffix_secs
  restart_valid_time=$(date -u +"%Y-%m-%d %H:%M:%S" -d "${SYEAR}${SMONTH}${SDAY} ${SHOUR} ${FHROT} hours")
  restart_secs=$(($(date -u -d "${restart_valid_time}" +%-H) * 3600))
  restart_file_suffix_date="$(date -u -d "${restart_valid_time}" +"%Y-%m-%d")"
  restart_file_suffix_secs="${restart_file_suffix_date}-$(printf "%05d" "${restart_secs}")"
  echo "${restart_file_suffix_secs}"
}

export_fv3_v16() {
  # Add support for v16 test cases. This section
  # will be removed once support for GFSv16 is
  # officially depricated.

  # Load in FV3 values
  export_fv3

  # Replace FV3 variable with old values as needed
  USE_MERRA2=.false.
  WRITE_NSFLIP=.false.

  DIAG_TABLE=diag_table_gfsv16.IN
  FIELD_TABLE=field_table_gfsv16
  FV3_RUN=control_run.IN
  INPUT_NML=control.nml.IN
  CCPP_SUITE=FV3_GFS_v16

  DOGP_CLDOPTICS_LUT=.false.
  DOGP_LWSCAT=.false.
  IAER=111
  ICLIQ_SW=1
  IOVR=1
  IMP_PHYSICS=11
  DNATS=1
  DO_SAT_ADJ=.true.
  LHEATSTRG=.true.
  LSEASPRAY=.false.
  GWD_OPT=1
  DO_UGWP_V0=.false.
  DO_GSL_DRAG_SS=.false.
  SATMEDMF=.false.
  ISATMEDMF=0
  LRADAR=.true.
  LTAEROSOL=.true.
  MRAEROSOL=.false.

  LSM=1
  LANDICE=.true.
  IALB=1
  IEMS=1

  NSTF_NAME=2,1,1,0,5
  FNALBC="'global_snowfree_albedo.bosu.t126.384.190.rg.grb'"
  FNVETC="'global_vegtype.igbp.t126.384.190.rg.grb'"
  FNSOTC="'global_soiltype.statsgo.t126.384.190.rg.grb'"
  FNSOCC="''"
  FNSMCC="'global_soilmgldas.t126.384.190.grb'"
  FNSMCC_control="'global_soilmgldas.statsgo.t1534.3072.1536.grb'"
  FNMSKH_control="'global_slmask.t1534.3072.1536.grb'"
  FNABSC="'global_mxsnoalb.uariz.t126.384.190.rg.grb'"

  RF_CUTOFF=30.0
  FAST_TAU_W_SEC=0.0

  TILEDFIX=.false.
  DO_CA=.false.
  CA_SGS=.false.
}

export_mpas() {
  export_gfs_physics
  # ufs.configure defaults
  UFS_CONFIGURE=ufs.configure.atm.IN
  atm_model=mpas

  #
  MPAS=true
  FV3=false
  S2S=false
  HAFS=false
  AQM=false
  FIRE_BEHAVIOR=false
  DATM_CDEPS=false
  DOCN_CDEPS=false
  DICE_CDEPS=false
  CICE_PRESCRIBED=false
  CDEPS_INLINE=false
  POSTAPP='global'
  USE_MERRA2=.true.
  NESTED=.false.
  BLOCKSIZE=32
  CHKSUM_DEBUG=.false.
  DYCORE_ONLY=.false

  # MPAS dynamical core defaults for RRFS
  MPAS_RESOLUTION=120

  ATM_compute_tasks=4

  #DJS2025 START: We don't need this for MPAS, but to setup the tests we do. CLEAN THIS UP!!!
  #Set defaults if ATMRES and DT_ATMOS are not set
  ATMRES=${ATMRES:-"C96"}
  DT_ATMOS=${DT_ATMOS:-"1800"}

  DAYS=1
  ENS_NUM=1
  SYEAR=2016
  SMONTH=10
  SDAY=03
  SHOUR=00
  SECS=$((10#${SHOUR} * 3600))
  FHMAX=$((DAYS * 24))
  FHCYC=0
  FHROT=0
  LDIAG3D=.false.
  QDIAG3D=.false.
  PRINT_DIFF_PGR=.false.
  MAX_OUTPUT_FIELDS=310
  UPDATE_FULL_OMEGA=.false.
  FHZERO=6
  FHCYC=0
  CPLWAV=.false.
  CPLCHM=.false.
  CPLWAV2ATM=.false.

  #DT_INNER=(Time step)/2
  DT_INNER_c96=360
  DT_INNER_c192=300
  DT_INNER_c384=150
  DT_INNER_c768=75

  if [[ ${DT_ATMOS} = 1800 ]]; then
    default_dt_atmos=1
    DT_INNER=${DT_INNER_c96}
  else
    default_dt_atmos=0
    DT_INNER=${DT_ATMOS}
  fi
  #DJS2025 END:

  # DJS2025: This is needed by rt_utils.sh, but not applicable to MPAS forecasts yet...
  NTILES=1
  QUILTING=.false.
  QUILTING_RESTART=.false.

  # stochastic phsyics (NOT USED in MPAS yet)
  DO_SPPT=.false.
  DO_SHUM=.false.
  DO_SKEB=.false.
  LNDP_TYPE=0
  N_VAR_LNDP=0

  INPES=${INPES_dflt}
  JNPES=${JNPES_dflt}

  # DJS2025: Needed for mpasatm_configure
  RESTART_INTERVAL=0
  ITASKS=1
  OUTPUT_HISTORY=.true.
  HISTORY_FILE_ON_NATIVE_GRID=.true.
  NUM_FILES=2
  FV3ATM_OUTPUT_DIR="./"
  FILENAME_BASE="'atm' 'sfc'"
  OUTPUT_GRID="'mpas'"
  OUTPUT_FILE="'netcdf'"
  ZSTANDARD_LEVEL=0

  DOMAINS_STACK_SIZE=3000000
}
export_mpas_rrfs() {
  # RRFS agnostic MPAS settings
  export_mpas

  # RRFS specific MPAS settings.
  DIAG_TABLE=diag_table_mpas
  FIELD_TABLE=field_table_rrfs_mpas
  FV3_RUN=rrfs_mpas_run.IN
  INPUT_NML=control_rrfs_mpas.nml.IN
  CCPP_SUITE=MPAS_RRFS

  MODEL_CONFIGURE=mpasrrfs_configure.IN
}

export_mpas_gfs() {
  # GFS agnostic MPAS settings.
  export_mpas

  # GFS specific MPAS setting
  DIAG_TABLE=diag_table_mpas
  FIELD_TABLE=field_table_gfsv17_mpas
  FV3_RUN=gfs_mpas_run.IN
  INPUT_NML=control_gfs_mpas.nml.IN
  # Use regional physics for now.
  CCPP_SUITE=MPAS_RRFS

  MODEL_CONFIGURE=mpasgfs_configure.IN
}

export_gfs_physics() {
  # Radiation
  ICLOUD=0
  ICLOUD_BL=1
  IAER=1011
  ICLIQ_SW=2
  IOVR=3
  LFNC_K=-999
  LFNC_P0=-999
  PDFCLD=.false.
  FHSWR=3600.
  FHLWR=3600.
  ICO2=2
  ISUBC_SW=2
  ISUBC_LW=2
  ISOL=2
  LWHTR=.true.
  SWHTR=.true.
  CNVGWD=.true.
  CAL_PRE=.false.
  REDRAG=.true.
  DSPHEAT=.true.
  HYBEDMF=.false.
  # RRTMGP
  DO_RRTMGP=.false.
  DOGP_CLDOPTICS_LUT=.true.
  DOGP_LWSCAT=.true.
  DOGP_SGS_CNV=.true.
  USE_LW_JACOBIAN=.false.
  DAMP_LW_FLUXADJ=.false.
  RRTMGP_LW_PHYS_BLKSZ=2
  EFFR_IN=.true.
  ACTIVE_GASES="'h2o_co2_o3_n2o_ch4_o2'"
  NGASES=6
  LW_FILE_GAS="'rrtmgp-data-lw-g128-210809.nc'"
  LW_FILE_CLOUDS="'rrtmgp-cloud-optics-coeffs-lw.nc'"
  SW_FILE_GAS="'rrtmgp-data-sw-g112-210809.nc'"
  SW_FILE_CLOUDS="'rrtmgp-cloud-optics-coeffs-sw.nc'"
  RRTMGP_NGPTSSW=112
  RRTMGP_NGPTSLW=128
  RRTMGP_NBANDSLW=16
  RRTMGP_NBANDSSW=14

  # Microphysics
  IMP_PHYSICS=8
  NWAT=6
  # GFDL MP
  DNATS=0
  DO_SAT_ADJ=.false.
  LHEATSTRG=.false.
  LSEASPRAY=.true.
  LGFDLMPRAD=.false.
  EFFR_IN=.false.
  # Thompson MP
  LRADAR=.false.
  LTAEROSOL=.false.
  EXT_DIAG_THOMPSON=.false.
  SEDI_SEMI=.true.
  DECFL=10
  # NSSL MP
  NSSL_CCCN=0.6e9
  NSSL_ALPHAH=0.0
  NSSL_ALPHAHL=1.0
  NSSL_HAIL_ON=.false.
  NSSL_CCN_ON=.true.
  NSSL_INVERTCCN=.true.

  # Smoke
  RRFS_SMOKE=.false.
  SMOKE_FORECAST=0
  RRFS_RESTART=NO
  SEAS_OPT=2

  # GWD
  LDIAG_UGWP=.false.
  DO_UGWP=.false.
  DO_TOFD=.false.
  GWD_OPT=2
  DO_UGWP_V0=.true.
  DO_UGWP_V1_W_GSLDRAG=.false.
  DO_UGWP_V0_OROG_ONLY=.false.
  DO_GSL_DRAG_LS_BL=.false.
  DO_GSL_DRAG_SS=.true.
  DO_GWD_OPT_PSL=.false.
  PSL_GWD_DX_FACTOR=6.0
  DO_GSL_DRAG_TOFD=.false.
  DO_UGWP_V1=.false.
  DO_UGWP_V1_OROG_ONLY=.false.
  KNOB_UGWP_SOLVER=2
  KNOB_UGWP_SOURCE=1,1,0,0
  KNOB_UGWP_WVSPEC=1,25,25,25
  KNOB_UGWP_AZDIR=2,4,4,4
  KNOB_UGWP_STOCH=0,0,0,0
  KNOB_UGWP_EFFAC=1,1,1,1
  KNOB_UGWP_DOAXYZ=1
  KNOB_UGWP_DOHEAT=1
  LAUNCH_LEVEL=54
  KNOB_UGWP_DOKDIS=1
  KNOB_UGWP_NDX4LH=1
  KNOB_UGWP_VERSION=0
  KNOB_UGWP_PALAUNCH=275.0e2
  KNOB_UGWP_NSLOPE=1
  KNOB_UGWP_LZMAX=15.750e3
  KNOB_UGWP_LZMIN=0.75e3
  KNOB_UGWP_LZSTAR=2.0e3
  KNOB_UGWP_TAUMIN=0.25e-3
  KNOB_UGWP_TAUAMP=3.0e-3
  KNOB_UGWP_LHMET=200.0e3
  KNOB_UGWP_OROSOLV="'pss-1986'"

  KNOB_UGWP_TAUAMP=3.0e-3
  DO_UGWP_V0_NST_ONLY=.false.

  # GWG resolution dependent settings
  CDMBGWD_c48='0.071,2.1,1.0,1.0'
  CDMBGWD_c96='0.14,1.8,1.0,1.0'
  CDMBGWD_c192='0.23,1.5,1.0,1.0'
  CDMBGWD_c384='1.1,0.72,1.0,1.0'
  CDMBGWD_c768='4.0,0.15,1.0,1.0'

  # set default
  CDMBGWD=${CDMBGWD_c96}

  # PBL
  ISATMEDMF=1
  TRANS_TRAC=.true.
  SATMEDMF=.true.
  HYBEDMF=.false.
  SHINHONG=.false.
  DO_YSU=.false.
  DO_MYNNEDMF=.false.
  HURR_PBL=.false.
  MONINQ_FAC=1.0
  SFCLAY_COMPUTE_FLUX=.false.

  # Shallow/deep convection
  DO_DEEP=.true.
  SHAL_CNV=.true.
  IMFSHALCNV=2
  HWRF_SAMFSHAL=.false.
  IMFDEEPCNV=2
  HWRF_SAMFDEEP=.false.
  RAS=.false.
  RANDOM_CLDS=.false.
  CNVCLD=.true.
  XR_CNVCLD=.false.
  PROGSIGMA=.false.
  BETASCU=8.0
  BETAMCU=1.0
  BETADCU=2.0

  # Aerosol convective scavenging
  FSCAV_AERO='"*:0.3","so2:0.0","msa:0.0","dms:0.0","nh3:0.4","nh4:0.6","bc1:0.6","bc2:0.6","oc1:0.4","oc2:0.4","dust1:0.6","dust2:0.6","dust3:0.6","dust4:0.6","dust5:0.6","seas1:0.5","seas2:0.5","seas3:0.5","seas4:0.5","seas5:0.5"'

  # SFC
  DO_MYJSFC=.false.
  DO_MYNNSFCLAY=.false.
  BL_MYNN_EDMF=1
  BL_MYNN_TKEADVECT=.true.
  BL_MYNN_EDMF_MOM=1

  # LSM
  PRSLRD0=0.
  IVEGSRC=1
  ISOT=1
  LSOIL=4
  LSM=2
  LSOIL_LSM=4
  LANDICE=.false.
  KICE=2
  IALB=2
  IEMS=2
  IOPT_DVEG=4
  IOPT_CRS=2
  IOPT_BTR=1
  IOPT_RUN=1
  IOPT_RAD=3
  IOPT_ALB=1
  IOPT_STC=3
  IOPT_FRZ=1
  IOPT_INF=1
  IOPT_SFC=3
  IOPT_TRS=2
  IOPT_DIAG=2
  IOPT_SNF=4
  IOPT_TBOT=2
  DEBUG=.false.
  NST_ANL=.true.
  PSAUTCO=0.0008,0.0005
  PRAUTCO=0.00015,0.00015

  D2_BG_K1=0.20
  D2_BG_K2=0.04
  PSM_BC=1

  DDDMP=0.1

  # Ozone / stratospheric H2O
  OZ_PHYS_OLD=.true.
  OZ_PHYS_NEW=.false.
  H2O_PHYS=.false.

  # Lake models
  LKM=0                   # 0=no lake, 1=run lake model, 2=run both lake and nsst on lake points
  IOPT_LAKE=2             # 1=flake, 2=clm lake
  LAKEFRAC_THRESHOLD=0.0  # lake fraction must be higher for lake model to run it
  LAKEDEPTH_THRESHOLD=1.0 # lake must be deeper (in meters) for a lake model to run it
  FRAC_ICE=.true.         # should be false for flake, true for clm_lake
}

export_fv3() {
  #Set defaults if ATMRES and DT_ATMOS are not set
  ATMRES=${ATMRES:-"C96"}
  DT_ATMOS=${DT_ATMOS:-"1800"}

  #DT_INNER=(Time step)/2
  DT_INNER_c96=360
  DT_INNER_c192=300
  DT_INNER_c384=150
  DT_INNER_c768=75

  if [[ ${DT_ATMOS} = 1800 ]]; then
    default_dt_atmos=1
    DT_INNER=${DT_INNER_c96}
  else
    default_dt_atmos=0
    DT_INNER=${DT_ATMOS}
  fi

  # ufs.configure defaults
  UFS_CONFIGURE=ufs.configure.atm.IN
  MODEL_CONFIGURE=model_configure.IN
  atm_model=fv3

  POST_ITAG=post_itag_gfs
  POSTXCONFIG=postxconfig-NT-gfs.txt
  POSTXCONFIG_FH00=postxconfig-NT-gfs_FH00.txt

  FV3=true
  S2S=false
  HAFS=false
  AQM=false
  DO_AQM_CANOPY=.false.
  aqm_rc_file=aqm.rc
  FIRE_BEHAVIOR=false
  DATM_CDEPS=false
  DOCN_CDEPS=false
  DICE_CDEPS=false
  CICE_PRESCRIBED=false
  CDEPS_INLINE=false
  POSTAPP='global'
  USE_MERRA2=.true.
  NESTED=.false.
  BLOCKSIZE=32
  CHKSUM_DEBUG=.false.
  DYCORE_ONLY=.false.

  IO_LAYOUT=1,1
  NTILES=6
  INPES=${INPES_dflt}
  JNPES=${JNPES_dflt}
  RESTART_INTERVAL=0
  USE_FV3_ROUTEHANDLES=.false.
  QUILTING=.true.
  QUILTING_RESTART=.true.
  WRITE_GROUP=1
  WRTTASK_PER_GROUP=6
  ITASKS=1
  OUTPUT_HISTORY=.true.
  HISTORY_FILE_ON_NATIVE_GRID=.false.
  WRITE_DOPOST=.false.
  NUM_FILES=2
  FV3ATM_OUTPUT_DIR="./"
  FILENAME_BASE="'atm' 'sfc'"
  OUTPUT_GRID="'cubed_sphere_grid'"
  OUTPUT_FILE="'netcdf'"
  ZSTANDARD_LEVEL=0
  IDEFLATE=0
  QUANTIZE_NSD=0
  ICHUNK2D=0
  JCHUNK2D=0
  ICHUNK3D=0
  JCHUNK3D=0
  KCHUNK3D=0
  IMO=384
  JMO=190
  WRITE_NSFLIP=.true.

  # New damping coefficients made the following
  #   dynamic based on resolution
  N_SPLIT=5
  K_SPLIT=2
  TAU=0.0
  RF_CUTOFF=10.
  FV_SG_ADJ=450

  DZ_MIN=6
  MIN_SEAICE=0.15
  FRAC_GRID=.true.
  MIN_LAKEICE=0.15

  #input file
  FV3_RUN=control_run.IN
  CCPP_SUITE=FV3_GFS_v17_p8
  FIELD_TABLE=field_table_thompson_noaero_tke
  DIAG_TABLE=diag_table_cpld.IN
  INPUT_NML=global_control.nml.IN

  DOMAINS_STACK_SIZE=3000000

  # Coldstart/warmstart
  #rt script for ICs
  MODEL_INITIALIZATION=false
  #namelist variable
  WARM_START=.false.
  READ_INCREMENT=.false.
  RES_LATLON_DYNAMICS="''"
  ATM_IGNORE_RST_CKSUM=.false.
  INCREMENT_FILE_ON_NATIVE_GRID=.false.
  NGGPS_IC=.true.
  EXTERNAL_IC=.true.
  MAKE_NH=.true.
  MOUNTAIN=.false.
  NA_INIT=1
  DO_VORT_DAMP=.true.
  N_SPONGE=42
  NUDGE_QV=.true.
  NUDGE_DZ=.false.
  HYDROSTATIC=.false.
  KORD_MT=9
  KORD_WZ=9
  KORD_TR=9
  KORD_TM=-9
  PHYS_HYDROSTATIC=.false.
  USE_HYDRO_PRESSURE=.false.
  NWAT=6
  NORD=2
  D4_BG=0.12
  VTDM4=0.02
  DELT_MAX=0.002
  EXTERNAL_ETA=.true.
  GFS_PHIL=.false.
  NCEP_IC=.false.
  D_CON=1.
  HORD_MT=5
  HORD_VT=5
  HORD_TM=5
  HORD_DP=-5
  HORD_TR=8
  ADJUST_DRY_MASS=.false.
  DRY_MASS=98320.0
  CONSV_TE=1.
  PRINT_FREQ=6
  NO_DYCORE=.false.

  FILTERED_TERRAIN=.true.
  GFS_DWINDS=.true.

  USE_UFO=.true.
  PRE_RAD=.false.
  TTENDLIM=-999

  # Radiation
  DO_RRTMGP=.false.
  DOGP_CLDOPTICS_LUT=.true.
  DOGP_LWSCAT=.true.
  DOGP_SGS_CNV=.true.
  USE_LW_JACOBIAN=.false.
  DAMP_LW_FLUXADJ=.false.
  RRTMGP_LW_PHYS_BLKSZ=2
  ICLOUD=0
  ICLOUD_BL=1
  IAER=1011
  ICLIQ_SW=2
  IOVR=3
  LFNC_K=-999
  LFNC_P0=-999
  PDFCLD=.false.
  FHSWR=3600.
  FHLWR=3600.

  ICO2=2
  ISUBC_SW=2
  ISUBC_LW=2
  ISOL=2
  LWHTR=.true.
  SWHTR=.true.
  CNVGWD=.true.
  CAL_PRE=.false.
  REDRAG=.true.
  DSPHEAT=.true.
  HYBEDMF=.false.

  # Microphysics
  IMP_PHYSICS=8
  NWAT=6
  # GFDL MP
  DNATS=0
  DO_SAT_ADJ=.false.
  LHEATSTRG=.false.
  LSEASPRAY=.true.
  LGFDLMPRAD=.false.
  EFFR_IN=.false.
  # Thompson MP
  LRADAR=.false.
  LTAEROSOL=.false.
  MRAEROSOL=.false.
  LTHAILAWARE=.false.
  EXT_DIAG_THOMPSON=.false.
  SEDI_SEMI=.true.
  DECFL=10
  # NSSL MP
  NSSL_CCCN=0.6e9
  NSSL_ALPHAH=0.0
  NSSL_ALPHAHL=1.0
  NSSL_HAIL_ON=.false.
  NSSL_CCN_ON=.true.
  NSSL_INVERTCCN=.true.

  # Smoke
  RRFS_SMOKE=.false.
  SMOKE_FORECAST=0
  RRFS_RESTART=NO
  SEAS_OPT=2

  # GWD
  DO_NGW_EC=.false.
  LDIAG_UGWP=.false.
  DO_UGWP=.false.
  DO_TOFD=.false.
  GWD_OPT=2
  DO_UGWP_V0=.true.
  DO_UGWP_V1_W_GSLDRAG=.false.
  DO_UGWP_V0_OROG_ONLY=.false.
  DO_GSL_DRAG_LS_BL=.false.
  DO_GSL_DRAG_SS=.true.
  DO_GWD_OPT_PSL=.false.
  PSL_GWD_DX_FACTOR=6.0
  DO_GSL_DRAG_TOFD=.false.
  DO_UGWP_V1=.false.
  DO_UGWP_V1_OROG_ONLY=.false.
  KNOB_UGWP_SOLVER=2
  KNOB_UGWP_SOURCE=1,1,0,0
  KNOB_UGWP_WVSPEC=1,25,25,25
  KNOB_UGWP_AZDIR=2,4,4,4
  KNOB_UGWP_STOCH=0,0,0,0
  KNOB_UGWP_EFFAC=1,1,1,1
  KNOB_UGWP_DOAXYZ=1
  KNOB_UGWP_DOHEAT=1
  LAUNCH_LEVEL=54
  KNOB_UGWP_DOKDIS=1
  KNOB_UGWP_NDX4LH=1
  KNOB_UGWP_VERSION=0
  KNOB_UGWP_PALAUNCH=275.0e2
  KNOB_UGWP_NSLOPE=1
  KNOB_UGWP_LZMAX=15.750e3
  KNOB_UGWP_LZMIN=0.75e3
  KNOB_UGWP_LZSTAR=2.0e3
  KNOB_UGWP_TAUMIN=0.25e-3
  KNOB_UGWP_TAUAMP=3.0e-3
  KNOB_UGWP_LHMET=200.0e3
  KNOB_UGWP_OROSOLV="'pss-1986'"

  KNOB_UGWP_TAUAMP=3.0e-3
  DO_UGWP_V0_NST_ONLY=.false.

  # resolution dependent settings
  CDMBGWD_c48='0.071,2.1,1.0,1.0'
  CDMBGWD_c96='0.14,1.8,1.0,1.0'
  CDMBGWD_c192='0.23,1.5,1.0,1.0'
  CDMBGWD_c384='1.1,0.72,1.0,1.0'
  CDMBGWD_c768='4.0,0.15,1.0,1.0'

  # set default
  CDMBGWD=${CDMBGWD_c96}

  if [[ ${default_dt_atmos} = 1 ]]; then
    DT_INNER=${DT_INNER_c96}
  else
    DT_INNER=${DT_ATMOS}
  fi

  ISATMEDMF=1
  TRANS_TRAC=.true.

  # PBL
  SATMEDMF=.true.
  HYBEDMF=.false.
  SHINHONG=.false.
  DO_YSU=.false.
  DO_MYNNEDMF=.false.
  HURR_PBL=.false.
  MONINQ_FAC=1.0
  SFCLAY_COMPUTE_FLUX=.false.
  TTE_EDMF=.false.
  CSCALE=1.0
  # Shallow/deep convection
  DO_DEEP=.true.
  SHAL_CNV=.true.
  IMFSHALCNV=2
  HWRF_SAMFSHAL=.false.
  IMFDEEPCNV=2
  HWRF_SAMFDEEP=.false.
  RAS=.false.
  RANDOM_CLDS=.false.
  CNVCLD=.true.
  XR_CNVCLD=.false.
  PROGSIGMA=.false.
  BETASCU=8.0
  BETAMCU=1.0
  BETADCU=2.0

  # Aerosol convective scavenging
  FSCAV_AERO='"*:0.3","so2:0.0","msa:0.0","dms:0.0","nh3:0.4","nh4:0.6","bc1:0.6","bc2:0.6","oc1:0.4","oc2:0.4","dust1:0.6","dust2:0.6","dust3:0.6","dust4:0.6","dust5:0.6","seas1:0.5","seas2:0.5","seas3:0.5","seas4:0.5","seas5:0.5"'

  # SFC
  DO_MYJSFC=.false.
  DO_MYNNSFCLAY=.false.
  BL_MYNN_EDMF=1
  BL_MYNN_TKEADVECT=.true.
  BL_MYNN_EDMF_MOM=1

  # LSM
  PRSLRD0=0.
  IVEGSRC=1
  ISOT=1
  LSOIL=4
  LSM=2
  LSOIL_LSM=4
  LANDICE=.false.
  KICE=2
  IALB=2
  IEMS=2
  IOPT_DVEG=4
  IOPT_CRS=2
  IOPT_BTR=1
  IOPT_RUN=1
  IOPT_RAD=3
  IOPT_ALB=1
  IOPT_STC=3
  IOPT_FRZ=1
  IOPT_INF=1
  IOPT_SFC=3
  IOPT_TRS=2
  IOPT_DIAG=2
  IOPT_SNF=4
  IOPT_TBOT=2
  DEBUG=.false.
  NST_ANL=.true.
  PSAUTCO=0.0008,0.0005
  PRAUTCO=0.00015,0.00015
  EFFR_IN=.true.
  ACTIVE_GASES="'h2o_co2_o3_n2o_ch4_o2'"
  NGASES=6
  LW_FILE_GAS="'rrtmgp-data-lw-g128-210809.nc'"
  LW_FILE_CLOUDS="'rrtmgp-cloud-optics-coeffs-lw.nc'"
  SW_FILE_GAS="'rrtmgp-data-sw-g112-210809.nc'"
  SW_FILE_CLOUDS="'rrtmgp-cloud-optics-coeffs-sw.nc'"
  RRTMGP_NGPTSSW=112
  RRTMGP_NGPTSLW=128
  RRTMGP_NBANDSLW=16
  RRTMGP_NBANDSSW=14

  D2_BG_K1=0.20
  D2_BG_K2=0.04
  PSM_BC=1

  DDDMP=0.1

  # Ozone / stratospheric H2O
  OZ_PHYS_OLD=.true.
  OZ_PHYS_NEW=.false.

  H2O_PHYS=.false.

  # Lake models
  LKM=0                   # 0=no lake, 1=run lake model, 2=run both lake and nsst on lake points
  IOPT_LAKE=2             # 1=flake, 2=clm lake
  LAKEFRAC_THRESHOLD=0.0  # lake fraction must be higher for lake model to run it
  LAKEDEPTH_THRESHOLD=1.0 # lake must be deeper (in meters) for a lake model to run it
  FRAC_ICE=.true.         # should be false for flake, true for clm_lake

  # Tiled Fix files
  TILEDFIX=.true.

  CPL=.false.
  CPLCHM=.false.
  CPLFLX=.false.
  CPLICE=.false.
  CPLWAV=.false.
  CPLWAV2ATM=.false.
  CPLLND=.false.
  CPLLND2ATM=.false.
  USE_MED_FLUX=.false.
  USE_OCEANUV=.false.
  DAYS=1
  NPX=97
  NPY=97
  NPZ=64
  NPZP=65
  NSTF_NAME=2,1,0,0,0
  OUTPUT_FH="12 -1"
  FHZERO=6
  FSICL=0
  FSICS=0

  # Dynamical core
  FV_CORE_TAU=0.
  FAST_TAU_W_SEC=0.2
  DRY_MASS=98320.0

  ENS_NUM=1
  SYEAR=2016
  SMONTH=10
  SDAY=03
  SHOUR=00
  SECS=$((10#${SHOUR} * 3600))
  FHMAX=$((DAYS * 24))
  FHCYC=24
  FHROT=0
  LDIAG3D=.false.
  QDIAG3D=.false.
  PRINT_DIFF_PGR=.false.
  MAX_OUTPUT_FIELDS=310
  UPDATE_FULL_OMEGA=.false.

  # Stochastic physics
  LCNORM=.false.
  PERT_MP=.false.
  PERT_RADTEND=.false.
  PERT_CLDS=.false.

  NEW_LSCALE=.false.
  STOCHINI=.false.
  DO_SPPT=.false.
  DO_SHUM=.false.
  DO_SKEB=.false.
  LNDP_TYPE=0
  N_VAR_LNDP=0
  SKEB=-999.
  SPPT=-999.
  SHUM=-999.
  LNDP_VAR_LIST="'XXX'"
  LNDP_PRT_LIST=-999
  LNDP_MODEL_TYPE=0
  LNDP_TAU=21600,
  LNDP_LSCALE=500000,
  ISEED_LNDP=2010,
  ISEED_SKEB=0
  SKEB_TAU=21600,
  SKEB_LSCALE=250000,
  SKEBNORM=0,
  SKEB_NPASS=30,
  SKEB_VDOF=5,
  ISEED_SHUM=1,
  SHUM_TAU=21600,
  SHUM_LSCALE=500000,
  ISEED_SPPT=20210325000103
  SPPT_TAU=2.16E4
  SPPT_LSCALE=500.E3
  SPPT_LOGIT=.true.,
  SPPT_SFCLIMIT=.true.,
  USE_ZMTNBLCK=.true.
  PBL_TAPER=0,0,0,0.125,0.25,0.5,0.75
  OCNSPPT=-999.
  OCNSPPT_LSCALE=500.E3
  OCNSPPT_TAU=2.16E4
  ISEED_OCNSPPT=20210325000108
  EPBL=-999.
  EPBL_LSCALE=500.E3
  EPBL_TAU=2.16E4
  ISEED_EPBL=20210325000113
  SKEBINT=0
  SHUMINT=0
  SPPTINT=0

  #IAU
  IAU_INC_FILES="''"
  IAU_DELTHRS=0
  IAUFHRS=-1
  IAU_OFFSET=0
  IAU_FILTER_INCREMENTS=.false.

  FH_DFI_RADAR='-2e10'

  #Cellular automata
  DO_CA=.true.
  CA_SGS=.true.
  CA_GLOBAL=.false.
  NCA=1
  NCELLS=5
  NLIVES=12
  NTHRESH=18
  NSEED=1
  NFRACSEED=0.5
  CA_TRIGGER=.true.
  NSPINUP=1
  ISEED_CA=12345

  #waves
  WW3_RSTDTHR=12
  WW3_DT_2_RST="$(printf "%02d" $((WW3_RSTDTHR * 3600)))"
  WW3_OUTDTHR=1
  WW3_DTFLD="$(printf "%02d" $((WW3_OUTDTHR * 3600)))"
  WW3_DTPNT="$(printf "%02d" $((WW3_OUTDTHR * 3600)))"
  WW3_GRD_OUTDIR='./'
  WW3_PNT_OUTDIR='./'
  WW3_RST_OUTDIR='./'
  DTRST=0
  RSTTYPE=T
  GOFILETYPE=1
  POFILETYPE=1
  WW3_OUTPARS="WND HS FP DP PHS PTP PDIR"
  CPLILINE='$'
  ICELINE='$'
  WINDLINE='$'
  CURRLINE='$'
  NFGRIDS=0
  NMGRIDS=1
  WW3GRIDLINE="'glo_1deg'  'no' 'no' 'CPL:native' 'no' 'no' 'no' 'no' 'no' 'no'   1  1  0.00 1.00  F"
  FUNIPNT=T
  IOSRV=1
  FPNTPROC=T
  FGRDPROC=T
  UNIPOINTS='points'
  FLAGMASKCOMP=' F'
  FLAGMASKOUT=' F'
  RUN_BEG="${SYEAR}${SMONTH}${SDAY} $(printf "%02d" $((SHOUR)))0000"
  RUN_END="2100${SMONTH}${SDAY} $(printf "%02d" $((SHOUR)))0000"
  OUT_BEG=${RUN_BEG}
  OUT_END=${RUN_END}
  RST_BEG=${RUN_BEG}
  RST_2_BEG=${RUN_BEG}
  RST_END=${RUN_END}
  RST_2_END=${RUN_END}
  WW3_WLEV='F'
  WW3_CUR='F'
  WW3_ICE='F'
  WW3_IC1='F'
  WW3_IC5='F'
  # ATMW
  WW3_MODDEF=mod_def.glo_1deg
  MESH_WAV=mesh.glo_1deg.nc
  WW3_RSTFLDS=" "
  # ATMA
  AOD_FRQ=060000

  # Regional
  WRITE_RESTART_WITH_BCS=.false.

  # Diagnostics
  PRINT_DIFF_PGR=.false.

  # Coupling
  coupling_interval_fast_sec=0
  CHOUR=06
  MOM6_OUTPUT_DIR=./MOM6_OUTPUT
  MOM6_RESTART_DIR=./RESTART/
  MOM6_RESTART_SETTING=n
  MOM6_HISTFREQ_N=6

  # Following not used for standalone
  USE_CICE_ALB=.false.

  # GFDL Cloud Microphysics
  FTSFS=90
  REIFLAG=2

  # NAM sfc
  FNGLAC="'global_glacier.2x2.grb'"
  FNMXIC="'global_maxice.2x2.grb'"
  FNTSFC="'RTGSST.1982.2012.monthly.clim.grb'"
  FNSNOC="'global_snoclim.1.875.grb'"
  FNZORC="'igbp'"
  FNAISC="'IMS-NIC.blended.ice.monthly.clim.grb'"
  LDEBUG=.false.

  # Land IAU defaults
  DO_LAND_IAU=.false.
  LAND_IAU_FHRS=3,6,9
  LAND_IAU_DELHRS=6
  LAND_IAU_INC_FILES="'sfc_inc',''"
  LSOIL_INCR=3
  LAND_IAU_FILTER_INC=.false.
  LAND_IAU_UPD_STC=.true.
  LAND_IAU_UPD_SLC=.true.
  LAND_IAU_DO_STCSMC_ADJ=.true.
  LAND_IAU_MIN_T_INC=0.0001
  LAND_IAU_MIN_SLC_INC=0.000001
}

# Add section for tiled grid namelist
export_tiled() {
  FNSMCC_control="'global_soilmgldas.statsgo.t1534.3072.1536.grb'"
  FNMSKH_control="'global_slmask.t1534.3072.1536.grb'"
  FNALBC="'${ATMRES}.snowfree_albedo.tileX.nc'"
  FNALBC2="'${ATMRES}.facsf.tileX.nc'"
  FNTG3C="'${ATMRES}.substrate_temperature.tileX.nc'"
  FNVEGC="'${ATMRES}.vegetation_greenness.tileX.nc'"
  FNVETC="'${ATMRES}.vegetation_type.tileX.nc'"
  FNSOTC="'${ATMRES}.soil_type.tileX.nc'"
  FNSOCC="'${ATMRES}.soil_color.tileX.nc'"
  FNSMCC=${FNSMCC_control}
  FNMSKH=${FNMSKH_control}
  FNVMNC="'${ATMRES}.vegetation_greenness.tileX.nc'"
  FNVMXC="'${ATMRES}.vegetation_greenness.tileX.nc'"
  FNSLPC="'${ATMRES}.slope_type.tileX.nc'"
  FNABSC="'${ATMRES}.maximum_snow_albedo.tileX.nc'"
  LSM=2
  LANDICE=".false."
}

export_ugwpv1() {
  DO_UGWP_V1=.true.
  DO_UGWP_V0=.false.
  GWD_OPT=2
  KNOB_UGWP_VERSION=1
  KNOB_UGWP_NSLOPE=1
  DO_GSL_DRAG_LS_BL=.true.
  DO_GSL_DRAG_SS=.true.
  DO_GSL_DRAG_TOFD=.true.
  DO_UGWP_V1_OROG_ONLY=.false.
  DO_UGWP_V0_NST_ONLY=.false.
  LDIAG_UGWP=.false.
  KNOB_UGWP_DOKDIS=2
  KNOB_UGWP_NDX4LH=4

  # Add updated damping and timestep variables
  case "${ATMRES}" in
  "C48")
    if [[ ${default_dt_atmos} = 1 ]]; then export DT_ATMOS=720; fi
    XR_CNVCLD=.false.
    CDMBGWD="0.071,2.1,1.0,1.0"
    CDMBGWD_GSL="40.0,1.77,1.0,1.0"
    KNOB_UGWP_TAUAMP=6.0e-3
    K_SPLIT=1
    N_SPLIT=4
    TAU=10.0
    RF_CUTOFF=100.0
    FV_SG_ADJ=3600
    ;;
  "C96")
    if [[ ${default_dt_atmos} = 1 ]]; then export DT_ATMOS=720; fi
    XR_CNVCLD=.false.
    CDMBGWD="0.14,1.8,1.0,1.0"
    CDMBGWD_GSL="20.0,2.5,1.0,1.0"
    KNOB_UGWP_TAUAMP=3.0e-3
    K_SPLIT=1
    N_SPLIT=4
    TAU=8.0
    RF_CUTOFF=100.0
    FV_SG_ADJ=1800
    ;;
  "C192")
    if [[ ${default_dt_atmos} = 1 ]]; then export DT_ATMOS=600; fi
    XR_CNVCLD=.true.
    CDMBGWD="0.23,1.5,1.0,1.0"
    CDMBGWD_GSL="5.0,5.0,1.0,1.0"
    KNOB_UGWP_TAUAMP=1.5e-3
    K_SPLIT=2
    N_SPLIT=5
    TAU=6.0
    RF_CUTOFF=100.0
    FV_SG_ADJ=1800
    ;;
  "C384")
    if [[ ${default_dt_atmos} = 1 ]]; then export DT_ATMOS=300; fi
    XR_CNVCLD=.true.
    CDMBGWD="1.1,0.72,1.0,1.0"
    CDMBGWD_GSL="5.0,5.0,1.0,1.0"
    KNOB_UGWP_TAUAMP=0.8e-3
    K_SPLIT=2
    N_SPLIT=4
    TAU=4.0
    RF_CUTOFF=100.0
    FV_SG_ADJ=900
    ;;
  "C768")
    if [[ ${default_dt_atmos} = 1 ]]; then export DT_ATMOS=150; fi
    XR_CNVCLD=.true.
    CDMBGWD="4.0,0.15,1.0,1.0"
    CDMBGWD_GSL="2.5,7.5,1.0,1.0"
    KNOB_UGWP_TAUAMP=0.5e-3
    K_SPLIT=2
    N_SPLIT=4
    TAU=3.0
    RF_CUTOFF=300.0
    FV_SG_ADJ=450
    ;;
  "C1152")
    if [[ ${default_dt_atmos} = 1 ]]; then export DT_ATMOS=150; fi
    XR_CNVCLD=.true.
    CDMBGWD="4.0,0.10,1.0,1.0"
    CDMBGWD_GSL="1.67,8.8,1.0,1.0"
    KNOB_UGWP_TAUAMP=0.35e-3
    K_SPLIT=2
    N_SPLIT=6
    TAU=2.5
    RF_CUTOFF=300.0
    FV_SG_ADJ=450
    ;;
  "C3072")
    if [[ ${default_dt_atmos} = 1 ]]; then export DT_ATMOS=90; fi
    XR_CNVCLD=.true.
    CDMBGWD="4.0,0.05,1.0,1.0"
    CDMBGWD_GSL="0.625,14.1,1.0,1.0"
    KNOB_UGWP_TAUAMP=0.13e-3
    K_SPLIT=4
    N_SPLIT=5
    TAU=0.5
    RF_CUTOFF=300.0
    FV_SG_ADJ=300
    ;;
  *)
    echo Invalid model resolution: "${ATMRES}". Please update specified variable ATMRES.
    exit 1
    ;;
  esac

  if [[ ${DO_GSL_DRAG_SS} = .true. ]]; then export CDMBGWD=${CDMBGWD_GSL}; fi
  if [[ ${SEDI_SEMI} = .false. ]]; then
    DT_INNER=$((DT_ATMOS / 2))
  else
    DT_INNER=${DT_ATMOS}
  fi
  default_dt_atmos=0
}

# Defaults for the CICE6 model namelist, mx100
export_cice6() {
  SECS=$((10#${SHOUR} * 3600))
  DT_CICE=${DT_ATMOS}
  CICE_NPT=999
  CICE_RUNTYPE=initial
  CICE_ICE_IC='cice_model.res.nc'
  CICE_RUNID=unknown
  CICE_USE_RESTART_TIME=.false.
  CICE_RESTART_DIR=./RESTART/
  CICE_RESTART_FILE=iced
  # CICE6 warmstarts
  OCNICE_WARMSTART=.false.

  CICE_RESTART_FORMAT='pnetcdf2'
  CICE_RESTART_IOTASKS=-99
  CICE_RESTART_REARR='box'
  CICE_RESTART_ROOT=-99
  CICE_RESTART_STRIDE=-99
  CICE_RESTART_CHUNK=0,0
  CICE_RESTART_DEFLATE=0

  CICE_HISTORY_FORMAT='pnetcdf2'
  CICE_HISTORY_IOTASKS=-99
  CICE_HISTORY_REARR='box'
  CICE_HISTORY_ROOT=-99
  CICE_HISTORY_STRIDE=-99
  CICE_HISTORY_CHUNK=0,0
  CICE_HISTORY_DEFLATE=0
  CICE_HISTORY_PREC=4

  CICE_DUMPFREQ=d
  CICE_DUMPFREQ_N=1000
  CICE_DIAGFREQ=$(((FHMAX * 3600) / DT_CICE))
  CICE_HISTFREQ_N="0, 0, 6, 0, 0"
  CICE_hist_suffix="'x','x','x','x','x'"
  CICE_HIST_AVG=.true.
  CICE_HISTORY_DIR=./history/
  CICE_INCOND_DIR=./history/
  CICE_GRID=grid_cice_NEMS_mx${OCNRES}.nc
  CICE_MASK=kmtu_cice_NEMS_mx${OCNRES}.nc
  CICE_GRIDATM=A
  CICE_GRIDOCN=A
  CICE_GRIDICE=B
  CICE_TR_POND_LVL=.true.
  CICE_RESTART_POND_LVL=.false.
  # setting to true will allow Frazil FW and Salt to be included in fluxes sent to ocean
  CICE_FRAZIL_FWSALT=.true.
  CICE_KTHERM=2
  CICE_TFREEZE_OPTION=mushy
  # SlenderX2
  CICE_NPROC=${ICE_tasks}
  np2=$((CICE_NPROC / 2))
  CICE_BLCKX=$((NX_GLB / np2))
  CICE_BLCKY=$((NY_GLB / 2))
  CICE_DECOMP=slenderX2

  #ds2s
  MESH_DICE=none
  stream_files_dice=none
  CICE_PRESCRIBED=false
  DICE_CDEPS=false

  #To modify aice on restart, "adjust_aice"
  CICE_RESTART_MOD='none'
}

# Defaults for the MOM6 model namelist, mx100
export_mom6() {
  DT_DYNAM_MOM6=1800
  DT_THERM_MOM6=3600
  MOM6_INPUT=MOM_input_100.IN
  MOM6_OUTPUT_DIR=./MOM6_OUTPUT
  MOM6_HISTFREQ_N=6
  MOM6_RESTART_DIR=./RESTART/
  MOM6_RESTART_SETTING=n
  MOM6_RIVER_RUNOFF=False
  MOM6_FRUNOFF=''
  MOM6_CHLCLIM=seawifs_1998-2006_smoothed_2X.nc
  MOM6_USE_LI2016=True
  MOM6_TOPOEDITS=''
  MOM6_HFREEZE=20.0
  MOM6_GUST_CONST=0.02
  MOM6_WRITE_GEOM=2
  # since CPL_SLOW is set to DT_THERM, this should be always be false
  MOM6_THERMO_SPAN=False
  MOM6_USE_WAVES=True
  MOM6_ALLOW_LANDMASK_CHANGES=False
  # MOM6 diag
  MOM6_DIAG_COORD_DEF_Z_FILE=interpolate_zgrid_40L.nc
  MOM6_DIAG_MISVAL='-1e34'
  # MOM6 IAU
  ODA_INCUPD=False
  ODA_INCUPD_NHOURS=6
  ODA_TEMPINC_VAR="'pt_inc'"
  ODA_SALTINC_VAR="'s_inc'"
  ODA_THK_VAR="'h_fg'"
  ODA_INCUPD_UV=False
  ODA_UINC_VAR="'u_inc'"
  ODA_VINC_VAR="'v_inc'"
  # MOM6 stochastics
  DO_OCN_SPPT=False
  PERT_EPBL=False
  OCN_SPPT=-999.
  EPBL=-999.
  # MOM6 warmstarts
  OCNICE_WARMSTART=.false.
  MOM6_INIT_FROM_Z=True
  MOM6_INIT_UV="zero"
  MOM6_WARMSTART_FILE="none"
}

# Defaults for the WW3 global model
export_ww3() {
  WW3_DOMAIN=mx${OCNRES}
  WW3_MODDEF=mod_def.mx${OCNRES}
  WW3_RSTDTHR=3
  WW3_DT_2_RST="$(printf "%02d" $((WW3_RSTDTHR * 3600)))"
  WW3_OUTDTHR=3
  WW3_DTFLD="$(printf "%02d" $((WW3_OUTDTHR * 3600)))"
  WW3_DTPNT="$(printf "%02d" $((WW3_OUTDTHR * 3600)))"
  WW3_WLEV='F'
  WW3_CUR='C'
  WW3_ICE='C'
  WW3_IC1='F'
  WW3_IC5='F'
  WW3_user_histname='false'
  WW3_historync='false'
  WW3_restartnc='true'
  WW3_restart_from_binary='false'
  # For default ufs_configure (fast loop), no added fields reqd
  WW3_RSTFLDS=" "
  # For either history_nc or restart_nc true
  WW3_PIO_FORMAT='pnetcdf'
  WW3_PIO_STRIDE=4
  WW3_PIO_IOTASKS=-99
  WW3_PIO_REARR='box'
  WW3_PIO_ROOT=-99
}

export_fire_behavior() {
  fbh_model=fire_behavior
  FIRE_BEHAVIOR=true
  FIRE_NML=namelist.fire.IN
  CPLFIRE=false
  DT_FIRE=${DT_ATMOS}
  OUTPUT_FS="$(printf "%02d" $((OUTPUT_FH * 3600)))"
  fire_atm_feedback=1.0
  fire_lsm_zcoupling=false
  fire_lsm_zcoupling_ref=60.0
  fire_num_ignitions=1
  fire_print_msg=0
  fire_upwinding=9
  fire_viscosity=0.4
  fire_wind_height=5.0
}

# Defaults for the global coupled
export_cmeps() {
  UFS_CONFIGURE=ufs.configure.s2sw_fast.IN
  med_model=cmeps
  atm_model=fv3
  ocn_model=mom6
  ice_model=cice6
  wav_model=ww3
  lnd_model=noahmp
  coupling_interval_slow_sec=${DT_THERM_MOM6}
  coupling_interval_fast_sec=${DT_ATMOS}
  MESH_OCN=mesh.mx${OCNRES}.nc
  MESH_ICE=mesh.mx${OCNRES}.nc
  MESH_WAV=mesh.${WW3_DOMAIN}.nc
  CPLMODE=ufs.frac
  CMEPS_PIO_FORMAT='pnetcdf'
  CMEPS_PIO_STRIDE=4
  CMEPS_PIO_IOTASKS=-99
  CMEPS_PIO_REARR='box'
  CMEPS_PIO_ROOT=-99
  RUNTYPE=startup
  RESTART_N=${FHMAX}
  RESTART_FH=" "
  CMEPS_RESTART_DIR=./RESTART/
  cap_dbug_flag=0
  WRITE_ENDOFRUN_RESTART=.false.
  # MOM6 attributes
  use_coldstart=false
  use_mommesh=true
  # CICE attributes
  eps_imesh=1.0e-1
  # mediator AO flux
  flux_convergence=0.0
  flux_iteration=2
  flux_scheme=0
  # mediator ocean albedo
  ocean_albedo_limit=0.06
  use_mean_albedos=.false.
  # vector remapping
  MAPUV3D=true
}

export_cpl() {
  FV3=true
  S2S=true
  HAFS=false
  AQM=false
  FIRE_BEHAVIOR=false
  DATM_CDEPS=false
  DOCN_CDEPS=false
  DICE_CDEPS=false
  CICE_PRESCRIBED=false
  CDEPS_INLINE=false
  ULTRALOW=.false.
  DAYS=1

  #model configure
  MODEL_CONFIGURE=model_configure.IN
  SYEAR=2021
  SMONTH=03
  SDAY=22
  SHOUR=06
  CHOUR=06
  FHMAX=24
  FHROT=0
  QUILTING_RESTART=.false.
  WRTTASK_PER_GROUP=${WPG_cpl_dflt}
  WRITE_NSFLIP=.true.
  OUTPUT_FH='6 -1'

  # default atm/ocn/ice resolution
  if [[ ${default_dt_atmos} = 1 ]]; then
    #If default DT_ATMOS is being used, set to 720 for RTs
    DT_ATMOS=720
    DT_INNER=${DT_ATMOS}
  fi
  if [[ -z ${OCNRES+x} || -z ${OCNRES} ]]; then
    OCNRES=100
  fi
  if [[ -z ${ICERES+x} || -z ${ICERES} ]]; then
    ICERES=1.00
  fi
  NX_GLB=360
  NY_GLB=320
  NPZ=127
  NPZP=128

  # Use updated omega calculations if
  #   hydrostatic is set to false
  if [[ "${HYDROSTATIC}" == .false. ]]; then
    UPDATE_FULL_OMEGA=.true.
  fi

  # default resources
  DOMAINS_STACK_SIZE=8000000
  INPES=${INPES_cpl_dflt}
  JNPES=${JNPES_cpl_dflt}
  THRD=${THRD_cpl_dflt}
  OCN_tasks=${OCN_tasks_cpl_dflt}
  ICE_tasks=${ICE_tasks_cpl_dflt}
  WAV_tasks=${WAV_tasks_cpl_dflt}

  # Set tiled file defaults
  export_tiled

  # Set CICE6 component defaults
  export_cice6

  # Set MOM6 component defaults
  export_mom6

  # Set WW3 component defaults
  export_ww3

  # Set CMEPS component defaults
  export_cmeps

  # FV3 defaults
  FRAC_GRID=.true.
  CCPP_SUITE=FV3_GFS_v17_coupled_p8
  INPUT_NML=global_control.nml.IN
  FIELD_TABLE=field_table_thompson_noaero_tke_GOCART
  DIAG_TABLE=diag_table_cpld.IN
  DIAG_TABLE_ADDITIONAL=''
  FIELD_TABLE_ADDITIONAL=''
  FV3_RUN=cpld_control_run.IN
  TILEDFIX=.true.

  FHZERO=6

  IALB=2
  IEMS=2
  LSM=2
  LANDICE=.false.
  IOPT_DVEG=4
  IOPT_CRS=2
  IOPT_RAD=3
  IOPT_ALB=1
  IOPT_STC=3

  IOPT_SFC=3
  IOPT_TRS=2
  IOPT_DIAG=2

  D2_BG_K1=0.20
  D2_BG_K2=0.04
  PSM_BC=1
  DDDMP=0.1

  DZ_MIN=6

  # Merra2 Aerosols & NSST
  USE_MERRA2=.true.
  IAER=1011
  NSTF_NAME=2,0,0,0,0

  LHEATSTRG=.false.
  LSEASPRAY=.true.

  # RRTMGP
  DO_RRTMGP=.false.
  DOGP_CLDOPTICS_LUT=.true.
  DOGP_LWSCAT=.true.
  DOGP_SGS_CNV=.true.

  # CA
  DO_CA=.true.
  CA_SGS=.true.
  CA_GLOBAL=.false.
  NCA=1
  NCELLS=5
  NLIVES=12
  NTHRESH=18
  NSEED=1
  NFRACSEED=0.5
  CA_TRIGGER=.true.
  NSPINUP=1
  ISEED_CA=12345

  FSICL=0
  FSICS=0

  USE_CICE_ALB=.true.
  MIN_SEAICE=1.0e-6
  DNATS=2
  IMP_PHYSICS=8
  LGFDLMPRAD=.false.
  DO_SAT_ADJ=.false.
  SATMEDMF=.true.

  CPLFLX=.true.
  CPLICE=.true.
  CPL=.true.
  CPLWAV=.true.
  CPLWAV2ATM=.true.
  USE_MED_FLUX=.false.
  CPLCHM=.false.
  CPLLND=.false.

  # for FV3: default values will be changed if doing a warm-warm restart
  WARM_START=.false.
  MAKE_NH=.true.
  NA_INIT=1
  EXTERNAL_IC=.true.
  NGGPS_IC=.true.
  MOUNTAIN=.false.
  # gocart inst_aod output; uses AERO_HIST.rc.IN from parm/gocart directory
  AOD_FRQ=060000

  # checkpoint restarts
  RESTART_FILE_PREFIX=''
  RESTART_FILE_SUFFIX_SECS=''
  RT35D=''

  #CDEPS ds2s
  MESH_DICE=none
  stream_files_dice=none
}
export_35d_run() {
  CNTL_DIR=""
  LIST_FILES=""
}
export_datm_cdeps() {
  FV3=false
  S2S=false
  HAFS=false
  AQM=false
  FIRE_BEHAVIOR=false
  DATM_CDEPS=true
  DOCN_CDEPS=false
  CDEPS_INLINE=false
  DAYS=1

  # model configure
  MODEL_CONFIGURE=datm_cdeps_configure.IN
  SYEAR=2011
  SMONTH=10
  SDAY=01
  SHOUR=00
  CHOUR=00
  FHMAX=24
  DT_ATMOS=900
  FHROT=0

  # required but unused
  WARM_START=.false.
  CPLWAV=.false.
  CPLCHM=.false.

  # atm/ocn/ice resolution
  IATM=1760
  JATM=880
  ATM_NX_GLB=${IATM}
  ATM_NY_GLB=${JATM}
  ATMRES="${IATM}x${JATM}"
  OCNRES=100
  ICERES=1.00
  NX_GLB=360
  NY_GLB=320

  # default resources
  ATM_compute_tasks=${ATM_compute_tasks_cdeps_100}
  OCN_tasks=${OCN_tasks_cdeps_100}
  ICE_tasks=${ICE_tasks_cdeps_100}

  # Set CICE6 component defaults
  export_cice6
  # default non-mushy thermo for CICE
  CICE_KTHERM=1
  CICE_TFREEZE_OPTION=linear_salt

  # Set MOM6 component defaults
  export_mom6
  # default no waves
  MOM6_USE_LI2016=False
  MOM6_USE_WAVES=False
  WW3_DOMAIN=''

  # Set CMEPS component defaults
  export_cmeps
  # vector remapping
  MAPUV3D=false
  # default configure
  UFS_CONFIGURE=ufs.configure.datm_cdeps.IN
  atm_model=datm
  CPLMODE=ufs.frac.aoflux

  # datm defaults
  INPUT_NML=input.mom6.nml.IN
  DIAG_TABLE=diag_table_cpld.IN
  DATM_SRC=CFSR
  FILEBASE_DATM=cfsr
  MESH_ATM=mesh.datm.1760x880.nc
  atm_datamode=GEFS
  stream_files=INPUT/${FILEBASE_DATM}.201110.nc
  EXPORT_ALL=.false.
  STREAM_OFFSET=0

  BL_SUFFIX=""
  RT_SUFFIX=""
}

export_hafs_datm_cdeps() {
  FV3=false
  S2S=false
  HAFS=true
  AQM=false
  FIRE_BEHAVIOR=false
  DATM_CDEPS=true
  DOCN_CDEPS=false
  CDEPS_INLINE=false
  INPES=${INPES_dflt}
  JNPES=${JNPES_dflt}
  NTILES=1

  atm_model=datm
  DATM_IN_CONFIGURE=datm_in.IN
  DATM_STREAM_CONFIGURE=hafs_datm.streams.era5.IN
  EXPORT_ALL=.false.
}

export_hafs_docn_cdeps() {
  FV3=true
  S2S=false
  HAFS=true
  AQM=false
  FIRE_BEHAVIOR=false
  DOCN_CDEPS=true
  CDEPS_INLINE=false
  INPES=${INPES_dflt}
  JNPES=${JNPES_dflt}
  NTILES=1

  ocn_model=docn
  ocn_datamode=sstdata
  CMEPS_PIO_FORMAT='pnetcdf'
  CMEPS_PIO_STRIDE=4
  CMEPS_PIO_IOTASKS=-99
  CMEPS_PIO_REARR='box'
  CMEPS_PIO_ROOT=-99
  DOCN_IN_CONFIGURE=docn_in.IN
  DOCN_STREAM_CONFIGURE=hafs_docn.streams.IN
}

export_hafs_regional() {
  FV3=true
  S2S=false
  HAFS=true
  AQM=false
  FIRE_BEHAVIOR=false
  DATM_CDEPS=false
  DOCN_CDEPS=false
  CDEPS_INLINE=false
  INPES=${INPES_dflt}
  JNPES=${JNPES_dflt}
  NTILES=1
  BLOCKSIZE=24

  # model_configure
  SYEAR=2019
  SMONTH=08
  SDAY=29
  SHOUR=00
  SECS=$((10#${SHOUR} * 3600))
  FHMAX=6
  ENS_NUM=1
  DT_ATMOS=900
  CPL=.true.
  RESTART_INTERVAL=0
  FHROT=0
  coupling_interval_fast_sec=0
  QUILTING=.true.
  WRITE_GROUP=1
  WRTTASK_PER_GROUP=6
  OUTPUT_HISTORY=.true.
  WRITE_DOPOST=.false.
  NUM_FILES=2
  FILENAME_BASE="'atm' 'sfc'"
  OUTPUT_GRID="'regional_latlon'"
  OUTPUT_FILE="'netcdf'"
  ZSTANDARD_LEVEL=0
  IDEFLATE=0
  QUANTIZE_NSD=0
  CEN_LON=-62.0
  CEN_LAT=25.0
  LON1=-114.5
  LAT1=-5.0
  LON2=-9.5
  LAT2=55.0
  DLON=0.03
  DLAT=0.03

  # shel.inp
  # input.nml
  CPL_IMP_MRG=.true.
  DIAG_TABLE=diag_table_hafs
  FIELD_TABLE=field_table_hafs

  OCNRES=''
  ICERES=''
  DT_THERM_MOM6=''

  # Set WW3 component defaults
  export_ww3
  # default hafs with no ice
  WW3_DOMAIN=natl_6m
  WW3_MODDEF=mod_def.${WW3_DOMAIN}
  WW3_WLEV='F'
  WW3_ICE='F'
  WW3_OUTPARS="WND HS T01 T02 DIR FP DP PHS PTP PDIR UST CHA USP"
  WW3_RSTFLDS=" "
  WW3_user_histname='false'
  WW3_historync='false'
  WW3_restartnc='true'
  WW3_restart_from_binary='false'
  # For either history_nc or restart_nc true
  WW3_PIO_FORMAT='pnetcdf'
  WW3_PIO_STRIDE=4
  WW3_PIO_IOTASKS=-99
  WW3_PIO_REARR='box'
  WW3_PIO_ROOT=-99

  # Set CMEPS component defaults
  export_cmeps
  # default hafs
  ocn_model=hycom
  CPLMODE=hafs
  MESH_WAV=mesh.hafs.nc
}

export_hafs() {
  export_fv3_v16
  FV3=true
  S2S=false
  HAFS=true
  AQM=false
  FIRE_BEHAVIOR=false
  DATM_CDEPS=false
  DOCN_CDEPS=false
  CDEPS_INLINE=false
  INPES=${INPES_dflt}
  JNPES=${JNPES_dflt}
  NTILES=1
  IMFSHALCNV=2
  IMFDEEPCNV=2
  HYBEDMF=.false.
  SATMEDMF=.true.
  MONINQ_FAC=-1.0
  HURR_PBL=.true.
  ISATMEDMF=1
  IOPT_SFC=1
  IOPT_DVEG=2
  IOPT_CRS=1
  IOPT_RAD=1
  IOPT_ALB=2
  IOPT_STC=1
  LSM=1
  LANDICE=.true.
  DO_GSL_DRAG_LS_BL=.true.
  DO_GSL_DRAG_SS=.true.
  DO_GSL_DRAG_TOFD=.true.
  IMP_PHYSICS=11
  IAER=111
  CNVGWD=.false.
  LTAEROSOL=.false.
  CDMBGWD=1.0,1.0,1.0,1.0
  MRAEROSOL=.false.
  LHEATSTRG=.false.
  LRADAR=.true.
  USE_OCEANUV=.false.

  FV_CORE_TAU=5.
  RF_CUTOFF=30.e2
  RF_CUTOFF_NEST=50.e2

  IS_MOVING_NEST=".false."
  VORTEX_TRACKER=0
  NTRACK=0
  MOVE_CD_X=0
  MOVE_CD_Y=0
  CPL_IMP_MRG=.true.

  OUTPUT_GRID=''
  IMO=''
  JMO=''
  CEN_LON=''
  CEN_LAT=''
  LON1=''
  LAT1=''
  LON2=''
  LAT2=''
  DLON=''
  DLAT=''
  STDLAT1=''
  STDLAT2=''
  NX=''
  NY=''
  DX=''
  DY=''

  OUTPUT_GRID_2=''
  IMO_2=''
  JMO_2=''
  CEN_LON_2=''
  CEN_LAT_2=''
  LON1_2=''
  LAT1_2=''
  LON2_2=''
  LAT2_2=''
  DLON_2=''
  DLAT_2=''
  STDLAT1_2=''
  STDLAT2_2=''
  NX_2=''
  NY_2=''
  DX_2=''
  DY_2=''

  OUTPUT_GRID_3=''
  IMO_3=''
  JMO_3=''
  CEN_LON_3=''
  CEN_LAT_3=''
  LON1_3=''
  LAT1_3=''
  LON2_3=''
  LAT2_3=''
  DLON_3=''
  DLAT_3=''
  STDLAT1_3=''
  STDLAT2_3=''
  NX_3=''
  NY_3=''
  DX_3=''
  DY_3=''

  OUTPUT_GRID_4=''
  IMO_4=''
  JMO_4=''
  CEN_LON_4=''
  CEN_LAT_4=''
  LON1_4=''
  LAT1_4=''
  LON2_4=''
  LAT2_4=''
  DLON_4=''
  DLAT_4=''
  STDLAT1_4=''
  STDLAT2_4=''
  NX_4=''
  NY_4=''
  DX_4=''
  DY_4=''

  OUTPUT_GRID_5=''
  IMO_5=''
  JMO_5=''
  CEN_LON_5=''
  CEN_LAT_5=''
  LON1_5=''
  LAT1_5=''
  LON2_5=''
  LAT2_5=''
  DLON_5=''
  DLAT_5=''
  STDLAT1_5=''
  STDLAT2_5=''
  NX_5=''
  NY_5=''
  DX_5=''
  DY_5=''

  OUTPUT_GRID_6=''
  IMO_6=''
  JMO_6=''
  CEN_LON_6=''
  CEN_LAT_6=''
  LON1_6=''
  LAT1_6=''
  LON2_6=''
  LAT2_6=''
  DLON_6=''
  DLAT_6=''
  STDLAT1_6=''
  STDLAT2_6=''
  NX_6=''
  NY_6=''
  DX_6=''
  DY_6=''

  OUTPUT_FH='3 -1'
}

export_hrrr() {
  export_fv3_v16
  NPZ=127
  NPZP=128
  DT_ATMOS=300
  SYEAR=2021
  SMONTH=03
  SDAY=22
  SHOUR=06
  OUTPUT_GRID='gaussian_grid'
  NSTF_NAME='2,0,0,0,0'
  WRITE_DOPOST=.true.
  IAER=5111
  FHMAX=12

  FRAC_GRID=.false.
  FRAC_ICE=.true.

  FV_CORE_TAU=10.
  RF_CUTOFF=7.5e2

  FV3_RUN=lake_control_run.IN
  CCPP_SUITE=FV3_HRRR
  INPUT_NML=rap.nml.IN
  FIELD_TABLE=field_table_thompson_aero_tke
  NEW_DIAGTABLE=diag_table_rap

  SFCLAY_COMPUTE_FLUX=.true.

  LKM=1
  IOPT_LAKE=2
  IMP_PHYSICS=8
  DNATS=0
  DO_SAT_ADJ=.false.
  LRADAR=.true.
  LTAEROSOL=.true.
  MRAEROSOL=.false.
  IALB=2
  IEMS=2
  HYBEDMF=.false.
  DO_MYNNEDMF=.true.
  DO_MYNNSFCLAY=.true.
  DO_DEEP=.false.
  SHAL_CNV=.false.
  IMFSHALCNV=-1
  IMFDEEPCNV=-1
  LHEATSTRG=.false.
  LSM=3
  LSOIL_LSM=9
  KICE=9

  GWD_OPT=3
  DO_UGWP_V0=.false.
  DO_UGWP_V0_OROG_ONLY=.false.
  DO_GSL_DRAG_LS_BL=.true.
  DO_GSL_DRAG_SS=.true.
  DO_GSL_DRAG_TOFD=.true.
  DO_UGWP_V1=.false.
  DO_UGWP_V1_OROG_ONLY=.false.
}

export_hrrr_conus13km() {
  export_fv3_v16
  SYEAR=2021
  SMONTH=05
  SDAY=12
  SHOUR=16
  FHMAX=2
  DT_ATMOS=120
  RESTART_INTERVAL=1
  QUILTING=.true.
  WRITE_GROUP=1
  WRTTASK_PER_GROUP=6
  NTILES=1
  WRITE_DOPOST=.false.
  OUTPUT_HISTORY=.true.
  OUTPUT_GRID=lambert_conformal
  OUTPUT_FILE="'netcdf'"

  # Revert these two to GFS_typedefs defaults to avoid a crash:
  SEDI_SEMI=.false.
  DECFL=8

  RRFS_SMOKE=.true.
  SEAS_OPT=0

  LKM=1
  SFCLAY_COMPUTE_FLUX=.true.
  IALB=2
  ICLIQ_SW=2
  IEMS=2
  IOVR=3
  KICE=9
  LSM=3
  LSOIL_LSM=9
  DO_MYNNSFCLAY=.true.
  DO_MYNNEDMF=.true.
  HYBEDMF=.false.
  SHAL_CNV=.false.
  DO_SAT_ADJ=.false.
  DO_DEEP=.false.
  CCPP_SUITE='FV3_HRRR'
  INPES=10
  JNPES=8
  NPX=397
  NPY=233
  NPZ=65
  MAKE_NH=.false.
  NA_INIT=0
  DNATS=0
  EXTERNAL_IC=.false.
  NGGPS_IC=.false.
  MOUNTAIN=.true.
  WARM_START=.true.
  READ_INCREMENT=.false.
  RES_LATLON_DYNAMICS="'fv3_increment.nc'"
  NPZP=66
  FHZERO=1.0
  IMP_PHYSICS=8
  LDIAG3D=.false.
  QDIAG3D=.false.
  PRINT_DIFF_PGR=.true.
  FHCYC=0.0
  IAER=1011
  LHEATSTRG=.false.
  RANDOM_CLDS=.false.
  CNVCLD=.false.
  IMFSHALCNV=-1
  IMFDEEPCNV=-1
  CDMBGWD='3.5,1.0'
  DO_SPPT=.false.
  DO_SHUM=.false.
  DO_SKEB=.false.
  LNDP_TYPE=0
  N_VAR_LNDP=0

  GWD_OPT=3
  DO_UGWP_V0=.false.
  DO_UGWP_V0_OROG_ONLY=.false.
  DO_GSL_DRAG_LS_BL=.true.
  DO_GSL_DRAG_SS=.true.
  DO_GSL_DRAG_TOFD=.true.
  DO_UGWP_V1=.false.
  DO_UGWP_V1_OROG_ONLY=.false.

  FV3_RUN=rrfs_warm_run.IN
  INPUT_NML=rrfs_conus13km_hrrr.nml.IN
  FIELD_TABLE=field_table_thompson_aero_tke_smoke
  DIAG_TABLE=diag_table_hrrr
  MODEL_CONFIGURE=model_configure_rrfs_conus13km.IN
  DIAG_TABLE_ADDITIONAL=diag_additional_rrfs_smoke
  FRAC_ICE=.true.
  USE_CDEPS_INLINE=.false.
}

export_rap_common() {
  export_fv3_v16
  NPZ=127
  NPZP=128
  DT_ATMOS=300
  SYEAR=2021
  SMONTH=03
  SDAY=22
  SHOUR=06
  OUTPUT_GRID='gaussian_grid'
  NSTF_NAME='2,0,0,0,0'
  WRITE_DOPOST=.true.
  IAER=5111

  FV_CORE_TAU=10.
  RF_CUTOFF=7.5e2

  FV3_RUN=control_run.IN
  INPUT_NML=rap.nml.IN
  FIELD_TABLE=field_table_thompson_aero_tke

  LHEATSTRG=.false.
  IMP_PHYSICS=8
  DNATS=0
  DO_SAT_ADJ=.false.
  LRADAR=.true.
  LTAEROSOL=.true.
  MRAEROSOL=.false.
  IALB=2
  IEMS=2
  HYBEDMF=.false.
  DO_MYNNEDMF=.true.
  DO_MYNNSFCLAY=.true.
}

export_rap() {
  export_rap_common

  FHMAX=12
  DIAG_TABLE=diag_table_rap
  CCPP_SUITE=FV3_RAP

  IMFSHALCNV=3
  IMFDEEPCNV=3
  LSM=3
  LSOIL_LSM=9
  KICE=9

  GWD_OPT=3
  DO_UGWP_V0=.false.
  DO_UGWP_V0_OROG_ONLY=.false.
  DO_GSL_DRAG_LS_BL=.true.
  DO_GSL_DRAG_SS=.true.
  DO_GSL_DRAG_TOFD=.true.
  DO_UGWP_V1=.false.
  DO_UGWP_V1_OROG_ONLY=.false.
}

export_rrfs_v1() {
  export_rap_common

  CCPP_SUITE=FV3_RRFS_v1beta
  DIAG_TABLE=diag_table_rap_noah

  DO_DEEP=.false.
  SHAL_CNV=.false.
  IMFSHALCNV=-1
  IMFDEEPCNV=-1
  LHEATSTRG=.false.
  LSM=2
  LSOIL_LSM=4
  LANDICE=.false.
}
