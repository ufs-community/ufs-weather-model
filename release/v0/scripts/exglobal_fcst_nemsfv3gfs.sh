#!/bin/ksh
################################################################################
# UNIX Script Documentation Block
# Script name:         exglobal_fcst_fv3gfs.sh.ecf
# Script description:  Runs a global FV3GFS model forecast
#
# Author:   Fanglin Yang       Org: NCEP/EMC       Date: 2016-11-15
# Abstract: This script runs a single GFS forecast with FV3 dynamical core.
#           This script is created based on a C-shell script that GFDL wrote
#           for the NGGPS Phase-II Dycore Comparison Project.
#
# Script history log:
# 2016-11-15  Fanglin Yang   First Version.
# 2017-02-09  Rahul Mahajan  Added warm start and restructured the code.
# 2017-03-10  Fanglin Yang   Updated for running forecast on Cray.
# 2017-03-24  Fanglin Yang   Updated to use NEMS FV3GFS with IPD4 
#
# Attributes:
#   Language: Portable Operating System Interface (POSIX) Shell
#   Machine: WCOSS-CRAY, Theia
################################################################################

#  Set environment.
export VERBOSE=${VERBOSE:-"YES"}
if [ $VERBOSE = YES ] ; then
  echo $(date) EXECUTING $0 $* >&2
  set -x
fi

# Cycling and forecast hour specific parameters
export PSLOT=${PSLOT:-fv3gfs}
export CASE=${CASE:-C768}
export CDATE=${CDATE:-2017032500}
export CDUMP=${CDUMP:-gfs}
export FHMIN=${FHMIN:-0}
export FHMAX=${FHMAX:-240}
export FHOUT=${FHOUT:-3}
export FHZER=${FHZER:-6}
export FHCYC=${FHCYC:-24}

# Directories.
export PTMP=${PTMP:-/gpfs/hps/ptmp}
export STMP=${STMP:-/gpfs/hps/stmp}
export NWPROD=${NWPROD:-${NWROOT:-/nwprod}}
export BASE_DATA=${BASE_DATA:-$NWPROD}
export FIX_DIR=${FIX_DIR:-$BASE_DATA/fix}
export FIX_AM=${FIX_AM:-$FIX_DIR/fix_am}
export FIX_FV3=${FIX_FV3:-$FIX_DIR/fix_fv3}
export DATA=${DATA:-$STMP/$LOGNAME/pr${PSLOT}${CASE}_${CDATE}}    #temporary running directory
export ROTDIR=${ROTDIR:-$PTMP/$LOGNAME/pr${PSLOT}}                #rorating archive directory
export IC_DIR=${IC_DIR:-$PTMP/$LOGNAME/ICs}                       #cold start initial conditions

# Model resolution specific parameters
export DELTIM=${DELTIM:-225}
export layout_x=${layout_x:-8}
export layout_y=${layout_y:-16}
export LEVS=${LEVS:-64}

# Utilities
export NCP=${NCP:-"/bin/cp -p"}
export NLN=${NLN:-"/bin/ln -sf"}
export SEND=${SEND:-"YES"}   #move final result to rotating directory
export ERRSCRIPT=${ERRSCRIPT:-'eval [[ $err = 0 ]]'}
export NDATE=${NDATE:-$NWPROD/util/exec/ndate}

# Other options
export MEMBER=${MEMBER:-"-1"} # -1: control, 0: ensemble mean, >0: ensemble member $MEMBER
export ENS_NUM=${ENS_NUM:-1}

# Model specific stuff
export FCSTEXECDIR=${FCSTEXECDIR:-${EXECDIR:-$BASE_DATA/sorc/fv3gfs.fd/BUILD/bin}}
export FCSTEXEC=${FCSTEXEC:-fv3_gfs.x}
export PARM_FV3DIAG=${PARM_FV3DIAG:-$FV3DIR_RELREASE/parm/parm_fv3diag}

# Model config options
export FCST_LAUNCHER=${FCST_LAUNCHER:-${APRUN:-""}}
export tasks=${tasks:-$((6*layout_x*layout_y))}
export nthreads=${nthreads:-${nth_f:-1}}
export cores_per_node=${cores_per_node:-${task_per_node:-24}}
export ntiles=${ntiles:-6}
export TYPE=${TYPE:-nh}                  # choices:  nh, hydro
export MONO=${MONO:-non-mono}            # choices:  mono, non-mono
export use_hyper_thread=${hyperthread:-".false."}

#-------------------------------------------------------
if [ ! -d $ROTDIR ]; then mkdir -p $ROTDIR; fi
if [ ! -d $DATA ]; then mkdir -p $DATA ;fi
mkdir -p $DATA/RESTART $DATA/INPUT
cd $DATA || exit 8

#-------------------------------------------------------
# member directory
if [ $MEMBER -lt 0 ]; then
  PREINP=$CDUMP
  MEMCHAR=""
else
  PREINP=enkf.$CDUMP
  MEMCHAR=mem`printf %03i $MEMBER`
fi
yyyymmdd=`echo $CDATE | cut -c1-8`
hh=`echo       $CDATE | cut -c9-10`
MEMDIR=$ROTDIR/${PREINP}.$yyyymmdd/$hh/$MEMCHAR
if [ ! -d $MEMDIR ]; then mkdir -p $MEMDIR; fi

#-------------------------------------------------------
# initial conditions
export warm_start=${warm_start:-".false."}
if [ $warm_start = ".false." ]; then
 if [ -d $IC_DIR/${CASE}_$CDATE ]; then
  $NCP $IC_DIR/${CASE}_$CDATE/* $DATA/INPUT/.
 else
  for file in $MEMDIR/INPUT/*.nc; do
    file2=$(echo $(basename $file))
    fsuf=`echo $file2 | cut -c1-3`
    if [ $fsuf = "gfs" -o $fsuf = "sfc" ]; then
      $NLN $file $DATA/INPUT/$file2
    fi
  done
 fi
else
  for file in $MEMDIR/RESTART/*.nc; do
    file2=$(echo $(basename $file))
    $NLN $file $DATA/INPUT/$file2
  done
  $NLN $MEMDIR/RESTART/coupler.res $DATA/INPUT/coupler.res
  export read_increment=${read_increment:-".false."}
  if [ $read_increment == ".true." ]; then
    if [ -f $MEMDIR/$increment_file ]; then
      $NLN $MEMDIR/$increment_file $DATA/INPUT/$increment_file
    else
      export read_increment=".false."
    fi
  fi
fi
nfiles=`ls -1 $DATA/INPUT/* | wc -l`
if [ $nfiles -lt 0 ]; then
  echo "Initial conditions must exist in $DATA/INPUT, ABORT!"
  exit 1
fi

#--------------------------------------------------------------------------
# Grid and orography data
for n in `seq 1 $ntiles`; do
  $NLN $FIX_FV3/$CASE/${CASE}_grid.tile${n}.nc     $DATA/INPUT/${CASE}_grid.tile${n}.nc
  $NLN $FIX_FV3/$CASE/${CASE}_oro_data.tile${n}.nc $DATA/INPUT/oro_data.tile${n}.nc
done
$NLN $FIX_FV3/$CASE/${CASE}_mosaic.nc  $DATA/INPUT/grid_spec.nc

# GFS standard input data
export iems=${iems:-1}
export isol=${isol:-2}
export iaer=${iaer:-111}
export ico2=${ico2:-2}

$NLN $FIX_AM/global_solarconstant_noaa_an.txt  $DATA/solarconstant_noaa_an.txt
$NLN $FIX_AM/global_o3prdlos.f77               $DATA/INPUT/global_o3prdlos.f77
$NLN $FIX_AM/global_sfc_emissivity_idx.txt     $DATA/sfc_emissivity_idx.txt

$NLN $FIX_AM/global_co2historicaldata_glob.txt $DATA/co2historicaldata_glob.txt
$NLN $FIX_AM/co2monthlycyc.txt                 $DATA/co2monthlycyc.txt
if [ $ico2 -gt 0 ]; then
  for file in `ls $FIX_AM/fix_co2_proj/global_co2historicaldata* ` ; do
    $NLN $file $DATA/$(echo $(basename $file) | sed -e "s/global_//g")
  done
fi

$NLN $FIX_AM/global_climaeropac_global.txt     $DATA/aerosol.dat
if [ $iaer -gt 0 ] ; then
  for file in `ls $FIX_AM/global_volcanic_aerosols* ` ; do
    $NLN $file $DATA/$(echo $(basename $file) | sed -e "s/global_//g")
  done
fi

export FNGLAC=${FNGLAC:-"$FIX_AM/global_glacier.2x2.grb"}
export FNMXIC=${FNMXIC:-"$FIX_AM/global_maxice.2x2.grb"}
export FNTSFC=${FNTSFC:-"$FIX_AM/RTGSST.1982.2012.monthly.clim.grb"}
export FNSNOC=${FNSNOC:-"$FIX_AM/global_snoclim.1.875.grb"}
export FNZORC=${FNZORC:-"igbp"}
export FNALBC=${FNALBC:-"$FIX_AM/global_snowfree_albedo.bosu.t1534.3072.1536.rg.grb"}
export FNALBC2=${FNALBC2:-"$FIX_AM/global_albedo4.1x1.grb"}
export FNAISC=${FNAISC:-"$FIX_AM/CFSR.SEAICE.1982.2012.monthly.clim.grb"}
export FNTG3C=${FNTG3C:-"$FIX_AM/global_tg3clim.2.6x1.5.grb"}
export FNVEGC=${FNVEGC:-"$FIX_AM/global_vegfrac.0.144.decpercent.grb"}
export FNVETC=${FNVETC:-"$FIX_AM/global_vegtype.igbp.t1534.3072.1536.rg.grb"}
export FNSOTC=${FNSOTC:-"$FIX_AM/global_soiltype.statsgo.t1534.3072.1536.rg.grb"}
export FNSMCC=${FNSMCC:-"$FIX_AM/global_soilmgldas.t1534.3072.1536.grb"}
export FNMSKH=${FNMSKH:-"$FIX_AM/seaice_newland.grb"}
export FNVMNC=${FNVMNC:-"$FIX_AM/global_shdmin.0.144x0.144.grb"}
export FNVMXC=${FNVMXC:-"$FIX_AM/global_shdmax.0.144x0.144.grb"}
export FNSLPC=${FNSLPC:-"$FIX_AM/global_slope.1x1.grb"}
export FNABSC=${FNABSC:-"$FIX_AM/global_mxsnoalb.uariz.t1534.3072.1536.rg.grb"}

# nstf_name contains the NSST related parameters
# nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled, 2 = NSSTM on and coupled
# nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
# nstf_name(3) : 1 = NSST analysis on, 0 = NSSTM analysis off
# nstf_name(4) : zsea1 in mm
# nstf_name(5) : zsea2 in mm
# nst_anl      : .true. or .false., NSST analysis over lake                       
export nstf_name=${nstf_name:-"2,0,1,0,5"}


#------------------------------------------------------------------
# changeable parameters
# dycore definitions
res=`echo $CASE |cut -c2-5`
resp=`expr $res + 1 `
export npx=$resp
export npy=$resp
export npz=`expr $LEVS - 1 `
export io_layout="1,1"
#export ncols=$(( (${npx}- 1)*(${npy}-1)*3/2 ))

# blocking factor used for threading and general physics performance
#export nyblocks=`expr \( $npy - 1 \) \/ $layout_y `
#export nxblocks=`expr \( $npx - 1 \) \/ $layout_x \/ 32`
#if [ $nxblocks -le 0 ]; then export nxblocks=1 ; fi
export blocksize=${blocksize:-32}

# export the pre-conditioning of the solution
# =0 implies no pre-conditioning
# >0 means new adiabatic pre-conditioning
# <0 means older adiabatic pre-conditioning
export na_init=${na_init:-1}

# variables for controlling initialization of NCEP/NGGPS ICs
export filtered_terrain=${filtered_terrain:-".true."}
export ncep_plevels=${ncep_plevels:-".true."}
export gfs_dwinds=${gfs_dwinds:-".true."}

# determines whether FV3 or GFS physics calculate geopotential
export gfs_phil=${gfs_phil:-".false."}

# determine whether ozone production occurs in GFS physics
export ozcalc=${ozcalc:-".true."}

# export various debug options
export no_dycore=${no_dycore:-".false."}
export dycore_only=${adiabatic:-".false."}
export chksum_debug=${chksum_debug:-".false."}
export print_freq=${print_freq:-6}

if [ ${TYPE} = "nh" ]; then
  # non-hydrostatic options
  export make_nh=".true."
  export hydrostatic=".false."
  export phys_hydrostatic=".false."     # can be tested
  export use_hydro_pressure=".false."   # can be tested
  export consv_te="1."
else
  # hydrostatic options
  export make_nh=".false."
  export hydrostatic=".true."
  export phys_hydrostatic=".false."     # will be ignored in hydro mode
  export use_hydro_pressure=".true."    # have to be .true. in hydro mode
  export consv_te="0."
fi

# time step parameters in FV3
export k_split=2
export n_split=6

if [ ${MONO} = "mono" -o ${MONO} = "monotonic" ];  then
  # monotonic options
  export d_con="1."
  export do_vort_damp=".false."
  if [ ${TYPE} = "nh" ]; then
    # non-hydrostatic
    export hord_mt="10"
    export hord_xx="10"
  else
    # hydrostatic
    export hord_mt="10"
    export hord_xx="10"
  fi
else
  # non-monotonic options
  export d_con="1."
  export do_vort_damp=".true."
  if [ ${TYPE} = "nh" ]; then
    # non-hydrostatic
    export hord_mt="6"
    export hord_xx="6"
  else
    # hydrostatic
    export hord_mt="10"
    export hord_xx="10"
  fi
fi

if [ ${MONO} = "non-mono" -a ${TYPE} = "nh" ]; then
  export vtdm4="0.02"
else
  export vtdm4="0.05"
fi

if [ $warm_start = ".false." ]; then # CHGRES'd GFS analyses
  export external_ic=".true."
  export mountain=".false."
  export read_increment=".false."
  export res_latlon_dynamics='""'
else # warm start from restart file
  export external_ic=".false."
  export mountain=".true."
  export make_nh=".false."
  export na_init=0                
  if [ $read_increment = ".true." ]; then # add increments on the fly to the restarts
    export res_latlon_dynamics="$increment_file"
  else
    export res_latlon_dynamics='""'
  fi
fi

# build the date for curr_date and diag_table from CDATE
export SYEAR=`echo $CDATE | cut -c1-4`
export SMONTH=`echo $CDATE | cut -c5-6`
export SDAY=`echo $CDATE | cut -c7-8`
export SHOUR=`echo $CDATE | cut -c9-10`
export curr_date="${SYEAR},${SMONTH},${SDAY},${SHOUR},0,0"
export restart_secs=${restart_secs:-0}

# copy over the tables
export DIAG_TABLE=${DIAG_TABLE:-$PARM_FV3DIAG/diag_table}
export DATA_TABLE=${DATA_TABLE:-$PARM_FV3DIAG/data_table}
export FIELD_TABLE=${FIELD_TABLE:-$PARM_FV3DIAG/field_table}

# build the diag_table with the experiment name and date stamp
cat > diag_table << EOF
FV3 Forecast
$SYEAR $SMONTH $SDAY $SHOUR 0 0
EOF
cat $DIAG_TABLE >> diag_table

$NCP $DATA_TABLE data_table
$NCP $FIELD_TABLE field_table

#------------------------------------------------------------------
cat > nems.configure <<EOF
 EARTH_component_list: ATM
 ATM_model:            fv3
 runSeq::
   ATM
 ::

EOF


cat > model_configure <<EOF
  total_member:            $ENS_NUM
  PE_MEMBER01:             $tasks
  start_year:              $SYEAR
  start_month:             $SMONTH
  start_day:               $SDAY
  start_hour:              $SHOUR
  start_minute:            0
  start_second:            0
  nhours_fcst:             $FHMAX
  RUN_CONTINUE:            ${RUN_CONTINUE:-".false."}
  ENS_SPS:                 ${ENS_SPS:-".false."}

  dt_atmos:                $DELTIM    
  calendar:                ${calendar:-'julian'}
  memuse_verbose:          ${memuse_verbose:-".false."}
  atmos_nthreads:          $nthreads
  use_hyper_thread:        ${hyperthread:-".false."}
  ncores_per_node:         $cores_per_node
  restart_interval:        ${restart_interval:-0}
EOF

#&coupler_nml
#  months = ${months:-0}
#  days = ${days:-$((FHMAX/24))}
#  hours = ${hours:-$((FHMAX-24*(FHMAX/24)))}
#  dt_atmos = $DELTIM
#  dt_ocean = $DELTIM
#  current_date = $curr_date
#  calendar = 'julian'
#  memuse_verbose = .false.
#  atmos_nthreads = $nthreads
#  use_hyper_thread = ${hyperthread:-".false."}
#  ncores_per_node = $cores_per_node
#  restart_secs = $restart_secs   ##DA
#/



cat > input.nml <<EOF
&amip_interp_nml
  interp_oi_sst = .true.
  use_ncep_sst = .true.
  use_ncep_ice = .false.
  no_anom_sst = .false.
  data_set = 'reynolds_oi'
  date_out_of_range = 'climo'
/

&atmos_model_nml
  blocksize = $blocksize
  chksum_debug = $chksum_debug
  dycore_only = $dycore_only
/

&diag_manager_nml
  prepend_date = .F.   
/

&fms_io_nml
  checksum_required = .false.
  max_files_r = 100
  max_files_w = 100
/

&fms_nml
  clock_grain = 'ROUTINE'
  domains_stack_size = ${domains_stack_size:-115200}
  print_memory_usage = ${print_memory_usage:-".false."}
/

&fv_grid_nml
  grid_file = 'INPUT/grid_spec.nc'
/

&fv_core_nml
  layout = $layout_x,$layout_y
  io_layout = $io_layout
  npx = $npx
  npy = $npy
  ntiles = $ntiles
  npz = $npz
  grid_type = -1
  make_nh = $make_nh
  fv_debug = ${fv_debug:-".false."}
  range_warn = ${range_warn:-".false."}
  reset_eta = .false.
  n_sponge = ${n_sponge:-24}
  nudge_qv = ${nudge_qv:-".true."}
  tau = 5.
  rf_cutoff = 7.5e2
  d2_bg_k1 = 0.15
  d2_bg_k2 = 0.02
  kord_tm = -9
  kord_mt = 9
  kord_wz = 9
  kord_tr = 9
  hydrostatic = $hydrostatic
  phys_hydrostatic = $phys_hydrostatic
  use_hydro_pressure = $use_hydro_pressure
  beta = 0.
  a_imp = 1.
  p_fac = 0.1
  k_split = $k_split
  n_split = $n_split
  nwat = 2
  na_init = $na_init
  d_ext = 0.
  dnats = 0
  fv_sg_adj = 450
  d2_bg = 0.
  nord = 2
  dddmp = 0.1
  d4_bg = 0.12
  vtdm4 = $vtdm4
  delt_max = 0.002
  ke_bg = 0.
  do_vort_damp = $do_vort_damp
  external_ic = $external_ic
  gfs_phil = $gfs_phil
  nggps_ic = ${nggps_ic:-".true."}
  mountain = $mountain
  ncep_ic = ${ncep_ic:-".false."}
  d_con = $d_con
  hord_mt = $hord_mt
  hord_vt = $hord_xx
  hord_tm = $hord_xx
  hord_dp = $hord_xx
  hord_tr = 8
  adjust_dry_mass = .false.
  consv_te = $consv_te
  consv_am = .false.
  fill = .true.
  dwind_2d = .false.
  print_freq = $print_freq
  warm_start = $warm_start
  no_dycore = $no_dycore
  z_tracer = .true.
/
##  res_latlon_dynamics = $res_latlon_dynamics    ###DA
##  read_increment = $read_increment              ###DA

&external_ic_nml
  filtered_terrain = $filtered_terrain
  ncep_plevels = $ncep_plevels
  levp = $LEVS
  gfs_dwinds = $gfs_dwinds
  checker_tr = .false.
  nt_checker = 0
/

##  ntoz        = ${ntoz:-2}
##  ntcw        = ${ntcw:-3}
&gfs_physics_nml
  fhzero      = $FHZER
  ldiag3d     = ${ldiag3d:-.false.}
  fhcyc       = $FHCYC
  nst_anl     = ${nst_anl:-".true."}
  use_ufo     = ${use_ufo:-".true."}
  pre_rad     = ${pre_rad:-".false."}
  ncld        = ${ncld:-1}
  zhao_mic    = ${zhao_mic:-".true."}
  pdfcld      = ${pdfcld:-".flase."}
  fhswr       = ${fhswr:-3600.}
  fhlwr       = ${fhlwr:-3600.}
  ialb        = ${ialb:-1}
  iems        = ${iems:-1}
  IAER        = ${iaer:-111}
  ico2        = ${ico2:-2}
  isubc_sw    = ${isubc_sw:-2}
  isubc_lw    = ${isubc_lw:-2}
  isol        = ${isol:-2}
  lwhtr       = ${lwhtr:-.true.}
  swhtr       = ${swhtr:-.true.}
  cnvgwd      = ${cnvgwd:-.true.}
  shal_cnv    = ${shal_cnv:-.true.}
  cal_pre     = ${cal_pre:-.true.}
  redrag      = ${redrag:-.true.}
  dspheat     = ${dspheat:-.true.}
  hybedmf     = ${hybedmf:-.true.}
  random_clds = ${random_clds:-.true.}
  trans_trac  = ${trans_trac:-.true.}
  cnvcld      = ${cnvcld:-.true.}
  imfshalcnv  = ${imfshalcnv:-2}
  imfdeepcnv  = ${imfdeepcnv:-2}
  cdmbgwd     = ${cdmbgwd:-"3.5,0.25"}
  prslrd0     = ${prslrd0:-0.}
  ivegsrc     = ${ivegsrc:-1}
  isot        = ${isot:-1}
  debug       = ${gfs_phys_debug:-".false."}
  nstf_name   = $nstf_name
/

&nggps_diag_nml
  fdiag = ${fdiag:-$FHOUT}
/

&interpolator_nml
  interp_method = 'conserve_great_circle'
/

&namsfc
  FNGLAC   = '${FNGLAC}' 
  FNMXIC   = '${FNMXIC}'
  FNTSFC   = '${FNTSFC}' 
  FNSNOC   = '${FNSNOC}' 
  FNZORC   = '${FNZORC}' 
  FNALBC   = '${FNALBC}' 
  FNALBC2  = '${FNALBC2}'
  FNAISC   = '${FNAISC}' 
  FNTG3C   = '${FNTG3C}' 
  FNVEGC   = '${FNVEGC}' 
  FNVETC   = '${FNVETC}' 
  FNSOTC   = '${FNSOTC}' 
  FNSMCC   = '${FNSMCC}' 
  FNMSKH   = '${FNMSKH}' 
  FNTSFA   = '${FNTSFA}' 
  FNACNA   = '${FNACNA}' 
  FNSNOA   = '${FNSNOA}' 
  FNVMNC   = '${FNVMNC}' 
  FNVMXC   = '${FNVMXC}'
  FNSLPC   = '${FNSLPC}'
  FNABSC   = '${FNABSC}'
  LDEBUG = .false.
  FSMCL(2) = 99999
  FSMCL(3) = 99999
  FSMCL(4) = 99999
  FTSFS = 90
  FAISS = 99999
  FSNOL = 99999
  FSICL = 99999
  FTSFL = 99999
  FAISL = 99999
  FVETL = 99999
  FSOTL = 99999
  FvmnL = 99999
  FvmxL = 99999
  FSLPL = 99999
  FABSL = 99999
  FSNOS = 99999
  FSICS = 99999
/

EOF

# Add namelist for stochastic physics options for ensemble member forecast
echo "" >> input.nml
if [ $MEMBER -gt 0 ]; then
  if [ "${SET_STP_SEED:-"NO"}" = "YES" ] ; then
      ISEED_SKEB=$((CDATE*1000 + MEMBER*10 + 1))
      ISEED_SHUM=$((CDATE*1000 + MEMBER*10 + 2))
      ISEED_SPPT=$((CDATE*1000 + MEMBER*10 + 3))
      ISEED_VC=$((CDATE*1000   + MEMBER*10 + 4))
  fi
  cat >> input.nml << EOF
&nam_stochy
  ntrunc = ${JCAP:-$((`echo $CASE | cut -c 2-`*2-2))}
  lon_s = ${LONB:-$((`echo $CASE | cut -c 2-`*4))}
  lat_s = ${LATB:-$((`echo $CASE | cut -c 2-`*2))}
  skeb = ${SKEB:-"-999."}
  shum = ${SHUM:-"-999."}
  sppt = ${SPPT:-"-999."}
  vcamp = ${VCAMP:-"-999."}
  iseed_skeb = ${ISEED_SKEB:-${ISEED:-"0"}}
  iseed_shum = ${ISEED_SHUM:-${ISEED:-"0"}}
  iseed_sppt = ${ISEED_SPPT:-${ISEED:-"0"}}
  iseed_vc = ${ISEED_VC:-${ISEED:-"0"}}
  skeb_tau = ${SKEB_TSCALE:-"-999."}
  shum_tau = ${SHUM_TSCALE:-"-999."}
  sppt_tau = ${SPPT_TSCALE:-"-999."}
  vc_tau = ${VC_TSCALE:-"-999."}
  skeb_lscale = ${SKEB_LSCALE:-"-999."}
  shum_lscale = ${SHUM_LSCALE:-"-999."}
  sppt_lscale = ${SPPT_LSCALE:-"-999."}
  vc_lscale = ${VC_LSCALE:-"-999."}
  sppt_logit = ${SPPT_LOGIT:-".true."}
  sppt_sfclimit = ${SPPT_SFCLIMIT:-".true."}
  vc = ${VC:-"0."}
/
EOF
else
cat >> input.nml << EOF
&nam_stochy
/
EOF
fi

#------------------------------------------------------------------
# run the executable
cd $DATA
$NCP $FCSTEXECDIR/$FCSTEXEC $DATA/.
$FCST_LAUNCHER ./$FCSTEXEC 1>&1 2>&2

export ERR=$?
export err=$ERR
$ERRSCRIPT || exit 2

#------------------------------------------------------------------
if [ $SEND = "YES" ]; then
  # Copy model output files
  cd $DATA
  for n in `seq 1 $ntiles`; do
    for file in $DATA/*.tile${n}.nc; do
      $NCP $file $MEMDIR/.
    done
  done

  # Copy model restart files
  cd $DATA/RESTART
  if [ $restart_secs -gt 0 ]; then
    restart_hrs=`echo "$restart_secs / 3600" | bc`
    RDATE=`$NDATE +$restart_hrs $CDATE`
    ryyyymmdd=`echo $RDATE | cut -c1-8`
    rhh=`echo       $RDATE | cut -c9-10`
    RMEMDIR=$ROTDIR/${PREINP}.$ryyyymmdd/$rhh/$MEMCHAR
    mkdir -p $RMEMDIR/RESTART
    for file in ${ryyyymmdd}.${rhh}0000.* ; do
      file2=$(echo $(basename $file) | sed -e "s/${ryyyymmdd}.${rhh}0000.//g")
      $NCP $file $RMEMDIR/RESTART/$file2
    done
  else
    mkdir -p $MEMDIR/RESTART
    for file in * ; do
      $NCP $file $MEMDIR/RESTART/$file
    done
  fi
fi

#------------------------------------------------------------------
# Clean up before leaving
if [ ${KEEPDATA:-YES} = NO ]; then rm -rf $DATA; fi

#------------------------------------------------------------------
set +x
if [ "$VERBOSE" = "YES" ] ; then
  echo $(date) EXITING $0 with return code $err >&2
fi
exit $err
