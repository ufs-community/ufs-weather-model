#!/bin/ksh
set -ax

#--------------------------------------
#-- remap FV3 6 tiles to global array
#-- Fanglin Yang, October 2016
#--------------------------------------

export CDATE=${CDATE:-"2016100300"}
export CASE=${CASE:-"C192"}           # C48 C96 C192 C384 C768 C1152 C3072
export GG=${master_grid:-"0p25deg"}   # 1deg 0p5deg 0p25deg 0p125deg

pwd=$(pwd)
export DATA=${DATA:-$pwd}
export NWPROD=${NWPROD:-$pwd}
export BASE_DATA=${BASE_DATA:-$NWPROD}
export REMAPEXE=${REMAPEXE:-$BASE_DATA/exec/fregrid_parallel}
export IPD4=${IPD4:-NO}

cycn=`echo $CDATE | cut -c 9-10`
export TCYC=${TCYC:-".t${cycn}z."}
export CDUMP=${CDUMP:-gfs}
export PREFIX=${PREFIX:-${CDUMP}${TCYC}}

#--------------------------------------------------
export FIX_FV3=${FIX_FV3:-$BASE_DATA/fix/fix_fv3}
export grid_loc=$FIX_FV3/$CASE/${CASE}_mosaic.nc
export weight_file=$FIX_FV3/$CASE/remap_weights_${CASE}_${GG}.nc

export REMAP_LAUNCHER=${REMAP_LAUNCHER:-${APRUN:-""}}
export NTHREADS_REMAP=${NTHREADS_REMAP:-${NTHREADS:-1}}

#--------------------------------------------------
if [ $GG = 1deg    ];  then  export nlon=360  ; export nlat=180  ; fi
if [ $GG = 0p5deg  ];  then  export nlon=720  ; export nlat=360  ; fi
if [ $GG = 0p25deg ];  then  export nlon=1440 ; export nlat=720  ; fi
if [ $GG = 0p125deg ]; then  export nlon=2880 ; export nlat=1440 ; fi

#--------------------------------------------------
hgt=h; if [ $IPD4 = YES ]; then hgt=z; fi

#--for non-hydrostatic case
export atmos_4xdaily_nh="slp, vort850, vort200,\
        us, u1000, u850, u700, u500, u200, u100, u50, u10,\
        vs, v1000, v850, v700, v500, v200, v100, v50, v10,\
        tm, t1000, t850, t700, t500, t200, t100, t50, t10,\
        ${hgt}1000, ${hgt}850, ${hgt}700, ${hgt}500, ${hgt}200, ${hgt}100, ${hgt}50, ${hgt}10,\
        q1000, q850, q700, q500, q200, q100, q50, q10,\
        rh1000, rh850, rh700, rh500, rh200,\
        omg1000, omg850, omg700, omg500, omg200, omg100, omg50, omg10,\
        w700,w850,w500, w200"

#--for hydrostatic case
export atmos_4xdaily_hy="slp, vort850, vort200,\
        us, u1000, u850, u700, u500, u200, u100, u50, u10,\
        vs, v1000, v850, v700, v500, v200, v100, v50, v10,\
        tm, t1000, t850, t700, t500, t200, t100, t50, t10,\
        ${hgt}1000, ${hgt}850, ${hgt}700, ${hgt}500, ${hgt}200, ${hgt}100, ${hgt}50, ${hgt}10,\
        q1000, q850, q700, q500, q200, q100, q50, q10,\
        rh1000, rh850, rh700, rh500, rh200,\
        omg1000, omg850, omg700, omg500, omg200, omg100, omg50, omg10,\
        w700"

export nggps2d_nh="ALBDOsfc, CPRATsfc, PRATEsfc, DLWRFsfc, ULWRFsfc,\
        DSWRFsfc, USWRFsfc, DSWRFtoa, USWRFtoa, ULWRFtoa,\
        GFLUXsfc, HGTsfc, HPBLsfc, ICECsfc, SLMSKsfc,\
        LHTFLsfc, SHTFLsfc, PRESsfc, PWATclm, SOILM,\
        SOILW1, SOILW2, SOILW3, SOILW4, SPFH2m,\
        TCDCclm, TCDChcl, TCDClcl, TCDCmcl,\
        SOILT1, SOILT2, SOILT3, SOILT4,\
        TMP2m, TMPsfc, UGWDsfc, VGWDsfc, UFLXsfc,\
        VFLXsfc, UGRD10m, VGRD10m, WEASDsfc, SNODsfc,\
        ZORLsfc, VFRACsfc, F10Msfc, VTYPEsfc, STYPEsfc"
export nggps2d_hy="$nggps2d_nh"

export nggps3d_nh="ucomp, vcomp, temp, delp, sphum, o3mr, clwmr, nhpres, w, delz"     #for non-hydrostatic case
export nggps3d_hy="ucomp, vcomp, temp, delp, sphum, o3mr, clwmr, hypres"              #for hydrostatic case

#--------------------------------------------------
cd $DATA || exit 8

testfile=${CDATE}0000.nggps3d.tile4.nc
if [ ! -s $testfile ]; then testfile=nggps3d.tile4.nc ; fi
nhrun=$(ncdump -c $testfile | grep nhpres)
nhrun=$?

export OMP_NUM_THREADS=$NTHREADS_REMAP

#--------------------------------------------------
err=0
for type in atmos_4xdaily nggps2d nggps3d ; do

  export in_file="${CDATE}0000.${type}"
  if [ ! -s $in_file ]; then export in_file="${type}"; fi
  export out_file=${PREFIX}${type}.${GG}.nc
  [[ -s $DATA/$out_file ]] && rm -f $DATA/$out_file
  if [ $nhrun -eq 0 ]; then
    export fld=$(eval echo \${${type}_nh})
  else
    export fld=$(eval echo \${${type}_hy})
  fi

  $REMAP_LAUNCHER $REMAPEXE --input_dir $DATA \
                         --input_file $in_file \
                         --output_dir $DATA \
                         --output_file $out_file \
                         --input_mosaic $grid_loc \
                         --scalar_field "$fld" \
                         --interp_method conserve_order1 \
                         --remap_file $weight_file \
                         --nlon $nlon \
                         --nlat $nlat
  rc=$?
  ((err+=$rc))

done

echo $(date) EXITING $0 with return code $err >&2
exit $err

