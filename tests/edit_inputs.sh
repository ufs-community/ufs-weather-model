#! /usr/bin/env bash
set -eu

function edit_ice_in {

  jday=$(date -d "${SYEAR}-${SMONTH}-${SDAY} ${SHOUR}:00:00" +%j)
  istep0=$(( ((10#$jday-1)*86400 + 10#$SHOUR*3600) / DT_CICE ))

  # assumes processor shape = "slenderX2"
  np2=$((NPROC_ICE/2))
  BLCKX=$((NX_GLB/$np2))
  BLCKY=$((NY_GLB/2))

  sed -e "s/YEAR_INIT/$SYEAR/g" \
      -e "s/ISTEP0/$istep0/g" \
      -e "s/DT_CICE/$DT_CICE/g" \
      -e "s/CICEGRID/$CICEGRID/g" \
      -e "s/CICEMASK/$CICEMASK/g" \
      -e "s/NPROC_ICE/$NPROC_ICE/g" \
      -e "s/NX_GLB/$NX_GLB/g" \
      -e "s/NY_GLB/$NY_GLB/g" \
      -e "s/BLCKX/$BLCKX/g" \
      -e "s/BLCKY/$BLCKY/g" \
      -e "s/RUNTYPE/$RUNTYPE/g" \
      -e "s/RUNID/$RUNID/g" \
      -e "s/CICE_HIST_AVG/$CICE_HIST_AVG/g" \
      -e "s/RESTART_EXT/$RESTART_EXT/g" \
      -e "s/USE_RESTART_TIME/$USE_RESTART_TIME/g" \
      -e "s/DUMPFREQ_N/$DUMPFREQ_N/g" \
      -e "s/DUMPFREQ/$DUMPFREQ/g" \
      -e "s/FRAZIL_FWSALT/$FRAZIL_FWSALT/g"
}

function edit_mom_input {

  sed -e "s/DT_THERM_MOM6/$DT_THERM_MOM6/g" \
      -e "s/DT_DYNAM_MOM6/$DT_DYNAM_MOM6/g" \
      -e "s/MOM6_RIVER_RUNOFF/$MOM6_RIVER_RUNOFF/g" \
      -e "s/MOM6_THERMO_SPAN/$MOM6_THERMO_SPAN/g" \
      -e "s/MOM6_REPRO_LA/$MOM6_REPRO_LA/g" \
      -e "s/MOM6_USE_WAVES/$MOM6_USE_WAVES/g" \
      -e "s/MOM6_ALLOW_LANDMASK_CHANGES/$MOM6_ALLOW_LANDMASK_CHANGES/g" \
      -e "s/NX_GLB/$NX_GLB/g" \
      -e "s/NY_GLB/$NY_GLB/g" \
      -e "s/CHLCLIM/$CHLCLIM/g"
}

function edit_data_table {
  sed -e "s/FRUNOFF/$FRUNOFF/g"
}

function edit_diag_table {
  sed -e "s/YMD/$SYEAR$SMONTH$SDAY/g" \
      -e "s/ATMRES/$ATMRES/g" \
      -e "s/SYEAR/$SYEAR/g" \
      -e "s/SMONTH/$SMONTH/g" \
      -e "s/SDAY/$SDAY/g"
}

function edit_ww3_input { 

  SDATEWW3="${SYEAR}${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
  EDATEWW3="2100${SMONTH}${SDAY} $(printf "%02d" $(( ${SHOUR}  )))0000"
  #note EDATE should be SDATE+FHMAX, but since this requires adding ndate 
  #a work around is to just put a date long in the future as the actual end time is
  #determined by the driver
  DT_2_RST_WAV="$(printf "%02d" $(( ${WW3RSTDTHR}*3600 )))"
  DTFLD_WAV="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"
  DTPNT_WAV="$(printf "%02d" $(( ${WW3OUTDTHR}*3600 )))"

  sed -e "s/NFGRIDS/$NFGRIDS/g" \
      -e "s/NMGRIDS/${NMGRIDS}/g" \
      -e "s/FUNIPNT/T/g" \
      -e "s/IOSRV/1/g" \
      -e "s/FPNTPROC/T/g" \
      -e "s/FGRDPROC/T/g" \
      -e "s/OUTPARS/${OUTPARS_WAV}/g" \
      -e "s/CPLILINE/${CPLILINE}/g" \
      -e "s/UNIPOINTS/'points'/g" \
      -e "s/GRIDLINE/${ww3gline}/g" \
      -e "s/ICELINE/$ICELINE/g" \
      -e "s/CURRLINE/$CURRLINE/g" \
      -e "s/WINDLINE/$WINDLINE/g" \
      -e "s/RUN_BEG/$SDATEWW3/g" \
      -e "s/RUN_END/$EDATEWW3/g" \
      -e "s/OUT_BEG/$SDATEWW3/g" \
      -e "s/OUT_END/$EDATEWW3/g" \
      -e "s/DTFLD/ $DTFLD_WAV/g" \
      -e "s/FLAGMASKCOMP/ F/g" \
      -e "s/FLAGMASKOUT/ F/g" \
      -e "s/GOFILETYPE/ $WW3OUTPUTTYPE/g" \
      -e "s/POFILETYPE/ $WW3OUTPUTTYPE/g" \
      -e "s/DTPNT/ $DTPNT_WAV/g" \
      -e "s/RST_BEG/$SDATEWW3/g" \
      -e "s/RSTTYPE/T/g" \
      -e "s/RST_2_BEG/$SDATEWW3/g" \
      -e "s/DTRST/0/g" \
      -e "s/DT_2_RST/$DT_2_RST_WAV/g" \
      -e "s/RST_END/$EDATEWW3/g" \
      -e "s/RST_2_END/$EDATEWW3/g" 
}
