#! /usr/bin/env bash
set -eu

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
