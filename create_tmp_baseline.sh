#!/bin/bash

set -e

if [ "$#" -ne 3 ]; then
    echo "Usage: ./create_tmp_baseline OLD_BASELINE_DIR TMP_BASELINE_DIR COPY_S2S_BASELINE"
    echo "Usage: ./create_tmp_baseline     Possible values: COPY_S2S_BASELINE=[Y|N]"
    exit 1
fi

OLD_BASELINE_DIR=$1
TMP_BASELINE_DIR=$2
COPY_S2S_BASELINE=$3

if [ ! -d $OLD_BASELINE_DIR ]; then
    echo "Existing baseline directory $OLD_BASELINE_DIR cannot be found - abort."
    exit 1
fi

if [ -d $TMP_BASELINE_DIR ]; then
    echo "New temporary baseline directory $TMP_BASELINE_DIR exists - abort."
    exit 1
fi

echo " "

if [[ $COPY_S2S_BASELINE =~ ^[Yy]$ ]]; then
    echo "Copy baseline from $OLD_BASELINE_DIR to $TMP_BASELINE_DIR, including S2S baseline"
elif [[ $COPY_S2S_BASELINE =~ ^[Nn]$ ]]; then
    echo "Copy baseline from $OLD_BASELINE_DIR to $TMP_BASELINE_DIR, excluding S2S baseline"
else
    echo "Invalid argument COPY_S2S_BASELINE, possible values: COPY_S2S_BASELINE=[Y|N]"
    exit 1
fi

mkdir -p $TMP_BASELINE_DIR
echo " "

################################################################################################################
#### THIS IS FROM rt.sh, uses RTPWD and NEW_BASELINE                                                        ####
################################################################################################################

RTPWD=$OLD_BASELINE_DIR
NEW_BASELINE=$TMP_BASELINE_DIR
echo "copy baseline inputs from: ${RTPWD}"
echo "                     to:   ${NEW_BASELINE}"

set -x

rsync -a "${RTPWD}"/FV3_* "${NEW_BASELINE}"/
rsync -a "${RTPWD}"/WW3_* "${NEW_BASELINE}"/
rsync -a "${RTPWD}"/DATM* "${NEW_BASELINE}"/

if [[ $COPY_S2S_BASELINE =~ ^[Yy]$ ]]; then
  rsync -a "${RTPWD}"/MOM6_* "${NEW_BASELINE}"/
  rsync -a "${RTPWD}"/CICE_* "${NEW_BASELINE}"/
  rsync -a "${RTPWD}"/CPL_* "${NEW_BASELINE}"/
  rsync -a "${RTPWD}"/BM_* "${NEW_BASELINE}"/
fi

# FIXME: move these namelist files to parm directory
rsync -a "${RTPWD}"/fv3_regional_control/input.nml "${NEW_BASELINE}"/fv3_regional_control/
rsync -a "${RTPWD}"/fv3_regional_quilt/input.nml   "${NEW_BASELINE}"/fv3_regional_quilt/
rsync -a "${RTPWD}"/fv3_regional_c768/input.nml    "${NEW_BASELINE}"/fv3_regional_c768/
rsync -a "${RTPWD}"/fv3_regional_restart/input.nml "${NEW_BASELINE}"/fv3_regional_restart/

rsync -a "${RTPWD}"/fv3_regional_control/model_configure "${NEW_BASELINE}"/fv3_regional_control/
rsync -a "${RTPWD}"/fv3_regional_quilt/model_configure   "${NEW_BASELINE}"/fv3_regional_quilt/
rsync -a "${RTPWD}"/fv3_regional_c768/model_configure    "${NEW_BASELINE}"/fv3_regional_c768/
rsync -a "${RTPWD}"/fv3_regional_restart/model_configure "${NEW_BASELINE}"/fv3_regional_restart/

rsync -a "${RTPWD}"/fv3_regional_control/INPUT     "${NEW_BASELINE}"/fv3_regional_control/
rsync -a "${RTPWD}"/fv3_regional_control/RESTART   "${NEW_BASELINE}"/fv3_regional_control/
rsync -a "${RTPWD}"/fv3_regional_quilt/INPUT       "${NEW_BASELINE}"/fv3_regional_quilt/
rsync -a "${RTPWD}"/fv3_regional_c768/INPUT        "${NEW_BASELINE}"/fv3_regional_c768/
rsync -a "${RTPWD}"/fv3_regional_restart/INPUT     "${NEW_BASELINE}"/fv3_regional_restart/
rsync -a "${RTPWD}"/fv3_stretched/INPUT            "${NEW_BASELINE}"/fv3_stretched/
rsync -a "${RTPWD}"/fv3_stretched_nest/INPUT       "${NEW_BASELINE}"/fv3_stretched_nest/
rsync -a "${RTPWD}"/fv3_stretched_nest_quilt/INPUT "${NEW_BASELINE}"/fv3_stretched_nest_quilt/

################################################################################################################

# Remove old high-resolution SAR input data and move low-resolution data in its place
cd $TMP_BASELINE_DIR
rm -fr FV3_input_data_sar
mv FV3_input_data_sar_25km FV3_input_data_sar

exit 0