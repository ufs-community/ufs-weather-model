#!/bin/bash
set -eux

export RT_COMPILER='intel'
source ../detect_machine.sh
echo "Machine ID: "+$MACHINE_ID
if [[ $MACHINE_ID = hera.* ]]; then
  WORKDIR=/scratch1/NCEPDEV/nems/Brian.Curtis/test
  export MODULEPATH="/apps/lmod/lmod/modulefiles/Core:/apps/modules/modulefiles/Linux:/apps/modules/modulefiles:/opt/cray/modulefiles:/opt/cray/craype/default/modulefiles"
  export PATH=/scratch1/NCEPDEV/nems/emc.nemspara/soft/miniconda3/bin:$PATH
  export PYTHONPATH=/scratch1/NCEPDEV/nems/emc.nemspara/soft/miniconda3/lib/python3.8/site-packages
elif [[ $MACHINE_ID = orion.* ]]; then
  WORKDIR=/work/noaa/nems/bcurtis/test
  export MODULEPATH="/apps/modulefiles/core"
  export PATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/bin:$PATH
  export PYTHONPATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/lib/python3.8/site-packages
elif [[ $MACHINE_ID = jet.* ]]; then
  WORKDIR=/lfs4/HFIP/h-nems/Brian.Curtis/test
  export MODULEPATH="/apps/lmod/lmod/modulefiles/Core:/apps/modules/modulefiles/Linux:/apps/modules/modulefiles"
  export ACCNR="h-nems"
  export PATH=/lfs4/HFIP/hfv3gfs/software/ecFlow-5.3.1/bin:$PATH
  export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/ecFlow-5.3.1/lib/python2.7/site-packages
elif [[ $MACHINE_ID = gaea.* ]]; then
  WORKDIR=/lustre/f2/pdata/ncep/Brian.Curtis/test
  export MODULEPATH="/usw/eslogin/modulefiles-c4:/sw/gaea-cle7/modulefiles/linux-sles15-x86_64:/opt/cray/pe/perftools/7.1.3/modulefiles:/opt/cray/ari/modulefiles:/opt/cray/pe/craype/2.6.3/modulefiles:/opt/cray/pe/modulefiles:/opt/cray/modulefiles:/opt/modulefiles:/sw/common/modulefiles"
  export LOADEDMODULES=$LOADEDMODULES
  export ACCNR="nggps_emc" # This applies to Brian.Curtis, may need change later
  export PATH=/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/envs/ufs-weather-model/bin:$PATH
  export PYTHONPATH=/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/lib/python3.8/site-packages
elif [[ $MACHINE_ID = cheyenne.* ]]; then
  #export PATH=/glade/p/ral/jntp/tools/ecFlow-5.3.1/bin:$PATH
  #export PYTHONPATH=/glade/p/ral/jntp/tools/ecFlow-5.3.1/lib/python2.7/site-packages
  echo "cheyenne not currently supported. automated RT not starting"
  exit 1
else
  echo "No Python Path for this machine. automated RT not starting"
  exit 1
fi

python rt_auto.py -m $MACHINE_ID -w $WORKDIR

exit 0
