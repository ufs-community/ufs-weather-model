#!/bin/bash --login
set -eux
if [ -f "accesstoken.sh" ]; then
  source ./accesstoken.sh
else
  echo "Please create accesstoken.sh (600) with the following content\n"
  echo "export ghapitoken=<GitHub API Token Here>"
  exit 1
fi

export RT_COMPILER='intel'
source ../detect_machine.sh
echo "Machine ID: "+$MACHINE_ID
if [[ $MACHINE_ID = hera.* ]]; then
  WORKDIR=/scratch1/NCEPDEV/nems/Brian.Curtis/test
  export PATH=/scratch1/NCEPDEV/nems/emc.nemspara/soft/miniconda3/bin:$PATH
  export PYTHONPATH=/scratch1/NCEPDEV/nems/emc.nemspara/soft/miniconda3/lib/python3.8/site-packages
elif [[ $MACHINE_ID = orion.* ]]; then
  WORKDIR=/work/noaa/nems/bcurtis/test
  export PATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/bin:$PATH
  export PYTHONPATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/lib/python3.8/site-packages
elif [[ $MACHINE_ID = jet.* ]]; then
  WORKDIR=/lfs4/HFIP/h-nems/Brian.Curtis/test
  export ACCNR="h-nems"
  export PATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/bin:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/bin:$PATH
  export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/lfs4/HFIP/hfv3gfs/software/miniconda3/4.8.3/lib/python3.8/site-packages
elif [[ $MACHINE_ID = gaea.* ]]; then
  WORKDIR=/lustre/f2/pdata/ncep/Brian.Curtis/test
  export LOADEDMODULES=$LOADEDMODULES
  export ACCNR="nggps_emc" # This applies to Brian.Curtis, may need change later
  export PATH=/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/envs/ufs-weather-model/bin:$PATH
  export PYTHONPATH=/lustre/f2/pdata/esrl/gsd/contrib/miniconda3/4.8.3/lib/python3.8/site-packages
elif [[ $MACHINE_ID = cheyenne.* ]]; then
  WORKDIR=/glade/work/heinzell/fv3/ufs-weather-model/auto-rt
  export ACCNR="P48503002"
  # DH* 20210226 temporary workaround to be able to use gnu
  export RT_COMPILER="gnu"
  export MACHINE_ID="cheyenne.gnu"
  # *DH 20210226
  export PATH=/glade/p/ral/jntp/tools/miniconda3/4.8.3/envs/ufs-weather-model/bin:/glade/p/ral/jntp/tools/miniconda3/4.8.3/bin:$PATH
  export PYTHONPATH=/glade/p/ral/jntp/tools/miniconda3/4.8.3/envs/ufs-weather-model/lib/python3.8/site-packages:/glade/p/ral/jntp/tools/miniconda3/4.8.3/lib/python3.8/site-packages
else
  echo "No Python Path for this machine. automated RT not starting"
  exit 1
fi

python rt_auto.py -m $MACHINE_ID -w $WORKDIR

exit 0
