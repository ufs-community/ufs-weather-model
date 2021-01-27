#!/bin/bash
set -eux

export RT_COMPILER='intel'
source ../detect_machine.sh
echo "Machine ID: "+$MACHINE_ID
export MACHINE_ID=$MACHINE_ID
if [[ $MACHINE_ID = hera.* ]]; then
  export PATH=/scratch1/NCEPDEV/nems/emc.nemspara/soft/miniconda3/bin:$PATH
  export PYTHONPATH=/scratch1/NCEPDEV/nems/emc.nemspara/soft/miniconda3/lib/python3.8/site-packages
  python rt_auto.py
elif [[ $MACHINE_ID = orion.* ]]; then
  export PATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/bin:$PATH
  export PYTHONPATH=/work/noaa/nems/emc.nemspara/soft/miniconda3/lib/python3.8/site-packages
  python rt_auto.py
elif [[ $MACHINE_ID = jet.* ]]; then
  export PATH=/lfs4/HFIP/hfv3gfs/software/ecFlow-5.3.1/bin:$PATH
  export PYTHONPATH=/lfs4/HFIP/hfv3gfs/software/ecFlow-5.3.1/lib/python2.7/site-packages
  python rt_auto.py
elif [[ $MACHINE_ID = gaea.* ]]; then
  module load cray-python/3.7.3.2
  export PATH=/lustre/f2/pdata/esrl/gsd/contrib/ecFlow-5.3.1/bin:$PATH
  export PYTHONPATH=/lustre/f2/pdata/esrl/gsd/contrib/ecFlow-5.3.1/lib/python3.7/site-packages
  python rt_auto.py
elif [[ $MACHINE_ID = cheyenne.* ]]; then
  export PATH=/glade/p/ral/jntp/tools/ecFlow-5.3.1/bin:$PATH
  export PYTHONPATH=/glade/p/ral/jntp/tools/ecFlow-5.3.1/lib/python2.7/site-packages
  python rt_auto.py
else
  echo "No Python Path for this machine. automated RT not starting"
fi
