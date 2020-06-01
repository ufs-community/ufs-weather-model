#!/bin/bash

#BSUB -oo out
#BSUB -eo err
#BSUB -J @[JBNME]
#BSUB -W 00:45
#BSUB -q @[QUEUE]
#BSUB -P GFS-DEV
#BSUB -n 1
#BSUB -R affinity[core(1)]
#BSUB -R rusage[mem=8192]

set -eux

echo "Compile started:  " `date`

@[PATHRT]/compile_cmake.sh @[PATHTR] @[MACHINE_ID] "@[MAKE_OPT]" @[COMPILE_NR]

echo "Compile ended:    " `date`
