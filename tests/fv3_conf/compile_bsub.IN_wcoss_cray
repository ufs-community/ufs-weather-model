#!/bin/bash

#BSUB -oo out
#BSUB -eo err
#BSUB -J @[JBNME]
#BSUB -W 00:30
#BSUB -q @[QUEUE]
#BSUB -P GFS-DEV
#BSUB -extsched "CRAYLINUX[]" -R "1*{select[craylinux && !vnode]} + 24*{select[craylinux && vnode] span [ptile=24]}"
#BSUB -M 500

set -eux

echo "Compile started:  " `date`

aprun -n 1 -j 1 -N 1 -d 24 @[PATHRT]/compile_cmake.sh @[PATHTR] @[MACHINE_ID] "@[MAKE_OPT]" @[COMPILE_NR]

echo "Compile ended:    " `date`
