#!/bin/bash

# This is awful, but quick to implement

mkdir -p ../Noah-Comp/ccpphys_files/
ccpp_files="funcphys.f90
machine.F
namelist_soilveg.f
physcons.F90
set_soilveg.f
sfc_diff.f
sfc_drv.f
sfc_drv_loop.F90
surface_perturbation.F90
"

for file in $ccpp_files; do
    echo $file
    ln -s ../../FV3/ccpp/physics/physics/$file ../Noah-Comp/ccpphys_files/$file
done
