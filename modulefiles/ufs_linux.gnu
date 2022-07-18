#!/bin/bash

source /usr/share/lmod/6.6/init/bash

module use /home/builder/opt/hpc-modules/modulefiles/stack

module load hpc/1.2.0
module load hpc-gnu/9.3.0
module load hpc-openmpi/4.0.1

module load jasper/2.0.25
module load zlib/1.2.11
module load png/1.6.35

module load hdf5/1.10.6
module load netcdf/4.7.4
module load pio/2.5.3
module load esmf/v8.3.0b09
module load fms/2022.01

module load bacio/2.4.1
module load crtm/2.3.0
module load g2/3.4.5
module load g2tmpl/1.10.0
module load ip/3.3.3
module load sp/2.3.3
module load w3nco/2.4.1

module load gftl-shared/v1.5.0
module load yafyaml/v0.5.1
module load mapl/2.22.0-esmf-v8.3.0b09

##
## load programming environment: compiler, flags, paths
##
export CC=${CC:-mpicc}
export CXX=${CXX:-mpicxx}
export F77=${F77:-mpif77}
export F90=${F90:-mpif90}
export FC=${FC:-mpif90}

##
## set up variables for ../cmake/configure_linux.gnu.cmake
##
export CMAKE_Platform=linux.gnu
