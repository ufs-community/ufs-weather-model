#%Module

proc ModulesHelp {} {
  puts stderr "\tcit - loads modules required for building and running UFS Model on Hera/GNU"
}

module-whatis "loads UFS Model prerequisites for Hera/GNU"

module use /contrib/sutils/modulefiles
module load sutils

module load cmake/3.16.1
setenv CMAKE_C_COMPILER mpicc
setenv CMAKE_CXX_COMPILER mpicxx
setenv CMAKE_Fortran_COMPILER mpif90
setenv CMAKE_Platform hera.gnu

module use /scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/modulefiles/stack
module load hpc/1.1.0
module load hpc-gnu/9.2.0
module load hpc-mpich/3.3.2

module load ufs_common
