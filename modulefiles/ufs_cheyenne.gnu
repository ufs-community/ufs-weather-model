#%Module

proc ModulesHelp {} {
  puts stderr "\tcit - loads modules required for building and running UFS Model on Cheyenne/GNU"
}

module-whatis "loads UFS Model prerequisites for Cheyenne/GNU"

module load cmake/3.16.4
setenv CMAKE_C_COMPILER mpicc
setenv CMAKE_CXX_COMPILER mpicxx
setenv CMAKE_Fortran_COMPILER mpif90
setenv CMAKE_Platform cheyenne.gnu

# load programming environment
module load ncarenv/1.3
module load gnu/9.1.0
module load mpt/2.22
module load ncarcompilers/0.5.0
module unload netcdf

module use /glade/p/ral/jntp/GMTB/tools/hpc-stack-v1.1.0/modulefiles/stack
module load hpc/1.1.0
module load hpc-gnu/9.1.0
module load hpc-mpt/2.22

module load ufs_common
