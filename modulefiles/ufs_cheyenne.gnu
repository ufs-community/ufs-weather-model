#%Module

proc ModulesHelp {} {
  puts stderr "\tcit - loads modules required for building and running UFS Model on Cheyenne/GNU"
}

module-whatis "loads UFS Model prerequisites for Cheyenne/GNU"

module load cmake/3.22.0
module load python/3.7.9

# load programming environment
module load ncarenv/1.3
module load gnu/10.1.0
module load mpt/2.22
module load ncarcompilers/0.5.0
module unload netcdf

module use /glade/work/epicufsrt/GMTB/tools/gnu/10.1.0/hpc-stack-v1.2.0/modulefiles/stack
module load hpc/1.2.0
module load hpc-gnu/10.1.0
module load hpc-mpt/2.22

module load ufs_common

setenv CC mpicc
setenv CXX mpicxx
setenv FC mpif90
setenv CMAKE_Platform cheyenne.gnu
