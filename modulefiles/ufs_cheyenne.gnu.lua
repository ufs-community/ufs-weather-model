help([[
loads UFS Model prerequisites for Cheyenne/GNU
]])

-- Fix module setup on Cheyenne
remove_path ("MODULEPATH", "/glade/u/apps/ch/modulefiles/default/compilers")
setenv("MODULEPATH_ROOT", "/glade/work/jedipara/cheyenne/spack-stack/modulefiles")

-- Load spack-stack modules
prepend_path("MODULEPATH", "/glade/work/jedipara/cheyenne/spack-stack/modulefiles/compilers")
prepend_path("MODULEPATH", "/glade/work/jedipara/cheyenne/spack-stack/modulefiles/misc")
prepend_path("MODULEPATH", "/glade/work/jedipara/cheyenne/spack-stack/spack-stack-v1/envs/skylab-2.0.0-plus-ufs-gnu-10.1.0/install/modulefiles/Core")

load("miniconda/3.9.12")
load("ecflow/5.8.4")

load("stack-gcc/10.1.0")
load("stack-openmpi/4.1.1")
load("stack-python/3.9.12")

load("ufs-weather-model-env/1.0.0")
load("w3emc/2.9.2")

setenv("CC", "mpicc")
setenv("CXX", "mpicxx")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "cheyenne.gnu")

whatis("Description: UFS build environment")
