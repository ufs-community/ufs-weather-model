help([[
loads UFS Model prerequisites for noaacloud/intel
]])


prepend_path("MODULEPATH", "/contrib/EPIC/spack-stack/spack-stack-1.3.0/envs/unified-dev/install/modulefiles/Core")

load("stack-intel/2021.3.0")
load("stack-intel-oneapi-mpi/2021.3.0")
load("stack-python/3.9.12")
load("cmake/3.23.1")

load("ufs_common_spack_debug")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "noaacloud.intel")

whatis("Description: UFS build environment")
