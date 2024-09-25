help([[
loads UFS Model prerequisites for NOAA Parallelworks/Intel
]])

prepend_path("MODULEPATH", "/contrib/spack-stack-rocky8/spack-stack-1.6.0/envs/fms-2024.01/install/modulefiles/Core")
prepend_path("MODULEPATH", "/apps/modules/modulefiles")
prepend_path("PATH", "/contrib/EPIC/bin")
load("gnu")
load("stack-intel")
load("stack-intel-oneapi-mpi")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.10.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
unload("gnu")

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "noaacloud.intel")

whatis("Description: UFS build environment")
