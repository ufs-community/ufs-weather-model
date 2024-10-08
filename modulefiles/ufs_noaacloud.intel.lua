help([[
loads UFS Model prerequisites for NOAA Parallelworks/Intel
]])

prepend_path("MODULEPATH", "/contrib/spack-stack-rocky8/spack-stack-1.8.0/envs/ue-intel-2021.10.0/install/modulefiles/Core")
prepend_path("MODULEPATH", "/apps/modules/modulefiles")
prepend_path("PATH", "/contrib/EPIC/bin")
load("gnu")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.10.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
unload("gnu")

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "noaacloud.intel")

whatis("Description: UFS build environment")
