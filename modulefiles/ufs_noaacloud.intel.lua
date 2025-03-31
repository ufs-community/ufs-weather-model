help([[
loads UFS Model prerequisites for NOAA Parallelworks/Intel
]])

prepend_path("MODULEPATH", "/contrib/spack-stack-rocky8/spack-stack-1.9.1/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")
prepend_path("MODULEPATH", "/apps/modules/modulefiles")

gnu_ver=os.getenv("gnu_ver") or ""
load(pathJoin("gnu", gnu_ver))

stack_intel_ver=os.getenv("stack_intel_ver") or ""
load(pathJoin("stack-intel", stack_intel_ver))

stack_intel_oneapi_mpi_ver=os.getenv("stack_intel_oneapi_mpi_ver") or ""
load(pathJoin("stack-intel-oneapi-mpi", stack_intel_oneapi_mpi_ver))

gnu_ver=os.getenv("gnu_ver") or ""
unload(pathJoin("gnu", gnu_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "noaacloud.intel")

whatis("Description: UFS build environment")
