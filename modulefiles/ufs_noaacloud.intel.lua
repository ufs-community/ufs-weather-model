help([[
loads UFS Model prerequisites for cloud/intel
]])

prepend_path("MODULE_PATH", "/contrib/spack-stack/envs/ufs-wm/install/modulefiles/Core")

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2021.3.0"
load(pathJoin("stack-intel", hpc_intel_ver))

hpc_impi_ver=os.getenv("hpc_impi_ver") or "2021.3.0"
load(pathJoin("stack-intel-oneapi-mpi", hpc_impi_ver))

load("ufs_common_spack")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "noaacloud.intel")

whatis("Description: UFS build environment")
