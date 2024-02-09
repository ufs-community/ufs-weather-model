help([[
loads UFS Model prerequisites for S4
]])

miniconda3_ver=os.getenv("miniconda3_ver") or "3.8-s4"
load(pathJoin("miniconda", miniconda3_ver))

license_ver=os.getenv("license_ver") or "S4"
load(pathJoin("license_intel",license_ver))


prepend_path("MODULEPATH", "/data/prod/hpc-stack/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.1.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2022.1"
load(pathJoin("hpc-intel", hpc_intel_ver))

hpc_impi_ver=os.getenv("hpc_impi_ver") or "2022.1"
load(pathJoin("hpc-impi", hpc_impi_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "s4.intel")

whatis("Description: UFS build environment")
