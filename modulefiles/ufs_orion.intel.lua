help([[
loads UFS Model prerequisites for Orion/Intel
]])

load("contrib")
load("noaatools")

cmake_ver=os.getenv("cmake_ver") or "3.22.1"
load(pathJoin("cmake", cmake_ver))

prepend_path("MODULEPATH", "/work/noaa/epic-ps/role-epic-ps/hpc-stack/libs/intel-2022.1.2/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2022.1.2"
load(pathJoin("hpc-intel", hpc_intel_ver))

hpc_impi_ver=os.getenv("hpc_impi_ver") or "2022.1.2"
load(pathJoin("hpc-impi", hpc_impi_ver))

scotch_ver=os.getenv("scotch_ver") or "7.0.3"
load(pathJoin("scotch", scotch_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "orion.intel")

whatis("Description: UFS build environment")
