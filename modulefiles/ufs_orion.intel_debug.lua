help([[
loads UFS Model prerequisites for Orion/Intel
]])

load("contrib")
load("noaatools")

cmake_ver=os.getenv("cmake_ver") or "3.22.1"
load(pathJoin("cmake", cmake_ver))

python_ver=os.getenv("python_ver") or "3.7.5"
load(pathJoin("python", python_ver))

--prepend_path("MODULEPATH", "/apps/contrib/NCEP/libs/hpc-stack/modulefiles/stack")
prepend_path("MODULEPATH", "/work/noaa/epic-ps/hpc-stack/libs/intel/2022.1.2/modulefiles/stack")

--hpc_ver=os.getenv("hpc_ver") or "1.1.0"
hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2022.1.2"
load(pathJoin("hpc-intel", hpc_intel_ver))

hpc_impi_ver=os.getenv("hpc_impi_ver") or "2022.1.2"
load(pathJoin("hpc-impi", hpc_impi_ver))

load("ufs_common_debug")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "orion.intel")

whatis("Description: UFS build environment")
