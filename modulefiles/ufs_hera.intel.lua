help([[
loads UFS Model prerequisites for Hera/Intel
]])

prepend_path("MODULEPATH", "/contrib/sutils/modulefiles")
load("sutils")

cmake_ver=os.getenv("cmake_ver") or "3.20.1"
load(pathJoin("cmake", cmake_ver))

intel_ver=os.getenv("intel_ver") or "2022.1.2"
load(pathJoin("intel", intel_ver))

impi_ver=os.getenv("impi_ver") or "2022.1.2"
load(pathJoin("impi", impi_ver))

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2_ncdf492/modulefiles/stack")

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
setenv("CMAKE_Platform", "hera.intel")

prepend_path("MODULEPATH", "/scratch2/BMC/ifi/modulefiles")
try_load("ifi/20230118-intel-2022.1.2")

whatis("Description: UFS build environment")
