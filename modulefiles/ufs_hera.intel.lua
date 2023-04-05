help([[
loads UFS Model prerequisites for Hera/Intel
]])

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/miniconda3/modulefiles")
miniconda3_ver=os.getenv("miniconda3_ver") or "4.12.0"
load(pathJoin("miniconda3", miniconda3_ver))

prepend_path("MODULEPATH", "/contrib/sutils/modulefiles")
load("sutils")

cmake_ver=os.getenv("cmake_ver") or "3.20.1"
load(pathJoin("cmake", cmake_ver))

intel_ver=os.getenv("intel_ver") or "2022.1.2"
load(pathJoin("intel", intel_ver))

impi_ver=os.getenv("impi_ver") or "2022.1.2"
load(pathJoin("impi", impi_ver))

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/intel-2022.1.2/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2022.1.2"
load(pathJoin("hpc-intel", hpc_intel_ver))

hpc_impi_ver=os.getenv("hpc_impi_ver") or "2022.1.2"
load(pathJoin("hpc-impi", hpc_impi_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "hera.intel")

prepend_path("MODULEPATH", "/scratch2/BMC/ifi/modulefiles")
try_load("ifi/20230118-intel-2022.1.2")

whatis("Description: UFS build environment")
