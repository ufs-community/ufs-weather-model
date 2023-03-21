help([[
  This module loads libraries required for building and running UFS Weather Model 
  on the NOAA RDHPC machine Gaea using Intel-2022.1.2
]])

whatis([===[Loads libraries needed for building the UFS Weather Model and debug on Gaea ]===])

prepend_path("MODULEPATH", "/lustre/f2/dev/role.epic/contrib/modulefiles")
load(pathJoin("miniconda3",os.getenv("miniconda_ver") or "4.12.0"))

load(pathJoin("cmake", os.getenv("cmake_ver") or "3.20.1"))

prepend_path("MODULEPATH","/lustre/f2/dev/role.epic/contrib/hpc-stack/intel-2021.3.0_noarch/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))
load(pathJoin("intel", os.getenv("intel_ver") or "2021.3.0"))
load(pathJoin("hpc-intel", os.getenv("hpc_intel_ver") or "2021.3.0"))
load(pathJoin("hpc-cray-mpich", os.getenv("hpc_cray_mpich_ver") or "7.7.11"))
load(pathJoin("gcc", os.getenv("gcc_ver") or "8.3.0"))
load(pathJoin("libpng", os.getenv("libpng_ver") or "1.6.37"))

-- needed for WW3 build
load(pathJoin("gcc", os.getenv("gcc_ver") or "8.3.0"))
-- Needed at runtime:
load("alps")
load("rocoto")

load("ufs_common_debug")

setenv("CC","cc")
setenv("FC","ftn")
setenv("CXX","CC")
setenv("CMAKE_C_COMPILER","cc")
setenv("CMAKE_CXX_COMPILER","CC")
setenv("CMAKE_Fortran_COMPILER","ftn")
setenv("CMAKE_Platform","gaea.intel")
