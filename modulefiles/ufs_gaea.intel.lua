help([[
  This module loads libraries required for building and running UFS Weather Model 
  on the NOAA RDHPC machine Gaea using Intel-2022.0.2
]])

whatis([===[Loads libraries needed for building the UFS Weather Model on Gaea ]===])

load_any(pathJoin("cmake", os.getenv("cmake_ver") or "3.20.1"),"cmake")

prepend_path("MODULEPATH","/lustre/f2/dev/role.epic/contrib/hpc-stack/intel-2022.0.2/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))

load(pathJoin("intel-classic", os.getenv("intel_ver") or "2022.0.2"))
load_any(pathJoin("cray-mpich", os.getenv("cray_mpich_ver") or "7.7.11"),"cray-mpich")
load(pathJoin("hpc-intel", os.getenv("hpc_intel_ver") or "2022.0.2"))
load(pathJoin("hpc-cray-mpich", os.getenv("hpc_cray_mpich_ver") or "7.7.11"))
load(pathJoin("libpng", os.getenv("libpng_ver") or "1.6.37"))

-- Needed at runtime:
load("alps")

load("ufs_common")

setenv("CC","cc")
setenv("FC","ftn")
setenv("CXX","CC")
setenv("CMAKE_Platform","gaea.intel")
