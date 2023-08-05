help([[
  This module loads libraries required for building and running UFS Weather Model 
  on the NOAA RDHPC machine Gaea using Intel-2022.0.2
]])

whatis([===[Loads libraries needed for building the UFS Weather Model on Gaea ]===])

unload("intel")
unload("cray-mpich")
unload("cray-python")
unload("darshan")

prepend_path("MODULEPATH", "/lustre/f2/dev/wpo/role.epic/contrib/spack-stack/spack-stack-1.4.1-c4/envs/ufs-pio-2.5.10/install/modulefiles/Core")
prepend_path("MODULEPATH", "/lustre/f2/pdata/esrl/gsd/spack-stack/modulefiles")

stack_intel_ver=os.getenv("stack_intel_ver") or "2022.0.2"
load(pathJoin("stack-intel", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "7.7.20"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.12"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

-- Needed at runtime:
load("alps")

load("ufs_common")

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")
setenv("CMAKE_Platform","gaea.intel")
