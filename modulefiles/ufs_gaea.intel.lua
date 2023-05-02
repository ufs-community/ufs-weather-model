help([[
  This module loads libraries required for building and running UFS Weather Model
  on the NOAA RDHPC machine Gaea using Intel-2021.3.0
]])

whatis([===[Loads libraries needed for building the UFS Weather Model on Gaea ]===])

prepend_path("MODULEPATH", "/lustre/f2/dev/wpo/role.epic/contrib/spack-stack/spack-stack-1.3.0/envs/unified-env/install/modulefiles/Core")
prepend_path("MODULEPATH", "/lustre/f2/pdata/esrl/gsd/spack-stack/modulefiles")
prepend_path("MODULEPATH", "/lustre/f2/dev/role.epic/contrib/modulefiles")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.3.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "7.7.11"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.12"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

-- needed for WW3 build
load(pathJoin("gcc", os.getenv("gcc_ver") or "8.3.0"))
-- Needed at runtime:
load("alps")
load("rocoto")

load("ufs_common")

setenv("CC","cc")
setenv("FC","ftn")
setenv("CXX","CC")
setenv("CMAKE_C_COMPILER","cc")
setenv("CMAKE_CXX_COMPILER","CC")
setenv("CMAKE_Fortran_COMPILER","ftn")
setenv("CMAKE_Platform","gaea.intel")
