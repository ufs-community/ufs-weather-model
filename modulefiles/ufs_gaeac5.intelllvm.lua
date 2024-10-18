help([[
  This module loads libraries required for building and running UFS Weather Model 
  on the NOAA RDHPC machine Gaea C5 using Intel-2023.2.0.
]])

whatis([===[Loads libraries needed for building the UFS Weather Model on Gaea ]===])

prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/spack-stack-1.6.0/envs/upp-addon-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2023.2.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.28"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.10.13"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

unload("darshan-runtime")
unload("cray-libsci")

unload("intel-classic/2023.2.0")
load("intel-oneapi/2023.2.0")

setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_F90", "ifx")

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")
setenv("CMAKE_Platform","gaeac5.intelllvm")
