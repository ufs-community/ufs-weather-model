help([[
  This module loads libraries required for building and running UFS Weather Model 
  on the NOAA RDHPC machine Gaea C7 using Intel-2025.3.1.
]])

whatis([===[Loads libraries needed for building the UFS Weather Model on Gaea C7]===])

unload("oneapi/2026.0.0")
unload("cray-mpich/9.2.0")

prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/c7/spack-stack-2.1.1/envs/ue-oneapi-2025.3.1/modules/Core/")

stack_intel_ver=os.getenv("stack_intel_ver") or "2025.3.1"
load(pathJoin("stack-intel-oneapi-compilers", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "9.2.0"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.31.11"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.11.0.0"
load(pathJoin("nccmp", nccmp_ver))

-- unload("darshan-runtime")
-- unload("cray-libsci")

remove_path("LD_LIBRARY_PATH","/autofs/ncrc-svm1_proj/epic/spack-stack/c7/spack-stack-2.1.1/envs/ue-oneapi-2025.3.1/install/intel-oneapi-compilers/2025.3.1/zstd-1.5.7-ogrl5fa/lib")

setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_FC", "ifx")

setenv("CC","icx")
setenv("CXX","icpx")
setenv("FC","ifx")
setenv("CMAKE_Platform","gaeac7.oneapi")
