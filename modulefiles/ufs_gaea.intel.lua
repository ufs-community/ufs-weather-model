help([[
loads UFS Model prerequisites for Gaea/Intel
]])

prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2023.1.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.25"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")
setenv("CMAKE_Platform","gaea.intel")
-- failed with:Could not find compiler set in environment variable CC:
--setenv("CC", "mpiicc")
--setenv("CXX", "mpiicpc")
--setenv("FC", "mpiifort")
--setenv("CMAKE_Platform", "gaea.intel")

whatis("Description: UFS build environment")
