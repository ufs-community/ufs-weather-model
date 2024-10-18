help([[
loads UFS Model prerequisites for Hera/GNU
]])

prepend_path("MODULEPATH", "/contrib/spack-stack/spack-stack-1.8.0/envs/ue-gcc-9.2.0/install/modulefiles/Core")

stack_gnu_ver=os.getenv("stack_gnu_ver") or "9.2.0"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_openmpi_ver=os.getenv("stack_openmpi_ver") or "4.1.6"
load(pathJoin("stack-openmpi", stack_openmpi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "hera.gnu")

whatis("Description: UFS build environment")
