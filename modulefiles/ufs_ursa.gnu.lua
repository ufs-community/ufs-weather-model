help([[
loads UFS Model prerequisites for Ursa/GNU
]])

prepend_path("MODULEPATH", "/collab1/data/Ratko.Vasic/spack-stack-1.6.0/envs/fms-2024.01-gnu/install/modulefiles/Core")

stack_gnu_ver=os.getenv("stack_gnu_ver") or "12.4.0"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_openmpi_ver=os.getenv("stack_openmpi_ver") or "4.1.6"
load(pathJoin("stack-openmpi", stack_openmpi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.1.0"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpif90")

setenv("CMAKE_Platform", "ursa.gnu")

whatis("Description: UFS build environment")
