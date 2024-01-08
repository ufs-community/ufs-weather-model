help([[
loads UFS Model prerequisites for Hercules/GNU
]])

prepend_path("MODULEPATH", "/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core")
-- for mvapich2, need:
prepend_path("MODULEPATH", "/work/noaa/epic/role-epic/spack-stack/hercules/modulefiles")

stack_gnu_ver=os.getenv("stack_gnu_ver") or "12.2.0"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_mvapich2_ver=os.getenv("stack_mvapich2_ver") or "2.3.7"
load(pathJoin("stack-mvapich2", stack_mvapich2_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "hercules.gnu")

whatis("Description: UFS build environment")
