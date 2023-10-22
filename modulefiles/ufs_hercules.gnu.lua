help([[
loads UFS Model prerequisites for Hercules/GNU
]])

prepend_path("MODULEPATH", "/work/noaa/epic/role-epic/spack-stack/hercules/spack-stack-1.5.0/envs/unified-env-gnu-mpich/install/modulefiles/Core")
-- for mpich, need:
prepend_path("MODULEPATH", "/work/noaa/epic/role-epic/spack-stack/hercules/modulefiles")

stack_gnu_ver=os.getenv("stack_gnu_ver") or "11.3.1"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_mpich_ver=os.getenv("stack_mpich_ver") or "4.1.2"
load(pathJoin("stack-mpich", stack_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "hercules.gnu")

whatis("Description: UFS build environment")
