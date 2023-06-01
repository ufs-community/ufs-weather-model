help([[
loads UFS Model prerequisites for Hera/GNU
]])

prepend_path("MODULEPATH", "/contrib/sutils/modulefiles")
load("sutils")

cmake_ver=os.getenv("cmake_ver") or "3.20.1"
load(pathJoin("cmake", cmake_ver))

gnu_ver=os.getenv("gnu_ver") or "9.2.0"
load(pathJoin("gnu", gnu_ver))

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/gnu-9.2/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_gnu_ver=os.getenv("hpc_gnu_ver") or "9.2"
load(pathJoin("hpc-gnu", hpc_gnu_ver))

hpc_mpich_ver=os.getenv("hpc_mpich_ver") or "3.3.2"
load(pathJoin("hpc-mpich", hpc_mpich_ver))

scotch_ver=os.getenv("scotch_ver") or "7.0.3"
load(pathJoin("scotch", scotch_ver))

load("ufs_common")

setenv("CC", "mpicc")
setenv("CXX", "mpicxx")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "hera.gnu")

whatis("Description: UFS build environment")
