help([[
Load environment to compile ufs-weather-model in a container using GNU compilers
]])

prepend_path("MODULEPATH", "/opt/modulefiles")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.9.2/envs/ufs-wm-env/install/modulefiles/Core")
stack_gcc_ver=os.getenv("stack_gcc_ver") or "13.3.1"

load(pathJoin("stack-gcc", stack_gcc_ver))
load("stack-openmpi")

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpif90")
setenv("I_MPI_CC", "gcc")
setenv("I_MPI_CXX", "g++")
setenv("I_MPI_F90", "gfortran")

setenv("CMAKE_Platform", "container.gnu")

whatis("Description: UFS build environment")
