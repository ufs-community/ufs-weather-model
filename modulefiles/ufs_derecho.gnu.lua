help([[
loads UFS Model prerequisites for Derecho/GNU
]])

setenv("LMOD_TMOD_FIND_FIRST","yes")
prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.9.2/envs/ue-gcc-12.4/install/modulefiles/Core")

unload("ncarcompilers")
stack_gnu_ver=os.getenv("stack_gnu_ver") or "12.4.0"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.29"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.11.7"
load(pathJoin("stack-python", stack_python_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpicc")
setenv("CXX", "mpicxx")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "derecho.gnu")

whatis("Description: UFS build environment")
