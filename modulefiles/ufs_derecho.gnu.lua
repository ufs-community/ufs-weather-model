help([[
loads UFS Model prerequisites for Derecho/GNU
]])

setenv("LMOD_TMOD_FIND_FIRST","yes")
prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/modulefiles")
load("ecflow/5.8.4")
load("mysql/8.0.33")

setenv("LMOD_TMOD_FIND_FIRST","yes")
prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/modulefiles_extra")
prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core")

unload("ncarcompilers")
stack_gnu_ver=os.getenv("stack_gnu_ver") or "12.2.0"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.25"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.26.3"
load(pathJoin("cmake", cmake_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.10.8"
load(pathJoin("stack-python", stack_python_ver))

setenv("CMAKE_Platform", "derecho.gnu")
load("ufs-weather-model-env")

whatis("Description: UFS build environment")
