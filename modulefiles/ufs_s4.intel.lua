help([[
loads UFS Model prerequisites for S4/Intel
]])

prepend_path("MODULEPATH", "/data/prod/jedi/spack-stack/spack-stack-1.6.0/envs/upp-addon-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.5.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.10.13"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "s4.intel")

whatis("Description: UFS build environment")
