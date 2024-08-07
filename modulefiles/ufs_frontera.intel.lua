help([[
loads UFS Model prerequisites for Frontera/Intel
]])

prepend_path("MODULEPATH", "/work2/06146/tg854455/frontera/spack-stack/modulefiles")
load("ecflow/5.8.4")

prepend_path("MODULEPATH", "/work2/01118/tg803972/frontera/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "19.1.1.217"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2020.4.304"
load(pathJoin("stack-intel-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.24.2"
load(pathJoin("cmake", cmake_ver))
--load("cmake/3.24.2")

load("ufs_common")

stack_python_ver=os.getenv("stack_python_ver") or "3.10.13"
load(pathJoin("stack-python", stack_python_ver))

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "frontera.intel")

whatis("Description: UFS build environment")
