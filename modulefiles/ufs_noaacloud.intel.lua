help([[
loads UFS Model prerequisites for NOAA Parallelworks/Intel
]])


prepend_path("MODULEPATH", "/contrib/EPIC/spack-stack/spack-stack-1.4.1/envs/ufs-pio-2.5.10/install/modulefiles/Core")
prepend_path("MODULEPATH", "/contrib/spack-stack/modulefiles/core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.3.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.3.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.12"
load(pathJoin("stack-python", stack_python_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "noaacloud.intel")

whatis("Description: UFS build environment")
