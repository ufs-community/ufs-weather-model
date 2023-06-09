help([[
loads UFS Model prerequisites for Hercules/Intel
]])

prepend_path("MODULEPATH", "/work/noaa/epic-ps/role-epic-ps/spack-stack/spack-stack-1.4.0-hercules/envs/unified-env-v2/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.7.1"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.7.1"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.14"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

setenv("CC", "icc")
setenv("CXX", "icpc")
setenv("FC", "ifort")
setenv("CMAKE_Platform", "hercules.intel")

whatis("Description: UFS build environment")
