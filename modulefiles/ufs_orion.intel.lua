help([[
loads UFS Model prerequisites for Orion/Intel
]])

prepend_path("MODULEPATH", "/work/noaa/epic-ps/role-epic-ps/spack-stack/spack-stack-1.3.1/envs/unified-env/install/modulefiles/Core")
prepend_path("MODULEPATH", "/work/noaa/da/role-da/spack-stack/modulefiles")

stack_intel_ver=os.getenv("stack_intel_ver") or "2022.0.2"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.1"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.7"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.22.1"
load(pathJoin("cmake", cmake_ver))

scotch_ver=os.getenv("scotch_ver") or "7.0.3"
load(pathJoin("scotch", scotch_ver))

load("ufs_common")

setenv("CC", "icc")
setenv("CXX", "icpc")
setenv("FC", "ifort")
setenv("CMAKE_Platform", "orion.intel")

whatis("Description: UFS build environment")
