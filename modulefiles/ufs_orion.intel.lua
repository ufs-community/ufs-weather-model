help([[
loads UFS Model prerequisites for Orion/Intel
]])

prepend_path("MODULEPATH", "/work/noaa/epic/role-epic/spack-stack/spack-stack-1.4.1/envs/ufs-pio-2.5.10/install/modulefiles/Core")
prepend_path("MODULEPATH", "/work/noaa/da/role-da/spack-stack/modulefiles")

stack_intel_ver=os.getenv("stack_intel_ver") or "2022.0.2"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.1"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.7"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "orion.intel")

whatis("Description: UFS build environment")
