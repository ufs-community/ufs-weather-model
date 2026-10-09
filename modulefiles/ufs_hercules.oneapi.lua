help([[
loads UFS Model prerequisites for Hercules/IntelLLVM (OneAPI)
]])

prepend_path("MODULEPATH", "/apps/contrib/spack-stack/spack-stack-2.1.1/envs/ue-oneapi-2025.3.1/modules/Core")
prepend_path("MODULEPATH", "/apps/contrib/spack-stack/modulefiles")

stack_intel_ver=os.getenv("stack_intel_ver") or "2025.3.1"
load(pathJoin("stack-intel-oneapi-compilers", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.17"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.31.8"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicx")
setenv("CXX", "mpiicpx")
setenv("FC", "mpiifx")
setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_FC", "ifx")

setenv("CMAKE_Platform", "hercules.oneapi")

whatis("Description: UFS build environment")
