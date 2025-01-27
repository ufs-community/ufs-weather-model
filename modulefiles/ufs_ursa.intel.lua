help([[
loads UFS Model prerequisites for Ursa/Intel
]])

--prepend_path("MODULEPATH", "/contrib/spack-stack/envs/1.8.0/ue-oneapi-ifort-2024.2.1/install/modulefiles/Core")
--prepend_path("MODULEPATH", "/contrib/spack-stack/envs/1.8.0/ue-oneapi-ifort-2024.2.1/install/modulefiles/intel-oneapi-mpi/2021.13-eaajhcw/oneapi/2024.2.1")

prepend_path("MODULEPATH", "/contrib/spack-stack/envs/1.6.0/fms-2024.01/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.1.0"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "icc")
setenv("CXX", "icpc")
setenv("FC", "ifort")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")

setenv("I_MPI_CC", "icc")
setenv("I_MPI_CXX", "icpc")
setenv("I_MPI_F90", "ifort")

setenv("CMAKE_Platform", "ursa.intel")

whatis("Description: UFS build environment")
