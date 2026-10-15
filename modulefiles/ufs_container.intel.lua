help([[
Load environment to compile ufs-weather-model in a container using Intel
]])

prepend_path("MODULEPATH", "/opt/intel/oneapi/2024.2/etc/modulefiles")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.9.2/envs/ufs-wm-env/install/modulefiles/Core")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.9.2/envs/ufs-wm-env/install/modulefiles/intel-oneapi-mpi/2021.13-ux7fmve/gcc/13.3.1")

load("tbb/2021.13")
load("compiler-rt/2024.2.1")
load("compiler/2024.2.1")
load("mkl/2024.2")
load("ifort/2024.2.1")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"

load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicx")
setenv("CXX", "mpiicpx")
setenv("FC", "mpiifort")
setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_F90", "ifort")

setenv("CMAKE_Platform", "container.intel")

whatis("Description: UFS build environment")
