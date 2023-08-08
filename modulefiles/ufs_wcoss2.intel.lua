help([[
loads UFS Model prerequisites on Cactus and Dogwood
]])

-- First, look for libraries in "prod" space
PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.1.0"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

intel_ver=os.getenv("intel_ver") or "19.1.3.304"
load(pathJoin("intel", intel_ver))

craype_ver=os.getenv("craype_ver") or "2.7.13"
load(pathJoin("craype", craype_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.12"
load(pathJoin("cray-mpich", cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.20.2"
load(pathJoin("cmake", cmake_ver))

prepend_path("MODULEPATH", "/apps/test/hpc-stack/i-19.1.3.304__m-8.1.12__h-1.14.0__n-4.9.2__p-2.5.10__e-8.4.2/modulefiles/compiler/intel/19.1.3.304")
prepend_path("MODULEPATH", "/apps/test/hpc-stack/i-19.1.3.304__m-8.1.12__h-1.14.0__n-4.9.2__p-2.5.10__e-8.4.2/modulefiles/mpi/intel/19.1.3.304/cray-mpich/8.1.12")

load("ufs_common")

setenv("CC", "cc")
setenv("CXX", "CC")
setenv("FC", "ftn")
setenv("CMAKE_Platform", "wcoss2")

whatis("Description: UFS build environment")
