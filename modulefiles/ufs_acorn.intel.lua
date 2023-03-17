help([[
Load environment to build UFS on Acorn with Intel compiler
]])

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.1.0"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

intel_ver=os.getenv("intel_ver") or "19.1.3.304"
load(pathJoin("intel", intel_ver))

craype_ver=os.getenv("craype_ver") or "2.7.13"
load(pathJoin("craype", craype_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.9"
load(pathJoin("cray-mpich", cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.20.2"
load(pathJoin("cmake", cmake_ver))

prepend_path("MODULEPATH", "/lfs/h1/emc/nceplibs/noscrub/hpc-stack/libs/hpc-stack/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
hpc_intel_ver=os.getenv("hpc_intel_ver") or "19.1.3.304"
hpc_cray_mpich_ver=os.getenv("hpc_cray_mpich_ver") or "8.1.9"
load(pathJoin("hpc", hpc_ver))
load(pathJoin("hpc-intel", hpc_intel_ver))
load(pathJoin("hpc-cray-mpich", hpc_cray_mpich_ver))

load("ufs_common")

prepend_path("MODULEPATH", "/lfs/h1/emc/nceplibs/noscrub/UPP_IFI/modulefiles")
load("ifi/20230118-intel-19.1.3.304")

setenv("CC", "cc")
setenv("CXX", "CC")
setenv("FC", "ftn")
setenv("CMAKE_Platform", "acorn")

whatis("Description: UFS build environment")
