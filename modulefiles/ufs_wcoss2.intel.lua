help([[
Load environment to build UFS on Cactus/Dogwood with Intel compiler
]])

prepend_path("MODULEPATH", "/apps/ops/test/spack-stack-nco-1.9/modulefiles/Core")

load("stack-oneapi")
load("stack-cray-mpich")
load("stack-python")

load("cmake")
load("glibc")
load("intel-oneapi-runtime")
load("zlib-ng")
load(pathJoin("python-venv", "1.0"))
load("py-pyyaml")
load("ufs_common")

setenv("INTEL_COMPILER_TYPE", "RECOMMENDED")
setenv("CC", "icx")
setenv("CXX", "icpx")
setenv("FC", "ifort")
setenv("CMAKE_Platform", "wcoss2")

whatis("Description: UFS build environment")
