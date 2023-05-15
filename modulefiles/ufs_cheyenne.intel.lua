help([[
loads UFS Model prerequisites for Cheyenne/Intel
]])

cmake_ver=os.getenv("cmake_ver") or "3.22.0"
load(pathJoin("cmake", cmake_ver))

python_ver=os.getenv("python_ver") or "3.7.9"
load(pathJoin("python", python))

--load programming environment
ncarenv_ver=os.getenv("ncarenv_ver") or "1.3"
load(pathJoin("ncarenv", ncarenv_ver))

intel_ver=os.getenv("intel_ver") or "2022.1"
load(pathJoin("intel", intel_ver))

mpt_ver=os.getenv("mpt_ver") or "2.25"
load(pathJoin("mpt", mpt_ver))

ncarcompilers_ver=os.getenv("ncarcompilers_ver") or "0.5.0"
load(pathJoin("ncarcompilers", ncarcompilers_ver))

unload("netcdf")

prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/hpc-stack/intel2022.1/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2022.1"
load(pathJoin("hpc-intel", hpc_intel_ver))

hpc_mpt_ver=os.getenv("hpc_mpt_ver") or "2.25"
load(pathJoin("hpc-mpt", hpc_mpt_ver))

scotch_ver=os.getenv("scotch_ver") or "7.0.3"
load(pathJoin("scotch", scotch_ver))

load("ufs_common")

setenv("CC", "mpicc")
setenv("CXX", "mpicxx")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "cheyenne.intel")

whatis("Description: UFS build environment")
