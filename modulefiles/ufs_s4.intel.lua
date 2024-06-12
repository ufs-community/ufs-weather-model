help([[
loads UFS Model prerequisites for S4/Intel
]])

prepend_path("MODULEPATH", "/data/users/dheinzeller/ufs-wm-chunked/spst151-i19/envs/unified-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2019.5.281"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2019.5.281"
load(pathJoin("stack-intel-mpi", stack_impi_ver))

load("ufs_common")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "s4.intel")

whatis("Description: UFS build environment")
