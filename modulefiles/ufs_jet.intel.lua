help([[
loads UFS Model prerequisites for Jet/Intel
]])

-- When Jet switches from CentOS to Rocky (June 2024), replace line with correct spack-stack
-- If you want to use Rocky OS now on Jet:
-- 1. login to one of Rocky Jet nodes: fe5 - fe8
-- 2. uncomment line to load correct spack-stack (comment/delete line to old spack-stack)
-- 3. in file  ../tests/rt.sh correct ecflow
-- prepend_path("MODULEPATH", "/mnt/lfs4/HFIP/hfv3gfs/role.epic/spack-stack/spack-stack-1.5.1/envs/unified-env-rocky8/install/modulefiles/Core")
prepend_path("MODULEPATH", "/mnt/lfs4/HFIP/hfv3gfs/role.epic/spack-stack/spack-stack-1.5.1/envs/unified-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.5.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.1"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
setenv("CMAKE_Platform", "jet.intel")

whatis("Description: UFS build environment")
