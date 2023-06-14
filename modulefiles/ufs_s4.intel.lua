help([[
loads UFS Model prerequisites for S4/Intel
]])

prepend_path("MODULEPATH", "/data/prod/jedi/spack-stack/spack-stack-1.4.0/envs/unified-env-v2/install/modulefiles/Core")
prepend_path("MODULEPATH", "/data/prod/jedi/spack-stack/modulefiles")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.5.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.12"
load(pathJoin("stack-python", stack_python_ver))

load("ufs_common")

setenv("CC", "icc")
setenv("CXX", "icpc")
setenv("FC", "ifort")
setenv("CMAKE_Platform", "s4.intel")

whatis("Description: UFS build environment")
