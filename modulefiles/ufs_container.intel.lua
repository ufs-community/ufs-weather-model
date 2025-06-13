help([[
Load environment to compile ufs-weather-model in a container using Intel
]])

prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.6.0/envs/fms-2024.01/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
stack_impi_ver=os.getenv("stack_impi_ver") or "2021.9.0"

load("gnu")
load(pathJoin("stack-intel", stack_intel_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
unload("gnu")

load("ufs_common")

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CMAKE_Platform", "container.intel")

whatis("Description: UFS build environment")
