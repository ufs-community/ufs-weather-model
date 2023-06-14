help([[
loads UFS Model prerequisites for Hercules/GNU
]])

prepend_path("MODULEPATH", "/work/noaa/epic-ps/role-epic-ps/spack-stack/spack-stack-1.4.0-hercules/envs/unified-env-v2/install/modulefiles/Core")

stack_gnu_ver=os.getenv("stack_gnu_ver") or "11.3.1"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_openmpi_ver=os.getenv("stack_openmpi_ver") or "4.1.4"
load(pathJoin("stack-openmpi", stack_openmpi_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.9.14"
load(pathJoin("stack-python", stack_python_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("ufs_common")

setenv("CC", "gcc")
setenv("CXX", "g++")
setenv("FC", "gfortran")
setenv("CMAKE_Platform", "hercules.gnu")

whatis("Description: UFS build environment")
