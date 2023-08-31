help([[
loads UFS Model prerequisites for NOAA Parallelworks/Intel
]])


setenv("LMOD_TMOD_FIND_FIRST","yes")
prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/modulefiles")
load("ecflow/5.8.4")
load("mysql/8.0.33")

setenv("LMOD_TMOD_FIND_FIRST","yes")
--prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/modulefiles_extra")
--prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/spack-stack/spack-stack-dev-20230825/envs/unified-env/install/modulefiles/Core")
--prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/spack-stack/spack-stack-dev-20230814/envs/unified-en2/install/modulefiles/Core")
--prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/spack-stack/spack-stack-dev-20230814/envs/unified-env/install/modulefiles/Core")
--prepend_path("MODULEPATH", "/lustre/desc1/scratch/mpotts/spack-stack/spack-stack-dev-20230814/envs/ufs-pio-2.5.10/install/modulefiles/Core")
prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/spack-stack/spack-stack-1.4.1/envs/unified-dev/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.8.0"
--stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.8.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.26.3"
load(pathJoin("cmake", cmake_ver))

stack_python_ver=os.getenv("stack_python_ver") or "3.10.10"
load(pathJoin("stack-python", stack_python_ver))


setenv("CMAKE_Platform", "derecho.intel")
load("ufs-weather-model-env")
--prepend_path("PATH","/opt/cray/pe/pals/1.2.11/bin")
--setenv("CC", "cc")
--setenv("CXX", "CC")
--setenv("FC", "ftn")

whatis("Description: UFS build environment")
