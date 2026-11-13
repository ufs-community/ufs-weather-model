help([[
Load environment to compile ufs-weather-model in a container using Intel Classic w/spack-stack v1.9.2. 
This modulefile will switch to OneAPI compilers in the near future. 
]])

prepend_path("MODULEPATH", "/opt/intel/oneapi/2024.2/etc/modulefiles")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.9.2/envs/ufs-wm-env/install/modulefiles/Core")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.9.2/envs/ufs-wm-env/install/modulefiles/intel-oneapi-mpi/2021.13-ux7fmve/gcc/13.3.1")

load("tbb/2021.13")
load("compiler-rt/2024.2.1")
load("compiler/2024.2.1")
load("mkl/2024.2")
load("ifort/2024.2.1")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"

load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

-- load("ufs_common") -- temporarily not use common list of module until moved to containers with spack-stack-2.1.1

local ufs_modules = {
  {["jasper"]          = "2.0.32"},
  {["libpng"]          = "1.6.37"},
  {["hdf5"]            = "1.14.3"},
  {["netcdf-c"]        = "4.9.2"},
  {["netcdf-fortran"]  = "4.6.1"},
  {["parallelio"]      = "2.6.2"},
  {["esmf"]            = "8.8.0"},
  {["fms"]             = "2024.02"},
  {["bacio"]           = "2.4.1"},
  {["crtm"]            = "2.4.0.1"},
  {["g2"]              = "3.5.1"},
  {["g2tmpl"]          = "1.13.0"},
  {["ip"]              = "5.1.0"},
  {["sp"]              = "2.5.0"},
  {["w3emc"]           = "2.10.0"},
  {["gftl-shared"]     = "1.9.0"},
  {["mapl"]            = "2.53.4-esmf-8.8.0"},
  {["scotch"]          = "7.0.4"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicx")
setenv("CXX", "mpiicpx")
setenv("FC", "mpiifort")
setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_F90", "ifort")

setenv("CMAKE_Platform", "container.oneapi")

whatis("Description: UFS build environment")
