help([[
Load environment to compile ufs-weather-model in a container using GNU compilers
]])

prepend_path("MODULEPATH", "/opt/modulefiles")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.9.2/envs/ufs-wm-env/install/modulefiles/Core")
stack_gcc_ver=os.getenv("stack_gcc_ver") or "13.3.1"

load(pathJoin("stack-gcc", stack_gcc_ver))
load("stack-openmpi")

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

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpif90")
setenv("I_MPI_CC", "gcc")
setenv("I_MPI_CXX", "g++")
setenv("I_MPI_F90", "gfortran")

setenv("CMAKE_Platform", "container.gnu")

whatis("Description: UFS build environment")
