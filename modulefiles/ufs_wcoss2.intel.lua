help([[
loads UFS Model prerequisites on Cactus and Dogwood
]])

-- First, look for libraries in "prod" space
PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.1.0"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

intel_ver=os.getenv("intel_ver") or "19.1.3.304"
load(pathJoin("intel", intel_ver))

craype_ver=os.getenv("craype_ver") or "2.7.13"
load(pathJoin("craype", craype_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.12"
load(pathJoin("cray-mpich", cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.20.2"
load(pathJoin("cmake", cmake_ver))

prepend_path("MODULEPATH", "/apps/test/hpc-stack/i-19.1.3.304__m-8.1.12__h-1.14.0__n-4.9.2__p-2.5.10__e-8.4.2/modulefiles/compiler/intel/19.1.3.304")
prepend_path("MODULEPATH", "/apps/test/hpc-stack/i-19.1.3.304__m-8.1.12__h-1.14.0__n-4.9.2__p-2.5.10__e-8.4.2/modulefiles/mpi/intel/19.1.3.304/cray-mpich/8.1.12")

local ufs_modules = {
  {["jasper"]          = "2.0.25"},
  {["zlib"]            = "1.2.11"},
  {["libpng"]          = "1.6.37"},
  {["hdf5"]            = "1.14.0"},
  {["netcdf"]          = "4.9.2"},
  {["pio"]             = "2.5.10"},
  {["esmf"]            = "8.4.2"},
  {["fms"]             = "2023.01"},
  {["bacio"]           = "2.4.1"},
  {["crtm"]            = "2.4.0"},
  {["g2"]              = "3.4.5"},
  {["g2tmpl"]          = "1.10.2"},
  {["ip"]              = "3.3.3"},
  {["sp"]              = "2.3.3"},
  {["w3emc"]           = "2.9.2"},
  {["gftl-shared"]     = "v1.5.0"},
  {["mapl"]            = "2.35.2-esmf-8.4.2"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

prepend_path("MODULEPATH", "/apps/prod/lmodules/INTEL_cray_mpich/19.1.3.304/cray-mpich/8.1.9")
scotch_ver=os.getenv("scotch_ver") or "7.0.4"
load(pathJoin("scotch",scotch_ver))

setenv("CC", "cc")
setenv("CXX", "CC")
setenv("FC", "ftn")
setenv("CMAKE_Platform", "wcoss2")

whatis("Description: UFS build environment")
