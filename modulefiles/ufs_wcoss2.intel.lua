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

local ufs_modules = {
  {["jasper"]      = "2.0.25"},
  {["zlib"]        = "1.2.11"},
  {["libpng"]      = "1.6.37"},
  {["hdf5-C"]      = "1.14.0"},
  {["netcdf-C"]    = "4.9.2"},
  {["pio-C"]       = "2.5.10"},
  {["esmf-C"]      = "8.6.0"},
  {["fms-C"]       = "2023.04"},
  {["bacio"]       = "2.4.1"},
  {["crtm"]        = "2.4.0"},
  {["g2"]          = "3.4.5"},
  {["g2tmpl"]      = "1.10.2"},
  {["ip"]          = "3.3.3"},
  {["sp"]          = "2.3.3"},
  {["w3emc"]       = "2.9.2"},
  {["gftl-shared"] = "1.6.1"},
  {["mapl-C"]      = "2.40.3"},
  {["pnetcdf-C"]   = "1.12.2"},
  {["scotch"]      = "7.0.4"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

setenv("CC", "cc")
setenv("CXX", "CC")
setenv("FC", "ftn")
setenv("CMAKE_Platform", "wcoss2")

whatis("Description: UFS build environment")
