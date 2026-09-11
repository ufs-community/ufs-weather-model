help([[
loads UFS Model prerequisites for Ursa/Intel LLVM
]])

source_sh("bash", "/opt/intel/oneapi/setvars.sh")
source_sh("bash", "/opt/intel/oneapi/compiler/2026.1/env/vars.sh")
source_sh("bash", "/opt/intel/oneapi/mpi/2021.18/env/vars.sh")
source_sh("bash", "/opt/intel/oneapi/mkl/2026.1/env/vars.sh")

prepend_path("MODULEPATH", "/work/alexander_richert/stack/spack-stack/envs/nco-sci-oneapi-2026.1.0/modules_flat/Core")

load("python/3.11.15")
load("python-venv/1.0")
load("openblas/0.3.33")

local ufs_modules = {
  {["jasper"]          = "4.2.8"},
  {["libpng"]          = "1.6.55"},
  {["hdf5"]            = "1.14.5"},
  {["netcdf-c"]        = "4.9.2"},
  {["netcdf-fortran"]  = "4.6.1"},
  {["py-pyyaml"]       = "6.0.3"},
  {["parallelio"]      = "2.6.2"},
  {["esmf"]            = "8.8.0"},
  {["fms"]             = "2024.03-gfs-constants"},
  {["bacio"]           = "2.4.1"},
  {["crtm"]            = "3.1.2"},
  {["g2"]              = "3.5.1"},
  {["g2tmpl"]          = "1.17.0"},
  {["ip"]              = "5.4.0"},
  {["sp"]              = "2.5.0"},
  {["w3emc"]           = "2.13.0"},
  {["gftl-shared"]     = "1.12.0"},
  {["mapl"]            = "2.53.4-esmf-8.8.0"},
  {["scotch"]          = "7.0.10"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

-- nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
-- load(pathJoin("nccmp", nccmp_ver))

setenv("CC", "mpiicx")
setenv("CXX", "mpiicpx")
setenv("FC", "mpiifx")
setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_FC", "ifx")

setenv("CMAKE_Platform", "ursa.oneapi")

whatis("Description: UFS build environment")
