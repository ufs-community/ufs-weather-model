whatis("Description: UFS build environment common libraries")

help([[Load UFS Model common libraries]])

local ufs_modules = {
  {["jasper"]          = "4.2.4"},
  {["libpng"]          = "1.6.37"},
  {["hdf5"]            = "1.14.5"},
  {["netcdf-c"]        = "4.9.2"},
  {["netcdf-fortran"]  = "4.6.1"},
  {["parallelio"]      = "2.6.2"},
  {["esmf"]            = "8.8.0"},
  {["fms"]             = "2025.03-gfs-constants"},
  {["bacio"]           = "2.6.0"},
  {["crtm"]            = "3.1.3"},
  {["g2"]              = "3.5.1"},
  {["g2tmpl"]          = "1.17.0"},
  {["ip"]              = "5.4.0"},
  {["w3emc"]           = "2.13.0"},
  {["gftl-shared"]     = "1.11.0"},
  {["mapl"]            = "2.53.4-esmf-8.8.0"},
  {["scotch"]          = "7.0.10"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end
