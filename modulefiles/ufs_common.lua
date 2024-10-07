whatis("Description: UFS build environment common libraries")

help([[Load UFS Model common libraries]])

local ufs_modules = {
  {["jasper"]          = "2.0.32"},
  {["zlib-ng"]         = "2.1.6"},
  {["libpng"]          = "1.6.37"},
  {["hdf5"]            = "1.14.3"},
  {["netcdf-c"]        = "4.9.2"},
  {["netcdf-fortran"]  = "4.6.1"},
  {["parallelio"]      = "2.6.2"},
  {["esmf"]            = "8.6.1"},
  {["fms"]             = "2024.02"},
  {["bacio"]           = "2.4.1"},
  {["crtm"]            = "2.4.0.1"},
  {["g2"]              = "3.5.1"},
  {["g2tmpl"]          = "1.13.0"},
  {["ip"]              = "5.0.0"},
  {["sp"]              = "2.5.0"},
  {["w3emc"]           = "2.10.0"},
  {["gftl-shared"]     = "1.9.0"},
  {["mapl"]            = "2.46.3-esmf-8.6.1"},
  {["scotch"]          = "7.0.4"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end
