whatis("Description: UFS build environment common libraries")

help([[Load UFS Model common libraries]])

local ufs_modules = {
  {["jasper"]          = "2.0.32"},
  {["zlib"]            = "1.2.13"},
  {["libpng"]          = "1.6.37"},
  {["hdf5"]            = "1.14.0"},
  {["netcdf-c"]        = "4.9.2"},
  {["netcdf-fortran"]  = "4.6.0"},
  {["parallelio"]      = "2.5.10"},
  {["esmf"]            = "8.5.0"},
  {["fms"]             = "2023.02.01"},
  {["bacio"]           = "2.4.1"},
  {["crtm"]            = "2.4.0"},
  {["g2"]              = "3.4.5"},
  {["g2tmpl"]          = "1.10.2"},
  {["ip"]              = "4.3.0"},
  {["sp"]              = "2.3.3"},
  {["w3emc"]           = "2.10.0"},
  {["gftl-shared"]     = "1.6.1"},
  {["mapl"]            = "2.40.3-esmf-8.5.0"},
  {["scotch"]          = "7.0.4"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end
