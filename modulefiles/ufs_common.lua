whatis("Description: UFS build environment common libraries")

help([[Load UFS Model common libraries]])

local ufs_modules = {
  {["jasper"]      = "2.0.25"},
  {["zlib"]        = "1.2.11"},
  {["libpng"]      = "1.6.37"},
  {["hdf5"]        = "1.10.6"},
  {["netcdf"]      = "4.7.4"},
  {["pio"]         = "2.5.7"},
  {["esmf"]        = "8.3.0b09"},
  {["fms"]         = "2022.04"},
  {["bacio"]       = "2.4.1"},
  {["crtm"]        = "2.4.0"},
  {["g2"]          = "3.4.5"},
  {["g2tmpl"]      = "1.10.2"},
  {["ip"]          = "3.3.3"},
  {["sp"]          = "2.3.3"},
  {["w3emc"]       = "2.9.2"},
  {["gftl-shared"] = "v1.5.0"},
  {["mapl"]        = "2.22.0-esmf-8.3.0b09"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end
