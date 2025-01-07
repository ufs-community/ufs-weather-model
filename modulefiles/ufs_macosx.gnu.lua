help([[
loads UFS Model prerequisites for MacOS clang/gcc ("gnu")
]])

prepend_path("MODULEPATH", "/Users/username/spack-stack/spack-stack-1.8.0/envs/ufs-srw-env/install/modulefiles/Core")

stack_gnu_ver=os.getenv("stack_apple_clang_ver") or "15.0.0"
load(pathJoin("stack-apple-clang", stack_gnu_ver))

stack_openmpi_ver=os.getenv("stack_openmpi_ver") or "5.0.3"
load(pathJoin("stack-openmpi", stack_openmpi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

local ufs_modules = {
  {["jasper"]          = "2.0.32"},
  {["zlib"]            = "1.2.13"},
  {["libpng"]          = "1.6.37"},
  {["hdf5"]            = "1.14.0"},
  {["netcdf-c"]        = "4.9.2"},
  {["netcdf-fortran"]  = "4.6.1"},
  {["parallelio"]      = "2.6.2"},
  {["esmf"]            = "8.6.0"},
  {["fms"]             = "2024.02"},
  {["bacio"]           = "2.4.1"},
  {["crtm"]            = "2.4.0.1"},
  {["g2"]              = "3.5.1"},
  {["g2tmpl"]          = "1.13.0"},
  {["ip"]              = "5.0.0"},
  {["sp"]              = "2.5.0"},
  {["w3emc"]           = "2.10.0"},
  {["gftl-shared"]     = "1.9.0"},
  {["mapl"]            = "2.40.3-esmf-8.6.0"},
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
setenv("CXX", "mpicxx")
setenv("FC", "mpifort")
setenv("CMAKE_Platform", "macosx.gnu")
setenv("CMAKE_Fortran_COMPILER_ID", "GNU")

osx_sysroot=os.getenv("OSX_SYSROOT")
setenv("CMAKE_OSX_SYSROOT","OSX_SYSROOT")

setenv("CFLAGS"," -Wno-implicit-function-declaration ")

if mode() == "load" then
  LmodMsgRaw([===[
   Please export these env. variables after the module is successfully loaded:
       > export LDFLAGS+=" -L${libjpeg_turbo_ROOT}/lib -ljpeg -Wl,-rpath,$libjpeg_turbo_ROOT}/lib -L${jasper_ROOT}/lib -ljasper -Wl,-rpath,${jasper_ROOT}/lib -L${libpng_ROOT}/lib -lpng -Wl,-rpath,${libpng_ROOT}/lib "
  ]===])
end
whatis("Description: UFS build environment")
