help([[
loads NEMS FV3 prerequisites on Cactus and Dogwood
]])

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.1.0"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

intel_ver=os.getenv("intel_ver") or "19.1.3.304"
load(pathJoin("intel", intel_ver))

craype_ver=os.getenv("craype_ver") or "2.7.10"
load(pathJoin("craype", craype_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.9"
load(pathJoin("cray-mpich", cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.20.2"
load(pathJoin("cmake", cmake_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

ip_ver=os.getenv("ip_ver") or "3.3.3"
load(pathJoin("ip", ip_ver))

sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))

w3nco_ver=os.getenv("w3nco_ver") or "2.4.1"
load(pathJoin("w3nco", w3nco_ver))

w3emc_ver=os.getenv("w3emc_ver") or "2.9.2"
load(pathJoin("w3emc", w3emc_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.2"
load(pathJoin("nemsio", nemsio_ver))

g2tmpl_ver=os.getenv("g2tmpl_ver") or "1.9.1"
load(pathJoin("g2tmpl", g2tmpl_ver))

upp_ver=os.getenv("upp_ver") or "8.3.0"
load(pathJoin("upp", upp_ver))

crtm_ver=os.getenv("crtm_ver") or "2.4.0.1"
load(pathJoin("crtm", crtm_ver))

g2_ver=os.getenv("g2_ver") or "3.4.5"
load(pathJoin("g2", g2_ver))

jasper_ver=os.getenv("jasper_ver") or "2.0.25"
load(pathJoin("jasper", jasper_ver))

zlib_ver=os.getenv("zlib_ver") or "1.2.11"
load(pathJoin("zlib", zlib_ver))

libpng_ver=os.getenv("libpng_ver") or "1.6.37"
load(pathJoin("libpng", libpng_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.10.6"
load(pathJoin("hdf5", hdf5_ver))

netcdf_ver=os.getenv("netcdf_ver") or "4.7.4"
load(pathJoin("netcdf", netcdf_ver))

esmf_ver=os.getenv("esmf_ver") or "8.0.1"
load(pathJoin("esmf", esmf_ver))

--## load cmake
setenv("CC", "cc")
setenv("CXX", "CC")
setenv("FC", "ftn")
setenv("CMAKE_Platform", "wcoss2")

whatis("Description: ufs-weather-model build environment")
