help([[
loads UFS Model prerequisites on Cactus and Dogwood
]])

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.1.0"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))
--module load PrgEnv-intel/8.1.0

intel_ver=os.getenv("intel_ver") or "19.1.3.304"
load(pathJoin("intel", intel_ver))
--module load intel/19.1.3.304

craype_ver=os.getenv("craype_ver") or "2.7.13"
load(pathJoin("craype", craype_ver))
--module load craype/2.7.13

cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.7"
load(pathJoin("cray-mpich", cray_mpich_ver))
--module load cray-mpich/8.1.7

cmake_ver=os.getenv("cmake_ver") or "3.20.2"
load(pathJoin("cmake", cmake_ver))
--module load cmake/3.20.2

--#module use /apps/ops/para/libs/modulefiles/stack
--#module load hpc/1.2.0
--#module load hpc-intel/19.1.3.304
--#module load hpc-cray-mpich/8.1.7
--#module load ufs_common

setenv("HPC_OPT", "/apps/ops/para/libs")
prepend_path("MODULEPATH", "/apps/ops/para/libs/modulefiles/compiler/intel/19.1.3.304")
prepend_path("MODULEPATH", "/apps/ops/para/libs/modulefiles/mpi/intel/19.1.3.304/cray-mpich/8.1.7")

jasper_ver=os.getenv("jasper_ver") or "2.0.25"
load(pathJoin("jasper", jasper_ver))
--module load jasper/2.0.25

zlib_ver=os.getenv("zlib_ver") or "1.2.11"
load(pathJoin("zlib", zlib_ver))
--module load zlib/1.2.11

libpng_ver=os.getenv("libpng_ver") or "1.6.37"
load(pathJoin("libpng", libpng_ver))
--module load libpng/1.6.37

hdf5_ver=os.getenv("hdf5_ver") or "1.10.6"
load(pathJoin("hdf5", hdf5_ver))
--module load hdf5/1.10.6

netcdf_ver=os.getenv("netcdf_ver") or "4.7.4"
load(pathJoin("netcdf", netcdf_ver))
--module load netcdf/4.7.4

pio_ver=os.getenv("pio_ver") or "2.5.2"
load(pathJoin("pio", pio_ver))
--module load pio/2.5.2

esmf_ver=os.getenv("esmf_ver") or "8.3.0b09"
load(pathJoin("esmf", esmf_ver))
--module load esmf/8.3.0b09

fms_ver=os.getenv("fms_ver") or "2022.01"
load(pathJoin("fms", fms_ver))
--module load fms/2022.01

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))
--module load bacio/2.4.1

crtm_ver=os.getenv("crtm_ver") or "2.3.0"
load(pathJoin("crtm", crtm_ver))
--module load crtm/2.3.0

g2_ver=os.getenv("g2_ver") or "3.4.5"
load(pathJoin("g2", g2_ver))
--module load g2/3.4.5

g2tmpl_ver=os.getenv("g2tmpl_ver") or "1.10.0"
load(pathJoin("g2tmpl", g2tmpl_ver))
--module load g2tmpl/1.10.0

ip_ver=os.getenv("ip_ver") or "3.3.3"
load(pathJoin("ip", ip_ver))
--module load ip/3.3.3

sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))
--module load sp/2.3.3

w3emc_ver=os.getenv("w3emc_ver") or "2.9.2"
load(pathJoin("w3emc", w3emc_ver))
--module load w3emc/2.9.2

--#module load gftl-shared/v1.3.3
--#module load yafyaml/v0.5.1
--#module load mapl/2.11.0-esmf-8.3.0b09

setenv("CC", "cc")
setenv("CXX", "CC")
setenv("FC", "ftn")
setenv("CMAKE_Platform", "wcoss2")

whatis("Description: UFS build environment")