help([[
loads UFS Model prerequisites for Cheyenne/GNU
]])

-- Fix module setup on Cheyenne
remove_path ("MODULEPATH", "/glade/u/apps/ch/modulefiles/default/compilers")
setenv("MODULEPATH_ROOT", "/glade/work/jedipara/cheyenne/spack-stack/modulefiles")

-- Load spack-stack modules
prepend_path("MODULEPATH", "/glade/work/jedipara/cheyenne/spack-stack/modulefiles/compilers")
prepend_path("MODULEPATH", "/glade/work/jedipara/cheyenne/spack-stack/modulefiles/misc")
prepend_path("MODULEPATH", "/glade/work/jedipara/cheyenne/spack-stack/spack-stack-v1/envs/skylab-2.0.0-plus-ufs-gnu-10.1.0/install/modulefiles/Core")

load("miniconda/3.9.12")
load("ecflow/5.8.4")

load("stack-gcc/10.1.0")
load("stack-openmpi/4.1.1")
load("stack-python/3.9.12")

mpt_ver=os.getenv("mpt_ver") or "2.25"
load(pathJoin("mpt", mpt_ver))

ncarcompilers_ver=os.getenv("ncarcompilers_ver") or "0.5.0"
load(pathJoin("ncarcompilers", ncarcompilers_ver))

unload("netcdf")

prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/hpc-stack/gnu10.1.0_ncdf492/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_gnu_ver=os.getenv("hpc_gnu_ver") or "10.1.0"
load(pathJoin("hpc-gnu", hpc_gnu_ver))

hpc_mpt_ver=os.getenv("hpc_mpt_ver") or "2.25"
load(pathJoin("hpc-mpt", hpc_mpt_ver))

load("jasper/2.0.25")
load("zlib/1.2.11")
load("hdf5/1.14.0")
load("netcdf/4.9.2")
load("pio/2.5.10")
load("esmf/8.4.2")
load("fms/2023.01")
load("bacio/2.4.1")
load("crtm/2.4.0")
load("g2/3.4.5")
load("g2tmpl/1.10.2")
load("ip/3.3.3")
load("sp/2.3.3")
load("w3emc/2.9.2")
load("gftl-shared/v1.5.0")
load("mapl/2.35.2-esmf-8.4.2")
load("scotch/7.0.3")


setenv("CC", "mpicc")
setenv("CXX", "mpicxx")
setenv("FC", "mpif90")
setenv("CMAKE_Platform", "cheyenne.gnu")

whatis("Description: UFS build environment")
