help([[
loads UFS Model prerequisites for WCOSS3/oneAPI
]])

-- Activate the Intel oneAPI compiler environment.
source_sh("bash", "/opt/intel/oneapi/setvars.sh")
source_sh("bash", "/opt/intel/oneapi/compiler/2026.1/env/vars.sh")
source_sh("bash", "/opt/intel/oneapi/mpi/2021.18/env/vars.sh")
source_sh("bash", "/opt/intel/oneapi/mkl/2026.1/env/vars.sh")

-- Use the WCOSS3 spack-stack module tree.
prepend_path("MODULEPATH", "/lfs/work/alexander_richert/stack/spack-stack/envs/nco-sci-oneapi-2026.1.0/modules_flat/Core")

-- Load the common UFS dependency stack.
load("python")
load("python-venv/1.0")
load("py-pyyaml/6.0.3")
load("ufs_common_wcoss3")

-- Additional WCOSS3 stack modules from the flat Core tree.
load("cmake")
load("python")
load("netcdf-c")
load("netcdf-fortran")
load("parallelio")
load("esmf")
load("fms")
load("bacio")
load("crtm")
load("g2")
load("g2tmpl")
load("ip")
load("w3emc")
load("gftl-shared")
load("mapl")
load("scotch")

setenv("CC", "mpiicx")
setenv("CXX", "mpiicpx")
setenv("FC", "mpiifx")
setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_FC", "ifx")
setenv("CMAKE_Platform", "wcoss3.oneapi")

whatis("Description: UFS build environment")
