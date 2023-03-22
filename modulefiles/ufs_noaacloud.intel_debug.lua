#%Module
module purge
module unuse /opt/cray/craype/default/modulefiles
module unuse /opt/cray/modulefiles
module use /contrib/spack-stack/miniconda/modulefiles/miniconda
module use /contrib/spack-stack/envs/ufs-wm/install/modulefiles/Core
module load stack-intel/2021.3.0
module load stack-intel/2021.3.0-oneapi-mpi/2021.3.0
module load stack-python/3.9.12
module load cmake/3.23.1

module load jasper/2.0.32 zlib/1.2.13 libpng/1.6.37 hdf5/1.10.6 netcdf-c/4.9.0 netcdf-fortran/4.6.0 parallelio/2.5.9 esmf/8.3.0b09-debug fms/2022.04 bacio/2.4.1 
module load crtm-fix/2.4.0_emc g2/3.4.5 g2tmpl/1.10.2 ip/3.3.3 sp/2.3.3 w3nco/2.4.1 gftl-shared/1.5.0 yafyaml/0.5.1 mapl/2.22.0-esmf-8.3.0b09-debug w3emc/2.9.2
module load crtm/2.4.0

export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
export CMAKE_Platform=noaacloud.intel
