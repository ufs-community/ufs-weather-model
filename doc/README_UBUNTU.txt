# Dom Heinzeller (dom.heinzeller@noaa.gov), 11/05/2018

In order to build and run the FV3 trunk (November 2018) with possible CCPP extensions by GMTB on Ubuntu Linux,
the following installation steps are recommended. The version numbers correspond to the default versions in
August 2018 and will change to newer versions in the future. Unless problems occur during the manual builds in
step 4, these differences can be ignored. It is also assumed that the bash shell is used in the following.

1. Install Ubuntu Linux 18.04.1 LTS or start up cloud instance (e.g. Amazon AWS) with Ubuntu Linux 18.04.1 LTS

2. Install standard packages as root (sudo su)

    apt-get update

    # Install gcc/gfortran/g++ 7.3.0
    apt-get install gfortran libgfortran-7-dev g++ libstdc++-7-dev gdb

    # Install make and cmake etc.
    apt-get install make cmake m4 ksh

    # Install Python
    apt-get install python

    # Install libxml2-dev
    apt-get install libxml2-dev

    # "Install" gmake
    ln -s /usr/bin/make /usr/bin/gmake

3. Install netCDF development libraries manually as root (with all the bells and whistles)

    export LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:$LD_LIBRARY_PATH"

    mkdir /home/ubuntu/src && cd /home/ubuntu/src

    # Download the following src files from the web to /home/ubuntu/src
    esmf_7_1_0r_src.tar.gz
    hdf5-1.8.21.tar.gz
    mpich-3.2.1.tar.gz
    netcdf-4.6.1.tar.gz
    netcdf-fortran-4.4.4.tar.gz
    parallel-netcdf-1.8.1.tar.gz
    szip-2.1.1.tar.gz
    zlib-1.2.11.tar.gz

    # zlib-1.2.11
    tar -xvf zlib-1.2.11.tar.gz 
    cd zlib-1.2.11/
    CC=mpicc \
    CFLAGS="-fPIC" \
    ./configure \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr zlib-1.2.11

    # szip-2.1.1
    gunzip szip-2.1.1.tar.gz 
    tar -xvf szip-2.1.1.tar 
    gzip szip-2.1.1.tar 
    cd szip-2.1.1/
    CC=mpicc \
    CFLAGS="-fPIC" \
    ./configure \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr szip-2.1.1

    # hdf5-1.8.21
    gunzip hdf5-1.8.21.tar.gz
    tar -xvzf hdf5-1.8.21.tar
    gzip hdf5-1.8.21.tar
    cd hdf5-1.8.21/
    CC=mpicc \
    CFLAGS="-fPIC" \
    FC=mpif90 \
    FCFLAGS="-fPIC" \
    CXX=mpicxx \
    CXXFLAGS="-fPIC" \
    ./configure \
    --enable-parallel \
    --enable-production \
    --enable-fortran \
    --enable-shared \
    --enable-static \
    --with-szlib=/usr/local \
    --with-zlib=/usr/local \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr hdf5-1.8.21

    # parallel-netcdf-1.8.1
    tar -xvf parallel-netcdf-1.8.1.tar.gz
    cd parallel-netcdf-1.8.1
    CC=gcc \
    CFLAGS="-g -O2 -fPIC" \
    CXX=g++ \
    CXXFLAGS="-g -O2 -fPIC" \
    FC=gfortran \
    FCFLAGS="-g -O2 -fPIC" \
    F77FLAGS="-fPIC -g -O2" \
    MPICC=mpicc \
    MPICXX=mpicxx \
    MPIF77=mpif77 \
    MPIF90=mpif90 \
    ./configure \
    --prefix=/usr/local \
    --enable-fortran \
    --enable-largefile \
    --disable-large-file-test 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr parallel-netcdf-1.8.1

    # netcdf-4.6.1
    tar -xvf netcdf-4.6.1.tar.gz 
    cd netcdf-4.6.1/
    CC=mpicc \
    CFLAGS="-fPIC" \
    ./configure \
    --prefix=/usr/local \
    --enable-netcdf-4 \
    --enable-parallel-tests \
    --disable-dap \
    --enable-cdf5 \
    --disable-large-file-tests \
    --enable-shared \
    --enable-static \
    --enable-parallel4 \
    --enable-pnetcdf 2>&1 | tee log.config
    make 2>&1 | tee log.make
    make install 2>&1 | tee log.install
    cd ..
    rm -fr netcdf-4.6.1

    # netcdf-fortran-4.4.4
    tar -xvf netcdf-fortran-4.4.4.tar.gz 
    cd netcdf-fortran-4.4.4/
    CC=mpicc \
    CFLAGS="-I/usr/local/include -fPIC" \
    FC=mpif90 \
    FCLAGS="-I/usr/local/include -fPIC" \
    F77=mpif77 \
    FFLAGS="-I/usr/local/include -fPIC" \
    LDFLAGS="-fPIC -L/usr/local/lib -lnetcdf -lhdf5 -lhdf5_hl -lsz -lz" \
    ./configure \
    --enable-parallel-tests \
    --prefix=/usr/local 2>&1 | tee log.config
    make 2>&1 | tee log.make
    patch -p0 < ../netcdf-fortran-4.4.4.patch
    make check 2>&1 | tee log.check
    make install 2>&1 | tee log.install
    cd ..
    rm -fr netcdf-fortran-4.4.4

    # NCEP libraries
    git clone https://github.com/climbfuji/NCEPlibs.git
    mv NCEPlibs NCEPlibs-20181105
    cd NCEPlibs-20181105
    mkdir /usr/local/NCEPlibs-gnu-20181105
    ./make_ncep_libs.sh -s linux -c gnu -d /usr/local/NCEPlibs-gnu-20181105 -o 1
    cd ..
    rm -fr NCEPlibs-20181105

    # Download esmf-7.1.0r to /usr/local/src
    tar -xvzf esmf_7_1_0r_src.tar.gz
    cd esmf
    export NETCDF=/usr/local
    export ESMF_DIR=`pwd`
    export ESMF_INSTALL_PREFIX=/usr/local/esmf-7.1.0r
    export ESMF_CXXCOMPILER=mpicxx
    export ESMF_CXXLINKER=mpicxx
    export ESMF_F90COMPILER=mpif90
    export ESMF_F90LINKER=mpif90
    export ESMF_COMM=mpich3
    export ESMF_MPIRUN=mpiexec
    export ESMF_NETCDF=1
    export ESMF_NETCDF_INCLUDE=$NETCDF/include
    export ESMF_NETCDF_LIBPATH=$NETCDF/lib
    export ESMF_NETCDF=split
    export ESMF_SL_LIBLIBS="-L/usr/local/lib -lmpichcxx -lmpichf90 -lmpich"
    make info 2>&1 | tee log.info
    make 2>&1 | tee log.make
    make check 2>&1 | tee log.check # this takes forever
    make install 2>&1 | tee log.install
    export -n NETCDF
    export -n ESMF_DIR
    export -n ESMF_INSTALL_PREFIX
    export -n ESMF_CXXCOMPILER
    export -n ESMF_CXXLINKER
    export -n ESMF_F90COMPILER
    export -n ESMF_F90LINKER
    export -n ESMF_COMM
    export -n ESMF_MPIRUN
    export -n ESMF_NETCDF
    export -n ESMF_NETCDF_INCLUDE
    export -n ESMF_NETCDF_LIBPATH
    export -n ESMF_NETCDF
    export -n ESMF_SL_LIBLIBS
    cd ..
    rm -fr esmf

4. Compile NEMSfv3gfs as normal user

    mkdir ~/scratch

    # rsync NEMSfv3gfs from GMTB's Github repository to ~/scratch/NEMSfv3gfs

    cd ~/scratch/NEMSfv3gfs/tests

    ./compile.sh $PWD/../FV3 linux.gnu 'CCPP=N'          2>&1 | tee log.compile # without CCPP
    ./compile.sh $PWD/../FV3 linux.gnu 'CCPP=Y'          2>&1 | tee log.compile # with CCPP, dynamic mode

5. Set up the run directory using the template on Theia or Cheyenne at some location on your machine:

    a) copy the contents of the run directory templates to where you want to run the model, change to this directory
       (these folders are read-only, i.e. users might have to add the write-flag after copying/rsyncing them)

        theia:    /scratch4/BMC/gmtb/Dom.Heinzeller/linux_rundirs/C96_trunk_20180831/gnu/
        cheyenne: /glade/p/ral/jntp/GMTB/NEMSfv3gfs/linux_rundirs/C96_trunk_20180831/gnu/

    b) edit run_linux_no_ccpp.sh/run_linux_ccpp.sh in change the variable FV3_BUILD_DIR to the top-level directory of your FV3-build

    c) edit model_configure to adjust the number of MPI tasks used (PE_MEMBER01 and ncores_per_node)

    d) edit input_no_ccpp.nml/input_ccpp.nml to adjust the splitting of the tile across the MPI tasks (parameter layout)

    e) start up the model using "OMP_NUM_THREADS=X ./run_linux_no_ccpp.sh 2>&1 | tee run_linux.log" (X=1,2,..; identical for  (run_linux_ccpp.sh)
